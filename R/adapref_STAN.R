#' Adaptive preferential sampling method_HMC using stan
#'
#' @param data list containing object or list containing vectors of coalescent
#'   times \code{coal_times}, sampling times \code{samp_times},  number
#'   sampled per sampling time \code{n_sampled}, \code{grid} grid points.\#' @param zeta1 precision parameter field on effective pop size
#' @param prior.c horseshoe or normal (gaussian) field on effective pop size
#' @param prior.s horseshoe or normal (gaussian) field on beta
#' @param order.s order field on eff pop
#' @param order.s order field on beta
#' @param zeta1 precision parameter field on effective pop size
#' @param zeta2 precision parameter field on beta
#' @param par vector parameter to save
#' @param chain number of chains to run
#' @param warmup burn in
#' @param thin thinning
#' @param iter  number of interations
#' @param control list control variables, typicall list(adapt_delta=0.98, max_treedepth=12)
#' @param save.loglik save log like if one wants to compute  = T
#' @param derivative logical whether to calculate estimates of the
#'   log-derivative.
#' @param forward logical whether to use the finite difference approximations of
#'   the log-derivative as a forward or backward derivative.
#'
#' @return Effective pop size and beta posterior- Stan output, to be passed into extract_out
#'   in order to have format for plotting
#' @export


adapref_sampling <- function(prior.c="normal", prior.s="horseshoe", order.c=1,order.s=1, zeta.c=0.01,zeta.s=0.01, fit=NA, data, pars=NA, chains=4, iter=2000, warmup=floor(iter/2), thin=1, control=list(adapt_delta=0.95, max_treedepth=12), save.loglik=FALSE, ...)  {

  # check for rstan
  if (!requireNamespace("rstan", quietly = TRUE)) {
    stop("Package 'rstan' needed for this function to work. Please install it.", call. = FALSE)
  }
  if (!(prior.c %in% c("horseshoe", "normal"))) stop("Must specify prior of type 'normal' or 'horseshoe'.")
  if (!(prior.s %in% c("horseshoe", "normal"))) stop("Must specify prior of type 'normal' or 'horseshoe'.")
  if (!(order.c %in% c(1,2))) stop("Model must be of order 1, 2, or 3.")
  if (zeta.c <= 0) stop("zeta must be > 0.")
  if (is.null(data) | class(data)!="list") stop("Must specify input data set as a list.")

  stlist <- list(...)


  # check if necessary elements are in coalescent data set
  if (is.null(data$y)) stop("Missing y: Coalescent data must be in proper format -- use make_coal_data function")
  if (is.null(data$N)) stop("Missing N: Coalescent data must be in proper format -- use make_coal_data function")
  if (is.null(data$J)) stop("Missing J: Coalescent data must be in proper format -- use make_coal_data function")
  if (is.null(data$gridrep)) stop("Missing gridrep: Coalescent data must be in proper format -- use make_coal_data function")
  if (is.null(data$Aik)) stop("Missing Aik: Coalescent data must be in proper format -- use make_coal_data function")
  if (is.null(data$dalpha)) stop("Missing dalpha: Coalescent data must be in proper format -- use make_coal_data function")
  if (is.null(data$log_mu)) stop("Missing log_mu: Coalescent data must be in proper format -- use make_coal_data function")


  tmp.dat <<- data
  mcode <- get_model_ada(prior.c=prior.c, prior.s=prior.s, order.c=order.c,order.s=order.s, zeta.c=zeta.c,zeta.s=zeta.s, save.loglik=save.loglik)
  finits <- get_init_ada(prior.c=prior.c, prior.s=prior.s, order.c=order.c,order.s=order.s)

  mname <- paste(prior.c, prior.s, order.c, order.s , sep="_")

  if ("init" %in% names(stlist) ) {
    if (!is.na(stlist$init)){
      stop("Do not specify an 'init' file - one will be automatically generated.")
    }
  }

  if ("file" %in% names(stlist) ) {
    if (!is.na(stlist$file)){
      stop("External file specified for model code. Use spmrf arguments to create model or use previously compiled spmrf fit object.")
    }
  }
  if ("model_code" %in% names(stlist)  ) {
    if (!is.na(stlist$file) & stlist$model_code!=""){
      stop("Argument 'model_code' specified. Use spmrf arguments to create model or use previously compiled spmrf fit object.")
    }
  }

  if (class(fit)=="stanfit") {
    #warning("Warning: Make sure 'prior','likelihood', and 'order' match those of 'fit' object!")
    sfit <- rstan::stan( model_name=mname, fit=fit, data=data, pars=pars, chains=chains, iter=iter, warmup=warmup,
                  thin=thin, init=finits, control=control, ...)
  }

  if (class(fit)!="stanfit"){
    if (class(fit)!="logical") stop("Must specify fit as 'stanfit' object or as NA.")
    if (class(fit)=="logical") {
      if (!is.na(fit)) stop("Must specify fit as 'stanfit' object or as NA.")
      if (is.na(fit)) {
        sfit <- rstan::stan(model_name=mname, model_code = mcode, fit=fit, data=data, chains=chains, iter=iter, warmup=warmup, thin=thin, init=finits, control=control, ...)
      }
    }
  }

  rm(tmp.dat, envir=.GlobalEnv)

  return(sfit)
}



get_model_ada <- function(prior.c="normal", prior.s="horseshoe", order.c=1, order.s=1, zeta.c=0.01,zeta.s=0.01, save.loglik=FALSE){


  if (zeta.c <= 0) stop("zeta.c must be > 0.")
  if (zeta.s <= 0) stop("zeta.s must be > 0.")

  ###  MODEL ORDER AND PRIOR CHOICE  #########
  # The first part is common to all of them ##
  general<-'
  functions {
  real coal_loglik_lp(vector ft, vector yy, int nc, vector aik, vector da) {
  vector [nc] ll;
  real sll;
  ll = -1.0*yy .* ft - da .* aik .* exp(-ft);
  sll = sum(ll);
  return sll ;
  }
  real sampling_loglik(vector t, vector ga, vector ss, int ns, vector da){
  vector [ns] ll;
  real sll;
  ll=(ga+t+log(da)) .* ss - exp(t) .* exp(ga) .* da;
  sll = sum(ll);
  return sll ;
  }
  vector coal_and_samp_loglik_vector(vector ft, vector yy, int nc, vector aik, vector da, int nci, int[] cstr, int[] cnd ,vector t, vector ga, vector ss, int ns, vector dg) {
  vector [nc] coal;
  vector [nci+ns] log_lik_V;
  coal = -1.0*yy .* ft - da .* aik .* exp(-ft);
  for (j in 1:nci){
  log_lik_V[j] = sum(coal[cstr[j]:cnd[j]]);
  }
  log_lik_V[(nci+1):(nci+ns)]=(ga+t+log(dg)) .* ss - exp(t) .* exp(ga) .* dg;
  return log_lik_V;
  }
  }

  data {
  int <lower=1> J; //number of grid points (theta params)
  int <lower=1> N; //number of grid subsections
  int <lower=1> S; //number of grid points up to which there are sampling events
  vector <lower=0>[N] y; //auxillary coal indicator
  int <lower=1> coalind[N]; // id variable for coalescent number
  int <lower=1> ncoal;  // number of coalescent events
  int <lower=1>  cstart[ncoal]; //indices for coal interval starts
  int <lower=1> cend[ncoal] ; //indices for coal interval ends
  int <lower=1>  gridrep [J]; //number of reps per theta
  vector <lower=0> [N] Aik; // active lineage combinatoric coefficient
  vector [N] dalpha; //delta alpha - subgrid widths
  vector [S] dgrid; //width of grid intervals in the sampling grid
  real log_mu; //mle for const Ne on log scale
  real log_beta; //mle for const samp_intensity on log scale
  vector <lower=0>[S] scount; //number of sampled for each grid point (up to )
  }
  '

  ####################
  ## BOTH ORDER 1 ####
  ####################


  ## Coalescent GMRF Order 1 --- Sampling HSMRF Order 1
  N.c_1_H.s_1_temp <- '
  parameters {
  vector [J-1] zdelta;
  real ztheta1;
  real <lower=0, upper=1> zgam;

  vector [S-1] zpsi;
  real zalpha1;
  vector <lower=0, upper=1> [S-1] zrho;
  real <lower=0, upper=1> zxi;
  }

  transformed parameters {
  vector [J] theta;
  vector [N] ftheta;
  real <lower=0> gam;

  vector [S] alpha;
  real <lower=0> xi;
  vector[S-1] rho;

  gam = ZETAVAL1*tan(zgam*pi()/2);
  theta[1] = 10*ztheta1 + log_mu;
  for (j in 1:(J-1)){
  theta[j+1] = gam*zdelta[j] + theta[j];
  }
  { int cnt;
  cnt = 0;
  for (j in 1:J){
  for (k in 1:gridrep[j]){
  cnt = cnt + 1;
  ftheta[cnt] = theta[j];
  }
  }
  }

  xi = ZETAVAL2*tan(zxi*pi()/2);
  alpha[1] = 10*zalpha1+log_beta;
  for (j in 1:(S-1)){
  rho[j] = xi*tan(zrho[j]*pi()/2);
  alpha[j+1] = zpsi[j]*rho[j] + alpha[j];
  }
  }

  model {
  zgam ~ uniform(0, 1);
  ztheta1 ~ normal(0, 1);
  zdelta ~ normal(0, 1);

  zxi ~ uniform(0, 1);
  zrho ~ uniform(0, 1);
  zalpha1 ~ normal(log_mu, 10);
  zpsi ~ normal(0, 1);

  target += coal_loglik_lp(ftheta, y, N, Aik, dalpha);
  target += sampling_loglik(theta[1:S], alpha, scount, S, dgrid);
  }

  LOGLIKCALC
  '

  ## Coalescent GMRF Order 1 --- Sampling GMRF Order 1
  N.c_1_N.s_1_temp <- '
  parameters {
  vector [J-1] zdelta;
  real ztheta1;
  real <lower=0, upper=1> zgam;

  vector [S-1] zpsi;
  real zalpha1;
  real <lower=0, upper=1> zxi;
  }

  transformed parameters {
  vector [J] theta;
  vector [N] ftheta;
  real <lower=0> gam;

  vector [S] alpha;
  real <lower=0> xi;

  gam = ZETAVAL1*tan(zgam*pi()/2);
  theta[1] = 10*ztheta1 + log_mu;
  for (j in 1:(J-1)){
  theta[j+1] = gam*zdelta[j] + theta[j];
  }
  { int cnt;
  cnt = 0;
  for (j in 1:J){
  for (k in 1:gridrep[j]){
  cnt = cnt + 1;
  ftheta[cnt] = theta[j];
  }
  }
  }

  xi = ZETAVAL2*tan(zxi*pi()/2);
  alpha[1] = 10*zalpha1+log_beta;
  for (j in 1:(S-1)){
  alpha[j+1] = zpsi[j]*xi + alpha[j];
  }
  }

  model {
  zgam ~ uniform(0, 1);
  ztheta1 ~ normal(0, 1);
  zdelta ~ normal(0, 1);

  zxi ~ uniform(0, 1);
  zalpha1 ~ normal(log_mu, 10);
  zpsi ~ normal(0, 1);

  target += coal_loglik_lp(ftheta, y, N, Aik, dalpha);
  target += sampling_loglik(theta[1:S], alpha, scount, S, dgrid);
  }

  LOGLIKCALC
  '


  ## Coalescent HSMRF Order 1 --- Sampling HSMRF Order 1
  H.c_1_H.s_1_temp <- '
  parameters {
  vector [J-1] zdelta;
  real ztheta1;
  vector <lower=0, upper=1> [J-1] ztau;
  real <lower=0, upper=1> zgam;

  vector [S-1] zpsi;
  real zalpha1;
  vector <lower=0, upper=1> [S-1] zrho;
  real <lower=0, upper=1> zxi;
  }

  transformed parameters {
  vector [J] theta;
  vector [N] ftheta;
  real <lower=0> gam;
  vector[J-1] tau;

  vector [S] alpha;
  real <lower=0> xi;
  vector[S-1] rho;

  gam = ZETAVAL1*tan(zgam*pi()/2);
  theta[1] = 10*ztheta1 + log_mu;
  for (j in 1:(J-1)){
  tau[j] = gam*tan(ztau[j]*pi()/2);
  theta[j+1] = zdelta[j]*tau[j] + theta[j];
  }
  { int cnt;
  cnt = 0;
  for (j in 1:J){
  for (k in 1:gridrep[j]){
  cnt = cnt + 1;
  ftheta[cnt] = theta[j];
  }
  }
  }

  xi = ZETAVAL2*tan(zxi*pi()/2);
  alpha[1] = 10*zalpha1+log_beta;
  for (j in 1:(S-1)){
  rho[j] = xi*tan(zrho[j]*pi()/2);
  alpha[j+1] = zpsi[j]*rho[j] + alpha[j];
  }
  }

  model {
  zgam ~ uniform(0, 1);
  ztau ~ uniform(0, 1);
  ztheta1 ~ normal(0, 1);
  zdelta ~ normal(0, 1);

  zxi ~ uniform(0, 1);
  zrho ~ uniform(0, 1);
  zalpha1 ~ normal(log_mu, 10);
  zpsi ~ normal(0, 1);

  target += coal_loglik_lp(ftheta, y, N, Aik, dalpha);
  target += sampling_loglik(theta[1:S], alpha, scount, S, dgrid);
  }

  LOGLIKCALC
  '


  ## Coalescent HSMRF Order 1 --- Sampling GMRF Order 1
  H.c_1_N.s_1_temp <- '
  parameters {
  vector [J-1] zdelta;
  real ztheta1;
  vector <lower=0, upper=1> [J-1] ztau;
  real <lower=0, upper=1> zgam;

  vector [S-1] zpsi;
  real zalpha1;
  real <lower=0, upper=1> zxi;
  }

  transformed parameters {
  vector [J] theta;
  vector [N] ftheta;
  real <lower=0> gam;
  vector[J-1] tau;

  vector [S] alpha;
  real <lower=0> xi;

  gam = ZETAVAL1*tan(zgam*pi()/2);
  theta[1] = 10*ztheta1 + log_mu;
  for (j in 1:(J-1)){
  tau[j] = gam*tan(ztau[j]*pi()/2);
  theta[j+1] = zdelta[j]*tau[j] + theta[j];
  }
  { int cnt;
  cnt = 0;
  for (j in 1:J){
  for (k in 1:gridrep[j]){
  cnt = cnt + 1;
  ftheta[cnt] = theta[j];
  }
  }
  }

  xi = ZETAVAL2*tan(zxi*pi()/2);
  alpha[1] = 10*zalpha1+log_beta;
  for (j in 1:(S-1)){
  alpha[j+1] = zpsi[j]*xi + alpha[j];
  }
  }

  model {
  zgam ~ uniform(0, 1);
  ztau ~ uniform(0, 1);
  ztheta1 ~ normal(0, 1);
  zdelta ~ normal(0, 1);

  zxi ~ uniform(0, 1);
  zalpha1 ~ normal(log_mu, 10);
  zpsi ~ normal(0, 1);

  target += coal_loglik_lp(ftheta, y, N, Aik, dalpha);
  target += sampling_loglik(theta[1:S], alpha, scount, S, dgrid);
  }

  LOGLIKCALC
  '

  #################################
  ## COAL ORDER 2  SAMP ORDER 1####
  #################################




  ## Coalescent GMRF Order 2 --- Sampling HSMRF Order 1
  N.c_2_H.s_1_temp <- '
  parameters {
  vector [J-1] zdelta;
  real ztheta1;
  real <lower=0, upper=1> zgam;

  vector [S-1] zpsi;
  real zalpha1;
  vector <lower=0, upper=1> [S-1] zrho;
  real <lower=0, upper=1> zxi;
  }

  transformed parameters {
  vector [J] theta;
  vector [N] ftheta;
  real <lower=0> gam;

  vector [S] alpha;
  real <lower=0> xi;
  vector[S-1] rho;

  gam = ZETAVAL1*tan(zgam*pi()/2);
  theta[1] = 10*ztheta1 + log_mu;
  theta[2] = sqrt(0.5)*gam*zdelta[1] + theta[1];
  for (j in 1:(J-2)){
  theta[j+2] = gam*zdelta[j+1] + 2*theta[j+1]-theta[j];
  }
  { int cnt;
  cnt = 0;
  for (j in 1:J){
  for (k in 1:gridrep[j]){
  cnt = cnt + 1;
  ftheta[cnt] = theta[j];
  }
  }
  }

  xi = ZETAVAL2*tan(zxi*pi()/2);
  alpha[1] = 10*zalpha1+log_beta;
  for (j in 1:(S-1)){
  rho[j] = xi*tan(zrho[j]*pi()/2);
  alpha[j+1] = zpsi[j]*rho[j] + alpha[j];
  }
  }

  model {
  zgam ~ uniform(0, 1);
  ztheta1 ~ normal(0, 1);
  zdelta ~ normal(0, 1);

  zxi ~ uniform(0, 1);
  zrho ~ uniform(0, 1);
  zalpha1 ~ normal(log_mu, 10);
  zpsi ~ normal(0, 1);

  target += coal_loglik_lp(ftheta, y, N, Aik, dalpha);
  target += sampling_loglik(theta[1:S], alpha, scount, S, dgrid);
  }

  LOGLIKCALC
  '

  ## Coalescent GMRF Order 2 --- Sampling GMRF Order 1
  N.c_2_N.s_1_temp <- '
  parameters {
  vector [J-1] zdelta;
  real ztheta1;
  real <lower=0, upper=1> zgam;

  vector [S-1] zpsi;
  real zalpha1;
  real <lower=0, upper=1> zxi;
  }

  transformed parameters {
  vector [J] theta;
  vector [N] ftheta;
  real <lower=0> gam;

  vector [S] alpha;
  real <lower=0> xi;

  gam = ZETAVAL1*tan(zgam*pi()/2);
  theta[1] = 10*ztheta1 + log_mu;
  theta[2] = sqrt(0.5)*gam*zdelta[1] + theta[1];
  for (j in 1:(J-2)){
  theta[j+2] = gam*zdelta[j+1] + 2*theta[j+1]-theta[j];
  }
  { int cnt;
  cnt = 0;
  for (j in 1:J){
  for (k in 1:gridrep[j]){
  cnt = cnt + 1;
  ftheta[cnt] = theta[j];
  }
  }
  }

  xi = ZETAVAL2*tan(zxi*pi()/2);
  alpha[1] = 10*zalpha1+log_beta;
  for (j in 1:(S-1)){
  alpha[j+1] = zpsi[j]*xi + alpha[j];
  }
  }

  model {
  zgam ~ uniform(0, 1);
  ztheta1 ~ normal(0, 1);
  zdelta ~ normal(0, 1);

  zxi ~ uniform(0, 1);
  zalpha1 ~ normal(log_mu, 10);
  zpsi ~ normal(0, 1);

  target += coal_loglik_lp(ftheta, y, N, Aik, dalpha);
  target += sampling_loglik(theta[1:S], alpha, scount, S, dgrid);
  }

  LOGLIKCALC
  '


  ## Coalescent HSMRF Order 2 --- Sampling HSMRF Order 1
  H.c_2_H.s_1_temp <- '
  parameters {
  vector [J-1] zdelta;
  real ztheta1;
  vector <lower=0, upper=1> [J-2] ztau;
  real <lower=0, upper=1> zgam;
  real <lower=0, upper=1> zptau2;

  vector [S-1] zpsi;
  real zalpha1;
  vector <lower=0, upper=1> [S-1] zrho;
  real <lower=0, upper=1> zxi;
  }

  transformed parameters {
  vector [J] theta;
  vector [N] ftheta;
  real <lower=0> gam;
  real <lower=0> ptau2;
  vector[J-2] tau;

  vector [S] alpha;
  real <lower=0> xi;
  vector[S-1] rho;

  gam = ZETAVAL1*tan(zgam*pi()/2);
  ptau2 = (gam/sqrt(2.0))*tan(zptau2*pi()/2);
  theta[1] = 10*ztheta1 + log_mu ;
  theta[2] = ptau2*zdelta[1] + theta[1];
  for (j in 1:(J-2)){
  tau[j] = gam*tan(ztau[j]*pi()/2);
  theta[j+2] = zdelta[j+1]*tau[j] + 2*theta[j+1]-theta[j];
  }
  { int cnt;
  cnt = 0;
  for (j in 1:J){
  for (k in 1:gridrep[j]){
  cnt = cnt + 1;
  ftheta[cnt] = theta[j];
  }
  }
  }

  xi = ZETAVAL2*tan(zxi*pi()/2);
  alpha[1] = 10*zalpha1+log_beta;
  for (j in 1:(S-1)){
  rho[j] = xi*tan(zrho[j]*pi()/2);
  alpha[j+1] = zpsi[j]*rho[j] + alpha[j];
  }
  }

  model {
  zgam ~ uniform(0, 1);
  ztau ~ uniform(0, 1);
  ztheta1 ~ normal(0, 1);
  zdelta ~ normal(0, 1);
  zptau2 ~ uniform(0, 1);

  zxi ~ uniform(0, 1);
  zrho ~ uniform(0, 1);
  zalpha1 ~ normal(log_mu, 10);
  zpsi ~ normal(0, 1);

  target += coal_loglik_lp(ftheta, y, N, Aik, dalpha);
  target += sampling_loglik(theta[1:S], alpha, scount, S, dgrid);
  }

  LOGLIKCALC
  '


  ## Coalescent HSMRF Order 2 --- Sampling GMRF Order 1
  H.c_2_N.s_1_temp <- '
  parameters {
  vector [J-1] zdelta;
  real ztheta1;
  vector <lower=0, upper=1> [J-2] ztau;
  real <lower=0, upper=1> zgam;
  real <lower=0, upper=1> zptau2;

  vector [S-1] zpsi;
  real zalpha1;
  real <lower=0, upper=1> zxi;
  }

  transformed parameters {
  vector [J] theta;
  vector [N] ftheta;
  real <lower=0> gam;
  real <lower=0> ptau2;
  vector[J-2] tau;

  vector [S] alpha;
  real <lower=0> xi;

  gam = ZETAVAL1*tan(zgam*pi()/2);
  ptau2 = (gam/sqrt(2.0))*tan(zptau2*pi()/2);
  theta[1] = 10*ztheta1 + log_mu ;
  theta[2] = ptau2*zdelta[1] + theta[1];
  for (j in 1:(J-2)){
  tau[j] = gam*tan(ztau[j]*pi()/2);
  theta[j+2] = zdelta[j+1]*tau[j] + 2*theta[j+1]-theta[j];
  }
  { int cnt;
  cnt = 0;
  for (j in 1:J){
  for (k in 1:gridrep[j]){
  cnt = cnt + 1;
  ftheta[cnt] = theta[j];
  }
  }
  }

  xi = ZETAVAL2*tan(zxi*pi()/2);
  alpha[1] = 10*zalpha1+log_beta;
  for (j in 1:(S-1)){
  alpha[j+1] = zpsi[j]*xi + alpha[j];
  }
  }

  model {
  zgam ~ uniform(0, 1);
  ztau ~ uniform(0, 1);
  zptau2 ~ uniform(0, 1);
  ztheta1 ~ normal(0, 1);
  zdelta ~ normal(0, 1);

  zxi ~ uniform(0, 1);
  zalpha1 ~ normal(log_mu, 10);
  zpsi ~ normal(0, 1);

  target += coal_loglik_lp(ftheta, y, N, Aik, dalpha);
  target += sampling_loglik(theta[1:S], alpha, scount, S, dgrid);
  }

  LOGLIKCALC
  '

  #################################
  ## COAL ORDER 1  SAMP ORDER 2####
  #################################


  ## Coalescent GMRF Order 1 --- Sampling HSMRF Order 2
  N.c_1_H.s_2_temp <- '
  parameters {
  vector [J-1] zdelta;
  real ztheta1;
  real <lower=0, upper=1> zgam;

  vector [S-1] zpsi;
  real zalpha1;
  vector <lower=0, upper=1> [S-2] zrho;
  real <lower=0, upper=1> zxi;
  real <lower=0, upper=1> zprho2;
  }

  transformed parameters {
  vector [J] theta;
  vector [N] ftheta;
  real <lower=0> gam;

  vector [S] alpha;
  real <lower=0> xi;
  real <lower=0> prho2;
  vector[S-2] rho;

  gam = ZETAVAL1*tan(zgam*pi()/2);
  theta[1] = 10*ztheta1 + log_mu;
  for (j in 1:(J-1)){
  theta[j+1] = gam*zdelta[j] + theta[j];
  }
  { int cnt;
  cnt = 0;
  for (j in 1:J){
  for (k in 1:gridrep[j]){
  cnt = cnt + 1;
  ftheta[cnt] = theta[j];
  }
  }
  }

  xi = ZETAVAL2*tan(zxi*pi()/2);
  prho2 = (xi/sqrt(2.0))*tan(zprho2*pi()/2);
  alpha[1] = 10*zalpha1+log_beta;
  alpha[2] = prho2*zpsi[1] + alpha[1];
  for (j in 1:(S-2)){
  rho[j] = xi*tan(zrho[j]*pi()/2);
  alpha[j+2] = zpsi[j+1]*rho[j] + 2*alpha[j+1]-alpha[j];
  }
  }

  model {
  zgam ~ uniform(0, 1);
  ztheta1 ~ normal(0, 1);
  zdelta ~ normal(0, 1);

  zxi ~ uniform(0, 1);
  zrho ~ uniform(0, 1);
  zprho2 ~ uniform(0, 1);
  zalpha1 ~ normal(log_mu, 10);
  zpsi ~ normal(0, 1);


  target += coal_loglik_lp(ftheta, y, N, Aik, dalpha);
  target += sampling_loglik(theta[1:S], alpha, scount, S, dgrid);
  }

  LOGLIKCALC
  '

  ## Coalescent GMRF Order 1 --- Sampling GMRF Order 2
  N.c_1_N.s_2_temp <- '
  parameters {
  vector [J-1] zdelta;
  real ztheta1;
  real <lower=0, upper=1> zgam;

  vector [S-1] zpsi;
  real zalpha1;
  real <lower=0, upper=1> zxi;
  }

  transformed parameters {
  vector [J] theta;
  vector [N] ftheta;
  real <lower=0> gam;

  vector [S] alpha;
  real <lower=0> xi;

  gam = ZETAVAL1*tan(zgam*pi()/2);
  theta[1] = 10*ztheta1 + log_mu;
  for (j in 1:(J-1)){
  theta[j+1] = gam*zdelta[j] + theta[j];
  }
  { int cnt;
  cnt = 0;
  for (j in 1:J){
  for (k in 1:gridrep[j]){
  cnt = cnt + 1;
  ftheta[cnt] = theta[j];
  }
  }
  }

  xi = ZETAVAL2*tan(zxi*pi()/2);
  alpha[1] = 10*zalpha1+log_beta;
  alpha[2] = sqrt(0.5)*gam*zpsi[1] + alpha[1];
  for (j in 1:(S-2)){
  alpha[j+2] = xi*zpsi[j+1] + 2*alpha[j+1]-alpha[j];
  }
  }

  model {
  zgam ~ uniform(0, 1);
  ztheta1 ~ normal(0, 1);
  zdelta ~ normal(0, 1);

  zxi ~ uniform(0, 1);
  zalpha1 ~ normal(log_mu, 10);
  zpsi ~ normal(0, 1);


  target += coal_loglik_lp(ftheta, y, N, Aik, dalpha);
  target += sampling_loglik(theta[1:S], alpha, scount, S, dgrid);
  }

  LOGLIKCALC
  '


  ## Coalescent HSMRF Order 1 --- Sampling HSMRF Order 2
  H.c_1_H.s_2_temp <- '
  parameters {
  vector [J-1] zdelta;
  real ztheta1;
  vector <lower=0, upper=1> [J-1] ztau;
  real <lower=0, upper=1> zgam;

  vector [S-1] zpsi;
  real zalpha1;
  vector <lower=0, upper=1> [S-2] zrho;
  real <lower=0, upper=1> zxi;
  real <lower=0, upper=1> zprho2;
  }

  transformed parameters {
  vector [J] theta;
  vector [N] ftheta;
  real <lower=0> gam;
  vector[J-1] tau;

  vector [S] alpha;
  real <lower=0> xi;
  real <lower=0> prho2;
  vector[S-2] rho;

  gam = ZETAVAL1*tan(zgam*pi()/2);
  theta[1] = 10*ztheta1 + log_mu;
  for (j in 1:(J-1)){
  tau[j] = gam*tan(ztau[j]*pi()/2);
  theta[j+1] = zdelta[j]*tau[j] + theta[j];
  }
  { int cnt;
  cnt = 0;
  for (j in 1:J){
  for (k in 1:gridrep[j]){
  cnt = cnt + 1;
  ftheta[cnt] = theta[j];
  }
  }
  }

  xi = ZETAVAL2*tan(zxi*pi()/2);
  prho2 = (xi/sqrt(2.0))*tan(zprho2*pi()/2);
  alpha[1] = 10*zalpha1+log_beta;
  alpha[2] = prho2*zpsi[1] + alpha[1];
  for (j in 1:(S-2)){
  rho[j] = xi*tan(zrho[j]*pi()/2);
  alpha[j+2] = zpsi[j+1]*rho[j] + 2*alpha[j+1]-alpha[j];
  }
  }

  model {
  zgam ~ uniform(0, 1);
  ztau ~ uniform(0, 1);
  ztheta1 ~ normal(0, 1);
  zdelta ~ normal(0, 1);

  zxi ~ uniform(0, 1);
  zrho ~ uniform(0, 1);
  zprho2 ~ uniform(0, 1);
  zalpha1 ~ normal(log_mu, 10);
  zpsi ~ normal(0, 1);

  target += coal_loglik_lp(ftheta, y, N, Aik, dalpha);
  target += sampling_loglik(theta[1:S], alpha, scount, S, dgrid);
  }

  LOGLIKCALC
  '


  ## Coalescent HSMRF Order 1 --- Sampling GMRF Order 2
  H.c_1_N.s_2_temp <- '
  parameters {
  vector [J-1] zdelta;
  real ztheta1;
  vector <lower=0, upper=1> [J-1] ztau;
  real <lower=0, upper=1> zgam;

  vector [S-1] zpsi;
  real zalpha1;
  real <lower=0, upper=1> zxi;
  }

  transformed parameters {
  vector [J] theta;
  vector [N] ftheta;
  real <lower=0> gam;
  vector[J-1] tau;

  vector [S] alpha;
  real <lower=0> xi;

  gam = ZETAVAL1*tan(zgam*pi()/2);
  theta[1] = 10*ztheta1 + log_mu;
  for (j in 1:(J-1)){
  tau[j] = gam*tan(ztau[j]*pi()/2);
  theta[j+1] = zdelta[j]*tau[j] + theta[j];
  }
  { int cnt;
  cnt = 0;
  for (j in 1:J){
  for (k in 1:gridrep[j]){
  cnt = cnt + 1;
  ftheta[cnt] = theta[j];
  }
  }
  }

  xi = ZETAVAL2*tan(zxi*pi()/2);
  alpha[1] = 10*zalpha1+log_beta;
  alpha[2] = sqrt(0.5)*gam*zpsi[1] + alpha[1];
  for (j in 1:(S-2)){
  alpha[j+2] = xi*zpsi[j+1] + 2*alpha[j+1]-alpha[j];
  }
  }

  model {
  zgam ~ uniform(0, 1);
  ztau ~ uniform(0, 1);
  ztheta1 ~ normal(0, 1);
  zdelta ~ normal(0, 1);

  zxi ~ uniform(0, 1);
  zalpha1 ~ normal(log_mu, 10);
  zpsi ~ normal(0, 1);

  target += coal_loglik_lp(ftheta, y, N, Aik, dalpha);
  target += sampling_loglik(theta[1:S], alpha, scount, S, dgrid);
  }

  LOGLIKCALC
  '

  #################################
  ## COAL ORDER 2  SAMP ORDER 2####
  #################################




  ## Coalescent GMRF Order 2 --- Sampling HSMRF Order 2
  N.c_2_H.s_2_temp <- '
  parameters {
  vector [J-1] zdelta;
  real ztheta1;
  real <lower=0, upper=1> zgam;

  vector [S-1] zpsi;
  real zalpha1;
  vector <lower=0, upper=1> [S-2] zrho;
  real <lower=0, upper=1> zxi;
  real <lower=0, upper=1> zprho2;
  }

  transformed parameters {
  vector [J] theta;
  vector [N] ftheta;
  real <lower=0> gam;

  vector [S] alpha;
  real <lower=0> xi;
  real <lower=0> prho2;
  vector[S-1] rho;

  gam = ZETAVAL1*tan(zgam*pi()/2);
  theta[1] = 10*ztheta1 + log_mu;
  theta[2] = sqrt(0.5)*gam*zdelta[1] + theta[1];
  for (j in 1:(J-2)){
  theta[j+2] = gam*zdelta[j+1] + 2*theta[j+1]-theta[j];
  }
  { int cnt;
  cnt = 0;
  for (j in 1:J){
  for (k in 1:gridrep[j]){
  cnt = cnt + 1;
  ftheta[cnt] = theta[j];
  }
  }
  }

  xi = ZETAVAL2*tan(zxi*pi()/2);
  prho2 = (xi/sqrt(2.0))*tan(zprho2*pi()/2);
  alpha[1] = 10*zalpha1+log_beta;
  alpha[2] = prho2*zpsi[1] + alpha[1];
  for (j in 1:(S-2)){
  rho[j] = xi*tan(zrho[j]*pi()/2);
  alpha[j+2] = zpsi[j+1]*rho[j] + 2*alpha[j+1]-alpha[j];
  }
  }

  model {
  zgam ~ uniform(0, 1);
  ztheta1 ~ normal(0, 1);
  zdelta ~ normal(0, 1);

  zxi ~ uniform(0, 1);
  zrho ~ uniform(0, 1);
  zprho2 ~ uniform(0, 1);
  zalpha1 ~ normal(log_mu, 10);
  zpsi ~ normal(0, 1);

  target += coal_loglik_lp(ftheta, y, N, Aik, dalpha);
  target += sampling_loglik(theta[1:S], alpha, scount, S, dgrid);
  }

  LOGLIKCALC
  '

  ## Coalescent GMRF Order 2 --- Sampling GMRF Order 2
  N.c_2_N.s_2_temp <- '
  parameters {
  vector [J-1] zdelta;
  real ztheta1;
  real <lower=0, upper=1> zgam;

  vector [S-1] zpsi;
  real zalpha1;
  real <lower=0, upper=1> zxi;
  }

  transformed parameters {
  vector [J] theta;
  vector [N] ftheta;
  real <lower=0> gam;

  vector [S] alpha;
  real <lower=0> xi;

  gam = ZETAVAL1*tan(zgam*pi()/2);
  theta[1] = 10*ztheta1 + log_mu;
  theta[2] = sqrt(0.5)*gam*zdelta[1] + theta[1];
  for (j in 1:(J-2)){
  theta[j+2] = gam*zdelta[j+1] + 2*theta[j+1]-theta[j];
  }
  { int cnt;
  cnt = 0;
  for (j in 1:J){
  for (k in 1:gridrep[j]){
  cnt = cnt + 1;
  ftheta[cnt] = theta[j];
  }
  }
  }

  xi = ZETAVAL2*tan(zxi*pi()/2);
  alpha[1] = 10*zalpha1+log_beta;
  alpha[2] = sqrt(0.5)*gam*zpsi[1] + alpha[1];
  for (j in 1:(S-2)){
  alpha[j+2] = xi*zpsi[j+1] + 2*alpha[j+1]-alpha[j];
  }
  }

  model {
  zgam ~ uniform(0, 1);
  ztheta1 ~ normal(0, 1);
  zdelta ~ normal(0, 1);

  zxi ~ uniform(0, 1);
  zalpha1 ~ normal(log_mu, 10);
  zpsi ~ normal(0, 1);

  target += coal_loglik_lp(ftheta, y, N, Aik, dalpha);
  target += sampling_loglik(theta[1:S], alpha, scount, S, dgrid);
  }

  LOGLIKCALC
  '


  ## Coalescent HSMRF Order 2 --- Sampling HSMRF Order 2
  H.c_2_H.s_2_temp <- '
  parameters {
  vector [J-1] zdelta;
  real ztheta1;
  vector <lower=0, upper=1> [J-2] ztau;
  real <lower=0, upper=1> zgam;
  real <lower=0, upper=1> zptau2;

  vector [S-1] zpsi;
  real zalpha1;
  vector <lower=0, upper=1> [S-2] zrho;
  real <lower=0, upper=1> zxi;
  real <lower=0, upper=1> zprho2;
  }

  transformed parameters {
  vector [J] theta;
  vector [N] ftheta;
  real <lower=0> gam;
  real <lower=0> ptau2;
  vector[J-2] tau;

  vector [S] alpha;
  real <lower=0> xi;
  real <lower=0> prho2;
  vector[S-1] rho;

  gam = ZETAVAL1*tan(zgam*pi()/2);
  ptau2 = (gam/sqrt(2.0))*tan(zptau2*pi()/2);
  theta[1] = 10*ztheta1 + log_mu ;
  theta[2] = ptau2*zdelta[1] + theta[1];
  for (j in 1:(J-2)){
  tau[j] = gam*tan(ztau[j]*pi()/2);
  theta[j+2] = zdelta[j+1]*tau[j] + 2*theta[j+1]-theta[j];
  }
  { int cnt;
  cnt = 0;
  for (j in 1:J){
  for (k in 1:gridrep[j]){
  cnt = cnt + 1;
  ftheta[cnt] = theta[j];
  }
  }
  }

  xi = ZETAVAL2*tan(zxi*pi()/2);
  prho2 = (xi/sqrt(2.0))*tan(zprho2*pi()/2);
  alpha[1] = 10*zalpha1+log_beta;
  alpha[2] = prho2*zpsi[1] + alpha[1];
  for (j in 1:(S-2)){
  rho[j] = xi*tan(zrho[j]*pi()/2);
  alpha[j+2] = zpsi[j+1]*rho[j] + 2*alpha[j+1]-alpha[j];
  }
  }

  model {
  zgam ~ uniform(0, 1);
  ztau ~ uniform(0, 1);
  ztheta1 ~ normal(0, 1);
  zdelta ~ normal(0, 1);
  zptau2 ~ uniform(0, 1);

  zxi ~ uniform(0, 1);
  zrho ~ uniform(0, 1);
  zprho2 ~ uniform(0, 1);
  zalpha1 ~ normal(log_mu, 10);
  zpsi ~ normal(0, 1);

  target += coal_loglik_lp(ftheta, y, N, Aik, dalpha);
  target += sampling_loglik(theta[1:S], alpha, scount, S, dgrid);
  }

  LOGLIKCALC
  '


  ## Coalescent HSMRF Order 2 --- Sampling GMRF Order 2
  H.c_2_N.s_2_temp <- '
  parameters {
  vector [J-1] zdelta;
  real ztheta1;
  vector <lower=0, upper=1> [J-2] ztau;
  real <lower=0, upper=1> zgam;
  real <lower=0, upper=1> zptau2;

  vector [S-1] zpsi;
  real zalpha1;
  real <lower=0, upper=1> zxi;
  }

  transformed parameters {
  vector [J] theta;
  vector [N] ftheta;
  real <lower=0> gam;
  real <lower=0> ptau2;
  vector[J-2] tau;

  vector [S] alpha;
  real <lower=0> xi;

  gam = ZETAVAL1*tan(zgam*pi()/2);
  ptau2 = (gam/sqrt(2.0))*tan(zptau2*pi()/2);
  theta[1] = 10*ztheta1 + log_mu ;
  theta[2] = ptau2*zdelta[1] + theta[1];
  for (j in 1:(J-2)){
  tau[j] = gam*tan(ztau[j]*pi()/2);
  theta[j+2] = zdelta[j+1]*tau[j] + 2*theta[j+1]-theta[j];
  }
  { int cnt;
  cnt = 0;
  for (j in 1:J){
  for (k in 1:gridrep[j]){
  cnt = cnt + 1;
  ftheta[cnt] = theta[j];
  }
  }
  }

  xi = ZETAVAL2*tan(zxi*pi()/2);
  alpha[1] = 10*zalpha1+log_beta;
  alpha[2] = sqrt(0.5)*gam*zpsi[1] + alpha[1];
  for (j in 1:(S-2)){
  alpha[j+2] = xi*zpsi[j+1] + 2*alpha[j+1]-alpha[j];
  }
  }

  model {
  zgam ~ uniform(0, 1);
  ztau ~ uniform(0, 1);
  zptau2 ~ uniform(0, 1);
  ztheta1 ~ normal(0, 1);
  zdelta ~ normal(0, 1);

  zxi ~ uniform(0, 1);
  zalpha1 ~ normal(log_mu, 10);
  zpsi ~ normal(0, 1);

  target += coal_loglik_lp(ftheta, y, N, Aik, dalpha);
  target += sampling_loglik(theta[1:S], alpha, scount, S, dgrid);
  }

  LOGLIKCALC
  '



  ## ###################
  #### Model choice ###
  ####################

  if (prior.c=="normal" & prior.s=="horseshoe" & order.c==1  & order.s==1) tmp.b <- paste(general,N.c_1_H.s_1_temp,sep=' ')
  if (prior.c=="normal" & prior.s=="normal" & order.c==1   & order.s==1) tmp.b <- paste(general,N.c_1_N.s_1_temp,sep=' ')
  if (prior.c=="horseshoe" & prior.s=="horseshoe" & order.c==1   & order.s==1) tmp.b <- paste(general,H.c_1_H.s_1_temp,sep=' ')
  if (prior.c=="horseshoe" & prior.s=="normal" & order.c==1   & order.s==1) tmp.b <- paste(general,H.c_1_N.s_1_temp,sep=' ')

  if (prior.c=="normal" & prior.s=="horseshoe" & order.c==2  & order.s==1) tmp.b <- paste(general,N.c_2_H.s_1_temp,sep=' ')
  if (prior.c=="normal" & prior.s=="normal" & order.c==2   & order.s==1) tmp.b <- paste(general,N.c_2_N.s_1_temp,sep=' ')
  if (prior.c=="horseshoe" & prior.s=="horseshoe" & order.c==2   & order.s==1) tmp.b <- paste(general,H.c_2_H.s_1_temp,sep=' ')
  if (prior.c=="horseshoe" & prior.s=="normal" & order.c==2   & order.s==1) tmp.b <- paste(general,H.c_2_N.s_1_temp,sep=' ')

  if (prior.c=="normal" & prior.s=="horseshoe" & order.c==1  & order.s==2) tmp.b <- paste(general,N.c_1_H.s_2_temp,sep=' ')
  if (prior.c=="normal" & prior.s=="normal" & order.c==1   & order.s==2) tmp.b <- paste(general,N.c_1_N.s_2_temp,sep=' ')
  if (prior.c=="horseshoe" & prior.s=="horseshoe" & order.c==1   & order.s==2) tmp.b <- paste(general,H.c_1_H.s_2_temp,sep=' ')
  if (prior.c=="horseshoe" & prior.s=="normal" & order.c==1   & order.s==2) tmp.b <- paste(general,H.c_1_N.s_2_temp,sep=' ')

  if (prior.c=="normal" & prior.s=="horseshoe" & order.c==2  & order.s==2) tmp.b <- paste(general,N.c_2_H.s_2_temp,sep=' ')
  if (prior.c=="normal" & prior.s=="normal" & order.c==2   & order.s==2) tmp.b <- paste(general,N.c_2_N.s_2_temp,sep=' ')
  if (prior.c=="horseshoe" & prior.s=="horseshoe" & order.c==2   & order.s==2) tmp.b <- paste(general,H.c_2_H.s_2_temp,sep=' ')
  if (prior.c=="horseshoe" & prior.s=="normal" & order.c==2   & order.s==2) tmp.b <- paste(general,H.c_2_N.s_2_temp,sep=' ')



  if (save.loglik==TRUE) {#for the time being this is going to be false always.
    tmp.txt <- ' generated quantities {
    vector[ncoal+S] log_lik;
    log_lik = coal_and_samp_loglik_vector(ftheta, y, N, Aik, dalpha, ncoal, cstart, cend,theta[1:S], alpha, scount, S, dgrid);
  }
    '
    tmp.b <- sub(pattern="LOGLIKCALC", replacement=tmp.txt, x=tmp.b)
}
  ## replace gamma #TO CHANGE TO ACCOMODATE TWO ZETAS
  tmp.c <- sub(pattern="ZETAVAL1", replacement=zeta.c, x=tmp.b)
  tmp.c <- sub(pattern="ZETAVAL2", replacement=zeta.s, x=tmp.c)




  if (save.loglik==FALSE) {
    tmp.txt <- '
    '
    tmp.c <- sub(pattern="LOGLIKCALC", replacement=tmp.txt, x=tmp.c)
  }

  return(tmp.c)
  }


get_init_ada <- function(prior.c="normal", prior.s="horseshoe", order.c=1,order.s=1) {
  # prior: ("horseshoe", "laplace", "normal")
  # likelihood: ("normal", "poisson", "binomial", "coalescent")
  # order: (1, 2, 3)



  #### ORder 1 both models

  if (prior.c=="normal" & prior.s=="horseshoe" & order.c==1 & order.s==1) {
    tmp.b <- 'function(chain_id=1, dat=tmp.dat) {
    zdel <- rnorm(dat$J-1, 0, 2)
    zth1 <- rnorm(1, 0, sd=0.5)
    zgm <- runif(1, 0.25, 0.75)

    zpsi <- rnorm(dat$S-1, 0, 2)
    zalph1 <- rnorm(1, 0, sd=0.5)
    zx <- runif(1, 0.25, 0.75)
    zrho <- runif(dat$S-1, 0.25, .75)
    list(zdelta =  zdel,
    ztheta1 = zth1,
    zgam = zgm,

    zpsi =  zpsi,
    zalpha1 = zalph1,
    zxi = zx,
    zrho = zrho
    )
  }'
  }


  if (prior.c=="normal" & prior.s=="normal" & order.c==1  & order.s==1) {
    tmp.b <- 'function(chain_id=1, dat=tmp.dat) {
    zdel <- rnorm(dat$J-1, 0, 2)
    zth1 <- rnorm(1, 0, sd=0.5)
    zgm <- runif(1, 0.25, 0.75)

    zpsi <- rnorm(dat$S-1, 0, 2)
    zalph1 <- rnorm(1, 0, sd=0.5)
    zx <- runif(1, 0.25, 0.75)

    list(zdelta =  zdel,
    ztheta1 = zth1,
    zgam = zgm,

    zpsi =  zpsi,
    zalpha1 = zalph1,
    zxi = zx
    )
  }'
  }


  if (prior.c=="horseshoe" & prior.s=="horseshoe" & order.c==1  & order.s==1) {
    tmp.b <- 'function(chain_id=1, dat=tmp.dat) {
    zdel <- rnorm(dat$J-1, 0, 2)
    zth1 <- rnorm(1, 0, sd=0.5)
    zgm <- runif(1, 0.2, 0.75)
    ztau <- runif(dat$J-1, 0.25, .75)

    zpsi <- rnorm(dat$S-1, 0, 2)
    zalph1 <- rnorm(1, 0, sd=0.5)
    zx <- runif(1, 0.25, 0.75)
    zrho <- runif(dat$S-1, 0.25, .75)
    list(zdelta =  zdel,
    ztheta1 = zth1,
    zgam = zgm,
    ztau = ztau,

    zpsi =  zpsi,
    zalpha1 = zalph1,
    zxi = zx,
    zrho = zrho
    )
  }'
  }


  if (prior.c=="horseshoe" & prior.s=="normal" & order.c==1 & order.s==1) {
    tmp.b <- 'function(chain_id=1, dat=tmp.dat) {
    zdel <- rnorm(dat$J-1, 0, 2)
    zth1 <- rnorm(1, 0, sd=0.5)
    zgm <- runif(1, 0.25, 0.75)
    ztau <- runif(dat$J-1, 0.25, .75)

    zpsi <- rnorm(dat$S-1, 0, 2)
    zalph1 <- rnorm(1, 0, sd=0.5)
    zx <- runif(1, 0.25, 0.75)
    list(zdelta =  zdel,
    ztheta1 = zth1,
    zgam = zgm,
    ztau = ztau,

    zpsi =  zpsi,
    zalpha1 = zalph1,
    zxi = zx
    )
  }'
  }


  #### Coalescent Order 2 Sampling Order 1 both models

  if (prior.c=="normal" & prior.s=="horseshoe" & order.c==2 & order.s==1) {
    tmp.b <- 'function(chain_id=1, dat=tmp.dat) {
    zdel <- rnorm(dat$J-1, 0, 4)
    zth1 <- rnorm(1, 0, sd=0.5)
    zgm <- runif(1, 0.4, 0.6)

    zpsi <- rnorm(dat$S-1, 0, 2)
    zalph1 <- rnorm(1, 0, sd=0.5)
    zx <- runif(1, 0.25, 0.75)
    zrho <- runif(dat$S-1, 0.25, .75)
    list(zdelta =  zdel,
    ztheta1 = zth1,
    zgam = zgm,

    zpsi =  zpsi,
    zalpha1 = zalph1,
    zxi = zx,
    zrho = zrho
    )
  }'
  }


  if (prior.c=="normal" & prior.s=="normal" & order.c==2 & order.s==1) {
    tmp.b <- 'function(chain_id=1, dat=tmp.dat) {
    zdel <- rnorm(dat$J-1, 0, 4)
    zth1 <- rnorm(1, 0, sd=0.5)
    zgm <- runif(1, 0.4, 0.6)

    zpsi <- rnorm(dat$S-1, 0, 2)
    zalph1 <- rnorm(1, 0, sd=0.5)
    zx <- runif(1, 0.25, 0.75)

    list(zdelta =  zdel,
    ztheta1 = zth1,
    zgam = zgm,

    zpsi =  zpsi,
    zalpha1 = zalph1,
    zxi = zx
    )
  }'
  }


  if (prior.c=="horseshoe" & prior.s=="horseshoe" & order.c==2 & order.s==1) {
    tmp.b <- 'function(chain_id=1, dat=tmp.dat) {
    zdel <- rnorm(dat$J-1, 0, 4)
    zth1 <- rnorm(1, 0, sd=0.5)
    zpt2 <- runif(1, 0.25, 0.75)
    zgm <- runif(1, 0.4, 0.6)
    ztau <- runif(dat$J-2, 0.4, .6)

    zpsi <- rnorm(dat$S-1, 0, 2)
    zalph1 <- rnorm(1, 0, sd=0.5)
    zx <- runif(1, 0.25, 0.75)
    zrho <- runif(dat$S-1, 0.25, .75)
    list(zdelta =  zdel,
    ztheta1 = zth1,
    zgam = zgm,
    ztau = ztau,
    zptau2 = zpt2,

    zpsi =  zpsi,
    zalpha1 = zalph1,
    zxi = zx,
    zrho = zrho
    )
  }'
  }


  if (prior.c=="horseshoe" & prior.s=="normal" & order.c==2 & order.s==1) {
    tmp.b <- 'function(chain_id=1, dat=tmp.dat) {
    zdel <- rnorm(dat$J-1, 0, 4)
    zth1 <- rnorm(1, 0, sd=0.5)
    zpt2 <- runif(1, 0.25, 0.75)
    zgm <- runif(1, 0.4, 0.6)
    ztau <- runif(dat$J-2, 0.4, .6)

    zpsi <- rnorm(dat$S-1, 0, 2)
    zalph1 <- rnorm(1, 0, sd=0.5)
    zx <- runif(1, 0.25, 0.75)

    list(zdelta =  zdel,
    ztheta1 = zth1,
    zgam = zgm,
    ztau = ztau,
    zptau2 = zpt2,

    zpsi =  zpsi,
    zalpha1 = zalph1,
    zxi = zx
    )
  }'
  }


  #### Coalescent Order 1 Sampling Order 2 both models

  if (prior.c=="normal" & prior.s=="horseshoe" & order.c==1 & order.s==2) {
    tmp.b <- 'function(chain_id=1, dat=tmp.dat) {
    zdel <- rnorm(dat$J-1, 0, 2)
    zth1 <- rnorm(1, 0, sd=0.5)
    zgm <- runif(1, 0.25, 0.75)

    zpsi <- rnorm(dat$S-1, 0, 4)
    zalph1 <- rnorm(1, 0, sd=0.5)
    zx <- runif(1, 0.4, 0.6)
    zrho <- runif(dat$S-2, 0.4, .6)
    zprho2 <- runif(1, 0.25, 0.75)

    list(zdelta =  zdel,
    ztheta1 = zth1,
    zgam = zgm,

    zpsi =  zpsi,
    zalpha1 = zalph1,
    zxi = zx,
    zprho2=zprho2,
    zrho = zrho
    )
  }'
  }


  if (prior.c=="normal" & prior.s=="normal" & order.c==1  & order.s==2) {
    tmp.b <- 'function(chain_id=1, dat=tmp.dat) {
    zdel <- rnorm(dat$J-1, 0, 2)
    zth1 <- rnorm(1, 0, sd=0.5)
    zgm <- runif(1, 0.25, 0.75)

    zpsi <- rnorm(dat$S-1, 0, 4)
    zalph1 <- rnorm(1, 0, sd=0.5)
    zx <- runif(1, 0.4, 0.6)

    list(zdelta =  zdel,
    ztheta1 = zth1,
    zgam = zgm,

    zpsi =  zpsi,
    zalpha1 = zalph1,
    zxi = zx
    )
  }'
  }


  if (prior.c=="horseshoe" & prior.s=="horseshoe" & order.c==1  & order.s==2) {
    tmp.b <- 'function(chain_id=1, dat=tmp.dat) {
    zdel <- rnorm(dat$J-1, 0, 2)
    zth1 <- rnorm(1, 0, sd=0.5)
    zgm <- runif(1, 0.2, 0.75)
    ztau <- runif(dat$J-1, 0.25, .75)

    zpsi <- rnorm(dat$S-1, 0, 4)
    zalph1 <- rnorm(1, 0, sd=0.5)
    zx <- runif(1, 0.4, 0.6)
    zrho <- runif(dat$S-2, 0.4, .6)
    zprho2 <- runif(1, 0.25, 0.75)
    list(zdelta =  zdel,
    ztheta1 = zth1,
    zgam = zgm,
    ztau = ztau,

    zpsi =  zpsi,
    zalpha1 = zalph1,
    zxi = zx,
    zprho2=zprho2,
    zrho = zrho
    )
  }'
  }


  if (prior.c=="horseshoe" & prior.s=="normal" & order.c==1 & order.s==2) {
    tmp.b <- 'function(chain_id=1, dat=tmp.dat) {
    zdel <- rnorm(dat$J-1, 0, 2)
    zth1 <- rnorm(1, 0, sd=0.5)
    zgm <- runif(1, 0.25, 0.75)
    ztau <- runif(dat$J-1, 0.25, .75)

    zpsi <- rnorm(dat$S-1, 0, 2)
    zalph1 <- rnorm(1, 0, sd=0.5)
    zx <- runif(1, 0.25, 0.75)
    list(zdelta =  zdel,
    ztheta1 = zth1,
    zgam = zgm,
    ztau = ztau,

    zpsi =  zpsi,
    zalpha1 = zalph1,
    zxi = zx
    )
  }'
  }

  ####  Order 2 both models

  if (prior.c=="normal" & prior.s=="horseshoe" & order.c==2 & order.s==2) {
    tmp.b <- 'function(chain_id=1, dat=tmp.dat) {
    zdel <- rnorm(dat$J-1, 0, 4)
    zth1 <- rnorm(1, 0, sd=0.5)
    zgm <- runif(1, 0.4, 0.6)

    zpsi <- rnorm(dat$S-1, 0, 4)
    zalph1 <- rnorm(1, 0, sd=0.5)
    zx <- runif(1, 0.4, 0.6)
    zrho <- runif(dat$S-2, 0.4, .6)
    zprho2 <- runif(1, 0.25, 0.75)
    list(zdelta =  zdel,
    ztheta1 = zth1,
    zgam = zgm,

    zpsi =  zpsi,
    zalpha1 = zalph1,
    zxi = zx,
    zprho2=zprho2,
    zrho = zrho
    )
  }'
  }


  if (prior.c=="normal" & prior.s=="normal" & order.c==2 & order.s==2) {
    tmp.b <- 'function(chain_id=1, dat=tmp.dat) {
    zdel <- rnorm(dat$J-1, 0, 4)
    zth1 <- rnorm(1, 0, sd=0.5)
    zgm <- runif(1, 0.4, 0.6)

    zpsi <- rnorm(dat$S-1, 0, 4)
    zalph1 <- rnorm(1, 0, sd=0.5)
    zx <- runif(1, 0.4, 0.6)

    list(zdelta =  zdel,
    ztheta1 = zth1,
    zgam = zgm,

    zpsi =  zpsi,
    zalpha1 = zalph1,
    zxi = zx
    )
  }'
  }


  if (prior.c=="horseshoe" & prior.s=="horseshoe" & order.c==2 & order.s==2) {
    tmp.b <- 'function(chain_id=1, dat=tmp.dat) {
    zdel <- rnorm(dat$J-1, 0, 4)
    zth1 <- rnorm(1, 0, sd=0.5)
    zpt2 <- runif(1, 0.25, 0.75)
    zgm <- runif(1, 0.4, 0.6)
    ztau <- runif(dat$J-2, 0.4, .6)

    zpsi <- rnorm(dat$S-1, 0, 4)
    zalph1 <- rnorm(1, 0, sd=0.5)
    zx <- runif(1, 0.4, 0.6)
    zrho <- runif(dat$S-2, 0.4, 0.6)
    zprho2 <- runif(1, 0.25, 0.75)
    list(zdelta =  zdel,
    ztheta1 = zth1,
    zgam = zgm,
    ztau = ztau,
    zptau2 = zpt2,

    zpsi =  zpsi,
    zalpha1 = zalph1,
    zxi = zx,
    zprho2=zprho2,
    zrho = zrho
    )
  }'
  }


  if (prior.c=="horseshoe" & prior.s=="normal" & order.c==2 & order.s==2) {
    tmp.b <- 'function(chain_id=1, dat=tmp.dat) {
    zdel <- rnorm(dat$J-1, 0, 4)
    zth1 <- rnorm(1, 0, sd=0.5)
    zpt2 <- runif(1, 0.25, 0.75)
    zgm <- runif(1, 0.4, 0.6)
    ztau <- runif(dat$J-2, 0.4, .6)

    zpsi <- rnorm(dat$S-1, 0, 4)
    zalph1 <- rnorm(1, 0, sd=0.5)
    zx <- runif(1, 0.4, 0.6)

    list(zdelta =  zdel,
    ztheta1 = zth1,
    zgam = zgm,
    ztau = ztau,
    zptau2 = zpt2,

    zpsi =  zpsi,
    zalpha1 = zalph1,
    zxi = zx
    )
  }'
  }




  ptmp <- parse(text=tmp.b)
  outf <- eval(ptmp)
  return(outf)

}  #end function


#' Prepare the format neeed for HMC
#'
#' @param samp_times vector of sampling times
#' @param n_sampled vector of number of samples per sampling time
#' @param coal_times vector of coalescent time
#' @param grid grid points
#'
#' @return all the information needed into the adapref_sampling function
#' @export

samp_and_coalescent_data <- function(samp_times, n_sampled, coal_times, grid)
{
  ns <- length(samp_times)
  nc <- length(coal_times)
  ng <- length(grid)-1

  if (sum(diff(coal_times)==0)>=1){
    id<-which(diff(coal_times)==0)
    corr<-min(diff(coal_times)[diff(coal_times)>0])*0.9
    coal_times[id+1]<-coal_times[id+1]+corr
    print("changed the coalescent times to avoid overlapping ones")
  }

  if (sum(diff(samp_times)==0)>=1){
    id<-which(diff(samp_times)==0)
    corr<-min(diff(samp_times)[diff(samp_times)>0])*0.9
    samp_times[id+1]<-samp_times[id+1]+samp_times
    print("changed the sampling times to avoid overlapping ones")
  }

  if (length(samp_times) != length(n_sampled))
    stop("samp_times vector of differing length than n_sampled vector.")

  if (length(coal_times) != sum(n_sampled) - 1)
    stop("Incorrect length of coal_times: should be sum(n_sampled) - 1.")


  if (max(samp_times, coal_times) > max(grid))
    stop("Grid does not envelop all sampling and/or coalescent times.")

  t <- sort(unique(c(samp_times, coal_times, grid)))  #combined times
  alin <- rep(0, length(t))   #number of active lineages

  for (i in 1:ns)
    alin[t >= samp_times[i]] <- alin[t >= samp_times[i]] + n_sampled[i]

  for (i in 1:nc)
    alin[t >= coal_times[i]] <- alin[t >= coal_times[i]] - 1

  #print(l)

  if (sum((alin < 1) & (t >= min(samp_times))) > 0)
    stop("Number of active lineages falls below 1 after the first sampling point.")

  mask <- alin > 0
  t <- t[mask]  # drop times with zero lineages
  alin <- head(alin[mask], -1) #drops first

  gridrep <- rep(0, ng)  #numbers of subintervals within a grid cell
  for (i in 1:ng)
    gridrep[i] <- sum(t > grid[i] & t <= grid[i+1])

  C <- 0.5 * alin * (alin-1)  #binomial coefficient for active lineages
  D <- diff(t)  #time step width

  y <- rep(0, length(D))  #indicator for coal event within sub interval
  y[t[-1] %in% coal_times] <- 1
  ny <- length(y)
  ncoal <- sum(y)

  #scount<-sum(n_sampled[1:max(which(samp_times<=sgrid$grid[1]))])
  scount<-c()
  for (i in 1:(length(sgrid$midpts))){
    id.low<-min(which(samp_times>=sgrid$grid[i]),length(samp_times)+1)
    if (id.low==length(samp_times)+1){break} #i.e. there are no more sampling times
    else{
      id.up<-max(which(samp_times<=sgrid$grid[i+1]))
      if (id.up<id.low){ #it happens when there are no sampling event in an interval
        scount<-c(scount,0)}  else {
          scount<-c(scount,sum(n_sampled[id.low:id.up]))}}
  }

  rep_idx <- cumsum(gridrep)
  rep_idx <- cbind(start=rep_idx-gridrep+1,end=rep_idx)

  coalind <- integer(ny)
  cnt <- 1
  for (j in 1:ny){
    coalind[j] <- cnt
    if (y[j]==1) cnt <- cnt+1
  }
  cstart <- integer(ncoal)
  cend <- integer(ncoal)
  for (k in 1:ncoal){
    tmpi <- which(coalind==k)
    cstart[k] <- min(tmpi)
    cend[k] <- max(tmpi)
  }
  mle_th <- -log(sum(y)/sum(C*D))
  mle_al <- log(ns/exp(mle_th)/sum(D))

  ncoalv <- numeric(ng)
  for (jj in 1:ng){
    ncoalv[jj] <- sum(y[rep_idx[jj,1]:rep_idx[jj,2]])
  }

  return(list(J=ng, N=ny, y=y, gridrep=gridrep, Aik=C, coalind=coalind,ncoal=ncoal,ncoalv=ncoalv, cstart=cstart,cend=cend, dalpha=D, rep.idx=rep_idx, log_mu=mle_th,log_beta=mle_al,scount=scount,dgrid=rep(grid[2],length(scount)),S=length(scount)))
}


