

#' adaPop: Inference for effective population size modeling two populations jointly (INLA approximation)
#'
#' @param tree1 \code{phylo} object or list containing vectors of coalescent
#'   times \code{coal_times}, sampling times \code{samp_times}, and number
#'   sampled per sampling time \code{n_sampled}.
#' @param tree2 Same as the previous parameters, now related to the second population.
#' @param samp_times1 vector of sampling times related to tree1. It necessary in the case that tips are not at zero
#' @param samp_times2 vector of sampling times related to tree2. It necessary in the case that tips are not at zero
#' @param lengthout number of grid points.
#' @param prec_alpha hyperparameter GMRF on Ne
#' @param prec_beta hyperparameter GMRF on gamma
#' @param beta1_prec hyperparameter GMRF on beta (pref sampling)
#' @param parSel Default is FALSE. If TRUE uses the parametric model, if FALSE, the adaptive one.
#' @param preferential. Default is FALSE. If TRUE includes also preferential sampling, if FALSE, preferential sampling is not included
#' @u.truncation It determines the upper bound the grid of the GMRF. Default is the minimum TMRCA between the two populations.
#' @l.truncation It determines the lower bound the grid of the GMRF. Default is 0
#' @beta0_remove It removes the intercept in the parameteric model if TRUE. Default is FALSE
#' @param simplify logical whether to fully bucket all Poisson points.
#' @param derivative logical whether to calculate estimates of the
#'   log-derivative.
#' @param forward logical whether to use the finite difference approximations of
#'   the log-derivative as a forward or backward derivative.
#'
#' @return Effective pop size posterior, and posterior distributions of additional parameters.
#'  \code{result} contains the INLA output, \code{data} contains the
#'   information passed to INLA, \code{grid}  grid boundaries,
#'   \code{x}  grid mid points \code{effpop} median effective population size estimates,
#'   \code{effpop025} and \code{effpop975}  2.5th and 97.5th
#'   posterior percentiles, \code{summary} contains a data.frame of the
#'   estimates
#'   If \code{parSel = FALSE}, \code{selInt} median gamma estimates,
#'   \code{selInt025} and \code{selInt975}  2.5th and 97.5th
#'   posterior percentiles, \code{selIntsummary} contains a data.frame of the
#'   estimates.
#'   If \code{parSel = TRUE}, \code{beta1} median beta estimates,
#'   \code{beta1post} beta posterior distribution and \code{beta1summ}.
#'   Same information are included also for the other parameter (if \code{beta0_remove = FALSE})
#'    and \code{derivative} (if \code{derivative = TRUE})  log-derivative,
#'   If \code{preferential = TRUE}, \code{sampInt} median beta, \code{sampIntmean} mean beta,
#'   \code{sampInt975} and \code{sampInt025}97.5th and 2.5th posterior percentile beta
#' @export


BNPR_pop <- function (tree1,tree2, samp_times1,samp_times2, lengthout = 100,
                      prec_alpha = 0.01,
                      prec_beta = 0.01,
                      beta1_prec = 0.001,
                      parSel=FALSE,
                      preferential=FALSE,
                      u.truncation=NA,
                      l.truncation=NA,
                      beta0_remove=FALSE,
                      simplify = TRUE,
                      derivative = FALSE,
                      forward = TRUE)
{

  if (min(samp_times1)>0){
    warning("Warning: group 1 was not sampled at 0, model may not be identifiable +
            the current implementation is suboptimal in that case",immediate. = T)
  }

  #If parSel is FALSE, I am using adaSel
  phy1 <- summarize_phylo(tree1)
  phy1$samp_times <- phy1$samp_times + min(samp_times1)
  phy1$coal_times <- phy1$coal_times + min(samp_times1)
  phy2 <- summarize_phylo(tree2)
  phy2$samp_times <- phy2$samp_times + min(samp_times2)
  phy2$coal_times <- phy2$coal_times + min(samp_times2)


  result <- infer_coal_samp_pop(phy1,phy2,lengthout = lengthout,
                                      prec_alpha = prec_alpha, prec_beta = prec_beta,
                                      beta1_prec = beta1_prec,
                                      parSel, preferential,
                                      u.truncation,l.truncation,
                                      simplify,beta0_remove)



  #result$samp_times <- phy$samp_times
  #result$n_sampled <- phy$n_sampled
  #result$coal_times <- phy$coal_times
  if (!preferential){
      result$effpop <- exp(-result$result$summary.random$time$`0.5quant`)
      result$effpopmean <- exp(-result$result$summary.random$time$mean)
      result$effpop975 <- exp(-result$result$summary.random$time$`0.025quant`)
      result$effpop025 <- exp(-result$result$summary.random$time$`0.975quant`)
      result$summary <- with(result$result$summary.random$time,
                             data.frame(time = ID, mean = exp(-mean), sd = sd * exp(-mean),
                                        quant0.025 = exp(-`0.975quant`), quant0.5 = exp(-`0.5quant`),
                                        quant0.975 = exp(-`0.025quant`)))
    } else {
    result$effpop <- exp(-result$result$summary.random$time1$`0.5quant`)
    result$effpopmean <- exp(-result$result$summary.random$time1$mean)
    result$effpop975 <- exp(-result$result$summary.random$time1$`0.025quant`)
    result$effpop025 <- exp(-result$result$summary.random$time1$`0.975quant`)
    result$summary <- with(result$result$summary.random$time1,
                           data.frame(time = ID, mean = exp(-mean), sd = sd * exp(-mean),
                                      quant0.025 = exp(-`0.975quant`), quant0.5 = exp(-`0.5quant`),
                                      quant0.975 = exp(-`0.025quant`)))
    result$sampInt <- exp(result$result$summary.random$time3b$`0.5quant`)
    result$sampIntmean <- exp(result$result$summary.random$time3b$mean)
    result$sampInt975 <- exp(result$result$summary.random$time3b$`0.975quant`)
    result$sampInt025 <- exp(result$result$summary.random$time3b$`0.025quant`)
    result$sampIntsummary <- with(result$result$summary.random$time3b,
                                 data.frame(time = ID, mean = exp(mean), sd = sd * exp(mean),
                                            quant0.025 = exp(`0.025quant`), quant0.5 = exp(`0.5quant`),
                                            quant0.975 = exp(`0.975quant`)))
    }

    if (parSel){
      if (beta0_remove==FALSE){
        result$beta0 <- result$result$summary.fixed["beta0", "0.5quant"]
        result$beta0summ <- result$result$summary.fixed["beta0",]
        rownames(result$beta0summ) <- "Beta0"
        result$beta0post <- result$result$marginals.fixed$beta0
      }
      result$beta1 <- result$result$summary.hyperpar[2, "0.5quant"]
      result$beta1summ <- result$result$summary.hyperpar[2,]
      rownames(result$beta1summ) <- "Beta1"
      result$beta1post <- result$result$marginals.hyperpar$`Beta for time2`
    } else {
      result$selInt <- exp(-result$result$summary.random$time2b$`0.5quant`)
      result$selIntmean <- exp(-result$result$summary.random$time2b$mean)
      result$selInt975 <- exp(-result$result$summary.random$time2b$`0.025quant`)
      result$selInt025 <- exp(-result$result$summary.random$time2b$`0.975quant`)
      result$selIntsummary <- with(result$result$summary.random$time2b,
                                   data.frame(time = ID, mean = exp(mean), sd = sd * exp(mean),
                                              quant0.025 = exp(`0.025quant`), quant0.5 = exp(`0.5quant`),
                                              quant0.975 = exp(`0.975quant`)))
    }

  result$preferential=preferential

  return(result)
}





infer_coal_samp_pop <- function(phy1,phy2, lengthout=100, prec_alpha=0.01, prec_beta=0.01,
                                      beta1_prec=0.001, parSel=TRUE,preferential=FALSE,
                                      u.truncation=NA,l.truncation=NA, simplify = TRUE,
                                      beta0_remove=FALSE)
{
  if (!requireNamespace("INLA", quietly = TRUE)) {
    stop('INLA needed for this function to work. Use install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE).',
         call. = FALSE)
  }


  coal_times1 <- phy1$coal_times
  coal_times2 <- phy2$coal_times
  n_sampled1 <- phy1$n_sampled
  n_sampled2 <- phy2$n_sampled
  samp_times1 <- phy1$samp_times
  samp_times2 <- phy2$samp_times

  grid <- seq(min(c(samp_times1,samp_times2)), max(c(coal_times1,coal_times2)), length.out = lengthout+1)

  #if (is.null(n_sampled))
  #  n_sampled <- rep(1, length(samp_times))

  coal_data1 <- coal_stats(grid = grid, samp_times = samp_times1, n_sampled = n_sampled1,
                           coal_times = coal_times1)
  if (simplify){coal_data1 <- with(coal_data1, condense_stats(time=time, event=event, E=E))}
  coal_data1 <- identify_off_grid(coal_data1,samp_times1,coal_times1)

  coal_data2 <- coal_stats(grid = grid, samp_times = samp_times2, n_sampled = n_sampled2,
                           coal_times = coal_times2)
  if (simplify){coal_data2 <- with(coal_data2, condense_stats(time=time, event=event, E=E))}
  coal_data2 <- identify_off_grid(coal_data2,samp_times2,coal_times2)


  if (preferential){
    samp_data1 <- samp_stats_adasel(grid = grid, samp_times = samp_times1,
                                n_sampled = n_sampled1)
    samp_data2 <- samp_stats_adasel(grid = grid, samp_times = samp_times2,
                                 n_sampled = n_sampled2)
  }




  #If I am using a truncation, this is the place that adjust for it
  if (!is.na(u.truncation)){
    id <- min(which(coal_data1$time>=u.truncation))
    coal_data1 <- coal_data1[1:id,]
    coal_data2 <- coal_data2[1:id,]
    grid <- grid[1:(id+1)]
  }
  if (!is.na(l.truncation)){
    id <- max(which(coal_data1$time<=l.truncation))
    coal_data1 <- coal_data1[id:length(coal_data1$time),]
    coal_data2 <- coal_data2[id:length(coal_data1$time),]
    grid <- grid[(id-1):length(grid)]
  }



  hyper <- list(prec = list(param = c(prec_alpha, prec_beta)))


  if (!preferential){
    data <- joint_coal_stats(coal_data1 = coal_data1, coal_data2 = coal_data2,
                             samp_times2=samp_times2,grid=grid)
  } else {
    data <- joint_coal_stats_adasel(coal_data1 = coal_data1, coal_data2 = coal_data2,
                                    samp_data1 = samp_data1, samp_data2 = samp_data2,
                                    grid=grid,samp_times2=samp_times2)

  }
  if (parSel){
    if (!preferential){
        if (beta0_remove){
          formula <- Y ~ -1 +
            f(time, model="rw1", hyper = hyper, constr = FALSE) +
            f(time2,copy="time", fixed=FALSE, param=c(0, beta1_prec))
        } else {
          formula <- Y ~ -1 + beta0 +
            f(time, model="rw1", hyper = hyper, constr = FALSE) +
            f(time2,copy="time", fixed=FALSE, param=c(0, beta1_prec))
        }
    } else {
      data$time3b=data$time3
      data$time4b=data$time4
      formula <- Y ~ -1 + beta0 +
        f(time1, model="rw1", hyper = hyper, constr = FALSE) +
        f(time2,copy="time1", fixed=FALSE, param=c(0, beta1_prec)) +
        f(time3,w,copy="time1") + f(time3b, model="rw1", hyper = hyper, constr=FALSE) +
        #f(time4,w,copy="time2") + f(time4b, copy="time3b")
        f(time4,w,copy="time2") + f(time4b, model="rw1", hyper = hyper, constr=FALSE)
    }
  } else {
    if (!preferential){
      data$beta0<-NULL
      data$time2b=data$time2
      formula <- Y ~ -1 +
        f(time, model="rw1", hyper = hyper, constr = FALSE) +
        f(time2, w, copy="time") + f(time2b, model="rw1",hyper = hyper, constr=FALSE)
    } else {
      data$beta0<-NULL
      data$time2b=data$time2
      data$time3b=data$time3
      data$time4b=data$time4
      data$time4c=data$time4
      formula <- Y ~ -1 +
        f(time1, model="rw1", hyper = hyper, constr = FALSE) +
        f(time2, copy="time1") + f(time2b, model="rw1", hyper = hyper, constr=FALSE) +
        f(time3, w, copy="time1") + f(time3b, model="rw1", hyper = hyper, constr=FALSE) +
        #f(time4, w, copy="time1") + f(time4b, w, copy="time2b") + f(time4c,copy="time3b")
        f(time4, w, copy="time1") + f(time4b, w, copy="time2b") +f(time4c, model="rw1", hyper = hyper, constr=FALSE)
    }

  }

  if(!preferential){
    family <- c("poisson", "poisson")
  } else {
    family <- c("poisson", "poisson","poisson","poisson")
  }

  lc_many <- NULL
  mod <- INLA::inla(formula, family = family, data = data,
                    lincomb = lc_many, offset = data$E_log,
                    #control.predictor = list(compute=TRUE),
                    control.compute = list(config=TRUE))
                    #control.inla = list(lincomb.derived.only=FALSE)
  #I am removing the control.predictor. Let's see if I need it.

  return(list(result = mod, data = data, grid = grid, x = coal_data1$time))
}



identify_off_grid <- function(coal_data,samp_times,coal_times)
{
  if (samp_times[1]>0){
    id <- which.min(abs(coal_data$time-samp_times[1]))
    coal_data$event[1:(id-1)] <- NA
  }
  if (coal_times[length(coal_times)]<coal_data$time[length(coal_data$time)]){
    id <- which.min(abs(coal_data$time-coal_times[length(coal_times)])) #Note that there are double points in $time, that's why id+1
    coal_data$event[(id+1):length(coal_data$event)] <- NA
  }
  return(coal_data)
}


joint_coal_stats <- function(coal_data1, coal_data2,samp_times2,grid)
{

  ##Here, I ensure that the field on the second population starts from the minimum sampling points
  id.field2 <- max(which(utils::tail(grid,-1) <=  min(samp_times2)),0) + 1
  coal_data2 <- coal_data2[id.field2:nrow(coal_data2),]

  n1 <- length(coal_data1$time) #Here are the midpts
  n2 <- length(coal_data2$time) #
  beta0 <- c(rep(0, n1), rep(1, n2))
  E_log <- c(coal_data1$E_log, coal_data2$E_log) #samp_data$E_log does not include NA in the pref_samp, so I am not modifying it here either
  Y <- matrix(c(coal_data1$event, rep(NA, n2), rep(NA, n1), coal_data2$event),
              nrow = n1 + n2, byrow = FALSE)
  #Need to correct Y for the parts of the grid that need to be eccluded
  #Note, there are some parts of Y where they are both NAs
  w <- c(rep(1, n1), rep(1, n2))  #what does the w do? I think that takes into account for the fact that one is at the numerator,
  #the other one at the numerator.

  #we are not correcting the time, i.e. there are times also to part that are not defined
  time  <- c(coal_data1$time, rep(NA, n2))
  time2 <- c(rep(NA, n1), coal_data2$time)

  return(list(Y = Y, beta0 = beta0, time = time, time2 = time2, w = w, E_log = E_log))
}


joint_coal_stats_adasel <- function(coal_data1, coal_data2,samp_data1,samp_data2,grid,samp_times2)
{

  id.field2 <- max(which(utils::tail(grid,-1) <=  min(samp_times2)),0) + 1
  coal_data2 <- coal_data2[id.field2:nrow(coal_data2),]
  samp_data2 <- samp_data2[id.field2:nrow(samp_data2),]


  n1 <- length(coal_data1$time) #Here are the midpts
  n2 <- length(coal_data2$time) #
  n3 <- length(samp_data1$time)
  n4 <- length(samp_data2$time)
  beta0 <- c(rep(0, n1), rep(1, n2),rep(0,n3),rep(1,n4))
  E_log <- c(coal_data1$E_log, coal_data2$E_log,samp_data1$E_log,samp_data2$E_log) #samp_data$E_log does not include NA in the pref_samp, so I am not modifying it here either
  c1 <- c(coal_data1$event, rep(NA, n2+n3+n4))
  c2 <- c(rep(NA, n1), coal_data2$event, rep(NA,n3+n4))
  c3 <- c(rep(NA,n1+n2),samp_data1$count, rep(NA,n4))
  c4 <- c(rep(NA,n1+n2+n3), samp_data2$count)
  Y <- cbind(c1,c2,c3,c4)
  #Need to correct Y for the parts of the grid that need to be eccluded
  #Note, there are some parts of Y where they are both NAs
  w <- c(rep(1, n1), rep(1, n2),rep(-1,n3),rep(-1,n4))  #what does the w do? I think that takes into account for the fact that one is at the numerator,
  #the other one at the numerator.

  #we are not correcting the time, i.e. there are times also to part that are not defined
  time1 <- c(coal_data1$time, rep(NA, n2+n3+n4))
  time2 <- c(rep(NA, n1), coal_data2$time, rep(NA,n3+n4))
  time3 <- c(rep(NA,n1+n2),samp_data1$time, rep(NA,n4))
  time4 <- c(rep(NA,n1+n2+n3), samp_data2$time)


  return(list(Y = Y, beta0 = beta0, time1 = time1, time2 = time2, time3=time3,time4=time4, w = w, E_log = E_log))
}




coal_stats <- function(grid, samp_times, coal_times, n_sampled = NULL,
                       log_zero = -100)
{
  lengthout <- length(grid) - 1
  field <- grid[-1] - diff(grid)/2 #It is practically the mid point of the grid

  #if (is.null(n_sampled))
  #  n_sampled <- rep(1, length(samp_times))
  args <- gen_INLA_args(samp_times = samp_times, n_sampled = n_sampled,
                        coal_times = coal_times) #Looks identical between, just shifted
  #args refer only to the period in which I have coal and samp events. It does not matter the grid. It simply sort them relatively


  coal_factor <- args$coal_factor
  s <- args$s
  event <- args$event

  grid_trimmed <- setdiff(x = grid, y = s) #elements in the grid that are not in the other vector
  sorting <- sort(c(grid_trimmed, s), index.return=TRUE)
  sgrid <- sorting$x #grid + samp_events + coal_events - duplicates (e.g. time=0)
  ordering <- sorting$ix # the numbers before the lenght of the grid are grid points, the rest no

  time_index <- cut(x = sgrid[-1], breaks = grid, labels = FALSE) # It tells me which point of sgrid belong to a given interval in the grid
  #it is of the same length of sgrid but without point 0
  time <- field[time_index] #repeats all the midpoint in the grid

  event_out <- c(rep(0, length(grid_trimmed)), event)[ordering] #spread out the location of the coal_event (1s) in the sgrid,
  #grid point and sampling event are 0s

  Cfun <- stats::stepfun(x = s, y = c(0, coal_factor, 0), right = TRUE)
  Cvec <- Cfun(sgrid[-1])
  E <- diff(sgrid)*Cvec

  E_log = log(E)
  E_log[E == 0] = log_zero

  return(data.frame(time = time, event = event_out[-1], E = E, E_log = E_log))
}

gen_INLA_args <- function(samp_times, n_sampled, coal_times)
{
  if (sum(n_sampled) != length(coal_times) + 1)
    stop("Number sampled not equal to number of coalescent events + 1.")

  if (length(intersect(coal_times, samp_times)) > 0)
    warning("Coincident sampling event and coalescent event: results may be unpredictable.")

  l <- length(samp_times)
  m <- length(coal_times)
  sorting <- sort(c(samp_times, coal_times), index.return=TRUE)

  lineage_change <- c(n_sampled, rep(-1, m))[sorting$ix]
  lineages <- utils::head(cumsum(lineage_change), -1) # remove entry for the post-final-coalescent-event open interval
  coal_factor <- lineages*(lineages-1)/2

  event <- c(rep(0, l), rep(1, m))[sorting$ix]

  return(list(coal_factor=coal_factor, s=sorting$x, event=event, lineages=lineages))
}



condense_stats <- function(time, event, E, log_zero = -100)
{
  result <- stats::aggregate(event ~ time, FUN = sum)
  result$E <- stats::aggregate(E ~ time, FUN = sum)$E

  E_log = log(result$E)
  E_log[result$E == 0] = log_zero
  result$E_log <- E_log

  return(result)
}








samp_stats_adasel <- function(grid, samp_times, n_sampled = NULL, trim_end = FALSE)
{
  lengthout <- length(grid) - 1
  field <- grid[-1] - diff(grid)/2
  E <- diff(grid)

  bins <- cut(x = samp_times, breaks = grid, include.lowest = TRUE)

  if (is.null(n_sampled))
    count <- as.vector(table(bins))
  else
  {
    tab <- stats::aggregate(n_sampled ~ bins, FUN = sum, labels = FALSE)
    count <- rep(0, lengthout)
    count[as.numeric(tab$bins)] <- tab$n_sampled
  }

  count[utils::head(grid, -1) >= max(samp_times)] <- NA
  count[utils::tail(grid,-1) <=  min(samp_times)] <- NA #Double Check
  result <- data.frame(time = field, count = count, E = E, E_log = log(E))

  if (trim_end) #I am not sure what's this is supposed to do
    result <- result[stats::complete.cases(result),]


  return(result)
}

###### For future: make the following publicly accessible



eff2_adapop <- function(res,preferential=F){

  #Sample population size of the second group from the joint inference
  result <- res$result
  l_time1 <- nrow(res$summary)
  l_time2 <- nrow(res$selIntsummary)
  if (!preferential){
    samp <- INLA::inla.posterior.sample(500, res$result,
                                        selection = list(time=seq((l_time1-l_time2)+1,l_time1),
                                                         time2b=seq(1,l_time2)))
    posterior <- sapply(samp,function(x) apply(exp(matrix(c(x$latent),ncol=2,byrow=F))^(-1),1,prod))
  } else {
    samp <- INLA::inla.posterior.sample(500, res$result,
                                        selection = list(time1=seq((l_time1-l_time2)+1,l_time1),
                                                         time2b=seq(1,l_time2)))
    posterior <- sapply(samp,function(x) apply(exp(matrix(c(x$latent),ncol=2,byrow=F))^(-1),1,prod))

  }
  median <- apply(posterior,1,quantile, probs=0.5)
  quant0.025 <- apply(posterior,1,quantile, probs=0.025)
  quant0.975 <- apply(posterior,1,quantile, probs=0.975)
  mean <- apply(posterior,1,mean)
  result$effpop2 <- median
  result$effpop2mean <- mean
  result$effpop2_975 <- quant0.975
  result$effpop2_025 <- quant0.025
  result$effpop2summary <- data.frame(time = res$selIntsummary$time, mean = mean,
                                      quant0.025 = quant0.025, quant0.5 = median,
                                      quant0.975 = quant0.975)

  return(result)
}




eff2_parametric <- function(res){

  #compute the ratio of the parameteric model

  l_time <- nrow(res$summary)
  samp1 <- INLA::inla.posterior.sample(500, res$result,selection = list(time=seq(1,l_time)))
  posterior1 <- sapply(samp1,function(x) exp(matrix(c(x$latent),ncol=1,byrow=F))^(-1))
  samp2 <- INLA::inla.posterior.sample(500, res$result,selection = list(time2=seq(1,l_time),beta0=1))
  posterior2 <- sapply(samp2,function(x) exp(matrix(c(x$latent),ncol=1,byrow=F))^(-1))
  posterior2 <- posterior2[1:l_time,]%*%diag(posterior2[l_time+1,])
  posterior <- posterior2
  median <- apply(posterior,1,quantile, probs=0.5)
  quant0.025 <- apply(posterior,1,quantile, probs=0.025)
  quant0.975 <- apply(posterior,1,quantile, probs=0.975)
  mean <- apply(posterior,1,mean)
  result <- list()
  result$effpop2 <- median
  result$effpop2mean <- mean
  result$effpop2_975 <- quant0.975
  result$effpop2_025 <- quant0.025
  result$effpop2summary <- data.frame(time = res$summary$time, mean = mean,
                                      quant0.025 = quant0.025, quant0.5 = median,
                                      quant0.975 = quant0.975)

  return(result)
}


