#From phylodyn

shade_band = function(x, ylo, yhi, xlim=NULL, col="gray")
{
  if (is.null(xlim))
    xlim = c(0, Inf)
  mask = x >= min(xlim) & x <= max(xlim)
  
  x = x[mask]
  ylo = ylo[mask]
  yhi = yhi[mask]
  
  graphics::polygon(c(x, rev(x)), c(yhi, rev(ylo)), col=col, border=NA)
}

float2gray = function(f)
{
  return(sprintf("#%02x%02x%02x", floor((1-f) * 255), floor((1-f) * 255), floor((1-f) * 255)))
}

#' Plot heatmap of a histogram
#' 
#' @param hist \code{histogram} object to be displayed.
#'
#' @param y numeric y-coordinate to display heatmap.
#' @param wd numeric width of heatmap in y-units.
#'
#' @export
hist2heat = function(hist, y, wd)
{
  breaks = hist$breaks
  counts = hist$counts
  upper  = max(counts)
  n = length(counts)
  cols = float2gray(counts / upper)
  graphics::segments(x0 = breaks[1:n], y0=y, x1 = breaks[2:(n+1)], y1=y, lwd=wd, col=cols, lend=1)
}


#From spmrf


varRef <- function(nn, kap, omg2, order=1){
  #nn is number of thetas, kap = 1/gam^2 (precision), omg2=omega^2 (variance of theta1)
  #returns vector of marginal variances
  nnv <- 1:nn
  if (order==1) out <- omg2 + (nnv-1)*(1/kap)
  if (order==2) out <- omg2 + (1/kap)*nnv*(nnv-1)*(2*nnv-1)/6
  if (order >= 3) out <- NULL
  return(out)
}

