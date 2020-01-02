#' @title Implement a random walk Metropolis sampler for generating the standard Laplace distribution
#' @description Implement a random walk Metropolis sampler for generating the standard Laplace distribution,For the increment, simulate from a normal distribution.
#' @param sigma the  shape paramaters of Laplace distribution
#' @param N the length of the series you wanted
#' @return a series satisfy Laplace distribution
#' @importFrom GeneralizedHyperbolic dskewlap
#' @examples
#' \dontrun{
#' N<-1e4
#' sigma<-4
#' laplace1<-Laplace.Metropolis(sigma,N)
#' }
#' @export

Laplace.Metropolis<-function(sigma,N){
  x<-numeric(N)
  x[1]<-rnorm(1,0,sigma)
  u<-runif(N)
  for (i in 2:N) {
    y<-rnorm(1,x[i-1],sigma)
    if(u[i]<=dskewlap(y)/dskewlap(x[i-1]))
      x[i]=y
    else{
      x[i]=x[i-1]
    }
  }
  return(x)
}
