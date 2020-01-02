#' @title A function to compute a Monte Carlo estimate of the Beta cdf with integer parameters
#' @description A function to compute a Monte Carlo estimate of the Beta cdf with integer parameters
#' @param y the parameter of cdf
#' @param alpha the  shape paramaters of beta distribution
#' @param beta the  shape paramaters of beta distribution
#' @param m the number of Monte Carlo experiment times
#' @return a Monte Carlo estimate of the beta cdf
#' @examples
#' \dontrun{
#' alpha<-3;beta<-3;m<-1e4
#' Pbeta1<-function(x){
#'  Pbeta(x,alpha,beta,m)
#' }
#' sapply((1:9)/10,Pbeta1)
#' pbeta1<-function(x){
#'  pbeta(x,alpha,beta)
#' }
#' sapply((1:9)/10,pbeta1)
#' }
#' @export
Pbeta<-function(y,alpha,beta,m){
  if(y<=0) value=0
  if(y>1) value=1
  if(y>=0&y<1){
    m <- 1e4; x <- runif(m, min=0, max=y)
    value <- mean(factorial(alpha+beta-1)*x^2*(1-x)^2/(factorial(alpha-1)*factorial(beta-1)))*(y-0)
  }
  value
}
