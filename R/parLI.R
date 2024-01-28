#######################################
#### Funcion para los limites de pareto
#######################################

#' .
#'
#' @param S .
#' @param x .
#' @param sim .
#' @param eval .
#' @param cum .
#'
#' @return
#'
#'
#' @examples

parLI <- function(S, x, sim, eval, cum){


  mu1 <- S[4]
  s1 <- S[2]
  s2 <- S[3]
  fS1 <- logf(s1)
  fS2 <- logf(s2)
  a1 <- log(1/(abs(s1 - mu1)))
  b1 <- log(1/(abs(s2 - mu1)))
  gamm1 <- (fS1 - fS2)/(a1 - b1)
  rho1 <- fS2 - (gamm1 * b1)

  m <- 1/((exp(rho1)/(gamm1 -1))*(((mu1 - S[2])^{1 - gamm1})  - ((mu1 - S[1])^{1 - gamm1})))

  if (eval){
    return(m*exp(rho1)*(abs(x - mu1))^{-gamm1})
  }

  if (sim){
    u <- runif(1)
    sim <- mu1 - (((gamm1 - 1)*u/(m*exp(rho1))) + (mu1 - S[1])^{1 - gamm1})^{1/(1-gamm1)}
    return(sim)
  }

  if (cum){
    return(exp(rho1)*(abs(x - mu1))^{-gamm1})
  }
}
