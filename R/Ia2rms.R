#######################################
#### Funcion para los limites de pareto
#######################################

#' Title
#'
#' @param S .
#' @param x .
#' @param sim .
#' @param eval .
#' @param cum .
#'
#' @return .
#'
#'
#' @examples .
#'
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



#' Title.
#'
#' @param S .
#' @param x .
#' @param eval .
#' @param sim .
#' @param cum .
#'
#' @return .
#'
#'
#' @examples .
parLS <- function(S, x, eval, sim, cum){

  j <- length(S)

  mu2 <- S[j-3]
  ns1 <- S[j-2]
  ns2 <- S[j-1]
  nfS1 <- logf(ns1)
  nfS2 <- logf(ns2)
  an <- log(1/(abs(ns1 - mu2)))
  bn <- log(1/(abs(ns2 - mu2)))
  gammn <- (nfS1 - nfS2)/(an - bn)
  rho2 <- nfS2 - (gammn * bn)

  m <- 1/((exp(rho2)/(1 - gammn))*((S[j] - mu2)^{1 - gammn}  - (ns2 - mu2)^{1 - gammn}))

  if (eval){
    return(m*exp(rho2)*(abs(x - mu2))^{-gammn})
  }

  if (sim){
    u <- runif(1)
    sim <- mu2 + ((S[j - 1] - mu2)^{1 - gammn} - ((gammn - 1)* u/(m * exp(rho2))))^{1/(1 - gammn)}
    return(sim)
  }

  if(cum){
    return(exp(rho2)*(abs(x - mu2))^{-gammn})
  }
}







####################################
###### Función que calcula F(x)
####################################

#' Title
#'
#' @param S .
#' @param logf .
#'
#' @return .
#'
#'
#' @examples .
#'
Fx <- function(S, logf){

  j <- length(S)

  a <- NULL

  intervalos <- NULL

  for (i in 2:j ){

    intervalos <- rbind(intervalos, c(S[i-1], S[i]))

    if ((i > 2) & (i < j) ) a <- c(a, max(logf(intervalos[(i-1),1]), logf(intervalos[(i-1),2])))
  }

  intI <- integrate(parLI, lower = S[1], upper = S[2], S = S,eval = F, sim = F,  cum = T)
  intS <- integrate(parLS, lower = S[j-1], upper = S[j], S = S,eval = F, sim = F,  cum = T)

  m <- sum(exp(a)* (intervalos[2:(j-2),2] - intervalos[2:(j-2),1])) + intI$value + intS$value


  Prob <- matrix(c(0, intI$value), ncol = 2)

  dim <- dim(intervalos)

  for (i in 2:(dim[1] - 1)){

    li <- Prob[(i-1),2]
    ls <- li + exp(a[i-1])*(S[i+1] - S[i])
    Prob <- rbind(Prob, c(li, ls))

    if(i == (dim[1] - 1)){
      li <- Prob[i,2]
      ls <- li + intS$value
      Prob <- (1/m)*rbind(Prob, c(li, ls))
    }
  }
  return(list(Prob, intervalos, a))
}



######################################
#FUNCION PARA SIMULAR DE H(X)
######################################

#' Title
#'
#' @param S .
#' @param logf .
#' @param Probs .
#' @param num.sim .
#'
#' @return .
#'
#'
#' @examples .
simH <- function(S, logf, Probs, num.sim){

  values <- NULL


  for (i in (1:num.sim)){
    u <- runif(1)
    bol <- c()
    for (j in (1:dim(Probs)[1])){
      bol <- c(bol, between(x = u, left = Probs[j ,1], right = Probs[j ,2]))
    }

    r <- which(bol)
    if (r == 1) {
      v <- parLI(S, sim = T, eval = F, cum = F)
      values <- c(values,v)
    }
    if (r == (length(S) - 1)){
      v <- parLS(S, sim = T, eval = F, cum = F)
      values <- c(values,v)
    }
    if (all(r != c(1, (length(S) - 1)))){
      v <-  runif(1, S[r], S[r+1])
      values <- c(values,v)
    }
  }
  return(values)
}


###############################################################
######## Funcion para evaluar h(x)
#######################################




#' Title
#'
#' @param x .
#' @param parLI .
#' @param parLS .
#' @param S .
#' @param intervalos .
#' @param a .
#'
#' @return .
#'
#'
#' @examples .
evalH <- function(x, parLI, parLS, S, intervalos, a){

  j <- length(S)

  bol <- c()

  for (i in (1:dim(intervalos)[1])){
    bol <- c(bol, between(x = x, left = intervalos[i ,1], right = intervalos[i ,2]))
  }

  bol <- which(bol)

  if (bol == 1) return(parLI(S, x = x, sim = F, eval = F, cum = T))
  if (bol == (j - 1)){ return(parLS(S, x = x, sim = F, eval = F, cum = T))
  }else{
    return(exp(a[bol-1]))
  }
}



######################################
####### IA2RMS  #######################
######################################

#' Independent Doubly Adaptive Rejection Metropolis Sampling
#'
#' @param y.initial Value inside the domain of the target density
#' @param logf log of the target density. (Class function)
#' @param indFunc Indicator function for the domain of the target density. (Class function)
#' @param n.sample Number of samples to draw from the taret density
#'
#' @return Sample of lenght n.sample

#' @importFrom stats runif integrate
#' @importFrom dplyr between
#' @importFrom dlm convex.bounds

#' @examples
#' y.initial <- 0.8
#' n.sample <- 5000
#' logf <- function(x) log(0.3 * dt(x, 3.2) + 0.35*dnorm(x,mean = 2.5, 0.3) + 0.5*dt(x, 1))
#' indFunc <- function(x) all(x > -5) * all(x < 6)
#' y <- Ia2rms(y.initial, logf, indFunc, n.sample)
#' hist(y[[1]],breaks = 50, probability = T, ylim = c(0,0.5))
#' curve(exp(Vectorize(logf)(x)), -5, 8, ylab = "f(x)", main = "Representación de f(x) y g(x)", add = T)
#' legend(legend = c("f(x)", "g(x)"), col = c("black", "blue"), x = "topleft", lty = c(1, 1))
#'
#' @export
#'
Ia2rms <- function(y.initial, logf, indFunc, n.sample){

  dim <- length(y.initial)


  if (dim == 1) {
    bounds <- y.initial + convex.bounds(y.initial, dir = 1, indFunc = indFunc)

    if (diff(bounds) < 1e-07){
      y.sample <- rep(y.initial, n.sample)}
  }

  n <- 15

  delta <- (bounds[2] - bounds[1]) / n

  S <- bounds[1] + (0:n) * delta

  path <- Fx(S, logf)

  Probs <- path[[1]]

  y.sample <- NULL


  x0 <- (S[3] + S[4])/2

  k <- 0

  while (length(y.sample) < n.sample){

    # Paso Aceptación Rechazo

    k <- k + 1
    x <- simH(S = S, logf = logf, Probs = Probs, num.sim = 1)
    pix <- evalH(x, parLI, parLS, S, path[[2]], path[[3]])

    if (runif(1) <= exp(logf(x))/pix){

      alpha <- min(1, (exp(logf(x)) * min(exp(logf(x0)), evalH(x0, parLI, parLS, S, path[[2]], path[[3]])))/ (exp(logf(x0)) * min(exp(logf(x)), pix)))


      if (runif(1) <= alpha){
        y.sample <- c(y.sample, x)
        y <- x0
        x0 <- x

      }else{
        x0 <- x0
        y <- x
      }

      if (runif(1) > evalH(y, parLI, parLS, S, path[[2]], path[[3]])/exp(logf(y)) & (length(S) < 80)){

        S <- sort(c(y,S))
        path <- Fx(S, logf)
        Probs <- path[[1]]
      }

    }else{
      if(length(S) < 80){
        S <- sort(c(x,S))
        path <- Fx(S, logf)
        Probs <- path[[1]]}
    }
  }
  return(y.sample)
}
