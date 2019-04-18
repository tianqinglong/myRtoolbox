#' Compute The MLEs of the Weibull Censored Data
#'
#' Find out the MLEs for censored data with \code{survival} package.
#'
#' @param censor_data The value of function \code{simulate_weibull_data()}.
#' @examples
#' censor_data <- simulate_weibull_data(1, 5, 0.2, 1, 1)
#' ( censor_data_with_mle <- add_Weibull_mle(censor_data) )
#' @export
#' @importFrom survival Surv survreg
add_Weibull_mle <- function(censor_data) {

  n_minus_r <- censor_data[[4]] - censor_data[[1]]



  sobj <- Surv(time = c( censor_data[[3]],

                         rep(censor_data[[2]], n_minus_r)

  ),

  event = c( rep(1, censor_data[[1]] ),

             rep(0, n_minus_r)

  ),

  type = 'right'

  )

  sfit <- survreg(sobj~1, dist = 'weibull')



  return( list( Shape = 1/sfit$scale,

                Scale = as.numeric( exp(sfit$coefficients) ),

                Number_of_Failures = censor_data[[1]],

                Censor_Time = censor_data[[2]],

                Failure_Times = censor_data[[3]],

                Total_Number = censor_data[[4]])

  )

}

#' A Wrapper Function of \code{simulate_weibull_data}
#'
#' This function only gives the MLEs, discarding the censoring data.
#'
#' @param censor_data The value of function\code{simulate_weibull_data()}.
#' @export
find_Weibull_MLEs <- function(censor_data){
  weibull_list <- add_Weibull_mle(censor_data)
  MLEs <- list(shape = weibull_list$Shape,
               scale = weibull_list$Scale)
  return(MLEs)
}

#' The Log-likelihood Function
#' @keywords internal
wbfuncFRWB <- function(para, r, n, weights, times){
  a = para[1]
  b = para[2]

  beta = exp(-a)
  eta = exp(b-exp(a)*log(-log(1-r/2/n)))

  value = sum(weights[1:r]*dweibull(times[1:r], beta, eta,log = 1))+
    sum(weights[(r+1):n]*pweibull(times[(r+1):n], beta, eta, lower.tail = 0, log.p = 1))

  return(-value)
}

#' A Robust Algorithm to Compute The MLEs of the Weighted Weibull Censored Data
#'
#' Compute the MLEs of censored Weibull data by reparameterization and using \code{optim()} function.
#' @param censor_data The value of function \code{simulate_weibull_data()}.
#' @param weight Specify the way to put a weight on the data. "no": do not put the weight on the data.
#' @examples
#' censor_data <- simulate_weibull_data(1, 5, 0.2, 1, 1)
#' ( censor_data_with_mle <- add_Weibull_mle_robust(censor_data) )
#' @export
add_Weibull_mle_robust <- function(censor_data, weight = "no")
{
  options(warn=-1)

  r <- censor_data$Number_of_Failures
  n <- censor_data$Total_Number

  times <- c(censor_data$Failure_Times,
            rep(censor_data$Censor_Time, times = n-r)
            )

  if(weight == "no"){
    weights <- rep(1, n)
  } else {
    weights <- rep(1, n) # I will add the weighted mle later
  }

  init_sigma <- (log(times[r+1])-log(times[1]))/(log(-log(1-r/n))-log(-log(1-0.5/n)))
  init_a <- log(init_sigma)
  init_eta <- exp(log(times[r+1])-init_sigma*log(-log(1-r/n)))
  init_b <- log(init_eta)+init_sigma*log(-log(1-r/2/n))

  op <- optim(c(init_a, init_b), wbfuncFRWB, r = r, n = n, weights = weights, times = times)
  a <- op$par[1]
  b <- op$par[2]

  beta = exp(-a)
  eta = exp(b-exp(a)*log(-log(1-r/2/n)))

  if(op$convergence != 0){
    print("not coverge")
  }

  return(
    c(
      list(shape = beta,
           scale = eta,
           censor_data)
    )
  )
}

#' A Wrapper Function of \code{add_Weibull_mle_robust}
#'
#' This function only gives the MLEs, discarding the censoring data.
#'
#' @param censor_data The value of function\code{add_Weibull_mle_robust()}.
#' @export
find_Weibull_MLEs_robust <- function(censor_data){
  weibull_list <- add_Weibull_mle_robust(censor_data)
  MLEs <- list(shape = weibull_list$Shape,
               scale = weibull_list$Scale)
  return(MLEs)
}
