#' Simulate Weibull censored data.
#'
#' @param type Censoring Type
#' @param r Expected Number of Failures
#' @param Pt Fraction of Failing
#' @param shape The Weibull shape parameter
#' @param scale The Weibull scale parameter
#' @examples
#' simulate_weibull_data(2, 10, 0.1, 1, 1)
#' @export
#' @importFrom stats dweibull optim pweibull runif
simulate_weibull_data <- function(type, r , Pt, shape, scale) {

  n <- round( r/Pt )

  if (type == 1){
    y_vector <- numeric(n)
    x <- 0

    for ( i in 0:(n-1) ){
      u <- runif(1)
      y <- 1 - (1-x) * ( 1-u )^ ( 1.0 / (n-i) )

      if (y > Pt) {break}

      y_vector[i+1] <- y
      x <- y
      }

    if( sum( y_vector > 0 ) < 2){

      output <- simulate_weibull_data(type, r, Pt, shape, scale)

      return(output)
    }

    num_failure <- sum( y_vector > 0 )

    failure_times <- scale * ( -log( 1 - y_vector[ 1:num_failure ] ) ) ^ ( 1/shape )

    censor_time <- scale * ( -log( 1 - Pt ) ) ^ ( 1/shape )

    output <- list(Number_of_Failures = num_failure, Censor_Time = censor_time,
                   Failure_Times = failure_times, Total_Number = n)

    return(output)
  }

  else if(type == 2)
  {
    y_vector <- numeric(r)
    x <- 0

    for(i in 0 : (r-1) )
    {
      u <- runif(1)

      y <- 1 - (1-x) * ( 1-u )^ ( 1.0 / (n-i) )

      y_vector[ i+1 ] <- y

      x <- y
    }

    failure_times <- scale * ( -log( 1 - y_vector ) ) ^ ( 1/shape )

    return(list(Number_of_Failures = r, Censor_Time = failure_times[r],
                Failure_Times = failure_times, Total_Number = n))
  }
}
