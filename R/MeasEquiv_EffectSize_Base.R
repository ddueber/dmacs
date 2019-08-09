#' Standard deviations of columns
#'
#' \code{colSD} computes standard deviations of columns.
#'
#' @param x is a matrix or data frame for which we want to obtain column sds
#' \@param ... are other arguments to be passed to \code{sd}, such as
#' \code{na.rm}
#'
#' @return A vector of standard deviations by column
#' @keywords internal
#' @importFrom stats sd
#'
colSD <- function(x, ...) {apply(X=x, MARGIN=2, FUN=sd, ...)}


#' dmacs measurement nonequivalence effect size
#'
#' \code{item_dmacs} computes the dmacs effect size for a single indicator
#' relative to a single factor in a single focal group
#'
#' \code{item_dmacs} is called by \code{dmacs_summary_single}, which
#' in turn is called by \code{\link{lavaan_dmacs}} and
#' \code{\link{mplus_dmacs}}, which are the only functions in this
#' package intended for casual users
#'
#' @param LambdaR is the factor loading of the item onto the factor of
#' interest for the reference group.
#' @param ThreshR is the indicator intercept (for continuous
#' indicators) or a vector of thresholds (for
#' categorical indicators) for the reference group.
#' @param LambdaF is the factor loading of the item onto the factor of
#' interest for the focal group.
#' @param ThreshF is the indicator intercept (for continuous
#' indicators) or a vector of thresholds (for
#' categorical indicators) for the focal group.
#' @param MeanF is the factor mean in the focal group
#' @param VarF is the factor variances in the focal group.
#' @param SD is the indicator standard deviations to be used as
#' the denominator of the dmacs effect size. This will usually either be
#' pooled standard deviation for the indicator or the standard deviation
#' for the indicator in the reference group.
#' @param categorical is a Boolean variable declaring whether the variables
#' in the model are ordered categorical. Models in which some variables are
#' categorical and others are continuous are not supported. If no value is
#' provided, categorical defaults to \code{FALSE}, although if a vector of
#' thresholds are provided, categorical will be forced to
#' \code{TRUE}. A graded response model with probit link (e.g., DWLS in
#' lavaan or WLSMV in Mplus) is used for categorical variables. If you desire
#' for other categorical models (e.g., IRT parameterization) to be supported,
#' e-mail the maintainer.
#' @param stepsize is the interval width for the Riemann sum used to estimate
#' the integral in equation 3 of Nye & Drasgow (2011). Default value is .001.
#' A larger value can be used for faster performance; accuracy is
#' excellent at \code{stepsize = .01} in my simulations.
#'
#' @return The dmacs effect size of equation 3 of Nye & Drasgow (2011).
#'
#' @examples
#' LambdaF <- 0.74
#' LambdaR <- 0.76
#' ThreshF <- 1.28
#' ThreshR <- 0.65
#' MeanF   <- 0.21
#' VarF    <- 1.76
#' SD      <- 1.85
#' item_dmacs(LambdaR, ThreshR, LambdaF, ThreshF, MeanF, VarF, SD)
#'
#' @section References:
#' Nye, C. & Drasgow, F. (2011). Effect size indices for analyses of
#' measurement equivalence: Understanding the practical importance of
#' differences between groups. \emph{Journal of Applied Psychology, 96}(5),
#' 966-980.
#'
#' @export
#' @importFrom stats dnorm

item_dmacs <- function (LambdaR, ThreshR,
                        LambdaF, ThreshF,
                        MeanF, VarF,
                        SD, categorical = FALSE,
                        stepsize = .001) {

  ## If item does not load on factor, return NA
  if(LambdaR == 0) {return(NA)}

  ## if more than one threshold, we must be in a categorical situation
  if (length(ThreshR) > 1) { categorical <- TRUE}

  ## integrate over z from -5 to 5 by .001, which gives us 10,000 integration points, which is enough!
  z <- seq(-5, 5, stepsize)
  ## Compute the integrand using the expected value function expected_value
  integrand <- (expected_value(LambdaF, ThreshF, MeanF+z*sqrt(VarF), categorical) -
                  expected_value(LambdaR, ThreshR, MeanF+z*sqrt(VarF), categorical))^2 * dnorm(z)
  ## Now, sum it to get the integral, and compute the effect size. Stepsize is in z units, not theta units!!
  sqrt(sum(integrand*stepsize*sqrt(VarF)))/SD

}



#' Expecte bias to item mean
#'
#' \code{delta_mean_item} computes the expected bias in item mean due to
#' measurement nonequivalence.
#'
#' \code{delta_mean_item} is called by \code{dmacs_summary_single}, which
#' in turn is called by \code{\link{lavaan_dmacs}} and
#' \code{\link{mplus_dmacs}}, which are the only functions in this
#' package intended for casual users
#'
#' @param LambdaR is the factor loading of the item onto the factor of
#' interest for the reference group.
#' @param ThreshR is the indicator intercept (for continuous
#' indicators) or a vector of thresholds (for
#' categorical indicators) for the reference group.
#' @param LambdaF is the factor loading of the item onto the factor of
#' interest for the focal group.
#' @param ThreshF is the indicator intercept (for continuous
#' indicators) or a vector of thresholds (for
#' categorical indicators) for the focal group.
#' @param MeanF is the factor mean in the focal group
#' @param VarF is the factor variances in the focal group.
#' @param categorical is a Boolean variable declaring whether the variables
#' in the model are ordered categorical. Models in which some variables are
#' categorical and others are continuous are not supported. If no value is
#' provided, categorical defaults to \code{FALSE}, although if a vector of
#' thresholds are provided, categorical will be forced to
#' \code{TRUE}. A graded response model with probit link (e.g., DWLS in
#' lavaan or WLSMV in Mplus) is used for categorical variables. If you desire
#' for other categorical models (e.g., IRT parameterization) to be supported,
#' e-mail the maintainer.
#' @param stepsize is the interval width for the Riemann sum used to estimate
#' the integral in equation 6 of Nye & Drasgow (2011). Default value is .001.
#' A larger value can be used for faster performance; accuracy is
#' excellent at \code{stepsize = .01} in my simulations.
#'
#' @return The expected bias in item mean due to
#' measurement nonequivalence in equation 4 of Nye & Drasgow (2011).
#'
#' @examples
#' LambdaF <- 0.74
#' LambdaR <- 0.76
#' ThreshF <- 1.28
#' ThreshR <- 0.65
#' MeanF   <- 0.21
#' VarF    <- 1.76
#' SD      <- 1.85
#' delta_mean_item(LambdaR, ThreshR, LambdaF, ThreshF, MeanF, VarF, SD)
#'
#' @section References:
#' Nye, C. & Drasgow, F. (2011). Effect size indices for analyses of
#' measurement equivalence: Understanding the practical importance of
#' differences between groups. \emph{Journal of Applied Psychology, 96}(5),
#' 966-980.
#'
#' @export
#' @importFrom stats dnorm

delta_mean_item <- function (LambdaR, ThreshR,
                             LambdaF, ThreshF,
                             MeanF, VarF,
                             categorical = FALSE, stepsize = .001) {

  ## If item does not load on factor, return NA
  if(LambdaR == 0) {return(NA)}

  ## if more than one threshold, we must be in a categorical situation
  if (length(ThreshR) > 1) { categorical <- TRUE}

  ## integrate over z from -5 to 5 by .001, which gives us 10,000 integration points, which is overkill
  z <- seq(-5, 5, stepsize)
  ## Compute the integrand using the expected value function expected_value
  integrand <- (expected_value(LambdaF, ThreshF, MeanF+z*sqrt(VarF), categorical) -
                  expected_value(LambdaR, ThreshR, MeanF+z*sqrt(VarF), categorical)) * dnorm(z)
  ## Now, sum it to get the integral. Stepsize is in z units, not theta units!!
  sum(integrand*stepsize*sqrt(VarF))

}


#' Expecte bias to total score variance
#'
#' \code{delta_var} computes the expected bias in total score variance due
#' to measurement nonequivalence. \code{delta_var} will only work for
#' unidimensional models.
#'
#' \code{delta_var} is called by \code{dmacs_summary_single}, which
#' in turn is called by \code{\link{lavaan_dmacs}} and
#' \code{\link{mplus_dmacs}}, which are the only functions in this
#' package intended for casual users
#'
#' @param LambdaR is the vector of factor loadings for the
#' reference group.
#' @param LambdaF is the vector of factor loadings for the
#' focal group.
#' @param VarF is the factor variance of the focal group.
#'
#' @return The expected bias in total score variance due to
#' measurement nonequivalence in equation 7, 8, and 9 of Nye & Drasgow (2011).
#'
#' @examples
#' LambdaF <- c(1.00, 0.74,  1.14, 0.92)
#' LambdaR <- c(1.00, 0.76,  1.31, 0.98)
#' VarF    <- 1.76
#' delta_var(LambdaR, LambdaF, VarF)
#'
#' @section References:
#' Nye, C. & Drasgow, F. (2011). Effect size indices for analyses of
#' measurement equivalence: Understanding the practical importance of
#' differences between groups. \emph{Journal of Applied Psychology, 96}(5),
#' 966-980.
#'
#' @export

delta_var <- function (LambdaR, LambdaF, VarF) {

  delta_cov_mat <- matrix(nrow=length(LambdaR), ncol=length(LambdaR))
  ## I know for loops are supposed to be bad, but this is SO CLEAN!
  for (i in 1:length(LambdaR)) {
    for (j in 1:length(LambdaR)) {
      delta_cov_mat[i,j] <- LambdaR[[j]]*(LambdaF[[i]]-LambdaR[[i]])*VarF
                          + LambdaR[[i]]*(LambdaF[[j]]-LambdaR[[j]])*VarF
      + (LambdaF[[i]]-LambdaR[[i]])*(LambdaF[[j]]-LambdaR[[j]])*VarF
    }
  }
  sum(delta_cov_mat)

}

#' Expected value of an indicator
#'
#' \code{expected_value} returns the expected value of an indicator given
#' item parameters and factor value.
#'
#' @param Lambda is the loading of the indicator on the factor.
#' @param Thresh is the intercept (continuous indicator) or vector of
#' thresholds (categorical indicator) for the indicator.
#' @param Theta is the value of the factor for which the expected value is
#' to be computed
#' \@param categorical is a Boolean variable declaring whether the variables
#' in the model are ordered categorical. Models in which some variables are
#' categorical and others are continuous are not supported. If no value is
#' provided, categorical defaults to \code{FALSE}, although if multiple
#' thresholds are provided for an item, categorical will be forced to
#' \code{TRUE}. A graded response model with probit link (e.g., DWLS in
#' lavaan or WLSMV in Mplus) is used for categorical variables. If you desire
#' for other categorical models (e.g., IRT parameterization) to be supported,
#' e-mail the maintainer.
#'
#' @return The expected value of the indicator when the factor score is \code{Theta}
#' @keywords internal
#' @importFrom stats pnorm

expected_value <- function (Lambda, Thresh, Theta, categorical = FALSE) {
  ## if more than one threshold, we must be in a categorical situation -- this line still needs to be tested with categorical variables
  if (length(Thresh[[1]]) > 1) { categorical <- TRUE }

  if (categorical) {
    ## Graded Response model with probit link.
    ## categories go from 0 to number of thresholds
    max <- length(Thresh)
    ## let's be lazy and make a max+1 category that is impossible to attain
    Thresh[max+1] <- 9999999

    ## I know for loops are bad, but using both current and next element of vector in sapply is hard for me!
    expected <- 0
    for (i in 1:max) {expected <- expected + i*(pnorm(Lambda*(Theta-Thresh[i]))-pnorm(Lambda*(Theta-Thresh[i+1])))}
    expected

  } else {
    ## The linear continuous world is so easy!!
    Thresh+Lambda*Theta
  }

}

