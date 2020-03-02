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
#' @param LambdaR is the factor loading of the indicator onto the factor of
#' interest for the reference group.
#' @param LambdaF is the factor loading of the indicator onto the factor of
#' interest for the focal group.
#' @param NuR is the indicator intercept for the reference group.
#' @param NuF is the indicator intercept for the focal group.
#' @param MeanF is the factor mean in the focal group
#' @param VarF is the factor variances in the focal group.
#' @param SD is the indicator standard deviations to be used as
#' the denominator of the dmacs effect size. This will usually either be
#' pooled standard deviation for the indicator or the standard deviation
#' for the indicator in the reference group.
#' @param ThreshR is a vector of thresholds (for categorical indicators)
#' for the reference group. Defaults to \code{NULL} for continuous
#' indicators.
#' @param ThreshF is a vector of thresholds (for categorical indicators)
#' for the focal group. Defaults to \code{NULL} for continuous
#' indicators.
#' @param ThetaR is the indicator residual variance in the
#' reference group. Defaults to \code{NULL} for continuous
#' indicators.
#' @param ThetaF is the indicator residual variance in the
#' focal group. Defaults to \code{NULL} for continuous
#' indicators.
#' @param categorical is a Boolean variable declaring whether the variables
#' in the model are ordered categorical. Models in which some variables are
#' categorical and others are continuous are not supported. If no value is
#' provided, categorical defaults to \code{FALSE}, although if a vector of
#' thresholds are provided, categorical will be forced to
#' \code{TRUE}. A graded response model with probit link (e.g., DWLS in
#' lavaan or WLSMV in Mplus) is used for categorical variables. If you desire
#' for other categorical models (e.g., IRT parameterization) to be supported,
#' e-mail the maintainer.
#'
#' @return The dmacs effect size of equation 3 of Nye & Drasgow (2011).
#'
#' @examples
#' LambdaF <- 0.74
#' LambdaR <- 0.76
#' NuF     <- 1.28
#' NuR     <- 0.65
#' MeanF   <- 0.21
#' VarF    <- 1.76
#' SD      <- 1.85
#' item_dmacs(LambdaR, LambdaF, NuR, NuF, MeanF, VarF, SD)
#'
#' @section References:
#' Nye, C. & Drasgow, F. (2011). Effect size indices for analyses of
#' measurement equivalence: Understanding the practical importance of
#' differences between groups. \emph{Journal of Applied Psychology, 96}(5),
#' 966-980.
#'
#' @export
#' @importFrom stats dnorm
#' @importFrom stats integrate

item_dmacs <- function (LambdaR, LambdaF,
                        NuR, NuF,
                        MeanF, VarF, SD,
                        ThreshR = NULL, ThreshF = NULL,
                        ThetaR = NULL, ThetaF = NULL,
                        categorical = FALSE) {

  # Use Thresholds as a check for categorical-ness
  if (!is.null(ThreshR)) {
    categorical <- TRUE
    ## If threshold vectors do not have the same length, throw an error
    if (length(ThreshR) != length(ThreshF)) stop("Item must have same number of thresholds in both reference and focal group")

  }

  ## If item does not load on factor, return NA
  if(LambdaR == 0) {return(NA)}

  ## Create a function for the integrand using the expected value function expected_value
  ## The sqrt(VarF) is there because we did a change of varianbles into the z metric

  integrand <- function(z, LambdaR, LambdaF,
                        NuR, NuF,
                        ThreshR, ThreshF,
                        ThetaR, ThetaF,
                        MeanF, VarF, categorical) {

    (expected_value(LambdaF, NuF, MeanF+z*sqrt(VarF), ThreshF, ThetaF, categorical) -
       expected_value(LambdaR, NuR, MeanF+z*sqrt(VarF), ThreshR, ThetaR, categorical))^2 * dnorm(z) * sqrt(VarF)

  }

  ## Now, sum it to get the integral, and compute the effect size. Stepsize is in z units, not theta units!!
  sqrt(integrate(integrand, -Inf, Inf,
                 LambdaR, LambdaF,
                 NuR, NuF,
                 ThreshR, ThreshF,
                 ThetaR, ThetaF,
                 MeanF, VarF,
                 categorical = categorical)$value)/SD

}



#' Expected bias to item mean
#'
#' \code{delta_mean_item} computes the expected bias in item mean due to
#' measurement nonequivalence.
#'
#' \code{delta_mean_item} is called by \code{dmacs_summary_single}, which
#' in turn is called by \code{\link{lavaan_dmacs}} and
#' \code{\link{mplus_dmacs}}, which are the only functions in this
#' package intended for casual users
#'
#' @param LambdaR is the factor loading of the indicator onto the factor of
#' interest for the reference group.
#' @param LambdaF is the factor loading of the indicator onto the factor of
#' interest for the focal group.
#' @param NuR is the indicator intercept for the reference group.
#' @param NuF is the indicator intercept for the focal group.
#' @param MeanF is the factor mean in the focal group
#' @param VarF is the factor variances in the focal group.
#' @param ThreshR is a vector of thresholds (for categorical indicators)
#' for the reference group. Defaults to \code{NULL} for continuous
#' indicators.
#' @param ThreshF is a vector of thresholds (for categorical indicators)
#' for the focal group. Defaults to \code{NULL} for continuous
#' indicators.
#' @param ThetaR is the indicator residual variance in the
#' reference group. Defaults to \code{NULL} for continuous
#' indicators.
#' @param ThetaF is the indicator residual variance in the
#' focal group. Defaults to \code{NULL} for continuous
#' indicators.
#' @param categorical is a Boolean variable declaring whether the variables
#' in the model are ordered categorical. Models in which some variables are
#' categorical and others are continuous are not supported. If no value is
#' provided, categorical defaults to \code{FALSE}, although if a vector of
#' thresholds are provided, categorical will be forced to
#' \code{TRUE}. A graded response model with probit link (e.g., DWLS in
#' lavaan or WLSMV in Mplus) is used for categorical variables. If you desire
#' for other categorical models (e.g., IRT parameterization) to be supported,
#' e-mail the maintainer.
#'
#' @return The expected bias in item mean due to
#' measurement nonequivalence in equation 4 of Nye & Drasgow (2011).
#'
#' @examples
#' LambdaF <- 0.74
#' LambdaR <- 0.76
#' NuF     <- 1.28
#' NuR     <- 0.65
#' MeanF   <- 0.21
#' VarF    <- 1.76
#' delta_mean_item(LambdaR, LambdaF, NuR, NuF, MeanF, VarF)
#'
#' @section References:
#' Nye, C. & Drasgow, F. (2011). Effect size indices for analyses of
#' measurement equivalence: Understanding the practical importance of
#' differences between groups. \emph{Journal of Applied Psychology, 96}(5),
#' 966-980.
#'
#' @export
#' @importFrom stats dnorm
#' @importFrom stats integrate

delta_mean_item <- function (LambdaR, LambdaF,
                             NuR, NuF,
                             MeanF, VarF,
                             ThreshR = NULL, ThreshF = NULL,
                             ThetaR = NULL, ThetaF = NULL,
                             categorical = FALSE) {
  # Use Thresholds as a check for categorical-ness
  if (categorical) {
    categorical <- TRUE
    ## If threshold vectors do not have the same length, throw an error
    if (length(ThreshR) != length(ThreshF)) stop("Item must have same number of thresholds in both reference and focal group")

  }

  ## If item does not load on factor, return NA
  if(LambdaR == 0) {return(NA)}

  ## Create a function for the integrand using the expected value function expected_value
  ## The sqrt(VarF) is there because we did a change of variables into the z metric.
  ## Probably would have been easier to fix the dnorm, but it's done.
  integrand <- function(z, LambdaR, LambdaF,
                           NuR, NuF,
                           ThreshR, ThreshF,
                           ThetaR, ThetaF,
                           MeanF, VarF, categorical) {

    (expected_value(LambdaF, NuF, MeanF+z*sqrt(VarF), ThreshF, ThetaF, categorical) -
       expected_value(LambdaR, NuR, MeanF+z*sqrt(VarF), ThreshR, ThetaR, categorical)) * dnorm(z) * sqrt(VarF)

  }
  ## Now, integrate
  integrate(integrand, -Inf, Inf,
            LambdaR, LambdaF,
            NuR, NuF,
            ThreshR, ThreshF,
            ThetaR, ThetaF,
            MeanF, VarF,
            categorical = categorical)$value

}


#' Expected bias to total score variance
#'
#' \code{delta_var} computes the expected bias in total score variance due
#' to measurement nonequivalence. \code{delta_var} will only work for
#' unidimensional linear models (not categorical).
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
#' @param categorical is a Boolean variable declaring whether the variables
#' in the model are ordered categorical. Categorical indicators are not supported
#' for this function.
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

delta_var <- function (LambdaR, LambdaF, VarF, categorical = FALSE) {

  if(categorical) {
    warning("At this time, delta variance can only be computed for linear models, not for categorical ones")
    return(NULL)
  }
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
#' @param Nu is the indicator intercept.
#' @param Eta is the value of the factor for which the expected value is
#' to be computed
#' @param Thresh is the vector of thresholds for a categorical indicator.
#' Defaults to NULL as \code{Thresh} is not needed for continuous indicators.
#' @param categorical is a Boolean variable declaring whether the variables
#' in the model are ordered categorical. Models in which some variables are
#' categorical and others are continuous are not supported. If no value is
#' provided, categorical defaults to \code{FALSE}. A graded response model
#' with probit link (e.g., DWLS in lavaan or WLSMV in Mplus) is used for
#' categorical variables. If you desire for other categorical models
#' (e.g., IRT parameterization) to be supported, e-mail the maintainer.
#'
#' @return The expected value of the indicator when the factor score is \code{Eta}
#' @keywords internal
#' @importFrom stats pnorm

expected_value <- function (Lambda, Nu, Eta, Thresh = NULL, Theta = NULL, categorical = FALSE) {

  if (categorical) {
    ## Graded Response model with probit link.

    ## Let's give ourselves a y* score
    Mu <- Nu + Lambda*Eta

    ## categories go from 0 to number of thresholds
    max <- length(Thresh)
    ## Make a max+1 category that is impossible to attain
    Thresh[max+1] <- Inf

    ## Start expected value at 0 and incremented with expectation value from each category
    expected <- 0
    for (i in 1:max) {
      expected <- expected + i*(pnorm(Mu - Thresh[i], sd = sqrt(Theta)) -
                                                  pnorm(Mu - Thresh[i+1], sd = sqrt(Theta)))
    }
    expected

  } else {
    ## The linear continuous world is so easy!!
    Nu+Lambda*Eta
  }

}








#' Expected bias in test score
#'
#' \code{DTF_graph} produces a graph of expected total score or of expected
#' bias in total score based on item parameters in the reference group and
#' based on item parameters in the focal group.
#'
#' \code{DTF_graph} is called by \code{dmacs_summary}, which
#' in turn is called by \code{\link{lavaan_dmacs}} and
#' \code{\link{mplus_dmacs}}, which are the only functions in this
#' package intended for casual users
#'
#' @param LambdaR is the factor loading of the indicator onto the factor of
#' interest for the reference group.
#' @param LambdaF is the factor loading of the indicator onto the factor of
#' interest for the focal group.
#' @param NuR is the indicator intercept for the reference group.
#' @param NuF is the indicator intercept for the focal group.
#' @param MeanF is the factor mean in the focal group
#' @param VarF is the factor variances in the focal group.
#' @param ThreshR is a vector of thresholds (for categorical indicators)
#' for the reference group. Defaults to \code{NULL} for continuous
#' indicators.
#' @param ThreshF is a vector of thresholds (for categorical indicators)
#' for the focal group. Defaults to \code{NULL} for continuous
#' indicators.
#' @param ThetaR is the indicator residual variance in the
#' reference group. Defaults to \code{NULL} for continuous
#' indicators.
#' @param ThetaF is the indicator residual variance in the
#' focal group. Defaults to \code{NULL} for continuous
#' indicators.
#' @param categorical is a Boolean variable declaring whether the variables
#' in the model are ordered categorical. Models in which some variables are
#' categorical and others are continuous are not supported. If no value is
#' provided, categorical defaults to \code{FALSE}, although if a vector of
#' thresholds are provided, categorical will be forced to
#' \code{TRUE}. A graded response model with probit link (e.g., DWLS in
#' lavaan or WLSMV in Mplus) is used for categorical variables. If you desire
#' for other categorical models (e.g., IRT parameterization) to be supported,
#' e-mail the maintainer.
#' @param minEta is the minimum factor score to be considered. Defaults to -3.
#' @param maxEta is the maximum factor score to be considered. Defaults to -3.
#'
#' @return A graph of expected total score or bias in expected total score due to
#' measurement nonequivalence, per Nye & Drasgow (2011).
#'
#' @examples
#' LambdaF <- matrix(c(1.00, 0.74,  1.14, 0.92), ncol = 1)
#' LambdaR <- matrix(c(1.00, 0.76,  1.31, 0.98), ncol = 1)
#' NuF     <- c(0.00, 1.28, -0.82, 0.44)
#' NuR     <- c(0.00, 0.65, -0.77, 0.47)
#' DTF_graph(LambdaR, LambdaF, NuR, NuF)
#'
#' @section References:
#' Nye, C. & Drasgow, F. (2011). Effect size indices for analyses of
#' measurement equivalence: Understanding the practical importance of
#' differences between groups. \emph{Journal of Applied Psychology, 96}(5),
#' 966-980.
#'
#' @export
#' @importFrom stats reshape

DTF_graph       <- function (LambdaR, LambdaF,
                             NuR, NuF,
                             RefGroup, FocGroup,
                             ThreshR = NULL, ThreshF = NULL,
                             ThetaR = NULL, ThetaF = NULL,
                             categorical = FALSE,
                             minEta = -3, maxEta = 3,
                             GraphType = "Scores") {

  ## Use Thresholds as a check for categorical-ness
  if (categorical) {
    categorical <- TRUE
    ## If threshold vectors do not have the same length, throw an error
    for (i in 1:length(ThreshR)) {
      if (length(ThreshR[[i]]) != length(ThreshF[[i]])) stop("Item must have same number of thresholds in both reference and focal group")
    }
  }

  ## Vector of Eta values to compute expected scores or biases
  Eta <- seq(minEta, maxEta, length.out = 100)

  ## iterate over the items to create expected scores
  # reference group
  ref_expected <- sapply(1:length(LambdaF), function (x) {
    mapply(expected_value, Eta = Eta,  MoreArgs = list(Lambda = LambdaR[x], Nu = NuR[x],
                                                       Thresh = ThreshR[[x]], Theta = ThetaR[[x]],
                                                       categorical = categorical))
  })
  # total score
  ref_expected_total <- rowSums(ref_expected)


  # focal group
  foc_expected <- sapply(1:length(LambdaF), function (x) {
    mapply(expected_value, Eta = Eta,  MoreArgs = list(Lambda = LambdaF[x], Nu = NuF[x],
                                                       Thresh = ThreshF[[x]], Theta = ThetaF[[x]],
                                                       categorical = categorical))
  })
  # total score
  foc_expected_total <- rowSums(foc_expected)

  #bias
  bias_expected <- foc_expected_total - ref_expected_total

  # merge into a dataframe for graphing with
  graph_data <- data.frame(list(Eta = Eta, Ref = ref_expected_total, Foc = foc_expected_total, Bias = bias_expected))


  if (GraphType == "Scores") {

    # make it long so that we can confrom to Hadley's ridiculous standards
    graph_data_long <- reshape(graph_data[,1:3], varying = list(2:3), direction = "long")
    colnames(graph_data_long) <- c("Eta", "Group", "Expected", "id")
    graph_data_long$Group[graph_data_long$Group == 1] <- RefGroup
    graph_data_long$Group[graph_data_long$Group == 2] <- FocGroup

    ## Plot both sets of scores
    ggplot2::ggplot(graph_data_long, aes(x = Eta)) +
                    geom_line(aes(y = Expected, group = Group, color = Group)) +
                    scale_colour_manual(values=c("red","blue")) +
                    ggtitle("Expected Total Scores") +
                    xlab("Factor Score") +
                    ylab("Expected Total Score") +
                    theme(plot.title = element_text(hjust = 0.5), legend.position =  c(0.9, 0.2)) +
                    scale_x_continuous(limits=c(minEta, maxEta), breaks=seq(minEta,maxEta,1))

  } else {
    # Plot the bias
    ggplot2::ggplot(graph_data, aes(x = Eta)) +
                    geom_line(aes(y = Bias), colour = "red") +
                    ggtitle(paste0("Expected Bias between ", FocGroup, " and ", RefGroup)) +
                    xlab("Factor Score") +
                    ylab("Expected Bias") +
                    theme(plot.title = element_text(hjust = 0.5)) +
                    scale_x_continuous(limits=c(minEta, maxEta), breaks=seq(minEta,maxEta,1))
  }


  expected_value <- function (Lambda, Nu, Eta, Thresh = NULL, Theta = NULL, categorical = FALSE) {}
}

