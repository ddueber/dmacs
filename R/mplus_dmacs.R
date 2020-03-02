
#' Summary of measurement nonequivalence effects
#'
#' \code{mplus_dmacs} returns a summary of measurement non-equivalence
#' effects given an Mplus .out file.
#'
#' @param fit is an Mplus .out file of a multigroup CFA analysis. The default
#' is to launch a window for choosing the file.
#' @param RefGroup can be the name of the reference group (as a string),
#' or the index of the reference group (as a number). RefGroup defaults to
#' the first group if no value is provided. It is strongly recommended to
#' provide the reference group as a string, since group names in data are
#' often ordered by their appearance in the data, not alphabetically.
#' @param dtype described the pooling of standard deviations for use in the
#' denominator of the dmacs effect size. Possibilities are "pooled" for
#' pooled standard deviations, or "glass" for always using the standard
#' deviation of the reference group.
#'
#' @return A list, indexed by group, of lists of measurement nonequivalence
#' effects from Nye and Drasgow (2011), including dmacs, expected bias in
#' the mean score by item,
#' expected bias in the mean total score, and expected bias in the variance
#' of the total score. Expected bias in the variance of the total score is
#' only supplied for unidimensional models in the current version of this
#' package
#'
#'
#' @section References:
#' Nye, C. & Drasgow, F. (2011). Effect size indices for analyses of
#' measurement equivalence: Understanding the practical importance of
#' differences between groups. \emph{Journal of Applied Psychology, 96}(5),
#' 966-980.
#' @export

mplus_dmacs <- function(fit = file.choose(),  RefGroup = 1, dtype = "pooled", MEtype = "Group") {
  ## We will need the raw text file later if continuous variables with pooled dtype
  fit0 <- fit

  ## Read in the model with MplusAutomation
  if (!any(class(fit) == "mplus.model")) {
    fit <- MplusAutomation::readModels(fit)
  }

  ## Things proceed completely differently for longitudinal vs group invariance
  if (grepl("ong", MEtype, fixed = TRUE)) {

    ## Because it will save a lot of space, and just make life easier, let's give fit$parameters$unstandardized a name
    Params <- fit$parameters$unstandardized

    ## Group names are now factor names
    Groups <- unique(Params[grepl(".BY",  Params$paramHeader), "paramHeader"])
    Groups <- gsub('.{3}$', '', Groups)
    names(Groups) <- Groups

    ## If RefGroup is a string, turn it into an index, else issue a warning that maybe it should be a string
    if (is.character(RefGroup)) {
      RefGroup <- match(RefGroup, Groups)
    } else {
      warning(paste("It is recommended that you provide the name of the factor for the reference time to RefGroup = as a string; see ?mplus_dmacs. The reference group being used is", Groups[RefGroup]))
    }

    ## Now, let's grab the loading vectors for each timepoint
    ## Warning! No item names are attached to these loadings
    LambdaList <- lapply (Groups, function (G) {
      as.matrix(Params[(Params$paramHeader == paste0(G, ".BY")), "est"])
    })

    ## Get an array of item names for the reference time. We will need this later to attach names to items
    ItemNames <- unique(Params[grepl(paste0(Groups[RefGroup],".BY"),  Params$paramHeader), "param"])
    names(ItemNames) <- ItemNames


    ## Make a list of factor means
    MeanList <- lapply(Groups, function (G) {
      Params[(Params$paramHeader == "Means") & (Params$param == G), "est"]
    })

    ## Make a list of factor variances
    VarList <- lapply(Groups, function (G) {
      Params[(Params$paramHeader == "Variances") & (Params$param == G), "est"]
    })

    ## Make a list of item intercepts (and thresholds and residual variances)
    ## THIS IS A HACK - IT ONLY CONSIDERS THE CONTINUOUS CASE. PLEASE FIX LATER TO ALSO CONSIDER THE CATEGORICAL CASE!!
    categorical <- FALSE
    ## Make a list of items names by timepoint
    GroupItems <- lapply(Groups, function (G) {
      ItemNames <- unique(Params[grepl(paste0(G,".BY"),  Params$paramHeader), "param"])
    })
    ## now fetch NuList
    NuList <- lapply(Groups, function (G) {
      sapply(GroupItems[[G]], function (I) {
        ## Intercept value for Item I, Group G
        Params[(Params$paramHeader == "Intercepts") & (Params$param == I), "est"]
      })
    })

    ## ThreshList and ThetaList do not exist here
    ThreshList <- NULL
    ThetaList  <- NULL

    ## Now fetch item SDs
    ## THIS IS A HACK - IT ONLY CONSIDERS THE CONTINUOUS CASE. PLEASE FIX LATER TO ALSO CONSIDER THE CATEGORICAL CASE!!

    ## Should probably throw an error if sampstat is not included
    AllSDs <- sqrt(diag(fit$sampstat$covariances))
    RawSDList <- lapply(Groups, function (G) {
      AllSDs[GroupItems[[G]]]
    })

    ## Do "pooled" and "glass" separately
    if (dtype == "pooled") {
      # Since sample sizes are the same, pooled just means "average"
      SDList <- lapply(Groups, function (G) {
        (RawSDList[[G]] + RawSDList[[RefGroup]]) / 2
      })
    } else if (dtype == "glass") {
      ItemSDs <- RawSDList[[RefGroup]]
      SDList <- lapply(Groups, function (x) {ItemSDs})
    } else {
      stop("Only \"pooled\" and \"glass\" SD types are supported")
    }


  } else { ## Group Invariance testing
    ## Get an array of group names.
    Groups <- unique(fit$parameters$unstandardized$Group)
    names(Groups) <- Groups

    ## If RefGroup is a string, turn it into an index, else issue a warning that maybe it should be a string
    if (is.character(RefGroup)) {
      RefGroup <- match(RefGroup, Groups)
    } else {
      warning(paste("It is recommended that you provide the name of the reference group as a string; see ?mplus_dmacs. The reference group being used is", Groups[RefGroup]))
    }

    ## Because it will save a lot of space, and just make life easier, let's give fit$parameters$unstandardized a name
    Params <- fit$parameters$unstandardized

    ## Get an array of factor names
    FactorNames = unique(Params[which(Params$paramHeader == "Means"),"param"])
    names(FactorNames) <- FactorNames

    ## Get an array of item names
    ItemNames <- unique(Params[grepl(".BY",  Params$paramHeader), "param"])
    names(ItemNames) <- ItemNames

    ## Make a list of loading matrices
    LambdaList <- lapply(Groups, function (G) {
      sapply(FactorNames, function (F) {
        sapply(ItemNames, function (I) {
          ## loading value for factor F, item I, Group G
          lambda <- Params[(Params$paramHeader == paste0(F, ".BY")) & (Params$param == I) & (Params$Group == G), "est"]
          ## If lambda is empty, return zero
          if (length(lambda) == 1) {
            lambda
          } else {
            0
          }
        })
      })
    })

    ## Make a list of factor means
    MeanList <- lapply(Groups, function (G) {
      sapply(FactorNames, function (F) {
        ## mean value for factor F, group G
        Params[(Params$paramHeader == "Means") & (Params$param == F) & (Params$Group == G), "est"]
      })
    })

    ## Make a list of factor variances
    VarList <- lapply(Groups, function (G) {
      sapply(FactorNames, function (F) {
        ## variance value for factor F, group G
        Params[(Params$paramHeader == "Variances") & (Params$param == F) & (Params$Group == G), "est"]
      })
    })

    ## Make a list of Intercepts/Thresholds
    ## Check to see if we are in a continuous or categorical context, because ThreshList works very differently in those two contexts
    if (length(fit$input$variable$categorical) == 0) {
      categorical <- FALSE
      ## now fetch NuList
      NuList <- lapply(Groups, function (G) {
        sapply(ItemNames, function (I) {
          ## Intercept value for Item I, Group G
          Params[(Params$paramHeader == "Intercepts") & (Params$param == I) & (Params$Group == G), "est"]
        })
      })

      ## ThreshList and ThetaList do not exist here
      ThreshList <- NULL
      ThetaList  <- NULL

    } else {
      categorical <- TRUE

      ## If "Scales" is in paramHeader, throw an error and demand THETA parameterization
      if ("Scales" %in% Params$paramHeader) stop(cat("Only THETA parameterization is supported for models fit by Mplus. Please insert \n Analysis: parameterization = THETA;"))

      ## Create NuList which is all zeros.
      NuList <- lapply(Groups, function (G) {
        sapply(ItemNames, function (I) {
          ## Intercept value for Item I, Group G
          0
        })
      })

      ## now fetch ThreshList
      ThreshList <- lapply(Groups, function (G) {
        lapply(ItemNames, function(I) {
          ## All threshold values (in a vector... they are in order in "Params") of item I, group G
          Params[(Params$paramHeader == "Thresholds") & (grepl(paste0("^", I, "\\$"), Params$param)) & (Params$Group == G), "est"]
        })
      })

      ## Now fetch ThetaList
      ThetaList <- lapply(Groups, function (G) {
        sapply(ItemNames, function (I) {
          # All residual variances in a vector
          Params[(Params$paramHeader == "Residual.Variances") & (Params$param == I) & (Params$Group == G), "est"]
        })
      })

    }

    ## All that's left is to compute the SDs... which doesn't sound fun at all
    if (categorical) {
      if (dtype == "pooled") {
        ## Make a list of (item SDs and item sample sizes) indexed by the groups
        ## We need the sample sizes so that we can pool
        ## This is faster, but harder to follow than making separate SD and SS lists
        RawSDList <- lapply (Groups, function (G) {
          ## MplusAutomation does a terrible job at organizing fit$sampstat - the groups are given the wrong names!!
          ## So we are going to take the Group name and turn it into an index
          GIndex <- which(Groups == G)
          ## This variable is just for readability
          GroupCounts <- fit$sampstat[[GIndex]]$proportions.counts
          ## Make a vector of item SDs
          sapply(ItemNames, function (I) {
            ## Get the counts for each category, make a vector with the right number in each category, then get the SD
            ItemCounts <- GroupCounts[grepl(I, GroupCounts$variable),]
            ItemData <- NULL
            for (i in 1:length(ItemCounts)) {
              ItemData <- c(ItemData, rep(ItemCounts[i, "category"], times = ItemCounts[i, "count"]))
            }
            c(sd(ItemData), length(ItemData))
          })
        })
        SDList <- lapply (Groups, function (G) {
          ## everything is vectorized here, so we can just have at it, I hope!!
          (RawSDList[[G]][1,]*(RawSDList[[G]][2,] - 1) + RawSDList[[RefGroup]][1,]*(RawSDList[[RefGroup]][2,] - 1)) /
            (RawSDList[[G]][2,] - 1 + RawSDList[[RefGroup]][2,] - 1)
        })

      } else if (dtype == "glass") {
        ## MplusAutomation organizes fit$sampstat strangely, so we need the index here
        GroupCounts <- fit$sampstat[[RefGroup]]$proportions.counts
        ## Make a vector of item SDs
        ItemSDs <- sapply(ItemNames, function (I) {
          ## Get the counts for each category, make a vector with the right number in each category, then get the SD
          ItemCounts <- GroupCounts[grepl(I, GroupCounts$variable),]
          ItemData <- NULL
          for (i in 1:length(ItemCounts)) {
            ItemData <- c(ItemData, rep(ItemCounts[i, "category"], times = ItemCounts[i, "count"]))
          }
          sd(ItemData)
        })
        ## Now, lapply that vector for each group
        SDList <- lapply(Groups, function (x) {ItemSDs})
      } else {
        stop("Only \"pooled\" and \"glass\" SD types are supported")
      }

    } else { ## We are in continuous case, and should extract SDs from sampstat if we can... estimate from model (with warning) if no sampstat

      if (dtype == "pooled") {
        ## Should probably throw an error if sampstat is not included
        RawSDList <- lapply(Groups, function (G) {
          ## MplusAutomation does a terrible job at organizing fit$sampstat - the groups are given the wrong names!!
          ## So we are going to take the Group name and turn it into an index
          GIndex <- which(Groups == G)
          sqrt(diag(fit$sampstat[[GIndex]]$covariances))
        })
        ## Now we need to fetch the group sample sizes... which MplusAutomation does not do for us!!!
        FitLines <- readLines(fit0)
        GroupCounts <- lapply (Groups, function (G) {
          ## The first time something like "Group GROUP1" shows us is the sample size line for that group
          GroupCountLine <- FitLines[grepl(paste("Group", G), FitLines)][1]
          JustGroupCount <- gsub(paste("Group", G), "", GroupCountLine)
          ## The only thing left in the string is spaces and the group sample size
          as.numeric(JustGroupCount)
        })

        ## Now, let's pool the results
        SDList <- lapply(Groups, function (G) {
          (RawSDList[[G]]*(GroupCounts[[G]]-1) + RawSDList[[RefGroup]]*(GroupCounts[[RefGroup]]-1)) /
            ((GroupCounts[[G]]-1) + (GroupCounts[[RefGroup]]-1))
        })
      } else if (dtype == "glass") {
        ItemSDs <- sqrt(diag(fit$sampstat[[RefGroup]]$covariances))
        SDList <- lapply(Groups, function (x) {ItemSDs})
      } else {
        stop("Only \"pooled\" and \"glass\" SD types are supported")
      }
    }
  }


  Results <- dmacs_summary(LambdaList, NuList,
                           MeanList, VarList, SDList,
                           Groups, RefGroup,
                           ThreshList, ThetaList,
                           categorical)

  ## Note to self - we may need to insert some names here!!

  Results

}
