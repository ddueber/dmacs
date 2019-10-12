# dmacs
dmacs R package for computing measurement nonequivalence indices

## Overview

dmacs provides indices related to the effects of measurement nonequivalence on observed scores, as described in Nye and Drasgow (2011).

## Installation

``` r
# Install the development version from GitHub:
# install.packages("devtools")
devtools::install_github("ddueber/dmacs")
```

## Usage

##### Mplus output files can be input directly into mplus_dmacs

``` r
# install.packages("MplusAutomation")
library(dmacs)
# Will launch a dialog box for opening an appropriate .out file
# Don't forget to request sampstat output
mplus_dmacs()
```

##### lavaan output can be input directly into lavaan_dmacs

``` r
# install.packages("lavaan")
library(dmacs)
HS.model <- '  visual  =~ x1 + x2 + x3
                textual =~ x4 + x5 + x6
                speed   =~ x7 + x8 + x9 '
fit <- lavaan::cfa(HS.model,
                   data = lavaan::HolzingerSwineford1939,
                   group = "school")
lavaan_dmacs(fit, RefGroup = "Pasteur")
```

## Output

dmacs_summary, mplus_dmacs, and lavaan_dmacs all return a list of measurement nonequivalence indices, including 
effect sizes (DMACS), expected bias in mean item score (ItemDeltaMean), expected bias in mean total scale score (MeanDiff), 
and expected bias in variance of total scale score (VarDiff) due to measurement nonequivalence.

## Supported and unsupported models

dmacs supports use of both continuous and categorical indicators, as well as both unidimensional and multiunidimensional CFA models 
for two or more groups. 

dmacs will produce output for models involving crossloading indicators. When an indicator loads on multiple correlated factors, 
the output will not be correct. When an indicator loads on multiple orthogonal factors (e.g., as in a bifactor model), expected 
bias in item mean will have to be summed across factors by the user. 

Correlated uniquenesses will not cause any problems with computation of dmacs indices.
