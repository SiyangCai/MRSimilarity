# MRSimilarity

This is the R package for a designed method to apply a heterogeneity statistic on a set of causal estimates in Mendelian Randomisation using summary statistics from GWAS studies. The aim is to test whether those estimates are similar.

For details of our methods, please refer to our paper.


## Installation
From Rstudio:

```
library(devtools)
install_github("SiyangCai/MRSimilarity")
```

## Usage
Using the main function `MR_similarity_analysis` generates causal estimates with scatterplot to initially visualise results from multiple selected methods. Test statistics, e.g. $\bar{F}$ statistics, $I^2$ statistics and MR-Egger intercept are provided to investigate potential dissimilarity due to bias.
Similarity statistics are then applied on overall hypothesis that all methods are similar. Users are allowed to further investigation through pairwise or group-wise copmparisons if dissimilarity is found in overall hypothesis.
An example of using the `MR_similarity_analysis` is demonstrated with a sample summary dataset as follow:

```
library(MRSimilarity)

# Load the sample dataset: LDL-C (exposure) vs. CAD (outcome)
data(LDL-CAD)

# Perform similarity analysis
# Recommend using non-parametric bootstrap with a large bootstrap sample
res = MR_similarity_analysis(data = LDL_CAD,
                               xlab = "Instrument-LDL-cholesterol effect",
                               ylab = "Instrument-CAD effect",
                               title_text = "Causal effect estimates of LDL on CAD",
                               parametric_bootstrap = FALSE,
                               n.boot = 10000)
print(res)
```
For more information, please refer to the help page in the R package.
