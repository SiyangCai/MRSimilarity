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
The main function `MR_similarity_analysis` generates causal estimates and a scatterplot to provide an initial visualisation of results obtained from multiple selected methods. Test statistics, including the $\bar{F}$ statistic, $I^2$ statistic, and MR-Egger intercept (p-value), are provided to investigate potential dissimilarity arising from weak instrument bias or horizontal pleiotropy.
\\
Similarity statistics are then applied to test the overall hypothesis that all selected methods produce similar causal estimates. If dissimilarity is detected in the overall hypothesis, users can further investigate the source of variation through pairwise or group-wise comparisons between methods.

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
