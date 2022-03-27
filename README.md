# RAREsim

The functions within the RAREsim R package are used to estimate the expected number of rare variants per minor allele count (MAC) bin via the *Number of Variants* and *AFS* functions. Users have the flexibility to estimate the number of variants with default parameters, user-defined parameters, or parameters estimated to match user-provided target data. The default parameters are provided for four ancestries and are expected to perform well for sample sizes up to 125,000.

The RAREsim R package is used to calculate the estimated number of variants per MAC bin. Throught the _Number of Variants_ and _AFS_ functions. The output of the RAREsim functions can be directly input into the pruning package for RAREsim, implemented within python: https://github.com/ryanlayer/raresim. The pruning process requires the expected number of variants per rare MAC bin and over-simulated haplotypes created from simulating from HAPGEN2 with information at all sequencing bases - including monomorphic bases. A complete example with all steps can be found on the Example page: https://github.com/meganmichelle/RAREsim_Example. 

# Installation

The RAREsim R package can be downloaded using devtools.

```{r}
library(devtools)
install_github('meganmichelle/RAREsim')
```

RAREsim also requires [HAPGEN2](https://mathgen.stats.ox.ac.uk/genetics_software/hapgen/hapgen2.html).

The pruning process is done with the [RAREsim python package](https://github.com/ryanlayer/raresim).

# More Information

RAREsim is a flexible and accurate rare variant simulation algorithm that emulates real data. Using parameters and haplotypes derived from real sequencing data, RAREsim efficiently simulates the expected variant distribution and enables real variant annotations. See the full article here: https://doi.org/10.1016/j.ajhg.2022.02.009.

# Beta

RAREsim is currently in beta.
