# RAREsim

The RAREsim R package is used to calculate the estimated number of variants per MAC bin. Throught the _Number of Variants_ and _AFS_ functions. Once the expected number of variants have been calculated, RAREsim can prune the over-simulated haplotypes created from simulating with HAPGEN2 and information at all  sequencing bases - including monomorphic bases. Finally, the simulated haplotypess are probabilistically pruned to match what is expected via the RAREsim package found here: https://github.com/ryanlayer/raresim. A complete example with all steps can be found on the [Example page](https://github.com/meganmichelle/RAREsim_Example)

# Installation

The RAREsim R package can be downloaded using devtools.

```{r}
library(devtools)
install_github('meganmichelle/RAREsim')
```

RAREsim also requires [HAPGEN2](https://mathgen.stats.ox.ac.uk/genetics_software/hapgen/hapgen2.html).

# More Information

RAREsim is a flexible and accurate rare variant simulation algorithm that emulates real data. Using parameters and haplotypes derived from real sequencing data, RAREsim efficiently simulates the expected variant distribution and enables real variant annotations. The RAREsim paper can be found at https://www.biorxiv.org/content/10.1101/2021.04.13.439644v2.

# Beta

RAREsim is currently in beta.
