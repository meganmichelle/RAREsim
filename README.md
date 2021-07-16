# RAREsim

RAREsim is a flexible and accurate rare variant simulation algorithm that emulates real data. Using parameters and haplotypes derived from real sequencing data, RAREsim efficiently simulates the expected variant distribution and enables real variant annotations. The RAREsim paper can be found at https://www.biorxiv.org/content/10.1101/2021.04.13.439644v2.

Here, we have the RAREsim R package. RAREsim is implemented using the R functions along with HAPGEN2 and a Bash script. For a complete implementation of RAREsim, including  generalizable example code, see https://github.com/meganmichelle/RAREsim_Example.

# Installation

The RAREsim R package can be downloaded using devtools.

```{r}
library(devtools)
install_github('meganmichelle/RAREsim')
```

RAREsim also requires [HAPGEN2](https://mathgen.stats.ox.ac.uk/genetics_software/hapgen/hapgen2.html) and sed/gsed.

# Beta

RAREsim is currently in beta.
