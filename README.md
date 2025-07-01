# Example data for flashfmZero

<!-- badges: start -->

[![DOI](https://zenodo.org/badge/1011933350.svg)](https://doi.org/10.5281/zenodo.15785324)
<!-- badges: end -->
 
Here, we provide summary-level [INTERVAL
study](https://doi.org/10.1186/1745-6215-15-363) data for blood cell traits in the *SMIM1* region (1:3441528-3959487). This includes GWAS summary statistics files for 99 blood cell traits (([Astle et
al. 2016](https://doi.org/10.1016/j.cell.2016.10.042), [Akbari et
al. 2023](https://doi.org/10.1038/s41467-023-40679-y))), their trait covariance matrix, and the linkage disequilibrium (LD) matrix for this region (1:3441528-3959487).

Latent factor methods boost power to detect and fine-map signals for many traits and flashfmZero improves fine-mapping resolution by leveraging multiple latent factors.
By using the tools in the [flashfmZero R package](https://jennasimit.github.io/flashfmZero/), latent factor GWAS summary statistics could be calculated from these files, followed by single and multi-trait fine-mapping.    
 
[FlashfmZero](https://doi.org/10.1016/j.xgen.2025.100847) is a computationally efficient approach to **jointly
fine-map signals in any number of uncorrelated quantitative traits**
(i.e. zero correlation), as may result from latent traits estimated from
factor analysis using a varimax rotation. For correlated traits, the
original [flashfm](https://www.nature.com/articles/s41467-021-26364-y)
multi-trait fine-mapping method should be used and cited, which is
available in the flashfmZero package for convenience. Flashfm and flashfmZero
output trait-specific results,leveraging information between traits; for
each trait, credible sets, SNP marginal posterior probabilities of
a causal association (MPP), and multi-SNP model posterior probabilities (PP) are
output.


For more details on flashfmZero, please see:

F Zhou, WJ Astle, AS Butterworth, JL Asimit. (2025). Improved genetic
discovery and fine-mapping resolution through multivariate latent factor
analysis of high-dimensional traits. *Cell Genomics*.
[link](https://doi.org/10.1016/j.xgen.2025.100847)

Website available at: <https://jennasimit.github.io/flashfmZero/>
 
 
 
