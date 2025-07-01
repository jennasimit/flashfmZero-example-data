# Example data for flashfmZero
 
Here, we provide summary-level [INTERVAL
study](https://doi.org/10.1186/1745-6215-15-363) data for fine-mapping associations with blood cell traits in the *SMIM1* region. This includes GWAS summary statistics files for 99 blood cell traits (([Astle et
al. 2016](https://doi.org/10.1016/j.cell.2016.10.042), [Akbari et
al. 2023](https://doi.org/10.1038/s41467-023-40679-y))), their trait covariance matrix, and the linkage disequilibrium (LD) matrix for this region.    
 
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

We also provide an approach to **estimate latent factor GWAS summary
statistics using only summary level data** - observed trait GWAS summary
statistics, observed trait covariance matrix.

For more details on flashfmZero, please see:

F Zhou, WJ Astle, AS Butterworth, JL Asimit. (2025). Improved genetic
discovery and fine-mapping resolution through multivariate latent factor
analysis of high-dimensional traits. *Cell Genomics*.
[link](https://doi.org/10.1016/j.xgen.2025.100847)

Website available at: <https://jennasimit.github.io/flashfmZero/>
 
 
 
