# FACETS-CH

FACETS-CH is an algorithm to detect mosaic chromosomal alterations (mCAs) from targeted sequencing data. It iteratively computes sample-specific noise profiles and identifies aberrant segments using a bivariate Wald test. Some of the preprocessing and segmentation steps are based on [FACETS] (Shen et al, NAR 2016). 

FACETS-CH is designed for detecting subclonal chromosomal alterations from normal tissue sequencing, where most of the genome should be in diploid state. It is not fit for analyzing tumor genomes that have complex karyotypes, for which we recommend using the original FACETS.

![](/samples.png)

## Install

Make sure you install first [FACETS], and its dependencies. You can install the current version (along with the vignette) using the command:

```R
devtools::install_github("mskcc/facets", build_vignettes = TRUE)
```

pctGCdata is a required package. So install that also (needs to be done only once)

```R
devtools::install_github("mskcc/pctGCdata")
```

If you get an error message about pctGCdata use

```R
devtools::install_github("veseshan/pctGCdata")
```

Finally, install this package inside an R session:

```R
devtools::install_local(build_vignettes = FALSE, upgrade="never")
```

## Test Installation

Run containerized tests with:

```bash
./test-container.sh
```

## Usage

The following steps are required to run analysis using FACETS-CH.
1. Build a panel of normal (PON) composed of unmatched normal samples sequenced using the same panel. We recommend >100 samples with both genders with optimal performance.

2. Run pileup using `snp-pileup` from the FACETS package, using a DBSNP VCF.
Example:
```
snp-pileup \\
    -g -q15 -Q20 -P100 -A \\
    -r20 \\
    -v dbsnp_138.b37.vcf.gz \\
    {outfile} \\
    {nbam}'
```
3. Run the FACETS-CH wrapper script `facets-ch.R`.
Example:
```
Rscript facets-ch.R \\
    --db {db} \\
    --pu {pu} \\
    --sex {sex} \\
    --outdir {outdir} \\
    --prefix {prefix} \\
    --panel {panel} \\
    --cval {cval} \\
    --bin {bin} \\
    --blacklist {blacklist} \\
    --snpdistance {snpdistance} \\
    --correctbaf {correctbaf} \\
    --gamma {gamma} \\
    --bestnormal {bestnormal}'")

}
```

<!-- References -->

[FACETS]: https://github.com/mskcc/facets
