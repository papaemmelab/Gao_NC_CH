# FACETS-CH

Clonal hematopoiesis (CH) calling in targeted sequencing data.

In the preprocessing step, FACETS-CH leverages a large panel of normal (PON) to reduce noise. It then iteratively computes sample-specific noise profiles and uses a bivariate Wald test to identify aberrant segments.

Some of the preprocessing and segmentation steps are based on [FACETS] (Shen et al, NAR 2016). FACETS-CH is designed for detecting subclonal chromosomal alterations from normal tissue sequencing, where most of the genome should be in diploid state. It is not fit for analyzing tumor genomes that have complex karyotypes.

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

FACETS-CH runs in 2 steps: pileup step and variant calling step.

### 1. Pileup step

```R
run_pileup = function(nid, nbam, outdir = '/work/isabl/home/gaot/ch_cnv_pileup', q = 15, Q = 20, r = 20) {

    outfile = glue('{outdir}/{nid}_pileup.tsv.gz')
    logfile = glue('{outdir}/{nid}.logs')

    glue(
        "bsub \\
        -M 10 \\
        -n 1 \\
        -We 59 \\
        -oo {logfile} \\
       'rm -f {outfile} && singularity exec --workdir $TMP_DIR/gaot_docker-facets_v0.1.1_`uuidgen` \\
            --bind /ifs:/ifs --bind /juno:/juno \\
            --bind /juno/work:/work \\
            --bind /juno/res:/res \\
            /juno/work/isabl/local/docker-facets/v0.1.1/docker-facets-v0.1.1.simg \\
            snp-pileup \\
                -g -q{q} -Q{Q} -P100 -A \\
                -r{r} \\
                -v /work/isabl/public/facets-ch/pon/dbsnp_138.b37.vcf.gz \\
                {outfile} \\
                {nbam}'")
}
```

### 2. Run FACETS-CH

```R
run_facets = function(prefix,
                      pu,
                      sex,
                      panel,
                      db,
                      cval = 30,
                      bin = 750,
                      memory = 4,
                      gamma = 0.25,
                      blacklist = TRUE,
                      bestnormal = TRUE,
                      snpdistance = TRUE,
                      correctbaf = TRUE,
                      outdir = '/work/isabl/home/gaot/ch_cnv/results') {

    glue(
        "bsub \\
        -n 1 \\
        -M {memory} \\
        -We 59 \\
        -oo {outdir}/{sid}.logs \\
       '/work/isabl/home/gaot/R/R-3.6.1/bin/Rscript /work/isabl/home/gaot/facets-ch/facets-ch.R \\
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
