#! /bin/bash

# see https://explainshell.com/explain?cmd=set+-euxo%20pipefail
set -eux pipefail

# Test data
PREFIX=tumor
PANEL=FAKE_PANEL
SEX=Male
DIR=/test

# INPUT files
INPUT_BAM=$DIR/tests/data/tumor/tumor.bam
DBSNP=$DIR/tests/data/dbsnp.vcf.gz
DB=$DIR/tests/data/db

# Outputs
OUTDIR=$DIR/tests/output
OUTPUT_PILEUP=$OUTDIR/${PREFIX}_pileup.tsv.gz
OUTPUT_PROFILE=$OUTDIR/${PREFIX}.png
OUTPUT_SEGMENTS=$OUTDIR/${PREFIX}_seg.tsv
OUTPUT_SIGNALS=$OUTDIR/${PREFIX}_sig.tsv

# Clean Output test directory
rm -f $OUTDIR/*

# Run Tumor dbsnps Pileup
# docker run -it -v /Users:/Users --entrypoint 'snp-pileup' papaemmelab/facets-ch \
snp-pileup -g -A -q15 -r20 -Q20 -P100 -v $DBSNP $OUTPUT_PILEUP $INPUT_BAM

# Run Facets CH
# docker run -it -v /Users:/Users papaemmelab/facets-ch \
Rscript $DIR/facets-ch.R \
  --db $DB \
  --pu $OUTPUT_PILEUP \
  --sex $SEX \
  --outdir $OUTDIR \
  --prefix $PREFIX \
  --panel $PANEL \
  --cval 30 \
  --bin 750 \
  --gamma 0.25 \
  --blacklist TRUE \
  --snpdistance TRUE \
  --correctbaf TRUE \
  --bestnormal TRUE

# Check output files exist
[[ -f $OUTPUT_PILEUP && -f $OUTPUT_PROFILE && -f $OUTPUT_SEGMENTS && -f $OUTPUT_SIGNALS ]] \
  && echo "Output files exist. TESTS PASSED" \
  || echo "Missing ouput files. TESTS FAILED."
