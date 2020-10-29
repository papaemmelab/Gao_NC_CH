suppressMessages(library("facetsCH", quietly = TRUE))
suppressMessages(library("argparse", quietly = TRUE))
suppressMessages(library("dplyr", quietly = TRUE))
suppressMessages(library("data.table", quietly = TRUE))
suppressMessages(library("ggplot2", quietly = TRUE))
suppressMessages(library("patchwork", quietly = TRUE))
suppressMessages(library('stringr', quietly = TRUE))
suppressMessages(library('glue', quietly = TRUE))

# R.utils::sourceDirectory('/work/isabl/home/gaot/facets-ch/R', pattern="*R")

parser <- ArgumentParser(description = "Run Facets for Clonal Hematopoiesis (CH)")

parser$add_argument(
    "--db",
    required = TRUE,
    help = "db directory")

parser$add_argument(
    "--pu",
    required = TRUE,
    help = "SNP PileUp output")

parser$add_argument(
    "--sex",
    required = TRUE,
    help = "sex")

parser$add_argument(
    "--prefix",
    required = TRUE,
    help = "OutFile")

parser$add_argument(
    "--outdir",
    required = TRUE,
    help = "OutFile")

parser$add_argument(
    "--panel",
    required = TRUE,
    help = "panel version")

parser$add_argument(
    "--bin",
    required = TRUE,
    type = 'integer',
    help = "SNP bin size")

parser$add_argument(
    "--cval",
    required = TRUE,
    type = 'double',
    help = "critical value for segmentation")

parser$add_argument(
    "--gamma",
    required = TRUE,
    type = 'double',
    help = "het SNP scaling factor")

parser$add_argument(
    "--blacklist",
    required = FALSE,
    default = TRUE,
    type = 'logical',
    help = "use SNP blacklist or not")

parser$add_argument(
    "--bestnormal",
    required = FALSE,
    default = TRUE,
    type = 'logical',
    help = "automatically use best unmatched normal in pool")

parser$add_argument(
    "--correctbaf",
    required = FALSE,
    default = TRUE,
    type = 'logical',
    help = "correct baf using PON")

parser$add_argument(
    "--snpdistance",
    required = FALSE,
    default = TRUE,
    type = 'logical',
    help = "plot using SNP distance instead of genomic")


opt = parser$parse_args()

panel = opt$panel
outdir = opt$outdir
prefix = opt$prefix
sample_pu = fread(opt$pu)

print(
    glue(
        'cval: {opt$cval} \\
        bin: {opt$bin} \\
        blacklist: {opt$blacklist} \\
        bestnormal: {opt$bestnormal} \\
        correctbaf: {opt$correctbaf} \\
        gamma: {opt$gamma} \\
        db: {opt$db}')
)

## centromere positions
acen = fread(paste0(opt$db, '/acen.tsv'))
colnames(acen) = c('chrom', 'start', 'end', 'arm', 'type')

acen = acen %>%
    mutate(chrom = str_remove(chrom, 'chr')) %>%
    mutate(
        chrom = ifelse(chrom == 'X', 23, chrom),
        chrom = ifelse(chrom == 'Y', 24, chrom),
        chrom = as.integer(chrom)
    ) %>%
    arrange(chrom) %>%
    group_by(chrom) %>%
    summarise(cen_start = min(start), cen_end = max(end)) %>%
    ungroup()

## Panel Reference Files
pon_samples = fread(paste0(opt$db, glue('/normals_{panel}.tsv')))
het_whitelist = fread(paste0(opt$db, glue('/het_whitelist_{panel}.tsv')), header = F) %>% pull(V1)
het_bias = fread(paste0(opt$db, glue('/het_bias_{panel}.tsv')))

## SNP lists
if (panel %in% c('IM3', 'IM5', 'IM6')) {
  snp_blacklist_IM6 = fread(paste0(opt$db, '/snp_blacklist_IM6.tsv'), header = F) %>% pull(V1)
  snp_blacklist_IM5 = fread(paste0(opt$db, '/snp_blacklist_IM5.tsv'), header = F) %>% pull(V1)
  snp_blacklist = union(snp_blacklist_IM6, snp_blacklist_IM5)
} else {
  snp_blacklist = fread(paste0(opt$db, glue('/snp_blacklist_{panel}.tsv')), header = F) %>% pull(V1)
}

if (opt$sex %in% c('MALE', 'Male', 'M')) {
  pon_samples = pon_samples %>% filter(Sex == 'Male')
} else {
  pon_samples = pon_samples %>% filter(Sex == 'Female')
}

pon_pileups = lapply(
    paste0(opt$db, '/pon_pileups/', pon_samples[['sample_id']], '_pileup.tsv.gz'),
    function(file) {
      fread(file)
    })

if (!opt$correctbaf) {
  het_bias = NULL
}

sigs = lapply(1:length(pon_pileups), function(i) {

  if (all(sample_pu$File1A == pon_pileups[[i]]$File1A)) {
    data.frame()
  } else {
    get_logr(
        sample_pu,
        pon_pileups[[i]],
        opt$bin,
        opt$cval,
        snp_blacklist = snp_blacklist,
        het_whitelist = het_whitelist,
        het_bias = het_bias,
        gamma = opt$gamma,
        blacklist = opt$blacklist
    ) %>% mutate(normal = i)
  }
}) %>%
    Reduce(rbind, .) %>%
    mutate(snp_index = paste0(chrom, '_', maploc))

errs = sigs %>%
    group_by(normal) %>%
    summarise(std.err = sum(cnlr ^ 2)) %>%
    arrange(std.err)

if (opt$bestnormal) {
  best_normal = errs %>% top_n(-1) %>% pull(normal)
} else {
  best_normal = 1
}

sig = sigs %>% filter(normal == best_normal)

sig = sig %>% mutate(snp_i = 1:n()) %>% mutate(aberrant = FALSE)

sig = sig %>% left_join(
        acen,
        by = 'chrom'
    ) %>%
    mutate(arm = ifelse(maploc < cen_start, 'p', 'q')) %>%
    select(-cen_start, - cen_end) %>%
    mutate(chrom_arm = paste0(chrom, arm))

print('calling abberrant segments')

seg = get_seg(sig)

seg = call_aberrant(sig, seg)

seg[c('type', 'phi', 'err')] = t(mapply(variant_type, seg$cnlr_adj, seg$valor, seg$aberrant))

p = facets_plot(sig, seg, title = prefix, snp_distance = opt$snpdistance)

options(bitmapType = 'cairo', device = 'png')
par(mar = c(0.2, 0.2, 0.2, 0.2))

## Outputs
ggsave(glue('{outdir}/{prefix}.png'), p, width = 10, height = 2.3, dpi = 300, units = "in")
fwrite(seg, glue('{outdir}/{prefix}_seg.tsv'), sep = '\t')
fwrite(sig, glue('{outdir}/{prefix}_sig.tsv'), sep = '\t')
