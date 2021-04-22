library(rtracklayer)
library(GenomicRanges)

readPQS <- function(file_path)
{
  tmp <- read.table(file_path, 
                    header = F, sep = '\t',
                    col.names = c('seqnames', 'start', 'end', 'name', 'length', 'strand', 'seq')
  )
  
  tmp$seqnames <- gsub('(.*)\\sN.*', '\\1', tmp$seqnames)
  tmp <- tmp[tmp$seqnames %in% paste0('chr', c(1:100, 'X', 'Y')), ]
  makeGRangesFromDataFrame(tmp)
}


bed_to_bw <- function(file_path, .name, genome) {
  if (length(file_path) > 1) {
    gr <- Reduce("c", sapply(file_path, readPQS))
  } else {
    gr <- readPQS(file_path)
  }
  
  export.bw(coverage(gr), paste0(.name, "_stack.bw"))
  
  gr <- reduce(gr, ignore.strand = T)
  gr$score <- 1
  if (genome == "hg19") {
    seqlengths(gr) <- seqlengths(BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19)[seqlevels(gr)]
  } 
  if (genome == "mm9") {
    seqlengths(gr) <- seqlengths(BSgenome.Mmusculus.UCSC.mm9::BSgenome.Mmusculus.UCSC.mm9)[seqlevels(gr)]
  }
  export.bw(gr, paste0(.name, "_binary.bw"))
}


bed_to_bw("data/hg19/canonical_PQS_hg19.bed", "canonical_PQS_hg19")
bed_to_bw("data/hg19/extended_canonical_PQS_hg19.bed", "extended_canonical_PQS_hg19")
bed_to_bw("data/hg19/two_tetrads_PQS_hg19.bed", "two_tetrads_PQS_hg19")

bed_to_bw("data/hg19/AAAA_hg19.bed", "AAAA_hg19", "hg19")
bed_to_bw("data/hg19/AAAB_hg19.bed", "AAAB_hg19", "hg19")
bed_to_bw("data/hg19/AABB_hg19.bed", "AABB_hg19", "hg19")
bed_to_bw("data/hg19/ABBB_hg19.bed", "ABBB_hg19", "hg19")
bed_to_bw("data/hg19/ABAB_hg19.bed", "ABAB_hg19", "hg19")
bed_to_bw("data/hg19/ABAA_hg19.bed", "ABAA_hg19", "hg19")
bed_to_bw("data/hg19/ABBA_hg19.bed", "ABBA_hg19", "hg19")
bed_to_bw("data/hg19/AABA_hg19.bed", "AABA_hg19", "hg19")

bed_to_bw(list.files("data/hg19", pattern = "(AA|AB)", full.names = T), "tranStrand_hg19", "hg19")


bed_to_bw("data/mm9/canonical_PQS_mm9.bed", "canonical_PQS_mm9", "mm9")
bed_to_bw("data/mm9/extended_canonical_PQS_mm9.bed", "extended_canonical_PQS_mm9", "mm9")
bed_to_bw("data/mm9/two_tetrads_PQS_mm9.bed", "two_tetrads_PQS_mm9", "mm9")

bed_to_bw("data/mm9/AAAA_mm9.bed", "AAAA_mm9", "mm9")
bed_to_bw("data/mm9/AAAB_mm9.bed", "AAAB_mm9", "mm9")
bed_to_bw("data/mm9/AABB_mm9.bed", "AABB_mm9", "mm9")
bed_to_bw("data/mm9/ABBB_mm9.bed", "ABBB_mm9", "mm9")
bed_to_bw("data/mm9/ABAB_mm9.bed", "ABAB_mm9", "mm9")
bed_to_bw("data/mm9/ABAA_mm9.bed", "ABAA_mm9", "mm9")
bed_to_bw("data/mm9/ABBA_mm9.bed", "ABBA_mm9", "mm9")
bed_to_bw("data/mm9/AABA_mm9.bed", "AABA_mm9", "mm9")

bed_to_bw(list.files("data/mm9", pattern = "(AA|AB)", full.names = T), "tranStrand_mm9", "mm9")