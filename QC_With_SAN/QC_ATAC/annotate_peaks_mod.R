### Annotate peaks w GenomicRanges ###
suppressPackageStartupMessages({
  library(argparse)
  library(Signac)
  library(GenomicRanges)
  library(data.table)
  library(EnsDb.Hsapiens.v86)
  library(tidyverse)
  library(EnsDb.Mmusculus.v79)
  })

# parser <- ArgumentParser()
# parser$add_argument("peak_file", type="character",
#                     help = "Path to tsv file contaning peak ids")
# parser$add_argument("outdir",
#                     type="character",
#                     help = "Path to output directory")
# parser$add_argument("--genome", default = "hg38",
#                     type="character",
#                     help = "Reference genome ID")
# args <- parser$parse_args()

# peak_file <- args$peak_file
# outdir <- args$outdir
# genome <- args$genome

peak_file <- '/nfs/team205/heart/cellatac/tic-1158/results200k/peak_matrix/peaks.txt'
outdir <- '/lustre/scratch117/cellgen/team205/dm19/kazu/notebooks/annotate_peaks'
genome <- 'hg38'

#' Get annotation overlap
#'
#' Computes overlap between GenomicRanges object and annotation coordinates
#'
#' @param gr GenomicRanges object
#' @param anno_gr GenomicRanges object of annotations (usually loaded from ensembldb package)
#' @param anno_name character of name of annotation (e.g. genes, exons, promoters ...)
#' @param minoverlap minimum no of overlapping bps to call an overlap
#'
#' @return GenomicRanges object with overlap stored in metadata column
#'
get_annotation_overlap <- function(gr, anno_gr, anno_name, minoverlap = 10){
  ## Clean annotation object
  genome(anno_gr) <- "hg38"
#  seqlevelsStyle(anno_gr) <- 'UCSC'
  # chr_names <- seqlevels(anno_gr)
  # seqlevels(anno_gr) <- paste0('chr', chr_names)
  anno_gr <- keepStandardChromosomes(anno_gr, pruning.mode = 'coarse')
  ## Find overlaps
  overlaps <- findOverlaps(gr, anno_gr, minoverlap = minoverlap)
  gr@elementMetadata[anno_name] <- 0
  gr@elementMetadata[anno_name][queryHits(overlaps),] <- 1
  return(gr)
}

#' Annotate GenomicRanges
#'
#' add width and overlap with annotations to metadata columns of GenomicRanges object
#'
#' @param gr GenomicRanges object to annotate
#' @param EnsDb.genome Ensembl annotation package to use for annotation. See details for installation (default = EnsDb.Hsapiens.v86, for hg38).
#' @param blacklist_gr GenomicRanges object of ENCODE blacklist regions (default: blacklist_hg38 from Signac package)
#' @param minoverlap minimum no of overlapping bps to call an overlap
#'
#' @details You can install EnsDb annotation packages by running
#' \code{BiocManager::install("EnsDb.Hsapiens.v86")}
#'
#' @import ensembldb
#' @import dplyr
#' @import GenomicRanges
#' @import GenomeInfoDb
#' @import Signac
#' @import ensembldb
#' @importFrom stats setNames
#'
#' @export
annotate_gr <- function(gr, EnsDb.genome = EnsDb.Hsapiens.v86, 
  blacklist_gr = Signac::blacklist_hg38, minoverlap=10){
  ## Compute width
  gr$peak_width <- width(gr)

  ## Compute overlap with annotations
  exon_coords <- ensembldb::exons(EnsDb.genome, filter = ~ gene_biotype == "protein_coding")
  genes_coords <- ensembldb::genes(EnsDb.genome, filter = ~ gene_biotype == "protein_coding")
  promoter_coords <- ensembldb::promoters(EnsDb.genome, filter = ~ gene_biotype == "protein_coding", upstream=2000, downstream = 100)
  anno_list <- list(exon=exon_coords, gene=genes_coords, promoter=promoter_coords)
  for (i in seq_along(anno_list)) {
    anno_gr = anno_list[[i]]
    chr_names <- seqlevels(anno_gr)
    seqlevels(anno_gr) <- paste0('chr', chr_names)
    gr <- get_annotation_overlap(gr, anno_gr = anno_gr, anno_name = names(anno_list)[i], minoverlap=minoverlap)
  }
  gr$annotation <- data.frame(gr@elementMetadata) %>%
    mutate(annotation = case_when(exon==1 ~ "exon",
                                  promoter==1 ~ "promoter",
                                  gene==1 & exon==0 ~ "intron",
                                  TRUE ~ "intergenic"
                                  )) %>%
    pull(annotation)

  ## Overlapping gene
  chr_names <- seqlevels(genes_coords)
  seqlevels(genes_coords) <- paste0('chr', chr_names)
  # seqlevelsStyle(genes_coords) <- 'UCSC'
  genes_coords <- keepStandardChromosomes(genes_coords, pruning.mode = 'coarse')
  genespromoters_coords <- Extend(genes_coords, upstream = 2000)
  overlap <- findOverlaps(gr, genespromoters_coords, minoverlap = minoverlap)
  gr$gene_name <- NA
  gr$gene_id <- NA
  gr$gene_name[queryHits(overlap)] <- genespromoters_coords$gene_name[subjectHits(overlap)]
  gr$gene_id[queryHits(overlap)] <- genespromoters_coords$gene_id[subjectHits(overlap)]

  ## Distance to nearest TSS
  transcript_coords <- ensembldb::transcripts(EnsDb.genome, filter = ~ gene_biotype == "protein_coding")
  chr_names <- seqlevels(transcript_coords)
  seqlevels(transcript_coords) <- paste0('chr', chr_names)
  # seqlevelsStyle(transcript_coords) <- 'UCSC'
  transcript_coords <- keepStandardChromosomes(transcript_coords, pruning.mode = 'coarse')
  tss_coords <- resize(unlist(range(split(transcript_coords, ~ tx_id))), width=1)
  tss_distance <- distanceToNearest(gr, tss_coords)
  gr$tss_distance <- tss_distance@elementMetadata$distance

  ## Overlap with blacklist
  gr <- get_annotation_overlap(gr, anno_gr = blacklist_gr, anno_name = "ENCODE_blacklist")

  return(gr)
  }

# prefix <- "GSE126074_P0_BrainCortex_SNAREseq_"
# atac.id <- "chromatin"
# peak_file <- paste0(data.dir, "ATAC_raw/", prefix, atac.id, ".peaks.tsv.gz")
# genome = "mm10"

if (genome=="hg38") {
  EnsDb.genome = EnsDb.Hsapiens.v86
  blacklist_gr = blacklist_hg38
} else if (genome=="mm10"){
  EnsDb.genome = EnsDb.Mmusculus.v79 
  blacklist_gr = blacklist_mm10
} else {
  stop("Unrecognized genome ID. Supported genomes: hg38, mm10")
}

peak_ids <- fread(peak_file, header = F)[["V1"]]
peaks_gr <- StringToGRanges(peak_ids, sep=c(":","-"))

# gr <- peaks_gr
peaks_gr <- annotate_gr(peaks_gr, EnsDb.genome, blacklist_gr)

peaks_df <- data.frame(peaks_gr@elementMetadata) 
peaks_df["peak_id"] <- peak_ids

write.csv(peaks_df, paste0(outdir, "/","ATACpeaks_annotation.csv" ))
