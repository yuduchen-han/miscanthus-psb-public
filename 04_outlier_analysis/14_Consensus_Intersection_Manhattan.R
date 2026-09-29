# 科学目的：严格取 pcadapt ∩ RF Habitat ∩ RDA 并为最终 21 SNP 加基因/功能注释。
# 输入：三个独立候选表、最终 Genlight、reference GFF 与功能注释表。
# 输出：21 SNP 清单、Venn、完整/精简注释、供 21-SNP 分析使用的 RDS。
# 流程：验证 SNP universe → 严格交集 → gene overlap/nearest → 功能注释；缺注释不删 SNP。

suppressPackageStartupMessages({
  library(adegenet); library(dplyr); library(GenomicRanges)
  library(rtracklayer); library(VennDiagram); library(grid)
})

# Load data: reference paths are relative to PROS_Pipeline; confirm these sources for release.
dir_out <- "04_Results_Archive/11_Consensus_Annotation_HabitatNew"
dir.create(dir_out, recursive=TRUE, showWarnings=FALSE)
gff_path <- paste0("../../2_data/referece_datasat/",
  "miscanthus_sinensis_v7.1_JGI/Phytozome/PhytozomeV12/",
  "early_release/Msinensis_497_v7.1/annotation/Msinensis_497_v7.1.gene.gff3.gz")
anno_path <- paste0("../../2_data/referece_datasat/",
  "miscanthus_sinensis_v7.1_JGI/Phytozome/PhytozomeV14/",
  "Msinensis/v7.1/annotation/Msinensis_497_v7.1.P14.annotation_info.txt.gz")
gl_main <- readRDS("00_rds_modules/02_final_genlight.rds")
stopifnot(nInd(gl_main)==359, nLoc(gl_main)==8625, !anyDuplicated(locNames(gl_main)))
snp_universe <- locNames(gl_main)
extract_ids <- function(path, preferred) {
  df <- read.csv(path, stringsAsFactors=FALSE)
  col <- intersect(preferred, names(df))[1]
  stopifnot(!is.na(col))
  ids <- as.character(df[[col]])
  unique(ids[!is.na(ids) & trimws(ids)!=""])
}
pcadapt_snps <- extract_ids("04_Results_Archive/09_pcadapt/Table_pcadapt_Outliers_Robust.csv",
                            c("SNP_ID","Candidate_SNP_ID","Outlier_SNP_ID"))
rf_snps <- extract_ids("03_Public_Code/data/handoff/RF_Habitat_final_candidates.csv",
                       c("SNP_ID","Candidate_SNP_ID"))
rda_snps <- extract_ids("04_Results_Archive/10_RDA_HabitatNew/Table_RDA_Candidate_SNPs_Pure_List.csv",
                        c("SNP_ID","Candidate_SNP_ID"))
stopifnot(all(pcadapt_snps %in% snp_universe), all(rf_snps %in% snp_universe),
          all(rda_snps %in% snp_universe))

# Prepare data: candidate membership is fixed before any annotation lookup.
pcadapt_rf <- intersect(pcadapt_snps, rf_snps)
pcadapt_rda <- intersect(pcadapt_snps, rda_snps)
rf_rda <- intersect(rf_snps, rda_snps)
triple_snps <- Reduce(intersect, list(pcadapt_snps, rf_snps, rda_snps))
stopifnot(length(triple_snps)==21)
overlap_summary <- data.frame(
  Comparison=c("pcadapt","RF","RDA","pcadapt_RF","pcadapt_RDA","RF_RDA","pcadapt_RF_RDA"),
  N_SNPs=c(length(pcadapt_snps),length(rf_snps),length(rda_snps),
           length(pcadapt_rf),length(pcadapt_rda),length(rf_rda),length(triple_snps)))

# Analysis: sequence-name harmonization and overlap/nearest rules match the raw method.
stopifnot(file.exists(gff_path))
gff <- rtracklayer::import(gff_path)
genes <- gff[gff$type=="gene"]
stopifnot(length(genes)>0)
idx <- match(triple_snps, locNames(gl_main))
chr <- as.character(gl_main@chromosome[idx])
numeric_chr <- grepl("^[0-9]+$", chr)
chr[numeric_chr] <- sprintf("Chr%02d", as.numeric(chr[numeric_chr]))
coords <- data.frame(SNP_ID=triple_snps, Chr=chr,
                     SNP_Pos=as.numeric(gl_main@position[idx]))
stopifnot(!anyNA(coords$SNP_Pos))
available <- coords$Chr %in% unique(as.character(seqnames(genes)))
annotatable <- coords[available, , drop=FALSE]
no_gff <- coords[!available, , drop=FALSE]
empty <- data.frame(SNP_ID=character(), Gene_ID=character(), Chr=character(),
                    SNP_Pos=numeric(), Gene_Start=numeric(), Gene_End=numeric(),
                    Distance=numeric(), Hit_Type=character())
annotated <- empty
if (nrow(annotatable)>0) {
  snp_gr <- GRanges(seqnames=annotatable$Chr,
                    ranges=IRanges(start=annotatable$SNP_Pos, end=annotatable$SNP_Pos),
                    SNP_ID=annotatable$SNP_ID)
  genes_clean <- genes[as.character(seqnames(genes)) %in%
                         unique(as.character(seqnames(snp_gr)))]
  ids_clean <- if (!is.null(genes_clean$Name)) as.character(genes_clean$Name) else
    as.character(genes_clean$ID)
  hits <- findOverlaps(snp_gr, genes_clean)
  inside <- empty
  if (length(hits)>0) {
    q <- queryHits(hits); s <- subjectHits(hits)
    inside <- data.frame(SNP_ID=as.character(snp_gr$SNP_ID[q]), Gene_ID=ids_clean[s],
                         Chr=as.character(seqnames(genes_clean[s])), SNP_Pos=start(snp_gr[q]),
                         Gene_Start=start(genes_clean[s]), Gene_End=end(genes_clean[s]),
                         Distance=0, Hit_Type="inside_gene") %>%
      distinct(SNP_ID, .keep_all=TRUE)
  }
  outside_ids <- setdiff(annotatable$SNP_ID, inside$SNP_ID)
  near <- empty
  if (length(outside_ids)>0) {
    outside_gr <- snp_gr[match(outside_ids, snp_gr$SNP_ID)]
    nearest_idx <- GenomicRanges::nearest(outside_gr, genes_clean, ignore.strand=TRUE)
    stopifnot(!anyNA(nearest_idx))
    nearest_genes <- genes_clean[nearest_idx]
    nearest_ids <- if (!is.null(nearest_genes$Name)) as.character(nearest_genes$Name) else
      as.character(nearest_genes$ID)
    near <- data.frame(SNP_ID=outside_ids, Gene_ID=nearest_ids,
                       Chr=as.character(seqnames(nearest_genes)), SNP_Pos=start(outside_gr),
                       Gene_Start=start(nearest_genes), Gene_End=end(nearest_genes),
                       Distance=as.numeric(GenomicRanges::distance(outside_gr, nearest_genes,
                                                                   ignore.strand=TRUE)),
                       Hit_Type="nearest_gene")
  }
  annotated <- bind_rows(inside, near)
  stopifnot(nrow(annotated)==nrow(annotatable),
            setequal(annotated$SNP_ID, annotatable$SNP_ID))
}
unannotated <- empty
if (nrow(no_gff)>0) {
  unannotated <- data.frame(SNP_ID=no_gff$SNP_ID, Gene_ID=NA_character_,
                            Chr=no_gff$Chr, SNP_Pos=no_gff$SNP_Pos,
                            Gene_Start=NA_real_, Gene_End=NA_real_, Distance=NA_real_,
                            Hit_Type="no_gene_annotation")
}
triple_anno <- bind_rows(annotated, unannotated)
stopifnot(nrow(triple_anno)==length(triple_snps),
          !anyDuplicated(triple_anno$SNP_ID), setequal(triple_anno$SNP_ID, triple_snps))
triple_anno <- triple_anno[match(triple_snps, triple_anno$SNP_ID), , drop=FALSE]
if (file.exists(anno_path)) {
  anno_df <- read.delim(anno_path, sep="\t", header=TRUE, stringsAsFactors=FALSE)
  stopifnot("locusName" %in% names(anno_df))
  triple_anno <- left_join(triple_anno, distinct(anno_df, locusName, .keep_all=TRUE),
                           by=c("Gene_ID"="locusName"))
  stopifnot(nrow(triple_anno)==length(triple_snps),
            identical(as.character(triple_anno$SNP_ID), as.character(triple_snps)))
} else warning("Functional annotation unavailable; consensus membership unchanged.")
compact_cols <- intersect(c("SNP_ID","Gene_ID","Chr","SNP_Pos","Gene_Start",
                            "Gene_End","Distance","Hit_Type","Best.hit.arabi.defline"),
                          names(triple_anno))
triple_compact <- triple_anno[, compact_cols, drop=FALSE]

# Save results: 21-SNP GST/RF consume triple_snps and triple_annotation_full.
write.csv(data.frame(SNP_ID=triple_snps), file.path(dir_out, "Table_Triple_Overlap_SNPs.csv"), row.names=FALSE)
write.csv(overlap_summary, file.path(dir_out, "Table_Consensus_Overlap_Summary.csv"), row.names=FALSE)
write.csv(triple_anno, file.path(dir_out, "Table_Triple_Overlap_Annotation_Full.csv"), row.names=FALSE)
write.csv(triple_compact, file.path(dir_out, "Table_Triple_Overlap_Annotation.csv"), row.names=FALSE)
venn_plot <- venn.diagram(x=list(pcadapt=pcadapt_snps, RF=rf_snps, RDA=rda_snps),
  filename=NULL, fill=c("#8DA0CB","#66C2A5","#FC8D62"), alpha=0.45,
  cex=1.35, cat.cex=1.15, cat.fontface="plain", margin=0.08, main=NULL)
pdf(file.path(dir_out, "Figure_pcadapt_RF_RDA_Consensus_Venn.pdf"), width=8, height=7)
grid.newpage(); grid.draw(venn_plot); dev.off()
saveRDS(list(pcadapt_snps=pcadapt_snps, rf_snps=rf_snps, rda_snps=rda_snps,
             pairwise_overlap=list(pcadapt_RF=pcadapt_rf, pcadapt_RDA=pcadapt_rda, RF_RDA=rf_rda),
             triple_snps=triple_snps, triple_annotation=triple_compact,
             triple_annotation_full=triple_anno),
        "00_rds_modules/11_Consensus_RF_pcadapt_RDA.rds")
