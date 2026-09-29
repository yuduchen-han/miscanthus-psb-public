# 科学目的：由最终 359 × 8625 VCF 构建 individual-level supplementary ML tree。
# 输入：最终 VCF、Genlight、外部 vcf2phylip.py、PATH 中的 iqtree2。
# 输出：variable-sites treefile 与供最终 plot 脚本使用的 midpoint-rooted tree RDS。
# 流程：按 Genlight 重排 VCF → PHYLIP → IQ-TREE 生成 variable sites → MFP+ASC → midpoint。

library(adegenet)
library(vcfR)
library(ape)
library(phytools)

dir_out <- "04_Results_Archive/03_ML_Tree"
dir.create(dir_out, recursive=TRUE, showWarnings=FALSE)
vcf2phylip_script <- "01_Data/vcf2phylip.py" # 外部 dependency，不在本层重写
stopifnot(file.exists(vcf2phylip_script), nzchar(Sys.which("iqtree2")))
clean_vcf_path <- file.path(dir_out, "Miscanthus_core_filtered.vcf.gz")
phylip_path <- file.path(dir_out, "Miscanthus_core_filtered.min4.phy")
qi_prefix <- file.path(dir_out, "Miscanthus_ML_final_8625")
varsites_path <- paste0(qi_prefix, ".varsites.phy")
tree_path <- paste0(qi_prefix, ".treefile")

# Load data
gl_main <- readRDS("00_rds_modules/02_final_genlight.rds")
vcf_raw <- read.vcfR("01_Data/vcf_260902_8625_fin/populations_final_359_8625.vcf",
                     verbose=FALSE)
stopifnot(nInd(gl_main)==359, nLoc(gl_main)==8625,
          ncol(vcf_raw@gt)-1==359, nrow(vcf_raw@fix)==8625,
          identical(as.character(gl_main@other$meta$SampleID), as.character(indNames(gl_main))))

# Prepare data: use the same suffix stripping as Core Genlight and its sample order.
stopifnot("FORMAT" %in% colnames(vcf_raw@gt))
clean_names <- sub("_.*$", "", setdiff(colnames(vcf_raw@gt), "FORMAT"))
stopifnot(!anyDuplicated(clean_names), setequal(indNames(gl_main), clean_names))
colnames(vcf_raw@gt) <- c("FORMAT", clean_names)
vcf_filtered <- vcf_raw[, c("FORMAT", indNames(gl_main))]
stopifnot(identical(colnames(vcf_filtered@gt)[-1], indNames(gl_main)))
vcfR::write.vcf(vcf_filtered, file=clean_vcf_path)

# Analysis: vcf2phylip writes its output into cwd. Keep the cwd change scoped.
script_abs <- normalizePath(vcf2phylip_script, mustWork=TRUE)
vcf_abs <- normalizePath(clean_vcf_path, mustWork=TRUE)
old_wd <- getwd()
python_status <- tryCatch({
  setwd(normalizePath(dir_out, mustWork=TRUE))
  system2("python3", c(shQuote(script_abs), "-i", shQuote(vcf_abs)))
}, finally=setwd(old_wd))
stopifnot(python_status==0, file.exists(phylip_path))
phylip_dims <- as.integer(strsplit(trimws(readLines(phylip_path, n=1)), "\\s+")[[1]])
stopifnot(identical(phylip_dims[1:2], c(359L, 8625L)))

# +ASC preflight may exit nonzero on invariant sites; its *.varsites.phy is the product.
if (file.exists(varsites_path)) file.remove(varsites_path)
system2("iqtree2", c("-s", shQuote(phylip_path), "-m", "MFP+ASC",
                      "-pre", shQuote(qi_prefix), "--redo"))
stopifnot(file.exists(varsites_path))
variable_dims <- as.integer(strsplit(trimws(readLines(varsites_path, n=1)), "\\s+")[[1]])
stopifnot(identical(variable_dims[1:2], c(359L, 3047L)))

# Final model: 1000 UFBoot + 1000 SH-aLRT, AUTO threads.
if (file.exists(tree_path)) file.remove(tree_path)
status <- system2("iqtree2", c("-s", shQuote(varsites_path), "-m", "MFP+ASC",
                              "-B", "1000", "--alrt", "1000", "-T", "AUTO",
                              "-pre", shQuote(qi_prefix), "--redo"))
stopifnot(status==0, file.exists(tree_path))
tree_raw <- read.tree(tree_path)
stopifnot(length(tree_raw$tip.label)==359,
          !anyDuplicated(tree_raw$tip.label),
          setequal(tree_raw$tip.label, gl_main@other$meta$SampleID))
tree_midpoint <- midpoint.root(tree_raw)

# Save results: 03b_ML_tree_plot_fin.r handles SH-aLRT >=80 AND UFBoot >=95.
saveRDS(list(tree_raw=tree_raw, tree_midpoint=tree_midpoint,
             source_vcf="01_Data/vcf_260902_8625_fin/populations_final_359_8625.vcf",
             n_individuals=359, n_snp_records=8625, n_variable_alignment_sites=3047,
             model="MFP+ASC", ufboot=1000, sh_alrt=1000,
             threads="AUTO", rooting="midpoint"),
        "00_rds_modules/03_ml_tree_models.rds")
