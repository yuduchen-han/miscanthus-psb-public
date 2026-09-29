# Publication scientific pipelines

Run all scripts from the `PROS_Pipeline/` root (parent of this directory). `00_rawcode/` is the read-only provenance layer. All paths below are relative to that root, except `../../2_data/` reference annotation (relative to `PROS_Pipeline/`) and the externally supplied Stacks input path. Never run these scripts against existing formal results without directing output to a separate location: static review did **not** rerun analyses or validate numerical equivalence.

## Workflow and raw → publication lines

| Raw | Publication script | Lines | Preserved scientific core |
| --- | --- | ---: | --- |
| `01_core_computation.R` | `01_data_preparation/01_core_computation.R` | 190 → 42 | 359 × 8625 VCF, PLINK suffix cleanup, explicit metadata alignment, Region_L2, downstream region/macro/habitat/variety palettes. |
| `01_pca_habitat.r` | `02_population_structure/01_pca_habitat.r` | 158 → 46 | glPca nf=3, parallel=TRUE/n.cores=4, PC variance, Habitat color/variety shape, manuscript PDF/PNG/CSV. |
| `04_structure.r` | `02_population_structure/04_structure.r` | 928 → 85 | Genlight→LEA .geno (missing=9), K=1:8, entropy=TRUE, 1000 runs/K, project="force", all cross-entropies, minimum CE run for K=2:5, barplots/data; K=3 minimum from established analysis. |
| `04b_Structure_Map_K3.R` | `02_population_structure/04b_Structure_Map_K3.R` | 765 → 107 | Existing project; 1000 K=3 runs, minimum CE Q, Region_L2×Habitat_new means and centroids, all original manual pie/label offsets, true coordinates. |
| `03b_PoMo_build_counts.py` | `02_population_structure/03b_PoMo_build_counts.py` | 132 → 60 | 359 VCF samples→46 Pop_L3, every all-sites record, missing GT skipped, A/C/G/T allele counts, NSITES=751805. |
| `03b_PoMo_iqtree.sh` | `02_population_structure/03b_PoMo_iqtree.sh` | 29 → 27 | GTR+P, 1000 UFBoot, nt AUTO; PoMo original defaults. |
| `03b_PoMo_plot.r` | `02_population_structure/03b_PoMo_plot.r` | 111 → 34 | Midpoint-rooted 46-tip tree, Region_L2 tip mapping/colors; no bootstrap-label addition. |
| `03_MLtree.r` | `02_population_structure/03_MLtree.r` | 1459 → 76 | Matched 359×8625 VCF/Genlight, vcf2phylip, IQ-TREE-generated 3047 variable sites, MFP+ASC, 1000 UFBoot/SH-aLRT, AUTO, midpoint rooting. Uses PATH `iqtree2`. |
| `03b_ML_tree_plot_fin.r` | `02_population_structure/03b_ML_tree_plot_fin.r` | 177 → 53 | Final fan plot; node support SH-aLRT≥80 **and** UFBoot≥95, original colors/shapes. |
| `06_IBD_sampsites.r` | `03_genetic_differentiation/06_IBD_sampsites.r` | 1444 → 131 | Pop_L3 n≥3, WC FST, negative FST zeroed only for distance, Haversine centroid, Pearson Mantel 9999, seed 2026. |
| `06b_IBD_plot_din.r` | `03_genetic_differentiation/06b_IBD_plot_din.r` | 126 → 41 | Figure 6 data/Mantel from IBD RDS, pair categories and original final plotting style. |
| `07_AMOVA_habitatmew.r` | `03_genetic_differentiation/07_AMOVA_habitatmew.r` | 349 → 69 | 359×8625, Habitat/Pop_L3, poppr.amova clonecorrect=FALSE/within=FALSE, randtest 999, seed 2026, original signed Sigma and percentage/Table 2. |
| `09_outlier_habitat.r` | `04_outlier_analysis/09_outlier_habitat.r` | 868 → 56 | pcadapt K=20 diagnostic/K=3 formal, min.maf=.05, BH q<.01, robust MAF≥.10. |
| `10_GEA.r` | `04_outlier_analysis/10_GEA.r` | 1184 → 127 | SNP mean imputation; environmental NA exclusion/VIF=10/scale; RDA permutation 999/seed 2026/parallel 1; RDA1 **or** RDA2 >3 SD. |
| `14_Consensus_Intersection_Manhattan.R` | `04_outlier_analysis/14_Consensus_Intersection_Manhattan.R` | 1747 → 145 | Strict pcadapt robust∩RF Habitat∩RDA, 8625-SNP universe, 21 consensus; GFF overlap/nearest/no_gene_annotation; functional join cannot change membership. |
| `14c_21snp_gst_mmod.r` | `04_outlier_analysis/14c_21snp_gst_mmod.r` | 466 → 59 | 359×8625, Species grouping, Genlight→Genind, one mmod::Gst_Nei genome-wide run, 21 SNP extraction, mean/median/90th/95th, annotation and manuscript table. |
| `08c_RF_Variety_21snp.r` | `04_outlier_analysis/08c_RF_Variety_21snp.r` | 789 → 70 | Fixed 21 SNPs (no selection), SNP mean imputation, Species, ranger 5000 trees×3, mtry=round(21/3)=7, permutation importance, seed=500000+r, threads=8, diagnostic PCA. |
| `08d_RF_Variety_21snp_plot_fin.r` | `04_outlier_analysis/08d_RF_Variety_21snp_plot_fin.r` | 152 → 40 | Final Figure 5: Habitat color, variety shape, PCA variance labels, PDF/PNG. |
| `13_4_0_make_table1_popmap.R` | `06_tables_figures/01_make_table1_popmap.R` | 187 → 32 | Final 359 SampleIDs; Region_L2 14 groups; Stacks two-column popmap. |
| `13_2_table1.sh` | `06_tables_figures/02_table1_stacks.sh` | 18 → 12 | Stacks -t 8 -R .80 --min-mac 3 --max-obs-het .60 --write-single-snp; original Stacks input supplied via STACKS_INPUT_DIR. |
| `13_4_table1_fin.r` | `06_tables_figures/03_table1_summary.R` | 505 → 91 | Ho/He/pi/FIS from Stacks **All positions (variant and fixed)**; AR from final **8625 LD-pruned SNP** Genlight via hierfstat::allelic.richness; final CSV. |
| `16a_upplementary_Population_Mapping.r` | `06_tables_figures/16a_upplementary_Population_Mapping.r` | 190 → 171 | User-confirmed sampling table source: 359 individuals, 46 Pop_L3, 14 Region_L2, population labels and true coordinates. **Not further reduced:** raw source does not produce the requested sampling figure. |

Typical order: core → PCA/sNMF/K3 map/individual ML tree/PoMo three-step → IBD/AMOVA/pcadapt/RDA → RF Habitat (separate pending producer) → strict consensus → 21-SNP GST/RF/final Figure 5 → Table 1 three-step/sampling information. Figure 6 uses the IBD RDS; the individual ML final plot uses the midpoint tree RDS; the 21-SNP analyses use the consensus RDS.

Removed engineering code: loggers, duplicated QC and console reports, temporary audit tables, unused FASTA/tree preview branches, old K4/K5 maps, superseded diagnostic figures, plot-object caches, Word preview, and redundant archive RDS. Final figures, data tables and RDS needed by downstream publication scripts remain. Nontrivial annotation/gene membership checks remain because deleting them could silently change the candidate set.

## Dependencies / unresolved work

- **RF Habitat:** `08_RF_habitat.r` remains untouched in `00_rawcode/`; Stage 1 percentage screening lacks a confirmed executable producer. Consensus reads the byte-identical archived handoff `03_Public_Code/data/handoff/RF_Habitat_final_candidates.csv` (534 SNP IDs). Do not claim end-to-end RF reproduction until Stage 1 is resolved.
- **Sampling figure:** confirmed `16a_upplementary_Population_Mapping.r` only writes a supplementary table. No code in this source generates the requested final map figure/manual positioning. Its publication copy was left unchanged beyond round 1; figure producer is **NEEDS VERIFICATION**, not substituted with `16_sampling_map.r`.
- **Table 1:** `02_table1_stacks.sh` directly reads the authoritative `03_Public_Code/00_upstream/final_359_samples.tsv`; set `STACKS_INPUT_DIR` to the original Stacks catalog before running. `01_make_table1_popmap.R` is only convenience regeneration from Genlight, **not historical popmap provenance**; it need not be run first. `populations` must be on PATH. The manuscript CSV does not conflate Stacks all-sites and AR LD-pruned universes. Publication population metadata now uses final Pop_L3 `Species` classifications; its separate `Table1_Variety` preserves the existing final Table 1's 14-region display labels without changing statistics or grouping (see `data/README.md`).
- **Tree dependencies:** PoMo requires final all-sites VCF, metadata and its original bundled IQ-TREE binary; `NSITES=751805` is preserved, not recounted. Individual ML requires external `01_Data/vcf2phylip.py`, `python3`, and `iqtree2` on PATH. Static review does not verify IQ-TREE/PoMo model compatibility or number of variable sites by executing the workflow.
- **Small data:** `03_Public_Code/data/metadata/` contains final-359 individual and final-46 population tables; `data/README.md` documents sources, location-release review, handoff checksum and external-data manifest. Publication scripts now read these small files from `PROS_Pipeline/` root; large VCF paths remain unchanged until archive download locations are confirmed.
- **Annotation:** reference GFF and functional table retain the source-defined paths relative to `PROS_Pipeline/` under `../../2_data/referece_datasat/`; their provenance/release copies need confirmation. Missing functional annotation cannot alter the 21 candidate SNPs, but the downstream GST manuscript table requires its expected annotation fields.
- Existing empty placeholder scripts are not confirmed publication workflows and were left untouched. Existing archived results were not overwritten; parse/syntax checks are not numerical-equivalence tests.

Static checks: all 32 publication R files `parse()` (19 nonempty scripts plus 13 empty placeholders); all three shell files `bash -n` (one empty placeholder); the PoMo Python script `ast.parse()`. No statistical analysis or official output regeneration was run.
