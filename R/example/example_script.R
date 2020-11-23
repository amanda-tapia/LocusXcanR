

# directory
dsdir <- "C:/Users/manda/OneDrive/Courses - UNC/Dissertation/Kaiser/RApp/"

# example header details
headernote <- "All genomic positions are from GRCh37.<p>TWAS results for 10 traits are presented in the figures and tables below. Traits and trait categories are listed and defined as follows:<p><ul><li>Platelet count (PLT)</li><li>Red blood cell indices (RBC): red blood cell count (RBC), hematocrit (HCT), hemoglobin (HGB), mean corpuscular volume (MCV), and red cell distribution width (RDW)</li><li>White blood cell indices (WBC): white blood cell count (WBC), monocyte (MONO), neutrophil (NEUTRO), and lymphocyte (LYMPH)</li></ul>"

# example methods
method <- source("R/example/methods.R")

# call the function
LocusXcanR(
         #twas_result = paste0(dsdir,"KaiserAnalysisDS.txt"),
         twas_result = "R/example/twas_ds.txt",
         
         #known_variants = paste0(dsdir,"RBC_knownsnp.match.txt"),
         known_variants = "R/example/gwas_sentinel.txt",
         
         #weight_tbl = paste0(dsdir,"DGN-weights.txt"),
         weight_tbl = "R/example/weight_tbl.txt",
         
         
         known_gwas = paste0(dsdir,"gwaspval_DGN.txt"),
         
         
         #db_genes = paste0(dsdir,"DGN_genes_Kaiser.txt"),
         db_genes = "R/example/db_genes.txt",

         
         all_gwas = paste0(dsdir,"kaiser.gwas_locusALL.txt"),
         pred_exp_corr = paste0(dsdir,"DGN_expression_dec_cor.Rdata"),
         ld_gwas = paste0(dsdir,"locusLD_topsub.ld"),
         
         
         study_name = "Genetic Epidemiology Research on Adult Health and Aging (GERA) Europeans",
         ref_expr_name = "PredictDB Depression Genes and Networks (DGN) weights",
         head_details = headernote,
         method_details = method,
         conditional_present = F,
         primary_tissue = "DGN",
         meta_present = F,
         multiple_tissues = F)
