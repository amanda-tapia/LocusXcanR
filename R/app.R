#' TWAS results integrated with other 'omics studies
#' R Shiny app to pull in relevant data for TWAS results interpretation
#'
#' created by: Amanda L. Tapia    Date: 06/08/2020
#'
#' @param twas_result File path to TWAS results, as a text string (required)
#' @param weight_tbl File path to TWAS weights database, as a text string (required)
#' @param known_variants File path to known GWAS variants, as a text string (required)
#' @param pred_exp_corr File path to predicted expression correlation, as a text string (required)
#' @param study_name A text string of the name of the study (optional, missing name is default)
#' @param conditional_present TRUE if conditional analysis results are available, FALSE if not (FALSE is default)
#' @param multiple_tissues TRUE if TWAS results are available from more than one tissue, FALSE if TWAS results are only available from a single tissue or a multi-tissue analysis (FALSE is default)
#' @param known_gwas File path to study GWAS data matching known variants
#' @export
#' @importFrom magrittr "%>%"
#' @importFrom tidyr "separate"
#' @importFrom dplyr "filter","select","distinct"
#' @importFrom plotly "renderPlotly","plotlyOutput","ggplotly"
#' @importFrom visNetwork "renderVisNetwork","visNetworkOutput","visNetwork"
#' @importFrom data.table ":=","data.table","as.data.table"
#' @importFrom shiny "fluidPage","h3","h4","HTML","tabPanel","tabsetPanel","br","hr","strong"
#' @importFrom DT "formatStyle","styleEqual","datatable"
#' @importFrom ggplot2 "scale_colour_manual","ggplot","aes","geom_point","geom_hline","theme_bw","geom_segment","annotate"
#' @importFrom grDevices "X11"
# @importFrom graphics "layout"
#' @importFrom utils "read.table"
#' @param db_genes File path to a list of genes in the database
#' @param all_gwas File path to study GWAS data
#' @param ld_gwas File path to the LD among study variants or an LD reference panel
####################################################################



shinyTWAS <- function(twas_result,weight_tbl,study_name="",pred_exp_corr,conditional_present=FALSE,multiple_tissues=FALSE,
                      known_variants,known_gwas,db_genes,all_gwas,ld_gwas){
  
  # load analysis dataset
  twas_ds <- data.table::fread(twas_result,stringsAsFactors = F, header=T, sep = '\t')
  
  
  # load known variants dataset
  rbcds <- read.table(known_variants,
                     stringsAsFactors = F, header=F, sep = '\t',fill=T,col.names=1:20, comment.char = "")
  #wbcds <- read.table('WBC_knownsnp.match.txt',
  #                   stringsAsFactors = F, header=F, sep = '\t',fill=T,col.names=1:20, comment.char = "")
  #pltds <- read.table('PLT_knownsnp.match.txt',
  #                   stringsAsFactors = F, header=F, sep = '\t',fill=T,col.names=1:20, comment.char = "")
  
  # create indicator for known variants in study data set
  rbcds$X21 <- ifelse(is.na(rbcds$X19), 0, 1)
  #wbcds$X21 <- ifelse(is.na(wbcds$X19), 0, 1)
  #pltds$X21 <- ifelse(is.na(pltds$X19), 0, 1)  
  
  # study gwas data (known variants at the locus)
  gwasds <- read.table(known_gwas,
                        stringsAsFactors = F, header = T, sep = '\t')
  gwasds$log10p <- -log10(gwasds$PVAL)
  gwasds$poslog10p <- log10(gwasds$PVAL)
  gwasdsfin <- gwasds %>% separate(SNP,c("chr","pos","all1","all2"),sep=':',remove=F)

   
  # genes included in PredictDB
  genes <- read.table(db_genes,
                     stringsAsFactors = T, header = T)
  
  
  # load all Kaiser GWAS data (all GWAS with pvalues < 0.05)
  gwasall <- data.table::fread(all_gwas, header = T,stringsAsFactors = F)
  gwasallfin <- gwasall %>% separate(SNP,c('chr','pos','all1','all2'),sep=':',remove=F,extra = 'merge')
  gwasallfin$neglog10p <- -log10(gwasallfin$PVAL)
  gwasallfin$poslog10p <- (-1)*gwasallfin$neglog10p
   
   
  # load correlated predicted expression data
  load(pred_exp_corr)
  
  
  # load the weights data set
  weight_ds <- data.table::fread(weight_tbl, header = T,stringsAsFactors = F, fill=T)

  
  # load the LD data
  LDds <- read.table(ld_gwas, header=T, stringsAsFactors = F)
  
   
  # color blind color palette
  cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")




# library(shiny)
# library(tidyverse)
# library(DT)
# library(plotly)
# library(RColorBrewer)
# library(data.table)
# library(visNetwork)
# 
# study_name="Genetic Epidemiology Research on Adult Health and Aging (GERA) Europeans"
# conditional_present=1
# metaanalysis_present=1
# twas_result="KaiserAnalysisDS.txt"
# weight_tbl="DGN-weights.txt"
# multiple_tissues=1
# 
# 
# # analysis dataset
# twas_ds <- fread(twas_result,stringsAsFactors = F, header=T, sep = '\t')
# 
# 
# # known variants datasets
# rbcds <- read.table('RBC_knownsnp.match.txt',
#                     stringsAsFactors = F, header=F, sep = '\t',fill=T,col.names=1:20, comment.char = "")
# wbcds <- read.table('WBC_knownsnp.match.txt',
#                     stringsAsFactors = F, header=F, sep = '\t',fill=T,col.names=1:20, comment.char = "")
# pltds <- read.table('PLT_knownsnp.match.txt',
#                     stringsAsFactors = F, header=F, sep = '\t',fill=T,col.names=1:20, comment.char = "")
# 
# # create indicator for variants in study data set
# rbcds$X21 <- ifelse(is.na(rbcds$X19), 0, 1)
# wbcds$X21 <- ifelse(is.na(wbcds$X19), 0, 1)
# pltds$X21 <- ifelse(is.na(pltds$X19), 0, 1)
# 

# # study gwas data (known variants at the locus)
# gwasds <- read.table('gwaspval_DGN.txt',
#                    stringsAsFactors = F, header = T, sep = '\t')
# gwasds$log10p <- -log10(gwasds$PVAL)
# gwasds$poslog10p <- log10(gwasds$PVAL)
# gwasdsfin <- gwasds %>% separate(SNP,c("chr","pos","all1","all2"),sep=':',remove=F)
# 
# 
# # genes included in PredictDB
# genes <- read.table('DGN_genes_Kaiser.txt',
#                     stringsAsFactors = T, header = T)
# 
# 
# # load all Kaiser GWAS data (all GWAS with pvalues < 0.05)
# #gwasall <- fread('kaiser.gwas_0005_OLD.txt', header = T,stringsAsFactors = F)
# gwasall <- fread('kaiser.gwas_locusALL.txt', header = T,stringsAsFactors = F)
# gwasallfin <- gwasall %>% separate(SNP,c('chr','pos','all1','all2'),sep=':',remove=F,extra = 'merge')
# gwasallfin$neglog10p <- -log10(gwasallfin$PVAL)
# gwasallfin$poslog10p <- (-1)*gwasallfin$neglog10p
# 
# 
# # load correlated predicted expression data
# #expds <- fread("DGN_expression_dec.txt", header=T, stringsAsFactors = F)
# load("DGN_expression_dec_cor.Rdata")
# 
# 
# # load the weights data set
# weight_ds <- fread(weight_tbl, header = T,stringsAsFactors = F, fill=T)
# 
# # load the LD data
# LDds <- read.table("locusLD_topsub.ld", header=T, stringsAsFactors = F)
# 
# # color blind color palette
# cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# 
# 
# 
# ####################################################################
# 
#
  # filter analysis dataset by tissue
  dgnds <- twas_ds %>% filter(tissue=='DGN')
  
  dgnds$locvar <- paste0("Locus ",dgnds$locus2," : ", "chr ",dgnds$chr," : ",dgnds$locstart,", ",
                         dgnds$locstop," : ",dgnds$index," : ",dgnds$pheno
                        )
  dgnds$genestartMB <- round(dgnds$genestart/1000000,4)
  dgnds$genestopMB <- round(dgnds$genestop/1000000,4)
  dgnds$log10pvalmeta <- -log10(dgnds$P_value_meta)
  
  
  sorteddgn <- dgnds[
    order( dgnds$locus2, dgnds$genestart, dgnds$genestop, dgnds$pheno ),
  ]
  
  # list of unique loci
  loclist <- sort(unique(dgnds$locus2[complete.cases(dgnds$locus2)]))
  loclistdet <- unique(sorteddgn$locvar[complete.cases(sorteddgn$locus2)])
  
  
  # p-value threshold for primary reference panel
  pthresh <- -log10(0.05/(nrow(dgnds)))
  
  # number of significant TWAS genes
  signifrow <- dgnds %>% filter(SignifGene==1 & is.na(HLARegion) & is.na(MHCRegion)
                               & SingleSNP!=1)
  numsignif <- nrow(signifrow)
  
  # meta-analysis threshold
  metathresh <- log10(0.05/numsignif)




############################################################################################
############################################################################################


  # # set up UI
  # ui <- shiny::fluidPage(
  #   h3("Transcriptome-wide association study (TWAS)"),
  #   shiny::h4(study_name),
  #   shiny::br(),
  #   shiny::br(),
  # 
  # 
  #   shiny::br(),
  #   shiny::fixedPanel(style="z-index:100;",
  #              top = 125, left = 10, right = 0, width = "50%", draggable = T,
  #              shiny::selectInput("locuslst","Select a locus to view (click and drag to reposition menu):",
  #                          choices = loclistdet, width = "550px")
  #              ),
  #   shiny::br(),
  #   shiny::br(),
  # 
  # 
  #   shiny::h4(shiny::strong("Figure 1. TWAS-GWAS mirror plot of genes and variants within the locus")),
  # 
  #   "Note: Top figure displays TWAS significant genes and any additional non-significant genes reported from GWAS, bottom figure displays GWAS variants. In the TWAS plot, \"reported in GWAS\" means that the GERA TWAS gene was reported in the GWAS catalog as the assigned gene for a single variant signal associated with the phenotype category, often based on physical proximity. In the GWAS plot, \"reported in GWAS\" means that the GERA GWAS variant was reported in the GWAS catalog as a single variant signal associated with the phenotype category. Marginal TWAS displays results of gene-trait associations. Conditional TWAS displays results of gene-trait associations, conditional on reported GWAS variants at the locus (conditional results only available for significant TWAS genes).",
  #   shiny::radioButtons("margcondradio", label = shiny::h5(shiny::strong("Select TWAS results to display:")),
  #                choices = list("Marginal" = 1, "Conditional" = 2),
  #                selected = 1, inline=T),
  #   plotly::plotlyOutput("TWASmirror", height = 550),
  #   shiny::br(),
  #   shiny::hr(),
  # 
  # 
  #   shiny::h4(shiny::strong("Figure 2. Predicted expression correlation between index TWAS gene and other TWAS significant and reported GWAS genes at the locus and their relationship to GWAS variants included in predicted gene expression model")),
  #   "Note: Color scale for genes denotes the degree of predicted expression correlation with the index gene",
  #   shiny::fluidRow(
  #     shiny::column(6,
  #            plotly::plotlyOutput("TWAScorr", height=600),
  #     ),
  #     shiny::column(6,
  #       visNetwork::visNetworkOutput("TWASnetwork", height=600),
  #     ),
  #   ),
  #   shiny::br(),
  #   shiny::hr(),
  # 
  # 
  #   shiny::h4(shiny::strong("Figure 3. Mirror plot of GERA TWAS and meta-analysis TWAS")),
  #   shiny::h5("Note: "),
  #   plotly::plotlyOutput("Tmetamirror", height=600),
  #   shiny::br(),
  #   shiny::hr(),
  # 
  #   shiny::h4(shiny::strong("Figure 4. Comparison of TWAS results from DGN reference panel to results from secondary reference panels")),
  #   shiny::h5("Note: DGN = Depression Genes and Networks, GWB = GTEx whole blood, GTL = GTEx EBV transformed lymphocytes, MSA = MESA monocytes; each represents a gene expression reference panel"),
  #   shiny::h5("The figure in each tab displays a mirror plot of GERA results using DGN reference panel versus GERA results using a secondary reference panel (GWB, GTL, or MSA)."),
  #   shiny::tabsetPanel(
  #     id = "twascompare",
  #     shiny::tabPanel("DGN vs. GWB",plotly::plotlyOutput("CompRefGWB", height=600)),
  #     shiny::tabPanel("DGN vs. GTL",plotly::plotlyOutput("CompRefGTL", height=600)),
  #     shiny::tabPanel("DGN vs. MSA",plotly::plotlyOutput("CompRefMSA", height=600))
  #   ),
  #   shiny::br(),
  #   shiny::hr(),
  # 
  # 
  #   shiny::h4(shiny::strong("Table 1. Overall TWAS results from primary and secondary reference panels within the locus")),
  #   shiny::h5("Note: DGN = Depression Genes and Networks, GWB = GTEx whole blood, GTL = GTEx EBV transformed lymphocytes, MSA = MESA monocytes; each represents a gene expression reference panel"),
  #   shiny::h5(shiny::HTML('<span style="background-color:lightgreen">Significant gene-trait associations highlighted in green, </span><span style="background-color:tomato">HLA genes / MHC regions highlighted in red</span>')),
  #   shiny::tabsetPanel(
  #     id = 'twasresult',
  #     shiny::tabPanel("DGN", DT::dataTableOutput("DGNtbl")),
  #     shiny::tabPanel("GWB", DT::dataTableOutput("GWBtbl")),
  #     shiny::tabPanel("GTL", DT::dataTableOutput("GTLtbl")),
  #     shiny::tabPanel("MSA", DT::dataTableOutput("MSAtbl"))
  #   ),
  #   shiny::br(),
  #   shiny::hr(),
  # 
  #   shiny::h4(shiny::strong("Table 2. Overall TWAS results for other traits in the trait category from primary and secondary reference panels within the locus")),
  #   shiny::h5("Note: DGN = Depression Genes and Networks, GWB = GTEx whole blood, GTL = GTEx EBV transformed lymphocytes, MSA = MESA monocytes; each represents a gene expression reference panel"),
  #   shiny::h5(shiny::HTML('<span style="background-color:lightgreen">Significant gene-trait associations highlighted in green, </span><span style="background-color:tomato">HLA genes / MHC regions highlighted in red</span>')),
  #   shiny::tabsetPanel(
  #     id = 'twasresulttrt',
  #     shiny::tabPanel("DGN", DT::dataTableOutput("DGNtbltrt")),
  #     shiny::tabPanel("GWB", DT::dataTableOutput("GWBtbltrt")),
  #     shiny::tabPanel("GTL", DT::dataTableOutput("GTLtbltrt")),
  #     shiny::tabPanel("MSA", DT::dataTableOutput("MSAtbltrt"))
  #   ),
  #   shiny::br(),
  #   shiny::hr(),
  # 
  # 
  #   shiny::h4(shiny::strong("Table 3. Reported GWAS variants within the locus (all traits in the category)")),
  #   shiny::h5(shiny::HTML('<span style="background-color:tomato">Table <strong>row</strong> is highlighted in red if the SNP is not in GERA imputation data</span>')),
  #   shiny::h5(shiny::HTML('<span style="background-color:lightgreen"><strong>RSIDs</strong> reported in Figure 1 highlighted in green</span>')),
  #   shiny::h5("GWAS ",shiny::strong("Genes"),":"),
  #   shiny::p(shiny::HTML('<ul>
  #     <li><span style="background-color:lightgreen">matching significant TWAS genes highlighted in green</span></li>
  #     <li><span style="background-color:yellow"    >matching non-significant TWAS genes highlighted in yellow</span></li>
  #     <li><span style="background-color:orange"    >included in DGN but not predicted in GERA highlighted in orange</span></li>
  #     <li><span style="background-color:tomato"    >included in DGN but predictive model does not contain SNPs highlighted in red</span></li>
  #     <li>not included in DGN (i.e. we have no info in this analysis) are not highlighted</li>
  #     </ul>
  #     ')),
  #   shiny::br(),
  #   DT::dataTableOutput("KnownSNPtbl"),
  #   shiny::br(),
  #   shiny:: hr(),
  # 
  # 
  #   shiny::h4(shiny::strong("Table 4. GERA GWAS results displayed in Figure 1 as variants \"Reported in GWAS\"")),
  #   shiny::h5("Note: Results listed below are from trait specific GWAS"),
  #   DT::dataTableOutput("GWASvars"),
  #   shiny::br()
  # )
  
  
  # set up UI
  ui <- 
    navbarPage(title="Shiny TWAS",
               tabPanel(title="Results",
                        fluidPage(
                          h4("Transcriptome-wide association study (TWAS) using PredictDB Depression Genes and Networks (DGN) weights"),
                          h5("Genetic Epidemiology Research on Adult Health and Aging (GERA) Europeans"),
                          "All genomic positions are from GRCh37. ",
                          "TWAS results for 10 traits are presented in the figures and tables below. Traits and trait categories are listed and defined as follows:",
                          HTML("<ul><li>Platelet count (PLT)</li><li>Red blood cell indices (RBC): red blood cell count (RBC), hematocrit (HCT), hemoglobin (HGB), mean corpuscular volume (MCV), and red cell distribution width (RDW)</li><li>White blood cell indices (WBC): white blood cell count (WBC), monocyte (MONO), neutrophil (NEUTRO), and lymphocyte (LYMPH)</li></ul>"),
                          br(),
                          br(),
                          
                          # selectInput("locuslst","Select a locus to view:",choices = loclist),
                          # br(),
                          
                          # selectInput("loclstdet","Select a locus to view: ", choices = loclistdet, width="500px"),
                          # br(),
                          
                          br(),
                          fixedPanel(style="z-index:100;",
                                     top = 200, left = 10, right = 0, width = "60%", draggable = T,
                                     selectInput("locuslst","Select a locus to view (click and drag to reposition menu):",
                                                 choices = loclistdet,width = "550px")
                          ),
                          br(),
                          br(),
                          br(),
                          
                          # h4("Figure 1. GERA TWAS gene-trait associations within the locus"),
                          # h5("Note: Figure displays TWAS significant genes and any additional TWAS genes known from GWAS"),
                          # plotlyOutput("TWASplt"),
                          # br(),
                          # 
                          # h4("Figure 2. GERA GWAS known variants within the locus"),
                          # plotlyOutput("GWASplt"),
                          # br(),
                          
                          # uncomment when memory available
                          # plotOutput("chrplt",height = 100),
                          # br(),
                          
                          h4(strong("Figure 1. TWAS-GWAS mirror plot of genes and variants within the locus")),
                          
                          "Note: Top figure displays TWAS significant genes and any additional non-significant genes reported from GWAS, bottom figure displays GWAS variants. In the TWAS plot, \"reported in GWAS\" means that the GERA TWAS gene was reported in the GWAS catalog as the assigned gene for a single variant signal associated with the phenotype category, often based on physical proximity. In the GWAS plot, \"reported in GWAS\" means that the GERA GWAS variant was reported in the GWAS catalog as a single variant signal associated with the phenotype category. Marginal TWAS displays results of gene-trait associations. Conditional TWAS displays results of gene-trait associations, conditional on reported GWAS variants at the locus (conditional results only available for significant TWAS genes).",
                          radioButtons("margcondradio", label = h5(strong("Select TWAS results to display:")),
                                       choices = list("Marginal" = 1, "Conditional" = 2), 
                                       selected = 1, inline=T),  
                          plotlyOutput("TWASmirror", height = 550),
                          br(),
                          hr(),
                          
                          
                          h4(strong("Figure 2a. TWAS-GWAS mirror locus-zoom plot")),
                          "Note: Top panel displays predicted expression correlation between index TWAS gene and other genes at the locus. Bottom panel displays LD between the index SNP and other SNPs at the locus. Lines connect genes to their predictive model variants. Color scale for genes denotes the degree of predicted expression correlation with the index gene. Color scale for SNPs and solid lines denotes the degree of LD with the index SNP. Dashed red line in the top panel denotes TWAS p-value threshold = 4.37e-7, and in bottom panel denotes GWAS p-value threshold = 5.0e-8",
                          br(),
                          # fluidRow(
                          #   column(6,
                          #     plotlyOutput("TWAScorr", height=600),
                          #   ),
                          #   column(6,
                          #     visNetworkOutput("TWASnetwork", height=600),
                          #   ),
                          # ),
                          plotlyOutput("TWAScorr", height=550),
                          br(),
                          h4(strong("Figure 2b. Network visualization of TWAS results")),
                          "Sentinel TWAS gene is indicated by a star, all other genes are squares. Sentinel GWAS variant is indicated by a triangle, all other variants are circles. Color scale of all lines and shapes is based on correlation with the index gene or index variant. Line thickness corresponds to the model weight. Solid line indicates a positive direction of effect and dashed line indicates a negative direction. Size of the shape corresponds to the size of the -log10(p-value).",
                          br(),
                          visNetworkOutput("TWASnetwork", height=600),
                          hr(),
                          
                          
                          h4(strong("Figure 3. Mirror plot of GERA TWAS and meta-analysis TWAS")),
                          "Note: Top figure displays TWAS significant genes and any additional non-significant genes reported from GWAS, bottom figure displays the same genes from TWAS meta-analysis of ARIC, WHI, and BioMe. In both plots, \"reported in GWAS\" means that the TWAS gene was reported in the GWAS catalog as the assigned gene for a single variant signal associated with the phenotype category, often based on physical proximity.",
                          plotlyOutput("Tmetamirror", height=600),
                          br(),
                          hr(),
                          
                          h4(strong("Figure 4. Comparison of TWAS results from DGN reference panel to results from secondary reference panels")),
                          "Note: DGN = Depression Genes and Networks, GWB = GTEx whole blood, GTL = GTEx EBV transformed lymphocytes, MSA = MESA monocytes; each represents a gene expression reference panel. The figure in each tab displays a mirror plot of GERA results using DGN reference panel versus GERA results using a secondary reference panel (GWB, GTL, or MSA).",
                          #plotlyOutput("CompRef", height = 600),
                          tabsetPanel(
                            id = "twascompare",
                            tabPanel("DGN vs. GWB",plotlyOutput("CompRefGWB", height=600)),
                            tabPanel("DGN vs. GTL",plotlyOutput("CompRefGTL", height=600)),
                            tabPanel("DGN vs. MSA",plotlyOutput("CompRefMSA", height=600))
                          ),
                          br(),
                          hr(),
                          
                          # h4("Table 1. Overall GERA TWAS results within the locus"),
                          # h5("Note: Significant gene-trait associations highlighted in green, HLA genes / MHC regions highlighted in red"),
                          # DT::dataTableOutput("TWAStbl"),
                          # br(),
                          
                          h4(strong("Table 1. Overall TWAS results from primary and secondary reference panels within the locus")),
                          "Note: DGN = Depression Genes and Networks, GWB = GTEx whole blood, GTL = GTEx EBV transformed lymphocytes, MSA = MESA monocytes; each represents a gene expression reference panel. ",
                          HTML('<span style="background-color:lightgreen"> Significant gene-trait associations highlighted in green. </span> <span style="background-color:tomato"> HLA genes / MHC regions / single SNP models highlighted in red. </span>'),
                          " MHC region is defined as GRCh37; chr6:28,477,797-33,448,354. Single SNP model indicates that the predictive expression model for the gene contained only a single SNP.",
                          tabsetPanel(
                            id = 'twasresult',
                            tabPanel("DGN", DT::dataTableOutput("DGNtbl")),
                            tabPanel("GWB", DT::dataTableOutput("GWBtbl")),
                            tabPanel("GTL", DT::dataTableOutput("GTLtbl")),
                            tabPanel("MSA", DT::dataTableOutput("MSAtbl"))
                          ),
                          br(),
                          hr(),
                          
                          h4(strong("Table 2. Overall TWAS results for other traits in the trait category from primary and secondary reference panels within the locus")),
                          "Note: DGN = Depression Genes and Networks, GWB = GTEx whole blood, GTL = GTEx EBV transformed lymphocytes, MSA = MESA monocytes; each represents a gene expression reference panel. ",
                          HTML('<span style="background-color:lightgreen">Significant gene-trait associations highlighted in green, </span><span style="background-color:tomato">HLA genes / MHC regions / single SNP models highlighted in red. </span>'),
                          " MHC region is defined as GRCh37; chr6:28,477,797-33,448,354. Single SNP model indicates that the predictive expression model for the gene contained only a single SNP.",
                          tabsetPanel(
                            id = 'twasresulttrt',
                            tabPanel("DGN", DT::dataTableOutput("DGNtbltrt")),
                            tabPanel("GWB", DT::dataTableOutput("GWBtbltrt")),
                            tabPanel("GTL", DT::dataTableOutput("GTLtbltrt")),
                            tabPanel("MSA", DT::dataTableOutput("MSAtbltrt"))
                          ),
                          br(),
                          hr(),
                          
                          
                          h4(strong("Table 3. Reported GWAS variants within the locus (all traits in the category)")),
                          h5(HTML('<span style="background-color:tomato">Table <strong>row</strong> is highlighted in red if the SNP is not in GERA imputation data</span>')),
                          h5(HTML('<span style="background-color:lightgreen"><strong>RSIDs</strong> reported in Figure 1 highlighted in green</span>')),
                          h5("GWAS ",strong("Genes"),":"),
                          p(HTML('<ul>
        <li><span style="background-color:lightgreen">matching significant TWAS genes highlighted in green</span></li>
        <li><span style="background-color:yellow"    >matching non-significant TWAS genes highlighted in yellow</span></li>
        <li><span style="background-color:orange"    >included in DGN but not predicted in GERA highlighted in orange</span></li>
        <li><span style="background-color:tomato"    >included in DGN but predictive model does not contain SNPs highlighted in red</span></li>
        <li>not included in DGN (i.e. we have no info in this analysis) are not highlighted</li>
        </ul>
        ')),
                          br(),
                          DT::dataTableOutput("KnownSNPtbl"),
                          br(),
                          hr(),
                          
                          
                          # h5(HTML('<span style="background-color:lightgreen">&nbsp;&nbsp; - matching significant TWAS genes highlighted in green</span>')),
                          # h5("   - matching non-significant TWAS genes highlighted in yellow",
                          #    style="color: yellow;"),
                          # h5(HTML('<span style="background-color:yellow">Yellow Text with header</span>')),
                          # p("yellow text with p()", style="background-color:yellow"),
                          # h5("   - included in DGN but not predicted in GERA highlighted in orange",
                          #    style="color: orange;"),
                          # h5("   - included in DGN but predictive model does not contain SNPs highlighted in red",
                          #    style="color: red;"),
                          # h5("   - not included in DGN (i.e. we have no info in this analysis) are not highlighted"),
                          
                          
                          
                          h4(strong("Table 4. GERA GWAS results displayed in Figure 1 as variants \"Reported in GWAS\"")),
                          h5("Note: Results listed below are from trait specific GWAS"),
                          DT::dataTableOutput("GWASvars"),
                          br()
                        )
               ),
               
               #########################################################################################    
               
               tabPanel("Methods",
                        h3(strong("Methods")),
                        "A current draft of the manuscript describing these results and detailing the methods can be found here:",
                        tags$a(href="https://docs.google.com/document/d/1TfqzBPGhL6TEqUuxN5lpy7cJVZcZhyiy/edit", "https://docs.google.com/document/d/1TfqzBPGhL6TEqUuxN5lpy7cJVZcZhyiy/edit",target="_blank"),
                        ". All cohorts included in the analysis are described individually below. We analyze 10 hematological phenotypes (platelet count, red blood cell count, hematocrit, hemoglobin, mean corpuscular volume, red cell distribution width, white blood cell count, monocyte count, neutrophil count, and lymphocyte count) across all cohorts.",
                        
                        hr(),
                        h4(strong("PrediXcan method:")),
                        "PrediXcan (26258848) is a gene-based association test that prioritizes genes which are likely to be causal for the phenotype. It implements an elastic net based method for selecting variants associated with gene expression in a given reference panel, and then uses those variants to predict gene expression in a cohort with only genotype data. We downloaded the PrediXcan software (see URLs) along with its prepackaged weights for gene expression data from PredictDB (see URLs). Weights for gene expression using RNA sequencing data were obtained from the Genotype-Tissue Expression project (23715323) (whole blood, genes= ; and EBV transformed lymphocytes, genes=; n=), Depression Genes and Networks (24092820) (whole blood, genes=11538, n=922), and Multi-Ethnic Study of Atherosclerosis (Europeans only, monocytes, genes=, n=) (23900078).  Imputed genotypes for all cohorts were filtered for imputation quality based on R2 > 0.3; variants not meeting this threshold were excluded from the analysis. We use DGN as our primary reference panel for all TWAS analyses as it is the largest single whole blood RNA-seq dataset.",
                        
                        hr(),
                        h4(strong("Included cohorts:")),
                        "These TWAS analyses were limited to self-reported white or European ancestry participants, for easy comparability with the DGN European ancestry eQTL panel, including input of LD information into the R Shiny application (see R Shiny Methods), and with the largest single-ancestry blood cell trait GWAS.",  
                        
                        h5(strong("Genetic Epidemiology Research on Adult Health and Aging (GERA)")),
                        "The GERA cohort includes over 100,000 adults who are members of the Kaiser Permanente Medical Care Plan, Northern California Region (KPNC) and consented to research on the genetic and environmental factors that affect health and disease, linking together clinical data from electronic health records, survey data on demographic and behavioral factors, and environmental data with genetic data.  The GERA cohort was formed by including all self-reported racial and ethnic minority participants with saliva samples (19%); the remaining participants were drawn sequentially and randomly from non-Hispanic White participants (81%). Genotyping was completed as previously described (26092718) using 4 different custom Affymetrix Axiom arrays with ethnic-specific content to increase genomic coverage. Genotype data were imputed to 1000 Genomes Phase 3. Principal components analysis was used to characterize genetic structure in this European sample (26092716). Hematological measures were extracted from medical records. In individuals with multiple measurements, the first visit with complete white blood cell differential (if any) was used for each participant. Otherwise, the first visit was used. In total, 54,542 participants with hematological measures were included in the analysis.",
                        
                        h5(strong("Women's Health Initiative (WHI)")),
                        "WHI originally enrolled 161,808 women aged 50-79 between 1993 and 1998 at 40 centers across the US, including both a clinical trial (including three trials for hormone therapy, dietary modification, and calcium/vitamin D) and an observational study arm (9492970). WHI recruited a socio-demographically diverse population representative of US women in this age range. Two WHI extension studies conducted additional follow-up on consenting women from 2005-2010 and 2010-2015. Genotyping was available on some WHI participants through the WHI SNP Health Association Resource (SHARe) resource, which used the Affymetrix 6.0 array (~906,600 SNPs, 946,000 copy number variation probes) and on other participants through the MEGA array (https://www.biorxiv.org/content/10.1101/188094v2). Imputation and association analysis was performed separately in individuals with Affymetrix only, MEGA only, and both Affymetrix and MEGA data. For variants with both Affymetrix and MEGA genotypes available, MEGA genotypes were used. In total, 18,100 self-reported white women with hematological phenotypes were included. Six sub-cohorts from the WHI study were included in the meta-analysis and phenotypes were not collected uniformly across the cohorts. Sample size information for each phenotype is contained in (Supplementary Table 5).",
                        
                        h5(strong("Atherosclerosis Risk in Communities (ARIC)")),
                        "The ARIC study was initiated in 1987 and recruited participants age 45-64 years from 4 field centers (Forsyth County, NC; Jackson, MS; northwestern suburbs of Minneapolis, MN; Washington County, MD) in order to study cardiovascular disease and its risk factors (2646917), including the participants of self-reported European ancestry included here. Standardized physical examinations and interviewer-administered questionnaires were conducted at baseline (1987-89), three triennial follow-up examinations, a fifth examination in 2011-13, and a sixth exam in 2016-2017. Genotyping was performed through the CARe consortium Affymetrix 6.0 array (20400780). ARIC EA genotype data were imputed to Haplotype Reference Consortium (HRC). In total, 9,345 European ancestry participants with hematological phenotypes were included in the analysis. All phenotypes were adjusted for study site, age, age squared, sex, and top ten PCs and were inverse normalized.",
                        
                        h5(strong("BioMe")),
                        "Details will be coming soon.",
                        
                        hr(),
                        h4(strong("Conditional analysis:")),
                        "For each statistically significant TWAS gene-trait association, the effect of predicted gene expression was conditioned on a set of previously reported GWAS sentinel variants from [32888494] meeting the following criteria: 1) the sentinel variant fell within a 1Mb region of the TWAS gene, 2) the trait with which the GWAS variant was associated matched the TWAS analytical trait or was within the same trait category as the analytical trait (platelets, red blood cell indices {hematocrit, hemoglobin, mean corpuscular volume, red blood cell count, red blood cell distribution width}, white blood cell indices {white blood cell count, neutrophils, monocytes, lymphocytes}), and 3) the GWAS variant met an imputation quality threshold of R2 > 0.3. We used a modified version of the cpgen R package (see cpgen Methods) to perform the conditional analysis, accounting for a PLINK KING-robust kinship matrix (20926424), which used only genotyped variants and excluded those with minor allele frequency less than 5% and those missing more than 1% of SNPs.",
                        
                        hr(),
                        h4(strong("Meta-analysis and replication with ARIC, WHI, BioMe:")),
                        "In order to replicate the conditionally significant gene-trait association, we tested each association via a meta-analysis of the ARIC, WHI, and BioMe cohorts. As described above, PrediXcan was used to facilitate gene expression imputation and association in each cohort separately, and the meta-analysis association test was conducted using METAL (20616382).",
                        br(),
                        br(),
                        "Replication of the GERA significant gene-trait associations was performed using meta-analyzed TWAS results from ARIC, WHI, and BioMe. Nine gene-trait associations remained statistically significant after conditional analysis; for this set of genes, we defined a Bonferroni-corrected statistically significant replication threshold at p-value < 5.56 X 10-3. For the fine-mapping analysis, statistical significance of replicated genes was qualified based on two different thresholds -- a stringent threshold Bonferroni-corrected for all 239 statistically significant TWAS gene-trait associations at p-value < 2.09 X 10-4, and a more lenient threshold at p-value < 0.05.",
                        
                        hr(),
                        h4(strong("FOCUS:")),
                        "We used the Fine-mapping Of CaUsal gene Sets (FOCUS) (30926970) software to fine-map TWAS statistics at genomic risk regions. As input, we used GWAS summary data from GERA along with eQTL weights from PredictDB Depression Genes and Networks whole blood data, and an LD reference panel from 1000 Genomes  phase 3. The software outputs a credible set of genes at each locus which can be used to explain observed genomic risk.",
                        
                        hr(),
                        h4(strong("Fine-mapping loci and locus categories:")),
                        "Fine-mapping loci refers to fine-mapping analysis of trait-specific genomic locations that contain, and are centered at, sentinel TWAS genes. That is, we take the set of trait-specific statistically significant TWAS genes, select the most significant gene in the set (the sentinel gene), and assign it to a locus along with any other statistically significant TWAS genes within a 1Mb window of the sentinel gene. We then select the next most significant TWAS gene which has not yet been assigned to a locus and continue in this fashion until all statistically significant TWAS genes have been assigned to a locus.
        We define locus categories based on whether the locus contains a single gene or multiple genes and whether the locus replicates in TWAS meta-analysis at either a lenient or strict threshold. Thus, locus categories are defined as follows: 1=single gene locus, strict replication (p < 2.09E-04); 2=single gene locus, replication (p < 0.05); 3=single gene locus, no replication; 4=multi gene locus, strict replication (p < 2.09E-04); 5=multi gene locus, replication (p < 0.05); 6=multi gene locus, no replication.",
                        
                        hr(),
                        h4(strong("R Shiny:")),
                        "We use R's convenient Shiny package to produce the web application which displays our GERA TWAS results. The IdeogramTrack (https://rdrr.io/bioc/Gviz/man/IdeogramTrack-class.html) uses Genome Reference Consortium Human Build 37 (GRCh37) and UCSC cytogentic bands from http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/. All GERA TWAS results were produced using PrediXcan as described above, all GERA GWAS results were produced using XXX, and GERA conditional analysis results were produced using cpgen as described below. Known GWAS sentinel variants were obtained from PMID: 32888494. Model weights and model variants were taken from our primary DGN reference panel from PredictDB (or secondary reference panels GWB, GTL, or MSA from PredictDB). Correlation of predicted expression among genes at the locus was calculated using R's cor() function, and LD among variants was computed using plink --r2 (https://zzz.bwh.harvard.edu/plink/ld.shtml). We used ggplot2() to produce all figures, except the network visualization used visNetwork(). Tables were produced using the DT package (https://www.rdocumentation.org/packages/DT/versions/0.16).",
                        
                        hr(),
                        h4(strong("cpgen:")),
                        "We used the R package cpgen to perform conditional analysis of TWAS-significant genes, while accounting for a KING kinship matrix. However, cpgen is designed in such a way that it performs eigenvalue decomposition on the cohort sample for every function call. Since we had 239 TWAS-significant associations, this would have required eigenvalue decomposition on a sample of N ~ 55,000 for each of those 239 associations, a computationally burdensome calculation. Thus, we slightly modified the cpgen script. Specifically, we computed the eigenvalue decomposition on the GERA sample outside of the cpgen script (for each phenotype), and then subsequently loaded the appropriate eigenvectors and eigenvalues into the program, modifying the script so that it could take these eigenvectors and eigenvalues as input.
        ",
                        br(),
                        br(),
               )
    )
  





############################################################################################


# # give instructions to the server
# server <- function(input,output){
# 
#   # define reactive variables for all plots
#   locds <- reactive({
#     filter(dgnds, locvar==input$locuslst)
#   })
# 
#   # get the locus number from the detailed drop down menu
#   locnum <- reactive({
#     tmp <- input$locuslst
#     tmp2 <- as.numeric(strsplit(strsplit(tmp,":")[[1]][1]," ")[[1]][2])
#     tmp2
#   })
# 
#   # define locus window
#   xhigh <- reactive({ unique(locds()$locstop) })
#   xlow <- reactive({ unique(locds()$locstart) })
#   xhighMB <- reactive({ round(xhigh()/1000000,4) })
#   xlowMB <- reactive({ round(xlow()/1000000,4) })
# 
#   # locus specific characteristics
#   locchr <-   reactive({ unique(locds()$chr) })
#   locpheno <- reactive({ unique(locds()$pheno) })
#   locphcat <- reactive({ unique(locds()$phenocat) })
#   loctitle <- reactive({ paste0("chr ",locchr(),": ",xlow(),", ",xhigh(),"; trait = ",locpheno(),
#                                 ", trait category = ",locphcat())
#     })
# 
#   # set the dataset to extract known variants at the locus
#     ds <- reactive({
#       if (locphcat()=="RBC") {
#         rbcds
#       } else if (locphcat()=="WBC") {
#         wbcds
#       } else {
#         pltds
#       }
#     })
# 
#   # select all known GWAS genes at locus
#   phenotbl <- reactive({
#     tmp <- filter(ds(),(X6)>=xlow() & xhigh()>=(X6) & X5==locchr()) %>% select(X8)
#     colnames(tmp) <- c("Gene")
#     tmp
#   })
# 
#   # select all TWAS genes at locus
#   dgntbl <-  reactive({
#     tmp <- filter(dgnds,genestart>=xlow() & genestop<=xhigh() & pheno==locpheno() & chr==locchr())
# 
#     if (nrow(phenotbl())==0){
#       tmp$kngene="Not Reported in GWAS"
#     } else{
#       tmp$kngene <- ifelse(grepl(paste(unique(phenotbl()$Gene), collapse="|"),tmp$genename),
#                                 "Reported in GWAS","Not Reported in GWAS")
#     }
#     tmp
#   })
# 
# 
#   ############################################################
# 
# 
#   # TWAS/GWAS mirror plot
#   output$TWASmirror <- renderPlotly({
#     locds <- dgnds %>% filter(locvar==input$locuslst)
#     yhigh <- max(locds$log10pval)+0.25*max(locds$log10pval)
#     pthresh <- -log10(0.05/(nrow(dgnds)))
#     locchr <- unique(locds$chr)
#     locpheno <- unique(locds$pheno)
#     locphcat <- unique(locds$phenocat)
# 
#     # select known variants at locus
#     if (locphcat=="RBC") {
#       ds=rbcds
#     } else if (locphcat=="WBC") {
#       ds=wbcds
#     } else {
#       ds=pltds
#     }
# 
#     # select all known GWAS genes at locus
#     phenotbl <- ds %>% filter(X6>=xlow() & xhigh()>=X6 & X5==locchr) %>% select(X8)
#     colnames(phenotbl) <- c("Gene")
# 
#     # select all TWAS genes at locus
#     dgntbl <- dgnds %>% filter(genestart>=xlow() & genestop<=xhigh() & pheno==locpheno & chr==locchr)
#     dgntbl$kngene <- ifelse(grepl(paste(unique(phenotbl$Gene), collapse="|"),dgntbl$genename),"Reported in GWAS","Not reported in GWAS")
# 
#     # select significant and known genes to plot
#     dgntblplt <- dgntbl %>% filter(SignifGene==1 | kngene=="Reported in GWAS")
# 
#     # select the GWAS results at the locus
#     gwasloc <- gwasallfin %>% filter(Locus==locnum())
#     gwasloc$knsnp <- ifelse(grepl(paste(unique(gwasdsfin$SNP[gwasdsfin$Locus==locnum()]),collapse="|"),
#                                   gwasloc$SNP),"Reported in GWAS","Not reported in GWAS")
# 
#     # not known variants at the locus pval < 0.1
#     gwaslocnotkn <- gwasloc %>% filter(knsnp=="Not reported in GWAS")
# 
#     # known variants at the locus, regardless of pvalue
#     knowngwas <- gwasdsfin %>% filter(Locus==locnum())
#     knowngwas$knsnp <- "Reported in GWAS"
# 
#     ylow <- min(gwaslocnotkn$poslog10p,log10(5*10^(-8)),knowngwas$poslog10p)+0.15*min(gwaslocnotkn$poslog10p,log10(5*10^(-8)),knowngwas$poslog10p)
#     nudgeval <- 0.05*max(abs(yhigh))
#     nudgevalg <- 0.05*ylow
# 
#     # define specific colors for known/not known genes
#     if(length(unique(dgntblplt$kngene))==1){
#       if(unique(dgntblplt$kngene)=="Not reported in GWAS"){
#         myColors <- setNames( c('#56B4E9','#000000'),
#                               c("Not reported in GWAS","Reported in GWAS") )
#         } else if(unique(dgntblplt$kngene)=="Reported in GWAS"){
#           myColors <- setNames( c('#000000','#56B4E9'),
#                                 c("Reported in GWAS","Not reported in GWAS") )
#         }
#       }
#     else {
#       myColors <- setNames( c('#56B4E9', '#000000', '#56B4E9'),
#                             c("Not reported in GWAS","Reported in GWAS",NA))
#     }
# 
#     colScale <- scale_colour_manual(values = myColors)
# 
#     # color scale for known/not known gwas variants
#     myColors <-
#       setNames( c('#56B4E9', '#000000', '#E69F00')
#                 , c("Reported in GWAS","Not reported in GWAS",NA))
#     colScaleg <- scale_colour_manual(values=myColors)
# 
#     if (input$margcondradio==1){
#       pvalplt <- dgntblplt$log10pval
#     }
#     else if(input$margcondradio==2){
#       pvalplt <- -log10(dgntblplt$p_final)
#     }
# 
#     # TWAS plot
#     p <- ggplotly(ggplot(data=dgntblplt, aes(x=genemid, y=pvalplt, color=as.factor(kngene))) +
#                     geom_point(pch=15) +
#                     geom_segment(aes(x = genestart, y = pvalplt, xend = genestop, yend = pvalplt, color=as.factor(kngene)), size=2) +
#                     geom_text(aes(x=genemid,y=pvalplt,label=genename,color=as.factor(kngene)), nudge_y = nudgeval, size=4) +
#                     geom_hline(aes(yintercept=pthresh), lty=2, color="red") +
#                     xlim(xlow(),xhigh()) + ylim(0,yhigh) +
#                     annotate(geom="text",x=xlow()+270000, y=pthresh-nudgeval, color='red',
#                              label=paste0("TWAS p-value: ",formatC(10^(-pthresh), format = "e", digits = 2)),size=4) +
#                     theme(legend.position = 'top', legend.title = element_blank()) +
#                     theme_bw() +
#                     colScale
#                     )
# 
#     p <- p %>% layout(title = list(text=loctitle(), x=0, xanchor='left', y=.99),
#                       yaxis = list(title = ' TWAS -log10(p)'),
#                       legend = list(orientation='h', x=0, y=1))
# 
#     # GWAS plot
#     m <- ggplotly(ggplot(data=gwaslocnotkn, aes(x=as.numeric(pos),y=poslog10p, color=as.factor(knsnp))) +
#                     geom_point() +
#                     geom_hline(aes(yintercept=(log10(5*10^(-8)))), lty=2, color="red") +
#                     xlim(xlow(),xhigh()) + ylim(ylow,0) +
#                     geom_point(data=knowngwas, aes(x=as.numeric(pos),y=poslog10p, color=as.factor(knsnp))) +
#                     annotate(geom="text",x=xlow()+270000, y=log10(5*10^(-8))+nudgevalg, color='red',
#                              label=paste0("GWAS p-value: ",formatC(5*10^(-8), format = "e", digits = 2)),size=4) +
#                     theme_bw() +
#                     colScaleg
#                   )
# 
#     m <- m %>% layout(xaxis = list(title="position"),
#                       yaxis=list(title="GWAS log10(p)"))
# 
#     # combine plots together
#     fig <- subplot(style(p,showlegend=F), m, shareX = T, nrows = 2, titleX = T,titleY = T,
#                    which_layout = 1, margin=0.005)
# 
#   })
# 
#   ############################################################
# 
# 
#   # TWAS correlation plots
#   output$TWAScorr <- renderPlotly({
# 
#     # select significant and known genes to plot
#     dgntblplt <- dgntbl() %>% filter(SignifGene==1 | kngene=="Reported in GWAS")
#     uniqgenes <- dgntblplt %>% select(gene)
#     indexgene <- dgntblplt$gene[dgntblplt$p == min(dgntblplt$p)]
# 
# 
#     # extract set of genes from correlation matrix
#     Mindex <- data.table::data.table(M[uniqgenes$gene,indexgene])
#     colnames(Mindex) <- c("corr")
#     Mindex$gene <- uniqgenes$gene
# 
# 
#     # categorize the correlation values into bins for plotting
#     Mindex[,corgroup:= cut(Mindex[,abs(corr)],
#                       c(0,.2,.4,.6,.8,1),
#                       labels=c("[0 - 0.2]","(0.2 - 0.4]","(0.4 - 0.6]","(0.6 - 0.8]","(0.8 - 1]"),
#                       right=T,include.lowest=T)]
# 
#     # color for correlation categories
#     Mindex[,corcol:= cut(Mindex[,abs(corr)],
#                          c(0,.2,.4,.6,.8,1),
#                          labels=c("#CCCCCC","#333FFF","#00FF00","#FF9900","#FF0000"),
#                          right=T,include.lowest = T)]
# 
# 
#     # match correlation values and categories back to overall table
#     corplt <- merge(dgntblplt,Mindex,by.x="gene",by.y="gene")
# 
#     # get the DGN weights for the specific genes at locus
#     weight_dsloc <- weight_ds %>% filter(gene %in% corplt$gene)
# 
#     # subset the GWAS variants for specific chr, and phenotype
#     gwasallfinloc <- gwasallfin %>% filter(Locus==locnum()) %>%
#       select(pos,all1,all2,poslog10p)
#     gwasallfinloc$posnum <- as.numeric(gwasallfinloc$pos)
# 
#     # get the matching GWAS variants for the weights
#     weight_dsgwasloc1 <- merge(weight_dsloc,gwasallfinloc,by.x=c('V7','ref_allele','eff_allele'),
#                            by.y=c('posnum','all1','all2'))
#     weight_dsgwasloc2 <- merge(weight_dsloc,gwasallfinloc,by.x=c('V7','ref_allele','eff_allele'),
#                            by.y=c('posnum','all2','all1'))
#     weight_dsgwasloc <- rbind(weight_dsgwasloc1,weight_dsgwasloc2)
# 
#     # merge GWAS results with TWAS info
#     TWASloc <- corplt %>% select(gene,genename,genestartMB,genemid,genestopMB,log10pval,corgroup)
#     TWASloc$genemidMB <- round(TWASloc$genemid/1000000,4)
#     weight_dsgwaslocfin <- merge(weight_dsgwasloc,TWASloc, by.x='gene',by.y='gene',all.x=T)
# 
#     # set y limit
#     yhigh <- max(corplt$log10pval)+0.25*max(corplt$log10pval)
#     ylow <- min(weight_dsgwasloc$poslog10p)+0.15*min(weight_dsgwasloc$poslog10p)
#     nudgeval <- 0.05*max(abs(yhigh))
# 
#     # set color scale values
#     colval <- c("[0 - 0.2]"="#CCCCCC","(0.2 - 0.4]"="#333FFF",
#                  "(0.4 - 0.6]"="#00FF00","(0.6 - 0.8]"="#FF9900","(0.8 - 1]"="#FF0000")
#     breakval <- c("[0 - 0.2]","(0.2 - 0.4]","(0.4 - 0.6]","(0.6 - 0.8]","(0.8 - 1]")
# 
# 
#     # select the LD for the given locus
#     LDdsloc <- as.data.table(filter(LDds,locus2==locnum()))
# 
#     # define correlation and color groups
#     LDdsloc[,ldcol:= cut(LDdsloc[,abs(corrab)],
#                          c(0,.2,.4,.6,.8,1),
#                          labels=c("#CCCCCC","#333FFF","#00FF00","#FF9900","#FF0000"),
#                          right=T,include.lowest = T)]
# 
#     LDdsloc[,ldgroup:= cut(LDdsloc[,abs(corrab)],
#                            c(0,.2,.4,.6,.8,1),
#                            labels=c("[0 - 0.2]","(0.2 - 0.4]","(0.4 - 0.6]","(0.6 - 0.8]","(0.8 - 1]"),
#                            right=T,include.lowest=T)]
# 
# 
# 
#     # match snp LD back with overall file
#     weight_dsgwaslocfin2 <- merge(weight_dsgwaslocfin,LDdsloc, by.x=c('V7'),by.y=c('posb'),all.x=T)
# 
#     weight_dsgwaslocfin2_comp <- weight_dsgwaslocfin2[complete.cases(weight_dsgwaslocfin2), ]
#     weight_dsgwaslocfin2_comp$line <- ifelse(weight_dsgwaslocfin2_comp$weight<0,2,1)
# 
#     xlowMB <- min(corplt$genestartMB,round(weight_dsgwaslocfin2_comp$V7/1000000,4))-.0005
#     xhighMB <- max(corplt$genestopMB,round(weight_dsgwaslocfin2_comp$V7/1000000,4))+.0005
# 
# 
#     # plot correlation categories at the locus, like locus zoom plot
#     loczoom <- ggplotly(ggplot(data=corplt, aes(x=round(genemid/1000000,4), y=log10pval)) +
#       geom_hline(aes(yintercept=pthresh), lty=2, color="red") +
#       geom_hline(aes(yintercept=log10(5*10^(-8))), lty=2, color="red") +
#       geom_hline(aes(yintercept=0),size=1,color='black') +
#       geom_segment(data=weight_dsgwaslocfin2_comp,aes(x=round(V7/1000000,4),y=poslog10p,xend=genemidMB,yend=log10pval
#                                             ,color=ldgroup)) +
#       geom_point(aes(color=corgroup), pch=15) +
#       geom_segment(aes(x=genestartMB,y=log10pval,xend=genestopMB,yend=log10pval, color=corgroup),size=2) +
#       geom_text(aes(label=genename,color=corgroup),nudge_y = nudgeval) +
#       scale_color_manual(breaks = breakval, values = colval) +
#       # scale_color_manual(name = NULL,
#       #    values = colval,
#       #    breaks = breakval,
#       #    guide = guide_legend(override.aes = list(shape = rep(15, 5),
#       #                                             color = c("#CCCCCC","#333FFF",
#       #                                                       "#00FF00","#FF9900","#FF0000")) ) ) +
#       theme_bw() +
#       #xlim(xlowMB(),xhighMB()) +
#       xlim(xlowMB,xhighMB) +
#       ylim(ylow,yhigh) +
#       theme(legend.position = 'top', legend.title = element_blank()) +
#       geom_point(data=weight_dsgwaslocfin2_comp, aes(x=round(V7/1000000,4),y=poslog10p,color=ldgroup)) #+
#     )
# 
#     loczoom <- loczoom %>% layout(title = list(text="(a)", x=0, xanchor='left', y=.99),
#                       yaxis = list(title = '<-- GWAS log10(p) ... 0 ... TWAS -log10(p) -->'),
#                       legend = list(orientation='h', x=0, y=1),
#                       xaxis = list(title="position (in Mb)"))
#   })
# 
# 
#   ############################################################
# 
# 
#   # TWAS network plots
#   output$TWASnetwork <- renderVisNetwork({
# 
#     # select significant and known genes to plot
#     dgntblplt <- dgntbl() %>% filter(SignifGene==1 | kngene=="Reported in GWAS")
#     uniqgenes <- dgntblplt %>% select(gene)
#     indexgene <- dgntblplt$gene[dgntblplt$p == min(dgntblplt$p)]
# 
#     # extract set of genes from correlation matrix
#     Mindex <- data.table::data.table(M[uniqgenes$gene,indexgene])
#     colnames(Mindex) <- c("corr")
#     Mindex$gene <- uniqgenes$gene
# 
#     # categorize the correlation values into bins for plotting
#     Mindex[,corgroup:= cut(Mindex[,abs(corr)],
#                            c(0,.2,.4,.6,.8,1),
#                            labels=c("[0 - 0.2]","(0.2 - 0.4]","(0.4 - 0.6]","(0.6 - 0.8]","(0.8 - 1]"),
#                            right=T,include.lowest=T)]
# 
#     # color for correlation categories
#     Mindex[,corcol:= cut(Mindex[,abs(corr)],
#                          c(0,.2,.4,.6,.8,1),
#                          labels=c("#CCCCCC","#333FFF","#00FF00","#FF9900","#FF0000"),
#                          right=T,include.lowest = T)]
# 
#     # match correlation values and categories back to overall table
#     corplt <- merge(dgntblplt,Mindex,by.x="gene",by.y="gene")
# 
#     # get the DGN weights for the specific genes at locus
#     weight_dsloc <- weight_ds %>% filter(gene %in% corplt$gene)
# 
#     # subset the GWAS variants for specific chr, and phenotype
#     gwasallfinloc <- gwasallfin %>% filter(Locus==locnum()) %>%
#       select(pos,all1,all2,poslog10p)
#     gwasallfinloc$posnum <- as.numeric(gwasallfinloc$pos)
# 
#     # get the matching GWAS variants for the weights
#     weight_dsgwasloc1 <- merge(weight_dsloc,gwasallfinloc,by.x=c('V7','ref_allele','eff_allele'),
#                            by.y=c('posnum','all1','all2'))
#     weight_dsgwasloc2 <- merge(weight_dsloc,gwasallfinloc,by.x=c('V7','ref_allele','eff_allele'),
#                            by.y=c('posnum','all2','all1'))
#     weight_dsgwasloc <- rbind(weight_dsgwasloc1,weight_dsgwasloc2)
# 
#     # merge GWAS results with TWAS info
#     TWASloc <- corplt %>% select(gene,genename,genestartMB,genemid,genestopMB,log10pval,corgroup)
#     TWASloc$genemidMB <- round(TWASloc$genemid/1000000,4)
#     weight_dsgwaslocfin <- merge(weight_dsgwasloc,TWASloc, by.x='gene',by.y='gene',all.x=T)
# 
#     # set y limit
#     yhigh <- max(corplt$log10pval)+0.15*max(corplt$log10pval)
#     ylow <- min(weight_dsgwasloc$poslog10p)+0.15*min(weight_dsgwasloc$poslog10p)
#     nudgeval <- 0.05*max(abs(yhigh))
# 
#     # set color scale values
#     colval <- c("[0 - 0.2]"="#CCCCCC","(0.2 - 0.4]"="#333FFF",
#                 "(0.4 - 0.6]"="#00FF00","(0.6 - 0.8]"="#FF9900","(0.8 - 1]"="#FF0000")
# 
#     # select the LD for the given locus
#     LDdsloc <- as.data.table(filter(LDds,locus2==locnum()))
# 
#     # define correlation and color groups
#     LDdsloc[,ldcol:= cut(LDdsloc[,abs(corrab)],
#                          c(0,.2,.4,.6,.8,1),
#                          labels=c("#CCCCCC","#333FFF","#00FF00","#FF9900","#FF0000"),
#                          right=T,include.lowest = T)]
# 
#     LDdsloc[,ldgroup:= cut(LDdsloc[,abs(corrab)],
#                            c(0,.2,.4,.6,.8,1),
#                            labels=c("[0 - 0.2]","(0.2 - 0.4]","(0.4 - 0.6]","(0.6 - 0.8]","(0.8 - 1]"),
#                            right=T,include.lowest=T)]
# 
#     # match snp LD back with overall file
#     weight_dsgwaslocfin2 <- merge(weight_dsgwaslocfin,LDdsloc, by.x=c('V7'),by.y=c('posb'),all.x=T)
# 
#     weight_dsgwaslocfin2_comp <- weight_dsgwaslocfin2[complete.cases(weight_dsgwaslocfin2), ]
#     weight_dsgwaslocfin2_comp$line <- ifelse(weight_dsgwaslocfin2_comp$weight<0,2,1)
# 
# 
#     # data for links
#     from <- weight_dsgwaslocfin2$rsid
#     to <- weight_dsgwaslocfin2$genename
#     weight <- weight_dsgwaslocfin2$weight
#     color <- weight_dsgwaslocfin2$ldcol
#     length <- 1/(weight_dsgwaslocfin2$log10pval-weight_dsgwaslocfin2$poslog10p)*1000
#     links <- data.frame(from,to,weight,color,length)
#     links$dashes <- ifelse(links$weight<0,TRUE,FALSE)
#     links$width <- abs(links$weight)*100
# 
#     links_complete <- links[complete.cases(links), ]
# 
# 
#     # data for gene nodes
#     nodesgene <- corplt %>% select(genename,log10pval,corcol)
#     colnames(nodesgene) <- c("id","size","color.background")
#     nodesgene$title=nodesgene$id
#     nodesgene$type="gene"
#     nodesgene$shape="square"
#     nodesgene$shadow=TRUE
#     nodesgene$color.border="black"
#     nodesgene$label=nodesgene$id
#     maxp <- max(nodesgene$size)
#     nodesgene$shape[nodesgene$size==maxp]="star"
#     # nodesgene$Gene = nodesgene$label
#     # topgene <- nodesgene$Gene[nodesgene$size==maxp]
#     # nodesgene$Gene[nodesgene$size==maxp]=paste0(topgene,"*")
# 
# 
#     nodessnp <- weight_dsgwaslocfin2 %>% select(rsid,poslog10p,ldcol)
#     colnames(nodessnp) <- c("id","poslog10p","color.background")
#     nodessnp$size <- abs(nodessnp$poslog10p)
#     nodessnp$title=nodessnp$id
#     nodessnp$type="snp"
#     nodessnp$shape="dot"
#     nodessnp$shadow=FALSE
#     nodessnp$color.border="black"
#     nodessnp$label=""
#     nodessnp$shape[nodessnp$size==max(nodessnp$size)]="triangle"
#     #nodessnp$Gene[nodessnp$Gene==topgene]=paste0(topgene,"*")
# 
# 
#     nodessnp <- select(nodessnp,-poslog10p)
#     nodessnpfin <- distinct(nodessnp)
# 
#     nodes <- rbind(nodesgene,nodessnpfin)
#     nodes_complete <- nodes[complete.cases(nodes), ]
# 
#     nodes_complete$color.highlight.background <- "purple"
#     nodes_complete$color.highlight.border <- "purple"
#     #nodes_complete$font.size <- max(nodessnp$size, nodesgene$size)*2
#     nodes_complete$font.size <- 20
# 
#     visnet <- visNetwork(nodes_complete,links_complete, main="(b)")
#     visOptions(visnet, highlightNearest = TRUE)
# 
# 
#   })
# 
# 
#   ############################################################
# 
#   # TWAS/GWAS mirror plot of meta-analysis
#   output$Tmetamirror <- renderPlotly({
# 
#     # select significant and known genes to plot
#     dgntblplt <- dgntbl() %>% filter(SignifGene==1 | kngene=="Reported in GWAS")
# 
# 
#     # define y limits and nudges for plotting
#     yhigh <- max(locds()$log10pval)+0.25*max(locds()$log10pval)
#     ylow <- min(-dgntblplt$log10pvalmeta,metathresh)+0.15*min(-dgntblplt$log10pvalmeta,metathresh)
#     nudgeval <- 0.05*max(abs(yhigh))
#     nudgevalg <- 0.05*ylow
# 
# 
#     # define specific colors for known/not known genes
#     myColors <- setNames( c('#56B4E9', '#000000', '#E69F00'),
#                           c("Not reported in GWAS","Reported in GWAS",NA))
# 
#     colScale <- scale_colour_manual(values = myColors)
# 
# 
#     # TWAS plot
#     p <- ggplotly(ggplot(data=dgntblplt, aes(x=genemid, y=log10pval, color=as.factor(kngene))) +
#                     geom_point(pch=15) +
#                     geom_segment(aes(x = genestart, y = log10pval, xend = genestop, yend = log10pval,
#                                       color=as.factor(kngene)), size=2) +
#                     geom_text(aes(x=genemid,y=log10pval,label=genename,color=as.factor(kngene)),
#                                nudge_y = nudgeval, size=4) +
#                     geom_hline(aes(yintercept=pthresh), lty=2, color="red") +
#                     xlim(xlow(),xhigh()) + ylim(0,yhigh) +
#                     #xlim(xlow,xhigh) + ylim(0,yhigh) +
#                     annotate(geom="text",x=xlow()+270000, y=pthresh-nudgeval, color='red',
#                             label=paste0("TWAS p-value: ",formatC(10^(-pthresh), format = "e", digits = 2)),size=4) +
#                     theme(legend.position = 'top', legend.title = element_blank()) +
#                     theme_bw() +
#                     colScale
#     )
# 
#     p <- p %>% layout(title = list(text=loctitle(), x=0, xanchor='left', y=.99),
#                       yaxis = list(title = ' TWAS -log10(p)'),
#                       legend = list(orientation='h', x=0, y=1))
# 
#     # Meta-analysis plot
#     m <- ggplotly(ggplot(data=dgntblplt, aes(x=genemid, y=-log10pvalmeta, color=as.factor(kngene))) +
#                      geom_point(pch=15) +
#                      geom_segment(aes(x = genestart, y = -log10pvalmeta, xend = genestop, yend = -log10pvalmeta,
#                                       color=as.factor(kngene)), size=2) +
#                      geom_text(aes(x=genemid,y=(-log10pvalmeta),label=genename,color=as.factor(kngene)),
#                                nudge_y = nudgevalg, size=4) +
#                      geom_hline(aes(yintercept=log10(0.05)), lty=2, color="orange") +
#                      geom_hline(aes(yintercept=metathresh), lty=2, color="red") +
#                      xlim(xlow(),xhigh()) + ylim(ylow,0) +
#                     # #xlim(xlow,xhigh) + ylim(ylow,0) +
#                     annotate(geom="text",x=xlow()+270000, y=log10(0.05)+nudgevalg, color='orange',
#                             label=paste0("Replication p-value: 0.05"),size=4) +
#                     annotate(geom="text",x=xlow()+370000, y=metathresh+nudgevalg, color='red',
#                             label=paste0("Strict replication p-value: ", formatC(10^(metathresh), format = "e", digits = 2)),size=4) +
#                     theme(legend.position = 'top', legend.title = element_blank()) +
#                     theme_bw() +
#                     colScale
#     )
# 
#     m <- m %>% layout(xaxis = list(title="position"),
#                       yaxis=list(title="Meta-analysis log10(p)"))
# 
#     # combine plots together
#     fig <- subplot(style(p,showlegend=F), m, shareX = T, nrows = 2, titleX = T,titleY = T,
#                    which_layout = 1, margin=0.005)
# 
#   })
# 
#   ############################################################
# 
# 
#   output$CompRefGWB <- renderPlotly({
# 
#     # select DGN results
#     locds <- dgnds %>% filter(locvar==input$locuslst)
#     xhigh <- max(locds$genestop)+1000000
#     xlow <- max(0,min(locds$genestart)-1000000)
#     locchr <- unique(locds$chr)
#     locpheno <- unique(locds$pheno)
# 
# 
#     #####################
# 
# 
#     # select reference panel specific results
#     dgntbl <- dgnds %>% filter(genestart>=xlow & genestop<=xhigh & pheno==locpheno & chr==locchr &
#                                  is.na(HLARegion) & is.na(MHCRegion)) %>%
#       select(genename,chr, pheno,genestart,genestop,genemid,log10pval)
# 
#     locgwb <- twas_ds %>% filter(tissue=="GWB" & genestart>=xlow & genestop<=xhigh & chr==locchr &
#                                     pheno==locpheno & is.na(HLARegion) & is.na(MHCRegion)) %>%
#       select(genename,chr, pheno,log10pval,genestart,genestop,genemid, SignifGene)
#     pthreshgwb <- -log10(0.05/(nrow(twas_ds[twas_ds$tissue=="GWB",])))
# 
#     locgtl <- twas_ds %>% filter(tissue=="GTL" & genestart>=xlow & genestop<=xhigh & chr==locchr &
#                                     pheno==locpheno & is.na(HLARegion) & is.na(MHCRegion)) %>%
#       select(genename,chr, pheno,log10pval,genestart,genestop,genemid, SignifGene)
#     pthreshgtl <- -log10(0.05/(nrow(twas_ds[twas_ds$tissue=="GTL",])))
# 
#     locmsa <- twas_ds %>% filter(tissue=="MSA" & genestart>=xlow & genestop<=xhigh & chr==locchr &
#                                     pheno==locpheno & is.na(HLARegion) & is.na(MHCRegion)) %>%
#       select(genename,chr, pheno,log10pval,genestart,genestop,genemid, SignifGene)
#     pthreshmsa <- -log10(0.05/(nrow(twas_ds[twas_ds$tissue=="MSA",])))
# 
# 
#     #####################
# 
# 
#     # merge secondary reference panel data sets with DGN
#     dgngwb <- merge(dgntbl,locgwb, by=c("genename","chr","pheno"), all = T)
#     dgngwb$inboth <- ifelse(complete.cases(dgngwb$log10pval.x,dgngwb$log10pval.y),"In Both","Not in Both")
# 
#     dgngtl <- merge(dgntbl,locgtl, by=c("genename","chr","pheno"), all = T)
#     dgngtl$inboth <- ifelse(complete.cases(dgngtl$log10pval.x,dgngtl$log10pval.y),"In Both","Not in Both")
# 
#     dgnmsa <- merge(dgntbl,locmsa, by=c("genename","chr","pheno"), all = T)
#     dgnmsa$inboth <- ifelse(complete.cases(dgnmsa$log10pval.x,dgnmsa$log10pval.y),"In Both","Not in Both")
# 
# 
#     #####################
# 
# 
#     # set y limit values for each plot
#     yhighdgn <- max(dgngwb$log10pval.x, na.rm = T)+0.15*max(dgngwb$log10pval.x, na.rm = T)
#     nudgevaldgn <- 0.05*yhighdgn
# 
#     yhighgtl <- max(dgngtl$log10pval.y,pthreshgtl, na.rm = T)+0.15*max(dgngtl$log10pval.y,pthreshgtl, na.rm = T)
#     nudgevalgtl <- 0.05*yhighgtl
# 
#     yhighgwb <- max(dgngwb$log10pval.y,pthreshgwb, na.rm = T)+0.15*max(dgngwb$log10pval.y,pthreshgwb, na.rm = T)
#     nudgevalgwb <- 0.05*yhighgwb
# 
#     yhighmsa <- max(dgnmsa$log10pval.y,pthreshmsa, na.rm = T)+0.15*max(dgnmsa$log10pval.y,pthreshmsa, na.rm = T)
#     nudgevalmsa <- 0.05*yhighmsa
# 
# 
#     myColors <- setNames( c('#000000', '#56B4E9', '#E69F00'),
#                           c("In Both","Not in Both",NA))
# 
#     colScale <- scale_colour_manual(values = myColors)
# 
# 
#     ####################
# 
#     # plot info for GWB
#     atop <- ggplotly(ggplot(data=dgngwb, aes(x=genemid.x,y=log10pval.x, color=inboth)) +
#                        geom_point(pch=15) +
#                        geom_segment(aes(x=genestart.x,y=log10pval.x, xend=genestop.x, yend=log10pval.x), size=2) +
#                        geom_text(aes(label=genename), nudge_y = nudgevaldgn) +
#                        geom_hline(aes(yintercept=pthresh), lty=2, color='red') +
#                        xlim(xlow,xhigh) +
#                        ylim(0,yhighdgn) +
#                        theme_bw() +
#                        theme(legend.position = 'top', legend.title = element_blank()) +
#                        colScale
#     )
# 
#     atop <- atop %>% layout(yaxis = list(title = ' DGN TWAS -log10(p)'),
#                             xaxis = list(range=c(xlow,xhigh)),
#                             annotations=list(x = 0.5 , y = 1.1, text = "(a) DGN vs. GWB", showarrow = F,
#                                              xref='paper', yref='paper',xanchor='center'),
#                             legend = list(orientation='h', x=0, y=1),
#                             title=list(text=paste0(loctitle(),"\n"),x=0,xanchor='left'), margin=list(t=60)
#     )
# 
#     #####
# 
#     abottom <- ggplotly(ggplot(data=dgngwb, aes(x=genemid.y,y=-log10pval.y, color=inboth)) +
#                           geom_hline(aes(yintercept=-pthreshgwb), lty=2, color='red') +
#                           geom_point(pch=15) +
#                           geom_segment(aes(x=genestart.y,y=-log10pval.y, xend=genestop.y, yend=-log10pval.y), size=2) +
#                           geom_text(aes(label=genename), nudge_y = -1.5*nudgevalgwb) +
#                           xlim(xlow,xhigh) +
#                           ylim(-yhighgwb,0) +
#                           theme_bw() +
#                           theme(legend.position = 'top', legend.title = element_blank()) +
#                           colScale
#     )
# 
#     abottom <- abottom %>% layout(yaxis = list(title = 'GWB TWAS log10(p)',
#                                                range=c(-yhighgwb,0.1)),
#                                   xaxis = list(range=c(xlow,xhigh), title="position"))
# 
#     #####
# 
#     afig <- subplot(atop,style(abottom,showlegend=F),nrows=2,shareX = T, titleX = T, titleY = T, which_layout = 1, margin=.01)
# 
#   })
# 
# 
#   output$CompRefGTL <- renderPlotly({
# 
#     # select DGN results
#     locds <- dgnds %>% filter(locvar==input$locuslst)
#     xhigh <- max(locds$genestop)+1000000
#     xlow <- max(0,min(locds$genestart)-1000000)
#     locchr <- unique(locds$chr)
#     locpheno <- unique(locds$pheno)
# 
# 
#     #####################
# 
# 
#     # select reference panel specific results
#     dgntbl <- dgnds %>% filter(genestart>=xlow & genestop<=xhigh & pheno==locpheno & chr==locchr &
#                                  is.na(HLARegion) & is.na(MHCRegion)) %>%
#       select(genename,chr, pheno,genestart,genestop,genemid,log10pval)
# 
#     locgwb <- twas_ds %>% filter(tissue=="GWB" & genestart>=xlow & genestop<=xhigh & chr==locchr &
#                                     pheno==locpheno & is.na(HLARegion) & is.na(MHCRegion)) %>%
#       select(genename,chr, pheno,log10pval,genestart,genestop,genemid, SignifGene)
#     pthreshgwb <- -log10(0.05/(nrow(twas_ds[twas_ds$tissue=="GWB",])))
# 
#     locgtl <- twas_ds %>% filter(tissue=="GTL" & genestart>=xlow & genestop<=xhigh & chr==locchr &
#                                     pheno==locpheno & is.na(HLARegion) & is.na(MHCRegion)) %>%
#       select(genename,chr, pheno,log10pval,genestart,genestop,genemid, SignifGene)
#     pthreshgtl <- -log10(0.05/(nrow(twas_ds[twas_ds$tissue=="GTL",])))
# 
#     locmsa <- twas_ds %>% filter(tissue=="MSA" & genestart>=xlow & genestop<=xhigh & chr==locchr &
#                                     pheno==locpheno & is.na(HLARegion) & is.na(MHCRegion)) %>%
#       select(genename,chr, pheno,log10pval,genestart,genestop,genemid, SignifGene)
#     pthreshmsa <- -log10(0.05/(nrow(twas_ds[twas_ds$tissue=="MSA",])))
# 
# 
#     #####################
# 
# 
#     # merge secondary reference panel data sets with DGN
#     dgngwb <- merge(dgntbl,locgwb, by=c("genename","chr","pheno"), all = T)
#     dgngwb$inboth <- ifelse(complete.cases(dgngwb$log10pval.x,dgngwb$log10pval.y),"In Both","Not in Both")
# 
#     dgngtl <- merge(dgntbl,locgtl, by=c("genename","chr","pheno"), all = T)
#     dgngtl$inboth <- ifelse(complete.cases(dgngtl$log10pval.x,dgngtl$log10pval.y),"In Both","Not in Both")
# 
#     dgnmsa <- merge(dgntbl,locmsa, by=c("genename","chr","pheno"), all = T)
#     dgnmsa$inboth <- ifelse(complete.cases(dgnmsa$log10pval.x,dgnmsa$log10pval.y),"In Both","Not in Both")
# 
# 
#     #####################
# 
# 
#     # set y limit values for each plot
#     yhighdgn <- max(dgngwb$log10pval.x, na.rm = T)+0.15*max(dgngwb$log10pval.x, na.rm = T)
#     nudgevaldgn <- 0.05*yhighdgn
# 
#     yhighgtl <- max(dgngtl$log10pval.y,pthreshgtl, na.rm = T)+0.15*max(dgngtl$log10pval.y,pthreshgtl, na.rm = T)
#     nudgevalgtl <- 0.05*yhighgtl
# 
#     yhighgwb <- max(dgngwb$log10pval.y,pthreshgwb, na.rm = T)+0.15*max(dgngwb$log10pval.y,pthreshgwb, na.rm = T)
#     nudgevalgwb <- 0.05*yhighgwb
# 
#     yhighmsa <- max(dgnmsa$log10pval.y,pthreshmsa, na.rm = T)+0.15*max(dgnmsa$log10pval.y,pthreshmsa, na.rm = T)
#     nudgevalmsa <- 0.05*yhighmsa
# 
# 
#     myColors <- setNames( c('#000000', '#56B4E9', '#E69F00'),
#                           c("In Both","Not in Both",NA))
# 
#     colScale <- scale_colour_manual(values = myColors)
# 
# 
#     ####################
# 
#     # plot info for GTL
#     btop <- ggplotly(ggplot(data=dgngtl, aes(x=genemid.x,y=log10pval.x, color=inboth)) +
#                        geom_point(pch=15) +
#                        geom_segment(aes(x=genestart.x,y=log10pval.x, xend=genestop.x, yend=log10pval.x), size=2) +
#                        geom_text(aes(label=genename), nudge_y = nudgevaldgn) +
#                        geom_hline(aes(yintercept=pthresh), lty=2, color='red') +
#                        xlim(xlow,xhigh) +
#                        ylim(0,yhighdgn) +
#                        theme_bw() +
#                        theme(legend.position = 'top', legend.title = element_blank()) +
#                        colScale
#     )
# 
#     btop <- btop %>% layout(yaxis = list(title = ' DGN TWAS -log10(p)'),
#                             xaxis = list(range=c(xlow,xhigh)),
#                             annotations=list(x = 0.5 , y = 1.1, text = "(b) DGN vs. GTL", showarrow = F,
#                                              xref='paper', yref='paper',xanchor='center'),
#                             legend = list(orientation='h', x=0, y=1),
#                             title=list(text=paste0(loctitle(),"\n"),x=0,xanchor='left'), margin=list(t=60))
# 
#     #####
# 
#     bbottom <- ggplotly(ggplot(data=dgngtl, aes(x=genemid.y,y=-log10pval.y, color=inboth)) +
#                           geom_hline(aes(yintercept=-pthreshgtl), lty=2, color='red') +
#                           geom_point(pch=15) +
#                           geom_segment(aes(x=genestart.y,y=-log10pval.y, xend=genestop.y, yend=-log10pval.y), size=2) +
#                           geom_text(aes(label=genename), nudge_y = -1.5*nudgevalgtl) +
#                           xlim(xlow,xhigh) +
#                           ylim(-yhighgtl,0) +
#                           theme_bw() +
#                           theme(legend.position = 'top', legend.title = element_blank()) +
#                           colScale
#     )
# 
#     bbottom <- bbottom %>% layout(yaxis = list(title = 'GTL TWAS log10(p)',
#                                                range=c(-yhighgtl,0.1)),
#                                   xaxis = list(range=c(xlow,xhigh), title="position"))
# 
#     #####
# 
#     bfig <- subplot(btop,style(bbottom,showlegend=F),nrows=2,shareX = T, titleX = T, titleY = T, which_layout = 1, margin=0.01)
# 
#   })
# 
# 
# 
#   output$CompRefMSA <- renderPlotly({
# 
#     # select DGN results
#     locds <- dgnds %>% filter(locvar==input$locuslst)
#     xhigh <- max(locds$genestop)+1000000
#     xlow <- max(0,min(locds$genestart)-1000000)
#     locchr <- unique(locds$chr)
#     locpheno <- unique(locds$pheno)
# 
# 
#     #####################
# 
# 
#     # select reference panel specific results
#     dgntbl <- dgnds %>% filter(genestart>=xlow & genestop<=xhigh & pheno==locpheno & chr==locchr &
#                                  is.na(HLARegion) & is.na(MHCRegion)) %>%
#       select(genename,chr, pheno,genestart,genestop,genemid,log10pval)
# 
#     locgwb <- twas_ds %>% filter(tissue=="GWB" & genestart>=xlow & genestop<=xhigh & chr==locchr &
#                                     pheno==locpheno & is.na(HLARegion) & is.na(MHCRegion)) %>%
#       select(genename,chr, pheno,log10pval,genestart,genestop,genemid, SignifGene)
#     pthreshgwb <- -log10(0.05/(nrow(twas_ds[twas_ds$tissue=="GWB",])))
# 
#     locgtl <- twas_ds %>% filter(tissue=="GTL" & genestart>=xlow & genestop<=xhigh & chr==locchr &
#                                     pheno==locpheno & is.na(HLARegion) & is.na(MHCRegion)) %>%
#       select(genename,chr, pheno,log10pval,genestart,genestop,genemid, SignifGene)
#     pthreshgtl <- -log10(0.05/(nrow(twas_ds[twas_ds$tissue=="GTL",])))
# 
#     locmsa <- twas_ds %>% filter(tissue=="MSA" & genestart>=xlow & genestop<=xhigh & chr==locchr &
#                                     pheno==locpheno & is.na(HLARegion) & is.na(MHCRegion)) %>%
#       select(genename,chr, pheno,log10pval,genestart,genestop,genemid, SignifGene)
#     pthreshmsa <- -log10(0.05/(nrow(twas_ds[twas_ds$tissue=="MSA",])))
# 
# 
#     #####################
# 
# 
#     # merge secondary reference panel data sets with DGN
#     dgngwb <- merge(dgntbl,locgwb, by=c("genename","chr","pheno"), all = T)
#     dgngwb$inboth <- ifelse(complete.cases(dgngwb$log10pval.x,dgngwb$log10pval.y),"In Both","Not in Both")
# 
#     dgngtl <- merge(dgntbl,locgtl, by=c("genename","chr","pheno"), all = T)
#     dgngtl$inboth <- ifelse(complete.cases(dgngtl$log10pval.x,dgngtl$log10pval.y),"In Both","Not in Both")
# 
#     dgnmsa <- merge(dgntbl,locmsa, by=c("genename","chr","pheno"), all = T)
#     dgnmsa$inboth <- ifelse(complete.cases(dgnmsa$log10pval.x,dgnmsa$log10pval.y),"In Both","Not in Both")
# 
# 
#     #####################
# 
# 
#     # set y limit values for each plot
#     yhighdgn <- max(dgngwb$log10pval.x, na.rm = T)+0.15*max(dgngwb$log10pval.x, na.rm = T)
#     nudgevaldgn <- 0.05*yhighdgn
# 
#     yhighgtl <- max(dgngtl$log10pval.y,pthreshgtl, na.rm = T)+0.15*max(dgngtl$log10pval.y,pthreshgtl, na.rm = T)
#     nudgevalgtl <- 0.05*yhighgtl
# 
#     yhighgwb <- max(dgngwb$log10pval.y,pthreshgwb, na.rm = T)+0.15*max(dgngwb$log10pval.y,pthreshgwb, na.rm = T)
#     nudgevalgwb <- 0.05*yhighgwb
# 
#     yhighmsa <- max(dgnmsa$log10pval.y,pthreshmsa, na.rm = T)+0.15*max(dgnmsa$log10pval.y,pthreshmsa, na.rm = T)
#     nudgevalmsa <- 0.05*yhighmsa
# 
# 
#     myColors <- setNames( c('#000000', '#56B4E9', '#E69F00'),
#                           c("In Both","Not in Both",NA))
# 
#     colScale <- scale_colour_manual(values = myColors)
# 
# 
#     ####################
# 
#     # plot info for MSA
#     ctop <- ggplotly(ggplot(data=dgnmsa, aes(x=genemid.x,y=log10pval.x, color=inboth)) +
#                        geom_point(pch=15) +
#                        geom_segment(aes(x=genestart.x,y=log10pval.x, xend=genestop.x, yend=log10pval.x), size=2) +
#                        geom_text(aes(label=genename), nudge_y = nudgevaldgn) +
#                        geom_hline(aes(yintercept=pthresh), lty=2, color='red') +
#                        xlim(xlow,xhigh) +
#                        ylim(0,yhighdgn) +
#                        theme_bw() +
#                        theme(legend.position = 'top', legend.title = element_blank()) +
#                        colScale
#     )
# 
#     ctop <- ctop %>% layout(yaxis = list(title = 'DGN TWAS -log10(p)'),
#                             xaxis = list(range=c(xlow,xhigh)),
#                             annotations=list(x = 0.5 , y = 1.1, text = "(c) DGN vs. MSA", showarrow = F,
#                                              xref='paper', yref='paper',xanchor='center'),
#                             legend = list(orientation='h', x=0, y=1),
#                             title=list(text=paste0(loctitle(),"\n"),x=0,xanchor='left'),margin=list(t=60))
# 
#     #####
# 
#     cbottom <- ggplotly(ggplot(data=dgnmsa, aes(x=genemid.y,y=-log10pval.y, color=inboth)) +
#                           geom_hline(aes(yintercept=-pthreshmsa), lty=2, color='red') +
#                           geom_point(pch=15) +
#                           geom_segment(aes(x=genestart.y,y=-log10pval.y, xend=genestop.y, yend=-log10pval.y), size=2) +
#                           geom_text(aes(label=genename), nudge_y = -1.5*nudgevalmsa) +
#                           xlim(xlow,xhigh) +
#                           theme_bw() +
#                           ylim(-yhighmsa-2,0) +
#                           colScale
#     )
# 
#     cbottom <- cbottom %>% layout(yaxis = list(title = 'MSA TWAS log10(p)',
#                                                range=c(-yhighmsa,0.1)),
#                                   xaxis = list(range=c(xlow,xhigh), title="position"))
# 
#     #####
# 
#     cfig <- subplot(ctop,style(cbottom,showlegend=F),nrows=2,shareX = T, titleX = T, titleY = T, which_layout = 1,
#                     margin = 0.01)
# 
# 
#     #####################
# 
#   })
# 
# 
# 
#   ############################################################
# 
# 
#   # TWAS results table with separate tabs for each reference panel
#   output$DGNtbl <- DT::renderDataTable({
# 
#     locds <- dgnds %>% filter(locvar==input$locuslst)
#     xhigh <- max(locds$genestop)+1000000
#     xlow <- max(0,min(locds$genestart)-1000000)
#     locpheno <- unique(locds$pheno) #locus phenotype
#     locchr<- unique(locds$chr) #locus chromosome
# 
#     # select all genes at locus
#     dgntbl <- dgnds %>% filter(genestart>=xlow & genestop <=xhigh & pheno==locpheno & chr==locchr) %>%
#       select(genename,chr, pheno,phenocat,se_beta_, p, SignifGene,HLARegion,MHCRegion)
#     colnames(dgntbl) <- c("gene","chr","trait","trait category","beta","p-value","TWAS significant",
#                           "HLA gene","MHC region")
# 
#     # print data table
#     DT::datatable(dgntbl, rownames=F, options = list(
#       order = list(list(5, 'asc')))) %>% DT::formatStyle(
#         "TWAS significant", target = "row",
#         backgroundColor = styleEqual(c(1), c('lightgreen'))
#       ) %>% DT::formatStyle(
#         "HLA gene", target = "row",
#         backgroundColor = styleEqual(c(1), c('red'))
#       ) %>% DT::formatStyle(
#         "MHC region", target = "row",
#         backgroundColor = styleEqual(c(1), c('red'))
#       )
#   })
# 
# 
#   #####################
# 
# 
#   output$GTLtbl <- DT::renderDataTable({
# 
#     locds <- dgnds %>% filter(locvar==input$locuslst)
#     xhigh <- max(locds$genestop)+1000000
#     xlow <- max(0,min(locds$genestart)-1000000)
#     locpheno <- unique(locds$pheno) #locus phenotype
#     locchr<- unique(locds$chr) #locus chromosome
# 
#     # select all genes at locus
#     gtltbl <- twas_ds %>% filter(tissue=='GTL',genestart>=xlow & genestop <=xhigh & pheno==locpheno & chr==locchr) %>%
#       select(genename,chr, pheno,phenocat,se_beta_, p, SignifGene,HLARegion,MHCRegion)
#     colnames(gtltbl) <- c("gene","chr","trait","trait category","beta","p-value","TWAS significant",
#                           "HLA gene","MHC region")
# 
#     # print data table
#     DT::datatable(gtltbl, rownames=F, options = list(
#       order = list(list(5, 'asc')))) %>% formatStyle(
#         "TWAS significant", target = "row",
#         backgroundColor = styleEqual(c(1), c('lightgreen'))
#       ) %>% formatStyle(
#         "HLA gene", target = "row",
#         backgroundColor = styleEqual(c(1), c('red'))
#       ) %>% formatStyle(
#         "MHC region", target = "row",
#         backgroundColor = styleEqual(c(1), c('red'))
#       )
#   })
# 
# 
#   #####################
# 
# 
#   output$GWBtbl <- DT::renderDataTable({
# 
#     locds <- dgnds %>% filter(locvar==input$locuslst)
#     xhigh <- max(locds$genestop)+1000000
#     xlow <- max(0,min(locds$genestart)-1000000)
#     locpheno <- unique(locds$pheno) #locus phenotype
#     locchr<- unique(locds$chr) #locus chromosome
# 
#     # select all genes at locus
#     gwbtbl <- twas_ds %>% filter(tissue=='GWB',genestart>=xlow & genestop <=xhigh & pheno==locpheno & chr==locchr) %>%
#       select(genename,chr, pheno,phenocat,se_beta_, p, SignifGene,HLARegion,MHCRegion)
#     colnames(gwbtbl) <- c("gene","chr","trait","trait category","beta","p-value","TWAS significant",
#                           "HLA gene","MHC region")
# 
#     # print data table
#     DT::datatable(gwbtbl, rownames=F, options = list(
#       order = list(list(5, 'asc')))) %>% DT::formatStyle(
#         "TWAS significant", target = "row",
#         backgroundColor = styleEqual(c(1), c('lightgreen'))
#       ) %>% DT::formatStyle(
#         "HLA gene", target = "row",
#         backgroundColor = styleEqual(c(1), c('red'))
#       ) %>% DT::formatStyle(
#         "MHC region", target = "row",
#         backgroundColor = styleEqual(c(1), c('red'))
#       )
#   })
# 
# 
#   #####################
# 
# 
#   output$MSAtbl <- DT::renderDataTable({
# 
#     locds <- dgnds %>% filter(locvar==input$locuslst)
#     xhigh <- max(locds$genestop)+1000000
#     xlow <- max(0,min(locds$genestart)-1000000)
#     locpheno <- unique(locds$pheno) #locus phenotype
#     locchr<- unique(locds$chr) #locus chromosome
# 
#     # select all genes at locus
#     msatbl <- twas_ds %>% filter(tissue=='MSA',genestart>=xlow & genestop <=xhigh & pheno==locpheno & chr==locchr) %>%
#       select(genename,chr, pheno,phenocat,se_beta_, p, SignifGene,HLARegion,MHCRegion)
#     colnames(msatbl) <- c("gene","chr","trait","trait category","beta","p-value","TWAS significant",
#                           "HLA gene","MHC region")
# 
#     # print data table
#     DT::datatable(msatbl, rownames=F, options = list(
#       order = list(list(5, 'asc')))) %>% DT::formatStyle(
#         "TWAS significant", target = "row",
#         backgroundColor = styleEqual(c(1), c('lightgreen'))
#       ) %>% DT::formatStyle(
#         "HLA gene", target = "row",
#         backgroundColor = styleEqual(c(1), c('red'))
#       ) %>% DT::formatStyle(
#         "MHC region", target = "row",
#         backgroundColor = styleEqual(c(1), c('red'))
#       )
#   })
# 
# 
#   ############################################################
# 
# 
#   # TWAS ohter traits results table with separate tabs for each reference panel
#   output$DGNtbltrt <- DT::renderDataTable({
# 
#     # select dgn all genes from other trait categories at locus
#     dgntbl <- filter(dgnds, genestart>=xlow() & genestop<=xhigh() & pheno!=locpheno()
#                      & chr==locchr() & phenocat==locphcat()) %>%
#         select(genename,chr, pheno,phenocat,se_beta_, p, SignifGene,HLARegion,MHCRegion)
#     colnames(dgntbl) <- c("Gene","Chr","Trait","Trait category","Beta","P-value","TWAS significant",
#                           "HLA gene","MHC region")
# 
#     # print data table
#     DT::datatable(dgntbl, rownames=F, options = list(
#       order = list(list(5, 'asc')))) %>% DT::formatStyle(
#         "TWAS significant", target = "row",
#         backgroundColor = styleEqual(c(1), c('lightgreen'))
#       ) %>% DT::formatStyle(
#         "HLA gene", target = "row",
#         backgroundColor = styleEqual(c(1), c('red'))
#       ) %>% DT::formatStyle(
#         "MHC region", target = "row",
#         backgroundColor = styleEqual(c(1), c('red'))
#       )
#   })
# 
# 
#   #####################
# 
# 
#   output$GTLtbltrt <- DT::renderDataTable({
# 
#     # select all genes at locus
#     gtltbl <- twas_ds %>% filter(tissue=='GTL',genestart>=xlow() & genestop<=xhigh()
#                                   & pheno!=locpheno() & chr==locchr() & phenocat==locphcat()) %>%
#       select(genename,chr, pheno,phenocat,se_beta_, p, SignifGene,HLARegion,MHCRegion)
#     colnames(gtltbl) <- c("Gene","Chr","Trait","Trait category","Beta","P-value","TWAS significant",
#                           "HLA gene","MHC region")
# 
#     # print data table
#     DT::datatable(gtltbl, rownames=F, options = list(
#       order = list(list(5, 'asc')))) %>% DT::formatStyle(
#         "TWAS significant", target = "row",
#         backgroundColor = styleEqual(c(1), c('lightgreen'))
#       ) %>% DT::formatStyle(
#         "HLA gene", target = "row",
#         backgroundColor = styleEqual(c(1), c('red'))
#       ) %>% DT::formatStyle(
#         "MHC region", target = "row",
#         backgroundColor = styleEqual(c(1), c('red'))
#       )
#   })
# 
# 
#   #####################
# 
# 
#   output$GWBtbltrt <- DT::renderDataTable({
# 
#     # select all genes at locus
#     gwbtbl <- twas_ds %>% filter(tissue=='GWB',genestart>=xlow() & genestop<=xhigh()
#                                   & pheno!=locpheno() & chr==locchr() & phenocat==locphcat()) %>%
#       select(genename,chr, pheno,phenocat,se_beta_, p, SignifGene,HLARegion,MHCRegion)
#     colnames(gwbtbl) <- c("Gene","Chr","Trait","Trait category","Beta","P-value","TWAS significant",
#                           "HLA gene","MHC region")
# 
#     # print data table
#     DT::datatable(gwbtbl, rownames=F, options = list(
#       order = list(list(5, 'asc')))) %>% formatStyle(
#         "TWAS significant", target = "row",
#         backgroundColor = styleEqual(c(1), c('lightgreen'))
#       ) %>% formatStyle(
#         "HLA gene", target = "row",
#         backgroundColor = styleEqual(c(1), c('red'))
#       ) %>% formatStyle(
#         "MHC region", target = "row",
#         backgroundColor = styleEqual(c(1), c('red'))
#       )
#   })
# 
# 
#   #####################
# 
# 
#   output$MSAtbltrt <- DT::renderDataTable({
# 
#     # select all genes at locus
#     msatbl <- twas_ds %>% filter(tissue=='MSA',genestart>=xlow() & genestop<=xhigh()
#                                   & pheno!=locpheno() & chr==locchr() & phenocat==locphcat()) %>%
#       select(genename,chr, pheno,phenocat,se_beta_, p, SignifGene,HLARegion,MHCRegion)
#     colnames(msatbl) <- c("Gene","Chr","Trait","Trait category","Beta","P-value","TWAS significant",
#                           "HLA gene","MHC region")
# 
#     # print data table
#     DT::datatable(msatbl, rownames=F, options = list(
#       order = list(list(5, 'asc')))) %>% formatStyle(
#         "TWAS significant", target = "row",
#         backgroundColor = styleEqual(c(1), c('lightgreen'))
#       ) %>% formatStyle(
#         "HLA gene", target = "row",
#         backgroundColor = styleEqual(c(1), c('red'))
#       ) %>% formatStyle(
#         "MHC region", target = "row",
#         backgroundColor = styleEqual(c(1), c('red'))
#       )
#   })
# 
# 
#   ############################################################
# 
#   # Table of known variants
#   output$KnownSNPtbl <- DT::renderDataTable({
#     locds <- dgnds %>% filter(locvar==input$locuslst)
#     xhigh <- max(locds$genestop)+1000000
#     xlow <- max(0,min(locds$genestart)-1000000)
#     locpheno<-unique(locds$phenocat) #locus phenotype category
#     locph<-unique(locds$pheno) #locus phenotype
#     locchr<- unique(locds$chr) #locus chromosome
# 
#     if (locpheno=="RBC") {
#       ds=rbcds
#     } else if (locpheno=="WBC") {
#       ds=wbcds
#     } else {
#       ds=pltds
#     }
# 
#     # select all known snps at locus
#     phenotbl <- ds %>% filter(X6>=xlow & xhigh>=X6 & X5==locchr) %>%
#       select(X1,X2,X3,X4,X5,X6,X8,X9,X10,X11,X12,X13,X14,X15,X21)
#     colnames(phenotbl) <- c("RSID", "Reference", "Ancestry", "Trait", "Chr", "Pos_b37", "Gene", "Effect Allele",
#                             "Other Allele", "EAF", "Beta", "SE", "P", "N","In Kaiser")
# 
# 	if (nrow(phenotbl)==0){
#       DT::datatable(phenotbl,
#                 options=list(columnDefs = list(list(visible=FALSE, targets=c(15,16,17,18,19,20,21,22,23)))))
#     } else {
# 
# 		# select set of TWAS insignificant genes
# 		dgntbl <- dgnds %>% filter(genestart>=xlow & genestop<=xhigh & pheno==locph & chr==locchr &
# 									SignifGene==0  & is.na(HLARegion) & is.na(MHCRegion)) %>% select(genename)
# 
# 		# indicator for significant TWAS gene pattern
# 		phenotbl$X16 <- ifelse(grepl(paste(unique(locds$genename), collapse="|"),phenotbl$Gene),1,0)
# 
# 		# indicator for non-significant TWAS gene pattern
# 		phenotbl$X17 <- ifelse(grepl(paste(unique(dgntbl$genename),collapse="|"),phenotbl$Gene),1,0)
# 
# 		# indicator for gene not predicted in Kaiser
# 		phenotbl$notpred <- ifelse(grepl(paste(genes$genename[genes$X0==0],collapse = "|"),phenotbl$Gene),1,0)
# 		phenotbl$notdgn <- ifelse(grepl(paste(genes$genename[genes$X0==1],collapse="|"),phenotbl$Gene),1,0)
# 
# 		phenotbl$chrposall <- paste0(phenotbl$Chr,":",phenotbl$Pos_b37,":",phenotbl$`Effect Allele`,":",phenotbl$`Other Allele`)
# 		phenotbl$chrposall2 <- paste0(phenotbl$Chr,":",phenotbl$Pos_b37,":",phenotbl$`Other Allele`,":",phenotbl$`Effect Allele`)
# 
# 		# indicator for snps included in GWAS results
# 		locgwas <- gwasdsfin %>% filter(Locus==locnum()) %>% select(SNP)
# 		phenotbl$X19 <- ifelse(grepl(paste(unique(locgwas$SNP),collapse="|"),phenotbl$chrposall),1,0)
# 		phenotbl$X20 <- ifelse(grepl(paste(unique(locgwas$SNP),collapse="|"),phenotbl$chrposall2),1,0)
# 
# 		# print data table
# 		DT::datatable(phenotbl, #rownames=FALSE,
# 				  options=list(columnDefs = list(list(visible=FALSE,
# 													  targets=c(15,16,17,18,19,20,21,22,23))))
# 				  ) %>%
# 		  formatStyle("Gene", "X16", backgroundColor = styleEqual(
# 			c(1), c("lightgreen"))
# 		  ) %>%
# 		  formatStyle("In Kaiser", target = "row",
# 			backgroundColor = styleEqual(c(0), c('red'))
# 		  ) %>%
# 		  formatStyle("Gene", "X17", backgroundColor = styleEqual(
# 			c(1), c("yellow"))
# 		  ) %>%
# 		  formatStyle("RSID", "X19", backgroundColor = styleEqual(
# 			c(1), c("lightgreen"))
# 		  ) %>%
# 		  formatStyle("RSID", "X20", backgroundColor = styleEqual(
# 			c(1), c("lightgreen"))
# 		  ) %>%
# 		  formatStyle("Gene", "notpred", backgroundColor = styleEqual(
# 			c(1), c("orange"))
# 		  ) %>%
# 		  formatStyle("Gene", "notdgn", backgroundColor = styleEqual(
# 			c(1), c("red"))
# 		  )
# 	}
#   })
# 
#   ############################################################
# 
#   # Table of known variants
#   output$GWASvars <- DT::renderDataTable({
#     locgwas <- gwasdsfin %>% filter(Locus==locnum())
# 
#     # print data table
#     DT::datatable(locgwas, #rownames=F,
#               options=list(columnDefs = list(list(visible=FALSE, targets=c(2,3,4,5,9,10,11)))))
#   })
# 
# }
  
  
  
  # give instructions to the server
  server <- function(input,output){
    
    # define reactive variables for all plots
    locds <- reactive({
      filter(dgnds, locvar==input$locuslst)
    })
    
    # get the locus number from the detailed drop down menu
    locnum <- reactive({ 
      tmp <- input$locuslst
      tmp2 <- as.numeric(strsplit(strsplit(tmp,":")[[1]][1]," ")[[1]][2])
      tmp2
    })
    
    # define locus window
    xhigh <- reactive({ unique(locds()$locstop) })
    xlow <- reactive({ unique(locds()$locstart) })
    xhighMB <- reactive({ round(xhigh()/1000000,4) })
    xlowMB <- reactive({ round(xlow()/1000000,4) })
    
    # locus specific characteristics
    locchr <-   reactive({ unique(locds()$chr) })
    locpheno <- reactive({ unique(locds()$phenoname) })
    locphcat <- reactive({ unique(locds()$phenocat) })
    loctitle <- reactive({ paste0("chr ",locchr(),": ",xlow(),", ",xhigh(),"; trait = ",locpheno(),
                                  ", trait category = ",locphcat())
    })
    
    # set the dataset to extract known variants at the locus
    ds <- reactive({
      if (locphcat()=="RBC") {
        rbcds
      } else if (locphcat()=="WBC") {
        wbcds
      } else {
        pltds
      }
    })
    
    # select all known GWAS genes at locus
    phenotbl <- reactive({
      tmp <- filter(ds(),(X6)>=xlow() & xhigh()>=(X6) & X5==locchr()) %>% select(X8)
      colnames(tmp) <- c("Gene")
      tmp
    })
    
    # select all TWAS genes at locus
    dgntbl <-  reactive({
      tmp <- filter(dgnds,genestart>=xlow() & genestop<=xhigh() & phenoname==locpheno() & chr==locchr())
      
      if (nrow(phenotbl())==0){
        tmp$kngene="Not Reported in GWAS"
      } else{
        tmp$kngene <- ifelse(grepl(paste(unique(phenotbl()$Gene), collapse="|"),tmp$genename),
                             "Reported in GWAS","Not Reported in GWAS")
      }
      tmp
    })
    
    
    ############################################################
    
    # uncomment when memory available
    # ideogram plot
    # output$chrplt <- renderPlot({
    #   
    #   itrack = IdeogramTrack(genome = "hg19", chromosome = paste0("chr",locchr()), bands = data)
    #   trackplot = c(itrack,GenomeAxisTrack());
    #   
    #   plotTracks(c(trackplot), from = xlow(),to=xhigh(),transcriptAnnotation="symbol",
    #              add53=TRUE,showBandID=TRUE,cex.bands=2,stackHeight=2,background.title = "white",
    #              col.axis="black",col.title="black",cex.title=2,
    #              cex.axis=4,just.group="below",cex.main = 1.5,cex=1.5)
    # })
    
    ############################################################
    
    
    # TWAS/GWAS mirror plot
    output$TWASmirror <- renderPlotly({
      locds <- dgnds %>% filter(locvar==input$locuslst)
      #xhigh <- max(locds$genestop)+1000000
      #xlow <- max(0,min(locds$genestart)-1000000)
      yhigh <- max(locds$log10pval)+0.25*max(locds$log10pval)
      pthresh <- -log10(0.05/(nrow(dgnds)))
      locchr <- unique(locds$chr)
      #locpheno <- unique(locds$phenoname)
      #locphcat <- unique(locds$phenocat)
      #loctitle <- paste0("chr ",locchr,": ",xlow,", ",xhigh,"; phenotype = ",locpheno,"; category = ",locphcat)
      
      # select known variants at locus
      if (locphcat()=="RBC") {
        ds=rbcds
      } else if (locphcat()=="WBC") {
        ds=wbcds
      } else {
        ds=pltds
      }
      
      # select all known GWAS genes at locus
      #phenotbl <- ds %>% filter(X6>=xlow() & xhigh()>=X6 & X5==locchr) %>% select(X8)
      #colnames(phenotbl) <- c("Gene")
      
      # select all TWAS genes at locus
      #dgntbl <- dgnds %>% filter(genestart>=xlow() & genestop<=xhigh() & phenoname==locpheno() & chr==locchr)
      #dgntbl$kngene <- ifelse(grepl(paste(unique(phenotbl$Gene), collapse="|"),dgntbl$genename),"Reported in GWAS","Not reported in GWAS")
      
      # select significant and known genes to plot
      dgntblplt <- dgntbl() %>% filter(SignifGene==1 | kngene=="Reported in GWAS")
      
      # select the GWAS results at the locus
      #gwasloc <- gwasallfin %>% filter(chr==locchr & as.numeric(pos)>=xlow() & xhigh()>=as.numeric(pos) & PHENO==locpheno)
      gwasloc <- gwasallfin %>% filter(Locus==locnum())
      #gwasloc$knsnp <- ifelse(grepl(paste(unique(gwasdsfin$SNP[gwasdsfin$Locus==locnum()]),collapse="|"),
      #                              gwasloc$SNP),"Reported in GWAS","Not reported in GWAS")
      
      if (nrow(gwasdsfin[gwasdsfin$Locus==locnum(),])==0){
        gwasloc$knsnp <- "Not reported in GWAS"
      } else {
        gwasloc$knsnp <- ifelse(grepl(paste(unique(gwasdsfin$SNP[gwasdsfin$Locus==locnum()]),collapse="|"),
                                      gwasloc$SNP),"Reported in GWAS","Not reported in GWAS")
      }
      
      # not known variants at the locus pval < 0.1
      gwaslocnotkn <- gwasloc %>% filter(knsnp=="Not reported in GWAS")
      
      # known variants at the locus, regardless of pvalue
      knowngwas <- gwasdsfin %>% filter(Locus==locnum())
      
      if (nrow(knowngwas)>0){
        knowngwas$knsnp <- "Reported in GWAS"
        ylow <- min(gwaslocnotkn$poslog10p,log10(5*10^(-8)),knowngwas$poslog10p)+0.15*min(gwaslocnotkn$poslog10p,log10(5*10^(-8)),knowngwas$poslog10p)
      } else{
        ylow <- min(gwaslocnotkn$poslog10p,log10(5*10^(-8)))+0.15*min(gwaslocnotkn$poslog10p,log10(5*10^(-8)))
      }
      
      nudgeval <- 0.05*max(abs(yhigh))
      nudgevalg <- 0.05*ylow
      
      # define specific colors for known/not known genes
      if(length(unique(dgntblplt$kngene))==1){
        if(unique(dgntblplt$kngene)=="Not reported in GWAS"){
          myColors <- setNames( c('#56B4E9','#000000'),
                                c("Not reported in GWAS","Reported in GWAS") )
        } else if(unique(dgntblplt$kngene)=="Reported in GWAS"){
          myColors <- setNames( c('#000000','#56B4E9'),
                                c("Reported in GWAS","Not reported in GWAS") )
        } else {
          myColors <- setNames( c('#56B4E9', '#000000', '#56B4E9'),
                                c("Not reported in GWAS","Reported in GWAS",NA))  
        }
      } else {
        myColors <- setNames( c('#56B4E9', '#000000', '#56B4E9'),
                              c("Not reported in GWAS","Reported in GWAS",NA))
      }
      
      colScale <- scale_colour_manual(values = myColors)
      
      # color scale for known/not known gwas variants
      myColorsg <-
        setNames( c('#56B4E9', '#000000', '#E69F00')
                  , c("Reported in GWAS","Not reported in GWAS",NA))
      colScaleg <- scale_colour_manual(values=myColorsg)
      
      if (input$margcondradio==1){
        pvalplt <- dgntblplt$log10pval
      }
      else if(input$margcondradio==2){
        pvalplt <- -log10(dgntblplt$p_final)
      }
      
      # TWAS plot
      p <- ggplotly(ggplot(data=dgntblplt, aes(x=round(genemid/1000000,4), y=pvalplt, color=as.factor(kngene))) + 
                      geom_point(pch=15) +
                      geom_segment(aes(x = genestartMB, y = pvalplt, xend = genestopMB, yend = pvalplt, color=as.factor(kngene)), size=2) +
                      geom_text(aes(x=round(genemid/1000000,4),y=pvalplt,label=genename,color=as.factor(kngene)), nudge_y = nudgeval, size=4) +
                      geom_hline(aes(yintercept=pthresh), lty=2, color="red") +
                      xlim(xlowMB(),xhighMB()) + ylim(0,yhigh) +
                      annotate(geom="text",x=round((xlow()+270000)/1000000,4), y=pthresh-nudgeval, color='red',
                               label=paste0("TWAS p-value: ",formatC(10^(-pthresh), format = "e", digits = 2)),size=4) +
                      theme(legend.position = 'top', legend.title = element_blank()) +
                      theme_bw() +
                      colScale
      )
      
      p <- p %>% plotly::layout(title = list(text=loctitle(), x=0, xanchor='left', y=.99),
                        yaxis = list(title = ' TWAS -log10(p)'),
                        legend = list(orientation='h', x=0, y=1))
      
      
      plta <- ggplot(data=gwaslocnotkn, aes(x=round(as.numeric(pos)/1000000,4),y=poslog10p, color=as.factor(knsnp))) + 
        geom_point() +
        geom_hline(aes(yintercept=(log10(5*10^(-8)))), lty=2, color="red") +
        xlim(xlowMB(),xhighMB()) + ylim(ylow,0) +
        theme_bw() +
        colScaleg  +
        annotate(geom="text",x=xlowMB()+.270000, y=log10(5*10^(-8))+nudgevalg, color='red',
                 label=paste0("GWAS p-value: ",formatC(5*10^(-8), format = "e", digits = 2)),size=4)
      
      
      if (nrow(phenotbl())==0 | nrow(knowngwas)==0){
        app_plt <- plta
      } else {
        app_plt <- plta + 
          geom_point(data=knowngwas, aes(x=round(as.numeric(pos)/1000000,4),y=poslog10p, color=as.factor(knsnp)))
      }
      
      # GWAS plot
      m <- ggplotly(app_plt)
      
      m <- m %>% plotly::layout(xaxis = list(title="position (in Mb)"), 
                        yaxis=list(title="GWAS log10(p)"))
      
      # combine plots together
      fig <- subplot(style(p,showlegend=F), m, shareX = T, nrows = 2, titleX = T,titleY = T, 
                     which_layout = 1, margin=0.005)
      
    })
    
    ############################################################
    
    
    # TWAS correlation plots
    output$TWAScorr <- renderPlotly({
      
      # select significant and known genes to plot
      dgntblplt <- dgntbl() %>% filter(SignifGene==1 | kngene=="Reported in GWAS")
      uniqgenes <- dgntblplt %>% select(gene)
      indexgene <- dgntblplt$gene[dgntblplt$p == min(dgntblplt$p)]
      
      # extract the set of genes from expression file
      # exprloc <- expds %>% select(uniqgenes$gene)
      # colnumindex <- as.numeric(which(colnames(exprloc)==indexgene))
      
      # extract set of genes from correlation matrix
      Mindex <- data.table(M[uniqgenes$gene,indexgene])
      colnames(Mindex) <- c("corr")
      Mindex$gene <- uniqgenes$gene
      
      # calculate the correlation between pairs of genes
      #M <- cor(exprloc)
      
      # select index gene
      # Mindex <- data.table(M[,colnumindex])
      # colnames(Mindex) <- c("corr")
      # Mindex$gene <- rownames(M)
      
      # categorize the correlation values into bins for plotting
      Mindex[,corgroup:= cut(Mindex[,abs(corr)],
                             c(0,.2,.4,.6,.8,1),
                             labels=c("[0 - 0.2]","(0.2 - 0.4]","(0.4 - 0.6]","(0.6 - 0.8]","(0.8 - 1]"),
                             right=T,include.lowest=T)]
      
      # color for correlation categories
      Mindex[,corcol:= cut(Mindex[,abs(corr)],
                           c(0,.2,.4,.6,.8,1),
                           labels=c("#CCCCCC","#333FFF","#00FF00","#FF9900","#FF0000"),
                           right=T,include.lowest = T)]
      
      #Mindex[,corgroup.raw:= cut(Mindex[,abs(corr)],c(0,.2,.4,.6,.8,1),right=T,include.lowest=T)]
      #Mindex[,table(corgroup,corgroup.raw)]
      
      # match correlation values and categories back to overall table
      corplt <- merge(dgntblplt,Mindex,by.x="gene",by.y="gene")
      
      # get the DGN weights for the specific genes at locus
      dgnwtloc <- weight_ds %>% filter(gene %in% corplt$gene)
      
      # subset the GWAS variants for specific chr, and phenotype
      #gwasallfinloc <- gwasallfin %>% filter(chr==locchr(),PHENO==locpheno()) %>%
      #  select(pos,all1,all2,poslog10p)
      gwasallfinloc <- gwasallfin %>% filter(Locus==locnum()) %>%
        select(pos,all1,all2,poslog10p)
      gwasallfinloc$posnum <- as.numeric(gwasallfinloc$pos)
      
      # get the matching GWAS variants for the weights
      dgnwtgwasloc1 <- merge(dgnwtloc,gwasallfinloc,by.x=c('V7','ref_allele','eff_allele'),
                             by.y=c('posnum','all1','all2'))
      dgnwtgwasloc2 <- merge(dgnwtloc,gwasallfinloc,by.x=c('V7','ref_allele','eff_allele'),
                             by.y=c('posnum','all2','all1'))
      dgnwtgwasloc <- rbind(dgnwtgwasloc1,dgnwtgwasloc2)
      
      # merge GWAS results with TWAS info
      TWASloc <- corplt %>% select(gene,genename,genestartMB,genemid,genestopMB,log10pval,corgroup)
      TWASloc$genemidMB <- round(TWASloc$genemid/1000000,4)
      dgnwtgwaslocfin <- merge(dgnwtgwasloc,TWASloc, by.x='gene',by.y='gene',all.x=T)
      
      # set y limit
      yhigh <- max(corplt$log10pval)+0.25*max(corplt$log10pval)
      ylow <- min(dgnwtgwasloc$poslog10p,log10(5*10^(-8)))+0.15*min(dgnwtgwasloc$poslog10p,log10(5*10^(-8)))
      nudgeval <- 0.05*max(abs(yhigh))
      
      # set color scale values
      colval <- c("[0 - 0.2]"="#CCCCCC","(0.2 - 0.4]"="#333FFF",
                  "(0.4 - 0.6]"="#00FF00","(0.6 - 0.8]"="#FF9900","(0.8 - 1]"="#FF0000")
      breakval <- c("[0 - 0.2]","(0.2 - 0.4]","(0.4 - 0.6]","(0.6 - 0.8]","(0.8 - 1]")
      
      
      # select the LD for the given locus
      LDdsloc <- as.data.table(filter(LDds,locus2==locnum()))
      
      # define correlation and color groups
      LDdsloc[,ldcol:= cut(LDdsloc[,abs(corrab)],
                           c(0,.2,.4,.6,.8,1),
                           labels=c("#CCCCCC","#333FFF","#00FF00","#FF9900","#FF0000"),
                           right=T,include.lowest = T)]
      
      LDdsloc[,ldgroup:= cut(LDdsloc[,abs(corrab)],
                             c(0,.2,.4,.6,.8,1),
                             labels=c("[0 - 0.2]","(0.2 - 0.4]","(0.4 - 0.6]","(0.6 - 0.8]","(0.8 - 1]"),
                             right=T,include.lowest=T)]
      
      
      
      # match snp LD back with overall file
      dgnwtgwaslocfin2 <- merge(dgnwtgwaslocfin,LDdsloc, by.x=c('V7'),by.y=c('posb'),all.x=T)
      
      dgnwtgwaslocfin2_comp <- dgnwtgwaslocfin2[complete.cases(dgnwtgwaslocfin2), ]
      dgnwtgwaslocfin2_comp$line <- ifelse(dgnwtgwaslocfin2_comp$weight<0,2,1)
      
      xlowMB <- min(corplt$genestartMB,round(dgnwtgwaslocfin2_comp$V7/1000000,4))-.0005
      xhighMB <- max(corplt$genestopMB,round(dgnwtgwaslocfin2_comp$V7/1000000,4))+.0005
      
      #cntuniq <- length(unique(corplt$corgroup,dgnwtgwaslocfin2_comp$ldgroup))
      
      # plot correlation categories at the locus, like locus zoom plot
      loczoom <- ggplotly(ggplot(data=corplt, aes(x=round(genemid/1000000,4), y=log10pval)) +
                            geom_hline(aes(yintercept=pthresh), lty=2, color="red") +
                            geom_hline(aes(yintercept=log10(5*10^(-8))), lty=2, color="red") +
                            geom_hline(aes(yintercept=0),size=1,color='black') +
                            geom_segment(data=dgnwtgwaslocfin2_comp,aes(x=round(V7/1000000,4),y=poslog10p,xend=genemidMB,yend=log10pval
                                                                        ,color=ldgroup)) +
                            geom_point(aes(color=corgroup), pch=15) +
                            geom_segment(aes(x=genestartMB,y=log10pval,xend=genestopMB,yend=log10pval, color=corgroup),size=2) +
                            geom_text(aes(label=genename,color=corgroup),nudge_y = nudgeval) +
                            scale_color_manual(breaks = breakval, values = colval) +
                            # scale_color_manual(name = NULL,
                            #    values = colval,
                            #    breaks = breakval,
                            #    guide = guide_legend(override.aes = list(shape = rep(15, 5),
                            #                                             color = c("#CCCCCC","#333FFF",
                            #                                                       "#00FF00","#FF9900","#FF0000")) ) ) +
                            theme_bw() +
                            #xlim(xlowMB(),xhighMB()) +
                            xlim(xlowMB,xhighMB) +
                            ylim(ylow,yhigh) +
                            theme(legend.position = 'top', legend.title = element_blank()) +
                            #annotate(geom="text",x=round((xlow()+270000)/1000000,4), y=pthresh-nudgeval, color='red',
                            #         label=paste0("TWAS p-value: ",formatC(10^(-pthresh), format = "e", digits = 2)),size=4) +
                            geom_point(data=dgnwtgwaslocfin2_comp, aes(x=round(V7/1000000,4),y=poslog10p,color=ldgroup)) #+
      )
      
      loczoom <- loczoom %>% plotly::layout(title = list(text=loctitle(), x=0, xanchor='left', y=.99),
                                    yaxis = list(title = '<-- GWAS log10(p) ... 0 ... TWAS -log10(p) -->'),
                                    legend = list(orientation='h', x=0, y=1),
                                    xaxis = list(title="position (in Mb)"))
    })
    
    
    ############################################################
    
    
    # TWAS network plots
    output$TWASnetwork <- renderVisNetwork({
      
      # select significant and known genes to plot
      dgntblplt <- dgntbl() %>% filter(SignifGene==1 | kngene=="Reported in GWAS")
      uniqgenes <- dgntblplt %>% select(gene)
      indexgene <- dgntblplt$gene[dgntblplt$p == min(dgntblplt$p)]
      
      # extract set of genes from correlation matrix
      Mindex <- data.table(M[uniqgenes$gene,indexgene])
      colnames(Mindex) <- c("corr")
      Mindex$gene <- uniqgenes$gene
      
      # categorize the correlation values into bins for plotting
      Mindex[,corgroup:= cut(Mindex[,abs(corr)],
                             c(0,.2,.4,.6,.8,1),
                             labels=c("[0 - 0.2]","(0.2 - 0.4]","(0.4 - 0.6]","(0.6 - 0.8]","(0.8 - 1]"),
                             right=T,include.lowest=T)]
      
      # color for correlation categories
      Mindex[,corcol:= cut(Mindex[,abs(corr)],
                           c(0,.2,.4,.6,.8,1),
                           labels=c("#CCCCCC","#333FFF","#00FF00","#FF9900","#FF0000"),
                           right=T,include.lowest = T)]
      
      # match correlation values and categories back to overall table
      corplt <- merge(dgntblplt,Mindex,by.x="gene",by.y="gene")
      
      # get the DGN weights for the specific genes at locus
      dgnwtloc <- weight_ds %>% filter(gene %in% corplt$gene)
      
      # subset the GWAS variants for specific chr, and phenotype
      #gwasallfinloc <- gwasallfin %>% filter(chr==locchr(),PHENO==locpheno()) %>%
      #  select(pos,all1,all2,poslog10p)
      gwasallfinloc <- gwasallfin %>% filter(Locus==locnum()) %>%
        select(pos,all1,all2,poslog10p)
      gwasallfinloc$posnum <- as.numeric(gwasallfinloc$pos)
      
      # get the matching GWAS variants for the weights
      dgnwtgwasloc1 <- merge(dgnwtloc,gwasallfinloc,by.x=c('V7','ref_allele','eff_allele'),
                             by.y=c('posnum','all1','all2'))
      dgnwtgwasloc2 <- merge(dgnwtloc,gwasallfinloc,by.x=c('V7','ref_allele','eff_allele'),
                             by.y=c('posnum','all2','all1'))
      dgnwtgwasloc <- rbind(dgnwtgwasloc1,dgnwtgwasloc2)
      
      # merge GWAS results with TWAS info
      TWASloc <- corplt %>% select(gene,genename,genestartMB,genemid,genestopMB,log10pval,corgroup)
      TWASloc$genemidMB <- round(TWASloc$genemid/1000000,4)
      dgnwtgwaslocfin <- merge(dgnwtgwasloc,TWASloc, by.x='gene',by.y='gene',all.x=T)
      
      # set y limit
      yhigh <- max(corplt$log10pval)+0.15*max(corplt$log10pval)
      ylow <- min(dgnwtgwasloc$poslog10p)+0.15*min(dgnwtgwasloc$poslog10p)
      nudgeval <- 0.05*max(abs(yhigh))
      
      # set color scale values
      colval <- c("[0 - 0.2]"="#CCCCCC","(0.2 - 0.4]"="#333FFF",
                  "(0.4 - 0.6]"="#00FF00","(0.6 - 0.8]"="#FF9900","(0.8 - 1]"="#FF0000")
      
      # select the LD for the given locus
      LDdsloc <- as.data.table(filter(LDds,locus2==locnum()))
      
      # define correlation and color groups
      LDdsloc[,ldcol:= cut(LDdsloc[,abs(corrab)],
                           c(0,.2,.4,.6,.8,1),
                           labels=c("#CCCCCC","#333FFF","#00FF00","#FF9900","#FF0000"),
                           right=T,include.lowest = T)]
      
      LDdsloc[,ldgroup:= cut(LDdsloc[,abs(corrab)],
                             c(0,.2,.4,.6,.8,1),
                             labels=c("[0 - 0.2]","(0.2 - 0.4]","(0.4 - 0.6]","(0.6 - 0.8]","(0.8 - 1]"),
                             right=T,include.lowest=T)]
      
      # match snp LD back with overall file
      dgnwtgwaslocfin2 <- merge(dgnwtgwaslocfin,LDdsloc, by.x=c('V7'),by.y=c('posb'),all.x=T)
      
      dgnwtgwaslocfin2_comp <- dgnwtgwaslocfin2[complete.cases(dgnwtgwaslocfin2), ]
      dgnwtgwaslocfin2_comp$line <- ifelse(dgnwtgwaslocfin2_comp$weight<0,2,1)
      
      
      # data for links
      from <- dgnwtgwaslocfin2$rsid
      to <- dgnwtgwaslocfin2$genename
      weight <- dgnwtgwaslocfin2$weight
      color <- dgnwtgwaslocfin2$ldcol
      length <- 1/(dgnwtgwaslocfin2$log10pval-dgnwtgwaslocfin2$poslog10p)*1000
      links <- data.frame(from,to,weight,color,length)
      links$dashes <- ifelse(links$weight<0,TRUE,FALSE)
      links$width <- abs(links$weight)*100
      
      links_complete <- links[complete.cases(links), ]
      
      
      # data for gene nodes
      nodesgene <- corplt %>% select(genename,log10pval,corcol)
      colnames(nodesgene) <- c("id","size","color.background")
      nodesgene$title=paste0(nodesgene$id,", -log10(p)=",round(nodesgene$size,2))
      nodesgene$type="gene"
      nodesgene$shape="square"
      nodesgene$shadow=TRUE
      nodesgene$color.border="black"
      nodesgene$label=nodesgene$id
      maxp <- max(nodesgene$size)
      nodesgene$shape[nodesgene$size==maxp]="star"
      # nodesgene$Gene = nodesgene$label
      # topgene <- nodesgene$Gene[nodesgene$size==maxp]
      # nodesgene$Gene[nodesgene$size==maxp]=paste0(topgene,"*")
      
      
      nodessnp <- dgnwtgwaslocfin2 %>% select(rsid,poslog10p,ldcol)
      colnames(nodessnp) <- c("id","poslog10p","color.background")
      nodessnp$size <- abs(nodessnp$poslog10p)
      nodessnp$title=paste0(nodessnp$id,", -log10(p)=",round(abs(nodessnp$poslog10p),2))
      nodessnp$type="snp"
      nodessnp$shape="dot"
      nodessnp$shadow=FALSE
      nodessnp$color.border="black"
      nodessnp$label=""
      nodessnp$shape[nodessnp$size==max(nodessnp$size)]="triangle"
      #nodessnp$Gene[nodessnp$Gene==topgene]=paste0(topgene,"*")
      
      
      nodessnp <- select(nodessnp,-poslog10p)
      nodessnpfin <- distinct(nodessnp)
      
      nodes <- rbind(nodesgene,nodessnpfin)
      nodes_complete <- nodes[complete.cases(nodes), ]
      
      nodes_complete$color.highlight.background <- "purple"
      nodes_complete$color.highlight.border <- "purple"
      #nodes_complete$font.size <- max(nodessnp$size, nodesgene$size)*2
      nodes_complete$font.size <- 20
      
      visnet <- visNetwork(nodes_complete,links_complete, main="")
      visOptions(visnet, highlightNearest = TRUE)
      
      
    })
    
    
    ############################################################
    
    # TWAS/GWAS mirror plot of meta-analysis
    output$Tmetamirror <- renderPlotly({
      
      # select significant and known genes to plot
      dgntblplt <- dgntbl() %>% filter(SignifGene==1 | kngene=="Reported in GWAS")
      
      
      # define y limits and nudges for plotting
      yhigh <- max(locds()$log10pval)+0.25*max(locds()$log10pval)
      ylow <- min(-dgntblplt$log10pvalmeta,metathresh)+0.15*min(-dgntblplt$log10pvalmeta,metathresh)
      nudgeval <- 0.05*max(abs(yhigh))
      nudgevalg <- 0.05*ylow
      
      
      # define specific colors for known/not known genes
      myColors <- setNames( c('#56B4E9', '#000000', '#E69F00'),
                            c("Not reported in GWAS","Reported in GWAS",NA))
      
      colScale <- scale_colour_manual(values = myColors)
      
      
      # TWAS plot
      p <- ggplotly(ggplot(data=dgntblplt, aes(x=round(genemid/1000000,4), y=log10pval, color=as.factor(kngene))) +
                      geom_point(pch=15) +
                      geom_segment(aes(x = genestartMB, y = log10pval, xend = genestopMB, yend = log10pval,
                                       color=as.factor(kngene)), size=2) +
                      geom_text(aes(x=round(genemid/1000000,4),y=log10pval,label=genename,color=as.factor(kngene)),
                                nudge_y = nudgeval, size=4) +
                      geom_hline(aes(yintercept=pthresh), lty=2, color="red") +
                      xlim(xlowMB(),xhighMB()) + ylim(0,yhigh) +
                      #xlim(xlow,xhigh) + ylim(0,yhigh) +
                      annotate(geom="text",x=xlowMB()+.270000, y=pthresh-nudgeval, color='red',
                               label=paste0("TWAS p-value: ",formatC(10^(-pthresh), format = "e", digits = 2)),size=4) +
                      theme(legend.position = 'top', legend.title = element_blank()) +
                      theme_bw() +
                      colScale
      )
      
      p <- p %>% plotly::layout(title = list(text=loctitle(), x=0, xanchor='left', y=.99),
                        yaxis = list(title = ' TWAS -log10(p)'),
                        legend = list(orientation='h', x=0, y=1))
      
      # Meta-analysis plot
      m <- ggplotly(ggplot(data=dgntblplt, aes(x=round(genemid/1000000,4), y=-log10pvalmeta, color=as.factor(kngene))) +
                      geom_point(pch=15) +
                      geom_segment(aes(x = genestartMB, y = -log10pvalmeta, xend = genestopMB, yend = -log10pvalmeta,
                                       color=as.factor(kngene)), size=2) +
                      geom_text(aes(x=round(genemid/1000000,4),y=(-log10pvalmeta),label=genename,color=as.factor(kngene)),
                                nudge_y = nudgevalg, size=4) +
                      geom_hline(aes(yintercept=log10(0.05)), lty=2, color="orange") +
                      geom_hline(aes(yintercept=metathresh), lty=2, color="red") +
                      xlim(xlowMB(),xhighMB()) + ylim(ylow,0) +
                      # #xlim(xlow,xhigh) + ylim(ylow,0) +
                      annotate(geom="text",x=xlowMB()+.270000, y=log10(0.05)+nudgevalg, color='orange',
                               label=paste0("Replication p-value: 0.05"),size=4) +
                      annotate(geom="text",x=xlowMB()+.370000, y=metathresh+nudgevalg, color='red',
                               label=paste0("Strict replication p-value: ", formatC(10^(metathresh), format = "e", digits = 2)),size=4) +
                      theme(legend.position = 'top', legend.title = element_blank()) +
                      theme_bw() +
                      colScale
      )
      
      m <- m %>% plotly::layout(xaxis = list(title="position (in Mb)"),
                        yaxis=list(title="Meta-analysis log10(p)"))
      
      # combine plots together
      fig <- subplot(style(p,showlegend=F), m, shareX = T, nrows = 2, titleX = T,titleY = T,
                     which_layout = 1, margin=0.005)
      
    })
    
    ############################################################
    
    
    
    # comparison of TWAS DGN results with results from other reference panels
    # output$CompRef <- renderPlotly({
    #   
    #   # select DGN results
    #   locds <- dgnds %>% filter(locus==input$locuslst)
    #   xhigh <- max(locds$genestop)+1000000
    #   xlow <- max(0,min(locds$genestart)-1000000)
    #   pthresh <- -log10(0.05/(nrow(dgnds)))
    #   locchr <- unique(locds$chr)
    #   locpheno <- unique(locds$pheno)
    #   locphcat <- unique(locds$phenocat)
    #   loctitle <- paste0("chr ",locchr,": ",xlow,", ",xhigh,"; phenotype = ",locpheno,"; category = ",locphcat)
    #   
    #   
    #   #####################
    #   
    #   
    #   # select reference panel specific results
    #   dgntbl <- dgnds %>% filter(genestart>=xlow & genestop<=xhigh & pheno==locpheno & chr==locchr & 
    #                                is.na(HLARegion) & is.na(MHCRegion)) %>%
    #     select(genename,chr, pheno,genestart,genestop,genemid,log10pval)
    #   
    #   locgwb <- kaiserds %>% filter(tissue=="GWB" & genestart>=xlow & genestop<=xhigh & chr==locchr & 
    #                                   pheno==locpheno & is.na(HLARegion) & is.na(MHCRegion)) %>%
    #     select(genename,chr, pheno,log10pval,genestart,genestop,genemid, SignifGene)
    #   pthreshgwb <- -log10(0.05/(nrow(kaiserds[kaiserds$tissue=="GWB",])))
    #   
    #   locgtl <- kaiserds %>% filter(tissue=="GTL" & genestart>=xlow & genestop<=xhigh & chr==locchr & 
    #                                   pheno==locpheno & is.na(HLARegion) & is.na(MHCRegion)) %>%
    #     select(genename,chr, pheno,log10pval,genestart,genestop,genemid, SignifGene)
    #   pthreshgtl <- -log10(0.05/(nrow(kaiserds[kaiserds$tissue=="GTL",])))
    #   
    #   locmsa <- kaiserds %>% filter(tissue=="MSA" & genestart>=xlow & genestop<=xhigh & chr==locchr & 
    #                                   pheno==locpheno & is.na(HLARegion) & is.na(MHCRegion)) %>%
    #     select(genename,chr, pheno,log10pval,genestart,genestop,genemid, SignifGene)
    #   pthreshmsa <- -log10(0.05/(nrow(kaiserds[kaiserds$tissue=="MSA",])))
    #   
    #   
    #   #####################
    #   
    #   
    #   # merge secondary reference panel data sets with DGN
    #   dgngwb <- merge(dgntbl,locgwb, by=c("genename","chr","pheno"), all = T)
    #   dgngwb$inboth <- ifelse(complete.cases(dgngwb$log10pval.x,dgngwb$log10pval.y),"In Both","Not in Both")
    #   
    #   dgngtl <- merge(dgntbl,locgtl, by=c("genename","chr","pheno"), all = T)
    #   dgngtl$inboth <- ifelse(complete.cases(dgngtl$log10pval.x,dgngtl$log10pval.y),"In Both","Not in Both")
    #   
    #   dgnmsa <- merge(dgntbl,locmsa, by=c("genename","chr","pheno"), all = T)
    #   dgnmsa$inboth <- ifelse(complete.cases(dgnmsa$log10pval.x,dgnmsa$log10pval.y),"In Both","Not in Both")
    #   
    #   
    #   #####################
    #   
    #   
    #   # set y limit values for each plot
    #   yhighdgn <- max(dgngwb$log10pval.x, na.rm = T)+0.15*max(dgngwb$log10pval.x, na.rm = T)
    #   nudgevaldgn <- 0.05*yhighdgn
    #   
    #   yhighgtl <- max(dgngtl$log10pval.y,pthreshgtl, na.rm = T)+0.15*max(dgngtl$log10pval.y,pthreshgtl, na.rm = T)
    #   nudgevalgtl <- 0.05*yhighgtl
    #   
    #   yhighgwb <- max(dgngwb$log10pval.y,pthreshgwb, na.rm = T)+0.15*max(dgngwb$log10pval.y,pthreshgwb, na.rm = T)
    #   nudgevalgwb <- 0.05*yhighgwb
    #   
    #   yhighmsa <- max(dgnmsa$log10pval.y,pthreshmsa, na.rm = T)+0.15*max(dgnmsa$log10pval.y,pthreshmsa, na.rm = T)
    #   nudgevalmsa <- 0.05*yhighmsa
    #   
    #   
    #   ####################
    #   
    #   # plot info for GWB
    #   atop <- ggplotly(ggplot(data=dgngwb, aes(x=genemid.x,y=log10pval.x, color=inboth)) +
    #     geom_point(pch=15) +
    #     geom_segment(aes(x=genestart.x,y=log10pval.x, xend=genestop.x, yend=log10pval.x), size=2) +
    #     geom_text(aes(label=genename), nudge_y = nudgevaldgn) +
    #     geom_hline(aes(yintercept=pthresh), lty=2, color='blue') +
    #       xlim(xlow,xhigh) +
    #       ylim(0,yhighdgn) +
    #       theme(legend.position = 'top', legend.title = element_blank())
    #   )
    #   
    #   atop <- atop %>% layout(yaxis = list(title = ' DGN TWAS -log10(p)'),
    #                           xaxis = list(range=c(xlow,xhigh)),
    #                           annotations=list(x = 0.5 , y = 1.1, text = "(a) DGN vs. GWB", showarrow = F, 
    #                                            xref='paper', yref='paper',xanchor='center'),
    #                           legend = list(orientation='h', x=0, y=1),
    #                           title=list(text=paste0(loctitle,"\n"),x=0,xanchor='left'), margin=list(t=60)
    #   )
    #   
    #   #####
    #   
    #   abottom <- ggplotly(ggplot(data=dgngwb, aes(x=genemid.y,y=-log10pval.y, color=inboth)) +
    #     geom_hline(aes(yintercept=-pthreshgwb), lty=2, color='blue') +
    #     geom_point(pch=15) +
    #     geom_segment(aes(x=genestart.y,y=-log10pval.y, xend=genestop.y, yend=-log10pval.y), size=2) +
    #     geom_text(aes(label=genename), nudge_y = -1.5*nudgevalgwb) +
    #     xlim(xlow,xhigh) +
    #     ylim(-yhighgwb,0) +
    #     theme(legend.position = 'top', legend.title = element_blank())
    #   )
    #   
    #   abottom <- abottom %>% layout(yaxis = list(title = 'GWB TWAS log10(p)', 
    #                                        range=c(-yhighgwb,0.1)),
    #                           xaxis = list(range=c(xlow,xhigh), title="position"))
    # 
    #   #####
    #   
    #   afig <- subplot(atop,style(abottom,showlegend=F),nrows=2,shareX = T, titleX = T, titleY = T, which_layout = 1, margin=.01)
    #   
    # 
    #   #####################
    #   
    #   
    #   # plot info for GTL
    #   btop <- ggplotly(ggplot(data=dgngtl, aes(x=genemid.x,y=log10pval.x, color=inboth)) +
    #                      geom_point(pch=15) +
    #                      geom_segment(aes(x=genestart.x,y=log10pval.x, xend=genestop.x, yend=log10pval.x), size=2) +
    #                      geom_text(aes(label=genename), nudge_y = nudgevaldgn) +
    #                      geom_hline(aes(yintercept=pthresh), lty=2, color='blue') +
    #                      xlim(xlow,xhigh) +
    #                      ylim(0,yhighdgn) +
    #                      theme(legend.position = 'top', legend.title = element_blank())
    #   )
    #   
    #   btop <- btop %>% layout(yaxis = list(title = ' DGN TWAS -log10(p)'),
    #                           xaxis = list(range=c(xlow,xhigh)),
    #                           annotations=list(x = 0.5 , y = 1.1, text = "(b) DGN vs. GTL", showarrow = F, 
    #                                            xref='paper', yref='paper',xanchor='center'),
    #                           legend = list(orientation='h', x=0, y=1))
    #   
    #   #####
    #   
    #   bbottom <- ggplotly(ggplot(data=dgngtl, aes(x=genemid.y,y=-log10pval.y, color=inboth)) +
    #                         geom_hline(aes(yintercept=-pthreshgtl), lty=2, color='blue') +
    #                         geom_point(pch=15) +
    #                         geom_segment(aes(x=genestart.y,y=-log10pval.y, xend=genestop.y, yend=-log10pval.y), size=2) +
    #                         geom_text(aes(label=genename), nudge_y = -1.5*nudgevalgtl) +
    #                         xlim(xlow,xhigh) +
    #                         ylim(-yhighgtl,0) +
    #                         theme(legend.position = 'top', legend.title = element_blank())
    #   )
    #   
    #   bbottom <- bbottom %>% layout(yaxis = list(title = 'GTL TWAS log10(p)', 
    #                                              range=c(-yhighgtl,0.1)),
    #                                 xaxis = list(range=c(xlow,xhigh), title="position"))
    #   
    #   #####
    #   
    #   bfig <- subplot(btop,style(bbottom,showlegend=F),nrows=2,shareX = T, titleX = T, titleY = T, which_layout = 1, margin=0.01)
    #   
    #   
    #   #####################
    #   
    #   
    #   # plot info for MSA
    #   ctop <- ggplotly(ggplot(data=dgnmsa, aes(x=genemid.x,y=log10pval.x, color=inboth)) +
    #                      geom_point(pch=15) +
    #                      geom_segment(aes(x=genestart.x,y=log10pval.x, xend=genestop.x, yend=log10pval.x), size=2) +
    #                      geom_text(aes(label=genename), nudge_y = nudgevaldgn) +
    #                      geom_hline(aes(yintercept=pthresh), lty=2, color='blue') +
    #                      xlim(xlow,xhigh) +
    #                      ylim(0,yhighdgn) +
    #                      theme(legend.position = 'top', legend.title = element_blank())
    #   )
    #   
    #   ctop <- ctop %>% layout(yaxis = list(title = 'DGN TWAS -log10(p)'),
    #                           xaxis = list(range=c(xlow,xhigh)),
    #                           annotations=list(x = 0.5 , y = 1.1, text = "(c) DGN vs. MSA", showarrow = F, 
    #                                            xref='paper', yref='paper',xanchor='center'),
    #                           legend = list(orientation='h', x=0, y=1))
    #   
    #   #####
    #   
    #   cbottom <- ggplotly(ggplot(data=dgnmsa, aes(x=genemid.y,y=-log10pval.y, color=inboth)) +
    #                         geom_hline(aes(yintercept=-pthreshmsa), lty=2, color='blue') +
    #                         geom_point(pch=15) +
    #                         geom_segment(aes(x=genestart.y,y=-log10pval.y, xend=genestop.y, yend=-log10pval.y), size=2) +
    #                         geom_text(aes(label=genename), nudge_y = -1.5*nudgevalmsa) +
    #                         xlim(xlow,xhigh) +
    #                         ylim(-yhighmsa-2,0)
    #   )
    #   
    #   cbottom <- cbottom %>% layout(yaxis = list(title = 'MSA TWAS log10(p)', 
    #                                              range=c(-yhighmsa,0.1)),
    #                                 xaxis = list(range=c(xlow,xhigh), title="position"))
    #   
    #   #####
    #   
    #   cfig <- subplot(ctop,style(cbottom,showlegend=F),nrows=2,shareX = T, titleX = T, titleY = T, which_layout = 1,
    #                   margin = 0.01)
    #   
    #   
    #   #####################
    #   
    #   
    #   # combine all figures
    #   allfig <- subplot(style(afig,showlegend=F),style(bfig,showlegend=F), cfig, 
    #                     nrows=1, shareY = F, titleX = T, titleY = T, which_layout = 1)
    # 
    # 
    # })
    
    
    output$CompRefGWB <- renderPlotly({
      
      # select DGN results
      locds <- dgnds %>% filter(locvar==input$locuslst)
      xhigh <- max(locds$genestop)+1000000
      xlow <- max(0,min(locds$genestart)-1000000)
      locchr <- unique(locds$chr)
      #locpheno <- unique(locds$pheno)
      #locphcat <- unique(locds$phenocat)
      #loctitle <- paste0("chr ",locchr,": ",xlow,", ",xhigh,"; phenotype = ",locpheno,"; category = ",locphcat)
      
      
      #####################
      
      
      # select reference panel specific results
      dgntbl <- dgnds %>% filter(genestart>=xlow & genestop<=xhigh & phenoname==locpheno() & chr==locchr & 
                                   is.na(HLARegion) & is.na(MHCRegion)) %>%
        select(genename,chr, phenoname,genestart,genestop,genemid,log10pval)
      
      locgwb <- twas_ds %>% filter(tissue=="GWB" & genestart>=xlow & genestop<=xhigh & chr==locchr & 
                                      phenoname==locpheno() & is.na(HLARegion) & is.na(MHCRegion)) %>%
        select(genename,chr, phenoname,log10pval,genestart,genestop,genemid, SignifGene)
      pthreshgwb <- -log10(0.05/(nrow(twas_ds[twas_ds$tissue=="GWB",])))
      
      locgtl <- twas_ds %>% filter(tissue=="GTL" & genestart>=xlow & genestop<=xhigh & chr==locchr & 
                                      phenoname==locpheno() & is.na(HLARegion) & is.na(MHCRegion)) %>%
        select(genename,chr, phenoname,log10pval,genestart,genestop,genemid, SignifGene)
      pthreshgtl <- -log10(0.05/(nrow(twas_ds[twas_ds$tissue=="GTL",])))
      
      locmsa <- twas_ds %>% filter(tissue=="MSA" & genestart>=xlow & genestop<=xhigh & chr==locchr & 
                                      phenoname==locpheno() & is.na(HLARegion) & is.na(MHCRegion)) %>%
        select(genename,chr, phenoname,log10pval,genestart,genestop,genemid, SignifGene)
      pthreshmsa <- -log10(0.05/(nrow(twas_ds[twas_ds$tissue=="MSA",])))
      
      
      #####################
      
      
      # merge secondary reference panel data sets with DGN
      dgngwb <- merge(dgntbl,locgwb, by=c("genename","chr","phenoname"), all = T)
      dgngwb$inboth <- ifelse(complete.cases(dgngwb$log10pval.x,dgngwb$log10pval.y),"In Both","Not in Both")
      
      dgngtl <- merge(dgntbl,locgtl, by=c("genename","chr","phenoname"), all = T)
      dgngtl$inboth <- ifelse(complete.cases(dgngtl$log10pval.x,dgngtl$log10pval.y),"In Both","Not in Both")
      
      dgnmsa <- merge(dgntbl,locmsa, by=c("genename","chr","phenoname"), all = T)
      dgnmsa$inboth <- ifelse(complete.cases(dgnmsa$log10pval.x,dgnmsa$log10pval.y),"In Both","Not in Both")
      
      
      #####################
      
      
      # set y limit values for each plot
      yhighdgn <- max(dgngwb$log10pval.x, na.rm = T)+0.15*max(dgngwb$log10pval.x, na.rm = T)
      nudgevaldgn <- 0.05*yhighdgn
      
      yhighgtl <- max(dgngtl$log10pval.y,pthreshgtl, na.rm = T)+0.15*max(dgngtl$log10pval.y,pthreshgtl, na.rm = T)
      nudgevalgtl <- 0.05*yhighgtl
      
      yhighgwb <- max(dgngwb$log10pval.y,pthreshgwb, na.rm = T)+0.15*max(dgngwb$log10pval.y,pthreshgwb, na.rm = T)
      nudgevalgwb <- 0.05*yhighgwb
      
      yhighmsa <- max(dgnmsa$log10pval.y,pthreshmsa, na.rm = T)+0.15*max(dgnmsa$log10pval.y,pthreshmsa, na.rm = T)
      nudgevalmsa <- 0.05*yhighmsa
      
      
      myColors <- setNames( c('#000000', '#56B4E9', '#E69F00'),
                            c("In Both","Not in Both",NA))
      
      colScale <- scale_colour_manual(values = myColors)
      
      
      ####################
      
      # plot info for GWB
      atop <- ggplotly(ggplot(data=dgngwb, aes(x=round(genemid.x/1000000,4),y=log10pval.x, color=inboth)) +
                         geom_point(pch=15) +
                         geom_segment(aes(x=round(genestart.x/1000000,4),y=log10pval.x, xend=round(genestop.x/1000000,4), yend=log10pval.x), size=2) +
                         geom_text(aes(label=genename), nudge_y = nudgevaldgn) +
                         geom_hline(aes(yintercept=pthresh), lty=2, color='red') +
                         xlim(round(xlow/1000000,4),round(xhigh/1000000,4)) +
                         ylim(0,yhighdgn) +
                         theme_bw() +
                         theme(legend.position = 'top', legend.title = element_blank()) +
                         annotate(geom="text",x=round(xlow/1000000,4)+.370000, y=pthresh+nudgevaldgn, color='red',
                                  label=paste0("DGN TWAS p-value: ", formatC(10^-(pthresh), format = "e", digits = 2)),size=4) +
                         
                         colScale
      )
      
      atop <- atop %>% plotly::layout(yaxis = list(title = ' DGN TWAS -log10(p)'),
                              xaxis = list(range=c(round(xlow/1000000,4),round(xhigh/1000000,4))),
                              annotations=list(x = 0.5 , y = 1.1, text = "(a) DGN vs. GWB", showarrow = F, 
                                               xref='paper', yref='paper',xanchor='center'),
                              legend = list(orientation='h', x=0, y=1),
                              title=list(text=paste0(loctitle(),"\n"),x=0,xanchor='left'), margin=list(t=60)
      )
      
      #####
      
      abottom <- ggplotly(ggplot(data=dgngwb, aes(x=round(genemid.y/1000000,4),y=-log10pval.y, color=inboth)) +
                            geom_hline(aes(yintercept=-pthreshgwb), lty=2, color='red') +
                            geom_point(pch=15) +
                            geom_segment(aes(x=round(genestart.y/1000000,4),y=-log10pval.y, xend=round(genestop.y/1000000,4), yend=-log10pval.y), size=2) +
                            geom_text(aes(label=genename), nudge_y = -1.5*nudgevalgwb) +
                            xlim(round(xlow/1000000,4),round(xhigh/1000000,4)) +
                            ylim(-yhighgwb,0) +
                            theme_bw() +
                            theme(legend.position = 'top', legend.title = element_blank()) +
                            annotate(geom="text",x=round(xlow/1000000,4)+.370000, y=-pthreshgwb-nudgevalgwb, color='red',
                                     label=paste0("GWB TWAS p-value: ", formatC(10^-(pthreshgwb), format = "e", digits = 2)),size=4) +
                            colScale
      )
      
      abottom <- abottom %>% plotly::layout(yaxis = list(title = 'GWB TWAS log10(p)', 
                                                 range=c(-yhighgwb,0.1)),
                                    xaxis = list(range=c(round(xlow/1000000,4),round(xhigh/1000000,4)), title="position (in Mb)"))
      
      #####
      
      afig <- subplot(atop,style(abottom,showlegend=F),nrows=2,shareX = T, titleX = T, titleY = T, which_layout = 1, margin=.01)
      
    })
    
    
    output$CompRefGTL <- renderPlotly({
      
      # select DGN results
      locds <- dgnds %>% filter(locvar==input$locuslst)
      xhigh <- max(locds$genestop)+1000000
      xlow <- max(0,min(locds$genestart)-1000000)
      locchr <- unique(locds$chr)
      #locpheno <- unique(locds$pheno)
      #locphcat <- unique(locds$phenocat)
      #loctitle <- paste0("chr ",locchr,": ",xlow,", ",xhigh,"; phenotype = ",locpheno,"; category = ",locphcat)
      
      
      #####################
      
      
      # select reference panel specific results
      dgntbl <- dgnds %>% filter(genestart>=xlow & genestop<=xhigh & phenoname==locpheno() & chr==locchr & 
                                   is.na(HLARegion) & is.na(MHCRegion)) %>%
        select(genename,chr, phenoname,genestart,genestop,genemid,log10pval)
      
      locgwb <- twas_ds %>% filter(tissue=="GWB" & genestart>=xlow & genestop<=xhigh & chr==locchr & 
                                      phenoname==locpheno() & is.na(HLARegion) & is.na(MHCRegion)) %>%
        select(genename,chr, phenoname,log10pval,genestart,genestop,genemid, SignifGene)
      pthreshgwb <- -log10(0.05/(nrow(twas_ds[twas_ds$tissue=="GWB",])))
      
      locgtl <- twas_ds %>% filter(tissue=="GTL" & genestart>=xlow & genestop<=xhigh & chr==locchr & 
                                      phenoname==locpheno() & is.na(HLARegion) & is.na(MHCRegion)) %>%
        select(genename,chr, phenoname,log10pval,genestart,genestop,genemid, SignifGene)
      pthreshgtl <- -log10(0.05/(nrow(twas_ds[twas_ds$tissue=="GTL",])))
      
      locmsa <- twas_ds %>% filter(tissue=="MSA" & genestart>=xlow & genestop<=xhigh & chr==locchr & 
                                      phenoname==locpheno() & is.na(HLARegion) & is.na(MHCRegion)) %>%
        select(genename,chr, phenoname,log10pval,genestart,genestop,genemid, SignifGene)
      pthreshmsa <- -log10(0.05/(nrow(twas_ds[twas_ds$tissue=="MSA",])))
      
      
      #####################
      
      
      # merge secondary reference panel data sets with DGN
      dgngwb <- merge(dgntbl,locgwb, by=c("genename","chr","phenoname"), all = T)
      dgngwb$inboth <- ifelse(complete.cases(dgngwb$log10pval.x,dgngwb$log10pval.y),"In Both","Not in Both")
      
      dgngtl <- merge(dgntbl,locgtl, by=c("genename","chr","phenoname"), all = T)
      dgngtl$inboth <- ifelse(complete.cases(dgngtl$log10pval.x,dgngtl$log10pval.y),"In Both","Not in Both")
      
      dgnmsa <- merge(dgntbl,locmsa, by=c("genename","chr","phenoname"), all = T)
      dgnmsa$inboth <- ifelse(complete.cases(dgnmsa$log10pval.x,dgnmsa$log10pval.y),"In Both","Not in Both")
      
      
      #####################
      
      
      # set y limit values for each plot
      yhighdgn <- max(dgngwb$log10pval.x, na.rm = T)+0.15*max(dgngwb$log10pval.x, na.rm = T)
      nudgevaldgn <- 0.05*yhighdgn
      
      yhighgtl <- max(dgngtl$log10pval.y,pthreshgtl, na.rm = T)+0.15*max(dgngtl$log10pval.y,pthreshgtl, na.rm = T)
      nudgevalgtl <- 0.05*yhighgtl
      
      yhighgwb <- max(dgngwb$log10pval.y,pthreshgwb, na.rm = T)+0.15*max(dgngwb$log10pval.y,pthreshgwb, na.rm = T)
      nudgevalgwb <- 0.05*yhighgwb
      
      yhighmsa <- max(dgnmsa$log10pval.y,pthreshmsa, na.rm = T)+0.15*max(dgnmsa$log10pval.y,pthreshmsa, na.rm = T)
      nudgevalmsa <- 0.05*yhighmsa
      
      
      myColors <- setNames( c('#000000', '#56B4E9', '#E69F00'),
                            c("In Both","Not in Both",NA))
      
      colScale <- scale_colour_manual(values = myColors)
      
      
      ####################
      
      # plot info for GTL
      btop <- ggplotly(ggplot(data=dgngtl, aes(x=round(genemid.x/1000000,4),y=log10pval.x, color=inboth)) +
                         geom_point(pch=15) +
                         geom_segment(aes(x=round(genestart.x/1000000,4),y=log10pval.x, xend=round(genestop.x/1000000,4), yend=log10pval.x), size=2) +
                         geom_text(aes(label=genename), nudge_y = nudgevaldgn) +
                         geom_hline(aes(yintercept=pthresh), lty=2, color='red') +
                         xlim(round(xlow/1000000,4),round(xhigh/1000000,4)) +
                         ylim(0,yhighdgn) +
                         theme_bw() +
                         theme(legend.position = 'top', legend.title = element_blank()) +
                         annotate(geom="text",x=round(xlow/1000000,4)+.370000, y=pthresh+nudgevaldgn, color='red',
                                  label=paste0("DGN TWAS p-value: ", formatC(10^-(pthresh), format = "e", digits = 2)),size=4) +
                         colScale
      )
      
      btop <- btop %>% plotly::layout(yaxis = list(title = ' DGN TWAS -log10(p)'),
                              xaxis = list(range=c(round(xlow/1000000,4),round(xhigh/1000000,4))),
                              annotations=list(x = 0.5 , y = 1.1, text = "(b) DGN vs. GTL", showarrow = F, 
                                               xref='paper', yref='paper',xanchor='center'),
                              legend = list(orientation='h', x=0, y=1),
                              title=list(text=paste0(loctitle(),"\n"),x=0,xanchor='left'), margin=list(t=60))
      
      #####
      
      bbottom <- ggplotly(ggplot(data=dgngtl, aes(x=round(genemid.y/1000000,4),y=-log10pval.y, color=inboth)) +
                            geom_hline(aes(yintercept=-pthreshgtl), lty=2, color='red') +
                            geom_point(pch=15) +
                            geom_segment(aes(x=round(genestart.y/1000000,4),y=-log10pval.y, xend=round(genestop.y/1000000,4), yend=-log10pval.y), size=2) +
                            geom_text(aes(label=genename), nudge_y = -1.5*nudgevalgtl) +
                            xlim(round(xlow/1000000,4),round(xhigh/1000000,4)) +
                            ylim(-yhighgtl,0) +
                            theme_bw() +
                            theme(legend.position = 'top', legend.title = element_blank()) +
                            annotate(geom="text",x=round(xlow/1000000,4)+.370000, y=-pthreshgtl-nudgevaldgn, color='red',
                                     label=paste0("GTL TWAS p-value: ", formatC(10^-(pthreshgtl), format = "e", digits = 2)),size=4) +
                            
                            colScale
      )
      
      bbottom <- bbottom %>% plotly::layout(yaxis = list(title = 'GTL TWAS log10(p)', 
                                                 range=c(-yhighgtl,0.1)),
                                    xaxis = list(range=c(round(xlow/1000000,4),round(xhigh/1000000,4)), title="position (in Mb)"))
      
      #####
      
      bfig <- subplot(btop,style(bbottom,showlegend=F),nrows=2,shareX = T, titleX = T, titleY = T, which_layout = 1, margin=0.01)
      
    })
    
    
    
    output$CompRefMSA <- renderPlotly({
      
      # select DGN results
      locds <- dgnds %>% filter(locvar==input$locuslst)
      xhigh <- max(locds$genestop)+1000000
      xlow <- max(0,min(locds$genestart)-1000000)
      locchr <- unique(locds$chr)
      #locpheno <- unique(locds$pheno)
      
      
      #####################
      
      
      # select reference panel specific results
      dgntbl <- dgnds %>% filter(genestart>=xlow & genestop<=xhigh & phenoname==locpheno() & chr==locchr & 
                                   is.na(HLARegion) & is.na(MHCRegion)) %>%
        select(genename,chr, phenoname,genestart,genestop,genemid,log10pval)
      
      locgwb <- twas_ds %>% filter(tissue=="GWB" & genestart>=xlow & genestop<=xhigh & chr==locchr & 
                                      phenoname==locpheno() & is.na(HLARegion) & is.na(MHCRegion)) %>%
        select(genename,chr, phenoname,log10pval,genestart,genestop,genemid, SignifGene)
      pthreshgwb <- -log10(0.05/(nrow(twas_ds[twas_ds$tissue=="GWB",])))
      
      locgtl <- twas_ds %>% filter(tissue=="GTL" & genestart>=xlow & genestop<=xhigh & chr==locchr & 
                                      phenoname==locpheno() & is.na(HLARegion) & is.na(MHCRegion)) %>%
        select(genename,chr, phenoname,log10pval,genestart,genestop,genemid, SignifGene)
      pthreshgtl <- -log10(0.05/(nrow(twas_ds[twas_ds$tissue=="GTL",])))
      
      locmsa <- twas_ds %>% filter(tissue=="MSA" & genestart>=xlow & genestop<=xhigh & chr==locchr & 
                                      phenoname==locpheno() & is.na(HLARegion) & is.na(MHCRegion)) %>%
        select(genename,chr, phenoname,log10pval,genestart,genestop,genemid, SignifGene)
      pthreshmsa <- -log10(0.05/(nrow(twas_ds[twas_ds$tissue=="MSA",])))
      
      
      #####################
      
      
      # merge secondary reference panel data sets with DGN
      dgngwb <- merge(dgntbl,locgwb, by=c("genename","chr","phenoname"), all = T)
      dgngwb$inboth <- ifelse(complete.cases(dgngwb$log10pval.x,dgngwb$log10pval.y),"In Both","Not in Both")
      
      dgngtl <- merge(dgntbl,locgtl, by=c("genename","chr","phenoname"), all = T)
      dgngtl$inboth <- ifelse(complete.cases(dgngtl$log10pval.x,dgngtl$log10pval.y),"In Both","Not in Both")
      
      dgnmsa <- merge(dgntbl,locmsa, by=c("genename","chr","phenoname"), all = T)
      dgnmsa$inboth <- ifelse(complete.cases(dgnmsa$log10pval.x,dgnmsa$log10pval.y),"In Both","Not in Both")
      
      
      #####################
      
      
      # set y limit values for each plot
      yhighdgn <- max(dgngwb$log10pval.x, na.rm = T)+0.15*max(dgngwb$log10pval.x, na.rm = T)
      nudgevaldgn <- 0.05*yhighdgn
      
      yhighgtl <- max(dgngtl$log10pval.y,pthreshgtl, na.rm = T)+0.15*max(dgngtl$log10pval.y,pthreshgtl, na.rm = T)
      nudgevalgtl <- 0.05*yhighgtl
      
      yhighgwb <- max(dgngwb$log10pval.y,pthreshgwb, na.rm = T)+0.15*max(dgngwb$log10pval.y,pthreshgwb, na.rm = T)
      nudgevalgwb <- 0.05*yhighgwb
      
      yhighmsa <- max(dgnmsa$log10pval.y,pthreshmsa, na.rm = T)+0.15*max(dgnmsa$log10pval.y,pthreshmsa, na.rm = T)
      nudgevalmsa <- 0.05*yhighmsa
      
      
      myColors <- setNames( c('#000000', '#56B4E9', '#E69F00'),
                            c("In Both","Not in Both",NA))
      
      colScale <- scale_colour_manual(values = myColors)
      
      
      ####################
      
      # plot info for MSA
      ctop <- ggplotly(ggplot(data=dgnmsa, aes(x=round(genemid.x/1000000,4),y=log10pval.x, color=inboth)) +
                         geom_point(pch=15) +
                         geom_segment(aes(x=round(genestart.x/1000000,4),y=log10pval.x, xend=round(genestop.x/1000000,4), yend=log10pval.x), size=2) +
                         geom_text(aes(label=genename), nudge_y = nudgevaldgn) +
                         geom_hline(aes(yintercept=pthresh), lty=2, color='red') +
                         xlim(round(xlow/1000000,4),round(xhigh/1000000,4)) +
                         ylim(0,yhighdgn) +
                         theme_bw() +
                         theme(legend.position = 'top', legend.title = element_blank()) +
                         annotate(geom="text",x=round(xlow/1000000,4)+.370000, y=pthresh+nudgevaldgn, color='red',
                                  label=paste0("DGN TWAS p-value: ", formatC(10^-(pthresh), format = "e", digits = 2)),size=4) +
                         colScale
      )
      
      ctop <- ctop %>% plotly::layout(yaxis = list(title = 'DGN TWAS -log10(p)'),
                              xaxis = list(range=c(round(xlow/1000000,4),round(xhigh/1000000,4))),
                              annotations=list(x = 0.5 , y = 1.1, text = "(c) DGN vs. MSA", showarrow = F, 
                                               xref='paper', yref='paper',xanchor='center'),
                              legend = list(orientation='h', x=0, y=1),
                              title=list(text=paste0(loctitle(),"\n"),x=0,xanchor='left'),margin=list(t=60))
      
      #####
      
      cbottom <- ggplotly(ggplot(data=dgnmsa, aes(x=round(genemid.y/1000000,4),y=-log10pval.y, color=inboth)) +
                            geom_hline(aes(yintercept=-pthreshmsa), lty=2, color='red') +
                            geom_point(pch=15) +
                            geom_segment(aes(x=round(genestart.y/1000000,4),y=-log10pval.y, xend=round(genestop.y/1000000,4), yend=-log10pval.y), size=2) +
                            geom_text(aes(label=genename), nudge_y = -1.5*nudgevalmsa) +
                            xlim(round(xlow/1000000,4),round(xhigh/1000000,4)) +
                            theme_bw() +
                            ylim(-yhighmsa-2,0) +
                            annotate(geom="text",x=round(xlow/1000000,4)+.370000, y=-pthreshmsa-nudgevaldgn, color='red',
                                     label=paste0("MSA TWAS p-value: ", formatC(10^-(pthreshmsa), format = "e", digits = 2)),size=4) +
                            
                            colScale
      )
      
      cbottom <- cbottom %>% plotly::layout(yaxis = list(title = 'MSA TWAS log10(p)', 
                                                 range=c(-yhighmsa,0.1)),
                                    xaxis = list(range=c(round(xlow/1000000,4),round(xhigh/1000000,4)), title="position (in Mb)"))
      
      #####
      
      cfig <- subplot(ctop,style(cbottom,showlegend=F),nrows=2,shareX = T, titleX = T, titleY = T, which_layout = 1,
                      margin = 0.01)
      
      
      #####################
      
    })
    
    ############################################################
    
    
    # # TWAS results table
    # output$TWAStbl <- DT::renderDataTable({
    #   locds <- dgnds %>% filter(locus==input$locuslst)
    #   xhigh <- max(locds$genestop)+1000000
    #   xlow <- max(0,min(locds$genestart)-1000000)
    #   locpheno<-unique(locds$pheno) #locus phenotype
    #   locchr<- unique(locds$chr) #locus chromosome
    #   
    #   # select all genes at locus
    #   dgntbl <- dgnds %>% filter(genestart>=xlow & genestop <=xhigh & pheno==locpheno & chr==locchr) %>%
    #     select(genename,chr, pheno,phenocat,se_beta_, p, SignifGene,HLARegion,MHCRegion)
    #   colnames(dgntbl) <- c("gene","chr","trait","trait category","beta","p-value","TWAS significant",
    #                         "HLA gene","MHC region")
    #   
    #   # print data table
    #   datatable(dgntbl, options = list(
    #     order = list(list(6, 'asc')))) %>% formatStyle(
    #     "TWAS significant", target = "row",
    #     backgroundColor = styleEqual(c(1), c('lightgreen'))
    #     ) %>% formatStyle(
    #       "HLA gene", target = "row",
    #       backgroundColor = styleEqual(c(1), c('red'))
    #     ) %>% formatStyle(
    #       "MHC region", target = "row",
    #       backgroundColor = styleEqual(c(1), c('red'))
    #     )
    # })
    
    
    ############################################################
    
    
    # TWAS results table with separate tabs for each reference panel
    output$DGNtbl <- DT::renderDataTable({
      
      locds <- dgnds %>% filter(locvar==input$locuslst)
      xhigh <- max(locds$genestop)+1000000
      xlow <- max(0,min(locds$genestart)-1000000)
      locpheno <- unique(locds$phenoname) #locus phenotype
      locchr<- unique(locds$chr) #locus chromosome
      
      # select all genes at locus
      dgntbl <- dgnds %>% filter(genestart>=xlow & genestop <=xhigh & phenoname==locpheno & chr==locchr) %>%
        select(genename,chr, genestart, genestop, phenoname,phenocat,se_beta_, p, 
               weightsnp, cohortsnp, SignifGene,HLARegion,MHCRegion,SingleSNP)
      colnames(dgntbl) <- c("gene","chr","gene start","gene stop","trait","trait category","beta","p-value",
                            "# model SNPs","# SNPs used","TWAS significant",
                            "HLA gene","MHC region","1 SNP model")
      
      # print data table
      datatable(dgntbl, rownames=F, options = list(
        order = list(list(7, 'asc')))) %>% formatStyle(
          "TWAS significant", target = "row",
          backgroundColor = styleEqual(c(1), c('lightgreen'))
        ) %>% formatStyle(
          "HLA gene", target = "row",
          backgroundColor = styleEqual(c(1), c('red'))
        ) %>% formatStyle(
          "MHC region", target = "row",
          backgroundColor = styleEqual(c(1), c('red'))
        ) %>% formatStyle(
          "1 SNP model", target = "row",
          backgroundColor = styleEqual(c(1), c('red'))
        )
    })
    
    
    #####################
    
    
    output$GTLtbl <- DT::renderDataTable({
      
      locds <- dgnds %>% filter(locvar==input$locuslst)
      xhigh <- max(locds$genestop)+1000000
      xlow <- max(0,min(locds$genestart)-1000000)
      locpheno <- unique(locds$phenoname) #locus phenotype
      locchr<- unique(locds$chr) #locus chromosome
      
      # select all genes at locus
      gtltbl <- twas_ds %>% filter(tissue=='GTL',genestart>=xlow & genestop <=xhigh & phenoname==locpheno & chr==locchr) %>%
        select(genename,chr, genestart, genestop, phenoname,phenocat,se_beta_, p, 
               weightsnp,cohortsnp, SignifGene,HLARegion,MHCRegion, SingleSNP)
      colnames(gtltbl) <- c("gene","chr","gene start","gene stop","trait","trait category","beta","p-value",
                            "# model SNPs","#SNPs used","TWAS significant",
                            "HLA gene","MHC region","1 SNP model")
      
      # print data table
      datatable(gtltbl, rownames=F, options = list(
        order = list(list(7, 'asc')))) %>% formatStyle(
          "TWAS significant", target = "row",
          backgroundColor = styleEqual(c(1), c('lightgreen'))
        ) %>% formatStyle(
          "HLA gene", target = "row",
          backgroundColor = styleEqual(c(1), c('red'))
        ) %>% formatStyle(
          "MHC region", target = "row",
          backgroundColor = styleEqual(c(1), c('red'))
        ) %>% formatStyle(
          "1 SNP model", target = "row",
          backgroundColor = styleEqual(c(1), c('red'))
        )
    })
    
    
    #####################
    
    
    output$GWBtbl <- DT::renderDataTable({
      
      locds <- dgnds %>% filter(locvar==input$locuslst)
      xhigh <- max(locds$genestop)+1000000
      xlow <- max(0,min(locds$genestart)-1000000)
      locpheno <- unique(locds$phenoname) #locus phenotype
      locchr<- unique(locds$chr) #locus chromosome
      
      # select all genes at locus
      gwbtbl <- twas_ds %>% filter(tissue=='GWB',genestart>=xlow & genestop <=xhigh & phenoname==locpheno & chr==locchr) %>%
        select(genename,chr, genestart,genestop,phenoname,phenocat,se_beta_, p,weightsnp,cohortsnp, 
               SignifGene,HLARegion,MHCRegion,SingleSNP)
      colnames(gwbtbl) <- c("gene","chr","gene start","gene stop","trait","trait category","beta","p-value",
                            "# model SNPs","# SNPs used","TWAS significant",
                            "HLA gene","MHC region","1 SNP model")
      
      # print data table
      datatable(gwbtbl, rownames=F, options = list(
        order = list(list(7, 'asc')))) %>% formatStyle(
          "TWAS significant", target = "row",
          backgroundColor = styleEqual(c(1), c('lightgreen'))
        ) %>% formatStyle(
          "HLA gene", target = "row",
          backgroundColor = styleEqual(c(1), c('red'))
        ) %>% formatStyle(
          "MHC region", target = "row",
          backgroundColor = styleEqual(c(1), c('red'))
        ) %>% formatStyle(
          "1 SNP model", target = "row",
          backgroundColor = styleEqual(c(1), c('red'))
        )
    })
    
    
    #####################
    
    
    output$MSAtbl <- DT::renderDataTable({
      
      locds <- dgnds %>% filter(locvar==input$locuslst)
      xhigh <- max(locds$genestop)+1000000
      xlow <- max(0,min(locds$genestart)-1000000)
      locpheno <- unique(locds$phenoname) #locus phenotype
      locchr<- unique(locds$chr) #locus chromosome
      
      # select all genes at locus
      msatbl <- twas_ds %>% filter(tissue=='MSA',genestart>=xlow & genestop <=xhigh & phenoname==locpheno & chr==locchr) %>%
        select(genename,chr, genestart,genestop,phenoname,phenocat,se_beta_, p,weightsnp,cohortsnp, 
               SignifGene,HLARegion,MHCRegion,SingleSNP)
      colnames(msatbl) <- c("gene","chr","gene start","gene stop","trait","trait category","beta","p-value",
                            "# model SNPs","#SNPs used","TWAS significant",
                            "HLA gene","MHC region","1 SNP model")
      
      # print data table
      datatable(msatbl, rownames=F, options = list(
        order = list(list(7, 'asc')))) %>% formatStyle(
          "TWAS significant", target = "row",
          backgroundColor = styleEqual(c(1), c('lightgreen'))
        ) %>% formatStyle(
          "HLA gene", target = "row",
          backgroundColor = styleEqual(c(1), c('red'))
        ) %>% formatStyle(
          "MHC region", target = "row",
          backgroundColor = styleEqual(c(1), c('red'))
        )  %>% formatStyle(
          "1 SNP model", target = "row",
          backgroundColor = styleEqual(c(1), c('red'))
        )
    })
    
    
    ############################################################
    
    
    # TWAS other traits results table with separate tabs for each reference panel
    output$DGNtbltrt <- DT::renderDataTable({
      
      # select dgn all genes from other trait categories at locus
      dgntbl <- filter(dgnds, genestart>=xlow() & genestop<=xhigh() & phenoname!=locpheno() 
                       & chr==locchr() & phenocat==locphcat()) %>%
        select(genename,chr, genestart,genestop,phenoname,phenocat,se_beta_, p,weightsnp,cohortsnp, 
               SignifGene,HLARegion,MHCRegion,SingleSNP)
      colnames(dgntbl) <- c("gene","chr","gene start","gene stop","trait","trait category","beta","p-value",
                            "# model SNPs","#SNPs used","TWAS significant",
                            "HLA gene","MHC region","1 SNP model")
      
      # print data table
      datatable(dgntbl, rownames=F, options = list(
        order = list(list(7, 'asc')))) %>% formatStyle(
          "TWAS significant", target = "row",
          backgroundColor = styleEqual(c(1), c('lightgreen'))
        ) %>% formatStyle(
          "HLA gene", target = "row",
          backgroundColor = styleEqual(c(1), c('red'))
        ) %>% formatStyle(
          "MHC region", target = "row",
          backgroundColor = styleEqual(c(1), c('red'))
        ) %>% formatStyle(
          "1 SNP model", target = "row",
          backgroundColor = styleEqual(c(1), c('red'))
        )
    })
    
    
    #####################
    
    
    output$GTLtbltrt <- DT::renderDataTable({
      
      #locds <- dgnds %>% filter(locvar==input$locuslst)
      #xhigh <- max(locds$genestop)+1000000
      #xlow <- max(0,min(locds$genestart)-1000000)
      #locpheno <- unique(locds$pheno) #locus phenotype
      #locchr<- unique(locds$chr) #locus chromosome
      
      # select all genes at locus
      gtltbl <- twas_ds %>% filter(tissue=='GTL',genestart>=xlow() & genestop<=xhigh() 
                                    & phenoname!=locpheno() & chr==locchr() & phenocat==locphcat()) %>%
        select(genename,chr, genestart,genestop,phenoname,phenocat,se_beta_, p,weightsnp,cohortsnp, 
               SignifGene,HLARegion,MHCRegion,SingleSNP)
      colnames(gtltbl) <- c("gene","chr","gene start","gene stop","trait","trait category","beta","p-value",
                            "# model SNPs","#SNPs used","TWAS significant",
                            "HLA gene","MHC region","1 SNP model")
      
      # print data table
      datatable(gtltbl, rownames=F, options = list(
        order = list(list(7, 'asc')))) %>% formatStyle(
          "TWAS significant", target = "row",
          backgroundColor = styleEqual(c(1), c('lightgreen'))
        ) %>% formatStyle(
          "HLA gene", target = "row",
          backgroundColor = styleEqual(c(1), c('red'))
        ) %>% formatStyle(
          "MHC region", target = "row",
          backgroundColor = styleEqual(c(1), c('red'))
        ) %>% formatStyle(
          "1 SNP model", target = "row",
          backgroundColor = styleEqual(c(1), c('red'))
        )
    })
    
    
    #####################
    
    
    output$GWBtbltrt <- DT::renderDataTable({
      
      # select all genes at locus
      gwbtbl <- twas_ds %>% filter(tissue=='GWB',genestart>=xlow() & genestop<=xhigh()
                                    & phenoname!=locpheno() & chr==locchr() & phenocat==locphcat()) %>%
        select(genename,chr, genestart,genestop,phenoname,phenocat,se_beta_, p,weightsnp,cohortsnp, 
               SignifGene,HLARegion,MHCRegion,SingleSNP)
      colnames(gwbtbl) <- c("gene","chr","gene start","gene stop","trait","trait category","beta","p-value",
                            "# model SNPs","#SNPs used","TWAS significant",
                            "HLA gene","MHC region","1 SNP model")
      
      # print data table
      datatable(gwbtbl, rownames=F, options = list(
        order = list(list(7, 'asc')))) %>% formatStyle(
          "TWAS significant", target = "row",
          backgroundColor = styleEqual(c(1), c('lightgreen'))
        ) %>% formatStyle(
          "HLA gene", target = "row",
          backgroundColor = styleEqual(c(1), c('red'))
        ) %>% formatStyle(
          "MHC region", target = "row",
          backgroundColor = styleEqual(c(1), c('red'))
        ) %>% formatStyle(
          "1 SNP model", target = "row",
          backgroundColor = styleEqual(c(1), c('red'))
        )
    })
    
    
    #####################
    
    
    output$MSAtbltrt <- DT::renderDataTable({
      
      #locds <- dgnds %>% filter(locvar==input$locuslst)
      #xhigh <- max(locds$genestop)+1000000
      #xlow <- max(0,min(locds$genestart)-1000000)
      #locpheno <- unique(locds$pheno) #locus phenotype
      #locchr<- unique(locds$chr) #locus chromosome
      
      # select all genes at locus
      msatbl <- twas_ds %>% filter(tissue=='MSA',genestart>=xlow() & genestop<=xhigh() 
                                    & phenoname!=locpheno() & chr==locchr() & phenocat==locphcat()) %>%
        select(genename,chr, genestart,genestop,phenoname,phenocat,se_beta_, p,weightsnp,cohortsnp, 
               SignifGene,HLARegion,MHCRegion,SingleSNP)
      colnames(msatbl) <- c("gene","chr","gene start","gene stop","trait","trait category","beta","p-value",
                            "# model SNPs","#SNPs used","TWAS significant",
                            "HLA gene","MHC region","1 SNP model")
      
      # print data table
      datatable(msatbl, rownames=F, options = list(
        order = list(list(7, 'asc')))) %>% formatStyle(
          "TWAS significant", target = "row",
          backgroundColor = styleEqual(c(1), c('lightgreen'))
        ) %>% formatStyle(
          "HLA gene", target = "row",
          backgroundColor = styleEqual(c(1), c('red'))
        ) %>% formatStyle(
          "MHC region", target = "row",
          backgroundColor = styleEqual(c(1), c('red'))
        ) %>% formatStyle(
          "1 SNP model", target = "row",
          backgroundColor = styleEqual(c(1), c('red'))
        )
    })
    
    
    ############################################################
    
    # Table of known variants
    output$KnownSNPtbl <- DT::renderDataTable({
      locds <- dgnds %>% filter(locvar==input$locuslst)
      xhigh <- max(locds$genestop)+1000000
      xlow <- max(0,min(locds$genestart)-1000000)
      locpheno<-unique(locds$phenocat) #locus phenotype category
      locph<-unique(locds$phenoname) #locus phenotype
      locchr<- unique(locds$chr) #locus chromosome
      
      if (locpheno=="RBC") {
        ds=rbcds
      } else if (locpheno=="WBC") {
        ds=wbcds
      } else {
        ds=pltds
      }
      
      # select all known snps at locus
      phenotbl <- ds %>% filter(X6>=xlow & xhigh>=X6 & X5==locchr) %>%
        select(X1,X2,X3,X4,X5,X6,X8,X9,X10,X11,X12,X13,X14,X15,X21)
      colnames(phenotbl) <- c("RSID", "Reference", "Ancestry", "Trait", "Chr", "Pos_b37", "Gene", "Effect Allele",
                              "Other Allele", "EAF", "Beta", "SE", "P", "N","In Kaiser")
      
      if (nrow(phenotbl)==0){
        #df <- data.frame(a="No reported GWAS variants at this locus")
        datatable(phenotbl,
                  options=list(columnDefs = list(list(visible=FALSE, targets=c(15,16,17,18,19,20,21,22,23)))))
      } else {
        
        # select set of TWAS insignificant genes
        dgntbl <- dgnds %>% filter(genestart>=xlow & genestop<=xhigh & phenoname==locph & chr==locchr &
                                     SignifGene==0  & is.na(HLARegion) & is.na(MHCRegion)) %>% select(genename)
        
        # indicator for significant TWAS gene pattern 
        phenotbl$X16 <- ifelse(grepl(paste(unique(locds$genename), collapse="|"),phenotbl$Gene),1,0)
        
        # indicator for non-significant TWAS gene pattern
        phenotbl$X17 <- ifelse(grepl(paste(unique(dgntbl$genename),collapse="|"),phenotbl$Gene),1,0)
        
        # indicator for gene not predicted in Kaiser
        phenotbl$notpred <- ifelse(grepl(paste(genes$genename[genes$X0==0],collapse = "|"),phenotbl$Gene),1,0)
        phenotbl$notdgn <- ifelse(grepl(paste(genes$genename[genes$X0==1],collapse="|"),phenotbl$Gene),1,0)
        
        phenotbl$chrposall <- paste0(phenotbl$Chr,":",phenotbl$Pos_b37,":",phenotbl$`Effect Allele`,":",phenotbl$`Other Allele`)
        phenotbl$chrposall2 <- paste0(phenotbl$Chr,":",phenotbl$Pos_b37,":",phenotbl$`Other Allele`,":",phenotbl$`Effect Allele`)
        
        # indicator for snps included in GWAS results
        locgwas <- gwasdsfin %>% filter(Locus==locnum()) %>% select(SNP)
        phenotbl$X19 <- ifelse(grepl(paste(unique(locgwas$SNP),collapse="|"),phenotbl$chrposall),1,0)
        phenotbl$X20 <- ifelse(grepl(paste(unique(locgwas$SNP),collapse="|"),phenotbl$chrposall2),1,0)
        
        # print data table
        datatable(phenotbl, #rownames=FALSE,
                  options=list(columnDefs = list(list(visible=FALSE, 
                                                      targets=c(15,16,17,18,19,20,21,22,23))))
        ) %>%
          formatStyle("Gene", "X16", backgroundColor = styleEqual(
            c(1), c("lightgreen"))
          ) %>% 
          formatStyle("In Kaiser", target = "row",
                      backgroundColor = styleEqual(c(0), c('red'))
          ) %>%
          formatStyle("Gene", "X17", backgroundColor = styleEqual(
            c(1), c("yellow"))
          ) %>%
          formatStyle("RSID", "X19", backgroundColor = styleEqual(
            c(1), c("lightgreen"))
          ) %>%
          formatStyle("RSID", "X20", backgroundColor = styleEqual(
            c(1), c("lightgreen"))
          ) %>%
          formatStyle("Gene", "notpred", backgroundColor = styleEqual(
            c(1), c("orange"))
          ) %>%
          formatStyle("Gene", "notdgn", backgroundColor = styleEqual(
            c(1), c("red"))
          )
      }
    })
    
    
    ############################################################
    
    # Table of known variants
    output$GWASvars <- DT::renderDataTable({
      locgwas <- gwasdsfin %>% filter(Locus==locnum())
      
      # print data table
      datatable(locgwas, #rownames=F, 
                options=list(columnDefs = list(list(visible=FALSE, targets=c(2,3,4,5,9,10,11)))))
    })
    
  }
  


############################################################################################


  
  # render the app
  shiny::shinyApp(ui=ui, server=server)

}
