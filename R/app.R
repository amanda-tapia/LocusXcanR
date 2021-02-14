#' LocusXcanR: An R Shiny App to Visualize TWAS Results Integrated with Other 'Omics Studies
#' 
#' \code{LocusXcanR} allows users to import TWAS and GWAS results (and additional important details) 
#' in order to visualize them together in a single web page
#' application via R Shiny and facilitates TWAS results interpretation in context.
#' 
#' Please reference the vignette here () for further details on required and optional variables for each input dataset
#' 
#' Created by: Amanda L. Tapia, Date: June 8, 2020, Updated: December 16, 2020
#'
#' @param twas_result character, file path to TWAS results (required). See Details for required column names.
#' @param pvalthresh numeric, -log10 p-value threshold for TWAS results (required).
#' @param weight_tbl character, file path to TWAS weights database (required). See Details for required column names.
#' @param known_variants character, file path to known GWAS variants (required). See Details for required column names.
#' @param pred_exp_corr character, file path to predicted expression correlation (required). See Details for required column names.
#' @param study_name character, the name of the study (optional, default is missing)
#' @param conditional_present logical, TRUE if conditional analysis results are available for plotting, FALSE otherwise (default is FALSE)
#' @param multiple_tissues logical, TRUE if TWAS results are available from more than one tissue, FALSE if TWAS results are only available from a single tissue or a multi-tissue analysis (default is FALSE)
#' @param known_gwas character, file path to study GWAS data matching known variants (required)
#' @param db_genes character, file path to a list of genes in the database (required)
#' @param all_gwas character, file path to study GWAS data (required)
#' @param ld_gwas character, file path to the LD among study variants or an LD reference panel (required)
#' @param ref_expr_name character, name of the reference expression data set used in the analysis (optional, default is missing)
#' @param head_details character, any additional header details to be included in the app (optional, default is no details). HTML formatting commands may be used.
#' @param method_details character, detailed methods section (optional, default is no details). HTLM formatting commands may be used.
#' @param primary_tissue character, if multiple tissues are available for analysis, list the name of the primary tissue (required if multiple_tissues = TRUE)
#' @param meta_present logical, TRUE if results from TWAS meta-analysis are present for comparison, FALSE otherwise (optional, default is FALSE)
#' @param meta_thresh numeric, p-value threshold for meta-analysis results (required if meta_present=TRUE)
#' @param ideogram_present logical, TRUE if an ideogramTrack should be plotted, FALSE otherwise (optional, default is FALSE)
#' @param genome_build character, the genome build for the ideogramTrack data (required if ideogram_present = TRUE)
#' @param cytoband_ds character, file path to the cytogenic bands data set (required if ideogram_present = TRUE)
#' @param add_UI character, additional elements to add to UI (optional, default is missing)
#' @param add_server character, additional elements to add to server (optional, default is missing)
#' 
#' @return An R Shiny application
#' 
#' @examples
#' 
#' # This is an example with only the required parameters
#' 
#' 
#' # This is an example with all parameters
#' 
#' 
#' @importFrom magrittr "%>%"
#' @importFrom tidyr "separate"
#' @importFrom dplyr "filter","select","distinct"
#' @importFrom plotly "renderPlotly","plotlyOutput","ggplotly","subplot","style"
#' @importFrom visNetwork "renderVisNetwork","visNetworkOutput","visNetwork","visOptions"
#' @importFrom data.table ":=","data.table","as.data.table"
#' @importFrom shiny "fluidPage","h3","h4","plotOutput","HTML","tabPanel","tabsetPanel","br","hr","strong","navbarPage","h5","fixedPanel","p"
#' @importFrom DT "formatStyle","styleEqual","datatable"
#' @importFrom ggplot2 "scale_colour_manual","ggplot","aes","theme","element_blank","geom_point","geom_hline","xlim","ylim","theme_bw","geom_segment","annotate","geom_text"
#' @importFrom utils "read.table"
#' @importFrom stats "complete.cases","setNames"
#' @importFrom rlang ".data"
#'
#
#
####################################################################


###
############## TO DO ###############################################
####


###
########## Begin Function #####################################################
###


#' @export

LocusXcanR <- function(twas_result,pvalthresh,weight_tbl,study_name="",pred_exp_corr,conditional_present=FALSE,multiple_tissues=FALSE,
                      known_variants,known_gwas,db_genes,all_gwas,ld_gwas,ref_expr_name="",head_details="",method_details="",
                      primary_tissue,meta_present=FALSE,meta_thresh,ideogram_present=FALSE,genome_build,
                      cytoband_ds,add_UI="",add_Server=""){
  
  # load analysis dataset
  twas_ds <- data.table::fread(twas_result,stringsAsFactors = F, header=T)
  
  # select variables for table displays
  #origvars <- names(twas_ds)
  #origvarstbl <- origvars[!origvars %in% c("locus2","locstart","locstop")]
  
  # define additional variables for plotting
  twas_ds$genestartMB <- round(twas_ds$genestart/1000000,4)
  twas_ds$genestopMB <- round(twas_ds$genestop/1000000,4)
  twas_ds$genemid <- (twas_ds$genestart+twas_ds$genestop)/2
  twas_ds$genemidMB <- round(twas_ds$genemid/1000000,4)
  twas_ds$log10pval <- (log10(twas_ds$p))*(-1)
  
  # calculate conditional p-value when conditional analysis results present
  # if (conditional_present==TRUE){
  #   if(is.na(twas_ds$p_conditional)==TRUE){
  #     twas_ds$p_final=twas_ds$p
  #   } else twas_ds$p_final=twas$p_conditional
  # }
  if (conditional_present==TRUE){
    twas_ds$p_final[is.na(twas_ds$p_conditional)] <- twas_ds$p[is.na(twas_ds$p_conditional)]
    twas_ds$p_final[!is.na(twas_ds$p_conditional)] <- twas_ds$p_conditional[!is.na(twas_ds$p_conditional)]
  } else{
    twas_ds$p_final <- twas_ds$p
  }
  
  
  # load known variants dataset
  GWAS_sentinel <- read.table(known_variants, stringsAsFactors = F, header=T, sep = '\t')
  
  # get GWAS_sentinel column names
  # GWAS_sentinel_col <- colnames(GWAS_sentinel)
  
  
  
  # study gwas data (known variants at the locus)
  cohort_gwas_known <- read.table(known_gwas,
                        stringsAsFactors = F, header = T, sep = '\t')
  cohort_gwas_known$log10p <- -log10(cohort_gwas_known$PVAL)
  cohort_gwas_known$poslog10p <- log10(cohort_gwas_known$PVAL)
  cohort_gwas_knownfin <- cohort_gwas_known %>% separate(.data$SNP,c("chr","pos","all1","all2"),sep=':',remove=F)

   
  # genes included in reference panel
  ref_panel_genes <- read.table(db_genes, stringsAsFactors = T, header = T)
  
  
  # load all Kaiser GWAS data (all GWAS with pvalues < 0.05)
  gwasall <- data.table::fread(all_gwas, header = T,stringsAsFactors = F)
  gwasallfin <- gwasall %>% separate(.data$SNP,c('chr','pos','all1','all2'),sep=':',remove=F,extra = 'merge')
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


  # data for ideogram
  if (ideogram_present == TRUE){
    data <- read.table(cytoband_ds, header = F, sep = "\t")
    colnames(data) <- c("chrom", "chromStart", "chromEnd", "name", "gieStain")
  }

   
 ###
 ########## Filter primary tissue ##########################################################
 ###

  
  # filter analysis dataset by primary tissue if multiple tissues present
  if (multiple_tissues==TRUE){
    
    # filter the analysis dataset
    primary_ref_ds <- twas_ds %>% filter(.data$tissue==primary_tissue)
    
    # set the list of tissues
    tissue_list <- unique(twas_ds$tissue[twas_ds$tissue!=primary_tissue])

  } else{
    primary_ref_ds <- twas_ds
    primary_tissue <- unique((twas_ds$tissue))
  }
  

  # set the locus specific variable for drop down list
  primary_ref_ds$locvar <- paste0("Locus ",primary_ref_ds$locus," : ", "chr ",primary_ref_ds$chr," : ",primary_ref_ds$locstart,", ",
                         primary_ref_ds$locstop," : ",primary_ref_ds$index," : ",primary_ref_ds$pheno
                        )

  # calculate -log10(meta-analysis p-value) if meta-analysis results are present
  if (meta_present==TRUE){
    primary_ref_ds$log10pvalmeta <- -log10(primary_ref_ds$p_meta)
  }
  
  #select(-locus2,-locstart,-locstop,-genestartMB,-genestopMB,-genemid,-genemidMB,-locvar,-log10pvalmeta,-log10pval)
  
  
  sortedprimary_ref_ <- primary_ref_ds[
    order( primary_ref_ds$locus, primary_ref_ds$genestart, primary_ref_ds$genestop, primary_ref_ds$pheno ),
  ]
  
  # list of unique loci
  loclist <- sort(unique(primary_ref_ds$locus[complete.cases(primary_ref_ds$locus)]))
  loclistdet <- unique(sortedprimary_ref_$locvar[complete.cases(sortedprimary_ref_$locus)])
  
  
  # p-value threshold for primary reference panel
  #pthresh <- -log10(0.05/(nrow(primary_ref_ds)))
  pthresh <- pvalthresh

  # number of significant TWAS genes
  # signifrow <- primary_ref_ds %>% filter(SignifGene==1 & is.na(HLARegion) & is.na(MHCRegion)
  #                              & SingleSNP!=1)
  # numsignif <- nrow(signifrow)
  # 
  # # meta-analysis threshold
  # metathresh <- log10(0.05/numsignif)
  
  
  # meta-analysis strict p-value threshold
  if (meta_present==TRUE){
    metathresh=log10(meta_thresh)
  }
  
  


####
####### Prep for UI #####################################################################################
####
  
  
  ############ Prep conditional UI ####################
  
  # Set the UI when conditional analysis results are present
  # radio button to toggle between marginal and conditional results, when conditional results available
  if (conditional_present==TRUE){
    condtext <- "Note: Top figure displays TWAS significant genes and any additional non-significant genes reported from GWAS, bottom figure displays GWAS variants. In the TWAS plot, \"reported in GWAS\" means that the study TWAS gene was reported in the GWAS catalog as the assigned gene for a single variant signal associated with the phenotype category, often based on physical proximity. In the GWAS plot, \"reported in GWAS\" means that the study GWAS variant was reported in the GWAS catalog as a single variant signal associated with the phenotype category. Marginal TWAS displays results of gene-trait associations. Conditional TWAS displays results of gene-trait associations, conditional on reported GWAS variants at the locus (conditional results only available for significant TWAS genes)."

    radios_cond <- shiny::radioButtons("margcondradio", label = h5(strong("Select TWAS results to display:")),
                                       choices = list("Marginal TWAS" = 1, "Conditional TWAS" = 2), 
                                       selected = 1, inline=T)
  } else if(conditional_present==FALSE) {
    condtext <- "Note: Top figure displays TWAS significant genes and any additional non-significant genes reported from GWAS, bottom figure displays GWAS variants. In the TWAS plot, \"reported in GWAS\" means that the study TWAS gene was reported in the GWAS catalog as the assigned gene for a single variant signal associated with the phenotype category, often based on physical proximity. In the GWAS plot, \"reported in GWAS\" means that the study GWAS variant was reported in the GWAS catalog as a single variant signal associated with the phenotype category."
    radios_cond <- ""
  }
  
  
  ############ Prep meta UI ####################
  
  # Set the UI when meta-analysis results are present for comparison
  # plot the figure if meta-analysis results are present, otherwise plot nothing
  if (meta_present==TRUE){
    meta_result1 <- shiny::h4(shiny::strong("Mirror plot of study TWAS and meta-analysis TWAS"))
    meta_result2 <- "Note: Top figure displays TWAS significant genes and any additional non-significant genes reported from GWAS, bottom figure displays the same genes from TWAS meta-analysis of meta-analysis cohort(s). In both plots, \"reported in GWAS\" means that the TWAS gene was reported in the GWAS catalog as the assigned gene for a single variant signal associated with the phenotype category, often based on physical proximity."
    meta_result3 <- plotly::plotlyOutput("Tmetamirror", height=600)
    meta_result4 <- shiny::br()
    meta_result5 <- shiny::hr()
  } else if(meta_present==FALSE) {
    meta_result1 <- ""
    meta_result2 <- ""
    meta_result3 <- ""
    meta_result4 <- ""
    meta_result5 <- ""
  }
  
  
  ############ Prep multi-tissue UI ####################
  
  # Set the UI when results from multiple tissues are available for comparison
  
 
  # plot a multi-tabbed figure if results are present, otherwise plot nothing
  if (multiple_tissues==TRUE){
    
    # multi-tissue plot set tabpanels functions as character strings
    numtabs <- length(tissue_list)
    a <- rep("tabPanel(paste0(primary_tissue,' vs. ",numtabs)
    b <- rep("'),plotlyOutput('CompRef",numtabs)
    c <- rep("', height=600)),",numtabs)
    alltab <- paste0(a,tissue_list,b,tissue_list,c)
    
    # multi-tissue plot cut comma from last entry
    alltab[numtabs] <- substring(alltab[numtabs],1,nchar(alltab[numtabs])-1)
    
    # multi-tissue plot write expression as character string
    tabexpr <- paste0("tabsetPanel( id =  'twascompare',",paste(alltab,collapse=""),")")
    
    
    # multi-tissue table, set tab panels functions as character strings
    tabpan1 <- rep("tabPanel('",numtabs+1)
    tabpan2 <- rep("', DT::dataTableOutput('",numtabs+1)
    tabpan3 <- rep("tbl')),",numtabs+1)
    alltabpantbl <- paste0(tabpan1,c(primary_tissue,tissue_list),tabpan2,c(primary_tissue,tissue_list),tabpan3)
    
    # multi-tissue table, cut comma from last entry
    alltabpantbl[numtabs+1] <- substring(alltabpantbl[numtabs+1],1,nchar(alltabpantbl[numtabs+1])-1)
    
    # multi-tissue table, write expression as character string
    tabexprtable <- paste0("  shiny::column(12, tabsetPanel(id = 'twasresult',",
                           paste(alltabpantbl,collapse=""),"))")
    
    
    # multi-tissue plot to print to the UI
    secondary_result1 <- h4(strong("Comparison of TWAS results from primary reference panel to results from secondary reference panel(s)"))
    secondary_result2 <- "Note: The figure in each tab displays a mirror plot of study results using the primary reference panel versus study results using a secondary reference panel."
    secondary_result3 <- eval(str2expression(tabexpr))
    secondary_result4 <- shiny::br()
    secondary_result5 <- shiny::hr()
    
    # multi-tissue table to print to the UI
    secondary_resulttab <- eval(str2expression(tabexprtable))
    
    
  } else if(multiple_tissues==FALSE) {
    
    # single-tissue plot to print to the UI (i.e. no plot)
    secondary_result1 <- ""
    secondary_result2 <- ""
    secondary_result3 <- ""
    secondary_result4 <- ""
    secondary_result5 <- ""
    
    
    # single tissue table to print to UI (i.e. just primary tissue table) 
    exprtable <-  "  shiny::column(12, DT::dataTableOutput('primarytissuetbl'))"
    secondary_resulttab <- eval(str2expression(exprtable))
  }
  
  # shiny::column(12,
  #               tabsetPanel(
  #                 id = 'twasresult',
  #                 tabPanel(primary_tissue, DT::dataTableOutput("DGNtbl")),
  #                 tabPanel("GWB", DT::dataTableOutput("GWBtbl")),
  #                 tabPanel("GTL", DT::dataTableOutput("GTLtbl")),
  #                 tabPanel("MSA", DT::dataTableOutput("MSAtbl"))
  #               )
  # )
  # 
  
  # ideogramtrack plot
  if (ideogram_present == TRUE){
    ideo1 <- plotOutput("chrplt",height = 100)
    ideo2 <- br()
  } else{
    ideo1 <- ""
    ideo2 <- ""
  }


  
  
  ####
  ############ UI setup ####################
  ####
  
  # define the UI
  ui <- #tagList(
    #useShinyjs(),
    navbarPage(title="LocusXcanR",
               tabPanel(title="Results",
                        fluidPage(
                          h4("Transcriptome-wide association study (TWAS)"),
                          h5(ref_expr_name),
                          h5(study_name),
                          HTML(head_details),
                          br(),
                          br(),
                          

                          br(),
                          fixedPanel(style="z-index:100;",
                                     top = 200, left = 10, right = 0, width = "60%", draggable = T,
                                     shiny::selectInput("locuslst","Select a locus to view (click and drag to reposition menu):",
                                                 choices = loclistdet,width = "550px")
                          ),
                          br(),
                          br(),
                          br(),
                          
                          
                          #### UI ideogram plot ####
                          ideo1,
                          ideo2,
                          
                          
                          #### UI mirror ####
                          h4(strong("TWAS-GWAS mirror plot of genes and variants within the locus")),
                          condtext,
                          radios_cond,
                          plotlyOutput("TWASmirror", height = 550),
                          br(),
                          hr(),
                          
                          
                          #### UI locus zoom ####
                          h4(strong("TWAS-GWAS mirror locus-zoom plot")),
                          paste0("Note: Top panel displays predicted expression correlation between index TWAS gene and other genes at the locus. Bottom panel displays LD between the index SNP and other SNPs at the locus. Lines connect genes to their predictive model variants. Color scale for genes denotes the degree of predicted expression correlation with the index gene. Color scale for SNPs and solid lines denotes the degree of LD with the index SNP. Dashed red line in the top panel denotes TWAS p-value threshold = ",formatC(10^(-pthresh), format = "e", digits = 2),", and in bottom panel denotes GWAS p-value threshold = 5.0e-8"),
                          br(),
                          plotlyOutput("TWAScorr", height=550),
                          br(),
                          
                          
                          #### UI network ####
                          h4(strong("Network visualization of TWAS results")),
                          "Sentinel TWAS gene is indicated by a star, all other genes are squares. Sentinel GWAS variant is indicated by a triangle, all other variants are circles. Color scale of all lines and shapes is based on correlation with the index gene or index variant. Line thickness corresponds to the model weight. Solid line indicates a positive direction of effect and dashed line indicates a negative direction. Size of the shape corresponds to the size of the -log10(p-value).",
                          br(),
                          visNetworkOutput("TWASnetwork", height=600),
                          hr(),
                          
                          
                          meta_result1,
                          meta_result2,
                          meta_result3,
                          meta_result4,
                          meta_result5,

                          
                          #### UI multi-tissue ####
                          secondary_result1,
                          secondary_result2,
                          secondary_result3,
                          secondary_result4,
                          secondary_result5,
                          
                          
                          #### UI TWAS tbl ####
                          h4(strong("Overall TWAS results from primary and secondary reference panel(s), if applicable, within the locus")),
                          "Note: Each tab represents TWAS results from a different gene expression reference panel. ",

                          shiny::fluidRow(
                            
                            #column(3,
                              # conditionalPanel(
                              #   'input.twasresult === "DGN"',
                              #   checkboxGroupInput("show_vars1", "Select DGN columns to view:",
                              #                      origvars, selected = origvars)
                              # ),
                              # conditionalPanel(
                              #   'input.twasresult === "GWB"',
                              #   checkboxGroupInput("show_vars2", "Select GWB columns to view:",
                              #                      names(twas_ds), selected = names(twas_ds))
                              # ),
                              # conditionalPanel(
                              #   'input.twasresult === "GTL"',
                              #   checkboxGroupInput("show_vars3", "Select GTL columns to view:",
                              #                      names(twas_ds), selected = names(twas_ds))
                              # ),
                              # conditionalPanel(
                              #   'input.twasresult === "MSA"',
                              #   checkboxGroupInput("show_vars4", "Select MSA columns to view:",
                              #                      names(twas_ds), selected = names(twas_ds))
                              # )
                            #  checkboxGroupInput("show_vars", "Select columns to view:",
                            #                     names(primary_ref_ds), selected = names(primary_ref_ds))
                            #),
                            
                            secondary_resulttab
                            
                          #   shiny::column(12,
                          #     tabsetPanel(
                          #       id = 'twasresult',
                          #       tabPanel(primary_tissue, DT::dataTableOutput("primary_ref_DT")),
                          #       tabPanel("GWB", DT::dataTableOutput("GWBtbl")),
                          #       tabPanel("GTL", DT::dataTableOutput("GTLtbl")),
                          #       tabPanel("MSA", DT::dataTableOutput("MSAtbl"))
                          #     )
                          #   )
                          ),
                          br(),
                          hr(),
                          
                          # h4(strong("Table 2. Overall TWAS results for other traits in the trait category from primary and secondary reference panels within the locus")),
                          # "Note: DGN = Depression Genes and Networks, GWB = GTEx whole blood, GTL = GTEx EBV transformed lymphocytes, MSA = MESA monocytes; each represents a gene expression reference panel. ",
                          # HTML('<span style="background-color:lightgreen">Significant gene-trait associations highlighted in green, </span><span style="background-color:tomato">HLA genes / MHC regions / single SNP models highlighted in red. </span>'),
                          # " MHC region is defined as GRCh37; chr6:28,477,797-33,448,354. Single SNP model indicates that the predictive expression model for the gene contained only a single SNP.",
                          # tabsetPanel(
                          #   id = 'twasresulttrt',
                          #   tabPanel("DGN", DT::dataTableOutput("primary_ref_tbltrt")),
                          #   tabPanel("GWB", DT::dataTableOutput("GWBtbltrt")),
                          #   tabPanel("GTL", DT::dataTableOutput("GTLtbltrt")),
                          #   tabPanel("MSA", DT::dataTableOutput("MSAtbltrt"))
                          # ),
                          # br(),
                          # hr(),
                          
                          
                          #### UI known GWAS tbl ####
                          h4(strong("Reported GWAS sentinel variants within the locus")),
                          h5(HTML('<span style="background-color:tomato">Table <strong>row</strong> is highlighted in red if the SNP is not in the study imputation data</span>')),
                          h5(HTML('<span style="background-color:lightgreen">Variants (<strong>chrposall</strong>) reported in figures are highlighted in green</span>')),
                          h5("GWAS ",strong("genename"),":"),
                          p(HTML('<ul>
                                  <li><span style="background-color:lightgreen">matching significant TWAS genes are highlighted in green</span></li>
                                  <li><span style="background-color:yellow"    >matching non-significant TWAS genes are highlighted in yellow</span></li>
                                  <li><span style="background-color:orange"    >included in primary reference panel but not predicted in study are highlighted in orange</span></li>
                                  <li><span style="background-color:tomato"    >included in primary reference panel but no predictive model was fit for the gene are highlighted in red</span></li>
                                  <li>not included in primary reference panel (i.e. no info is available in this analysis) are not highlighted</li>
                                  </ul>
                                  ')),
                              br(),
                          DT::dataTableOutput("KnownSNPtbl"),
                          br(),
                          hr(),
                          
                          
                          #### UI GWAS variants ####
                          h4(strong("Study GWAS results displayed in figures as variants \"Reported in GWAS\"")),
                          h5("Note: Results listed below are from trait specific GWAS"),
                          DT::dataTableOutput("GWASvars"),
                          br(),
                          
                          
                          #### UI weight tbl ####
                          h4(strong("Gene expression prediction model weights")),
                          h5("Note: \"LD\" refers to the LD between the most significant rsid at the locus with the variant listed in the row. \"Corr\" refers to the Pearson correlation between the most significant gene at the locus and the gene listed in the row."),
                          DT::dataTableOutput("weighttbl"),
                          br(),
                          
                          add_UI
                        )
               ),
               
               ######### Methods UI ################################################################################    
               
               shiny::tabPanel("Methods",
                        method_details
               )
    )
#  )



####
####### Server setup #####################################################################################
####


  
  # give instructions to the server
  server <- function(input,output){
    
    # define reactive variables for all plots
    locds <- reactive({
      filter(primary_ref_ds, locvar==input$locuslst)
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
    locpheno <- reactive({ unique(locds()$pheno) })
    #locphcat <- reactive({ unique(locds()$phenocat) })
    loctitle <- reactive({ paste0("chr ",locchr(),": ",xlow()," - ",xhigh(),"; trait = ",locpheno())
      })
    
    # set the dataset to extract known variants at the locus
    # ds <- reactive({
    #   if (locphcat()=="RBC") {
    #     GWAS_sentinel
    #   } else if (locphcat()=="WBC") {
    #     wbcds
    #   } else {
    #     pltds
    #   }
    # })
    
    ds <- GWAS_sentinel
    
    # select all known GWAS genes at locus
    # phenotbl <- reactive({
    #   tmp <- filter(ds(),(X6)>=xlow() & xhigh()>=(X6) & X5==locchr()) %>% select(X8)
    #   colnames(tmp) <- c("Gene")
    #   tmp
    # })
    
    phenotbl <- reactive({
      tmp <- filter(ds,(position)>=xlow() & xhigh()>=(position) & chr==locchr()) %>% select(genename)
      #colnames(tmp) <- c("Gene")
      tmp
    })
    
    
    # select all TWAS genes at locus
    primary_ref_tbl <-  reactive({
      tmp <- filter(primary_ref_ds,genestart>=xlow() & genestop<=xhigh() & pheno==locpheno() & chr==locchr())
      
      if (nrow(phenotbl())==0){
        tmp$kngene="Not Reported in GWAS"
      } else{
        tmp$kngene <- ifelse(grepl(paste(unique(phenotbl()$genename), collapse="|"),tmp$genename),
                             "Reported in GWAS","Not Reported in GWAS")
      }
      tmp
    })
    

    ############# Ideogram ###############################################

    # ideogram plot
    output$chrplt <- renderPlot({

      itrack = Gviz::IdeogramTrack(genome = genome_build, chromosome = paste0("chr",locchr()), bands = data)
      trackplot = c(itrack,Gviz::GenomeAxisTrack());

      Gviz::plotTracks(c(trackplot), from = xlow(),to=xhigh(),transcriptAnnotation="symbol",
                 add53=TRUE,showBandID=TRUE,cex.bands=2,stackHeight=2,background.title = "white",
                 col.axis="black",col.title="black",cex.title=2,
                 cex.axis=4,just.group="below",cex.main = 1.5,cex=1.5)
    })
    
    ############# TWAS/GWAS mirror ###############################################
    
    
    
    # TWAS/GWAS mirror plot
    output$TWASmirror <- renderPlotly({
      #locds <- primary_ref_ds %>% filter(locvar==input$locuslst)
      yhigh <- max(locds()$log10pval)+0.25*max(locds()$log10pval)


      # select significant and known genes to plot
      primary_ref_tblplt <- primary_ref_tbl() %>% filter(SignifGene==1 | kngene=="Reported in GWAS")
      
      # select the GWAS results at the locus
      gwasloc <- gwasallfin %>% filter(LOCUS==locnum())

      if (nrow(cohort_gwas_knownfin[cohort_gwas_knownfin$Locus==locnum(),])==0){
        gwasloc$knsnp <- "Not reported in GWAS"
      } else {
        gwasloc$knsnp <- ifelse(grepl(paste(unique(cohort_gwas_knownfin$SNP[cohort_gwas_knownfin$Locus==locnum()]),collapse="|"),
                                      gwasloc$SNP),"Reported in GWAS","Not reported in GWAS")
      }
      
      # not known variants at the locus pval < 0.1
      gwaslocnotkn <- gwasloc %>% filter(knsnp=="Not reported in GWAS")
      
      # known variants at the locus, regardless of pvalue
      knowngwas <- cohort_gwas_knownfin %>% filter(Locus==locnum())
      
      if (nrow(knowngwas)>0){
        knowngwas$knsnp <- "Reported in GWAS"
        ylow <- min(gwaslocnotkn$poslog10p,log10(5*10^(-8)),knowngwas$poslog10p)+0.15*min(gwaslocnotkn$poslog10p,log10(5*10^(-8)),knowngwas$poslog10p)
      } else{
        ylow <- min(gwaslocnotkn$poslog10p,log10(5*10^(-8)))+0.15*min(gwaslocnotkn$poslog10p,log10(5*10^(-8)))
      }
      
      nudgeval <- 0.05*max(abs(yhigh))
      nudgevalg <- 0.05*ylow
      
      # define specific colors for known/not known genes
      if(length(unique(primary_ref_tblplt$kngene))==1){
        if(unique(primary_ref_tblplt$kngene)=="Not reported in GWAS"){
          myColors <- setNames( c('#56B4E9','#000000'),
                                c("Not reported in GWAS","Reported in GWAS") )
        } else if(unique(primary_ref_tblplt$kngene)=="Reported in GWAS"){
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
      
      if (conditional_present==TRUE){
        if (input$margcondradio==1){
          pvalplt <- primary_ref_tblplt$log10pval
        }
        else if(input$margcondradio==2){
          pvalplt <- -log10(primary_ref_tblplt$p_final)
        }
      } else{
          pvalplt <- primary_ref_tblplt$log10pval
      }
      
      #colnames(primary_ref_tblplt)[colnames(primary_ref_tblplt) == "genemidMB"] <- "gene midpoint"
      
      # TWAS plot
      p <- ggplotly(ggplot(data=primary_ref_tblplt, aes(x=genemidMB, y=pvalplt, 
                                                        color=as.factor(kngene),
                                                        text = paste("Gene: ",genename,
                                                                     "<br>-log10 p-value: ",round(pvalplt,3)))) + 
                      geom_point(pch=15) +
                      geom_segment(aes(x = genestartMB, y = pvalplt, xend = genestopMB, yend = pvalplt, color=as.factor(kngene)), size=2) +
                      geom_text(aes(x=genemidMB,y=pvalplt,label=genename,color=as.factor(kngene)), nudge_y = nudgeval, size=4) +
                      geom_hline(aes(yintercept=pthresh), lty=2, color="red") +
                      xlim(xlowMB(),xhighMB()) + ylim(0,yhigh) +
                      annotate(geom="text",x=round((xlow()+270000)/1000000,4), y=pthresh-nudgeval, color='red',
                               label=paste0("TWAS p-value: ",formatC(10^(-pthresh), format = "e", digits = 2)),size=4) +
                      theme(legend.position = 'top', legend.title = element_blank()) +
                      theme_bw() +
                      colScale
      , tooltip=c("text"))
      
      p <- p %>% plotly::layout(title = list(text=loctitle(), x=0, xanchor='left', y=.99),
                        yaxis = list(title = ' TWAS -log10(p)'),
                        legend = list(orientation='h', x=0, y=1))
      
      
      plta <- ggplot(data=gwaslocnotkn, aes(x=round(as.numeric(pos)/1000000,4),y=poslog10p, 
                                            color=as.factor(knsnp),
                                            text = paste("position: ",as.numeric(pos),
                                                         "<br>log10 p-value: ",round(poslog10p,3)))) + 
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
      m <- ggplotly(app_plt, tooltip=c("text"))
      
      m <- m %>% plotly::layout(xaxis = list(title="position (in Mb)"), 
                        yaxis=list(title="GWAS log10(p)"))
      
      # combine plots together
      fig <- subplot(style(p,showlegend=F), m, shareX = T, nrows = 2, titleX = T,titleY = T, 
                     which_layout = 1, margin=0.005)
      
    })
    
    
    ############# Locus zoom ###############################################
    
    
    # TWAS correlation plots
    output$TWAScorr <- renderPlotly({
      
      # select significant and known genes to plot
      primary_ref_tblplt <- primary_ref_tbl() %>% filter(SignifGene==1 | kngene=="Reported in GWAS")
      uniqgenes <- primary_ref_tblplt %>% select(genename)
      indexgene <- primary_ref_tblplt$genename[primary_ref_tblplt$p == min(primary_ref_tblplt$p)]
      
      
      # extract set of genes from correlation matrix
      Mindex <- data.table(M[uniqgenes$genename,indexgene])
      colnames(Mindex) <- c("corr")
      Mindex$genename <- uniqgenes$genename
      
      
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
      corplt <- merge(primary_ref_tblplt,Mindex,by.x="genename",by.y="genename")
      
      # get the primary_ref weights for the specific genes at locus
      primary_ref_wtloc <- weight_ds %>% filter(genename %in% corplt$genename)
      
      # subset the GWAS variants for specific chr, and phenotype
      gwasallfinloc <- gwasallfin %>% filter(LOCUS==locnum()) %>%
        select(pos,all1,all2,poslog10p)
      gwasallfinloc$posnum <- as.numeric(gwasallfinloc$pos)
      
      # get the matching GWAS variants for the weights
      primary_ref_wtgwasloc1 <- merge(primary_ref_wtloc,gwasallfinloc,by.x=c('position','ref_allele','eff_allele'),
                             by.y=c('posnum','all1','all2'))
      primary_ref_wtgwasloc2 <- merge(primary_ref_wtloc,gwasallfinloc,by.x=c('position','ref_allele','eff_allele'),
                             by.y=c('posnum','all2','all1'))
      primary_ref_wtgwasloc <- rbind(primary_ref_wtgwasloc1,primary_ref_wtgwasloc2)
      
      # merge GWAS results with TWAS info
      TWASloc <- corplt %>% select(genename,genestartMB,genemid,genestopMB,log10pval,corgroup,corr)
      TWASloc$genemidMB <- round(TWASloc$genemid/1000000,4)
      primary_ref_wtgwaslocfin <- merge(primary_ref_wtgwasloc,TWASloc, by.x='genename',by.y='genename',all.x=T)
      
      # set y limit
      yhigh <- max(corplt$log10pval)+0.25*max(corplt$log10pval)
      ylow <- min(primary_ref_wtgwasloc$poslog10p,log10(5*10^(-8)))+0.15*min(primary_ref_wtgwasloc$poslog10p,log10(5*10^(-8)))
      nudgeval <- 0.05*max(abs(yhigh))
      
      # set color scale values
      colval <- c("[0 - 0.2]"="#CCCCCC","(0.2 - 0.4]"="#333FFF",
                  "(0.4 - 0.6]"="#00FF00","(0.6 - 0.8]"="#FF9900","(0.8 - 1]"="#FF0000")
      breakval <- c("[0 - 0.2]","(0.2 - 0.4]","(0.4 - 0.6]","(0.6 - 0.8]","(0.8 - 1]")
      
      
      # select the LD for the given locus
      LDdsloc <- as.data.table(filter(LDds,locus==locnum()))
      
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
      primary_ref_wtgwaslocfin2 <- merge(primary_ref_wtgwaslocfin,LDdsloc, by.x=c('position'),by.y=c('posb'),all.x=T)
      
      primary_ref_wtgwaslocfin2_comp <- primary_ref_wtgwaslocfin2[complete.cases(primary_ref_wtgwaslocfin2), ]
      primary_ref_wtgwaslocfin2_comp$line <- ifelse(primary_ref_wtgwaslocfin2_comp$weight<0,2,1)
      
      xlowMB <- min(corplt$genestartMB,round(primary_ref_wtgwaslocfin2_comp$position/1000000,4))-.0005
      xhighMB <- max(corplt$genestopMB,round(primary_ref_wtgwaslocfin2_comp$position/1000000,4))+.0005
      

      # plot correlation categories at the locus, like locus zoom plot
      loczoom <- ggplotly(ggplot(data=corplt, aes(x=round(genemid/1000000,4), y=log10pval,color=corgroup,
                                                  text = paste("Gene: ",genename,
                                                               "<br>-log10 p-value: ",round(log10pval,3),
                                                               "<br>correlation: ",round(corr,3)))) +
                            geom_hline(aes(yintercept=pthresh), lty=2, color="red") +
                            geom_hline(aes(yintercept=log10(5*10^(-8))), lty=2, color="red") +
                            geom_hline(aes(yintercept=0),size=1,color='black') +
                            
                            geom_segment(data=primary_ref_wtgwaslocfin2_comp,aes(x=round(position/1000000,4),y=poslog10p,xend=genemidMB,yend=log10pval
                                                                        ,color=ldgroup)) +
                            geom_point(data=primary_ref_wtgwaslocfin2_comp,
                                       aes(color=corgroup), pch=15) +
                            geom_segment(aes(x=genestartMB,y=log10pval,xend=genestopMB,yend=log10pval, color=corgroup),size=2) +
                            geom_text(aes(label=genename,color=corgroup),nudge_y = nudgeval) +
                            scale_colour_manual(breaks = breakval, values = colval) +
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
                            geom_point(data=primary_ref_wtgwaslocfin2_comp, aes(x=round(position/1000000,4),y=poslog10p,color=ldgroup,
                                                                                text = paste("Position: ",position,
                                                                                             "<br>log10 p-value: ",round(poslog10p,3),
                                                                                             "<br>LD: ",round(corrab,3)))) #+
      , tooltip=c("text"))
      
      loczoom <- loczoom %>% plotly::layout(title = list(text=loctitle(), x=0, xanchor='left', y=.99),
                                    yaxis = list(title = '<-- GWAS log10(p) ... 0 ... TWAS -log10(p) -->'),
                                    legend = list(orientation='h', x=0, y=1),
                                    xaxis = list(title="position (in Mb)"))
    })
    
    
    ############# Network ###############################################
    
    
    # TWAS network plots
    output$TWASnetwork <- renderVisNetwork({
      
      # select significant and known genes to plot
      primary_ref_tblplt <- primary_ref_tbl() %>% filter(SignifGene==1 | kngene=="Reported in GWAS")
      uniqgenes <- primary_ref_tblplt %>% select(genename)
      indexgene <- primary_ref_tblplt$genename[primary_ref_tblplt$p == min(primary_ref_tblplt$p)]
      
      # extract set of genes from correlation matrix
      Mindex <- data.table(M[uniqgenes$genename,indexgene])
      colnames(Mindex) <- c("corr")
      Mindex$genename <- uniqgenes$genename
      
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
      corplt <- merge(primary_ref_tblplt,Mindex,by.x="genename",by.y="genename")
      
      # get the primary_ref weights for the specific genes at locus
      primary_ref_wtloc <- weight_ds %>% filter(genename %in% corplt$genename)
      
      # subset the GWAS variants for specific chr, and phenotype
      gwasallfinloc <- gwasallfin %>% filter(LOCUS==locnum()) %>%
        select(pos,all1,all2,poslog10p)
      gwasallfinloc$posnum <- as.numeric(gwasallfinloc$pos)
      
      # get the matching GWAS variants for the weights
      primary_ref_wtgwasloc1 <- merge(primary_ref_wtloc,gwasallfinloc,by.x=c('position','ref_allele','eff_allele'),
                             by.y=c('posnum','all1','all2'))
      primary_ref_wtgwasloc2 <- merge(primary_ref_wtloc,gwasallfinloc,by.x=c('position','ref_allele','eff_allele'),
                             by.y=c('posnum','all2','all1'))
      primary_ref_wtgwasloc <- rbind(primary_ref_wtgwasloc1,primary_ref_wtgwasloc2)
      
      # merge GWAS results with TWAS info
      TWASloc <- corplt %>% select(genename,genestartMB,genemid,genestopMB,log10pval,corgroup)
      TWASloc$genemidMB <- round(TWASloc$genemid/1000000,4)
      primary_ref_wtgwaslocfin <- merge(primary_ref_wtgwasloc,TWASloc, by.x='genename',by.y='genename',all.x=T)
      
      # set y limit
      yhigh <- max(corplt$log10pval)+0.15*max(corplt$log10pval)
      ylow <- min(primary_ref_wtgwasloc$poslog10p)+0.15*min(primary_ref_wtgwasloc$poslog10p)
      nudgeval <- 0.05*max(abs(yhigh))
      
      # set color scale values
      colval <- c("[0 - 0.2]"="#CCCCCC","(0.2 - 0.4]"="#333FFF",
                  "(0.4 - 0.6]"="#00FF00","(0.6 - 0.8]"="#FF9900","(0.8 - 1]"="#FF0000")
      
      # select the LD for the given locus
      LDdsloc <- as.data.table(filter(LDds,locus==locnum()))
      
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
      primary_ref_wtgwaslocfin2 <- merge(primary_ref_wtgwaslocfin,LDdsloc, by.x=c('position'),by.y=c('posb'),all.x=T)
      
      primary_ref_wtgwaslocfin2_comp <- primary_ref_wtgwaslocfin2[complete.cases(primary_ref_wtgwaslocfin2), ]
      primary_ref_wtgwaslocfin2_comp$line <- ifelse(primary_ref_wtgwaslocfin2_comp$weight<0,2,1)
      
      
      # data for links
      from <- primary_ref_wtgwaslocfin2$rsid
      to <- primary_ref_wtgwaslocfin2$genename
      weight <- primary_ref_wtgwaslocfin2$weight
      color <- primary_ref_wtgwaslocfin2$ldcol
      length <- 1/(primary_ref_wtgwaslocfin2$log10pval-primary_ref_wtgwaslocfin2$poslog10p)*1000
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
      
      
      nodessnp <- primary_ref_wtgwaslocfin2 %>% select(rsid,poslog10p,ldcol)
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
    
    
    ############# Meta analysis ###############################################
    
    # TWAS/GWAS mirror plot of meta-analysis
    output$Tmetamirror <- renderPlotly({
      
      # select significant and known genes to plot
      primary_ref_tblplt <- primary_ref_tbl() %>% filter(SignifGene==1 | kngene=="Reported in GWAS")
      
      
      # define y limits and nudges for plotting
      yhigh <- max(locds()$log10pval)+0.25*max(locds()$log10pval)
      ylow <- min(-primary_ref_tblplt$log10pvalmeta,metathresh)+0.15*min(-primary_ref_tblplt$log10pvalmeta,metathresh)
      nudgeval <- 0.05*max(abs(yhigh))
      nudgevalg <- 0.05*ylow
      
      
      # define specific colors for known/not known genes
      myColors <- setNames( c('#56B4E9', '#000000', '#E69F00'),
                            c("Not reported in GWAS","Reported in GWAS",NA))
      
      colScale <- scale_colour_manual(values = myColors)
      
      
      # TWAS plot
      p <- ggplotly(ggplot(data=primary_ref_tblplt, 
                           aes(x=round(genemid/1000000,4), y=log10pval, color=as.factor(kngene),
                               text = paste("Gene: ",genename,
                                            "<br>-log10 p-value: ",round(log10pval,3)))) +
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
      , tooltip = c("text"))
      
      p <- p %>% plotly::layout(title = list(text=loctitle(), x=0, xanchor='left', y=.99),
                        yaxis = list(title = ' TWAS -log10(p)'),
                        legend = list(orientation='h', x=0, y=1))
      
      # Meta-analysis plot
      m <- ggplotly(ggplot(data=primary_ref_tblplt, 
                           aes(x=round(genemid/1000000,4), y=-log10pvalmeta, 
                               color=as.factor(kngene),
                               text = paste("Gene: ",genename,
                                            "<br>log10 p-value: ",round(-log10pvalmeta,3)))) +
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
      , tooltip = c("text"))
      
      m <- m %>% plotly::layout(xaxis = list(title="position (in Mb)"),
                        yaxis=list(title="Meta-analysis log10(p)"))
      
      # combine plots together
      fig <- subplot(style(p,showlegend=F), m, shareX = T, nrows = 2, titleX = T,titleY = T,
                     which_layout = 1, margin=0.005)
      
    })
    
    
    ############# Ref panel compare ###############################################
    
    
    ######### NEED TO EDIT CODE HERE FOR comparing reference panels #####
    if (multiple_tissues == TRUE){
      
      n <- length(tissue_list)
      #varnames <- vector(mode='list', length=n)
      #varnames  <- paste("output$CompRef", refpanelist[1:n], sep="")
      varnames <- as.list(paste0("output$CompRef",tissue_list))
      
      #refunc <- paste("renderfunc", refpanelist[1:n], sep="")
      refunc <- vector(mode='list', length=n)
      
      for (i in seq(1,n)){
        refunc[[i]] <- renderPlotly({
          
          # select reference panel specific results
          primary_ref_tbl <- primary_ref_ds %>% filter(genestart>=xlow() & genestop<=xhigh() & pheno==locpheno() & 
                                                         chr==locchr() ) %>%
            select(genename,chr,pheno,genestart,genestop,genemid,log10pval)
          
          locrefpanel <- twas_ds %>% filter(.data$tissue==tissue_list[i] & genestart>=xlow() & genestop<=xhigh() & chr==locchr() & 
                                              pheno==locpheno()) %>%
            select(genename,chr, pheno,log10pval,genestart,genestop,genemid, SignifGene)
          
          pthreshrefpanel <- -log10(0.05/(nrow(twas_ds[twas_ds$tissue==tissue_list[i],])))
          
          #####################
          
          # merge secondary reference panel data sets with primary_ref
          primary_ref_refpanel <- merge(primary_ref_tbl,locrefpanel, by=c("genename","chr","pheno"), all = T)
          primary_ref_refpanel$inboth <- ifelse(complete.cases(primary_ref_refpanel$log10pval.x,
                                                               primary_ref_refpanel$log10pval.y),"In Both","Not in Both")
          
          #####################
          
          
          # set y limit values for each plot
          yhighprimary_ref_ <- max(primary_ref_refpanel$log10pval.x, na.rm = T)+0.15*max(primary_ref_refpanel$log10pval.x, na.rm = T)
          nudgevalprimary_ref_ <- 0.05*yhighprimary_ref_
          
          yhighrefpanel <- max(primary_ref_refpanel$log10pval.y,pthreshrefpanel, na.rm = T)+0.15*max(primary_ref_refpanel$log10pval.y,
                                                                                                     pthreshrefpanel,na.rm = T)
          nudgevalrefpanel <- 0.05*yhighrefpanel
          
          myColors <- setNames( c('#000000', '#56B4E9', '#E69F00'),
                                c("In Both","Not in Both",NA))
          colScale <- scale_colour_manual(values = myColors)
          
          
          ####################
          
          # plot info for ref panel
          atop <- ggplotly(ggplot(data=primary_ref_refpanel, 
                                  aes(x=round(genemid.x/1000000,4),
                                      y=log10pval.x, color=inboth,
                                      text = paste("Gene: ",genename,
                                                   "<br>-log10 p-value: ",round(log10pval.x,3)))) +
                             geom_point(pch=15) +
                             geom_segment(aes(x=round(genestart.x/1000000,4),y=log10pval.x, xend=round(genestop.x/1000000,4), 
                                              yend=log10pval.x), size=2) +
                             geom_text(aes(label=genename), nudge_y = nudgevalprimary_ref_) +
                             geom_hline(aes(yintercept=pthresh), lty=2, color='red') +
                             xlim(round(xlow()/1000000,4),round(xhigh()/1000000,4)) +
                             ylim(0,yhighprimary_ref_) +
                             theme_bw() +
                             theme(legend.position = 'top', legend.title = element_blank()) +
                             annotate(geom="text",x=round(xlow()/1000000,4)+.370000, y=pthresh+nudgevalprimary_ref_, 
                                      color='red',
                                      label=paste0(primary_tissue," TWAS p-value: ", formatC(10^-(pthresh), format = "e",
                                                                                             digits = 2)),size=4) +
                             
                             colScale
                           , tooltip = "text")
          
          atop <- atop %>% plotly::layout(yaxis = list(title = paste0(primary_tissue,' TWAS -log10(p)')),
                                          xaxis = list(range=c(round(xlow()/1000000,4),round(xhigh()/1000000,4))),
                                          annotations=list(x = 0.5 , y = 1.1, text = paste0("(a) ",primary_tissue," vs. ",tissue_list[i]), 
                                                           showarrow = F, 
                                                           xref='paper', yref='paper',xanchor='center'),
                                          legend = list(orientation='h', x=0, y=1),
                                          title=list(text=paste0(loctitle(),"\n"),x=0,xanchor='left'), margin=list(t=60)
          )
          
          #####
          
          abottom <- ggplotly(ggplot(data=primary_ref_refpanel, 
                                     aes(x=round(genemid.y/1000000,4),
                                         y=-log10pval.y, color=inboth,
                                         text = paste("Gene: ",genename,
                                                      "<br>log10 p-value: ",round(-log10pval.y,3)))) +
                                geom_hline(aes(yintercept=-pthreshrefpanel), lty=2, color='red') +
                                geom_point(pch=15) +
                                geom_segment(aes(x=round(genestart.y/1000000,4),y=-log10pval.y, 
                                                 xend=round(genestop.y/1000000,4), yend=-log10pval.y), size=2) +
                                geom_text(aes(label=genename), nudge_y = -1.5*nudgevalrefpanel) +
                                xlim(round(xlow()/1000000,4),round(xhigh()/1000000,4)) +
                                ylim(-yhighrefpanel,0) +
                                theme_bw() +
                                theme(legend.position = 'top', legend.title = element_blank()) +
                                annotate(geom="text",x=round(xlow()/1000000,4)+.370000, y=-pthreshrefpanel-nudgevalrefpanel, 
                                         color='red',
                                         label=paste0(tissue_list[i]," TWAS p-value: ", formatC(10^-(pthreshrefpanel), format = "e", 
                                                                                                digits = 2)),size=4) +
                                colScale
                              , tooltip = "text")
          
          abottom <- abottom %>% plotly::layout(yaxis = list(title = paste0(tissue_list[i],' TWAS log10(p)'), 
                                                             range=c(-yhighrefpanel,0.1)),
                                                xaxis = list(range=c(round(xlow()/1000000,4),round(xhigh()/1000000,4)), 
                                                             title=paste0("position (in Mb)",i)))
          
          #####
          
          afig <- subplot(atop,style(abottom,showlegend=F),nrows=2,shareX = T, titleX = T, titleY = T, which_layout = 1, 
                          margin=.01)
          
        })
        
        eval(parse(text=paste0("output$CompRefGTL <- refunc[[",i,"]]")))
        exprtext <- "output$CompRefGWB <- refunc[[2]]"
        eval(str2expression(exprtext))
        eval(parse(text="output$CompRefMSA <- refunc[[3]]"))
  
        
      }
    }
    
    
    ######### begin Commenting out here for testing ####################
    # 
    # # comparison of TWAS primary_ref results with results from other reference panels
    # output$CompRefGWB <- renderPlotly({
    #   
    #   # select primary_ref results
    #   #locds <- primary_ref_ds %>% filter(locvar==input$locuslst)
    #   #xhigh <- max(locds$genestop)+1000000
    #   #xlow <- max(0,min(locds$genestart)-1000000)
    #   #locchr <- unique(locds$chr)
    # 
    #   
    #   #####################
    #   
    #   
    #   # select reference panel specific results
    #   primary_ref_tbl <- primary_ref_ds %>% filter(genestart>=xlow() & genestop<=xhigh() & pheno==locpheno() & 
    #                                                  chr==locchr() #& 
    #                                #is.na(HLARegion) & is.na(MHCRegion)
    #                                ) %>%
    #     select(genename,chr,pheno,genestart,genestop,genemid,log10pval)
    #   
    #   locgwb <- twas_ds %>% filter(tissue=="GWB" & genestart>=xlow() & genestop<=xhigh() & chr==locchr() & 
    #                                   pheno==locpheno() #& is.na(HLARegion) & is.na(MHCRegion)
    #                                ) %>%
    #     select(genename,chr, pheno,log10pval,genestart,genestop,genemid, SignifGene)
    #   pthreshgwb <- -log10(0.05/(nrow(twas_ds[twas_ds$tissue=="GWB",])))
    #   
    #   # locgtl <- twas_ds %>% filter(tissue=="GTL" & genestart>=xlow() & genestop<=xhigh() & chr==locchr() & 
    #   #                                 pheno==locpheno() #& is.na(HLARegion) & is.na(MHCRegion)
    #   #                              ) %>%
    #   #   select(genename,chr, pheno,log10pval,genestart,genestop,genemid, SignifGene)
    #   # pthreshgtl <- -log10(0.05/(nrow(twas_ds[twas_ds$tissue=="GTL",])))
    #   # 
    #   # locmsa <- twas_ds %>% filter(tissue=="MSA" & genestart>=xlow() & genestop<=xhigh() & chr==locchr() & 
    #   #                                 pheno==locpheno() #& is.na(HLARegion) & is.na(MHCRegion)
    #   #                              ) %>%
    #   #   select(genename,chr, pheno,log10pval,genestart,genestop,genemid, SignifGene)
    #   # pthreshmsa <- -log10(0.05/(nrow(twas_ds[twas_ds$tissue=="MSA",])))
    #   
    #   
    #   #####################
    #   
    #   
    #   # merge secondary reference panel data sets with primary_ref
    #   primary_ref_gwb <- merge(primary_ref_tbl,locgwb, by=c("genename","chr","pheno"), all = T)
    #   primary_ref_gwb$inboth <- ifelse(complete.cases(primary_ref_gwb$log10pval.x,primary_ref_gwb$log10pval.y),"In Both","Not in Both")
    #   
    #   # primary_ref_gtl <- merge(primary_ref_tbl,locgtl, by=c("genename","chr","pheno"), all = T)
    #   # primary_ref_gtl$inboth <- ifelse(complete.cases(primary_ref_gtl$log10pval.x,primary_ref_gtl$log10pval.y),"In Both","Not in Both")
    #   # 
    #   # primary_ref_msa <- merge(primary_ref_tbl,locmsa, by=c("genename","chr","pheno"), all = T)
    #   # primary_ref_msa$inboth <- ifelse(complete.cases(primary_ref_msa$log10pval.x,primary_ref_msa$log10pval.y),"In Both","Not in Both")
    #   
    #   
    #   #####################
    #   
    #   
    #   # set y limit values for each plot
    #   yhighprimary_ref_ <- max(primary_ref_gwb$log10pval.x, na.rm = T)+0.15*max(primary_ref_gwb$log10pval.x, na.rm = T)
    #   nudgevalprimary_ref_ <- 0.05*yhighprimary_ref_
    #   
    #   # yhighgtl <- max(primary_ref_gtl$log10pval.y,pthreshgtl, na.rm = T)+0.15*max(primary_ref_gtl$log10pval.y,pthreshgtl, na.rm = T)
    #   # nudgevalgtl <- 0.05*yhighgtl
    #   
    #   yhighgwb <- max(primary_ref_gwb$log10pval.y,pthreshgwb, na.rm = T)+0.15*max(primary_ref_gwb$log10pval.y,pthreshgwb,
    #                                                                               na.rm = T)
    #   nudgevalgwb <- 0.05*yhighgwb
    #   
    #   # yhighmsa <- max(primary_ref_msa$log10pval.y,pthreshmsa, na.rm = T)+0.15*max(primary_ref_msa$log10pval.y,pthreshmsa, na.rm = T)
    #   # nudgevalmsa <- 0.05*yhighmsa
    #   
    #   
    #   myColors <- setNames( c('#000000', '#56B4E9', '#E69F00'),
    #                         c("In Both","Not in Both",NA))
    #   
    #   colScale <- scale_colour_manual(values = myColors)
    #   
    #   
    #   ####################
    #   
    #   # plot info for GWB
    #   atop <- ggplotly(ggplot(data=primary_ref_gwb, 
    #                           aes(x=round(genemid.x/1000000,4),
    #                               y=log10pval.x, color=inboth,
    #                               text = paste("Gene: ",genename,
    #                                            "<br>-log10 p-value: ",round(log10pval.x,3)))) +
    #                      geom_point(pch=15) +
    #                      geom_segment(aes(x=round(genestart.x/1000000,4),y=log10pval.x, xend=round(genestop.x/1000000,4), 
    #                                       yend=log10pval.x), size=2) +
    #                      geom_text(aes(label=genename), nudge_y = nudgevalprimary_ref_) +
    #                      geom_hline(aes(yintercept=pthresh), lty=2, color='red') +
    #                      xlim(round(xlow()/1000000,4),round(xhigh()/1000000,4)) +
    #                      ylim(0,yhighprimary_ref_) +
    #                      theme_bw() +
    #                      theme(legend.position = 'top', legend.title = element_blank()) +
    #                      annotate(geom="text",x=round(xlow()/1000000,4)+.370000, y=pthresh+nudgevalprimary_ref_, 
    #                               color='red',
    #                               label=paste0("Primary ref TWAS p-value: ", formatC(10^-(pthresh), format = "e",
    #                                                                                  digits = 2)),size=4) +
    #                      
    #                      colScale
    #   , tooltip = "text")
    #   
    #   atop <- atop %>% plotly::layout(yaxis = list(title = 'Primary ref TWAS -log10(p)'),
    #                           xaxis = list(range=c(round(xlow()/1000000,4),round(xhigh()/1000000,4))),
    #                           annotations=list(x = 0.5 , y = 1.1, text = paste0("(a) ",primary_tissue," vs. GWB"), showarrow = F, 
    #                                            xref='paper', yref='paper',xanchor='center'),
    #                           legend = list(orientation='h', x=0, y=1),
    #                           title=list(text=paste0(loctitle(),"\n"),x=0,xanchor='left'), margin=list(t=60)
    #   )
    #   
    #   #####
    #   
    #   abottom <- ggplotly(ggplot(data=primary_ref_gwb, 
    #                              aes(x=round(genemid.y/1000000,4),
    #                                  y=-log10pval.y, color=inboth,
    #                                  text = paste("Gene: ",genename,
    #                                               "<br>log10 p-value: ",round(-log10pval.y,3)))) +
    #                         geom_hline(aes(yintercept=-pthreshgwb), lty=2, color='red') +
    #                         geom_point(pch=15) +
    #                         geom_segment(aes(x=round(genestart.y/1000000,4),y=-log10pval.y, 
    #                                          xend=round(genestop.y/1000000,4), yend=-log10pval.y), size=2) +
    #                         geom_text(aes(label=genename), nudge_y = -1.5*nudgevalgwb) +
    #                         xlim(round(xlow()/1000000,4),round(xhigh()/1000000,4)) +
    #                         ylim(-yhighgwb,0) +
    #                         theme_bw() +
    #                         theme(legend.position = 'top', legend.title = element_blank()) +
    #                         annotate(geom="text",x=round(xlow()/1000000,4)+.370000, y=-pthreshgwb-nudgevalgwb, 
    #                                  color='red',
    #                                  label=paste0("GWB TWAS p-value: ", formatC(10^-(pthreshgwb), format = "e", 
    #                                                                             digits = 2)),size=4) +
    #                         colScale
    #   , tooltip = "text")
    #   
    #   abottom <- abottom %>% plotly::layout(yaxis = list(title = 'GWB TWAS log10(p)', 
    #                                              range=c(-yhighgwb,0.1)),
    #                                 xaxis = list(range=c(round(xlow()/1000000,4),round(xhigh()/1000000,4)), 
    #                                              title="position (in Mb)"))
    #   
    #   #####
    #   
    #   afig <- subplot(atop,style(abottom,showlegend=F),nrows=2,shareX = T, titleX = T, titleY = T, which_layout = 1, 
    #                   margin=.01)
    #   
    # })
    # 
    # 
    # output$CompRefGTL <- renderPlotly({
    #   
    #   # select primary_ref results
    #   #locds <- primary_ref_ds %>% filter(locvar==input$locuslst)
    #   #xhigh <- max(locds$genestop)+1000000
    #   #xlow <- max(0,min(locds$genestart)-1000000)
    #   #locchr <- unique(locds$chr)
    #  
    #   
    #   #####################
    #   
    #   
    #   # select reference panel specific results
    #   primary_ref_tbl <- primary_ref_ds %>% filter(genestart>=xlow() & genestop<=xhigh() & pheno==locpheno() & 
    #                                                  chr==locchr() 
    #                                                #& is.na(HLARegion) & is.na(MHCRegion)
    #                                                ) %>%
    #     select(genename,chr, pheno,genestart,genestop,genemid,log10pval)
    #   
    #   # locgwb <- twas_ds %>% filter(tissue=="GWB" & genestart>=xlow() & genestop<=xhigh() & chr==locchr() & 
    #   #                                pheno==locpheno() #& is.na(HLARegion) & is.na(MHCRegion)
    #   #                              ) %>%
    #   #   select(genename,chr, pheno,log10pval,genestart,genestop,genemid, SignifGene)
    #   # pthreshgwb <- -log10(0.05/(nrow(twas_ds[twas_ds$tissue=="GWB",])))
    #   
    #   locgtl <- twas_ds %>% filter(tissue=="GTL" & genestart>=xlow() & genestop<=xhigh() & chr==locchr() & 
    #                                  pheno==locpheno() #& is.na(HLARegion) & is.na(MHCRegion)
    #                                ) %>%
    #     select(genename,chr,pheno,log10pval,genestart,genestop,genemid, SignifGene)
    #   pthreshgtl <- -log10(0.05/(nrow(twas_ds[twas_ds$tissue=="GTL",])))
    #   
    #   # locmsa <- twas_ds %>% filter(tissue=="MSA" & genestart>=xlow() & genestop<=xhigh() & chr==locchr() & 
    #   #                                pheno==locpheno() #& is.na(HLARegion) & is.na(MHCRegion)
    #   #                              ) %>%
    #   #   select(genename,chr, pheno,log10pval,genestart,genestop,genemid, SignifGene)
    #   # pthreshmsa <- -log10(0.05/(nrow(twas_ds[twas_ds$tissue=="MSA",])))
    #   
    #   
    #   #####################
    #   
    #   
    #   # merge secondary reference panel data sets with DGN
    #   # primary_ref_gwb <- merge(primary_ref_tbl,locgwb, by=c("genename","chr","pheno"), all = T)
    #   # primary_ref_gwb$inboth <- ifelse(complete.cases(primary_ref_gwb$log10pval.x,primary_ref_gwb$log10pval.y),"In Both","Not in Both")
    #   
    #   primary_ref_gtl <- merge(primary_ref_tbl,locgtl, by=c("genename","chr","pheno"), all = T)
    #   primary_ref_gtl$inboth <- ifelse(complete.cases(primary_ref_gtl$log10pval.x,primary_ref_gtl$log10pval.y),"In Both","Not in Both")
    #   
    #   # primary_ref_msa <- merge(primary_ref_tbl,locmsa, by=c("genename","chr","pheno"), all = T)
    #   # primary_ref_msa$inboth <- ifelse(complete.cases(primary_ref_msa$log10pval.x,primary_ref_msa$log10pval.y),"In Both","Not in Both")
    #   
    #   
    #   #####################
    #   
    #   
    #   # set y limit values for each plot
    #   yhighprimary_ref_ <- max(primary_ref_gtl$log10pval.x, na.rm = T)+0.15*max(primary_ref_gtl$log10pval.x, na.rm = T)
    #   nudgevalprimary_ref_ <- 0.05*yhighprimary_ref_
    #   
    #   yhighgtl <- max(primary_ref_gtl$log10pval.y,pthreshgtl, na.rm = T)+0.15*max(primary_ref_gtl$log10pval.y,pthreshgtl, 
    #                                                                               na.rm = T)
    #   nudgevalgtl <- 0.05*yhighgtl
    #   
    #   # yhighgwb <- max(primary_ref_gwb$log10pval.y,pthreshgwb, na.rm = T)+0.15*max(primary_ref_gwb$log10pval.y,pthreshgwb, na.rm = T)
    #   # nudgevalgwb <- 0.05*yhighgwb
    #   # 
    #   # yhighmsa <- max(primary_ref_msa$log10pval.y,pthreshmsa, na.rm = T)+0.15*max(primary_ref_msa$log10pval.y,pthreshmsa, na.rm = T)
    #   # nudgevalmsa <- 0.05*yhighmsa
    #   
    #   
    #   myColors <- setNames( c('#000000', '#56B4E9', '#E69F00'),
    #                         c("In Both","Not in Both",NA))
    #   
    #   colScale <- scale_colour_manual(values = myColors)
    #   
    #   
    #   ####################
    #   
    #   # plot info for GTL
    #   btop <- ggplotly(ggplot(data=primary_ref_gtl, 
    #                           aes(x=round(genemid.x/1000000,4),
    #                               y=log10pval.x, color=inboth,
    #                               text = paste("Gene: ",genename,
    #                                            "<br>-log10 p-value: ",round(log10pval.x,3)))) +
    #                      geom_point(pch=15) +
    #                      geom_segment(aes(x=round(genestart.x/1000000,4),y=log10pval.x, xend=round(genestop.x/1000000,4), 
    #                                       yend=log10pval.x), size=2) +
    #                      geom_text(aes(label=genename), nudge_y = nudgevalprimary_ref_) +
    #                      geom_hline(aes(yintercept=pthresh), lty=2, color='red') +
    #                      xlim(round(xlow()/1000000,4),round(xhigh()/1000000,4)) +
    #                      ylim(0,yhighprimary_ref_) +
    #                      theme_bw() +
    #                      theme(legend.position = 'top', legend.title = element_blank()) +
    #                      annotate(geom="text",x=round(xlow()/1000000,4)+.370000, y=pthresh+nudgevalprimary_ref_, 
    #                               color='red',
    #                               label=paste0("Primary ref TWAS p-value: ", formatC(10^-(pthresh), format = "e", 
    #                                                                                  digits = 2)),size=4) +
    #                      colScale
    #   , tooltip = "text")
    #   
    #   btop <- btop %>% plotly::layout(yaxis = list(title = paste0(primary_tissue,' TWAS -log10(p)')),
    #                           xaxis = list(range=c(round(xlow()/1000000,4),round(xhigh()/1000000,4))),
    #                           annotations=list(x = 0.5 , y = 1.1, text = paste0("(b) ",primary_tissue," vs. GTL"), showarrow = F, 
    #                                            xref='paper', yref='paper',xanchor='center'),
    #                           legend = list(orientation='h', x=0, y=1),
    #                           title=list(text=paste0(loctitle(),"\n"),x=0,xanchor='left'), margin=list(t=60))
    #   
    #   #####
    #   
    #   bbottom <- ggplotly(ggplot(data=primary_ref_gtl, 
    #                              aes(x=round(genemid.y/1000000,4),
    #                                  y=-log10pval.y, color=inboth,
    #                                  text = paste("Gene: ",genename,
    #                                               "<br>log10 p-value: ",round(-log10pval.y,3)))) +
    #                         geom_hline(aes(yintercept=-pthreshgtl), lty=2, color='red') +
    #                         geom_point(pch=15) +
    #                         geom_segment(aes(x=round(genestart.y/1000000,4),y=-log10pval.y, 
    #                                          xend=round(genestop.y/1000000,4), yend=-log10pval.y), size=2) +
    #                         geom_text(aes(label=genename), nudge_y = -1.5*nudgevalgtl) +
    #                         xlim(round(xlow()/1000000,4),round(xhigh()/1000000,4)) +
    #                         ylim(-yhighgtl,0) +
    #                         theme_bw() +
    #                         theme(legend.position = 'top', legend.title = element_blank()) +
    #                         annotate(geom="text",x=round(xlow()/1000000,4)+.370000, y=-pthreshgtl-nudgevalprimary_ref_, 
    #                                  color='red',
    #                                  label=paste0("GTL TWAS p-value: ", formatC(10^-(pthreshgtl), format = "e", 
    #                                                                             digits = 2)),size=4) +
    #                         
    #                         colScale
    #                       , tooltip = "text")
    #   
    #   bbottom <- bbottom %>% plotly::layout(yaxis = list(title = 'GTL TWAS log10(p)', 
    #                                              range=c(-yhighgtl,0.1)),
    #                                 xaxis = list(range=c(round(xlow()/1000000,4),round(xhigh()/1000000,4)), 
    #                                              title="position (in Mb)"))
    #   
    #   #####
    #   
    #   bfig <- subplot(btop,style(bbottom,showlegend=F),nrows=2,shareX = T, titleX = T, titleY = T, 
    #                   which_layout = 1, margin=0.01)
    #   
    # })
    # 
    # 
    # 
    # ######## compare Ref #############
    # 
    # 
    # output$CompRefMSA <- renderPlotly({
    #   
    #   # select primary_ref results
    #   #locds <- primary_ref_ds %>% filter(locvar==input$locuslst)
    #   #xhigh <- max(locds$genestop)+1000000
    #   #xlow <- max(0,min(locds$genestart)-1000000)
    #   #locchr <- unique(locds$chr)
    # 
    #   
    #   ####
    #   
    #   
    #   # select reference panel specific results
    #   primary_ref_tbl <- primary_ref_ds %>% filter(genestart>=xlow() & genestop<=xhigh() & pheno==locpheno() & chr==locchr() 
    #                                                #& is.na(HLARegion) & is.na(MHCRegion)
    #                                                ) %>%
    #     select(genename,chr, pheno,genestart,genestop,genemid,log10pval)
    #   
    #   # locgwb <- twas_ds %>% filter(tissue=="GWB" & genestart>=xlow() & genestop<=xhigh() & chr==locchr() & 
    #   #                                pheno==locpheno() #& is.na(HLARegion) & is.na(MHCRegion)
    #   #                              ) %>%
    #   #   select(genename,chr, pheno,log10pval,genestart,genestop,genemid, SignifGene)
    #   # pthreshgwb <- -log10(0.05/(nrow(twas_ds[twas_ds$tissue=="GWB",])))
    #   # 
    #   # locgtl <- twas_ds %>% filter(tissue=="GTL" & genestart>=xlow() & genestop<=xhigh() & chr==locchr() & 
    #   #                                pheno==locpheno() #& is.na(HLARegion) & is.na(MHCRegion)
    #   #                              ) %>%
    #   #   select(genename,chr, pheno,log10pval,genestart,genestop,genemid, SignifGene)
    #   # pthreshgtl <- -log10(0.05/(nrow(twas_ds[twas_ds$tissue=="GTL",])))
    #   
    #   locmsa <- twas_ds %>% filter(tissue=="MSA" & genestart>=xlow() & genestop<=xhigh() & chr==locchr() & 
    #                                  pheno==locpheno() #& is.na(HLARegion) & is.na(MHCRegion)
    #                                ) %>%
    #     select(genename,chr, pheno,log10pval,genestart,genestop,genemid, SignifGene)
    #   pthreshmsa <- -log10(0.05/(nrow(twas_ds[twas_ds$tissue=="MSA",])))
    #   
    #   
    #   ####
    #   
    #   
    #   # merge secondary reference panel data sets with primary ref
    #   # primary_ref_gwb <- merge(primary_ref_tbl,locgwb, by=c("genename","chr","pheno"), all = T)
    #   # primary_ref_gwb$inboth <- ifelse(complete.cases(primary_ref_gwb$log10pval.x,primary_ref_gwb$log10pval.y),"In Both","Not in Both")
    #   # 
    #   # primary_ref_gtl <- merge(primary_ref_tbl,locgtl, by=c("genename","chr","pheno"), all = T)
    #   # primary_ref_gtl$inboth <- ifelse(complete.cases(primary_ref_gtl$log10pval.x,primary_ref_gtl$log10pval.y),"In Both","Not in Both")
    #   
    #   primary_ref_msa <- merge(primary_ref_tbl,locmsa, by=c("genename","chr","pheno"), all = T)
    #   primary_ref_msa$inboth <- ifelse(complete.cases(primary_ref_msa$log10pval.x,primary_ref_msa$log10pval.y),"In Both","Not in Both")
    #   
    #   
    #   ####
    #   
    #   
    #   # set y limit values for each plot
    #   yhighprimary_ref_ <- max(primary_ref_msa$log10pval.x, na.rm = T)+0.15*max(primary_ref_msa$log10pval.x, na.rm = T)
    #   nudgevalprimary_ref_ <- 0.05*yhighprimary_ref_
    #   
    #   # yhighgtl <- max(primary_ref_gtl$log10pval.y,pthreshgtl, na.rm = T)+0.15*max(primary_ref_gtl$log10pval.y,pthreshgtl, na.rm = T)
    #   # nudgevalgtl <- 0.05*yhighgtl
    #   # 
    #   # yhighgwb <- max(primary_ref_gwb$log10pval.y,pthreshgwb, na.rm = T)+0.15*max(primary_ref_gwb$log10pval.y,pthreshgwb, na.rm = T)
    #   # nudgevalgwb <- 0.05*yhighgwb
    #   
    #   yhighmsa <- max(primary_ref_msa$log10pval.y,pthreshmsa, na.rm = T)+0.15*max(primary_ref_msa$log10pval.y,pthreshmsa, 
    #                                                                               na.rm = T)
    #   nudgevalmsa <- 0.05*yhighmsa
    #   
    #   
    #   myColors <- setNames( c('#000000', '#56B4E9', '#E69F00'),
    #                         c("In Both","Not in Both",NA))
    #   
    #   colScale <- scale_colour_manual(values = myColors)
    #   
    #   
    #   ####################
    #   
    #   # plot info for MSA
    #   ctop <- ggplotly(ggplot(data=primary_ref_msa, 
    #                           aes(x=round(genemid.x/1000000,4),
    #                               y=log10pval.x, color=inboth,
    #                               text = paste("Gene: ",genename,
    #                                            "<br>-log10 p-value: ",round(log10pval.x,3)))) +
    #                      geom_point(pch=15) +
    #                      geom_segment(aes(x=round(genestart.x/1000000,4),y=log10pval.x, xend=round(genestop.x/1000000,4), 
    #                                       yend=log10pval.x), size=2) +
    #                      geom_text(aes(label=genename), nudge_y = nudgevalprimary_ref_) +
    #                      geom_hline(aes(yintercept=pthresh), lty=2, color='red') +
    #                      xlim(round(xlow()/1000000,4),round(xhigh()/1000000,4)) +
    #                      ylim(0,yhighprimary_ref_) +
    #                      theme_bw() +
    #                      theme(legend.position = 'top', legend.title = element_blank()) +
    #                      annotate(geom="text",x=round(xlow()/1000000,4)+.370000, y=pthresh+nudgevalprimary_ref_, 
    #                               color='red',
    #                               label=paste0("Primary ref TWAS p-value: ", formatC(10^-(pthresh), format = "e",
    #                                                                                  digits = 2)),size=4) +
    #                      colScale
    #   , tooltip = "text")
    #   
    #   ctop <- ctop %>% plotly::layout(yaxis = list(title = paste0(primary_tissue, 'TWAS -log10(p)')),
    #                           xaxis = list(range=c(round(xlow()/1000000,4),round(xhigh()/1000000,4))),
    #                           annotations=list(x = 0.5 , y = 1.1, text = paste0("(c) ",primary_tissue," vs. MSA"), showarrow = F, 
    #                                            xref='paper', yref='paper',xanchor='center'),
    #                           legend = list(orientation='h', x=0, y=1),
    #                           title=list(text=paste0(loctitle(),"\n"),x=0,xanchor='left'),margin=list(t=60))
    #   
    #   #####
    #   
    #   cbottom <- ggplotly(ggplot(data=primary_ref_msa, 
    #                              aes(x=round(genemid.y/1000000,4),
    #                                  y=-log10pval.y, color=inboth,
    #                                  text = paste("Gene: ",genename,
    #                                               "<br>log10 p-value: ",round(-log10pval.y,3)))) +
    #                         geom_hline(aes(yintercept=-pthreshmsa), lty=2, color='red') +
    #                         geom_point(pch=15) +
    #                         geom_segment(aes(x=round(genestart.y/1000000,4),y=-log10pval.y, 
    #                                          xend=round(genestop.y/1000000,4), yend=-log10pval.y), size=2) +
    #                         geom_text(aes(label=genename), nudge_y = -1.5*nudgevalmsa) +
    #                         xlim(round(xlow()/1000000,4),round(xhigh()/1000000,4)) +
    #                         theme_bw() +
    #                         ylim(-yhighmsa-2,0) +
    #                         annotate(geom="text",x=round(xlow()/1000000,4)+.370000, y=-pthreshmsa-nudgevalprimary_ref_, 
    #                                  color='red',
    #                                  label=paste0("MSA TWAS p-value: ", formatC(10^-(pthreshmsa), format = "e", 
    #                                                                             digits = 2)),size=4) +
    #                         
    #                         colScale
    #   , tooltip = "text")
    #   
    #   cbottom <- cbottom %>% plotly::layout(yaxis = list(title = 'MSA TWAS log10(p)', 
    #                                              range=c(-yhighmsa,0.1)),
    #                                 xaxis = list(range=c(round(xlow()/1000000,4),round(xhigh()/1000000,4)), 
    #                                              title="position (in Mb)"))
    #   
    #   #####
    #   
    #   cfig <- subplot(ctop,style(cbottom,showlegend=F),nrows=2,shareX = T, titleX = T, titleY = T, which_layout = 1,
    #                   margin = 0.01)
    #   
    #   
    #   #####################
    #   
    # })
    # 
    # 
    ############# End commenting out here for testing
    
    
    ############# Table TWAS result ###############################################
    
    if (multiple_tissues == FALSE){
      output$primarytissuetbl <- DT::renderDataTable({

        # select all genes at locus
        primary_ref_tbl <- primary_ref_ds %>% filter(genestart>=xlow() & genestop<=xhigh() & pheno==locpheno() & chr==locchr()) %>%
          select(-locus,-locstart,-locstop,-genestartMB,-genestopMB,-genemid,-genemidMB,-locvar,-log10pvalmeta,-log10pval)

        # print the data table
        DT::datatable(primary_ref_tbl, rownames=F)

      })
    } else{ 
    
    # TWAS results table with separate tabs for each reference panel
    output$DGNtbl <- DT::renderDataTable({
      
      #locds <- primary_ref_ds %>% filter(locvar==input$locuslst)
      #xhigh <- max(locds$genestop)+1000000
      #xlow <- max(0,min(locds$genestart)-1000000)
      #locpheno <- unique(locds$phenoname) #locus phenotype
      #locchr<- unique(locds$chr) #locus chromosome
      
      # select all genes at locus
      primary_ref_tbl <- primary_ref_ds %>% filter(genestart>=xlow() & genestop<=xhigh() & pheno==locpheno() & chr==locchr()) %>%
        #select(genename,chr, genestart, genestop, phenoname,phenocat,se_beta_, p, 
        #       weightsnp, cohortsnp, SignifGene,HLARegion,MHCRegion,SingleSNP)
        select(-locus,-locstart,-locstop,-genestartMB,-genestopMB,-genemid,-genemidMB,-locvar,-log10pvalmeta,-log10pval)
      #colnames(primary_ref_tbl) <- c("gene","chr","gene start","gene stop","trait","trait category","beta","p-value",
      #                      "# model SNPs","# SNPs used","TWAS significant",
      #                      "HLA gene","MHC region","1 SNP model")
      
      # print data table
      ## commenting this out and testing another
      # datatable(primary_ref_tbl, rownames=F, options = list(
      #   order = list(list(7, 'asc')))) %>% formatStyle(
      #     "TWAS significant", target = "row",
      #     backgroundColor = styleEqual(c(1), c('lightgreen'))
      #   ) %>% formatStyle(
      #     "HLA gene", target = "row",
      #     backgroundColor = styleEqual(c(1), c('red'))
      #   ) %>% formatStyle(
      #     "MHC region", target = "row",
      #     backgroundColor = styleEqual(c(1), c('red'))
      #   ) %>% formatStyle(
      #     "1 SNP model", target = "row",
      #     backgroundColor = styleEqual(c(1), c('red'))
      #   )
      
      #primary_ref_tblfin <- primary_ref_tbl[,c(input$show_vars)]
      DT::datatable(primary_ref_tbl, rownames=F)
      #primary_ref_tblplt <- primary_ref_tbl[, input$show_vars]
      #datatable(primary_ref_tblfin, rownames=F)
      #df <- as.data.frame(primary_ref_tbl)
      #DT::datatable(df)
      #DT::datatable(primary_ref_tbl())
      #DT::datatable(primary_ref_tbl[,input$show_vars, drop=FALSE])
    })
    
    
    #####################
    
    
    output$GTLtbl <- DT::renderDataTable({
      
      # locds <- primary_ref_ds %>% filter(locvar==input$locuslst)
      # xhigh <- max(locds$genestop)+1000000
      # xlow <- max(0,min(locds$genestart)-1000000)
      # locpheno <- unique(locds$phenoname) #locus phenotype
      # locchr<- unique(locds$chr) #locus chromosome
      
      # select all genes at locus
      gtltbl <- twas_ds %>% filter(.data$tissue=='GTL',genestart>=xlow() & genestop <=xhigh() & pheno==locpheno() & chr==locchr()) %>%
        #select(genename,chr, genestart, genestop, phenoname,phenocat,se_beta_, p, 
        #       weightsnp,cohortsnp, SignifGene,HLARegion,MHCRegion, SingleSNP)
        select(-locus,-locstart,-locstop,-genestartMB,-genestopMB,-genemid,-genemidMB,-log10pval)
      
      #colnames(gtltbl) <- c("gene","chr","gene start","gene stop","trait","trait category","beta","p-value",
      #                      "# model SNPs","#SNPs used","TWAS significant",
      #                      "HLA gene","MHC region","1 SNP model")
      
      # print data table
      # datatable(gtltbl, rownames=F, options = list(
      #   order = list(list(7, 'asc')))) %>% formatStyle(
      #     "TWAS significant", target = "row",
      #     backgroundColor = styleEqual(c(1), c('lightgreen'))
      #   ) %>% formatStyle(
      #     "HLA gene", target = "row",
      #     backgroundColor = styleEqual(c(1), c('red'))
      #   ) %>% formatStyle(
      #     "MHC region", target = "row",
      #     backgroundColor = styleEqual(c(1), c('red'))
      #   ) %>% formatStyle(
      #     "1 SNP model", target = "row",
      #     backgroundColor = styleEqual(c(1), c('red'))
      #   )
      DT::datatable(gtltbl, rownames=F)
    })
    
    
    #####################
    
    
    output$GWBtbl <- DT::renderDataTable({
      
      # locds <- primary_ref_ds %>% filter(locvar==input$locuslst)
      # xhigh <- max(locds$genestop)+1000000
      # xlow <- max(0,min(locds$genestart)-1000000)
      # locpheno <- unique(locds$phenoname) #locus phenotype
      # locchr<- unique(locds$chr) #locus chromosome
      
      # select all genes at locus
      gwbtbl <- twas_ds %>% filter(.data$tissue=='GWB',genestart>=xlow() & genestop<=xhigh() & pheno==locpheno() & chr==locchr()) %>%
        # select(genename,chr, genestart,genestop,phenoname,phenocat,se_beta_, p,weightsnp,cohortsnp, 
        #        SignifGene,HLARegion,MHCRegion,SingleSNP)
        select(-locus,-locstart,-locstop,-genestartMB,-genestopMB,-genemid,-genemidMB,-log10pval)
      # colnames(gwbtbl) <- c("gene","chr","gene start","gene stop","trait","trait category","beta","p-value",
      #                       "# model SNPs","# SNPs used","TWAS significant",
      #                       "HLA gene","MHC region","1 SNP model")
      
      # print data table
      # datatable(gwbtbl, rownames=F, options = list(
      #   order = list(list(7, 'asc')))) %>% formatStyle(
      #     "TWAS significant", target = "row",
      #     backgroundColor = styleEqual(c(1), c('lightgreen'))
      #   ) %>% formatStyle(
      #     "HLA gene", target = "row",
      #     backgroundColor = styleEqual(c(1), c('red'))
      #   ) %>% formatStyle(
      #     "MHC region", target = "row",
      #     backgroundColor = styleEqual(c(1), c('red'))
      #   ) %>% formatStyle(
      #     "1 SNP model", target = "row",
      #     backgroundColor = styleEqual(c(1), c('red'))
      #   )
      
      DT::datatable(gwbtbl, rownames=F) 
    })
    
    
    #####################
    
    
    output$MSAtbl <- DT::renderDataTable({
      
      # locds <- primary_ref_ds %>% filter(locvar==input$locuslst)
      # xhigh <- max(locds$genestop)+1000000
      # xlow <- max(0,min(locds$genestart)-1000000)
      # locpheno <- unique(locds$phenoname) #locus phenotype
      # locchr<- unique(locds$chr) #locus chromosome
      
      # select all genes at locus
      msatbl <- twas_ds %>% filter(.data$tissue=='MSA',genestart>=xlow() & genestop<=xhigh() & pheno==locpheno() & chr==locchr()) %>%
        # select(genename,chr, genestart,genestop,phenoname,phenocat,se_beta_, p,weightsnp,cohortsnp, 
        #        SignifGene,HLARegion,MHCRegion,SingleSNP)
        select(-locus,-locstart,-locstop,-genestartMB,-genestopMB,-genemid,-genemidMB,-log10pval)
      # colnames(msatbl) <- c("gene","chr","gene start","gene stop","trait","trait category","beta","p-value",
      #                       "# model SNPs","#SNPs used","TWAS significant",
      #                       "HLA gene","MHC region","1 SNP model")
      
      # print data table
      # datatable(msatbl, rownames=F, options = list(
      #   order = list(list(7, 'asc')))) %>% formatStyle(
      #     "TWAS significant", target = "row",
      #     backgroundColor = styleEqual(c(1), c('lightgreen'))
      #   ) %>% formatStyle(
      #     "HLA gene", target = "row",
      #     backgroundColor = styleEqual(c(1), c('red'))
      #   ) %>% formatStyle(
      #     "MHC region", target = "row",
      #     backgroundColor = styleEqual(c(1), c('red'))
      #   )  %>% formatStyle(
      #     "1 SNP model", target = "row",
      #     backgroundColor = styleEqual(c(1), c('red'))
      #   )
      DT::datatable(msatbl, rownames=F)
    })
    
    }
    ############################################################
    
    
    # TWAS other traits results table with separate tabs for each reference panel
    # output$primary_ref_tbltrt <- DT::renderDataTable({
    #   
    #   # select primary_ref_ all genes from other trait categories at locus
    #   primary_ref_tbl <- filter(primary_ref_ds, genestart>=xlow() & genestop<=xhigh() & phenoname!=locpheno() 
    #                    & chr==locchr() & phenocat==locphcat()) %>%
    #     select(genename,chr, genestart,genestop,phenoname,phenocat,se_beta_, p,weightsnp,cohortsnp, 
    #            SignifGene,HLARegion,MHCRegion,SingleSNP)
    #   colnames(primary_ref_tbl) <- c("gene","chr","gene start","gene stop","trait","trait category","beta","p-value",
    #                         "# model SNPs","#SNPs used","TWAS significant",
    #                         "HLA gene","MHC region","1 SNP model")
    #   
    #   # print data table
    #   datatable(primary_ref_tbl, rownames=F, options = list(
    #     order = list(list(7, 'asc')))) %>% formatStyle(
    #       "TWAS significant", target = "row",
    #       backgroundColor = styleEqual(c(1), c('lightgreen'))
    #     ) %>% formatStyle(
    #       "HLA gene", target = "row",
    #       backgroundColor = styleEqual(c(1), c('red'))
    #     ) %>% formatStyle(
    #       "MHC region", target = "row",
    #       backgroundColor = styleEqual(c(1), c('red'))
    #     ) %>% formatStyle(
    #       "1 SNP model", target = "row",
    #       backgroundColor = styleEqual(c(1), c('red'))
    #     )
    # })
    # 
    # 
    # #####################
    # 
    # 
    # output$GTLtbltrt <- DT::renderDataTable({
    #   
    # 
    #   # select all genes at locus
    #   gtltbl <- twas_ds %>% filter(tissue=='GTL',genestart>=xlow() & genestop<=xhigh() 
    #                                 & phenoname!=locpheno() & chr==locchr() & phenocat==locphcat()) %>%
    #     select(genename,chr, genestart,genestop,phenoname,phenocat,se_beta_, p,weightsnp,cohortsnp, 
    #            SignifGene,HLARegion,MHCRegion,SingleSNP)
    #   colnames(gtltbl) <- c("gene","chr","gene start","gene stop","trait","trait category","beta","p-value",
    #                         "# model SNPs","#SNPs used","TWAS significant",
    #                         "HLA gene","MHC region","1 SNP model")
    #   
    #   # print data table
    #   datatable(gtltbl, rownames=F, options = list(
    #     order = list(list(7, 'asc')))) %>% formatStyle(
    #       "TWAS significant", target = "row",
    #       backgroundColor = styleEqual(c(1), c('lightgreen'))
    #     ) %>% formatStyle(
    #       "HLA gene", target = "row",
    #       backgroundColor = styleEqual(c(1), c('red'))
    #     ) %>% formatStyle(
    #       "MHC region", target = "row",
    #       backgroundColor = styleEqual(c(1), c('red'))
    #     ) %>% formatStyle(
    #       "1 SNP model", target = "row",
    #       backgroundColor = styleEqual(c(1), c('red'))
    #     )
    # })
    # 
    # 
    # #####################
    # 
    # 
    # output$GWBtbltrt <- DT::renderDataTable({
    #   
    #   # select all genes at locus
    #   gwbtbl <- twas_ds %>% filter(tissue=='GWB',genestart>=xlow() & genestop<=xhigh()
    #                                 & phenoname!=locpheno() & chr==locchr() & phenocat==locphcat()) %>%
    #     select(genename,chr, genestart,genestop,phenoname,phenocat,se_beta_, p,weightsnp,cohortsnp, 
    #            SignifGene,HLARegion,MHCRegion,SingleSNP)
    #   colnames(gwbtbl) <- c("gene","chr","gene start","gene stop","trait","trait category","beta","p-value",
    #                         "# model SNPs","#SNPs used","TWAS significant",
    #                         "HLA gene","MHC region","1 SNP model")
    #   
    #   # print data table
    #   datatable(gwbtbl, rownames=F, options = list(
    #     order = list(list(7, 'asc')))) %>% formatStyle(
    #       "TWAS significant", target = "row",
    #       backgroundColor = styleEqual(c(1), c('lightgreen'))
    #     ) %>% formatStyle(
    #       "HLA gene", target = "row",
    #       backgroundColor = styleEqual(c(1), c('red'))
    #     ) %>% formatStyle(
    #       "MHC region", target = "row",
    #       backgroundColor = styleEqual(c(1), c('red'))
    #     ) %>% formatStyle(
    #       "1 SNP model", target = "row",
    #       backgroundColor = styleEqual(c(1), c('red'))
    #     )
    # })
    # 
    # 
    # #####################
    # 
    # 
    # output$MSAtbltrt <- DT::renderDataTable({
    #   
    #   # select all genes at locus
    #   msatbl <- twas_ds %>% filter(tissue=='MSA',genestart>=xlow() & genestop<=xhigh() 
    #                                 & phenoname!=locpheno() & chr==locchr() & phenocat==locphcat()) %>%
    #     select(genename,chr, genestart,genestop,phenoname,phenocat,se_beta_, p,weightsnp,cohortsnp, 
    #            SignifGene,HLARegion,MHCRegion,SingleSNP)
    #   colnames(msatbl) <- c("gene","chr","gene start","gene stop","trait","trait category","beta","p-value",
    #                         "# model SNPs","#SNPs used","TWAS significant",
    #                         "HLA gene","MHC region","1 SNP model")
    #   
    #   # print data table
    #   datatable(msatbl, rownames=F, options = list(
    #     order = list(list(7, 'asc')))) %>% formatStyle(
    #       "TWAS significant", target = "row",
    #       backgroundColor = styleEqual(c(1), c('lightgreen'))
    #     ) %>% formatStyle(
    #       "HLA gene", target = "row",
    #       backgroundColor = styleEqual(c(1), c('red'))
    #     ) %>% formatStyle(
    #       "MHC region", target = "row",
    #       backgroundColor = styleEqual(c(1), c('red'))
    #     ) %>% formatStyle(
    #       "1 SNP model", target = "row",
    #       backgroundColor = styleEqual(c(1), c('red'))
    #     )
    # })
    # 
    
    ############# Table known variants ###############################################
    
    # Table of known variants
    output$KnownSNPtbl <- DT::renderDataTable({
      #locds <- primary_ref_ds %>% filter(locvar==input$locuslst)
      #xhigh <- max(locds$genestop)+1000000
      #xlow <- max(0,min(locds$genestart)-1000000)
      locpheno<-unique(locds()$pheno) #locus phenotype category
      locph<-unique(locds()$pheno) #locus phenotype
      #locchr<- unique(locds$chr) #locus chromosome
      

      ds=GWAS_sentinel
      
      # select all known snps at locus
      phenotbl <- ds %>% filter(position>=xlow() & xhigh()>=position & chr==locchr()) #%>%
        #select(X1,X2,X3,X4,X5,X6,X8,X9,X10,X11,X12,X13,X14,X15,X21)
      #colnames(phenotbl) <- c("RSID", "Reference", "Ancestry", "Trait", "Chr", "Pos_b37", "Gene", "Effect Allele",
      #                        "Other Allele", "EAF", "Beta", "SE", "P", "N","In Kaiser")
      
      if (nrow(phenotbl)==0){
        #df <- data.frame(a="No reported GWAS variants at this locus")
        datatable(phenotbl)#,
                  #options=list(columnDefs = list(list(visible=FALSE, targets=c(15,16,17,18,19,20,21,22,23)))))
      } else {
        
        # select set of TWAS insignificant genes
        primary_ref_tbl <- primary_ref_ds %>% filter(genestart>=xlow() & genestop<=xhigh() & pheno==locph & chr==locchr() &
                                     SignifGene==0  #& is.na(HLARegion) & is.na(MHCRegion)
                                     ) %>% select(genename)
        
        # indicator for significant TWAS gene pattern 
        phenotbl$signiftwas <- ifelse(grepl(paste(unique(locds()$genename), collapse="|"),phenotbl$genename),1,0)
        
        # indicator for non-significant TWAS gene pattern
        phenotbl$nonsigniftwas <- ifelse(grepl(paste(unique(primary_ref_tbl$genename),collapse="|"),phenotbl$genename),1,0)
        
        # indicator for gene not predicted in primary tissue
        phenotbl$notpred <- ifelse(grepl(paste(ref_panel_genes$genename[ref_panel_genes$genestatus==0],collapse = "|"),phenotbl$genename),1,0)
        phenotbl$notprimary_ref <- ifelse(grepl(paste(ref_panel_genes$genename[ref_panel_genes$genestatus==1],collapse="|"),phenotbl$genename),1,0)
        
        phenotbl$chrposall <- paste0(phenotbl$chr,":",phenotbl$position,":",phenotbl$`effectallele`,":",phenotbl$`otherallele`)
        phenotbl$chrposall2 <- paste0(phenotbl$chr,":",phenotbl$position,":",phenotbl$`otherallele`,":",phenotbl$`effectallele`)
        
        # indicator for snps included in GWAS results
        locgwas <- cohort_gwas_knownfin %>% filter(Locus==locnum()) %>% select(.data$SNP)
        phenotbl$X19 <- ifelse(grepl(paste(unique(locgwas$SNP),collapse="|"),phenotbl$chrposall),1,0)
        phenotbl$X20 <- ifelse(grepl(paste(unique(locgwas$SNP),collapse="|"),phenotbl$chrposall2),1,0)
        
        # find target numbers for variable numbers we don't want to print
        phenotblcols <- ncol(phenotbl)
        phenotbltarg <- phenotblcols-7
        targetlst <- seq(to=phenotblcols,from=phenotbltarg,by=1)
        targetlstfin <- targetlst[targetlst != (phenotblcols-3)]
        
        # print data table
        datatable(phenotbl , #rownames=FALSE,
                  options=list(columnDefs = list(list(visible=FALSE, 
                  #                                    targets=c(15,16,17,18,19,20,21,22,23))))
                  #                                   targets=c("X19","X20"))))
                   targets=c(targetlstfin))))
                  ) %>%
          formatStyle("genename", "signiftwas", backgroundColor = styleEqual(
            c(1), c("lightgreen"))
          ) %>% 
          formatStyle("incohort", target = "row",
                      backgroundColor = styleEqual(c(0), c('red'))
          ) %>%
          formatStyle("genename", "nonsigniftwas", backgroundColor = styleEqual(
            c(1), c("yellow"))
          ) %>%
          formatStyle("chrposall", "X19", backgroundColor = styleEqual(
            c(1), c("lightgreen"))
          ) %>%
          formatStyle("chrposall", "X20", backgroundColor = styleEqual(
            c(1), c("lightgreen"))
          ) %>%
          formatStyle("genename", "notpred", backgroundColor = styleEqual(
            c(1), c("orange"))
          ) %>%
          formatStyle("genename", "notprimary_ref", backgroundColor = styleEqual(
            c(1), c("red"))
          )
      }
    })
    
    
    ############# Table GWAS variants ###############################################
    
    # Table of known variants
    output$GWASvars <- DT::renderDataTable({
      locgwas <- cohort_gwas_knownfin %>% filter(Locus==locnum())
      
      # find target numbers for variable numbers we don't want to print
      snpcolnum <- match("SNP",names(cohort_gwas_knownfin)) #col num of SNP variable
      cohort_gwas_cols <- ncol(cohort_gwas_knownfin) #total num of cols
      targetlst1 <- seq(to=snpcolnum+4,from=snpcolnum+1,by=1)
      targetlst2 <- seq(to=cohort_gwas_cols, from=cohort_gwas_cols-1,by=1)
      targetlstfin <- c(targetlst1,targetlst2)
      
      # print data table
      datatable(locgwas, #rownames=F, 
                options=list(columnDefs = list(list(visible=FALSE, 
                                                    targets=targetlstfin)))
                )
    })
    
    
    
    ############# Weights table ###############################################
    
    
    # Primary reference panel weights table
    output$weighttbl <- DT::renderDataTable({
      
      # select significant and known genes to plot
      primary_ref_tblplt <- primary_ref_tbl() %>% filter(SignifGene==1 | kngene=="Reported in GWAS")
      uniqgenes <- primary_ref_tblplt %>% select(genename)
      indexgene <- primary_ref_tblplt$genename[primary_ref_tblplt$p == min(primary_ref_tblplt$p)]
      
      
      # extract set of genes from correlation matrix
      Mindex <- data.table(M[uniqgenes$genename,indexgene])
      colnames(Mindex) <- c("corr")
      Mindex$genename <- uniqgenes$genename
      
      
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
      corplt <- merge(primary_ref_tblplt,Mindex,by.x="genename",by.y="genename")
      
      # get the primary_ref weights for the specific genes at locus
      primary_ref_wtloc <- weight_ds %>% filter(genename %in% corplt$genename)
      
      # subset the GWAS variants for specific chr, and phenotype
      gwasallfinloc <- gwasallfin %>% filter(LOCUS==locnum()) %>%
        select(pos,all1,all2,poslog10p)
      gwasallfinloc$posnum <- as.numeric(gwasallfinloc$pos)
      
      # get the matching GWAS variants for the weights
      primary_ref_wtgwasloc1 <- merge(primary_ref_wtloc,gwasallfinloc,by.x=c('position','ref_allele','eff_allele'),
                                      by.y=c('posnum','all1','all2'))
      primary_ref_wtgwasloc2 <- merge(primary_ref_wtloc,gwasallfinloc,by.x=c('position','ref_allele','eff_allele'),
                                      by.y=c('posnum','all2','all1'))
      primary_ref_wtgwasloc <- rbind(primary_ref_wtgwasloc1,primary_ref_wtgwasloc2)
      
      # merge GWAS results with TWAS info
      TWASloc <- corplt %>% select(genename,genestartMB,genemid,genestopMB,log10pval,corgroup,corr)
      TWASloc$genemidMB <- round(TWASloc$genemid/1000000,4)
      primary_ref_wtgwaslocfin <- merge(primary_ref_wtgwasloc,TWASloc, by.x='genename',by.y='genename',all.x=T)
      
      # set y limit
      yhigh <- max(corplt$log10pval)+0.25*max(corplt$log10pval)
      ylow <- min(primary_ref_wtgwasloc$poslog10p,log10(5*10^(-8)))+0.15*min(primary_ref_wtgwasloc$poslog10p,log10(5*10^(-8)))
      nudgeval <- 0.05*max(abs(yhigh))
      
      # set color scale values
      colval <- c("[0 - 0.2]"="#CCCCCC","(0.2 - 0.4]"="#333FFF",
                  "(0.4 - 0.6]"="#00FF00","(0.6 - 0.8]"="#FF9900","(0.8 - 1]"="#FF0000")
      breakval <- c("[0 - 0.2]","(0.2 - 0.4]","(0.4 - 0.6]","(0.6 - 0.8]","(0.8 - 1]")
      
      
      # select the LD for the given locus
      LDdsloc <- as.data.table(filter(LDds,locus==locnum()))
      
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
      primary_ref_wtgwaslocfin2 <- merge(primary_ref_wtgwaslocfin,LDdsloc, by.x=c('position'),by.y=c('posb'),all.x=T)
      
      primary_ref_wtgwaslocfin2_comp <- primary_ref_wtgwaslocfin2[complete.cases(primary_ref_wtgwaslocfin2), ] %>%
        select(rsid,weight,corrab,genename,corr,ref_allele,eff_allele,chr,position)
      colnames(primary_ref_wtgwaslocfin2_comp) <-
        c("rsid","weight","LD","gene","corr","ref allele","eff allele","chr","position")
    
      # print data table
      datatable(primary_ref_wtgwaslocfin2_comp
                )
      
    })
    
  if (add_UI == ""){
    newdata <- reactive({
      req(add_server!="")
      
      add_server
    })
  } else{
    add_server
  }

  }
  


######### Render app ##################################################################


  
  # render the app
  shiny::shinyApp(ui=ui, server=server)

}
