

##################
# load packages and data
library("shiny")
library("readxl")
library("DT")

df <- read_excel("data/OXTR_table.xlsx")
df <- df[!is.na(df$Effect_CTQ), ]
df$`chr3 position` <- as.character(df$`chr3 position`)


##################
# layout


ui <- fluidPage(
  
  tabsetPanel(
    
    ############
    # Tutorial page
    tabPanel("Tutorial",
             
             h3("Background"),
             
             p("DNA methylation of the OXTR gene is a popular research target 
               in relation to stress, early adversity, and psychosocial functioning."),
             
             p("This shiny app was designed to guide researchers which CpG site to investigate in their research.
               It allows selecting CpG sites based on different filters, e.g. their 
               location on the gene, statistical features, and associations with three external outcomes
               (Childhood Trauma Questionnaire [CTQ], dichotomous trauma group membership, OXTR mRNA expression).
               "),
             
             p("It is based on an effective sample of 110 women (67 with mRNA expression data), covering 183 CpG sites (ADD PREPRINT)"),
             
             
             h3("Getting Started"),
             
             p("Switch to the 'Use App' tab and select desired filters to choose CpG sites. The app will
             display a filtered table which is both search- and downloadable. On the bottom of the page you 
             can see how many CpG sites are selected. Clicking on the column names enables sorting the table
             for the values in this column. More information on individual filters and variables
             in the output table can be found below"),
             
             
             h3("Referencing this app"),
             
             p("When using this app for a publication, please consider citing [...]"),
             
             
             h3("Details on Filters"),
             
             p("CpG sites can be filtered in respect to three external outcomes (CTQ, Trauma group, mRNA expression).
               All three outcomes come with the same filters. The p-value filter selects the maximum allowed value,
               while the effect size/Bayes factor filter selects the minimum allowed value. Effect size filters are applied
               to the absolute effect size, i.e. without respect to positive/negative sign"),
             p("The Signal-to-Noise filter is the ratio between empirical standard deviation and the sensitivity of the assay"),
             
             p("Cluster Assignment was based on a cluster analysis reported in [...]"),
             
             
             h3("Variable Notes"),
             
             p("Effect_XXX: Association between methylation and outcome. Correlation for CTQ/mRNA-expression, cohen's d for trauma group"),
             
             p("Effect_XXX_abs: absolute value of effect sizes"),
             
             p("P_XXX: p-value of bivariate association"),
             
             p("FDR_XXX: FDR corrected p-value at FDR < .05"),
             
             p("BF_XXX: Bayes factor for bivariate associations using the Bayes factor package and default priors"),
             
             p("Diff_group: Difference in methylation between groups in original unstandardized metric (percent methylated)"),
             
             p("pInsuffVar: P-value from the chi²-test of empirical variance against squared assay sensitivity")
             
             
             
             ),
    

    
    ############
    # App page
    tabPanel("Use App",
      
      #########
      # filters
      fluidRow(
        
        column(2, 
               wellPanel(
                 
                 h3("Gene section", align = "center"),
                 
                 div(checkboxGroupInput("segment",
                                        label = "Section",
                                        choices = c("Exon 1" = "Exon 1",
                                                    "Intron 1" = "Intron 1",
                                                    "Exon 2" = "Exon 2",
                                                    "Intron 2" = "Intron 2",
                                                    "Exon 3" = "Exon 3"),
                                        selected = c("Exon 1", "Intron 1", "Exon 2", "Intron 2", "Exon 3")),
                     style="font-size:90%;"),
                 
                 div(checkboxInput("MT2",
                                   label = "MT2 only",
                                   value = F),
                     style="font-size:90%;"),
                 
                 div(checkboxInput("Exon3_transl",
                                   label = "Protein-coding region only",
                                   value = F),
                     style="font-size:90%;"),
                 
                 div(checkboxInput("Infin450k",
                                   label = "Infinium 450k only",
                                   value = F),
                     style="font-size:90%;"),
                 
                 div(checkboxInput("InfinEPIC",
                                   label = "Infinium EPIC only",
                                   value = F),
                     style="font-size:90%;")
                 
                 
               )
        ),
        
        #####
        # Association with CTQ
        column(2,
          wellPanel(
            
            h3("CTQ", align = "center"),
          
            div(sliderInput("Effect_CTQ", 
                        label = "Correlation",
                        min = 0, max = 1, value = 0, step = 0.01,),
                style="font-size:90%;"),
            
            div(sliderInput("P_CTQ", 
                        label = "P-value",
                        min = 0, max = 1, value = 1, step = 0.001),
                style="font-size:90%;"),
            
            div(sliderInput("FDR_CTQ", 
                            label = "P-value (FDR corrected)",
                            min = 0, max = 1, value = 1, step = 0.001),
                style="font-size:90%;"),
            
            div(sliderInput("BF_CTQ", 
                            label = "Bayes Factor",
                            min = 0, max = 10, value = 0, step = 0.1),
                style="font-size:90%;")
          )
        ),
        
        #####
        # Association with group
        column(2,
               wellPanel(
                 
                 h3("Trauma group", align = "center"),
                 
                 div(sliderInput("Effect_group", 
                                 label = "Cohen's d",
                                 min = 0, max = 1, value = 0, step = 0.01,),
                     style="font-size:90%;"),
                 
                 div(sliderInput("P_group", 
                                 label = "P-value",
                                 min = 0, max = 1, value = 1, step = 0.001),
                     style="font-size:90%;"),
                 
                 div(sliderInput("FDR_group", 
                                 label = "P-value (FDR corrected)",
                                 min = 0, max = 1, value = 1, step = 0.001),
                     style="font-size:90%;"),
                 
                 div(sliderInput("BF_group", 
                                 label = "Bayes Factor",
                                 min = 0, max = 10, value = 0, step = 0.1),
                     style="font-size:90%;")
               )
        ),
        
        #####
        # Association with mRNA expression
        column(2,
               wellPanel(
                 
                 h3("mRNA expression", align = "center"),
                 
                 div(sliderInput("Effect_mRNA_expr", 
                                 label = "Correlation",
                                 min = 0, max = 1, value = 0, step = 0.01,),
                     style="font-size:90%;"),
                 
                 div(sliderInput("P_mRNA_expr", 
                                 label = "P-value",
                                 min = 0, max = 1, value = 1, step = 0.001),
                     style="font-size:90%;"),
                 
                 div(sliderInput("FDR_mRNA_expr", 
                                 label = "P-value (FDR corrected)",
                                 min = 0, max = 1, value = 1, step = 0.001),
                     style="font-size:90%;"),
                 
                 div(sliderInput("BF_mRNA_expr", 
                                 label = "Bayes Factor",
                                 min = 0, max = 10, value = 0, step = 0.1),
                     style="font-size:90%;")
               )
        ),
        
        
        #####
        # other options
        column(2,
               wellPanel(
                 
                 h3("Other options", align = "center"),
                 
                 div(checkboxGroupInput("clusterAssignment",
                                      label = "Cluster Assignment",
                                      choices = c("No Cluster" = 0,
                                                  "Cluster 1" = 1,
                                                  "Cluster 2" = 2),
                                      selected = c(0,1,2)),
                                    style="font-size:90%;"),
                 
                 div(sliderInput("percOutlier", 
                                 label = "Max. Percent Outliers",
                                 min = 0, max = 100, value = 100, step = 0.1),
                     style="font-size:90%;"),
                 
                 div(sliderInput("Signal2Noise", 
                                 label = "Signal-to-Noise",
                                 min = 0, max = 10, value = 0, step = 0.1),
                     style="font-size:90%;"),
                 
                 div(sliderInput("pInsuffVar", 
                                 label = "P-Value Variance-Test",
                                 min = 0, max = 1, value = 1, step = 0.001),
                     style="font-size:90%;")
                     )
               
                  )
        
      ),
      
      
      #########
      # table output
      dataTableOutput("view"),
      
    )
  )
)




server <- function(input, output){
  
  
  output$view <- renderDataTable({
    df[
      df$Effect_CTQ_abs >= input$Effect_CTQ &
      df$P_CTQ <= input$P_CTQ &
      df$FDR_CTQ <= input$FDR_CTQ &
      df$BF_CTQ >= input$BF_CTQ &
         
      df$Effect_group_abs >= input$Effect_group &
      df$P_group <= input$P_group &
      df$FDR_group <= input$FDR_group &
      df$BF_group >= input$BF_group & 
      
      df$Effect_mRNA_expr_abs >= input$Effect_mRNA_expr &
      df$P_mRNA_expr <= input$P_mRNA_expr &
      df$FDR_mRNA_expr <= input$FDR_mRNA_expr &
      df$BF_mRNA_expr >= input$BF_mRNA_expr &
      
      df$clusterAssignment %in% input$clusterAssignment &
      df$percOutlier <= input$percOutlier &
      df$Signal_to_Noise >= input$Signal2Noise &
      df$pInsuffVar <= input$pInsuffVar &
      
      df$Segment %in% input$segment & 
        
        ifelse(input$MT2 == T & df$MT2 == T, 
               T, 
               ifelse(
                 input$Exon3_transl == T & df$Exon3_Translational == T,
                 T,
                 ifelse(
                   input$MT2 == F & input$Exon3_transl == F,
                   T,
                   F
                 )
               )) &
        
        ifelse(input$Infin450k == T & !is.na(df$`Infinium 450K`), 
               T, 
               ifelse(
                 input$InfinEPIC == T & !is.na(df$`Infinium EPIC`),
                 T,
                 ifelse(
                   input$Infin450k == F & input$InfinEPIC == F,
                   T,
                   F
                 )
               ))
        
  
      , ]
      
  },
  extensions = c("Buttons", "Scroller"),
  options = list(
    dom = 'Bfrtip',
    deferRender = TRUE,
    scrollY = 400,
    scrollX = 400,
    scroller = TRUE,
    buttons = c('copy', 'csv', 'excel'))
  )
  
  

}

shinyApp(ui = ui, server = server)








