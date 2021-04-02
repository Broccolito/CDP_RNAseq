library(shiny)
library(shinydashboard)
library(shinydashboardPlus)
library(shinyWidgets)
library(dashboardthemes)

library(dplyr)
library(ggplot2)
library(plotly)
library(ggfortify)

ui = dashboardPage(
  header = dashboardHeader(),
  sidebar = dashboardSidebar(
    textAreaInput(inputId = "genes", 
                  label = "Gene Number List", height = 500,
                  placeholder = "ENSG00000000419\nENSG00000000457\nENSG00000000460",
                  value = "ENSG00000000419\nENSG00000000457\nENSG00000000460"),
    switchInput(inputId = "load", label = "PC", value = TRUE),
    actionBttn(inputId = "drawPCA", label = "Run", style = "simple", block = FALSE)
  ),
  body = dashboardBody(
    
    shinyDashboardThemes(
      theme = "grey_dark"
    ),
    
    boxPlus(plotlyOutput("PCA_Plot12"), 
            solidHeader = TRUE, collapsible = TRUE, width = 6, title = "PC1 and PC2"),
    boxPlus(plotlyOutput("PCA_Plot23"), 
            solidHeader = TRUE, collapsible = TRUE, width = 6, title = "PC2 and PC3"),
    boxPlus(plotlyOutput("PCA_Plot34"), 
            solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, 
            width = 6, title = "PC3 and PC4"),
    boxPlus(plotlyOutput("PCA_Plot45"), 
            solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, 
            width = 6, title = "PC4 and PC5")
    
  )
)

server = function(input, output, session) {
  load("shiny_data.RData")
  cdp$Health_Status = as.factor(ifelse(cdp$CMS == 1, "CMS", "NON-CMS"))
  cdp$Exercise_Status = as.factor(ifelse(cdp$Exercise == 1, "EXERCISE", "REST"))
  cdp$Fasting_Status = as.factor(ifelse(cdp$Fasting == 1, "FASTING", "NON-FASTING"))
  cdp$Treatment = as.factor(ifelse(cdp$Hemodilution == 1, "HEMODILUTION", "NON-HEMODILUTION"))
  
  gene_list_to_PCA = function(gene_list, first_pc = 1, second_pc = 2){
    df = select(cdp, all_of(gene_list))
    plt = ggplotly(autoplot(prcomp(df), data = cdp, x = first_pc, y = second_pc,
                            colour = "Health_Status", shape = "Treatment",
                            loadings = input$load, loadings.label = input$load,
                            loadings.colour = 'black',
                            loadings.label.size = 3))
    return(plt)
  }
  
  observeEvent(input$drawPCA, {
    gene_list = unlist(strsplit(input$genes, "\n"))
    gene_list = unlist(strsplit(gene_list, "\t"))
    output$PCA_Plot12 = renderPlotly({
      plt12 = gene_list_to_PCA(gene_list, first_pc = 1, second_pc = 2)
      plt12
    })
    
    output$PCA_Plot23 = renderPlotly({
      plt23 = gene_list_to_PCA(gene_list, first_pc = 2, second_pc = 3)
      plt23
    })
    
    output$PCA_Plot34 = renderPlotly({
      plt34 = gene_list_to_PCA(gene_list, first_pc = 3, second_pc = 4)
      plt34
    })
    
    output$PCA_Plot45 = renderPlotly({
      plt45 = gene_list_to_PCA(gene_list, first_pc = 4, second_pc = 5)
      plt45
    })
  })
  
  onSessionEnded(function(){stopApp()})
}

shinyApp(ui, server)