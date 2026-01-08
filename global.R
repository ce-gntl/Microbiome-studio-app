library(shiny)
library(DT)
library(bslib)
library(microbiome)
library(microbiomeutilities)
library(ggpubr)
library(phyloseq)
library(ggvenn)
library(ANCOMBC)
library(tidyr)
library(reticulate)
library(hrbrthemes)
library(ggplot2)
library(vegan)
library(htmltools)
library(rlang)
library(bsicons)
library(spsComps)
library(grid)
library(ggrepel)
library(forcats)
library(scales)
library(shinydashboard)

###------Datatable options-------

options(DT.options = list(
  
  paging = TRUE,    ## paginate the output
  pageLength = 15,  ## number of rows to output for each page
  scrollX = TRUE,   ## enable scrolling on X axis
  scrollY = TRUE,   ## enable scrolling on Y axis
  autoWidth = FALSE))

tabPanel<- function(...) {
  shiny::tabPanel(..., class = "p-3")
}