evap <- function(function_call) {
  call <- substitute(function_call)
  reactiveVals$call_list <<- append(reactiveVals$call_list, call) # in the shiny app, we don't have to make this global call, I think
  return_value <- eval(call, envir = parent.frame())
  return(return_value)
}

reactiveVals <- list()
input <- list()
reactiveVals$sce <- readRDS('data/pbmc/sce.rds')
source("functions/vis_functions.R")
library(CATALYST)
library(shiny)

input$dr <- 'UMAP'
input$ncells <- 100
reactiveVals$sce <- evap(reactiveVals$sce <- runDR(reactiveVals$sce, input$dr, ncells = input$ncells), list(dr=input$dr, ncells=input$ncells))