library(rsconnect)
library(shiny)
library(tidyverse)
library(shinyjs)
library(DT)
library(waiter)
library(shinyalert)
library(shinyFeedback)
library(ggrepel)
library(shinydashboard)
library(directlabels)

source("maximin design.R")
source("Numerical derivation optimal proportions.R")
source("Optimal values fixed N.R")
source("Maximin Design N fixed.R")



ui <- dashboardPage(
  dashboardHeader(title = "Optimal allocation to treatments in a SMART", titleWidth = "500px"),
  dashboardSidebar(width = "350px",
                   sidebarMenu(
                     menuItem("Introduction and overview", tabName = "introduction", icon = icon("home")),
                     menuItem("Locally optimal design", tabName = "optimalproportions", icon = icon("list-alt"),
                              menuSubItem("Fixed sample size", tabName = "optimalnfixed"),
                              menuSubItem("Fixed budget", tabName = "optimalcost")),
                     menuItem("Maximin Design", tabName = "maximin", icon = icon("list-alt"),
                              menuSubItem("Fixed sample size", tabName = "maximinnfixed"),
                              menuSubItem("Fixed budget", tabName = "maximincost")),
                     br(),
                     br(),
                     br(),
                     br(),
                     tags$li(a(href = "https://www.uu.nl/en/",
                               img(src = "UU-logo.png", height = "80px",width = "300px")),
                            class = "dropdown")
                   )
  ),
  dashboardBody(useShinyjs(),
                useShinyFeedback(),
                use_waiter(),
                useShinyalert(),
                
                tabItems(
                  tabItem(tabName = "introduction",
                          fluidRow(column(10,offset = 1,
                                          h3("Introduction"),
                                          
                                          span("Version 0.3.1., created by Andrea Morciano."),
                                          br(),
                                          br(),
                                          span("This Shiny App is designed to help users derive optimal proportions of patients assigned
                                               to treatments in a SMART design and thus obtain the optimal design."),
                                          br(),
                                          span("Using the prototypical SMART design as shown in the figure below, users are asked 
                                               to specify the response rates for the first-stage
                                               treatment options and the weights used in the multiple objective optimal design. 
                                               The formula for the multiple objective optimal design is given by:"),
                                          withMathJax("$$\\Phi = Var(\\bar{Y}_{d_{1}} - \\bar{Y}_{d_{3}})\\lambda_{13} +
                                                      Var(\\bar{Y}_{d_{1}} - \\bar{Y}_{d_{4}})\\lambda_{14} +
                                                      Var(\\bar{Y}_{d_{2}} - \\bar{Y}_{d_{3}})\\lambda_{23} +
                                                      Var(\\bar{Y}_{d_{2}} - \\bar{Y}_{d_{4}})\\lambda_{24}$$"),
                                          h3("Locally Optimal Design: Overview"),
                                          div("The first tab allows calculation of the locally optimal design. 
                                              This optimal design is compared to the balanced design in terms of 
                                              relative efficiency and also contour plots for relative efficiency 
                                              are shown. Two scenarios are considered: one where sample size is 
                                              fixed and cost of treatments is not taken into account, and another 
                                              in which cost of treatments is considered and sample size depends 
                                              on a budget constraint. In both scenarios, the response rates and 
                                              the weights for the multiple-objective optimal design should be 
                                              specified. Furthermore, the total sample size should be specified 
                                              in the first scenario, while the cost per  treatment and the budget 
                                              should be specified in the second scenario.", style = "text-align: justify;"),
                                          br(),
                                          h3("Maximin Design: Overview"),
                                          div("The second tab allows calculation of the maximin design, which 
                                              provides a solution to the problem of the unknown response rates 
                                              for the first-stage treatment options. The maximin design is 
                                              compared to the balanced design in terms of minimum relative 
                                              efficiency and contour plots for the relative efficiency are shown.
                                              As for the locally optimal design, a fixed total sample size or
                                              a budgetary constraint is used.",style = "text-align: justify;"),
                                          br(),
                                      
                                          h3("Example of a SMART design"),
                                          br(),
                                         div(imageOutput("SMART"), style="text-align: center;"),
                                         br(),
                                          div(tags$b("A SMART design example:"), tags$i("A scheme of the prototypical SMART design. Circled 
                                                                                        “R” denotes randomization at each stage. p1 and (1-p1)
                                                                                        are, respectively, the probabilities of receiving 
                                                                                        first-stage treatments A and B. p2 and (1-p2) are, 
                                                                                        respectively, the probabilities of receiving 
                                                                                        second-stage treatments D and E for non-responders 
                                                                                        starting with first-stage treatment A. p3 and (1-p3)
                                                                                        are, respectively, the probabilities of receiving 
                                                                                        second-stage treatments G and H for non-responders 
                                                                                        starting with first-stage treatment B. γ1 and γ2 
                                                                                        indicate, respectively, response rates for the 
                                                                                        first-stage treatments A and B."),style="text-align: justify;"),
                                          br(),
                                          br(),
                                          br(),
                                          br()
                                          ))),
                  tabItem(tabName = "optimalnfixed",
                              fluidRow(box(title = "Fixed sample size", width =12, status = "primary", solidHeader = T,
                                          dataTableOutput(outputId = "opt"),br(),br())),
                              fluidRow(column(width = 4,
                                              box(title = "", width = 12, status = "primary", solidHeader = T,
                                              div(
                                                id = "side-panel3",
                                                h4("Enter response rates:"),
                                                
                                                numericInput(inputId = "g1a", label = paste0(intToUtf8(947),1),
                                                             value = 0.50,
                                                             min = 0, max = 1, step = 0.05),
                                                
                                                numericInput(inputId = "g2a", label = paste0(intToUtf8(947),2),
                                                             value = 0.50,
                                                             min = 0, max = 1, step = 0.05),
                                                h4("Enter N:"),
                                                
                                                numericInput(inputId = "N", label = "Sample size",
                                                             value = 1000,
                                                             min = 10, max = 10000, step = 10),
                                                
                                                h4("Enter weights:"),
                                                tags$i(h5("Weights must be different from zero and must sum to 1")),
                                                
                                                numericInput(inputId = "w13a", label = paste0(intToUtf8(955),"13"), 
                                                             value = 0.70, 
                                                             min = 0, max = 1, step = 0.05),
                                                
                                                numericInput(inputId = "w14a", label = paste0(intToUtf8(955),"14"), 
                                                             value = 0.10, 
                                                             min = 0, max = 1, step = 0.05),
                                                
                                                numericInput(inputId = "w23a", label = paste0(intToUtf8(955),"23"), 
                                                             value = 0.10, 
                                                             min = 0, max = 1, step = 0.05),
                                                
                                                numericInput(inputId = "w24a", label = paste0(intToUtf8(955),"24"), 
                                                             value = 0.10, 
                                                             min = 0, max = 1, step = 0.05)),
                                              tags$i( h4("Submit")),
                                              
                                              actionButton("nfixed", "Submit", icon = icon("paper-plane", lib = "font-awesome"), width = "250px"),
                                              br(),
                                             
                                             tags$i( h4("Reset fields to default")),
                                              actionButton("reset3", "Reset fields", icon = icon("sync-alt", lib = "font-awesome"), width = "250px"))),
                                       column(width = 8,
                                              box(title = "Relative efficiency", width = 12, status = "primary", solidHeader = T,
                                                  uiOutput(outputId = "r_eff_nf"),br()),
                                              box(title = "Contour plots for Relative Efficiency", width = 12, status = "primary", solidHeader = T,
                                                  fluidRow(
                                                    splitLayout(cellWidths = c("33%","33%","33%"),
                                                                plotOutput("contourp1n"), plotOutput("contourp2n"),plotOutput("contourp3n"))),
                                                  br(),
                                                  br()),
                                              box(title = "Elapsed time", width = 12, status = "primary", solidHeader = T,
                                                  uiOutput(outputId = "time3") )))
                              
                            
                          ),
                  tabItem(tabName = "optimalcost",
                        
                              fluidRow(
                              box(title = "Fixed budget", width = 12, status = "primary", solidHeader = T,
                                  dataTableOutput(outputId = "optimals"),
                                  br(),
                                  br())),
                              fluidRow(column(width = 4,
                                              box(title = "",width = 12, status = "primary", solidHeader = T,
                                              div(
                                                id = "side-panel1",
                                                h4("Enter response rates:"),
                                                
                                                numericInput(inputId = "g1", label = paste0(intToUtf8(947),1),
                                                             value = 0.50,
                                                             min = 0, max = 1, step = 0.05),
                                                
                                                numericInput(inputId = "g2", label = paste0(intToUtf8(947),2),
                                                             value = 0.50,
                                                             min = 0, max = 1, step = 0.05),
                                                
                                                h4("Enter cost of treatments:"),
                                                
                                                numericInput(inputId = "Ca", label = "First-stage treatment A:",
                                                             value = 50,
                                                             min = 0, max = 10000, step = 10),
                                                
                                                numericInput(inputId = "Cb", label = "First-stage treatment B:",
                                                             value = 50,
                                                             min = 0, max = 10000, step = 10),
                                                
                                                numericInput(inputId = "Cc", label = "Second-stage treatment responders C:",
                                                             value = 50,
                                                             min = 0, max = 10000, step = 10),
                                                
                                                numericInput(inputId = "Cf", label = "Second-stage treatment responders F:",
                                                             value = 50,
                                                             min = 0, max = 10000, step = 10),
                                                
                                                numericInput(inputId = "Cd", label = "Second-stage treatment non-responders D:",
                                                             value = 50,
                                                             min = 0, max = 10000, step = 10),
                                                
                                                numericInput(inputId = "Ce", label = "Second-stage treatment non-responders E:",
                                                             value = 50,
                                                             min = 0, max = 10000, step = 10),
                                                
                                                numericInput(inputId = "Cg", label = "Second-stage treatment non-responders G:",
                                                             value = 50,
                                                             min = 0, max = 10000, step = 10),
                                                
                                                numericInput(inputId = "Ch", label = "Second-stage treatment non-responders H:",
                                                             value = 50,
                                                             min = 0, max = 10000, step = 10),
                                                
                                                numericInput(inputId = "C", label = "Budget:",
                                                             value = 100000,
                                                             min = 0, max = 10000000, step = 1000), 
                                                
                                                h4("Enter weights:"),
                                                tags$i(h5("Weights must be different from zero and must sum to 1")),
                                                
                                                numericInput(inputId = "w13", label = paste0(intToUtf8(955),"13"), 
                                                             value = 0.70, 
                                                             min = 0, max = 1, step = 0.05),
                                                
                                                numericInput(inputId = "w14", label = paste0(intToUtf8(955),"14"), 
                                                             value = 0.10, 
                                                             min = 0, max = 1, step = 0.05),
                                                
                                                numericInput(inputId = "w23", label = paste0(intToUtf8(955),"23"), 
                                                             value = 0.10, 
                                                             min = 0, max = 1, step = 0.05),
                                                
                                                numericInput(inputId = "w24", label = paste0(intToUtf8(955),"24"), 
                                                             value = 0.10, 
                                                             min = 0, max = 1, step = 0.05)),
                                              tags$i( h4("Submit")),
                                              
                                              actionButton("optval", "Submit", icon = icon("paper-plane", lib = "font-awesome"), width = "250px"),
                                              br(),
                                              tags$i(h4("Reset fields to default")),
                                              actionButton("reset1", "Reset fields", icon = icon("sync-alt", lib = "font-awesome"), width = "250px"))
                              ),
                                       column(width = 8,
                                              box(title = "Relative efficiency",width = 12,status = "primary", solidHeader = T,
                                                  uiOutput(outputId = "r_eff"),
                                                  br()),
                                              br(),
                                              box(title = "Contour plots for Relative Efficiency",width = 12, status = "primary", solidHeader = T,
                                                  fluidRow(
                                                    splitLayout(cellWidths = c("33%","33%","33%"),
                                                                plotOutput("contourp1"), plotOutput("contourp2"),plotOutput("contourp3"))),
                                                  br()),
                                              br(),
                                              box(title = "Elapsed time", width = 12,status = "primary", solidHeader = T,
                                                  uiOutput(outputId = "time1") )))
                              
                            
                            
                           ),
                  tabItem(tabName = "maximinnfixed",
        
                              fluidRow(box(title = "Fixed sample size", width = 12, status = "primary", solidHeader = T,
                                           dataTableOutput(outputId = "mmdnf"),
                                           br(),
                                           br())),
                              fluidRow(column(width = 4,
                                              box(title = "", width = 12, status = "primary", solidHeader = T,
                                                  div(
                                                    id = "side-panel4",
                                                    h4("Enter response rates:"),
                                                    
                                                    numericInput(inputId = "g1b", label = paste0(intToUtf8(947),1),
                                                                 value = 0.50,
                                                                 min = 0, max = 1, step = 0.05),
                                                    
                                                    numericInput(inputId = "g2b", label = paste0(intToUtf8(947),2),
                                                                 value = 0.50,
                                                                 min = 0, max = 1, step = 0.05),
                                                    h4("Enter N:"),
                                                    
                                                    numericInput(inputId = "Nb", label = "Sample size",
                                                                 value = 1000,
                                                                 min = 10, max = 10000, step = 10),
                                                    
                                                    h4("Enter weights:"),
                                                    tags$i(h5("Weights must be different from zero and must sum to 1")),
                                                    
                                                    numericInput(inputId = "w13b", label = paste0(intToUtf8(955),"13"), 
                                                                 value = 0.70, 
                                                                 min = 0, max = 1, step = 0.05),
                                                    
                                                    numericInput(inputId = "w14b", label = paste0(intToUtf8(955),"14"), 
                                                                 value = 0.10, 
                                                                 min = 0, max = 1, step = 0.05),
                                                    
                                                    numericInput(inputId = "w23b", label = paste0(intToUtf8(955),"23"), 
                                                                 value = 0.10, 
                                                                 min = 0, max = 1, step = 0.05),
                                                    
                                                    numericInput(inputId = "w24b", label = paste0(intToUtf8(955),"24"), 
                                                                 value = 0.10, 
                                                                 min = 0, max = 1, step = 0.05)),
                                                  tags$i( h4("Submit")),
                                                  
                                                  actionButton("mmdbutton", "Submit", icon = icon("paper-plane", lib = "font-awesome"), width = "250px"),
                                                  br(),
                                                  tags$i(h4("Reset fields to default")),
                                                  actionButton("reset4", "Reset fields", icon = icon("sync-alt", lib = "font-awesome"), width = "250px"))),
                                       column(width = 8,
                                              box(title = "Relative efficiency", width = 12, status = "primary", solidHeader = T,
                                                  uiOutput(outputId = "r_eff_mmd_nf"),
                                                  br()),
                                              br(),
                                              box(title = "Contour plots for Relative Efficiency", width = 12, status = "primary", solidHeader = T,
                                                  fluidRow(
                                                    splitLayout(cellWidths = c("33%","33%","33%"),
                                                                plotOutput("contourp1mmdnf"), plotOutput("contourp2mmdnf"),plotOutput("contourp3mmdnf"))),
                                                  br()),
                                              br(),
                                              box(title = "Elapsed time", width = 12, status = "primary", solidHeader = T,
                                                  uiOutput(outputId = "time4"))
                                       )
                          )),
                          
                  tabItem(tabName = "maximincost",
                     
                              fluidRow(
                                box(title = "Fixed budget", width = 12, status = "primary", solidHeader = T,
                                    dataTableOutput(outputId = "mmd"),
                                    br(),
                                    br())
                              ),
                              fluidRow(column(width = 4,
                                              box(title = "", width = 12, status = "primary", solidHeader = T,
                                                  div(
                                                    id = "side-panel2",
                                                    h4("Enter response rates:"),
                                                    
                                                    numericInput(inputId = "g11", label = paste0(intToUtf8(947),1),
                                                                 value = 0.50,
                                                                 min = 0, max = 1, step = 0.05),
                                                    
                                                    numericInput(inputId = "g21", label = paste0(intToUtf8(947),2),
                                                                 value = 0.50,
                                                                 min = 0, max = 1, step = 0.05),
                                                    
                                                    h4("Enter cost of treatments:"),
                                                    
                                                    
                                                    numericInput(inputId = "Ca1", label = "First-stage treatment A:",
                                                                 value = 50,
                                                                 min = 0, max = 10000, step = 10),
                                                    
                                                    numericInput(inputId = "Cb1", label = "First-stage treatment B:",
                                                                 value = 50,
                                                                 min = 0, max = 10000, step = 10),
                                                    
                                                    numericInput(inputId = "Cc1", label = "Second-stage treatment responders C:",
                                                                 value = 50,
                                                                 min = 0, max = 10000, step = 10),
                                                    
                                                    numericInput(inputId = "Cf1", label = "Second-stage treatment responders F:",
                                                                 value = 50,
                                                                 min = 0, max = 10000, step = 10),
                                                    
                                                    numericInput(inputId = "Cd1", label = "Second-stage treatment non-responders D:",
                                                                 value = 50,
                                                                 min = 0, max = 10000, step = 10),
                                                    
                                                    numericInput(inputId = "Ce1", label = "Second-stage treatment non-responders E:",
                                                                 value = 50,
                                                                 min = 0, max = 10000, step = 10),
                                                    
                                                    numericInput(inputId = "Cg1", label = "Second-stage treatment non-responders G:",
                                                                 value = 50,
                                                                 min = 0, max = 10000, step = 10),
                                                    
                                                    numericInput(inputId = "Ch1", label = "Second-stage treatment non-responders H:",
                                                                 value = 50,
                                                                 min = 0, max = 10000, step = 10),
                                                    
                                                    numericInput(inputId = "C1", label = "Budget:",
                                                                 value = 100000,
                                                                 min = 0, max = 10000000, step = 1000), 
                                                    
                                                    h4("Enter weights:"),
                                                    tags$i(h5("Weights must be different from zero and must sum to 1.")),
                                                    
                                                    numericInput(inputId = "w131", label = paste0(intToUtf8(955),"13"), 
                                                                 value = 0.70, 
                                                                 min = 0, max = 1, step = 0.05),
                                                    
                                                    numericInput(inputId = "w141", label = paste0(intToUtf8(955),"14"), 
                                                                 value = 0.10, 
                                                                 min = 0, max = 1, step = 0.05),
                                                    
                                                    numericInput(inputId = "w231", label = paste0(intToUtf8(955),"23"), 
                                                                 value = 0.10, 
                                                                 min = 0, max = 1, step = 0.05),
                                                    
                                                    numericInput(inputId = "w241", label = paste0(intToUtf8(955),"24"), 
                                                                 value = 0.10, 
                                                                 min = 0, max = 1, step = 0.05)),
                                                  tags$i( h4("Submit")),
                                                  
                                                  actionButton("mmdopt", "Submit", icon = icon("paper-plane", lib = "font-awesome"), width = "250px"),
                                                  br(),
                                                  tags$i(h4("Reset fields to default")),
                                                  actionButton("reset2", "Reset fields", icon = icon("sync-alt", lib = "font-awesome"), width = "250px"))),
                                       column(width = 8,
                                              box(title = "Relative efficiency", width = 12, status = "primary", solidHeader = T,
                                                  uiOutput(outputId = "r_eff_mmd"),
                                                  br()),
                                              br(),
                                              box(title = "Contour plots for Minimum Relative Efficiency", width = 12, status = "primary", solidHeader = T,
                                                  fluidRow(
                                                    splitLayout(cellWidths = c("33%","33%","33%"),
                                                                plotOutput("contourp1mmd"), plotOutput("contourp2mmd"),plotOutput("contourp3mmd"))),
                                                  br()),
                                              br(),
                                              box(title = "Elapsed time", width = 12, status = "primary", solidHeader = T,
                                                  uiOutput(outputId = "time2"))
                                              ))
                         
                            
                          )
                                          )
                          )
                )


server <- function(input, output, session){
  
 output$SMART <- renderImage({
    list(src = "www/SMARTexample.png", alt = paste("SMART example"), width = "700px", height = "400px")
  }, deleteFile = F)
  #######################################Analytical output######################################
  anal_values <- reactiveValues()
  
  table_cache <- matrix(nrow = 0, ncol = 20)
  colnames(table_cache) <- c(paste0(intToUtf8(947),c(1:2)), "Ca", "Cb", "Cc", "Cd", "Ce", "Cf", "Cg", "Ch",
                             "p1.opt", "p2.opt", "p3.opt", "N.opt",
                             paste(intToUtf8(981), "value"), "Cost", paste0(intToUtf8(955),c("13","14","23","24")))
  table_cache <- DT::datatable(table_cache)
  
  
  observeEvent(input$optval, {
    
    show_waiter(
      shiny::tagList(
        spin_folding_cube(),
        span(h3("Loading ..."), style = "color:white;"),
        br(),
        br(),
        br(),
        span("NOTE: This may take some minutes.", style = "color:white;")
      )
    )
    
    anal_values$all <- analytic.prop(input$g1, input$g2, input$Ca, input$Cb, input$Cc,
                                     input$Cd, input$Ce, input$Cf, input$Cg, input$Ch, input$w13,
                                     input$w14, input$w23, input$w24, input$C)
    
    
    anal_values$final <- datatable(rbind(anal_values$final$x$data,anal_values$all$final), rownames = F)
    
    
    anal_values$r.eff <- anal_values$all$r.eff
    anal_values$C.eff <- anal_values$all$cost.eff
    
    anal_values$time <- anal_values$all$time
    
    anal_values$cont1 <- anal_values$all$cont_plot1
    anal_values$cont2 <- anal_values$all$cont_plot2
    anal_values$cont3 <- anal_values$all$cont_plot3
    
    
    hide_waiter()    
    shinyalert("Computation completed!", type = "success")
    
  })
  
  output$optimals <- renderDataTable({
    wgts <- input$w13 + input$w14 + input$w23 + input$w24
    
    validate(need(near(wgts,1),"Sum of weights must be equal to 1."))
    
    anal_values$final
  })
  
  output$r_eff <- renderUI({
    wgts <- input$w13 + input$w14 + input$w23 + input$w24
    
    validate(need(near(wgts,1),"Sum of weights must be equal to 1."))
    
    paste("Relative efficiency of the balanced design with respect to the optimal design is:", 
          paste0(round(as.numeric(anal_values$r.eff), digits = 4),"."),"The balanced design would cost",
          format(round(as.numeric(anal_values$C.eff),2), nsmall = 2, big.mark = " "), 
          "€ to perform as the optimal design, given an initial budget of:", 
          format(as.integer(anal_values$final$x$data[,16]), big.mark = " "),
          "Note: this refers to the last row in the above table.")
  })
  
  
  
  output$contourp1 <- renderPlot({
    wgts <- input$w13 + input$w14 + input$w23 + input$w24
    
    validate(need(near(wgts,1),"Sum of weights must be equal to 1."))
    anal_values$cont1
  })
  
  output$contourp2 <- renderPlot({
    wgts <- input$w13 + input$w14 + input$w23 + input$w24
    
    validate(need(near(wgts,1),"Sum of weights must be equal to 1."))
    anal_values$cont2
  })
  
  output$contourp3 <- renderPlot({
    wgts <- input$w13 + input$w14 + input$w23 + input$w24
    
    validate(need(near(wgts,1),"Sum of weights must be equal to 1."))
    anal_values$cont3
  })
  
  output$time1 <- renderUI({
    paste("The function took", round(as.numeric(anal_values$time),2), "seconds.")
  })
  
  observeEvent(input$reset1, {
    shinyjs::reset("side-panel1")
  })
  ######################################Optimal proportions N fixed Output######################################
  opt_val <- reactiveValues()
  
  table_cache2 <- matrix(nrow = 0, ncol = 11)
  colnames(table_cache2) <- c(paste0(intToUtf8(947),c(1:2)), "p1.opt", "p2.opt", "p3.opt", 
                              paste(intToUtf8(981), "value"), "N", paste0(intToUtf8(955),c("13","14","23","24")))
  
  table_cache2 <- DT::datatable(table_cache2)
  
  observeEvent(input$nfixed ,{
    show_waiter(
      shiny::tagList(
        spin_folding_cube(),
        span(h3("Loading ..."), style = "color:white;"),
        br(),
        br(),
        br(),
        span("NOTE: This may take some minutes.", style = "color:white;")
      )
    )
    
    opt_val$all <- optimal(input$g1a, input$g2a, input$w13a, input$w14a, input$w23a, input$w24a, input$N)
    
    opt_val$opt <- datatable(rbind(opt_val$opt$x$data,opt_val$all$opt), rownames = F)
    
    opt_val$r.eff <- opt_val$all$r.eff
    
    opt_val$n.eff <- opt_val$all$N.eff
    
    opt_val$time <- opt_val$all$time
    
    opt_val$cont1 <- opt_val$all$cont_plot1
    opt_val$cont2 <- opt_val$all$cont_plot2
    opt_val$cont3 <- opt_val$all$cont_plot3
    
    hide_waiter()
    
    shinyalert("Computation completed!", type = "success")
  })
  
  output$opt <- renderDataTable({
    wgtsa <- input$w13a + input$w14a + input$w23a + input$w24a
    
    validate(need(near(wgtsa,1),"Sum of weights must be equal to 1."))
    opt_val$opt
  })
  
  output$r_eff_nf <- renderUI({
    wgtsa <- input$w13a + input$w14a + input$w23a + input$w24a
    
    validate(need(near(wgtsa,1),"Sum of weights must be equal to 1."))
    
    paste("Relative efficiency of the balanced design with respect to the optimal design is:", 
          paste0(round(as.numeric(opt_val$r.eff),digits = 4),"."),"The balanced design would need",
          format(round(as.numeric(opt_val$n.eff) ,2), nsmall = 2, big.mark = " "), 
          "subjects to perform as the optimal design, given a sample size of:", 
          format(as.integer(opt_val$opt$x$data[,7]), big.mark = " "), ".",
          "Note: this refers to the last row in the above table.")
  })
  
  output$contourp1n <- renderPlot({
    wgts <- input$w13 + input$w14 + input$w23 + input$w24
    
    validate(need(near(wgts,1),"Sum of weights must be equal to 1."))
    opt_val$cont1
  })
  
  output$contourp2n <- renderPlot({
    wgts <- input$w13 + input$w14 + input$w23 + input$w24
    
    validate(need(near(wgts,1),"Sum of weights must be equal to 1."))
    opt_val$cont2
  })
  
  output$contourp3n <- renderPlot({
    wgts <- input$w13 + input$w14 + input$w23 + input$w24
    
    validate(need(near(wgts,1),"Sum of weights must be equal to 1."))
    opt_val$cont3
  })
  
  output$time3 <- renderUI({
    
    paste("The function took", round(as.numeric(opt_val$time),2), "seconds.")
  })
  
  observeEvent(input$reset3, {
    shinyjs::reset("side-panel3")
  })
  ###########################################MMD Output######################################
  mmd_values <- reactiveValues()
  
  table_cache1 <- matrix(nrow = 0, ncol = 21)
  colnames(table_cache1) <- c(paste0(intToUtf8(947),c(1:2)," Range"), "Ca", "Cb", "Cc", "Cd", "Ce", "Cf", "Cg", "Ch",
                              "p1.mmd", "p2.mmd", "p3.mmd", 
                              "N.mmd", paste(intToUtf8(981), "value"), "Cost", "MMV", 
                              paste0(intToUtf8(955),c("13","14","23","24")))
  
  table_cache1 <- DT::datatable(table_cache1)
  
  
  observeEvent(input$mmdopt, {
    
    show_waiter(
      shiny::tagList(
        spin_folding_cube(),
        span(h3("Loading ..."), style = "color:white;"),
        br(),
        br(),
        br(),
        span("NOTE: This may take some minutes.", style = "color:white;")
      )
    )
    
    mmd_values$mmm <- maximin_design(gamma1 = input$g11, gamma2 = input$g21, input$Ca1, input$Cb1, 
                                     input$Cc1, input$Cd1, input$Ce1, input$Cf1, input$Cg1, input$Ch1,
                                     input$w131,
                                     input$w141, input$w231, input$w241, input$C1, stepsize = 0.05)
    
    
    mmd_values$mmd_value <- datatable(rbind(mmd_values$mmd_value$x$data,mmd_values$mmm$mmd_value), rownames = F)
    
    #RE
    
    mmd_values$r.eff <- mmd_values$mmm$r.eff
    mmd_values$C.eff <- mmd_values$mmm$cost.eff
    
    #Contour plots
    
    mmd_values$cont_plot1 <- mmd_values$mmm$cont_plot1
    mmd_values$cont_plot2 <- mmd_values$mmm$cont_plot2
    mmd_values$cont_plot3 <- mmd_values$mmm$cont_plot3
    
    mmd_values$time <- mmd_values$mmm$time
    
    hide_waiter()
    
    shinyalert("Computation completed!", type = "success")
  })
  
  
  output$mmd <- renderDataTable({
    wgts1 <- input$w131 + input$w141 + input$w231 + input$w241
    
    validate(need(near(wgts1,1),"Sum of weights must be equal to 1."))
    mmd_values$mmd_value
  })
  
  output$r_eff_mmd <- renderUI({
    wgts1 <- input$w131 + input$w141 + input$w231 + input$w241
    
    validate(need(near(wgts1,1),"Sum of weights must be equal to 1."))
    
    paste("Relative efficiency of the balanced design with respect to the maximin design is:", 
          paste0(round(as.numeric(mmd_values$r.eff),digits = 4),"."),"The balanced design would cost",
          format(round(as.numeric(mmd_values$C.eff) ,2), nsmall = 2, big.mark = " "), 
          "€ to perform as the maximin design, given an initial budget of:", 
          format(as.integer(mmd_values$mmd_value$x$data[,16]), big.mark = " "),
          "Note: this refers to the last row in the above table.")
  })
  
  #Contour plots
  output$contourp1mmd <- renderPlot({
    wgts <- input$w13 + input$w14 + input$w23 + input$w24
    
    validate(need(near(wgts,1),"Sum of weights must be equal to 1."))
    mmd_values$cont_plot1
  })
  
  output$contourp2mmd <- renderPlot({
    wgts <- input$w13 + input$w14 + input$w23 + input$w24
    
    validate(need(near(wgts,1),"Sum of weights must be equal to 1."))
    mmd_values$cont_plot2
  })
  
  output$contourp3mmd <- renderPlot({
    wgts <- input$w13 + input$w14 + input$w23 + input$w24
    
    validate(need(near(wgts,1),"Sum of weights must be equal to 1."))
    mmd_values$cont_plot3
  })
  
  output$time2 <- renderUI({
    
    paste("The function took", round(as.numeric(mmd_values$time),2), "minutes.")
  })
  
  observeEvent(input$reset2, {
    shinyjs::reset("side-panel2")
  })
  
  ###########################################MMD N fixed Output######################################
  mmd_nf <- reactiveValues()
  
  table_cache3 <- matrix(nrow = 0, ncol = 12)
  colnames(table_cache3) <- c(paste0(intToUtf8(947),c(1:2)," Range"), "p1.mmd", "p2.mmd", "p3.mmd", 
                              "N", paste(intToUtf8(981), "value"), "MMV", 
                              paste0(intToUtf8(955),c("13","14","23","24")))
  
  table_cache3 <- DT::datatable(table_cache3)
  
  
  observeEvent(input$mmdbutton, {
    
    show_waiter(
      shiny::tagList(
        spin_folding_cube(),
        span(h3("Loading ..."), style = "color:white;"),
        br(),
        br(),
        br(),
        span("NOTE: This may take some minutes.", style = "color:white;")
      )
    )
    
    mmd_nf$all <- maximin_design_nfixed(gamma1 = input$g1b, gamma2 = input$g2b, input$w13b,
                                        input$w14b, input$w23b, input$w24b, input$Nb, stepsize = 0.05)
    
    
    mmd_nf$mmd_value <- datatable(rbind(mmd_nf$mmd_value$x$data,mmd_nf$all$mmd_value), rownames = F)
    
    #RE
    
    mmd_nf$r.eff <- mmd_nf$all$r.eff
    mmd_nf$N.eff <- mmd_nf$all$N.eff
    
    #Contour plots
    
    mmd_nf$cont_plot1 <- mmd_nf$all$cont_plot1
    mmd_nf$cont_plot2 <- mmd_nf$all$cont_plot2
    mmd_nf$cont_plot3 <- mmd_nf$all$cont_plot3
    
    mmd_nf$time <- mmd_nf$all$time
    
    hide_waiter()
    
    shinyalert("Computation completed!", type = "success")
  })
  
  
  output$mmdnf <- renderDataTable({
    wgtsb <- input$w13b + input$w14b + input$w23b + input$w24b
    
    validate(need(near(wgtsb,1),"Sum of weights must be equal to 1."))
    
    mmd_nf$mmd_value
  })
  
  output$r_eff_mmd_nf <- renderUI({
    wgtsb <- input$w13b + input$w14b + input$w23b + input$w24b
    
    validate(need(near(wgtsb,1),"Sum of weights must be equal to 1."))
    
    paste("Relative efficiency of the balanced design with respect to the maximin design is:", 
          paste0(round(as.numeric(mmd_nf$r.eff),digits = 4),"."),"The balanced design would need",
          format(round(as.numeric(mmd_nf$N.eff) ,2), nsmall = 2, big.mark = " "), 
          "subjects to perform as the maximin design, given a sample size of:", 
          format(as.integer(mmd_nf$mmd_value$x$data[,6]), big.mark = " "), ".",
          "Note: this refers to the last row in the above table.")
  })
  
  #Contour plots
  output$contourp1mmdnf <- renderPlot({
    wgts <- input$w13 + input$w14 + input$w23 + input$w24
    
    validate(need(near(wgts,1),"Sum of weights must be equal to 1."))
    mmd_nf$cont_plot1
  })
  
  output$contourp2mmdnf <- renderPlot({
    wgts <- input$w13 + input$w14 + input$w23 + input$w24
    
    validate(need(near(wgts,1),"Sum of weights must be equal to 1."))
    mmd_nf$cont_plot2
  })
  
  output$contourp3mmdnf <- renderPlot({
    wgts <- input$w13 + input$w14 + input$w23 + input$w24
    
    validate(need(near(wgts,1),"Sum of weights must be equal to 1."))
    mmd_nf$cont_plot3
  })
  
  output$time4 <- renderUI({
    
    paste("The function took", round(as.numeric(mmd_nf$time),2), "minutes.")
  })
  
  observeEvent(input$reset4, {
    shinyjs::reset("side-panel4")
  })
  
  session$onSessionEnded(stopApp)
}

shinyApp(ui = ui, server = server)

