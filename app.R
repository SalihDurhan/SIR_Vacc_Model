# 
# Shiny app for SIR model with vaccinations
# Info, Style, Localization changes and improvements over 
# https://github.com/ekstroem/Shiny-Vaccine
# 
# Linkteki shiny uygulamasi uzerine yapilan degisiklik ve eklerle olusturulmustur
#
# Salih Durhan

library("shiny")
library("deSolve")
library("ggplot2")
library("tidyverse")
library("ggrepel")
library("shinydashboard")
library("ggthemes")
library("shinyWidgets")

# Differential equations defining the SIR model
# SIR modelini tanimlayan diferansiyel denklemler

sir <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS <- -beta * S * I
    dI <-  beta * S * I - gamma * I
    dR <-                 gamma * I
    dV <- 0
    return(list(c(dS, dI, dR, dV)))
  })
}


#
# Shiny User Interface 
# Shiny Arayuzu
#

ui <- dashboardPage(
  dashboardHeader(title = "SIR Salgın Modelleri"),
  dashboardSidebar(
    
    
    checkboxGroupInput("add_var", "Diğer Değişkenler:",
                       c("Aşılı" = "V",
                         "Enfeksiyona Açık" = "S",
                         "İyileşmiş" = "R")),
    chooseSliderSkin("HTML5", color = "#DC322F"),
    sliderInput("popsize",
                "Nüfus (Milyon):",
                min = 1, max = 300, value = 83
    ),
    sliderInput("connum",
                "Temel Üreme Katsayısı (R0):",
                min = .5, max = 20, value = 3
    ),
    
    sliderInput("infper",
                "Enfeksiyon Süresi (Gün):",
                min = 1, max = 30, value = 14
    ),
    sliderInput("pinf",
                "Başlangıç Enfeksiyon Sayısı:",
                min = 1, max = 5000, value = 2
    ),
    sliderInput("pvac",
                "Aşılı Nüfus Oranı (%):",
                min = 0, max = 100, value = 0
    ),
    sliderInput("vaceff",
                "Aşı Verimliliği (%):",
                min = 0, max = 100, value = 85
    ),
    
    sliderInput("timeframe",
                "Zaman Aralığı (Gün):",
                min = 1, max = 400, value = 300
    )
    
  ),
  dashboardBody(
    tags$head(tags$style(HTML('
                              /* logo */
                                .skin-blue .main-header .logo {
                              background-color: #073642;
                              }
                              
                              /* logo when hovered */
                              .skin-blue .main-header .logo:hover {
                              background-color: #073642;
                              }
                              
                              /* navbar (rest of the header) */
                              .skin-blue .main-header .navbar {
                              background-color: #073642;
                              }
                              
                              /* main sidebar */
                              .skin-blue .main-sidebar {
                              background-color: #073642;
                              }
                              
                              /* active selected tab in the sidebarmenu */
                              .skin-blue .main-sidebar .sidebar .sidebar-menu .active a{
                              background-color: #073642;
                              }
                              
                              /* other links in the sidebarmenu */
                              .skin-blue .main-sidebar .sidebar .sidebar-menu a{
                              background-color: #073642;
                              color: #000000;
                              }
                              
                              /* other links in the sidebarmenu when hovered */
                              .skin-blue .main-sidebar .sidebar .sidebar-menu a:hover{
                              background-color: #073642;
                              }
                              /* toggle button when hovered  */
                              .skin-blue .main-header .navbar .sidebar-toggle:hover{
                              background-color: #073642;
                              }
                              
                              /* body */
                              .content-wrapper, .right-side {
                              background-color: #002B36;
                              }



                              '))),
    
    #    mainPanel(
    fluidRow(plotOutput("distPlot")),
    br(),
    fluidRow(
      # Dynamic valueBoxes 
      valueBoxOutput("Box1", width = 6),
      valueBoxOutput("Box2", width = 6),
      valueBoxOutput("Box3", width = 6),
      valueBoxOutput("Box4", width = 6)
    ),
    br(),
    br()
    )
    )

#
# Shiny server definition 
# Shiny sunucu tanimi
# 
server <- function(input, output) {
  # Create reactive input
  dataInput <- reactive({
    init       <-
      c(
        S = 1 - input$pinf / (input$popsize*1000000) - input$pvac / 100 * input$vaceff / 100,
        I = input$pinf /  (input$popsize*1000000),
        R = 0,
        V = input$pvac / 100 * input$vaceff / 100
      )
    ## beta: infection parameter; gamma: recovery parameter
    ## beta: enfeksiyon parametresi; gamma: iyilesme parametresi
    parameters <-
      c(beta = input$connum * 1 / input$infper,
        # * (1 - input$pvac/100*input$vaceff/100),
        gamma = 1 / input$infper)
    ## Time frame
    ## Zaman
    times      <- seq(0, input$timeframe, by = .2)
    
    ## Solve using ode (General Solver for Ordinary Differential Equations)
    ## SIR model denklemlerini ode kutuphanesi ile coz
    out <- ode(
      y = init,
      times = times,
      func = sir,
      parms = parameters
    )   
    #    out
    as.data.frame(out)
  })
  
  # Data keys
  # Veri tanimlari
  output$distPlot <- renderPlot({
    out <-
      dataInput() %>%
      gather(key, value, -time) %>%
      mutate(
        id = row_number(),
        key2 = recode(
          key,
          S = "Enfeksiyona Açık (S)",
          I = "Hastalanmış (I)",
          R = "İyileşmiş (R)",
          V = "Aşılı (V)"
        ),

        keyright = recode(
          key,
          S = "",
          I = "Hastalanmış (I)",
          R = "",
          V = ""
        )
      )
    
    # Graph
    # Grafik
    
    ggplot(data =  out[out$key %in% c("I", input$add_var),], #out[out$key == 'I' ,],
           aes(
             x = time,
             y = value,
             col = key2
           )) + 
      ylab("Toplam Nüfusa Oran") + xlab("Zaman (Gün)") +
      geom_line(size = 1, key_glyph = "rect") +
      geom_text_repel(
        data = subset(out, time == max(time)),
        aes(label = keyright),
        size = 6,
        segment.size  = 0.2,
        segment.color = "grey50",
        nudge_x = 0,
        hjust = 1,
        direction = "y"
      ) +
      
      scale_colour_manual(values = c("#DC322F", "#B58900", "#268BD2", "#2AA198" )) +
      scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
      theme_solarized(light = FALSE) + 
      theme(
        legend.position= c(0.9, 0.8),
        legend.title=element_blank(),
        legend.text = element_text(colour="white", size=10, face="bold"),
        axis.title.x = element_text(colour = "white"),
        axis.title.y = element_text(colour = "white"),
        axis.text.x = element_text(colour="white"), 
        axis.text.y = element_text(colour="white")
      )
    
  })
  
  
  output$Box1 <- renderValueBox({
    valueBox(
      tags$p("Salgının Tepe Noktası", style = "font-size: 80%;"),
      tags$p(dataInput() %>% filter(I == max(I)) %>% select(time) %>% mutate(time = round(time, 0)) %>% paste("Gün "), style = "font-size: 220%;"), 
      icon = icon("burn"),
      color = "olive"  
      #red, yellow, aqua, blue, light-blue, green, navy, teal, olive, lime, orange, fuchsia, purple, maroon, black.
    )
  })
  
  output$Box2 <- renderValueBox({
    valueBox(
      tags$p("Tepe Noktasındaki Hasta Sayısı", style = "font-size: 80%;"), 
      tags$p(dataInput() %>% filter(I == max(I)) %>% select(I) %>% mutate(I = round(I * input$popsize, 2)) %>% paste("Milyon"), style = "font-size: 220%;"), 
      icon = icon("medkit"),
      color = "orange"
    )
  })
  
  output$Box3 <- renderValueBox({
    valueBox(
      tags$p("Sürü Bağışıklığı için Gerekli Aşılı Nüfus Oranı", style = "font-size: 80%;"),
      tags$p(paste0(round(100 * (1 - 1 / (input$connum)), 2), "%"), style = "font-size: 220%;"),
      icon = icon("syringe"),
      color = "olive"
    )
  })
  
  output$Box4 <- renderValueBox({
    valueBox(
      tags$p("Salgın Süresince Toplam Hasta Sayısı", style = "font-size: 80%;"), 
      tags$p(dataInput() %>% filter(time == max(time)) %>% select(R) %>% mutate(R = round(R * input$popsize, 2)) %>% paste("Milyon"), style = "font-size: 220%;"),
      icon = icon("medkit"),
      color = "orange"
    )
  })  
}

shinyApp(ui = ui, server = server)