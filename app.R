library(shiny)
library(readxl)
library(dplyr)
library(ggplot2)
library(grid)
library(png)

# --- Load global variables ---
source("globals.R")  # must define: my_cols, my_alpha, size_lims, size_breaks, col_low, col_high, col_lims, meninges_width, meninges_height, meninges_bottom

options(shiny.launch.browser = TRUE)

# --- Upload module UI ---
uploadUI <- function(id) {
  ns <- NS(id)
  tagList(
    fileInput(ns("file1"), "Upload first file (structures.csv)", accept = c(".csv", ".xlsx", ".xls")),
    fileInput(ns("file2"), "Upload second file (angles.csv)", accept = c(".csv", ".xlsx", ".xls")),
    actionButton(ns("submit"), "Submit", class = "btn-primary"),
    textOutput(ns("status"))
  )
}

# --- Upload module server ---
uploadServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    
    data <- reactiveValues(df1 = NULL, df2 = NULL, df = NULL, text_df = NULL, ready = FALSE)
    
    observeEvent(input$submit, {
      req(input$file1, input$file2)
      
      read_file <- function(file) {
        df <- tryCatch(read.csv(file$datapath, stringsAsFactors = FALSE),
                       error = function(e) NULL)
        if(is.null(df)) {
          df <- tryCatch(read_excel(file$datapath),
                         error = function(e) stop("Cannot read file as CSV or Excel."))
        }
        df
      }
      
      data$df1 <- read_file(input$file1)
      data$df2 <- read_file(input$file2)
      
      structures <- data$df1
      angles <- data$df2
      
      # --- Compute offsets ---
      angles$Offset0 <- -angles$True0
      angles$Offset90 <- 90 - angles$True90
      angles$Offset180 <- 180 - angles$True180
      
      # Join structures and angles 
      df <- structures %>% 
        left_join(angles) %>% 
        select(-starts_with("True"))
      
      # Apply offset by sinus
      df$Offset[df$Sinus == 0] <- df$Offset0[df$Sinus == 0]
      df$Offset[df$Sinus == 90] <- df$Offset90[df$Sinus == 90]
      df$Offset[df$Sinus == 180] <- df$Offset180[df$Sinus == 180]
      
      # Compute coordinates
      df <- df %>%
        mutate(Angle = Angle + Offset,
               Length = Length / 1000,
               Size = Size / 1000) %>%
        mutate(X = Length * cos(Angle * pi / 180),
               Y = Length * sin(Angle * pi / 180)) %>%
        select(-starts_with("Offset"))
      
      df$Phenotype <- as.factor(df$Phenotype)
      
      # Compute summary for text labels
      data$text_df <- df %>%
        group_by(Group) %>%
        summarise(Mice = n_distinct(TissueID),
                  Count = n()) %>%
        mutate(Label = paste(Count, "ELS from", Mice, "samples"))
      
      data$df <- df
      data$ready <- TRUE
    })
    
    return(data)
  })
}

# --- Plotting function ---
make_tls_plot <- function(df_group, text_group, figure,
                          side_lim, upper_lim, lower_lim) {
  
  gg <- ggplot(df_group, aes(x = X, y = Y)) +
    lims(x = c(-side_lim, side_lim), y = c(lower_lim, upper_lim)) +
    geom_bin2d(binwidth = c(2.3, 2.3), boundary = c(-1, -1)) +
    coord_fixed() +
    annotation_custom(
      rasterGrob(figure, width = unit(1,"npc"), height = unit(1,"npc")),
      -meninges_width, meninges_width,
      -meninges_bottom, meninges_height
    ) +
    geom_point(aes(colour = Phenotype, size = Size), alpha = my_alpha) +
    scale_colour_manual(values = my_cols, name = "ELS type") +
    scale_size_continuous(limits = size_lims, breaks = size_breaks, name = "Size (\u00B5m²)") +
    scale_fill_gradient2(low = col_low, high = col_high, limits = col_lims, name = "ELS count") +
    labs(x = "Distance from Lambda (mm)", y = "Distance from Lambda (mm)", fill = "Count") +
    guides(
      colour = guide_legend(order = 1),
      size   = guide_legend(order = 2),
      fill   = guide_legend(order = 3)
    ) +
    theme_classic()
  
  if(nrow(text_group) > 0) {
    gg <- gg + geom_text(
      data = text_group,
      aes(label = Label),
      x = -side_lim + 2.4,
      y = upper_lim + 0.3,
      size = 3
    )
  }
  
  return(gg)
}

# --- Shiny UI ---
ui <- fluidPage(
  h2("Density map ELS"),
  uiOutput("main_ui"),
  hr(),
  uiOutput("group_ui"),
  uiOutput("return_button_ui"),
  uiOutput("download_buttons"),
  plotOutput("plot_tls", height = "600px")
)

# --- Shiny server ---
server <- function(input, output, session) {
  
  step <- reactiveVal(1)
  uploaded <- uploadServer("upload")
  
  # --- Main UI ---
  output$main_ui <- renderUI({
    if(step() == 1) uploadUI("upload")
    else tagList(h4("Processing complete"), p("Select a group to view the plot below."))
  })
  
  # --- Observe upload completion ---
  observeEvent(uploaded$ready, {
    req(uploaded$ready)
    step(2)
  })
  
  # --- Group selection (step 2 only) ---
  output$group_ui <- renderUI({
    req(step() == 2, uploaded$ready)
    selectInput("group", "Select Group", choices = unique(uploaded$df$Group))
  })
  
  # --- Return to upload button (step 2 only) ---
  output$return_button_ui <- renderUI({
    req(step() == 2)
    actionButton("return_upload", "Return to Upload", class = "btn-warning")
  })
  
  observeEvent(input$return_upload, {
    step(1)
    uploaded$ready <- FALSE
    uploaded$df <- NULL
    uploaded$text_df <- NULL
    output$group_ui <- renderUI({})
    output$download_buttons <- renderUI({})
  })
  
  # --- Download buttons ---
  output$download_buttons <- renderUI({
    req(step() == 2, uploaded$ready, !is.null(input$group))
    tagList(
      downloadButton("download_plot", "Download Selected Plot as PDF"),
      downloadButton("download_all_plots", "Download All Groups as PDF")
    )
  })
  
  # --- Render plot ---
  output$plot_tls <- renderPlot({
    req(step() == 2, uploaded$ready, !is.null(input$group))
    
    df_group <- uploaded$df %>% filter(Group == input$group)
    text_group <- uploaded$text_df %>% filter(Group == input$group)
    figure <- png::readPNG("data/opac_meninges.png")
    
    buffer <- 2
    upper_lim <- max(meninges_height, max(df_group$Y)) + buffer
    lower_lim <- min(-meninges_bottom, min(df_group$Y)) - buffer
    side_lim <- max(meninges_width, max(abs(df_group$X))) + buffer
    
    make_tls_plot(df_group, text_group, figure, side_lim, upper_lim, lower_lim)
  })
  
  # --- Download handlers ---
  output$download_plot <- downloadHandler(
    filename = function() {
      paste0("TLS_plot_", input$group, "_", Sys.Date(), ".pdf")
    },
    content = function(file) {
      req(uploaded$ready, input$group)
      
      df_group <- uploaded$df %>% filter(Group == input$group)
      text_group <- uploaded$text_df %>% filter(Group == input$group)
      figure <- png::readPNG("data/opac_meninges.png")
      
      buffer <- 2
      upper_lim <- max(df_group$Y, meninges_height) + buffer
      lower_lim <- min(df_group$Y, -meninges_bottom) - buffer
      side_lim  <- max(abs(df_group$X), meninges_width) + buffer
      
      p <- make_tls_plot(
        df_group, text_group, figure,
        side_lim, upper_lim, lower_lim
      )
      
      ggsave(
        filename = file,
        plot = p,
        device = cairo_pdf,
        width = 8,
        height = 6,
        units = "in"
      )
    }
  )
  
  
  output$download_all_plots <- downloadHandler(
    filename = function() {
      paste0("TLS_all_groups_", Sys.Date(), ".pdf")
    },
    content = function(file) {
      req(uploaded$ready)
      
      figure <- png::readPNG("data/opac_meninges.png")
      
      # open ONE cairo PDF device
      cairo_pdf(
        filename = file,
        width = 8,
        height = 6
      )
      
      for (group in unique(uploaded$df$Group)) {
        
        df_group <- uploaded$df %>% filter(Group == group)
        text_group <- uploaded$text_df %>% filter(Group == group)
        
        buffer <- 2
        upper_lim <- max(df_group$Y, meninges_height) + buffer
        lower_lim <- min(df_group$Y, -meninges_bottom) - buffer
        side_lim  <- max(abs(df_group$X), meninges_width) + buffer
        
        p <- make_tls_plot(
          df_group,
          text_group,
          figure,
          side_lim,
          upper_lim,
          lower_lim
        ) +
          ggtitle(paste("Group:", group))
        
        print(p)
      }
      
      dev.off()
    }
  )
  
}

# --- Run app ---
shinyApp(ui, server)
