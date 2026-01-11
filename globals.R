library(ggplot2)
library(patchwork)
library(dplyr)
library(png)
library(grid)
library(purrr)

outlier_threshold <- 13


my_cols <- c("#008ECE", "orange", "red")
figure <- readPNG("data/opac_meninges.png")

my_alpha <- 2/3
my_label_size <- 2
my_label_colour <- "black"

size_lims <- c(0, 330)
size_breaks <- seq(0, 300, by = 100)
col_low <- "white"
col_high <- "black"
col_lims <- c(0, 50)
bin_num <- 4
bin_width <- 3

my_palette <- 8
my_dir <- 1

asp_ratio <- 0.92
meninges_height <- 10.1
meninges_bottom <- meninges_height / 4
meninges_width <- (meninges_height + meninges_bottom) * asp_ratio / 2

get_inset <- function(df, my_cols){
  
  sub_df <- df %>%
    group_by(Phenotype) %>%
    summarise(Amount = 100 * round(n() / nrow(.), 3))
  
  p <- ggplot(sub_df, aes(x = "", y = Amount, fill = Phenotype)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y", start = 0) +
    scale_fill_manual(values = my_cols) +
    geom_text(aes(label = paste0(Amount, "%")),
              position = position_stack(vjust = 0.5), size = 2) +
    theme_classic() +
    theme(panel.background = element_rect(fill = "transparent"),
          plot.background = element_rect(fill = "transparent", colour = NA),
          axis.title = element_blank(),
          axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none")
  
  return(p)
}

annotation_custom2 <- function(grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, data) {
  
  layer(data = data, stat = StatIdentity, position = PositionIdentity, 
        geom = ggplot2:::GeomCustomAnn,
        inherit.aes = TRUE, params = list(grob = grob, 
                                          xmin = xmin, xmax = xmax, 
                                          ymin = ymin, ymax = ymax))
}