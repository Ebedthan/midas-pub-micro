source("renv/activate.R")

# My functions and theme -------------------------------------------------------
theme_Publication <-
  function(base_size = 14,
           base_family = "helvetica") {
    require(grid)
    require(ggthemes)
    require(ggplot2)
    (
      ggthemes::theme_foundation(base_size = base_size, base_family = base_family) +
        ggplot2::theme(
          plot.title = element_text(
            face = "bold",
            size = rel(1.2),
            hjust = 0.5
          ),
          text = element_text(),
          panel.background = element_rect(colour = NA),
          plot.background = element_rect(colour = NA),
          panel.border = element_rect(colour = NA),
          axis.title = element_text(face = "bold", size = rel(1)),
          axis.title.y = element_text(angle = 90, vjust = 2),
          axis.title.x = element_text(vjust = -0.2),
          axis.text = element_text(),
          axis.line = element_line(colour = "black"),
          axis.ticks = element_line(),
          panel.grid.major = element_line(colour = "#f0f0f0"),
          panel.grid.minor = element_blank(),
          legend.key = element_rect(colour = NA),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.key.size = unit(0.2, "cm"),
          legend.margin = margin(t = 0, unit = "cm"),
          legend.title = element_text(face = "italic"),
          plot.margin = unit(c(10, 5, 5, 5), "mm"),
          strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
          strip.text = element_text(face = "bold")
        )
    )
  }

scale_fill_Publication <- function(...) {
  library(scales)
  ggplot2::discrete_scale("fill",
                          "Publication",
                          manual_pal(
                            values = c(
                              "#386cb0",
                              "#fdb462",
                              "#7fc97f",
                              "#ef3b2c",
                              "#662506",
                              "#a6cee3",
                              "#fb9a99",
                              "#984ea3",
                              "#ffff33"
                            )
                          ), ...)
}

scale_colour_Publication <- function(...) {
  library(scales)
  ggplot2::discrete_scale("colour",
                          "Publication",
                          manual_pal(
                            values = c(
                              "#386cb0",
                              "#fdb462",
                              "#7fc97f",
                              "#ef3b2c",
                              "#662506",
                              "#a6cee3",
                              "#fb9a99",
                              "#984ea3",
                              "#ffff33"
                            )
                          ), ...)
}

export_plot <- function(myplot, file, width = 8, height = 6, units = "in", type = "tiff", bg = "white") {
  require(Cairo)
  require(cli)
  
  # Export plot using Cairo
  Cairo::Cairo(width = width,
               height = height,
               file = file,
               type = type,
               bg = bg,
               units = units,
               dpi = 300)
  plot(myplot)
  dev.off()
  
  cli::cli_alert_success(paste("Exporting plot to", file, "succeeded"))
  
  # Show plot
  myplot
}
