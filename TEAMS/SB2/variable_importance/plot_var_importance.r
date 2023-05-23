require(ggplot2)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

refs <- c("ref_v1", "ref_v2", "ref_v3", "ref_v4")

lapply(refs, function(ref) {
  
  models <- readRDS(file.path("../", ref, "models.rds"))
  
  x <- as.data.frame(unlist(models))
  
  xx <- data.frame(
    model = sapply(strsplit(rownames(x), ".", fixed = TRUE), "[", 1),
    variable = sapply(strsplit(rownames(x), ".", fixed = TRUE), "[", 2),
    beta = abs(x[, 1])
  )
  
  xx$type <- "covariates"
  xx$type[grep("cluster_", xx$variable)] <- "microbiome"
  xx$type[grep("p__", xx$variable)] <- "microbiome"
  
  plot <- 
    ggplot(
      xx,
      aes(
        x = reorder(variable, beta),
        y = beta,
        color = type)) +
    geom_point() +
    geom_segment(aes(xend = variable, yend = 0)) +
    coord_flip() + 
    facet_grid(~model) +
    theme_bw() +
    theme(
      axis.text.y = element_text(size = 7),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.position = "bottom"
    )
  plot
  
  ggsave(
    plot = plot,
    filename = "plot_VI.pdf",
    device = "pdf",
    path = file.path("../", ref),
    width = 5,
    height = 6
  )
  
})


