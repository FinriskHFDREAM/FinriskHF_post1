library(dplyr)
library(ggplot2)
library(forcats)

args=(commandArgs(TRUE))
args <- "~/Documents/other_project/DREAM/Rerun_2025"
PARAM <- list()
#this argument is to specify the path of input folder
PARAM$folder.R <- paste0(args[1])
PARAM$folder.data <- paste0(PARAM$folder.R, "/")

# load all the files
cox_data <-read.csv(paste0(PARAM$folder.data, "input/Cox_features.tsv"),sep="\t")
rsf_data <- read.csv(paste0(PARAM$folder.data, "input/RSF_features.tsv"),sep="\t")

cox_plot <- cox_data %>%
  mutate(type = "CoxPH",
         value = HR,
         lower = CI_lower,
         upper = CI_upper,
         pval=p)

rsf_plot <- rsf_data %>%
  mutate(type = "RSF",
         value = Importance,
         lower = NA,
         upper = NA)

plot_df <- bind_rows(cox_plot, rsf_plot)

# Set consistent feature order based on average value
plot_df$features <- fct_reorder(plot_df$features, plot_df$value, .fun = mean, .desc = TRUE)

# === Plot ===
plot_df <- bind_rows(cox_plot, rsf_plot)

# Set consistent feature order based on average value
plot_df$features <- fct_reorder(plot_df$features, plot_df$value, .fun = mean, .desc = TRUE)
# Create zebra stripe data for background shading
feature_levels <- levels(plot_df$features)
# Add significance stars to CoxPH data
plot_df <- plot_df %>%
  mutate(signif_label = case_when(
    type == "CoxPH" & p < 0.001 ~ "***",
    type == "CoxPH" & p < 0.01  ~ "**",
    type == "CoxPH" & p < 0.05  ~ "*",
    type == "CoxPH"             ~ "",
    TRUE                        ~ NA_character_
  ))

plot_df <- plot_df %>%
  mutate(alpha = ifelse(type == "CoxPH" & p >= 0.05, 0.4, 1))

zebra_df <- data.frame(
  features = feature_levels[seq(1, length(feature_levels), by = 2)],
  ymin = seq_along(feature_levels)[seq(1, length(feature_levels), by = 2)] - 0.5,
  ymax = seq_along(feature_levels)[seq(1, length(feature_levels), by = 2)] + 0.5
)
# === Plot ===

model_colors <- c(
  "SB2_v2_enetCoxPh" = "#018571",   # orange
  "SB2_v2_lassoCoxPh"  = "#d7191c",   # light blue
  "SB2_v3_enetCoxPh" = "#2c7bb6",   # green
  "SB2_v3_lassoCoxPh"  = "#fd8861",   # reddish orange
  "SB2_v2_enetRSF" = "#018571",   # purple/pink
  "SB2_v2_lassoRSF"  = "#d7191c",   # strong blue
  "SB2_v3_enetRSF"  = "#2c7bb6",   # yellow
  "SB2_v3_lassoRSF" = "#fd8861",   # gray
  "SB2_v2_allRSF"      = "#000000"    # black
)
ggplot(plot_df, aes(x = value, y = features, color = Model)) +
  geom_rect(data = zebra_df, 
            aes(ymin = ymin, ymax = ymax), 
            xmin = -Inf, xmax = Inf, 
            fill = "gray95", inherit.aes = FALSE)+
  # HR dots with error bars
  geom_point(data = filter(plot_df, type == "CoxPH"),
             size = 1, position = position_dodge(width = 0.7)) +
  
  geom_errorbarh(data = filter(plot_df, type == "CoxPH"),
                 aes(xmin = lower, xmax = upper),
                 height = 0.2,
                 position = position_dodge(width = 0.7)) +
  
  # RSF bars
  geom_col(data = filter(plot_df, type == "RSF"),
           aes(fill = Model),
           color=NA,
           width = 0.5,
           position = position_dodge(width = 0.7),
           alpha = 0.7) +
  scale_fill_manual(values = model_colors) +
  scale_color_manual(values = model_colors) +
  facet_grid(. ~ type, scales = "free_x", space = "free_x") +
  theme_classic() +
  labs(
    x = "Effect Size (HR or RSF Importance)",
    y = NULL,
  ) +
  theme(
    axis.text.y = element_text(margin = margin(r = 10),size = 12),  # Extra space between labels
    axis.text.x = element_text(size = 12),
    strip.text = element_text(size = 12),
    panel.spacing.y = unit(0.5, "lines"),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 11),
    legend.spacing.x = unit(0.6, 'cm'),
    legend.box.just = "left"
  ) +
  guides(
    color = guide_legend(nrow = 2),
    fill = guide_legend(nrow = 2)
  ) +
  geom_vline(data = data.frame(type = "CoxPH"), aes(xintercept = 1),
             linetype = "dashed", color = "gray40", inherit.aes = FALSE)
ggsave(filename = paste0(PARAM$folder.data, "/output/",
                                    "Plot_features_SB2.pdf"),
       width = 12, height = 10) # Try 8–10 height
