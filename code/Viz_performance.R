library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(reshape2)
library(scales)
library(patchwork)
library(ggrepel)


args=(commandArgs(TRUE))
args <- "~/Documents/other_project/DREAM/Rerun_2025"
PARAM <- list()
#this argument is to specify the path of input folder
PARAM$folder.R <- paste0(args[1])
PARAM$folder.data <- paste0(PARAM$folder.R, "/")

# load all the files
bootstrap_Harrell_c <- read.csv(paste0(PARAM$folder.data, "input/bootstrap_Harrell_c_baseline1020.csv"))[,-1]
bootstrap_hoslem <- read.csv(paste0(PARAM$folder.data,"input/bootstrap_hoslem_baseline1020.csv"))[,-1]
team_list <- read.csv(paste0(PARAM$folder.data,"input/team_list.csv"))
Final_leaderboard1 <- as.data.frame(read.csv(paste0(PARAM$folder.data, "input/Final_leaderboard.csv")))
Final_leaderboard <- merge(team_list, Final_leaderboard1, by = "id", all.x = TRUE) 
# prepare dataset
Final_leaderboard$`Teams_Participants` <- c("Baseline Age-Sex", "Baseline Covariates","Ensemble",
                                            "Baseline All","Metformin-121", "YG-HZ team", 
                                            "DFH","TristanF","SB2", "Pteam", 
                                            "UTKteam"
                                            )
Final_leaderboard$harrell_c.r <- round(Final_leaderboard$harrell_c, digits = 3)
Final_leaderboard$hoslem_sci<- scientific_format(digits=2)(Final_leaderboard$hoslem_test)
Final_leaderboard_arr <- Final_leaderboard %>% arrange(-harrell_c.r)
limits <- Final_leaderboard_arr$Teams_Participants

df.harrelc=as.data.frame(t(bootstrap_Harrell_c))
df.harrelc$folder_name <- rownames(df.harrelc)
df.harrelc=melt(df.harrelc)
colnames(df.harrelc)  <- c("folder_name", 
                           "variable", 
                           "Bootstrapped Harrell's C Index (N=1000)")
df.harrelc <- df.harrelc %>%
  left_join(Final_leaderboard[, c("folder_name", "Teams_Participants")], by = "folder_name")
df.harrelc$Teams_Participants <- factor(df.harrelc$Teams_Participants, levels = limits)


df.hoslem = as.data.frame(t(bootstrap_hoslem))
df.hoslem$folder_name <- rownames(df.hoslem)
df.hoslem=melt(df.hoslem)

colnames(df.hoslem)  <- c("folder_name", 
                          "variable", 
                          "Bootstrapped Hosmer-Lemeshow Test (N=1000)")
df.hoslem <- df.hoslem %>%
  left_join(Final_leaderboard[, c("folder_name", "Teams_Participants")], by = "folder_name")
df.hoslem$Teams_Participants <- factor(df.hoslem$Teams_Participants, levels = limits)

# Step 1: Summarise bootstrap distributions into medians + percentiles

summary_harrelc <- df.harrelc %>%
  group_by(Teams_Participants) %>%
  summarise(
    C_index_median = median(`Bootstrapped Harrell's C Index (N=1000)`, na.rm = TRUE),
    C_index_lower = quantile(`Bootstrapped Harrell's C Index (N=1000)`, 0.05, na.rm = TRUE),
    C_index_upper = quantile(`Bootstrapped Harrell's C Index (N=1000)`, 0.95, na.rm = TRUE),
    .groups = "drop"
  )

summary_hoslem <- df.hoslem %>%
  group_by(Teams_Participants) %>%
  summarise(
    HL_median = median(`Bootstrapped Hosmer-Lemeshow Test (N=1000)`, na.rm = TRUE),
    HL_lower = quantile(`Bootstrapped Hosmer-Lemeshow Test (N=1000)`, 0.05, na.rm = TRUE),
    HL_upper = quantile(`Bootstrapped Hosmer-Lemeshow Test (N=1000)`, 0.95, na.rm = TRUE),
    .groups = "drop"
  )

# Step 2: Merge summaries
df_plot <- summary_harrelc %>%
  left_join(summary_hoslem, by = "Teams_Participants") %>%
  mutate(Teams_Participants = factor(Teams_Participants, levels = limits))

df_plot <- df_plot %>%
  left_join(Final_leaderboard_arr[, c("Teams_Participants", "harrell_c.r", "hoslem_test")],
            by = "Teams_Participants") %>%
  mutate(Teams_Participants = factor(Teams_Participants, levels = limits))

# Plot
plot_all <- ggplot(df_plot, aes(x = harrell_c.r, y = hoslem_test, label = Teams_Participants)) +
  geom_point(size = 3, color = "steelblue") +  # Point from leaderboard
  geom_errorbarh(aes(xmin = C_index_lower, xmax = C_index_upper), height = 0.002) +
  geom_errorbar(aes(ymin = HL_lower, ymax = HL_upper), width = 0.01) +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "gray40") +
  theme_minimal(base_size = 14) +
  labs(
    x = "Harrell's C-index",
    y = "Hosmer-Lemeshow p-value (log scale)",
    title = "Model Discrimination vs Calibration (N=1000 Bootstraps)"
  ) +
  geom_text(vjust = -1.2, size = 3.5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(plot_all,filename = paste0(PARAM$folder.data, "/output/",
                                  "Plot_harrel_hoslem_all.pdf"),
       width = 8, height = 4)

#plot top team and baseline

selected_teams <- c("Baseline All", 
                    "Baseline Covariates", 
                    "Baseline Age-Sex", 
                    "DFH", 
                    "SB2", 
                    "Ensemble")

df_plot_subset <- df_plot %>%
  filter(Teams_Participants %in% selected_teams) %>%
  mutate(Teams_Participants = factor(Teams_Participants, levels = selected_teams))
df_plot_subset <- df_plot_subset %>%
  mutate(color_group = case_when(
    grepl("SB2", Teams_Participants) ~ "SB2",
    grepl("DFH", Teams_Participants) ~ "DFH",
    grepl("Baseline All", Teams_Participants) ~ "Baseline",
    grepl("Baseline Age-Sex", Teams_Participants) ~ "Baseline",
    grepl("Baseline Covariates", Teams_Participants) ~ "Baseline",
    grepl("Ensemble", Teams_Participants) ~ "Ensemble",
    TRUE ~ "Other"
  ))

color_palette <- c(
  "SB2" = "#2ca02c",       # Green
  "DFH" = "#ff7f0e",       # Orange
  "Baseline" = "#9467bd",  # Purple
  "Ensemble" = "#1f77b4",  # Dark blue
  "Other" = "gray60"
)

plot_subset <- ggplot(df_plot_subset, aes(x = harrell_c.r, y = hoslem_test, label = Teams_Participants)) +
  # Error bars first, so they appear behind the dots
  geom_errorbarh(aes(xmin = C_index_lower, xmax = C_index_upper, color = color_group),
                 height = 0.5, size = 0.4, linewidth = 0.4,width=0.3) +
  geom_errorbar(aes(ymin = HL_lower, ymax = HL_upper, color = color_group),
                width = 0.002, size = 0.4, linewidth = 0.4) +
  # Main dots
  geom_point(aes(color = color_group), size = 4) +
  # Add labels
  #geom_text(vjust = -0.2, hjust=-0.1, size = 4,face = "bold") +
  ggrepel::geom_text_repel(
    aes(label = Teams_Participants),
    size = 4,
    fontface = "bold",
    max.overlaps = 1000,        # allows more flexibility
    box.padding = 0.3,          # padding around labels
    point.padding = 0.2,        # padding around points
    segment.color = "gray60",   # color of connecting lines
    segment.size = 0.2          # thin connecting lines
  )+
  # Formatting axes
  scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_color_manual(values = color_palette) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "gray40", linewidth = 0.5) +
  # Aesthetics
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    panel.grid.minor = element_blank()
  ) +
  #theme(
  #  legend.position = "none",
  #  axis.text.x = element_text(angle = 45, hjust = 1),
  #  panel.grid.minor = element_blank()
  #) +
  labs(
    x = "Harrell's C-index",
    y = "Hosmer-Lemeshow p-value (log scale)"
  )

ggsave(plot_subset,filename = paste0(PARAM$folder.data, "/output/",
                                  "Plot_harrel_hoslem_subset.pdf"),
       width = 10, height = 10)


#Plot the refined model
#load refined model file
bootstrap_Harrell_c_ref <- read.csv(paste0(PARAM$folder.data, "input/bootstrap_Harrell_c_refined.csv"))[,-1]
bootstrap_hoslem_ref <- read.csv(paste0(PARAM$folder.data,"input/bootstrap_hoslem_refined.csv"))[,-1]
team_list_ref <- read.csv(paste0(PARAM$folder.data,"input/mapping_name.csv"))
Final_leaderboard_ref1 <- as.data.frame(read.csv(paste0(PARAM$folder.data, "input/refined_leaderboard.csv")))
Final_leaderboard_ref <- merge(team_list_ref, Final_leaderboard_ref1, by = "Folder_name", all.x = TRUE) 
Final_leaderboard_ref$harrell_c.r <- round(Final_leaderboard_ref$harrell_c, digits = 3)
Final_leaderboard_ref$hoslem_sci<- scientific_format(digits=2)(Final_leaderboard_ref$hoslem_test)
Final_leaderboard_ref_arr <- Final_leaderboard_ref %>% arrange(-harrell_c.r)
limitsR <- Final_leaderboard_ref_arr$Model_name

df.harrelc_ref=as.data.frame(t(bootstrap_Harrell_c_ref))
df.harrelc_ref$Model_name <- rownames(df.harrelc_ref)
df.harrelc_ref=melt(df.harrelc_ref)
colnames(df.harrelc_ref)  <- c("Model_name", 
                           "variable", 
                           "Bootstrapped Harrell's C Index (N=1000)")
df.harrelc_ref <- df.harrelc_ref %>%
  left_join(Final_leaderboard_ref[, c("Folder_name", "Model_name")], by = "Model_name")
df.harrelc_ref$Model_name <- factor(df.harrelc_ref$Model_name, levels = limitsR)


df.hoslem_ref = as.data.frame(t(bootstrap_hoslem_ref))
df.hoslem_ref$Model_name <- rownames(df.hoslem_ref)
df.hoslem_ref=melt(df.hoslem_ref)
colnames(df.hoslem_ref)  <- c("Model_name", 
                          "variable", 
                          "Bootstrapped Hosmer-Lemeshow Test (N=1000)")
df.hoslem_ref <- df.hoslem_ref %>%
  left_join(Final_leaderboard_ref[, c("Folder_name", "Model_name")], by = "Model_name")
df.hoslem_ref$Model_name <- factor(df.hoslem_ref$Model_name, levels = limitsR)

# Step 1: Summarise bootstrap distributions into medians + percentiles

summary_harrelc <- df.harrelc_ref %>%
  group_by(Model_name) %>%
  summarise(
    C_index_median = median(`Bootstrapped Harrell's C Index (N=1000)`, na.rm = TRUE),
    C_index_lower = quantile(`Bootstrapped Harrell's C Index (N=1000)`, 0.05, na.rm = TRUE),
    C_index_upper = quantile(`Bootstrapped Harrell's C Index (N=1000)`, 0.95, na.rm = TRUE),
    .groups = "drop"
  )

summary_hoslem <- df.hoslem_ref %>%
  group_by(Model_name) %>%
  summarise(
    HL_median = median(`Bootstrapped Hosmer-Lemeshow Test (N=1000)`, na.rm = TRUE),
    HL_lower = quantile(`Bootstrapped Hosmer-Lemeshow Test (N=1000)`, 0.05, na.rm = TRUE),
    HL_upper = quantile(`Bootstrapped Hosmer-Lemeshow Test (N=1000)`, 0.95, na.rm = TRUE),
    .groups = "drop"
  )

# Step 2: Merge summaries
df_plot_ref <- summary_harrelc %>%
  left_join(summary_hoslem, by = "Model_name") %>%
  mutate(Model_name = factor(Model_name, levels = limitsR))

df_plot_ref <- df_plot_ref %>%
  left_join(Final_leaderboard_ref_arr[, c("Model_name", "harrell_c","harrell_c.r", "hoslem_test")],
            by = "Model_name") %>%
  mutate(Model_name = factor(Model_name, levels = limitsR))


df_plot_ref_final <- merge(team_list_ref, df_plot_ref, by = "Model_name", all.x = TRUE) 

# Duplicate SB2_15 for both groups
sb2_extra <- df_plot_ref_final %>%
  filter(Folder_name == "SB2_15") %>%
  slice(rep(1, 2)) %>%  # duplicate row
  mutate(Group = c("SB2_v2", "SB2_v3"))  # assign to two groups

# Bind back into the original dataframe
df_plot_ref_final <- df_plot_ref_final %>%
  filter(Folder_name != "SB2_15") %>%
  bind_rows(sb2_extra)
write.csv(df_plot_ref_final, paste0(PARAM$folder.data, "/output/",
                                    "refined_model_df.csv"))
# Make sure factor levels are clean for consistent shape/color legends
df_plot_ref_final$Strategy <- factor(df_plot_ref_final$Strategy)
df_plot_ref_final$Features_selection <- factor(df_plot_ref_final$Features_selection)
df_plot_ref_final$Group <- factor(df_plot_ref_final$Group)

# Find best performing model (highest C-index) per group
df_labels <- df_plot_ref_final %>%
  group_by(Group) %>%
  slice_max(harrell_c, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(label_text = paste0("C=", round(harrell_c.r, 3), "\nHL=", formatC(hoslem_test, format = "e", digits = 2)))

# Plot
plot_facet <- ggplot(df_plot_ref_final, aes(x = harrell_c.r, y = hoslem_test,
                    color = Strategy, shape = Features_selection)) +
  
  # Error bars
  geom_errorbarh(aes(xmin = C_index_lower, xmax = C_index_upper), height = 0.5, size = 0.4, linewidth = 0.4,width=0.3) +
  geom_errorbar(aes(ymin = HL_lower, ymax = HL_upper), width = 0.002, size = 0.4, linewidth = 0.4) +
  
  # Main dots
  geom_point(size = 3) +
  # Best model labels
  geom_text(
    data = df_labels,
    aes(x = harrell_c.r, y = hoslem_test, label = label_text),
    color = "black",
    size = 3.3,
    hjust = -0.2,
    vjust = -0.1,
    inherit.aes = FALSE
  )+
  
  # Add HL p-value threshold line
  geom_hline(yintercept = 0.04, linetype = "dashed", color = "gray40", linewidth = 0.5) +
  
  # Scales and facets
  scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_color_brewer(palette = "Dark2") +
  facet_wrap(~ Group, ncol = 3, scales = "free_x") +
  
  # Aesthetics
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "Harrell's C-index",
    y = "Hosmer-Lemeshow p-value (log scale)"
  )
ggsave(plot_facet,filename = paste0(PARAM$folder.data, "/output/",
                                     "Plot_harrel_hoslem_facet_refined.pdf"),
       width = 10, height = 10)
final_plot <- plot_subset / plot_facet + plot_layout(heights = c(0.8, 0.8)) +plot_annotation(tag_levels = "A")
ggsave(final_plot,filename = paste0(PARAM$folder.data, "/output/",
                                    "Plot_harrel_hoslem_final2.pdf"),
       width = 10, height = 13)
