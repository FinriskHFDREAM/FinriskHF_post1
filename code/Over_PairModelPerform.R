# load data
library(readr)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(ggupset)
library(tidyverse)
library(patchwork)
library(viridis)
library(ggpattern)
library(scales)
library(ggbreak)

# set parameters
mainDir <- getwd()
PARAM <- list()
PARAM$folder.R <- paste0(getwd(), "/")
PARAM$folder.result <- paste0(PARAM$folder.R, "output/")
PARAM$folder.data <- paste0(PARAM$folder.R, "final_results")
PARAM$folder.synthetic <- paste0(PARAM$folder.R) #change here accordingly
PARAM$folder.finrisk <- paste0(PARAM$folder.R) #change here accordingly


TimepointColor=c( "#E69F00","#8dd9c9","#fb9a99","#984ea3", "#1f78b4", "#beaed4", "#33a02c",  "#e31a1c")

# set parameters
values_color=c("#9970ab", "#9970ab", "#9970ab","#fb9a99", "#999999", 
               "#999999","#73baae", "#999999", "#999999", "#999999")


limits=c("Pteam", "UTKteam","Metformin-121","YG-HZ team","TristanF",
         "Baseline All", "Baseline Covariates",
         "Baseline Age-Sex", "DFH", "SB2")
labels=c("Pteam", "UTKteam","Metformin-121",
         "YG-HZ team","TristanF", "Baseline All", 
         "Baseline Covariates", "Baseline Age-Sex", "DFH", 
         "SB2")
# load all the files
bootstrap_Harrell_c <- read_csv(paste0(PARAM$folder.data, "/bootstrap_Harrell_c_baseline3.csv"))[,-1]
bootstrap_hoslem <- read_csv(paste0(PARAM$folder.data,"/bootstrap_hoslem_baseline3.csv"))[,-1]
Final_leaderboard <- as.data.frame(read_csv(paste0(PARAM$folder.data, "/Final_leaderboard.csv")))

# prepare dataset
Final_leaderboard$`Teams/Participants` <- c("Metformin-121", "YG-HZ team", 
                                            "DFH","TristanF","SB2", "Pteam", 
                                            "UTKteam","Baseline Age-Sex",
                                            "Baseline Covariates", "Baseline All")
Final_leaderboard$harrell_c.r <- round(Final_leaderboard$harrell_c, digits = 3)
Final_leaderboard$hoslem_sci<- scientific_format(digits=2)(Final_leaderboard$hoslem_test)

df.harrelc=as.data.frame(t(bootstrap_Harrell_c))
df.harrelc$Teams_Participants <- rownames(df.harrelc)
df.harrelc=melt(df.harrelc)
colnames(df.harrelc)  <- c("Teams/Participants", 
                           "variable", 
                           "Bootstrapped Harrell's C Index (N=1000)")
df.harrelc$`Teams/Participants` <- factor(df.harrelc$`Teams/Participants`, levels = limits)

df.hoslem = as.data.frame(t(bootstrap_hoslem))
df.hoslem$Teams_Participants <- rownames(df.hoslem)
df.hoslem=melt(df.hoslem)
colnames(df.hoslem)  <- c("Teams/Participants", 
                          "variable", 
                          "Bootstrapped Hosmer-Lemeshow Test (N=1000)")

df.hoslem$`Teams/Participants` <- factor(df.hoslem$`Teams/Participants`, levels = limits)


#################################################################################
# plot hoslem and harrelc bootstrep
pharrelc <- ggplot(df.harrelc, aes(x=`Teams/Participants`, 
                                   y=`Bootstrapped Harrell's C Index (N=1000)`,  
                                   color=`Teams/Participants`,
                                   fill=`Teams/Participants`)) +
  geom_boxplot(color="grey50")   +
  #scale_color_manual(values=values_color) +
  scale_fill_manual(values=c("#999999", "#999999","#999999", "#999999", "#999999",
                             "#9970ab", "#9970ab", "#9970ab","#fb9a99" ,"#73baae")) +
  theme_bw() +  
  theme(axis.text.x = element_text(color = "grey20", size = 12, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 12,  hjust = 1, vjust = 0, face = "plain"),
        axis.title.x = element_text(color = "grey20", size = 12, hjust = .5, vjust = 0, face = "bold"),
        axis.title.y = element_text(color = "grey20", size = 12, hjust = .5, vjust = .5, face = "bold"),
        legend.position="none") +    
  scale_x_discrete(labels=as.factor(labels), limits=as.factor(limits)) +
  #scale_x_discrete(limits=limits, labels=labels) +
  scale_y_continuous(limits=c(0, 1)) +
  #scale_y_continuous(limits=c(0, 1), breaks=c(0, 0.1, 0.2, 0.3, 0.7, 0.8, 0.9, 1)) +
  #scale_y_break(c(0.35, 0.7)) + 
  #scale_y_break(c(0.4, 0.7), scales="free") + 
  # ylim(0.2,0.9) + 
  coord_flip() 
pharrelc
#################################################################################

phoslem <- ggplot(df.hoslem, aes(x=as.factor(`Teams/Participants`), 
                                 y=`Bootstrapped Hosmer-Lemeshow Test (N=1000)`,  
                                 color=`Teams/Participants`,
                                 fill=`Teams/Participants`)) +
  geom_boxplot(color="grey50")   +
  scale_fill_manual(values=c("#999999", "#999999","#999999", "#999999", "#999999",
                             "#9970ab", "#9970ab", "#9970ab","#fb9a99", "#73baae")) +
  theme_bw() + 
  rremove("x.grid") +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x = element_text(color = "grey20", size = 12, hjust = .5, vjust = .5, face = "plain"),
        # axis.text.y = element_text(color = "grey20", size = 12,  hjust = 1, vjust = 0, face = "plain"),
        axis.title.x = element_text(color = "grey20", size =12, hjust = .5, vjust = 0, face = "bold"),
        # axis.title.y = element_text(color = "grey20", size =12, hjust = .5, vjust = .5, face = "bold"),
        legend.position="none") +  scale_y_log10() +
  
  labs(y ="Bootstrapped Hosmer-Lemeshow Test (N=1000)") + 
  
  coord_flip() 
phoslem 

# combine plots
figure2 <- (pharrelc.bar | phoslem.bar) / (pharrelc| phoslem) + 
  plot_annotation(tag_levels = 'A')
figure2
ggsave(figure2, filename = paste0(PARAM$folder.result, 
                                  "Figure 2 Overall model performance.pdf"),
       width = 13, height = 11)


#################################################################################
# Figure 5: pair comparisons
#################################################################################

pair_model <- as.data.frame(read_csv(paste0(PARAM$folder.data, "/pair_average_eval.csv"))) #updated
pair_model$harrell_c.r <- round(pair_model$`Harrell's C index`, digits = 3)
pair_model=as_tibble(pair_model)
#pair_model$`Hosmer lemeshow`[pair_model$`Hosmer lemeshow`< 0.0005] <- "< e-04"
pair_model$`Hosmer-Lemeshow` = round(pair_model$`Hosmer-Lemeshow`, digits=4)

# add matches
list_name <-list(c("SB2", "DFH"), 
                 c("SB2", "DFH","TristanF"),
                 c("SB2", "DFH", "YG-HZ team"),
                 c("SB2", "DFH", "Metformin-121"), 
                 c("SB2", "DFH", "YG-HZ team", "TristanF"),
                 c("SB2", "DFH", "YG-HZ team", "Metformin-121", "TristanF"),
                 c("SB2", "DFH", "YG-HZ team", "Metformin-121", "TristanF", "UTKteam"),
                 c("SB2", "DFH", "YG-HZ team", "Metformin-121", "TristanF", "UTKteam", "PTeam"),
                 c("SB2"), c("DFH"))

pair_model$Teams <- list_name
values=c("#73baae", "#fb9a99", "#beaed4", "#beaed4", "#beaed4", 
         "#beaed4", "#beaed4", "#beaed4", "#beaed4")

# plot
p.pairs1 <- pair_model %>%
  distinct(`Harrell's C index`, .keep_all=TRUE) %>%
  arrange(order) %>%
  ggplot(aes(x=Teams, y=`Harrell's C index`, color=`Team Pairs`)) +
  geom_path(group=1, color="grey20", alpha=0.7) +
  geom_point(alpha=1, color=values, size=3)+
  scale_x_upset(n_intersections = 20, order_by = "degree") +
  labs(y="Harrell's C Index", x="Combinations of ensemble models") + 
  ylim(0.8, 0.9)+
  theme_bw() +
  theme(axis.text.x = element_text(color = "grey20", size = 10, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 10,  face = "plain"),
        axis.title.x = element_text(color = "grey20", size = 10, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 10, face = "plain")) +  
  theme_combmatrix(combmatrix.panel.point.color.fill = "grey30",
                   combmatrix.panel.line.size = 0.5,
                   combmatrix.label.make_space = FALSE,
                   combmatrix.panel.point.size = 3) +
  theme(plot.margin = margin(1,1,1,1, "cm")) + ggtitle("A")

p.pairs2 <- pair_model %>%
  distinct(`Harrell's C index`, .keep_all=TRUE) %>%
  arrange(order) %>%
  ggplot(aes(x=Teams, y=log2(`Hosmer-Lemeshow`), color=`Team Pairs`)) +
  geom_path(group=1, color="grey20", alpha=0.7)+
  geom_point(alpha=1, color=values, size=3) +
  geom_hline(yintercept=log2(0.05), linetype="dashed") +
  scale_x_upset(n_intersections = 20, order_by = "degree") +
  labs(y="log2(p-value)", x="Combinations of ensemble models") +
  #scale_y_break(c(1e-06,1e-30),scales = "free") + 
  theme_bw() +
  theme(axis.text.x = element_text(color = "grey20", size = 10, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 10,  face = "plain"),
        axis.title.x = element_text(color = "grey20", size = 10, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 10, face = "plain"),
        legend.position="right") +  
  theme(plot.margin = margin(1,1,1,1, "cm")) + ggtitle("B") +
  theme_combmatrix(combmatrix.panel.point.color.fill = "grey30",
                   combmatrix.panel.line.size = 0.5,
                   combmatrix.label.make_space = FALSE,
                   combmatrix.panel.point.size = 3,
                   combmatrix.label.text = element_blank())

p.pairs2
#combine
figure5 <- p.pairs1 + p.pairs2 +
  plot_layout()
#save
ggsave(figure5, filename = paste0(PARAM$folder.result, 
                                  "Figure 5 Aggregated models.pdf"),
       width = 8, height = 4)


