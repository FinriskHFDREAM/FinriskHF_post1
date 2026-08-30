Suppfigure1B <- df.counts %>%
  mutate(value=factor(value, levels = c( "No","Yes", "Not available"))) %>%
  #mutate(Group=factor(Group, levels = c("Train", "Test", "Score"))) %>%
  mutate(data=factor(data, levels = c("FINRISK", "Synthetic"))) %>%
  
  mutate(Groupset = paste0(data, " ", Group)) %>%
  mutate(Groupset=factor(Groupset, 
                         levels = c("FINRISK Train", "FINRISK Test", "FINRISK Score", 
                                    "Synthetic Score", "Synthetic Test", "Synthetic Train"))) %>%
  group_by(variable, Groupset) %>% #do calculations by siteID
  mutate(percent = n / sum(n) * 100) %>%
  ggplot() +
  geom_col_pattern(aes(x = as.factor(Groupset),
                       y = percent,
                       pattern =  value,
                       fill = as.factor(Groupset)),
                   position = position_dodge(width=0.8),
                   alpha=0.8,
                   colour  = 'black',
                   pattern_density=0.08,
                   pattern_fill="black",
                   pattern_colour="black") +  
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  ylim(0, 100)+
  facet_wrap(. ~ variable,  ncol=4) + 
  scale_color_manual(values=values) +
  scale_pattern_manual(
    values = c(
      "No" = "none",
      "Yes" = "stripe",
      "Not available" = "crosshatch"
    )
  ) +
  scale_fill_manual(values=values) +
  
  labs(y = "Percent of individuals (%)", 
       x="", 
       pattern = "Condition",
       fill = "Dataset") +
  theme_classic() +  
  theme(panel.border=element_blank(), 
        axis.line=element_line(),
        legend.position=c(0.9, 0.2)) +  ggtitle("B") +
  guides(
    pattern = guide_legend(
      override.aes = list(
        fill = "white",
        colour = "black",
        pattern_fill = "black",
        pattern_density = 0.15,
        pattern_spacing = 0.03
      )
    )
  )
