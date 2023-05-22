cran_packages <- c("tidyverse", "BiocManager", "knitr", "ggrepel", "pROC", "vegan",
                   "reshape2", "ggplot2", "ggpubr", "car",  "dplyr", "plyr")
for (i in cran_packages) {
  if (!require(i, character.only = TRUE))
    install.packages(i)
}
#load library
library(knitr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(car)
library(vegan)
​
# load file
folder <- gsub("/scp", "", getwd())
folder.results <- paste0(folder, "/results/")
file.path(folder, 'data/mobi.Rdata')
​
# load data
load(file=file.path(folder, 'data/mobi.Rdata'))
meta$ID <- rownames(meta)
​
# motu.abs is the absolue counts of species
# head(motus.abs)
# MMPC35551931ST MMPC41376408ST
#'Candidatus Kapabacteria' thiocyanatum [r_10354]              0              0
#Abiotrophia defectiva [r_04788]                               0              0
​
# head meta
# environment_material status
#MMPC35551931ST	feces [ENVO:00002003]	PC	1	0	
#MMPC41376408ST	feces [ENVO:00002003]	PC	1	0	
#MMPC59659730ST	feces [ENVO:00002003]
​
# apply filter and check if we loose so muc, we can adjust
motu.abs.fil <- motu.abs[rowSums(motu.abs >= 10^-5) >= 2,colSums(motu.abs) > 1200 ]
​
dim(motu.abs)
dim(motu.abs.fil)
​
# calculate bray curtis
motu.t = t(motu.abs.fil)
min_seq = min(rowSums(motu.t))
​
beta_dist.rar <-t(motu.abs.fil) %>%
  vegan::avgdist(dmethod = "bray", sample = min_seq)
nmds.rar <- metaMDS(beta_dist.rar) %>% scores(display=c("sites")) %>% as.tibble(rownames="ID")
​
# combine metadata and betadiv
meta_nmds.rar <- dplyr::inner_join(meta, nmds.rar)
​
p.betadiv.ord.rar <- ggplot(meta_nmds.rar, 
                            aes(x = NMDS1, 
                                y = NMDS2, 
                                color = status)) + #change color here accordingly such as age groups, gender
  geom_point() +
  stat_ellipse() +
  scale_colour_manual(values = c("#e41a1c",  # add colors according to number of categories
                                # "#984ea3", 
                                 #"#377eb8",
                                 #"#4daf4a",
                                 "#ff7f00")) +
  theme_classic()
​
p.betadiv.ord.rar
ggsave(p.betadiv.ord.rar, file=paste0(folder.results,
                                      "betadiv_nmds_ordination.pdf"), 
       width = 5, height=4)
​
​
# prepare pcoa
dist_tbl <- as_tibble(as.matrix(beta_dist.rar), 
                      rownames="samples")
​
dist_matrix <- dist_tbl %>%
  pivot_longer(cols=-samples, names_to="b", values_to="distances") %>%
  select(samples, b, distances) %>%
  pivot_wider(names_from="b", values_from="distances") %>%
  select(-samples) %>%
  as.dist()
​
pcoa <- cmdscale(dist_matrix, eig=TRUE, add=TRUE)
positions <- pcoa$points
colnames(positions) <- c("pcoa1", "pcoa2")
percent_explained <- 100 * pcoa$eig / sum(pcoa$eig)
​
pretty_pe <- format(round(percent_explained, digits =1), 
                    nsmall=1, 
                    trim=TRUE)
​
labels <- c(glue("PCoA 1 ({pretty_pe[1]}%)"),
            glue("PCoA 2 ({pretty_pe[2]}%)"))
​
# plot percent explained
tibble(pe = percent_explained,
       axis = 1:length(percent_explained)) %>%
  ggplot(aes(x=axis, y=pe)) +
  geom_line() +
  coord_cartesian(xlim = c(1, 10), ylim=c(0, 10)) +
  scale_x_continuous(breaks=1:10) +
  theme_classic()
​
# combine metadata and betadiv
positions <- positions %>% 
  as.data.frame() %>% 
  rownames_to_column( var = "ID")
​
pcoa.rar <- dplyr::inner_join(meta, p)
​
#please adapt the code to finrisk meta covariates
#plot pcoa
p.pcoa <- pcoa.rar %>%
ggplot(aes(x=pcoa1, y=pcoa2, 
           color= as.factor(status), #replace with Event
           size= stage,#replace with Age group
           shape = as.factor(gender)))+ 
  geom_point(alpha=0.5) + 
  labs(x=labels[1], y=labels[2]) +
  theme_classic() + 
  #scale_shape_manual(values=c(17, 1), name="Gender") + #keep
  scale_size_continuous(name="Stage") + #change to age group
  scale_colour_manual(values = c("#e41a1c",  # add colors according to number of categories 
                                 #"#ff7f00", 
                                 #"#984ea3",
                                 #"#4daf4a",
                                 "#377eb8"))+
  ggtitle("Bray Curtis")  +
  stat_ellipse(inherit.aes = F, #in your example should be able to just use stat_ellipse()
               aes(color=status, x=pcoa1, y=pcoa2)) 
​
p.pcoa
ggsave(p.pcoa, file=paste0(folder.results,
                                      "betadiv_pcoa.pdf"), 
       width = 5, height=5)