#### Rarefaction ####
library(vegan)
library(reshape2)
library(tidyr)
library(ggplot2)
library(dplyr)

ps <- readRDS('./Output/ps_itsxpress_assigned.RDS')
ps_fungi <- ps %>% subset_taxa(Kingdom == 'k__Fungi') %>% subset_samples(Site != 'NA')

# Selecting rarefaction threshold
fungi_otu <- otu_table(ps_fungi) %>% as.data.frame()

fungi_otu$n_seqs <-  fungi_otu %>% rowSums()
fungi_otu <- arrange(fungi_otu, n_seqs)


### rarefaction curves
rare_curve <- rarecurve(fungi_otu %>% select(-n_seqs), step = 100, label = F, tidy = TRUE)
levels = rare_curve %>% group_by(Site) %>%
  summarise(Max = max(Species)) %>%
  arrange(desc(Max)) %>%
  .$Site

rare_curve$Site <- factor(rare_curve$Site, levels = levels)

ggplot(data = rare_curve, aes(x = Sample, y = Species, colour = Site)) + 
  geom_line(linewidth = 1) +
  theme_bw() +
  theme(legend.position = 'none')


ggsave('./Plots/rarecurve_all.pdf')

# for bottom 20 samples with the least reads
rare_curve20 <- rarecurve(fungi_otu[1:5,] %>% select(-n_seqs), step = 20, label = F, tidy = TRUE)
levels = rare_curve20 %>% group_by(Site) %>%
  summarise(Max = max(Species)) %>%
  arrange(desc(Max)) %>%
  .$Site

rare_curve20$Site <- factor(rare_curve20$Site, levels = levels)

ggplot(data = rare_curve20, aes(x = Sample, y = Species, colour = Site)) + 
  geom_line(linewidth = 1) +
  theme_bw()


ggsave('./Plots/rarecurve_bottom20.pdf')

### visualizing reads distribution
fungi_otu %>% ggplot(aes(x = n_seqs)) +
  geom_histogram(binwidth = 200) +
  coord_cartesian(xlim=c(0,5000))

meta <- microbiome::meta(ps_fungi)
fungi_otu$SampleType = meta$Compartment

# jittering number of sequences to identify potential outliers
fungi_otu %>%
  ggplot(aes(x=1, y=n_seqs, colour = SampleType)) +
  geom_jitter(size = 2) +
  scale_y_log10() +
  theme_bw() +
  labs(y = 'Number of sequences') +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

ggsave('./Plots/number_of_seqs_jitter.pdf')

# identifying samples with significantly low number of reads
fungi_otu %>%
  arrange(n_seqs) %>%
  ggplot(aes(x=1:nrow(.), y = n_seqs)) +
  geom_line() +
  coord_cartesian(xlim = c(0,500), ylim=c(0,10000))

fungi_otu %>%
  arrange(n_seqs) %>%
  select(n_seqs) # minimum number of reads = 598

# Investigate goods coverage based on different thresholds
fungi_otu$Sample <- rownames(fungi_otu) 
goods_stats <- fungi_otu %>%
  select(-c(n_seqs,SampleType)) %>%
  pivot_longer(-Sample) %>%
  group_by(Sample) %>%
  summarize(n_seqs = sum(value),
            n_sing = sum(value == 2),
            Goods = 1- n_sing/n_seqs) %>%
  filter(n_seqs > 0)

goods_stats$SampleType = fungi_otu$SampleType

ggplot(goods_stats, aes(x = n_seqs, y = Goods, colour = SampleType)) + 
  geom_point(size = 2) +
  coord_cartesian(xlim=c(0,5000), ylim=c(0,1)) +
  coord_cartesian(ylim=c(0.98,1)) +
  theme_bw()

ggsave('./Plots/goods_coverage_doubletons.pdf')

# Investigate if there are any significant richness vs number of sequences trend
richness_reads <- fungi_otu %>%
  select(-c(n_seqs,SampleType)) %>%
  pivot_longer(-Sample) %>% 
  group_by(Sample) %>% 
  summarize(richness = sum(value > 0),
            n_seqs = sum(value)) 

richness_reads$SampleType <- fungi_otu$SampleType

ggplot(richness_reads, aes(x = n_seqs, y = richness)) + 
  geom_point(aes(colour = SampleType), size = 2) +
  geom_line(stat = 'smooth') +
  theme_bw() +
  coord_cartesian(xlim=c(0,50000)) +
  labs(x='Number of sequences', y='ASV richness')

ggsave('./Plots/richness_w_sequences.pdf')


# Perform rarefaction
fungi_otu <- otu_table(ps_fungi) %>% as.data.frame()

rarerified_otu <- matrix(0, nrow=nrow(fungi_otu), ncol=ncol(fungi_otu))
rownames(rarerified_otu) <- rownames(fungi_otu)
colnames(rarerified_otu) <- colnames(fungi_otu)

for (i in 1:999){
  temp <- rrarefy(fungi_otu, sample = fungi_otu %>% rowSums() %>% min())
  rarerified_otu = rarerified_otu + temp
  print(paste("Iteration:", i))
}
beepr::beep()

final_rarerified_otu = rarerified_otu/999
rownames(final_rarerified_otu) <- rownames(rarerified_otu)
colnames(final_rarerified_otu) <- colnames(rarerified_otu)

final_rarerified_otu %>% rowSums()

ps_fungi_rare <- ps_fungi
ps_fungi_rare@otu_table <- otu_table(final_rarerified_otu, taxa_are_rows = F)

saveRDS(ps_fungi_rare, './Output/ps_fungi_rare.RDS')
