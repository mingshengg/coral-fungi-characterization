# Load packages ####
library(ggplot2)
library(phyloseq)
library(tidyverse)
library(Biostrings)
library(vegan)
library(speedyseq)
library(emmeans)

# Standard error formula #
se <- function (x) {
  sd(x)/sqrt(length(x))
}

# Read dataframes
ps <- readRDS('./Output/ps_itsxpress_assigned.RDS')
ps_coral_skeleton <- ps %>% subset_samples(Compartment == 'Skeleton')
ps_coral_skeleton@sam_data


##### removing contamination from PS1_JO_S_F1 as there are high number of random ASVs found from mangroves only exclusively found here #####
cont <- ps %>% subset_taxa(Kingdom == 'k__Fungi') %>% subset_samples(SampleID=='PS1_JO_S_F1')
ps_wo_cont <- ps %>% subset_taxa(Kingdom == 'k__Fungi') %>% subset_samples(SampleID != 'PS1_JO_S_F1')

A <- cont %>% prune_taxa(taxa_sums(.)>0, .) %>% tax_table() %>% row.names()
B <- ps_wo_cont %>% prune_taxa(taxa_sums(.)>0, .) %>% tax_table() %>% row.names()

contam = setdiff(A,B)

ps_primer <- ps_coral_skeleton %>% subset_taxa(!taxa_names(ps_coral_skeleton) %in% contam)

# Rarefaction curves; Overall
otu_tab <- otu_table(ps_primer %>% subset_taxa(Kingdom == 'k__Fungi')) %>% data.frame()
rare_curve <- rarecurve(otu_tab, step = 30, label = F, tidy = TRUE)
levels = rare_curve %>% group_by(Site) %>%
  summarise(Max = max(Species)) %>%
  arrange(desc(Max)) %>%
  .$Site

rare_curve$Site <- factor(rare_curve$Site, levels = levels)
rare_curve$Coral = substr(rare_curve$Site, 1, 2)
rare_curve$Primer = substr(rare_curve$Site, 10, 11)
rare_curve$Primer_val = as.numeric(substr(rare_curve$Site, 11, 11))

ggplot(data = rare_curve, aes(x = Sample, y = Species, colour = Site)) + 
  geom_line(aes(size = Primer_val)) +
  theme_bw()

ggsave('./Plots/Coral_skeleton/skeleton_rarefaction.pdf')

# Rarefaction curve; lowest 5 samples
lowest5 <- otu_tab %>% rowSums() %>% sort() %>% names() %>% .[1:5]
rare_curve_l5 <- rarecurve(otu_tab[lowest5,], step = 5, label = F, tidy = TRUE)

rare_curve_l5$Site <- factor(rare_curve_l5$Site, levels = levels)
rare_curve_l5$Coral = substr(rare_curve_l5$Site, 1, 2)
rare_curve_l5$Primer = substr(rare_curve_l5$Site, 10, 11)

ggplot(data = rare_curve_l5, aes(x = Sample, y = Species, colour = Site)) + 
  geom_line(linewidth = 1, aes(linetype = Primer)) +
  theme_bw() + 
  scale_y_continuous(limits = c(0, 20), breaks = seq(0, 20, 2)) +
  scale_x_continuous(limits= c(0,1500))

ggsave('./Plots/Coral_skeleton/skeleton_rarefaction_zoomed.pdf')

#### re-do rarefaction after removal of contamination ####
fungi_otu <- otu_table(ps_primer %>% subset_taxa(Kingdom == 'k__Fungi')) %>% as.data.frame()

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

ps_fungi_rare <- ps_primer
ps_fungi_rare@otu_table <- otu_table(final_rarerified_otu, taxa_are_rows = F)

saveRDS(ps_fungi_rare, './Output/ps_itsxpress_assigned_skeleton_rare.RDS')

# 1. investigate whether percentage fungal reads were significantly affected by primer ####
ps_reads <- ps_primer
meta_reads <- microbiome::meta(ps_reads)

meta_reads$Fungi_reads <- otu_table(ps_reads %>% subset_taxa(Kingdom == 'k__Fungi')) %>% rowSums()
meta_reads$Primer_1 <- ifelse(meta_reads$Primer_set == 'F1'|meta_reads$Primer_set == 'F2',
                              'F1/F2',
                              'F3')

# add fungal reads to meta
meta_reads$Percentage <- ps_reads %>% microbiome::transform('compositional') %>% psmelt() %>% filter(Kingdom == 'k__Fungi') %>% group_by(Sample) %>%
  summarise(temp = sum(Abundance)*100) %>% as.data.frame() %>% .[,2]

mod1 <- lm(Ratio ~ Primer*Species, data = meta_reads)
shapiro.test(resid(mod1))
plot(mod1)
summary(mod1)
anova(mod1)

# pairwise testing
emmeans(mod1, specs=pairwise~Primer|Species)
write.csv(emmeans(mod1, specs=pairwise~Primer|Species)$contrast, './Output/Coral_skeleton/fungi_percentage_emmeans.csv')

# Reads averaged across primer
meta_reads %>% 
  group_by(Species, Primer_1) %>%
  summarise(mean = mean(Anthozoa_Percentage),
            se = se(Anthozoa_Percentage))

# Plot % fungal reads
meta_reads %>% 
  ggplot(data = ., aes(x=Primer, y=Percentage, fill = Primer)) +
  stat_summary(geom = "bar", fun = mean, position = "dodge", color="black") +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = 'dodge', width=.2, size = 1) +
  ylim(c(0,100)) + 
  facet_wrap(~Species) +
  theme_bw()

ggsave('./Plots/Coral_skeleton/fungi_percentage.pdf')

# 2. compare alpha-diversity indices ####
ps_primer_fungi_only <- readRDS('./Output/ps_itsxpress_assigned_skeleton_rare.RDS')
ps_primer_fungi_only <- prune_taxa(taxa_sums(ps_primer_fungi_only) > 0, ps_primer_fungi_only)

# Extract metadata and calculate alpha-diversity indices
meta <- microbiome::meta(ps_primer_fungi_only)                
meta$Shannon <- vegan::diversity(otu_table(ps_primer_fungi_only),index = "shannon")
meta$Richness <- vegan::specnumber(otu_table(ps_primer_fungi_only))
meta$Evenness <- microbiome::evenness(otu_table(ps_primer_fungi_only), index = "pielou") %>% .$pielou

meta$Primer_1 <- ifelse(meta$Primer_set == 'F1'|meta$Primer_set == 'F2',
                        'F1/F2',
                        'F3')

# add to ps object
ps_primer_fungi_only@sam_data$Richness <- meta$Richness
ps_primer_fungi_only@sam_data$Shannon <- meta$Shannon
ps_primer_fungi_only@sam_data$Evenness <- meta$Evenness

# Shannon ANOVA model
shannon_mod = lm((Shannon) ~ Primer*Species, data = meta)
plot(shannon_mod)
shapiro.test(resid(shannon_mod))
summary(shannon_mod)
anova(shannon_mod)

# Richness ANOVA model
richness_mod = lm(log(Richness) ~ Primer_1*Species, data = meta)
plot(richness_mod)
shapiro.test(resid(richness_mod))
summary(richness_mod)
anova(richness_mod)

emmeans(richness_mod, specs = pairwise~Primer_1|Species)
write.csv(emmeans(richness_mod, specs = pairwise~Primer|Species)$contrast,
          './Output/Coral_skeleton/emmeans_skeleton_richness.csv')

# Evenness ANOVA model
evenness_mod = lm((Evenness) ~ Primer*Species, data = meta)
plot(evenness_mod)
shapiro.test(resid(evenness_mod))
summary(evenness_mod)
anova(evenness_mod)

### alpha-diversity overall
meta %>% group_by(Primer) %>% 
  summarise(richness_mean = mean(Richness), richness_se = se(Richness),
            shannon_mean = mean(Shannon), shannon_se = se(Shannon),
            evenness_mean = mean(Evenness), evenness_se = se(Evenness))

meta %>%
  summarise(richness_mean = mean(Richness), richness_se = se(Richness),
            shannon_mean = mean(Shannon), shannon_se = se(Shannon),
            evenness_mean = mean(Evenness), evenness_se = se(Evenness))

meta %>% group_by(Primer) %>%
  summarise(richness_mean = mean(Richness), richness_se = se(Richness),
            shannon_mean = mean(Shannon), shannon_se = se(Shannon),
            evenness_mean = mean(Evenness), evenness_se = se(Evenness))


# 3. compositional differences across methods ####
# calculate bray-curtis dissimilarity
mdf_com <- ps_primer_fungi_only  %>% otu_table() %>% vegdist()

# NMDS
otu_res <- metaMDS(mdf_com, center = T)

data.scores = as.data.frame(scores(otu_res, display =  'sites'))
data.scores$Primer = meta$Primer
data.scores$Species = meta$Species

data.scores <- data.scores[order(rownames(data.scores)),]

gg = ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(data = data.scores, aes(colour = Primer, shape = Species), size = 3, alpha = 1) +
  geom_polygon(stat = "ellipse", aes(fill = Primer), alpha = 0.3) +
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
        legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 9, colour = "grey30")) + 
  labs(title = 'Stress = 0.1566') +
  xlab('NMDS1') +
  ylab('NMDS2')


ggsave('./Plots/Coral_skeleton/skeleton_NMDS_primer.pdf')

# PERMANOVA to investigate difference in composition across species*primer
comp_mod <- adonis2(mdf_com ~ Primer*Species, data = meta, by = 'terms')
pairwise.adonis(mdf_com, factors = meta$Primer)

write.csv(pairwise.adonis(mdf_com, factors = meta$Primer), './Output/Coral_skeleton/permanova_primer_pairwise.csv')

# betadisper to assess homogeneity
bd <- betadisper(vegdist(mdf_com), meta$Primer)
plot(bd)
boxplot(bd)
anova(bd)
permutest(bd, pairwise = T)

# Compositional plots
ps_primer_fungi_only <- prune_taxa(taxa_sums(ps_primer_fungi_only)>0, ps_primer_fungi_only)
ps_merge <- ps_primer_fungi_only %>% merge_samples2('Primer', fun_otu = mean, funs = list(mean))
ps_merge_comp <- ps_merge %>% microbiome::transform('compositional')

compositional_plot <- ps_merge_comp %>%
  psmelt() %>%
  ggplot(aes(y = Abundance, x = Sample, fill = Class)) + #change fill = Taxonomic rank
  geom_bar(stat = 'identity') +
  theme_minimal()

# venn diagramm ####
ps_venn <- ps_primer_fungi_only %>% microbiome::transform('compositional')

A <- ps_venn %>% subset_samples(Primer == "fITS7+ITS4") %>%
  subset_taxa(taxa_sums(ps_venn %>% subset_samples(Primer == "fITS7+ITS4"))>0) %>% tax_table() %>% row.names()
B <- ps_venn %>% subset_samples(Primer == "ITS86+ITS4") %>%
  subset_taxa(taxa_sums(ps_venn %>% subset_samples(Primer == "ITS86+ITS4"))>0) %>% tax_table() %>% row.names()
C <- ps_venn %>% subset_samples(Primer == "MS7+ITS4") %>%
  subset_taxa(taxa_sums(ps_venn %>% subset_samples(Primer == "MS7+ITS4"))>0) %>% tax_table() %>% row.names()

full <- unique(c(A,B,C))
n12 <- sum(full %in% unique(A) & full %in% unique(B))
n13 <- sum(full %in% unique(A) & full %in% unique(C))
n23 <- sum(full %in% unique(B) & full %in% unique(C))

n123 <- sum(full %in% unique(A) & full %in% unique(B) & full %in% unique(C))

venn.plot1 <- VennDiagram::draw.triple.venn(length((A)), length((B)), length((C)),
                                            n12, n23, n13,
                                            n123,
                                            category = c('fITS7+ITS4', 'ITS86+ITS4', 'MS7+ITS4'))

dev.off()

# ASVs that are shared across the three primer pairs
shared = A[full %in% unique(A) & full %in% unique(B) & full %in% unique(C)]
ps_venn@otu_table[,shared] %>% as.data.frame() %>% rowSums()
ps_venn@tax_table[shared] %>% as.data.frame() %>% dplyr::select(Family)

# Adding % of shared ASVs to metadata
meta$Shared <- ps_venn@otu_table[,shared] %>% as.data.frame() %>% rowSums()

# Assessing if % of shared ASVs are statistically different across different primers
mod_shared <- lm((Shared) ~ Primer*Species, data = meta)
plot(mod_shared)
shapiro.test(resid(mod_shared))
summary(mod_shared)
anova(mod_shared)

meta %>% summarise(mean_perc = mean(Shared), se_perc = se(Shared))

# composition of shared
my_subset <- subset(tax_table(ps_venn), rownames(tax_table(ps_venn)) %in% shared)
ps_shared <- merge_phyloseq(my_subset, otu_table(ps_venn), sample_data(ps_venn))

ps_shared_comp <- ps_shared %>% speedyseq::merge_samples2('Primer', fun_otu = mean, funs = list(mean)) %>%
  microbiome::transform('compositional')

shared_compositional_plot <- ps_shared_comp %>%
  psmelt() %>%
  ggplot(aes(y = Abundance, x = Primer, fill = Family)) + #change fill = Taxonomic rank
  geom_bar(stat = 'identity') +
  theme_minimal()
