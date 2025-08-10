# Load packages ####
library(ggplot2)
library(phyloseq)
library(tidyverse)
library(dplyr)
library(Biostrings)
library(vegan)
library(speedyseq)
library(emmeans)
library(pairwiseAdonis)

# Standard error formula #
se <- function(x){
  sd(x)/sqrt(length(x))
}

# Read dataframes
raw_reads <- readxl::read_xlsx('./raw_reads.csv') # dataframe of reads
ps_raw <- readRDS('./Output/ps_itsxpress_assigned.RDS') # overall phyloseq object
ps_raw <- ps_raw %>% subset_samples(Compartment == 'Tissue')


# Rarefaction curve; Overall
otu_tab <- otu_table(ps_raw %>% subset_taxa(Kingdom == 'k__Fungi')) %>% data.frame()
rare_curve <- rarecurve(otu_tab, step = 5, label = F, tidy = TRUE)
levels = rare_curve %>% group_by(Site) %>%
  summarise(Max = max(Species)) %>%
  arrange(desc(Max)) %>%
  .$Site

rare_curve$Site <- factor(rare_curve$Site, levels = levels)
rare_curve$Coral = substr(rare_curve$Site, 1, 2)
rare_curve$Primer = substr(rare_curve$Site, 10, 11)

ggplot(data = rare_curve, aes(x = (Sample), y = Species, colour = Site)) + 
  geom_line() +
  theme_bw() +
  theme(legend.position = "none")

ggsave('./Plots/Coral_tissue/tissue_rarefaction_lowest5.pdf')

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

ggsave('./Plots/tissue_rarefaction_zoomed.pdf')

# 1. investigate whether percentage fungal reads were significantly affected by dna extraction method + primer ####
raw_reads$Method <- paste0(raw_reads$DNA_extraction, raw_reads$Primer)
mod1 <- lm(log(Fungi_percentage) ~ DNA_extraction*Primer*Species, data = raw_reads)
shapiro.test(resid(mod1))
plot(mod1)
mod1
summary(mod1)
anova(mod1)
write.csv(anova(mod1), './Output/Coral_tissue/fungi_percentage_aov.csv')

# pairwise testing
emmeans_primer <- emmeans(mod1, specs=pairwise~Primer)
emmeans_pair <- emmeans(mod1, specs=pairwise~DNA_extraction:Primer|Species)
write.csv(emmeans_pair$contrasts, './Output/Coral_tissue/fungi_percentage_emmeans.csv')

# Reads averaged across primer and species
raw_reads %>%  group_by(Primer, Species) %>% 
  summarise(Fungal_Mean = mean(Fungi_reads), Fungal_se = se(Fungi_reads),
            Percentage_Mean = mean(Fungi_percentage), Percentage_se = se(Fungi_percentage))

# Reads averaged across primer and DNA extraction
raw_reads %>%  group_by(Primer, DNA_extraction) %>% 
  summarise(Fungal_Mean = mean(Fungi_reads), Fungal_se = se(Fungi_reads),
            Percentage_Mean = mean(Fungi_percentage), Percentage_se = se(Fungi_percentage))


# plot % fungal reads graph
# combining F1 and F2 since there are no statistical difference in generated % fungal reads
raw_reads$Primer_1 <- ifelse(raw_reads$Primer == 'F1'|raw_reads$Primer == 'F2', 'F1/F2', raw_reads$Primer)
raw_reads$Method_1 <- paste0(raw_reads$DNA_extraction, raw_reads$Primer_1)

raw_reads %>% 
  ggplot(data = ., aes(x=DNA_extraction, y=Fungi_percentage, fill = Primer_1)) +
  stat_summary(geom = "bar", fun = mean, position = "dodge", color="black") +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = 'dodge', linewidth = 1) +
  ylim(c(0,100)) + 
  facet_wrap(~Species) +
  theme_bw()

ggsave('./Plots/Coral_tissue/fungi_percentage.pdf')

# Reads averaged across library construction methods
raw_reads %>%  group_by(Method_1) %>% 
  summarise(Fungal_Mean = mean(Fungi_reads), Fungal_se = se(Fungi_reads),
            Percentage_Mean = mean(Fungi_percentage), Percentage_se = se(Fungi_percentage))

# 2. compare alpha-diversity indices ####
ps <- readRDS('./Output/ps_fungi_rare.RDS')

# subset to tissue samples
ps_tissue <- ps %>% subset_samples(Compartment  == 'Tissue')
ps_tissue <- prune_taxa(taxa_sums(ps_tissue) > 0, ps_tissue)

# Extract metadata and calculate alpha-diversity indices
meta <- microbiome::meta(ps_tissue)                
meta$Shannon <- vegan::diversity(otu_table(ps_tissue),index = "shannon")
meta$Richness <- vegan::specnumber(otu_table(ps_tissue))
meta$Evenness <- microbiome::evenness(otu_table(ps_tissue), index = "pielou") %>% .$pielou

#add to ps object
ps_tissue@sam_data$Richness <- meta$Richness
ps_tissue@sam_data$Shannon <- meta$Shannon
ps_tissue@sam_data$Evenness <- meta$Evenness

# Create new factor for overall library construction method
meta$Method = paste(meta$DNA_extraction, meta$Primer, sep='-')
ps_tissue@sam_data$Method = paste(meta$DNA_extraction, meta$Primer, sep='-')

# Shannon ANOVA model
shannon_mod = lm((Shannon) ~ Method*Species, data = meta)
plot(shannon_mod)
shapiro.test(resid(shannon_mod))
summary(shannon_mod)
anova(shannon_mod)

emmeans(shannon_mod, specs=pairwise~Primer|DNA_extraction)

meta %>% group_by(Species) %>% summarise(Mean = mean(Shannon),  SE = se(Shannon))

# Richness ANOVA model
richness_mod = lm(log(Richness) ~ Method*Species, data = meta)
plot(richness_mod)
shapiro.test(resid(richness_mod))
summary(richness_mod)
anova(richness_mod)

emmeans(richness_mod, specs=pairwise~Method|Species)

meta %>% filter(Species == 'Pocillopora_acuta') %>% group_by(Method) %>% summarise(Mean = mean(Richness),  SE = se(Richness))

write.csv(emmeans::emmeans(richness_mod, specs = pairwise ~ Method|Species)$contrasts,'./Output/Coral_tissue/emmeans_tissue_richness.csv')

# Evenness ANOVA model
evenness_mod = lm(Evenness ~ Method*Species, data = meta)
plot(evenness_mod)
shapiro.test(resid(evenness_mod))
summary(evenness_mod)
anova(evenness_mod)

emmeans(evenness_mod, specs=pairwise~Method|Species)

write.csv(emmeans::emmeans(evenness_mod, specs = pairwise ~ Method|Species)$contrasts,'./Output/Coral_tissue/emmeans_tissue_evenness.csv')

meta %>% filter(Species == 'Pachyseris_speciosa') %>%
  group_by(Method) %>%
  summarise(evenness_mean = mean(Evenness), evenness_se = se(Evenness))

### alpha-diversity overall
meta %>% group_by(Method) %>% 
  summarise(richness_mean = mean(Richness), richness_se = se(Richness),
            shannon_mean = mean(Shannon), shannon_se = se(Shannon),
            evenness_mean = mean(Evenness), evenness_se = se(Evenness))

meta %>% 
  summarise(richness_mean = mean(Richness), richness_se = se(Richness),
            shannon_mean = mean(Shannon), shannon_se = se(Shannon),
            evenness_mean = mean(Evenness), evenness_se = se(Evenness))

# 3. compositional differences across methods ####
# calculate bray-curtis dissimilarity
mdf_com <- ps_tissue %>% otu_table() %>% vegdist()

# NMDS
otu_res <- metaMDS(mdf_com)
otu_res

data.scores = as.data.frame(scores(otu_res, display =  'sites'))
rownames(data.scores) == rownames(meta)

data.scores$Method = meta$Method
data.scores$Primer = meta$Primer
data.scores$DNA_extraction = meta$DNA_extraction
data.scores$Species = meta$Species
data.scores$Method_1 = meta$Method_1

# Plot NMDS
gg = ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(data = data.scores, aes(colour = Method, shape = Species), size = 3, alpha = 1) +
  geom_polygon(stat = 'ellipse', aes(fill = Method), alpha = 0.3)  +
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
        legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 9, colour = "grey30")) + 
  labs(title = 'Stress = 0.1626') +
  xlab('NMDS1') +
  ylab('NMDS2')

gg

ggsave('./Plots/Coral_tissue/tissue_NMDS_Method.pdf')

# PERMANOVA to investigate difference in composition across species*method
comp_mod <- adonis2(mdf_com ~ Species*Method, data = meta, by = 'terms')

meta$Primer_1 <- ifelse(meta$Primer_set=='F1'|meta$Primer_set=='F2', 'F1/F2', 'F3') 
meta$Method_1 <- paste0(meta$DNA_extraction, meta$Primer_1)

method_pairwise <- pairwise.adonis(mdf_com, factors = meta$Method)
write.csv(method_pairwise, './Output/Coral_tissue/method_pairwise.csv')

# betadisper to assess homogeneity
bd_mod <- betadisper(vegdist(mdf_com), meta$Method)
anova(bd_mod)
boxplot(bd_mod)
permutest(bd_mod, pairwise = T)

# Compositional plots
ps_tissue@sam_data$Method = meta$Method
ps_merge <- ps_tissue %>% merge_samples2('Method', fun_otu = mean, funs = list(mean))
ps_merge_comp <- ps_merge %>% microbiome::transform('compositional')

tax_tab <- data.frame(ps_merge_comp@tax_table)
tax_tab$Class <- ifelse(str_detect(tax_tab$Class, 'Incertae_sedis'), 'NA', tax_tab$Class)
ps_merge_comp@tax_table <- tax_table(as.matrix(tax_tab))

compositional_plot <- ps_merge_comp %>%
  psmelt() %>% 
  ggplot(aes(y = Abundance, x = Sample, fill = Class)) + #change fill = Taxonomic rank
  geom_bar(stat = 'identity') +
  theme_minimal()

ggsave('./Plots/Coral_tissue/compositional_method_class.pdf')

#### venn diagram ####
ps_tissue <- readRDS('./Output/ps_fungi_rare.RDS') %>% subset_samples(Compartment  == 'Tissue')

meta <- microbiome::meta(ps_tissue)
ps_tissue@sam_data$Method = paste(meta$DNA_extraction, meta$Primer, sep='-')
meta$Method <- ps_tissue@sam_data$Method

ps_venn <- ps_tissue %>% tax_glom('Family') %>% microbiome::transform('compositional')
ps_venn <- prune_taxa(taxa_sums(ps_venn) > 0, ps_venn)

A <- ps_venn %>% subset_samples(Method == "Zymo_HostZero-fITS7+ITS4") %>% 
  subset_taxa(taxa_sums(ps_venn %>% subset_samples(Method == "Zymo_HostZero-fITS7+ITS4"))>0) %>% tax_table() %>% row.names()
B <- ps_venn %>% subset_samples(Method == "Zymo_HostZero-ITS86+ITS4") %>% 
  subset_taxa(taxa_sums(ps_venn %>% subset_samples(Method == "Zymo_HostZero-ITS86+ITS4"))>0) %>% tax_table() %>% row.names()
C <- ps_venn %>% subset_samples(Method == "Zymo_HostZero-MS7+ITS4") %>% 
  subset_taxa(taxa_sums(ps_venn %>% subset_samples(Method == "Zymo_HostZero-MS7+ITS4"))>0) %>% tax_table() %>% row.names()
D <- ps_venn %>% subset_samples(Method == "Qiagen_PowerSoil-MS7+ITS4") %>% 
  subset_taxa(taxa_sums(ps_venn %>% subset_samples(Method == "Qiagen_PowerSoil-MS7+ITS4"))>0) %>% tax_table() %>% row.names()


full <- unique(c(A,B,C,D))
n12 <- sum(full %in% unique(A) & full %in% unique(B))
n13 <- sum(full %in% unique(A) & full %in% unique(C))
n14 <- sum(full %in% unique(A) & full %in% unique(D))
n23 <- sum(full %in% unique(B) & full %in% unique(C))
n24 <- sum(full %in% unique(B) & full %in% unique(D))
n34 <- sum(full %in% unique(C) & full %in% unique(D))

n123 <- sum(full %in% unique(A) & full %in% unique(B) & full %in% unique(C))
n124 <- sum(full %in% unique(A) & full %in% unique(B) & full %in% unique(D))
n134 <- sum(full %in% unique(A) & full %in% unique(C) & full %in% unique(D))
n234 <- sum(full %in% unique(B) & full %in% unique(C) & full %in% unique(D))

n1234 <- sum(full %in% unique(A) & full %in% unique(B) & full %in% unique(C) & full %in% unique(D))

venn.plot1 <- VennDiagram::draw.quad.venn(length((A)),length((B)),length((C)), length((D)),
                                          n12,n13,n14,n23,n24,n34,
                                          n123,n124,n134,n234,n1234,
                                          category = c('HostZero-fITS7+ITS4','HostZero-ITS86+ITS4',
                                                       'HostZero-MS7+ITS4', 'PowerSoil-MS7+ITS4'))


dev.off()

# ASVs that are shared across the four library construction methods
shared = full[full %in% unique(A) & full %in% unique(B) & full %in% unique(C) & full %in% unique(D)]
ps_venn@otu_table[,shared] %>% as.data.frame() %>% rowSums()
ps_venn@tax_table[shared] %>% as.data.frame() %>% dplyr::select(Family)

# Adding % of shared ASVs to metadata
meta$Shared <- ps_venn@otu_table[,shared] %>% as.data.frame() %>% rowSums()

# Assessing if % of shared ASVs are statistically different across library construction methods
mod_shared <- lm((Shared) ~ Method*Species, data = meta)
plot(mod_shared)
shapiro.test(resid(mod_shared))
summary(mod_shared)
anova(mod_shared)
emmeans(mod_shared, specs = pairwise ~ Method)

meta %>% group_by(Method) %>% summarise(mean_perc = mean(Shared), se_perc = se(Shared))
meta %>% group_by(Species) %>% summarise(mean_perc = mean(Shared), se_perc = se(Shared))

# composition of shared 
my_subset <- subset(tax_table(ps_venn), rownames(tax_table(ps_venn)) %in% shared)
ps_shared <- merge_phyloseq(my_subset, otu_table(ps_venn), sample_data(ps_venn))

ps_shared_comp <- ps_shared %>% speedyseq::merge_samples2('Method', fun_otu = mean, funs = list(mean)) %>%
  microbiome::transform('compositional')

shared_compositional_plot <- ps_shared_comp %>%
  psmelt() %>%
  ggplot(aes(y = Abundance, x = Method, fill = Family)) + #change fill = Taxonomic rank
  geom_bar(stat = 'identity') +
  theme_minimal()