# Load packages ####
library(phyloseq)
library(emmeans)
library(dplyr)
library(ggplot2)
library(speedyseq)
library(vegan)

# Standard error formula #
se <- function(x) {
  sd(x)/sqrt(length(x))
}

# Read ps object
ps_rare <- readRDS('./Output/ps_fungi_rare.RDS') %>% 
  subset_samples(Species == 'Sediment')

# 1. compare alpha-diversity indices ####
meta = microbiome::meta(ps_rare)
meta$Shannon <- vegan::diversity(otu_table(ps_rare),index = "shannon")
meta$Richness <- vegan::specnumber(otu_table(ps_rare))
meta$Evenness <- microbiome::evenness(otu_table(ps_rare)) %>% .$pielou

# add variables
meta$Method <- paste0(meta$Primer,'_', meta$DNA_extraction)
meta$Primer_1 <- ifelse(meta$Primer == 'MS7+ITS4','F3','F1/F2')

ps_rare@sam_data$Method <- meta$Method

# Richness ANOVA model
richness_mod <- lm(data = meta, (Richness) ~ Primer*DNA_extraction + Site)
plot(richness_mod)
shapiro.test(resid(richness_mod))
summary(richness_mod)
anova(richness_mod)

emmeans(richness_mod, specs = pairwise ~ Primer|DNA_extraction)$contrasts %>% write.csv('./Output/Sediment/richness_nested_emmeans_sed.csv')
emmeans(richness_mod, specs = pairwise ~ DNA_extraction|Primer)
emmeans(richness_mod, specs = pairwise ~ Primer*DNA_extraction)$contrasts %>% write.csv('./Output/Sediment/richness_emmeans_sed.csv')

meta %>% group_by(Primer_1) %>%
  summarise(Mean = mean(Richness), SE = se(Richness))

# Shannon ANOVA model
shannon_mod <- lm(data = meta, (Shannon) ~ Primer*DNA_extraction + Site)

plot(shannon_mod)
shapiro.test(resid(shannon_mod))
summary(shannon_mod)
anova(shannon_mod)

emmeans(shannon_mod, specs = pairwise ~ Primer)
emmeans(shannon_mod, specs = pairwise ~ Primer|DNA_extraction)
emmeans(shannon_mod, specs = pairwise ~ DNA_extraction|Primer)
emmeans(shannon_mod, specs = pairwise ~ Primer*DNA_extraction)$contrasts %>% write.csv('./Output/Sediment/shannon_emmeans_sed.csv')

meta %>% group_by(Primer_1) %>%
  summarise(Mean = mean(Shannon), SE = se(Shannon))

# Evenness ANOVA model
evenness_mod <- lm(data = meta, (Evenness) ~ Primer*DNA_extraction + Site)
plot(evenness_mod)
shapiro.test(resid(evenness_mod))
summary(evenness_mod)
anova(evenness_mod)

emmeans(evenness_mod, specs=pairwise~Primer|DNA_extraction)
emmeans(evenness_mod, specs=pairwise~DNA_extraction|Primer)

# Average alpha-diversity indicies
meta %>% group_by(Primer, DNA_extraction) %>%
  summarise(richness_mean = mean(Richness), richness_se = se(Richness),
            shannon_mean = mean(Shannon), shannon_se = se(Shannon),
            evenness_mean = mean(Evenness), evenness_se = se(Evenness))

meta %>% group_by(DNA_extraction) %>%
  summarise(richness_mean = mean(Richness), richness_se = se(Richness),
            shannon_mean = mean(Shannon), shannon_se = se(Shannon),
            evenness_mean = mean(Evenness), evenness_se = se(Evenness))

meta %>% group_by(Primer_1) %>%
  summarise(richness_mean = mean(Richness), richness_se = se(Richness),
            shannon_mean = mean(Shannon), shannon_se = se(Shannon),
            evenness_mean = mean(Evenness), evenness_se = se(Evenness))

# 3. compositional differences across methods ####
# calculate bray-curtis dissimilarity
mdf_com <- ps_rare %>% otu_table() %>% vegdist()

# NMDS
otu_res <- metaMDS(mdf_com, center = T, k = 3)
otu_res

data.scores = as.data.frame(scores(otu_res, display =  'sites'))
data.scores$Primer = meta$Primer
data.scores$Site = meta$Site
data.scores$DNA_extraction = meta$DNA_extraction
meta$Method = paste0(meta$Primer,'_',meta$DNA_extraction)
data.scores$Method = meta$Method

data.scores <- data.scores[order(rownames(data.scores)),]

# Plot NMDS
gg = ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(data = data.scores, aes(colour = DNA_extraction), size = 3, alpha = 1) +
  geom_polygon(stat = "ellipse", aes(fill = DNA_extraction), alpha = 0.3) +
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
        legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 9, colour = "grey30")) + 
  labs(title = 'Stress = 0.1525575') +
  xlab('NMDS1') +
  ylab('NMDS2')

gg

# PERMANOVA to investigate difference in composition across library construction methods
perm_model <- adonis2(mdf_com ~ meta$DNA_extraction*meta$Primer + meta$Site, by = 'terms')
method_pair <- pairwiseAdonis::pairwise.adonis(mdf_com, factors = meta$Method)
write.csv(method_pair, './Output/Sediment/pairwise_adonis_sed.csv')

# betadisper to assess homogeneity
bd <- betadisper(mdf_com, meta$DNA_extraction)
anova(bd)
permutest(bd, pairwise = T)
plot(bd)

# Compositional plots
ps_merge <- ps_rare %>% speedyseq::merge_samples2('Primer', fun_otu = mean, funs = list(mean))
ps_merge_comp <- ps_merge %>% microbiome::transform('compositional')

compositional_plot <- ps_merge_comp %>%
  psmelt() %>%
  ggplot(aes(y = Abundance, x = Sample, fill = Phylum)) + #change fill = Taxonomic rank
  geom_bar(stat = 'identity') +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Differential abundance using ANCOM-BC2 ####
library(ANCOMBC)

# ANCOM-BC2 takes count data
ps <- readRDS('./Output/ps_itsxpress_assigned.RDS') %>% 
  subset_samples(Species == 'Sediment') %>%
  subset_taxa(Kingdom == 'k__Fungi')

# Reformatting ps-object
tax <- as.data.frame(tax_table(ps))
tax$Phylum <- gsub("p__Mucuromycota", "p__Mucoromycota", tax$Phylum)
tax_table(ps) <- tax_table(as.matrix(tax))

ps_da <- prune_taxa(taxa_sums(ps) > 0, ps) %>%
  tax_glom('Order') # analysis at order level

ancom_result_primer <- ancombc2(data=ps_da,
                        fix_formula = "Primer",
                        p_adj_method = "BH",
                        prv_cut = 0, #prevalence filter
                        lib_cut = 0, #filtering samples based on library sizes
                        group = "Primer", #for detection of structural zeros and performing global tests
                        struc_zero = TRUE, #to detect structural zeros
                        neg_lb = TRUE) #whether to classify taxon as a structural zero using its asymptotic lower bound, recommended for n>30

ancom_result_dna <- ancombc2(data=ps_da,
                                fix_formula = "DNA_extraction",
                                p_adj_method = "BH",
                                prv_cut = 0, #prevalence filter
                                lib_cut = 0, #filtering samples based on library sizes
                                group = "DNA_extraction", #for detection of structural zeros and performing global tests
                                struc_zero = TRUE, #to detect structural zeros
                                neg_lb = TRUE)

res_df_dna <- ancom_result_dna$res
res_df_primer <- ancom_result_primer$res

# visualizing DA results
library(tidyr)

# Visualizing DA taxa across DNA extraction methods
# filtering only significant results
res_df_dna_filt <- res_df_dna %>% filter(`diff_(Intercept)` == TRUE|
                                              `diff_DNA_extractionZymo_HostZero` == TRUE)

# creating new dataframe for plotting
res_df_dna_melt <- 
  pivot_longer(
    res_df_dna_filt,
    cols = starts_with('lfc_'),
    names_to = 'group',
    names_prefix = 'lfc_',
    values_to = 'lfc'
  ) %>%
  left_join(
    pivot_longer(
      res_df_dna_filt,
      cols = starts_with('se_'),
      names_to = 'group',
      names_prefix = 'se_',
      values_to = 'se'
    ),
    by = c('taxon','group')
  )

# Re-ordering taxa so that they are in descending order
taxa_order_dna <- res_df_dna_melt %>%
  filter(group == "DNA_extractionZymo_HostZero") %>%
  arrange(desc(lfc)) %>%
  pull(taxon)

res_df_dna_melt <- res_df_dna_melt %>%
  mutate(taxon = factor(taxon, levels = taxa_order_dna))

# Plot differential abundant plot
ggplot(data = res_df_dna_melt, aes(x = lfc, y = taxon, fill = group)) +
  geom_bar(stat = 'identity',position= 'dodge') +
  geom_errorbar(aes(y = taxon, xmin = lfc-se, xmax = lfc+se), position = 'dodge', size = 0.01, alpha = 0.5) +
  theme_bw()

res_df_dna$taxa_phylum = ps_da@tax_table[res_df_dna$taxon] %>% data.frame() %>% .$Phylum
res_df_dna$taxa_order = ps_da@tax_table[res_df_dna$taxon] %>% data.frame() %>% .$Order

write.csv(res_df_dna, './Output/Sediment/DA_extraction.csv')


# Visualizing DA taxa across primers
# filtering only significant results
res_df_dna_filt <- res_df_dna %>% filter(`diff_(Intercept)` == TRUE|
                                           `diff_DNA_extractionZymo_HostZero` == TRUE)

# creating new dataframe for plotting
res_df_primer_melt <- 
  pivot_longer(
    res_df_primer_filt,
    cols = starts_with('lfc_'),
    names_to = 'group',
    names_prefix = 'lfc_',
    values_to = 'lfc'
  ) %>%
  left_join(
    pivot_longer(
      res_df_primer_filt,
      cols = starts_with('se_'),
      names_to = 'group',
      names_prefix = 'se_',
      values_to = 'se'
    ),
    by = c('taxon','group')
  )

# Re-ordering taxa so that they are in descending order
taxa_order_primer <- res_df_primer_melt %>%
  filter(group == "(Intercept)") %>%
  arrange(lfc) %>%
  pull(taxon)

res_df_primer_melt <- res_df_primer_melt %>%
  mutate(taxon = factor(taxon, levels = taxa_order_primer))

# Plot differential abundant plot
ggplot(data = res_df_primer_melt, aes(x = lfc, y = taxon, fill = group)) +
  geom_bar(stat = 'identity',position= 'dodge') +
  geom_errorbar(aes(y = taxon, xmin = lfc-se, xmax = lfc+se), position = 'dodge', size = 0.01, alpha = 0.5) +
  theme_bw()

res_df_primer$taxa_phylum = ps_da@tax_table[res_df_primer$taxon] %>% data.frame() %>% .$Phylum
res_df_primer$taxa_order = ps_da@tax_table[res_df_primer$taxon] %>% data.frame() %>% .$Order

write.csv(res_df_primer, './Output/Sediment/DA_primer.csv')

# venn diagram ####
ps_venn <- ps %>% 
  microbiome::transform('compositional')

ps_venn@sam_data$Method <- paste0(ps_venn@sam_data$Primer, '_', ps_venn@sam_data$DNA_extraction)

A <- ps_venn %>% subset_samples(Method == "fITS7+ITS4_Qiagen_PowerSoil") %>% 
  prune_taxa(taxa_sums(.)>0, .) %>% tax_table() %>% row.names()
B <- ps_venn %>% subset_samples(Method == "MS7+ITS4_Qiagen_PowerSoil") %>% 
  prune_taxa(taxa_sums(.)>0, .) %>% tax_table() %>% row.names()
C <- ps_venn %>% subset_samples(Method == "fITS7+ITS4_Zymo_HostZero") %>% 
  prune_taxa(taxa_sums(.)>0, .) %>% tax_table() %>% row.names()
D <- ps_venn %>% subset_samples(Method == "MS7+ITS4_Zymo_HostZero") %>% 
  prune_taxa(taxa_sums(.)>0, .) %>% tax_table() %>% row.names()

full <- unique(c(A,B,C,D))
n12 <- length(intersect(A,B))
n13 <- length(intersect(A,C))
n23 <- length(intersect(B,C))

n123 <- length(Reduce(intersect, list(A, B, C)))

VennDiagram::venn.diagram(list(A,B,C,D),
                          category.names = c('fITS7+PowerSoil','ITS3-Coral+PowerSoil',
                          'fITS7+Zymo','ITS3-Coral+Zymo'),
                          filename='./Plots/sed_venn_1.tiff',
                          output = TRUE)


# ASVs that are shared across the four library construction methods
shared = A[full %in% unique(A) & full %in% unique(B) & full %in% unique(C) & full %in% unique(D)]
ps_venn@otu_table[,shared] %>% as.data.frame() %>% rowSums()
ps_venn@tax_table[shared] %>% as.data.frame() %>% select(Genus)

# Adding % of shared ASvs to metadata
ps_genus <- ps_rare
meta_qiagen <- microbiome::meta(ps_genus)
meta_qiagen$Shared <- ps_genus@otu_table[,shared] %>% as.data.frame() %>% rowSums()

# Assessing if % of shared ASVs are statistically different across library construction methods
mod_shared <- lm((Shared) ~ Primer + Site, data = meta_qiagen)
plot(mod_shared)
shapiro.test(resid(mod_shared))
summary(mod_shared)
anova(mod_shared)
emmeans(mod_shared, specs = pairwise~Primer)

mean(meta_qiagen$Shared)
se(meta_qiagen$Shared)

meta_qiagen %>% group_by(Primer) %>% 
  summarise(mean = mean(Shared), se = se(Shared))

# venn diagram across extraction method
ps_venn <- prune_taxa(taxa_sums(ps_venn) > 0, ps_venn)
ps_venn_F1 <- ps_venn %>% subset_samples(Primer_set =='F1')

A <- ps_venn %>% subset_samples(DNA_extraction == "Qiagen_PowerSoil") %>% 
  prune_taxa(taxa_sums(.)>0, .) %>% tax_table() %>% row.names()
B <- ps_venn %>% subset_samples(DNA_extraction == "Zymo_HostZero") %>% 
  prune_taxa(taxa_sums(.)>0, .) %>% tax_table() %>% row.names()

VennDiagram::venn.diagram(list(A,B),
                          category.names = c('Qiagen_PowerSoil','Zymo_HostZero'),
                          filename = './Plots/sed_venn_extraction.png',
                          output=TRUE)

ps_venn@tax_table[setdiff(B,A)] %>% data.frame() %>% dplyr::count(Phylum)
ps_venn@tax_table[setdiff(A,B)] %>% data.frame() %>% dplyr::count(Phylum)

ps_venn@tax_table[intersect(A,B)] %>% data.frame() %>% dplyr::count(Phylum)
ps_venn@tax_table[B] %>% data.frame() %>% dplyr::count(Phylum)
ps_venn@tax_table[A] %>% data.frame() %>% dplyr::count(Phylum)

# venn diagram across primers
A <- ps_venn %>% subset_samples(Primer == "fITS7+ITS4") %>% 
  prune_taxa(taxa_sums(.)>0, .) %>% tax_table() %>% row.names()
B <- ps_venn %>% subset_samples(Primer == "MS7+ITS4") %>% 
  prune_taxa(taxa_sums(.)>0, .) %>% tax_table() %>% row.names()

VennDiagram::venn.diagram(list(A,B),
                          category.names = c("fITS7+ITS4","MS7+ITS4"),
                          filename = './Plots/sed_venn_primer.png',
                          output=TRUE)

ps_venn@tax_table[setdiff(B,A)] %>% data.frame() %>% dplyr::count(Phylum)
ps_venn@tax_table[intersect(A,B)] %>% data.frame() %>% dplyr::count(Phylum)
ps_venn@tax_table[A] %>% data.frame() %>% dplyr::count(Phylum)

# Scatter plots of taxa across different methods ####
# Investigate congruence of taxa relative abundance across different methods: R^2 = 1 : Perfect congruence
# Not used in final analysis
ps@sam_data$Primer_2 <- ifelse(ps@sam_data$Primer_set == 'F3','New', 'Old')
ps@sam_data$Method_overall <- paste0(ps@sam_data$DNA_extraction, '_', ps@sam_data$Primer_2)
ps_scatter <- ps %>% speedyseq::merge_samples2('Method_overall', fun_otu = mean, funs = list(mean))
ps_scatter_d <- ps_scatter %>% tax_glom('Family')

otu_old_powersoil <- ps_scatter_d %>% subset_samples(Method_overall == 'Qiagen_PowerSoil_Old') %>% microbiome::transform('compositional') %>% otu_table() %>% as.data.frame()
otu_new_powersoil <- ps_scatter_d %>% subset_samples(Method_overall == 'Qiagen_PowerSoil_New') %>% microbiome::transform('compositional') %>% otu_table() %>% as.data.frame()
otu_old_hostzero <- ps_scatter_d %>% subset_samples(Method_overall == 'Zymo_HostZero_Old') %>% microbiome::transform('compositional') %>% otu_table() %>% as.data.frame()
otu_new_hostzero <- ps_scatter_d %>% subset_samples(Method_overall == 'Zymo_HostZero_New') %>% microbiome::transform('compositional') %>% otu_table() %>% as.data.frame()

compare_primer <- rbind(otu_old_powersoil, otu_new_powersoil) %>% t()
ggplot(data = compare_primer, aes(x = Qiagen_PowerSoil_Old, y = Qiagen_PowerSoil_New)) +
  geom_point()

summary(lm(Qiagen_PowerSoil_Old ~ Qiagen_PowerSoil_New, data = data.frame(compare_primer)))

compare_extraction <- rbind(otu_old_hostzero, otu_old_powersoil) %>% t()
ggplot(data = compare_extraction, aes(x = Zymo_HostZero_Old , y = Qiagen_PowerSoil_Old)) +
  geom_point()

summary(lm(Zymo_HostZero_Old ~ Qiagen_PowerSoil_Old, data = data.frame(compare_extraction)))

compare_both <- rbind(otu_new_hostzero, otu_old_powersoil) %>% t()
ggplot(data = compare_both, aes(x = Zymo_HostZero_New , y = Qiagen_PowerSoil_Old)) +
  geom_point()

summary(lm(Zymo_HostZero_New ~ Qiagen_PowerSoil_Old, data = data.frame(compare_both)))
