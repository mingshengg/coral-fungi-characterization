# Load packages ####
library(ggplot2)
library(phyloseq)
library(tidyverse); packageVersion("tidyverse")
library(Biostrings)
library(vegan)

# Read ps object
ps_sp <- readRDS("./Output/ps_itsxpress_assigned.RDS")
meta <- ps_sp@sam_data %>% as.data.frame

# Visualize successful taxonomic assignments at each Rank ####
# Restructure data-frame
full_ps <- ps_sp
kin <- tax_table(ps_sp)[,1] %>% as.data.frame() != 'Unassigned'
phy <- tax_table(ps_sp)[,2] %>% as.data.frame() != 'NA'
cla <- tax_table(ps_sp)[,3] %>% as.data.frame() != 'NA'
ord <- tax_table(ps_sp)[,4] %>% as.data.frame() != 'NA'
fam <- tax_table(ps_sp)[,5] %>% as.data.frame() != 'NA'
gen <- tax_table(ps_sp)[,6] %>% as.data.frame() != 'NA'
spp <- tax_table(ps_sp)[,7] %>% as.data.frame() != 'NA'
assignments <- data.frame(Kingdom = kin, Phylum=phy, Class=cla,Order=ord,Family=fam,Genus=gen,Species=spp)

# plot assignments at each rank
assignments %>% pivot_longer(1:7) %>% mutate(name=factor(name,levels = c("Kingdom", "Phylum","Class","Order","Family","Genus","Species"))) %>%
  ggplot(aes(x=name,fill=value)) + geom_bar() + scale_fill_manual(values=c("Gray","Black")) +
  labs(x="Taxonomic level",y="Count",fill="Unambiguous\nassignment")

ggsave("./Plots/Taxonomic_Assignment_Efficiency_at_Each_Taxonomic_Rank.png",dpi=300)
rm(phy,cla,ord,fam,gen,spp,assignments,ps_sp)

# Assessing success in amplifying mock
mock <- ps_sp %>% subset_taxa(Kingdom == 'k__Fungi') %>% subset_samples(Species == 'Mock')
mock <- prune_taxa(taxa_sums(mock) >0, mock)

mock@otu_table
mock@tax_table %>% data.frame() %>% dplyr::select(Order)

mock@refseq[c('ASV9','ASV54')] %>% Biostrings::writeXStringSet('./ASVs/cryptococcus_neoformans.fa')
mock@refseq[c('ASV508','ASV1660')] %>% Biostrings::writeXStringSet('./ASVs/Malassezia_globosa.fa')

ps_mock <- mock %>% merge_samples2('Primer', fun_otu = mean, funs = list(mean))
ps_mock_comp <- ps_mock %>% microbiome::transform('compositional')

compositional_plot <- ps_mock_comp %>%
  psmelt() %>% 
  ggplot(aes(y = Abundance, x = Sample, fill = Genus)) + #change fill = Taxonomic rank
  geom_bar(stat = 'identity') +
  theme_minimal()

ggsave('./Plots/mock_composition.pdf')