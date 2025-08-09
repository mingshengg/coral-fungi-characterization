### taxonomy update from unassigned ASVs
taxa <- readxl::read_xlsx('./Output/BLAST results/processed/fungi_only_ASVs.xlsx') %>% as.data.frame()
taxa <- readxl::read_xlsx('./Output/BLAST results/processed/unassigned_ASVs.xlsx') %>% as.data.frame() #repeat twice for unassigned and those only with kingdom fungi assignment

taxa_dat <- taxa[,2]
n_tab <- length(taxa_dat)
taxa_table <- matrix(, nrow = n_tab, ncol = 8)

for (i in 1:n_tab){
  temp = strsplit(taxa_dat[i], ';')
  taxa_table[i,] = c(temp[[1]], rep('NA', 8 - length(temp[[1]])))
}

rownames(taxa_table) <- taxa[,1]

ps <- readRDS('./Output/ps_itsxpress.RDS')
tax <- ps@tax_table %>% as.data.frame()

tax %>% filter(rownames(tax) %in% rownames(taxa_table))

taxa_table -> tax[match(rownames(taxa_table), rownames(tax)), ]

## checking to ensure that taxonomy has been updated
tax %>% filter(rownames(tax) %in% rownames(taxa_table))

ps_new <- ps
ps_new@tax_table <- tax_table(as.matrix(tax))
temp <- ps_new@tax_table %>% as.data.frame()

saveRDS(ps_new, './Output/ps_itsxpress_assigned.RDS')
