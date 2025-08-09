# Load packages ####
library(dada2); packageVersion("dada2")
library(purrr); packageVersion("purrr")
library(tidyverse); packageVersion("tidyverse")
library(readxl); packageVersion("readxl")
library(ShortRead); packageVersion("ShortRead")

# Load metadata ####
full_meta <- readxl::read_xlsx("./metadata.xlsx")

# Find raw fastq files and prepare workspace ####
path <- "./ITSxpress"
fqs <- list.files(path, full.names = TRUE, recursive = FALSE, pattern = "_trimmed_r1.fastq.gz$|_trimmed_r2.fastq.gz$")

# Parse fwd and rev reads
fnFs <- fqs[grep("_trimmed_r1.fastq.gz",fqs)]

# Get Sample Names
`%notin%` = Negate(`%in%`)
sample.names <- str_remove(basename(fnFs),"_trimmed_r1.fastq.gz")
full_meta %>% filter(SampleID %notin% sample.names) %>% select(SampleID)# should be 0, to indicate all samples are accounted for
length(sample.names) == nrow(full_meta) # should also be true to indicate all accounted for
meta <- full_meta %>% filter(SampleID %in% sample.names)

# subset metadata to this run's samples
meta <- full_meta %>% filter(SampleID %in% sample.names)

# Make filtered outfile names
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))

# make new directory for filtered files
if(!dir.exists(file.path(path,"filtered"))){
  dir.create(file.path(path,"filtered"))
}

# check for duplicated sample names
sum(duplicated(sample.names))

# prefilter reads, no need for length truncation due to ITSxpress ITS extraction
out <- filterAndTrim(fnFs, filtFs,
                     maxN=0, truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE,verbose = TRUE)

saveRDS(out,"./Output/out.RDS") # save object

# reload point
out <- readRDS("./Output/out.RDS") 
filtpath <- file.path(path,"filtered")

## get sample names again
sample.names <- str_remove(basename(fnFs),"_trimmed_r1.fastq.gz")

# reassign filts for any potentially lost samples
filtFs <- list.files(filtpath, pattern = "_F_filt", full.names = TRUE)

# Learn error rates ####
errF <- learnErrors(filtFs, multithread=TRUE)
saveRDS(errF,"./Output/errF.RDS") # save object

errF <- readRDS("./Output/errF.RDS") # reload point
out <- readRDS("./Output/out.RDS")

# plot error rates for sanity
plotErrors(errF, nominalQ=TRUE)

# Infer sequence variants ####

# add names to filts
names(filtFs) <- sample.names

# Sample inference ####
dadaFs <- dada(filtFs, err = errF, multithread = TRUE)

saveRDS(dadaFs,"./Output/dadaFs.RDS")

# Inspect the merger data.frame from the first sample
head(dadaFs[[1]])

dadaFs <- readRDS('./Output/dadaFs.RDS')

# Construct sequence table ####
seqtab <- makeSequenceTable(dadaFs)
saveRDS(seqtab,"./Output/seqtab.RDS") # save object
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

# Remove chimeras ####
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

# Save progress
saveRDS(seqtab.nochim,"./Output/seqtab.nochim.RDS")
seqtab.nochim = readRDS("./Output/seqtab.nochim.RDS")

#Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), rowSums(seqtab.nochim))

# If processing a single sample, remove the sapply calls: e.g. replace
# sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "DADA", "nonchim")

rownames(track) <- sample.names
head(track)

write.csv(track, './Output/track_reads.csv')

# Assign taxonomy
# export to qiime
uniquesToFasta(seqtab.nochim, fout='./Output/rep-seqs.fna', ids=colnames(seqtab.nochim))

# save as phyloseq object
library(phyloseq)

seqtab.nochim <- readRDS('./Output/seqtab.nochim.RDS')

taxa <- readRDS('./Output/taxa.RDS')

seqtab.df <- as.data.frame(seqtab.nochim)
row.names(seqtab.df)

otu <- otu_table(seqtab.nochim,taxa_are_rows = FALSE)
met <- sample_data(full_meta)
tax <- tax_table(taxa)

sample_names(met) <- met$SampleID

names <- rownames(otu)
rownames(otu) <- sapply(strsplit(names,"_trim"), `[`, 1)

sample_names(met) <- met$SampleID
taxa_names(tax) <- colnames(otu)

ps <- phyloseq(otu,met,tax)

sequences <- Biostrings::DNAStringSet(taxa_names(ps))
names(sequences) <- taxa_names(ps)
ps <- merge_phyloseq(ps, sequences)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

saveRDS(ps, './Output/ps_itsxpress.RDS')

#### Export to BLAST for further taxonomic assignment ####
ps <- readRDS('./Output/ps_itsxpress.RDS')
ps_unassigned <- ps %>% subset_taxa(Kingdom == "Unassigned")
ps_unassigned@refseq %>% Biostrings::writeXStringSet('./Output/unassigned_ASVs.fa')

ps_fungi_only <- ps %>% subset_taxa(Kingdom == 'k__Fungi' & Class == 'NA')
ps_fungi_only@refseq %>% Biostrings::writeXStringSet('./Output/fungi_only_ASVs.fa')

# submit to 2.2 blastn.sh
