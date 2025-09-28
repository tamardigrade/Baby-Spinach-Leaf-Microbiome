#_16S CIDA Code
#Loading the environment
#load('/local1/workdir1/tw488/RData_CIDA/Basic_16S_results.RData')
load('/local1/workdir1/tw488/RData_CIDA/x16S_analysiscode_forgithub.RData')
#does not include mock ps object
#VIEW THE SAVED FILE IN A THE DIRECTORY
dir()
#SAVE ENVIRONMENT
save.image(file='/local1/workdir1/tw488/RData_CIDA/x16S_analysiscode_forgithub.RData')

#PACKAGES NEEDED____________________________________________________________________
install.packages("BiocManager", repos = "https://cloud.r-project.org")
BiocManager::install("phyloseq")
library(tibble)
library(ggplot2)
library(tidyr)
library(tibble)
library(stringr)
library(maditr)
install.packages("devtools")
library("devtools")
install_github("biobakery/Maaslin2")
library(Maaslin2)
library("theseus")
library(dplyr)
library(phyloseq)
library(vegan)
library(corrplot)
library(devtools)
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
install.packages("ggsignif")
library(ggsignif)
library(pbkrtest)
library(lme4)
library(emmeans)
library(multcompView)
library(rcompanion)
#if (!requireNamespace("BiocManager", quietly = TRUE))y
#install.packages("BiocManager")
#BiocManager::install(version = '3.9')
#BiocManager::install("dada2", version = "3.9")
#remove.packages("rlang")
#install.packages("rlang")
library(dada2); packageVersion("dada2")
library(Biostrings); packageVersion("Biostrings")
#BiocManager::install("decontam")
library(decontam); packageVersion("decontam")
#Just for reference:for testing normality
library(ggpubr)
library(ShortRead)
library("tools")

#Most of sample processing workflow is from https://benjjneb.github.io/dada2/ITS_workflow.html

#map the path to the sequencing files on my drive
path <- "/local1/workdir1/tw488/CIDA_Sequencing_Runs/16S_RUNS1_7_allinonev2/fastq.gz_files_forprocessing/"
list.files(path)

#generate matched lists of the forward and reverse read files

fnFs <- sort(list.files(path, pattern = "_R1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2.fastq.gz", full.names = TRUE))

#identify the primers used to make these amplicons. Inorder to get it to actually work,
# you need to (i) exclude the Illumina adaptor part of the primer (generally the 5' end of each primer)
#and (ii) put the reverse and forward primers in below at 5' to 3' 

FWD <- "CADACTCCTACGGGAGGC"  
REV <- "ATCCTGTTTGMTMCCCVCRC"  

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
#check the orientations, make sure theyre right
FWD.orients
REV.orients
#both 5'-3'
#FWD <- "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGCADACTCCTACGGGAGGC"  
#REV <- "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGATCCTGTTTGMTMCCCVCRC" 
#_________________
#remove poly Ns from reads, create output filtN
fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filtered files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))

#************this step takes a long time
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)

#counting primers in F and R reads
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]),multithread = TRUE)

#the above step took awhile too
#what are the results? all F reads in F position, all R reads in R position?

#remove primers with cutadapt- 
cutadapt <-"/workdir1/tw488/CIDA_Sequencing_Runs/cutadapt-venv/bin/cutadapt"
system2(cutadapt, args = "--version")

#create output filenames for the cutadapt-ed files,
path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 
# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads 
                             "-m", 1, "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}

#double check for primers, output should all be 0
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]),multithread = TRUE)

# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "_R1.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_R2.fastq.gz", full.names = TRUE))

q<-basename(file_path_sans_ext(cutFs))
#this only got rid of .gz, need to get rid of .fastq as well
sample.names <-(file_path_sans_ext(q))
#.fastq.gz should be removed
head(sample.names)
list(sample.names)
type(sample.names)

#inspect quality reads
plotQualityProfile(cutFs[1:2])
plotQualityProfile(cutRs[1:2])

#_____start following this workflow: https://benjjneb.github.io/dada2/tutorial.html
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))

out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, truncLen=c(230,230), maxN = 0, maxEE = c(2, 2), 
                     truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE)  # on windows, set multithread = FALSE

head(out)
list(out)
plotQualityProfile(filtFs[1:2])
plotQualityProfile(filtRs[1:2])
class(out)

list(filtFs)
class(filtFs)

list(filtFs)
list(filtRs)
list(filtRs)
list(sample.names)

#Learn the Error Rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)

#run DADA
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

#Merge paired reads
#Most of your reads should successfully merge. If that is not the case upstream parameters may need to be revisited: Did you trim away the overlap between your reads?
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

#Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
#Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)


#Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
list(track)

#code to check # of reads to make box plots to look for outlier in read counts_________________
track<-as.data.frame(track)
track <- tibble::rownames_to_column(track, "samplename")
track <- track %>%mutate(newcol=samplename)
track <- track %>%separate(newcol, c("16s","month", "lot","day","rep","loc","r1"), "_",remove = TRUE)

#Assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim, "/local1/workdir1/tw488/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
data.frame(taxa.print)

#___________________________________
## 1. ASV Table
# Prep the ASV table! 
samples_out <- rownames(seqtab.nochim)
# Pull out sample names from the fastq file name 
#sample_names_reformatted <- sapply(strsplit(samples_out, split = "_"), `[`, 1)
l<-basename(file_path_sans_ext(samples_out))
#this only got rid of .gz, need to get rid of .fastq as well
sample_names_reformatted <-(file_path_sans_ext(l))
#.fastq.gz should be removed
head(sample_names_reformatted)
# Replace the names in our seqtable 
rownames(seqtab.nochim) <- sample_names_reformatted
### intuition check : result shoul be 'NULL'
stopifn=stopifnot(rownames(seqtab.nochim) == sample_names_reformatted)
stopifn

############## Modify the ASV names 
# Give headers more manageable names
# First pull the ASV sequences
asv_seqs <- colnames(seqtab.nochim)
# make headers for our ASV seq fasta file, which will be our asv names
asv_headers <- vector(dim(seqtab.nochim)[2], mode = "character")
# loop through vector and fill it in with ASV names 
for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep = "_")
}
# intitution check
asv_headers
##### Rename ASVs in table then write out our ASV fasta file! 
View(seqtab.nochim)
asv_tab <- t(seqtab.nochim)
#View(asv_tab)
## Rename our asvs! 
row.names(asv_tab) <- sub(">", "", asv_headers)
View(asv_tab)

## 2. Taxonomy Table 
#```{r prepare-tax-table}
View(taxa)
##### Prepare tax table 
# Add the ASV sequences from the rownames to a column 
new_tax_tab <- taxa %>%
  as.data.frame() %>%
  rownames_to_column(var = "ASVseqs") 
head(new_tax_tab)
# intution check 
h<-stopifnot(new_tax_tab$ASVseqs == colnames(seqtab.nochim))
h
list(asv_tab)
type(asv_tab)

# Now let's add the ASV names 
rownames(new_tax_tab) <- rownames(asv_tab)
View(new_tax_tab)
### Final prep of tax table. Add new column with ASV names 
asv_tax <- 
  new_tax_tab %>%
  # add rownames from count table for phyloseq handoff
  mutate(ASV = rownames(asv_tab)) %>%
  # Resort the columns with select
  dplyr::select(Kingdom, Phylum, Class, Order, Family, Genus, ASV, ASVseqs)
View(asv_tax)
# Intution check
j<-stopifnot(asv_tax$ASV == rownames(asv_tax), rownames(asv_tax) == rownames(asv_tab))
j

theme_set(theme_bw())
samples.out <- rownames(seqtab.nochim)
# Pull out sample names from the fastq file name 
#sample_names_reformatted <- sapply(strsplit(samples_out, split = "_"), `[`, 1)
l<-basename(file_path_sans_ext(samples.out))
#this only got rid of .gz, need to get rid of .fastq as well
sample_names_reformatted <-(file_path_sans_ext(l))
#.fastq.gz should be removed
head(sample_names_reformatted)
#function(fname) strsplit(basename(fname), "_")[[1]][1]
SampleName=sample.names
#SampleName=sample_names
Month <- sapply(strsplit(samples.out, "_"), `[`, 2)
Lot <- sapply(strsplit(samples.out, "_"), `[`, 3)
Day <- sapply(strsplit(samples.out, "_"), `[`, 4)
Rep <- sapply(strsplit(samples.out, "_"), `[`, 5)
Loc<- sapply(strsplit(samples.out, "_"), `[`, 6)
samdf <- data.frame(SampleName=SampleName,Month=Month, Lot=Lot, Day=Day, Rep=Rep,Loc=Loc)
rownames(samdf) <- samples.out
data.frame(samples.out)
#need to update 1122B files that have CA instead of GA-
rownames(samdf)[258] <- "16S_1122_B_D12_2_GA_R1"
rownames(samdf)[261] <- "16S_1122_B_D22_2_GA_R1"

#MAKE SAMPLE DATA FILE WITH METADATA FOR USE OF THIS
#We now construct a phyloseq object directly from the dada2 outputs.
asv_tab<-as.matrix(asv_tab)
asv_tax<-as.matrix(asv_tax)
ps <- phyloseq(otu_table(asv_tab, taxa_are_rows=TRUE), 
               sample_data(samdf), 
               tax_table(as.matrix(asv_tax)))
sample_variables(ps)
ps

#____________________________________________________________________________________
#subset mock samples
mock <- subset_samples(ps, Month== "mock")

#_____________________________________________________________________
#remove contaminant OTUs with decontam, make new ps object without them_______________
#read in merge samdf and decondf
decondf <- read.csv("/local1/workdir1/tw488/CIDA_Sequencing_Runs/16S_RUNS1_7_allinone/fastq.gz_files_forprocessing/decondf_16S_finalrun.csv",row.names = 1)
samdf2 <- read.csv("/local1/workdir1/tw488/CIDA_Sequencing_Runs/16S_RUNS1_7_allinone/fastq.gz_files_forprocessing/samdf16S_withactualday.csv",row.names = 1)
decondf <- merge(decondf, samdf2, by = 'row.names', all = TRUE)
decondf<-decondf %>% column_to_rownames(var="Row.names")
decondf <- decondf %>%mutate_at("SampleName", str_replace, "_R1", "")
str(decondf)

## Handoff to phyloseq 
ps1 <- phyloseq(otu_table(asv_tab, taxa_are_rows = TRUE),
                sample_data(decondf_nomock),
                tax_table(as.matrix(asv_tax)))

df_samps1<- as.data.frame(ps1@sam_data)

ps1 <- subset_samples(ps1, Month!= "mock")


with(sample_data(ps1), table(SampleName))
with(sample_data(ps1), table(quant_reading))

#run decontam wth combined method
#cannot run cont with 0 ng/uL value for quant_reading, so take half of smallest value and replace 0 with that
cont<-isContaminant(ps1, conc = "quant_reading", neg = "is.neg", method = "combined",
                    batch = NULL, threshold = 0.05,normalize = TRUE,detailed = TRUE)
cont

#check out score dist on ggplot histogram
p<-ggplot(cont, aes(x=p.freq)) + geom_histogram()
p

o<-ggplot(cont, aes(x=p.prev)) + geom_histogram()
o

i<-ggplot(cont, aes(x=p)) + geom_histogram()
i

#filter out the contaminants
cont1<-cont%>%
  filter(contaminant=="TRUE")

table(cont$contaminant)
head(which(cont$contaminant))
#match the cont1 to tax table and make sure the contaminant ASVs look right
contamTaxtbl <- merge(cont1, asv_tax, by = 'row.names', all = TRUE)
contamTaxtbl <-contamTaxtbl %>%filter(contaminant=="TRUE")

contamTaxtbl_1 <-contamTaxtbl[ -c(8:12) ]

#Make new PS without contamns
ps.noncontam <- prune_taxa(!cont$contaminant, ps1)
ps.noncontam

#check if the contam ASVs have actually been filtered out
df_psnoncontam_taxtbl<- as.data.frame(ps.noncontam@tax_table)

#code to remove ASVs==1 and chloroplast
ps.noncontam <- ps.noncontam %>% subset_taxa( Family!= "mitochondria" | is.na(Family) & Class!="Chloroplast" | is.na(Class) )
# remove taxa with 1 read count, from https://deneflab.github.io/Diversity_Productivity/analysis/OTU_Removal_Analysis.html
ps.noncontam<-prune_taxa(taxa_sums(ps.noncontam) > 1, ps.noncontam) 

ps.noncontam <- subset_samples(ps.noncontam, Month!= "mock")

# check to see if above worked by summing rows and columns, as well as 
smdta<-as.data.frame(ps.noncontam@sam_data)
otutbl<-as.data.frame(ps.noncontam@otu_table)
txtbl<-as.data.frame(ps.noncontam@tax_table)

#check the sum of rows/ASVs
otu_tbl_rowsums<-as.data.frame(rowSums(otutbl))
#no row sums under 2
#check the sum of rows/samples
otu_tbl_colsums<-as.data.frame(colSums(otutbl))
#what is the smallest column sum, i.e. which sample has the smallest # of reads?

#____16S MOCK COMMUNITY ANALYSIS BETWEEN SEQUENCING RUNS____________________

#count_seqs function for use in calculating relative abundance
count_seqs <- function(pool){
  count_dat <- sample_sums(pool)
  names <- names(count_dat)
  values <- unname(count_dat)
  counts_1 <- cbind(names, values) %>% as.data.frame()
  colnames(counts_1) = c("sample", "counts") 
  
  counts_1 <- counts_1 %>%
    arrange(as.numeric(counts)) %>%
    mutate(counts = as.numeric(counts))
  counts_1$n <- seq.int(nrow(counts_1))
  
  plot <- ggplot(counts_1, aes(n, counts)) +
    geom_line()+
    theme_classic()
  
  return(list("counts"=counts_1, "plot"=plot))
  
}

#### Calculate relative abundance ####

mock.counts <- count_seqs(mock@otu_table)

abund.MOCK <- mock@otu_table %>% as.matrix() %>% as.data.frame() %>%
  rownames_to_column("ASV_") %>%
  pivot_longer((-ASV_), names_to = "sample") %>%
  merge(mock.counts$counts, by="sample") %>%
  dplyr::select(-n) %>%
  mutate(newcol=sample)%>%
  separate(newcol, c("mock","x1","x2", "x3","x4","x5"), "_",remove = TRUE)%>%
  group_by(sample) %>%
  mutate(rel_abund = value/counts) %>% #relative abundance by treatment/day/rep
  ungroup() 

#check that relative abundance is calculated correctly- everything should add up to one
abund.MOCK %>% group_by(sample) %>% summarise(sum = sum(rel_abund))

#get taxonomy
mocktax <- mock@tax_table %>%
  as.data.frame() %>%
  rownames_to_column("ASV_") 
head(mocktax)

#merge taxonomy with rarefied relative abundance
count_tax.MOCK <- merge(abund.MOCK, mocktax, by="ASV_")
head(count_tax.MOCK)

# Get summary abundance of each taxa
taxa.meta.sumMOCK <- count_tax.MOCK%>%
  dplyr::select(-c(value, counts)) %>%
  pivot_longer(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "ASV_"),
               names_to = "level", 
               values_to = "taxon")

#format phylum names
genus_rel_abundMOCK <- taxa.meta.sumMOCK %>%
  filter(level == "Genus") %>%
  group_by(taxon, sample,mock,x1, x2,x3,x4,x5) %>% 
  summarise(rel_abund = 100*sum(rel_abund), .groups="drop") %>%
  mutate(taxon = str_replace(taxon,
                             "^unclassified (.*)", "Unclassified *\\1*"),
         taxon = str_replace(taxon, "^(\\S*)$", "\\1")) %>%
  ungroup()
genus_rel_abundMOCK$taxon <-genus_rel_abundMOCK$taxon%>%replace_na('Unclassified')
genus_rel_abundMOCK <- genus_rel_abundMOCK %>% filter(rel_abund!="0") 

#make a chart to see the rel ab of each genera
genus_rel_abundMOCK<-dcast(genus_rel_abundMOCK,taxon~sample,value.var="rel_abund")

# Create variables for the Defined Composition of the 16S mock community
Genus <- c("Listeria","Pseudomonas","Bacillus","Escherichia-Shigella","Salmonella","Lactobacillus","Enterococcus","Staphylococcus")
Original_Composition <- c(89.1, 8.9, 0.89, 0.089, 0.089, 0.0089,0.00089,0.000089)

# Join the variables to create a data frame
original <- data.frame(Genus,Original_Composition)
colnames(genus_rel_abundMOCK)[1] <-"Genus"
totalmock <- merge(original,genus_rel_abundMOCK,by=c("Genus"), all=T)
totalmock <-totalmock%>%arrange(desc(Original_Composition))%>%replace(is.na(.), 0)
colnames(totalmock)[3] <-"Sequencing Run 1 Relative Abundance"
colnames(totalmock)[4] <-"Sequencing Run 2 Relative Abundance"
colnames(totalmock)[5] <-"Sequencing Run 3 Relative Abundance"
colnames(totalmock)[6] <-"Sequencing Run 4 Relative Abundance"
colnames(totalmock)[7] <-"Sequencing Run 5 Relative Abundance"
colnames(totalmock)[8] <-"Sequencing Run 6 Relative Abundance"
colnames(totalmock)[9] <-"Sequencing Run 7 Relative Abundance"
colnames(totalmock)[2] <-"Defined Composition"

write.csv(totalmock, "/local1/workdir1/tw488/CIDA/For_Publication/mock16Ssampletable_forsupp.csv", row.names=FALSE)

#___MAKE Phyloseq OBJECTS FOR EACH ANALYSIS__________________________________________

#Import phyloseq object ps.noncontam, which is from is from 16S_Runs1_7.R
ps.noncontam16SOTU<-ps.noncontam@otu_table
ps.noncontam16SOTU_colsums<-as.data.frame(colSums(ps.noncontam16SOTU))
ps.noncontam16SOTU_rowsums<-as.data.frame(rowSums(ps.noncontam16SOTU))

#REMOVE NC SAMPLES, samples being excluded AND rarefy ps.noncontam
ps.noncontam_noNC <- subset_samples(ps.noncontam, Rep!= "NC")

#remove 0722_B
ps.noncontam_noNC_exclsamples <- subset_samples(ps.noncontam_noNC, SampleName!= "16S_0722_B_H_3_CA")

#remove 0922_B shelf life samples (Harvest is ok)
ps.noncontam_noNC_exclsamples <- subset_samples(ps.noncontam_noNC_exclsamples, SampleName!= "16S_0922_B_DI_2_CA")
ps.noncontam_noNC_exclsamples <- subset_samples(ps.noncontam_noNC_exclsamples, SampleName!= "16S_0922_B_DI_3_CA")
ps.noncontam_noNC_exclsamples <- subset_samples(ps.noncontam_noNC_exclsamples, SampleName!= "16S_0922_B_D7_3_CA")
ps.noncontam_noNC_exclsamples <- subset_samples(ps.noncontam_noNC_exclsamples, SampleName!= "16S_0922_B_D7_2_CA")
ps.noncontam_noNC_exclsamples <- subset_samples(ps.noncontam_noNC_exclsamples, SampleName!= "16S_0922_B_D22_2_CA")
ps.noncontam_noNC_exclsamples <- subset_samples(ps.noncontam_noNC_exclsamples, SampleName!= "16S_0922_B_D22_1_CA")
ps.noncontam_noNC_exclsamples <- subset_samples(ps.noncontam_noNC_exclsamples, SampleName!= "16S_0922_B_D17_3_CA")
ps.noncontam_noNC_exclsamples <- subset_samples(ps.noncontam_noNC_exclsamples, SampleName!= "16S_0922_B_D17_1_CA")
ps.noncontam_noNC_exclsamples <- subset_samples(ps.noncontam_noNC_exclsamples, SampleName!= "16S_0922_B_D12_2_CA")
ps.noncontam_noNC_exclsamples <- subset_samples(ps.noncontam_noNC_exclsamples, SampleName!= "16S_0922_B_D12_1_CA")

#remove the third 0422_B sample
ps.noncontam_noNC_exclsamples <- subset_samples(ps.noncontam_noNC_exclsamples, SampleName!= "16S_0422_B_H_1_CA")

#remove extra 16S_0522_A_H_3_CA_R1_R1
ps.noncontam_noNC_exclsamples <- subset_samples(ps.noncontam_noNC_exclsamples, SampleName!= "16S_0522_A_H_3_CA_R1")

#rarefy
ps.noncontam_noNC.rarefied = rarefy_even_depth(ps.noncontam_noNC_exclsamples, rngseed=1, sample.size=min(sample_sums(ps.noncontam_noNC_exclsamples)), replace=F)
ps.noncontam_noNC.rarefiedotu<-ps.noncontam_noNC.rarefied@otu_table
ps.noncontam_noNC.rarefiedotu_colsums<-as.data.frame(colSums(ps.noncontam_noNC.rarefiedotu))
ps.noncontam_noNC.rarefiedotu_rowsums<-as.data.frame(rowSums(ps.noncontam_noNC.rarefiedotu))

#HARVEST AND DAY INITIAL ONLY
HDI <- subset_samples(ps.noncontam_noNC.rarefied,Day!= "D22")
HDI <- subset_samples(HDI,Day!= "D17")
HDI <- subset_samples(HDI,Day!= "D12")
HDI <- subset_samples(HDI,Day!= "D7")

#HARVEST ONLY
H <- subset_samples(HDI,Day!= "DI")

#SHELF LIFE SAMPLES
DI_D22_D28 <- subset_samples(ps.noncontam_noNC.rarefied,Day!= "H")

#SEVEN DAY INTERVAL SHELF LIFE SAMPLES
sevenday <- subset_samples(DI_D22_D28,Month!= "0122")
sevenday <- subset_samples(sevenday,Month!= "0222")
sevenday <- subset_samples(sevenday,Month!= "1221")
sevenday <- subset_samples(sevenday,Month!= "0322")
sevenday <- subset_samples(sevenday,Month!= "0722")
sevenday <- subset_samples(sevenday,Month!= "0822")
sevenday <- subset_samples(sevenday,Month!= "0922")
sevenday <- subset_samples(sevenday,Month!= "1022")
sevenday <- subset_samples(sevenday,Month!= "1122")
sevenday <- subset_samples(sevenday,Month!= "1222")
sevenday <- subset_samples(sevenday,Actual_Day!= "D22")
sevenday <- subset_samples(sevenday,Actual_Day!= "D12")
sevenday <- subset_samples(sevenday,Actual_Day!= "D17")
sevenday <- subset_samples(sevenday,SampleName!= "16S_0622_B_D7_2_CA")
sevenday <- subset_samples(sevenday,SampleName!= "16S_0622_B_D7_3_CA")
sevenday <- subset_samples(sevenday,SampleName!= "16S_0622_B_DI_1_CA")
sevenday <- subset_samples(sevenday,SampleName!= "16S_0622_B_DI_3_CA")
sevenday <- subset_samples(sevenday,SampleName!= "16S_0622_B_H_2_CA")
sevenday <- subset_samples(sevenday,SampleName!= "16S_0622_B_H_3_CA")

#FIVE DAY INTERVAL SHELF LIFE SAMPLES
fiveday <- subset_samples(DI_D22_D28,Month!= "0522")
fiveday <- subset_samples(fiveday,Month!= "0422")
fiveday <- subset_samples(fiveday,SampleName!= "16S_0622_A_D12_1_CA")
fiveday <- subset_samples(fiveday,SampleName!= "16S_0622_A_D12_2_CA")
fiveday <- subset_samples(fiveday,SampleName!= "16S_0622_A_D17_1_CA")
fiveday <- subset_samples(fiveday,SampleName!= "16S_0622_A_D17_2_CA")
fiveday <- subset_samples(fiveday,SampleName!= "16S_0622_A_D22_2_CA")
fiveday <- subset_samples(fiveday,SampleName!= "16S_0622_A_D22_3_CA")
fiveday <- subset_samples(fiveday,SampleName!= "16S_0622_A_D7_3_CA")
fiveday <- subset_samples(fiveday,SampleName!= "16S_0622_A_D7_1_CA")
fiveday <- subset_samples(fiveday,SampleName!= "16S_0622_A_DI_1_CA")
fiveday <- subset_samples(fiveday,SampleName!= "16S_0622_A_DI_3_CA")

#MOVING ON TO HARVEST SAMPLE ANALYSIS

#___CALCULATE RELATIVE ABUNDANCE FOR HARVEST ONLY___________

count_seqs <- function(pool){
  count_dat <- sample_sums(pool)
  names <- names(count_dat)
  values <- unname(count_dat)
  counts_1 <- cbind(names, values) %>% as.data.frame()
  colnames(counts_1) = c("sample", "counts") 
  
  counts_1 <- counts_1 %>%
    arrange(as.numeric(counts)) %>%
    mutate(counts = as.numeric(counts))
  counts_1$n <- seq.int(nrow(counts_1))
  
  plot <- ggplot(counts_1, aes(n, counts)) +
    geom_line()+
    theme_classic()
  
  return(list("counts"=counts_1, "plot"=plot))
  
}

#### Calculate relative abundance ####

H.counts <- count_seqs(H@otu_table)

abund.H <- H@otu_table %>% as.matrix() %>% as.data.frame() %>%
  rownames_to_column("ASV_") %>%
  pivot_longer(-ASV_, names_to = "sample") %>%
  merge(H.counts$counts, by="sample") %>%
  dplyr::select(-n) %>%
  mutate(newcol=sample)%>%
  separate(newcol, c("X16S","Month", "Lot","Day","Rep","Loc","Actual_Day"), "_",remove = TRUE)%>%
  group_by(sample) %>%
  mutate(rel_abund = value/counts) %>% #relative abundance by treatment/day/rep
  ungroup() 

#check that relative abundance is calculated correctly
abund.H %>% group_by(sample) %>% summarise(sum = sum(rel_abund))

#check abund.H for top ASVs w/ the highest counts, use taxH to identify the ASVs

#get taxonomy
taxH <- H@tax_table %>%
  as.data.frame() %>%
  rownames_to_column("ASV_") 
head(taxH)

#merge taxonomy with rarefied relative abundance
count_tax.H <- merge(abund.H, taxH, by="ASV_")
head(count_tax.H)

# Get summary abundance of each taxa
taxa.meta.sumH <- count_tax.H%>%
  dplyr::select(-c(value, counts)) %>%
  pivot_longer(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "ASV"),
               names_to = "level", 
               values_to = "taxon") 
#format phylum names
genus_rel_abundH <- taxa.meta.sumH %>%
  filter(level == "Genus") %>%
  group_by(taxon, sample,Month, Lot,Day,Rep,Loc,Actual_Day) %>% 
  summarise(rel_abund = 100*sum(rel_abund), .groups="drop") %>%
  mutate(taxon = str_replace(taxon,
                             "^unclassified (.*)", "Unclassified *\\1*"),
         taxon = str_replace(taxon, "^(\\S*)$", "\\1")) %>%
  ungroup()
genus_rel_abundH$taxon <-genus_rel_abundH$taxon%>%replace_na('Unclassified')

#check how many taxa there are total 
genus_rel_abundH_total<-genus_rel_abundH%>%group_by(taxon)%>%summarise(mean = mean(rel_abund))%>%filter(mean!=0)

#pool taxon <5% relative abundance into "other" category
taxon_poolH <- genus_rel_abundH %>%
  group_by(sample, taxon) %>%
  summarise(mean = mean(rel_abund), .groups="drop") %>%
  group_by(taxon) %>%
  summarise(pool = max(mean) <5,
            mean=mean(mean), 
            .groups="drop") 

#arrange data frame to group samples by relative abundance within each taxon; 
#this will aid visual comparison between the barplots
join_taxon_rel_abundH <- inner_join(genus_rel_abundH, taxon_poolH, by="taxon") %>%
  mutate(taxon = ifelse(pool, "< 5 %", taxon)) %>%
  group_by(sample, taxon,Month, Lot,Day,Rep,Loc,Actual_Day) %>%
  summarise(rel_abund = sum(rel_abund),
            mean = mean(mean), .groups="drop") %>%
  mutate(taxon = factor(taxon),
         taxon = fct_reorder(taxon, rel_abund, .desc=TRUE),
         taxon = fct_shift(taxon, n=1)) %>%
  ungroup() %>%
  unique() %>%
  group_by(taxon) %>%
  arrange(desc(rel_abund), .by_group=TRUE) 

levels(join_taxon_rel_abundH$taxon)

#MAKE A LIST OF ALL GENERA IN HARVEST SAMPLES >5% RELATIVE ABUNDANCE
totalHgenera_over5<-join_taxon_rel_abundH%>%group_by(taxon)%>%filter(rel_abund>5)%>%distinct(taxon)

#make different columns: monthlot, month-lot-day-rep (MLDR) and day-rep (DR)
join_taxon_rel_abundH$Monthlot = paste(join_taxon_rel_abundH$Month, join_taxon_rel_abundH$Lot, sep="_")
join_taxon_rel_abundH$MLDR = paste(join_taxon_rel_abundH$Month,join_taxon_rel_abundH$Lot,join_taxon_rel_abundH$Day,join_taxon_rel_abundH$Rep,sep="_")
join_taxon_rel_abundH$DR = paste(join_taxon_rel_abundH$Day,join_taxon_rel_abundH$Rep,sep="_")
join_taxon_rel_abundH$Month<-factor(join_taxon_rel_abundH$Month,
                                    levels=c("1221","0122","0222","0322","0422","0522","0622","0722","0822","0922","1022","1122","1222"))

join_taxon_rel_abundH$DR<-factor(join_taxon_rel_abundH$DR,
                                 levels=c("H_1","H_2","H_3", "DI_1","DI_2","DI_3", "D7_1","D7_2","D7_3", "D12_1","D12_2","D12_3", "D17_1","D17_2","D17_3","D22_1","D22_2","D22_3"))


color_map <- c( "< 5 %"="black","Unclassified" = "red", "Bacillus" = "blue",
                "Aerococcus"="bisque4", 
                "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium"= "chartreuse",
                "Massilia"="darkgoldenrod1", "Kosakonia"="darkorange", 
                "Clostridium sensu stricto 1"= "darkslategray2", 
                "Planomicrobium"="brown","Salinicoccus"="magenta", 
                "Brevundimonas"= "gray28", "Flavobacterium"= "lightsalmon2", 
                "Buchnera"="darkorchid4", "Romboutsia"= "gray64", "Pseudomonas"="green4",
                "Stenotrophomonas"="papayawhip", "Psychrobacter"= "yellow3", 
                "Staphylococcus"= "thistle3", "Rheinheimera"="navy", 
                "Vagococcus"="purple",  "Parvibaculum"= "palevioletred",
                "Halomonas"="skyblue", "	Prosthecobacter"= "sienna", 
                "Stenotrophomonas"="darkolivegreen2", "Prosthecomicrobium"="chocolate", 
                "Pseudolabrys"= "deeppink", "Pantoea"="aquamarine", 
                "Ramlibacter"= "grey99", "Exiguobacterium"="yellow",
                "Paracoccus"="pink","Rhodobacter"="mediumspringgreen",
                "Runella"="firebrick","Prosthecobacter"= "orange",
                "Sphingomonas"= "plum1","Sediminibacterium"= "olivedrab4", 
                "Shewanella"= "darkcyan","Planococcus"= "lightblue",
                "Chryseobacterium" ="tomato", "Methylobacterium-Methylorubrum"="grey86")

#### Relative abundance plot ####
join_taxon_rel_abundH %>%
  mutate(Month = factor(Month, c("1221","0122","0222","0322","0422","0522","0622","0722","0822","0922","1022","1122","1222"))) %>%
  ggplot()+
  geom_col(mapping = aes(x = MLDR, y = rel_abund, fill = taxon), color = "black", position = "fill", show.legend = TRUE, width=.4)+
  facet_grid(cols = vars(Loc),scales="free",space = "free")+
  ylab("Proportion of Community") +
  xlab("Samples")+
  scale_fill_manual(values = color_map) +
  theme_linedraw()+
  theme(axis.text.y = element_text(size = 8, color = "black"),
        axis.title.y = element_text(size = 10, color = "black"),
        axis.text.x = element_text(size = 8, angle = 90, vjust = 1, hjust = 0.5, color = "black"),
        legend.text = element_text(size = 10,face = "italic"),
        legend.position = "bottom",
        legend.spacing.x = unit(0.1, 'mm'),
        legend.spacing.y = unit(0.05, 'mm'),
        plot.margin=grid::unit(c(0.1,0.1,0.1,0.1), "mm"),
        strip.text = element_text(size = 18, face = "bold", angle = 0),
        legend.title = element_text(face="bold", size = 10))+
  guides(fill=guide_legend(ncol=6,byrow=TRUE))+ labs(fill='Genus') 

ggsave("Bact_RelAB_Barplot_05282025.tiff", units="in", width=20, height=15, dpi=500, compression = 'lzw')

#Calculate average relative abundance of unclassified bacteria
unclrelab<-join_taxon_rel_abundH%>%filter(taxon== "Unclassified")%>% summarise(AvgbyLoc= mean(rel_abund))
  
#____16S Harvest permanova_______________________________________________________________________________________

#PERMANOVA WITH HOW FUNCTION
Hotu<-pstoveg_otu(H)
Hsd<-pstoveg_sd(H)
Hsd$Monthlot = paste(Hsd$Month, Hsd$Lot, sep="_")
Hvegdist<-vegdist(Hotu, method="bray") 
q<-as.matrix(Hvegdist)
#make Monthlot in to a plot

#need to randomly get rid of 24-16=8 CA samples and both FL samples so the # of samples are even

H.rarefied_noFL_4permanova <- subset_samples(H,SampleName!= "16S_1221_A_H_1_FL")
H.rarefied_noFL_4permanova <- subset_samples(H.rarefied_noFL_4permanova,SampleName!= "16S_1221_A_H_3_FL")
H.rarefied_noFL_4permanova <- subset_samples(H.rarefied_noFL_4permanova,SampleName!= "16S_0522_C_H_2_CA")
H.rarefied_noFL_4permanova <- subset_samples(H.rarefied_noFL_4permanova,SampleName!= "16S_0522_C_H_3_CA")
H.rarefied_noFL_4permanova <- subset_samples(H.rarefied_noFL_4permanova,SampleName!= "16S_0622_A_H_2_CA")
H.rarefied_noFL_4permanova <- subset_samples(H.rarefied_noFL_4permanova,SampleName!= "16S_0622_A_H_3_CA")
H.rarefied_noFL_4permanova <- subset_samples(H.rarefied_noFL_4permanova,SampleName!= "16S_0922_B_H_1_CA")
H.rarefied_noFL_4permanova <- subset_samples(H.rarefied_noFL_4permanova,SampleName!= "16S_0922_B_H_2_CA")
H.rarefied_noFL_4permanova <- subset_samples(H.rarefied_noFL_4permanova,SampleName!= "16S_1022_A_H_1_CA")
H.rarefied_noFL_4permanova <- subset_samples(H.rarefied_noFL_4permanova,SampleName!= "16S_1022_A_H_2_CA")
H.rarefied_noFL_4permanovaotu<-pstoveg_otu(H.rarefied_noFL_4permanova)
H.rarefied_noFL_4permanovasd<-pstoveg_sd(H.rarefied_noFL_4permanova)
H.rarefied_noFL_4permanovavegdist<-vegdist(H.rarefied_noFL_4permanovaotu, method="bray") 

H.rarefied_noFL_4permanovasd$Monthlot = paste(H.rarefied_noFL_4permanovasd$Month, H.rarefied_noFL_4permanovasd$Lot, sep="_")
#run permanova
h <- with(H.rarefied_noFL_4permanovasd, how(within=Within(type="free"),plots=Plots(strata=Monthlot, type="free"),nperm=999))
h
table(H.rarefied_noFL_4permanovasd$Monthlot,H.rarefied_noFL_4permanovasd$Loc)
H_perm <- adonis2(H.rarefied_noFL_4permanovavegdist ~Loc, data=H.rarefied_noFL_4permanovasd, permutations = h, method = "bray", 
                  sqrt.dist = FALSE, add = FALSE, by = "terms",
                  parallel = getOption("mc.cores"), na.action = na.fail)
H_perm


#___Non-metric Multidimensional Scaling (NMDS)______________________________________________________

H_NMDS <- metaMDS(Hvegdist)
H_NMDS
stressplot(H_NMDS)

HsdNMDS<-Hsd%>%rownames_to_column("sample")
plot_dfH <- scores(H_NMDS, display = "sites") %>% 
  as.data.frame() %>% 
  rownames_to_column("sample") %>% 
  full_join(HsdNMDS, by = "sample")
palNMDS <- c("#A6CEE3", "#B2DF8A","brown", "#1F78B4","#FB9A99", "#33A02C","#E31A1C" ,"#FDBF6F" ,"#FF7F00", "#CAB2D6","darkgoldenrod","purple","grey","black","darkorchid4")

plot_dfH$Month<-gsub("0122","Jan2022",as.character(plot_dfH$Month))
plot_dfH$Month<-gsub("0222","Feb2022",as.character(plot_dfH$Month))
plot_dfH$Month<-gsub("0322","March2022",as.character(plot_dfH$Month))
plot_dfH$Month<-gsub("0422","April2022",as.character(plot_dfH$Month))
plot_dfH$Month<-gsub("0522","May2022",as.character(plot_dfH$Month))
plot_dfH$Month<-gsub("0622","June2022",as.character(plot_dfH$Month))
plot_dfH$Month<-gsub("0722","July2022",as.character(plot_dfH$Month))
plot_dfH$Month<-gsub("0822","Aug2022",as.character(plot_dfH$Month))
plot_dfH$Month<-gsub("0922","Sept2022",as.character(plot_dfH$Month))
plot_dfH$Month<-gsub("1022","Oct2022",as.character(plot_dfH$Month))
plot_dfH$Month<-gsub("1122","Nov2022",as.character(plot_dfH$Month))
plot_dfH$Month<-gsub("1222","Dec2022",as.character(plot_dfH$Month))
plot_dfH$Month<-gsub("1221","Dec2021",as.character(plot_dfH$Month))
plot_dfH$Month<-factor(plot_dfH$Month,
                       levels=c("Dec2021","Jan2022","Feb2022","March2022","April2022","May2022","June2022","July2022","Aug2022","Sept2022","Oct2022","Nov2022","Dec2022"))

plot_H_nmds <- ggplot(plot_dfH, aes(x = NMDS1, y = NMDS2, color = Month, shape = Loc)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = palNMDS) +
  stat_ellipse(aes(color=Loc), linetype=2, level=0.6)+
  theme_classic()
plot_H_nmds
ggsave("H_BactL_NMDS_04162025.tiff", units="in", width=8, height=6, dpi=300, compression = 'lzw')

#_____DIFFERENTIAL ABUNDANCE AT THE GENUS LEVEL____________________________________
#### Calculate relative abundance from rarefied phyloseq object
H.rarefied_noFL <- subset_samples(H,Loc!= "FL")
H.rarefied_noFL.counts <- count_seqs(H.rarefied_noFL@otu_table)

abund.H.rarefied_noFL <- H.rarefied_noFL@otu_table %>% as.matrix() %>% as.data.frame() %>%
  rownames_to_column("ASV_") %>%
  pivot_longer(-ASV_, names_to = "sample") %>%
  merge(H.rarefied_noFL.counts$counts, by="sample") %>%
  dplyr::select(-n) %>%
  mutate(newcol=sample)%>%
  separate(newcol, c("X16S","Month", "Lot","Day","Rep","Loc","Actual_Day"), "_",remove = TRUE)%>%
  group_by(sample) %>%
  mutate(rel_abund = value/counts) %>% #relative abundance by treatment/day/rep
  ungroup() 

#check that relative abundance is calculated correctly
abund.H.rarefied_noFL %>% group_by(sample) %>% summarise(sum = sum(rel_abund))

#get taxonomy
taxH.rarefied_noFL <- H.rarefied_noFL@tax_table %>%
  as.data.frame() %>%
  rownames_to_column("ASV_") 
head(taxH.rarefied_noFL)

#merge taxonomy with rarefied relative abundance
count_tax.H.rarefied_noFL <- merge(abund.H.rarefied_noFL, taxH, by="ASV_")
head(count_tax.H.rarefied_noFL)

# Get summary abundance of each taxa
taxa.meta.sumH.rarefied_noFL <- count_tax.H.rarefied_noFL%>%
  dplyr::select(-c(value, counts)) %>%
  pivot_longer(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "ASV"),
               names_to = "level", 
               values_to = "taxon") 
#format phylum names
genus_rel_abundH.rarefied_noFL <- taxa.meta.sumH.rarefied_noFL %>%
  filter(level == "Genus") %>%
  group_by(taxon, sample,Month, Lot,Day,Rep,Loc,Actual_Day) %>% 
  summarise(rel_abund = 100*sum(rel_abund), .groups="drop") %>%
  mutate(taxon = str_replace(taxon,
                             "^unclassified (.*)", "Unclassified *\\1*"),
         taxon = str_replace(taxon, "^(\\S*)$", "\\1")) %>%
  ungroup()
#REPLACE 'NA' OF UNCLASSIFIED TAXA WITH 'UNCLASSIFIED'
genus_rel_abundH.rarefied_noFL$taxon <-genus_rel_abundH.rarefied_noFL$taxon%>%replace_na('Unclassified')
#DCAST TO CONVERT THE DATA FRAME
genus_rel_abundH.rarefied_noFL_4maaslin<-dcast(genus_rel_abundH.rarefied_noFL,sample~taxon,value.var="rel_abund")
genus_rel_abundH.rarefied_noFL_4maaslin<-genus_rel_abundH.rarefied_noFL_4maaslin%>% column_to_rownames(var="sample")

#make taxa and sample data tables from phyloseq object
H.rarefied_noFL_TAX<-data.frame(H.rarefied_noFL@tax_table)
H.rarefied_noFL_sampledata<-data.frame(H.rarefied_noFL@sam_data)
H.rarefied_noFL_sampledata$sampling_code = paste(H.rarefied_noFL_sampledata$Month, H.rarefied_noFL_sampledata$Lot, sep="_")

#check to make sure all are actually data.frames
class(genus_rel_abundH.rarefied_noFL_4maaslin)
class(H.rarefied_noFL_sampledata)
class(H.rarefied_noFL_TAX)

#run maaslin2 with relative abundance at the genus level

Hrare_noFL_fit_2 = Maaslin2(input_data = genus_rel_abundH.rarefied_noFL_4maaslin, 
                          input_metadata = H.rarefied_noFL_sampledata,
                          analysis_method = "CPLM",
                          normalization  = "NONE", 
                          transform="NONE",
                          min_prevalence = 0,
                          min_abundance = 0,
                          save_models="TRUE",
                          save_scatter="TRUE",
                          random_effects = c("sampling_code"),
                          fixed_effects  = "Loc",
                          reference      = c("Loc,CA"),
                          output         = "/local1/workdir1/tw488/CIDA/nus_output_11172024",
                          cores=2) 

#Import results in and filter by a <0.05 FDR
H.rare_noFL_sig_results<-read.table("/local1/workdir1/tw488/CIDA/16S_Maaslin2_Output_11172024/H_only_genus_output_11172024/all_results.tsv",header=TRUE)
H.rare_noFL_sig_results0.05<-H.rare_noFL_sig_results%>%filter(qval<0.05)%>%
  select("feature", "value", "coef","qval")%>% mutate(log2_fold_change = log2(exp(coef)))%>%select("feature", "value","qval","log2_fold_change")

#sort out which genera are differential abundance in AZ vs CA
#differential abundance in AZ 
H.rare_noFL_sig_results0.05_pos<-H.rare_noFL_sig_results%>%filter("log2_fold_change">0)
#differential abundance in  CA
H.rare_noFL_sig_results0.05_neg<-H.rare_noFL_sig_results%>%filter("log2_fold_change"<0)

#get genera out in to a list
H.rare_noFL_sig_results0.05_comb<-c(H.rare_noFL_sig_results0.05$feature)
print(H.rare_noFL_sig_results0.05_comb)

#______________________________________________________________________________________-
#use maaslin2 at ASV level for H samples
H.rarefied_noFL_OTU<-data.frame(H.rarefied_noFL@otu_table)
Hrarefit_noFL_ASVs_11172024 = Maaslin2(input_data = H.rarefied_noFL_OTU, 
                                       input_metadata = H.rarefied_noFL_sampledata,
                                       analysis_method = "CPLM",
                                       min_prevalence = 0.1,
                                       min_abundance = 0.01,
                                       transform="NONE",
                                       normalization  = "NONE",
                                       random_effects = c("sampling_code"),
                                       fixed_effects  = "Loc",
                                       reference      = c("Loc,CA"),
                                       output         = "/local1/workdir1/tw488/CIDA/16S_Maaslin2_Output_11172024/H_only_ASV_output_11172024",
                                       cores=2) 

#Import results in and filter by a <0.05 FDR
H_ASV_sig_results<-read.table("/local1/workdir1/tw488/CIDA/16S_Maaslin2_Output_11172024/H_only_ASV_output_11172024/all_results.tsv",header=TRUE)
H_ASV_sig_results0.05<-H_ASV_sig_results%>%filter(qval<0.05)%>%
  select("feature", "value", "coef","qval")%>% mutate(log2_fold_change = log2(exp(coef)))%>%select("feature", "value","qval","log2_fold_change")
names(H_ASV_sig_results0.05)[1]<-paste("ASV") 
#merge ASV table with significant results
H_ASV_sig_results0.05<-H_ASV_sig_results0.05%>% full_join(H.rarefied_noFL_TAX, by = "ASV")%>%
  filter(!is.na(log2_fold_change))%>%
  select("ASV", "value", "log2_fold_change","qval","Family","Genus")


#CODE TO MAKE Harvest ASV and Genus results tables________________________________________________________

#get avg rel abund by genus and location for Genus differential abundance results

genus_rel_abundH.rarefied_noFL_no0_taxaMEANbyLOC<-genus_rel_abundH.rarefied_noFL %>% group_by(taxon,Loc) %>% summarise(AvgbyLoc= mean(rel_abund))
genus_rel_abundH.rarefied_noFL_no0_taxaMEANbyLOC1<-genus_rel_abundH.rarefied_noFL_no0_taxaMEANbyLOC%>%filter(taxon %in% H.rare_noFL_sig_results0.05_comb)
genus_rel_abundH.rarefied_noFL_no0_taxaMEANbyLOC_allclostridiumandANPR<-genus_rel_abundH.rarefied_noFL_no0_taxaMEANbyLOC%>%filter(taxon %in% c('Clostridium sensu stricto 8', 'Clostridium sensu stricto 1', 'Clostridium sensu stricto 10'))
genus_rel_abundH.rarefied_noFL_no0_taxaMEANbyLOC1<-dcast(genus_rel_abundH.rarefied_noFL_no0_taxaMEANbyLOC1,taxon~Loc,value.var="AvgbyLoc")
genus_rel_abundH.rarefied_noFL_no0_taxaMEANbyLOC_allclostridiumandANPR<-dcast(genus_rel_abundH.rarefied_noFL_no0_taxaMEANbyLOC_allclostridiumandANPR,taxon~Loc,value.var="AvgbyLoc")
genus_rel_abundH.rarefied_noFL_no0_taxaMEANbyLOC4export<-rbind(genus_rel_abundH.rarefied_noFL_no0_taxaMEANbyLOC1,genus_rel_abundH.rarefied_noFL_no0_taxaMEANbyLOC_allclostridiumandANPR)

#merge them all
names(H.rare_noFL_sig_results0.05)[1]<-paste("Genus") 
names(genus_rel_abundH.rarefied_noFL_no0_taxaMEANbyLOC4export)[1]<-paste("Genus") 
genus_rel_abundH.rarefied_noFL_no0_taxaMEANbyLOC4export<-merge(genus_rel_abundH.rarefied_noFL_no0_taxaMEANbyLOC4export,H.rare_noFL_sig_results0.05, by="Genus",all.x = TRUE,all.y = TRUE)
names(genus_rel_abundH.rarefied_noFL_no0_taxaMEANbyLOC4export)[2]<-paste("Average Relative Abundance (%) in Yuma, AZ area")
names(genus_rel_abundH.rarefied_noFL_no0_taxaMEANbyLOC4export)[3]<-paste("Salinas, CA area")
write.csv(genus_rel_abundH.rarefied_noFL_no0_taxaMEANbyLOC4export,"/local1/workdir1/tw488/CIDA/16S_Maaslin2_Output_11172024/H_only_genus_output_11172024/H_avg_rel_abund_by_genus_and_location", row.names = TRUE)



#_____16S analysis by days from start of the planting season and between areas___________________________________________________________________

#TFSOS= day from start of season (PLANTING SEASON)
#calculate days from start of planting season (ie when the fist seeds are to be planted in a growing area)
TFSOSv2<-read.csv("/local1/workdir1/tw488/CIDA_Sequencing_Runs/DFSOSv2.csv", header=TRUE)
TFSOSv2$first_planting_date_from_TF<-gsub("_", "-", TFSOSv2$first_planting_date_from_TF)
TFSOSv2$date_of_harvest_yyyy_mm_dd<-gsub("_", "-", TFSOSv2$date_of_harvest_yyyy_mm_dd)
TFSOSv2$first_planting_date_from_TF<-as.Date(TFSOSv2$first_planting_date_from_TF)
TFSOSv2$date_of_harvest_yyyy_mm_dd<-as.Date(TFSOSv2$date_of_harvest_yyyy_mm_dd)
TFSOSv2$time_from_first_planting<-difftime(TFSOSv2$date_of_harvest_yyyy_mm_dd,TFSOSv2$first_planting_date_from_TF,units = "days")
TFSOSv2$time_from_first_planting<-as.integer(TFSOSv2$time_from_first_planting)
TFSOSv2<-TFSOSv2%>%filter(location!="FL"&location!="GA")%>%filter(sample_type=="r")
TFSOSv2<-TFSOSv2[c("sampling_code","time_from_first_planting")]

#ALPHA DIVERSITY INDICES AND TFSOS
Halphadiv<-estimate_richness(H, split = TRUE, measures = NULL)
Halphadiv$sample <- row.names(Halphadiv)  
Halphadiv <- Halphadiv %>% separate(sample, c("x16s", "month", "lot", "day", "rep", "loc","r1"), "_")
Halphadiv$sample <- row.names(Halphadiv) 
Halphadiv$sampling_code = paste(Halphadiv$month, Halphadiv$lot, sep="_")
Halphadiv <- Halphadiv %>% 
  as.data.frame() %>% 
  full_join(TFSOSv2, by = "sampling_code")
Halphadiv1 <- Halphadiv%>% 
  filter(complete.cases(.))

#calculate Pielou's Eveness
Halphadiv1 <- Halphadiv1 %>% mutate (Pielou=(Halphadiv1$Shannon)/log(Halphadiv1$Observed))

#Visualize
Halphadiv1 %>% ggplot(aes(x=time_from_first_planting, y=Observed)) +
  geom_jitter(size=0.25) +
  labs(x="Day from Start of Season",
       y="Number of Observed ASVs") +
  scale_x_continuous() +
  guides(color = guide_legend(override.aes = list(size=1))) +
  theme_classic()

Halphadiv1 %>% ggplot(aes(x=time_from_first_planting, y=Pielou)) +
  geom_jitter(size=0.25) +
  labs(x="Day from Start of Season",
       y="Pielou's Evenness") +
  scale_x_continuous() +
  guides(color = guide_legend(override.aes = list(size=1))) +
  theme_classic()

Halphadiv1 %>% ggplot(aes(x=time_from_first_planting, y=Shannon)) +
  geom_jitter(size=0.25) +
  labs(x="Day from Start of Season",
       y="Shannon Index") +
  scale_x_continuous() +
  guides(color = guide_legend(override.aes = list(size=1))) +
  theme_classic()

#SPEARMAN CORRELATION TEST TO CORRELATE TIME FROM THE BEGINING OF THE PLANTING SEASON AND SHANNON VALUES
mtrxcorrShn<-Halphadiv1%>% 
  select(time_from_first_planting, Shannon)
#CORRELATION TESTING WITH P VALUES: need to average replicates because the calculation needs one observation per unique metadata label
mtrxcorrShn1<-mtrxcorrShn%>%group_by(time_from_first_planting)%>%mutate(avgSh=mean(Shannon))%>%distinct(time_from_first_planting,avgSh)
cor.test(mtrxcorrShn1$time_from_first_planting,mtrxcorrShn1$avgSh,method=("spearman"))

#SPEARMAN CORRELATION TEST TO CORRELATE TIME FROM THE BEGINING OF THE PLANTING SEASON AND PIELOU VALUES
mtrxcorrPielou<-Halphadiv1%>% 
  select(time_from_first_planting, Pielou) 
#CORRELATION TESTING WITH P VALUES: need to average replicates because the calculation needs one observation per unique metadata label
mtrxcorrPielou1<-mtrxcorrPielou%>%group_by(time_from_first_planting)%>%mutate(avgPielou=mean(Pielou))%>%distinct(time_from_first_planting,avgPielou)
cor.test(mtrxcorrPielou1$time_from_first_planting,mtrxcorrPielou1$avgPielou,method=("spearman"))

#SPEARMAN CORRELATION TEST TO CORRELATE TIME FROM THE BEGINING OF THE PLANTING SEASON AND RICHNESS
mtrxcorrObs<-Halphadiv1%>% 
  select(time_from_first_planting,Observed) 
#CORRELATION TESTING WITH P VALUES: need to average replicates because the calculation needs one observation per unique metadata label
mtrxcorrObs1<-mtrxcorrObs%>%group_by(time_from_first_planting)%>%mutate(avgObs=mean(Observed))%>%distinct(time_from_first_planting,avgObs)
cor.test(mtrxcorrObs1$time_from_first_planting,mtrxcorrObs1$avgObs,method=("spearman"))

#alpha diversity by area

ggqqplot(Halphadiv1$Pielou)
ggqqplot(Halphadiv1$Shannon)
ggqqplot(Halphadiv1$Observed)
ggdensity(Halphadiv1$Shannon)
ggdensity(Halphadiv1$Observed)
ggdensity(Halphadiv1$Pielou)

H_loc_shan_ttest<-t.test(Shannon ~ loc, data = Halphadiv1)
H_loc_shan_ttest
H_loc_pielou_ttest<-t.test(Pielou ~ loc, data = Halphadiv1)
H_loc_pielou_ttest
H_loc_obs_ttest<-t.test(Observed ~ loc, data = Halphadiv1)
H_loc_obs_ttest

#____which genera are over-represented in a specific time of planting season? (maybe specific by region)___________________________________________________________________________________________

#add DFSOS to data frames for Masslin2 analysis

H.rarefied_noFLsampledata<-data.frame(H.rarefied_noFL@sam_data)
H.rarefied_noFLsampledata$sampling_code = paste(H.rarefied_noFLsampledata$Month, H.rarefied_noFLsampledata$Lot, sep="_")
H.rarefied_noFLsampledata_dfsos <- H.rarefied_noFLsampledata %>%
tibble::rownames_to_column("rowsample") %>% 
inner_join(TFSOSv2, by = "sampling_code")
H.rarefied_noFLsampledata_dfsos <-H.rarefied_noFLsampledata_dfsos%>%column_to_rownames( var = "rowsample")

#run maaslin2 w/ updated data frames

H_LOC_TFSOS_11172024 = Maaslin2(input_data = genus_rel_abundH.rarefied_noFL_4maaslin, 
                                   input_metadata = H.rarefied_noFLsampledata_dfsos,
                                   analysis_method = "CPLM",
                                   transform= "NONE",
                                   min_prevalence = 0,
                                   min_abundance = 0,
                                   normalization  = "NONE",
                                   random_effects = c("sampling_code"),
                                   output         = "/local1/workdir1/tw488/CIDA/16S_Maaslin2_Output_11172024/H_LOC_DFSOS_11172024", 
                                   fixed_effects  = c("Loc","time_from_first_planting"),
                                   reference      = c("Loc,CA"),
                                max_pngs = 200,
                                save_models="TRUE",
                                save_scatter="TRUE",
                                cores=2)

#import and filter results
Hrarefit_output_LocTFSOSresults<-read.table("/local1/workdir1/tw488/CIDA/16S_Maaslin2_Output_11172024/H_LOC_DFSOS_11172024/all_results.tsv",header=TRUE)
Hrarefit_output_LocTFSOSresults0.05<-Hrarefit_output_LocTFSOSresults%>%filter(qval<0.05)%>%
  select("feature", "value", "coef","qval")%>% mutate(log2_fold_change = log2(exp(coef)))%>%select(feature,value,qval,log2_fold_change)

#merge w/ relative abundance data to make a table
#time data, split time in to 1st half and second half of season by dividing month associated with a location
tffh_pg4<-Hrarefit_output_LocTFSOSresults0.05%>%filter(value=="time_from_first_planting")
colnames(tffh_pg4)[3] <- "FDR"
tffh_taxon<-c(tffh_pg4$feature)
colnames(genus_rel_abundH.rarefied_noFL)[2] <- "SampleName"
#filter out unneeded taxa
genus_rel_abundH.rarefied_noFL_sampName<-genus_rel_abundH.rarefied_noFL%>%filter(taxon %in% tffh_taxon)
#filter out Clostridium.sensu.stricto.1 and Clostridium.sensu.stricto.13 and OM27.clade
genus_rel_abundH.rarefied_noFL_Clostridium.sensu.stricto.1<- genus_rel_abundH.rarefied_noFL%>%filter(taxon=="Clostridium sensu stricto 1")
genus_rel_abundH.rarefied_noFL_Clostridium.sensu.stricto.13<- genus_rel_abundH.rarefied_noFL%>%filter(taxon=="Clostridium sensu stricto 13")
genus_rel_abundH.rarefied_noFL_Clostridium.sensu.stricto.8<- genus_rel_abundH.rarefied_noFL%>%filter(taxon=="Clostridium sensu stricto 8")
genus_rel_abundH.rarefied_noFL_OM27.clade<- genus_rel_abundH.rarefied_noFL%>%filter(taxon=="OM27 clade")
genus_rel_abundH.rarefied_noFL_Clostridium.sensu.stricto<-rbind(genus_rel_abundH.rarefied_noFL_Clostridium.sensu.stricto.1,genus_rel_abundH.rarefied_noFL_Clostridium.sensu.stricto.13)
genus_rel_abundH.rarefied_noFL_Clostridium.sensu.stricto<-rbind(genus_rel_abundH.rarefied_noFL_Clostridium.sensu.stricto,genus_rel_abundH.rarefied_noFL_OM27.clade)
genus_rel_abundH.rarefied_noFL_Clostridium.sensu.stricto<-rbind(genus_rel_abundH.rarefied_noFL_Clostridium.sensu.stricto,genus_rel_abundH.rarefied_noFL_Clostridium.sensu.stricto.8)
genus_rel_abundH.rarefied_noFL_sampName<-rbind(genus_rel_abundH.rarefied_noFL_sampName,genus_rel_abundH.rarefied_noFL_Clostridium.sensu.stricto)

#separate in to CA and AZ
genus_rel_abundH.rarefied_noFL_sampNameCA<-genus_rel_abundH.rarefied_noFL_sampName%>%filter(Loc=="CA")%>%mutate(whichhalfofszn="First")
genus_rel_abundH.rarefied_noFL_sampNameAZ<-genus_rel_abundH.rarefied_noFL_sampName%>%filter(Loc=="AZ")%>%mutate(whichhalfofszn="First")

#make monthlot
genus_rel_abundH.rarefied_noFL_sampNameCA$sampling_code = paste(genus_rel_abundH.rarefied_noFL_sampNameCA$Month,genus_rel_abundH.rarefied_noFL_sampNameCA$Lot, sep="_")
genus_rel_abundH.rarefied_noFL_sampNameAZ$sampling_code = paste(genus_rel_abundH.rarefied_noFL_sampNameAZ$Month,genus_rel_abundH.rarefied_noFL_sampNameAZ$Lot, sep="_")

#add in first or second in whichhalfofszn column
genus_rel_abundH.rarefied_noFL_sampNameCA$whichhalfofszn[genus_rel_abundH.rarefied_noFL_sampNameCA$sampling_code == "0722_B"] <- "Second"
genus_rel_abundH.rarefied_noFL_sampNameCA$whichhalfofszn[genus_rel_abundH.rarefied_noFL_sampNameCA$Month == "0822"] <- "Second"
genus_rel_abundH.rarefied_noFL_sampNameCA$whichhalfofszn[genus_rel_abundH.rarefied_noFL_sampNameCA$Month == "0922"] <- "Second"
genus_rel_abundH.rarefied_noFL_sampNameCA$whichhalfofszn[genus_rel_abundH.rarefied_noFL_sampNameCA$Month == "1022"] <- "Second"

#get average relative abundance bywhichhalfofszn
#genus_rel_abundH.rarefied_noFL_sampNameCA1<-genus_rel_abundH.rarefied_noFL_sampNameCA%>%group_by(taxon,whichhalfofszn)%>%mutate(AvybyWhos=mean(rel_abund))%>%distinct(taxon,whichhalfofszn,AvybyWhos)
#genus_rel_abundH.rarefied_noFL_sampNameCA1<-dcast(genus_rel_abundH.rarefied_noFL_sampNameCA1,taxon~whichhalfofszn,value.var="AvybyWhos")

#now AZ 
genus_rel_abundH.rarefied_noFL_sampNameAZ$whichhalfofszn[genus_rel_abundH.rarefied_noFL_sampNameAZ$sampling_code == "0222_B"] <- "Second"
genus_rel_abundH.rarefied_noFL_sampNameAZ$whichhalfofszn[genus_rel_abundH.rarefied_noFL_sampNameAZ$Month == "0322"] <- "Second"
genus_rel_abundH.rarefied_noFL_sampNameAZ$whichhalfofszn[genus_rel_abundH.rarefied_noFL_sampNameAZ$Month == "0422"] <- "Second"

#get average relative abundance bywhichhalfofszn
#genus_rel_abundH.rarefied_noFL_sampNameAZ1<-genus_rel_abundH.rarefied_noFL_sampNameAZ%>%group_by(taxon,whichhalfofszn)%>%mutate(AvybyWhos=mean(rel_abund))%>%distinct(taxon,whichhalfofszn,AvybyWhos)
#genus_rel_abundH.rarefied_noFL_sampNameAZ1<-dcast(genus_rel_abundH.rarefied_noFL_sampNameAZ1,taxon~whichhalfofszn,value.var="AvybyWhos")

#average AZ and CA together
genus_rel_abundH.rarefied_noFL_sampNameAZCA<-rbind(genus_rel_abundH.rarefied_noFL_sampNameCA,genus_rel_abundH.rarefied_noFL_sampNameAZ)
colnames(genus_rel_abundH.rarefied_noFL_sampNameAZCA)[1] <- "Genus"
genus_rel_abundH.rarefied_noFL_sampNameAZCA<-genus_rel_abundH.rarefied_noFL_sampNameAZCA%>%group_by(Genus,whichhalfofszn)%>%mutate(AvgRelAb=mean(rel_abund))%>%distinct(Genus,whichhalfofszn,AvgRelAb)
genus_rel_abundH.rarefied_noFL_sampNameAZCA<-dcast(genus_rel_abundH.rarefied_noFL_sampNameAZCA,Genus~whichhalfofszn,value.var="AvgRelAb")

#merge em all
colnames(tffh_pg4)[1] <- "Genus"
tffh_pg4_1<-merge(genus_rel_abundH.rarefied_noFL_sampNameAZCA,tffh_pg4, by="Genus",all.x = TRUE,all.y = TRUE)#%>%mutate(across(where(is.numeric), ~replace_na(., 0)))
colnames(tffh_pg4_1)[2] <- "Relative Abundance(%) for First Half of Planting Season"
colnames(tffh_pg4_1)[3] <- "Relative Abundance(%) for Second Half of Planting Season"
write.csv(tffh_pg4_1,"/local1/workdir1/tw488/CIDA/16S_Maaslin2_Output_11172024/H_LOC_DFSOS_11172024/Time_from_start_of_season_Maaslin_results_with_RelAb", row.names = TRUE)

#Area data, split time in to 1st half and second half of season___________________________________
Loc_tfsosandloc<-Hrarefit_output_LocTFSOSresults0.05%>%filter(value=="AZ")
colnames(Loc_tfsosandloc)[1] <- "taxon"
colnames(Loc_tfsosandloc)[3] <- "FDR"
genus_rel_abundH.rarefied_noFL_no0_taxaMEANbyLOC_loc1<-dcast(genus_rel_abundH.rarefied_noFL_no0_taxaMEANbyLOC,taxon~Loc,value.var="AvgbyLoc")
Loc_tfsosandloc$taxon<-gsub("Clostridium.sensu.stricto.1", "Clostridium sensu stricto 1", Loc_tfsosandloc$taxon)
Loc_tfsosandloc$taxon<-gsub("Clostridium.sensu.stricto.10", "Clostridium sensu stricto 10", Loc_tfsosandloc$taxon)
Loc_tfsosandloc$taxon<-gsub("Clostridium.sensu.stricto.8", "Clostridium sensu stricto 8", Loc_tfsosandloc$taxon)
Loc_tfsosandloc$taxon<-gsub("Clostridium.sensu.stricto.13", "Clostridium sensu stricto 13", Loc_tfsosandloc$taxon)
Loc_tfsosandloc_v2<-Loc_tfsosandloc%>%right_join(genus_rel_abundH.rarefied_noFL_no0_taxaMEANbyLOC_loc1, by='taxon')%>%na.omit()
colnames(Loc_tfsosandloc_v2)[5] <- "Relative Abundance (%) in Yuma, AZ area"
colnames(Loc_tfsosandloc_v2)[6] <- "Salinas, CA area"
Loc_tfsosandloc_v2 <- Loc_tfsosandloc_v2[, -2] 
colnames(Loc_tfsosandloc_v2)[1] <- "Genus"
colnames(Loc_tfsosandloc_v2)[2] <- "Significance Level"
write.csv(Loc_tfsosandloc_v2,"/local1/workdir1/tw488/CIDA/16S_Maaslin2_Output_11172024/H_LOC_DFSOS_11172024/Location_Maaslin_results_with_RelAb", row.names = TRUE)

#___HARVEST VS DAY INITIAL SECTION___________________________________________________________________________________________
#how many D samples?
DIonly<-data.frame(HDI@sam_data)
DIonly<-DIonly%>%filter(Day=="DI")

#PERMANOVA WITH HOW FUNCTION
HDIotu<-pstoveg_otu(HDI)
HDIsd<-pstoveg_sd(HDI)
HDIsd$Monthlot = paste(HDIsd$Month, HDIsd$Lot, sep="_")
HDIvegdist<-vegdist(HDIotu, method="bray") 
table(HDIsd$Monthlot)
#remove 0822_B,0922_B, 1122_A, 1222_B and C
HDIsd1<-HDIsd%>%filter(Monthlot!="0822_B")
HDIsd1<-HDIsd1%>%filter(Monthlot!="0922_B")
HDIsd1<-HDIsd1%>%filter(Monthlot!="1122_A")
HDIsd1<-HDIsd1%>%filter(Monthlot!="1222_B")
HDIsd1<-HDIsd1%>%filter(Monthlot!="1222_C")
HDIsd1<-HDIsd1%>%filter(Monthlot!="0122_A")
table(HDIsd1$Monthlot)

#filter out above monthlots from otu table
HDIotu1<-data.frame(HDIotu)%>%tibble::rownames_to_column("samples")%>%mutate(SampleNames=samples)
HDIotu1$SampleNames<-gsub("_R1", "", HDIotu1$SampleNames)
str(HDIotu1)
HDIotu1 <- HDIotu1 %>%mutate(sample=samples)%>% separate(sample, c("x16s", "month", "lot", "day", "rep", "loc","r1"), "_")
HDIotu1$sampling_code = paste(HDIotu1$month,HDIotu1$lot, sep="_")
HDIotu1<-HDIotu1%>%filter(sampling_code!="0822_B")
HDIotu1<-HDIotu1%>%filter(sampling_code!="0922_B")
HDIotu1<-HDIotu1%>%filter(sampling_code!="1122_A")
HDIotu1<-HDIotu1%>%filter(sampling_code!="1222_B")
HDIotu1<-HDIotu1%>%filter(sampling_code!="1222_C")
HDIotu1<-HDIotu1%>%filter(sampling_code!="0122_A")
HDIotu1<-HDIotu1%>%column_to_rownames(var="samples") %>%data.matrix
str(HDIotu1)
HDIvegdist1<-vegdist(HDIotu1, method="bray")

#permanova
table(HDIsd1$Monthlot)
h4 <- with(HDIsd1, how(within=Within(type="free"),plots=Plots(strata=Monthlot, type="free"),nperm=999))
h4

HDIvegdist_perm <- adonis2(HDIvegdist1~Day*Loc, data=HDIsd1, permutations = h4, method = "bray", 
                           sqrt.dist = FALSE, add = FALSE, by = "terms",
                           parallel = getOption("mc.cores"), na.action = na.fail)
HDIvegdist_perm

#NMDS FOR H-DI SAMPLES
HDI_NMDS <- metaMDS(HDIvegdist)
stressplot(HDI_NMDS)
HDIsdNMDS<-HDIsd%>%rownames_to_column("sample")
plot_HDI_NMDS <- scores(HDI_NMDS, display = "sites") %>% 
  as.data.frame() %>% 
  rownames_to_column("sample") %>% 
  full_join(HDIsdNMDS, by = "sample")


plot_HDI_NMDS$Day<-factor(plot_HDI_NMDS$Day,
                       levels=c("H","DI"))


palNMDS_HDI <- c( "#1F78B4", "#B2DF8A","black","grey","grey39")

plot_HDI_NMDS_ <- ggplot(plot_HDI_NMDS, aes(x = NMDS1, y = NMDS2, color = Day, shape = Loc)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = palNMDS_HDI) +
  stat_ellipse(aes(color=Loc), linetype=2, level=0.6)+
  theme_classic()
plot_HDI_NMDS_
ggsave("HDI_Bact_NMDS_04262025.tiff", units="in", width=8, height=6, dpi=300, compression = 'lzw')

# alpha diversity for H-DI______________________________________________________

#calculate alpha diversity stats
HDIrichnessstats<-estimate_richness(HDI, split = TRUE, measures = NULL)
HDIrichnessstats$sample <- row.names(HDIrichnessstats)  
HDIrichnessstats <- HDIrichnessstats %>% separate(sample, c("x16s", "month", "lot", "day", "rep", "loc","R2","r1"), "_")
#calculate Pielou
HDIrichnessstats <- HDIrichnessstats %>% mutate (Pielou=(HDIrichnessstats$Shannon)/log(HDIrichnessstats$Observed))
#check for normality
ggqqplot(HDIrichnessstats$Shannon)
ggqqplot(HDIrichnessstats$Pielou)
ggqqplot(HDIrichnessstats$Observed)
ggdensity(HDIrichnessstats$Shannon)
ggdensity(HDIrichnessstats$Pielou)
ggdensity(HDIrichnessstats$Observed)

#ITS P NORMAL SO LETS USE A T TEST
HDI_shan_ttest<-t.test(Shannon ~ day, data = HDIrichnessstats)
HDI_shan_ttest
HDI_Observed_ttest<-t.test(Observed ~ day, data = HDIrichnessstats)
HDI_Observed_ttest
HDI_Pielou_ttest<-t.test(Pielou ~ day, data = HDIrichnessstats)
HDI_Pielou_ttest

#graph them
#shannon
HDIrichnessstats$day<-gsub("DI", "Day Initial", HDIrichnessstats$day)
HDIrichnessstats$day<-gsub("H", "Harvest", HDIrichnessstats$day)

HDIrichnessstats$day<-factor(HDIrichnessstats$day,
                             levels=c("Harvest","Day Initial"))
#visualize
HDIrichnessstats%>% 
  ggplot(aes(x=day, y=Shannon, color=day)) +
  geom_boxplot(size=0.25) +
  labs(x="Day",
       y="Shannon Index") +
  theme_classic()+geom_signif(comparisons=list(c("Harvest","Day Initial")),color="black")

#observed
HDIrichnessstats%>% 
  ggplot(aes(x=day, y=Observed,color=day)) +
  geom_boxplot(size=0.25) +
  labs(x="Day",
       y="Richness") +
  theme_classic()+geom_signif(comparisons=list(c("Harvest","Day Initial")),color="black")

#Pielou
HDIrichnessstats%>% 
  ggplot(aes(x=day, y=Pielou,color=day)) +
  geom_boxplot(size=0.25) +
  labs(x="Day",
       y="Pielou's Evenness") +
  theme_classic()+geom_signif(comparisons=list(c("Harvest","Day Initial")),color="black")

# Find most common genera (TOP 5) in DI samples

#### Calculate relative abundance ####
DI <- subset_samples(HDI, Day!= "H")
DI.counts <- count_seqs(DI@otu_table)

abund.DI <- DI@otu_table %>% as.matrix() %>% as.data.frame() %>%
  rownames_to_column("OTU") %>%
  pivot_longer(-OTU, names_to = "sample") %>%
  merge(DI.counts$counts, by="sample") %>%
  dplyr::select(-n) %>%
  mutate(newcol=sample)%>%
  separate(newcol, c("X16S","Month", "Lot","Day","Rep","Loc","Actual_Day"), "_",remove = TRUE)%>%
  group_by(sample) %>%
  mutate(rel_abund = value/counts) %>% #relative abundance by treatment/day/rep
  ungroup() 

#check that relative abundance is calculated correctly
abund.DI %>% group_by(sample) %>% summarise(sum = sum(rel_abund))

#get taxonomy
taxDI <- DI@tax_table %>%
  as.data.frame() %>%
  rownames_to_column("OTU") 
head(taxDI)

#how many otus are there total?
nrow(taxDI)

#merge taxonomy with rarefied relative abundance
count_tax.DI <- merge(abund.DI, taxDI, by="OTU")
head(count_tax.DI)

#### Relative abundance plot ####

# Get summary abundance of each taxa
taxa.meta.sumDI <- count_tax.DI%>%
  dplyr::select(-c(value, counts)) %>%
  pivot_longer(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "OTU"),
               names_to = "level", 
               values_to = "taxon") 
#format phylum names
genus_rel_abundDI <- taxa.meta.sumDI %>%
  filter(level == "Genus") %>%
  group_by(taxon, sample,Month, Lot,Day,Rep,Loc,Actual_Day) %>% 
  summarise(rel_abund = 100*sum(rel_abund), .groups="drop") %>%
  mutate(taxon = str_replace(taxon,
                             "^unclassified (.*)", "Unclassified *\\1*"),
         taxon = str_replace(taxon, "^(\\S*)$", "\\1")) %>%
  ungroup()
#how many DI genera are there total?
totalDIgenera_before5<-genus_rel_abundDI%>%group_by(taxon)%>%filter(rel_abund!=0)%>%distinct(taxon)

#how many DI genera are there >5% rel ab?
genus_rel_abundDI_AVG <- genus_rel_abundDI%>%group_by(taxon) %>%mutate(avgrelab=mean(rel_abund))%>%distinct(taxon,avgrelab)%>%filter(avgrelab!=0)%>%filter(avgrelab>5)

#____DIFFERNTIAL ABUNDANCE BETWEEN HARVEST AND DAY INITIAL SAMPLES_________________________________________________
#remove FL sampels for maaslin2
HDI.rarefied_4maaslin <- subset_samples(HDI, Month!= "1221")

#make tax and sample data files
HDI.rarefied_4maaslin_TAX<-data.frame(HDI.rarefied_4maaslin@tax_table)
HDI.rarefied_4maaslin_sampledata<-data.frame(HDI.rarefied_4maaslin@sam_data)
HDI.rarefied_4maaslin_sampledata$sampling_code = paste(HDI.rarefied_4maaslin_sampledata$Month, HDI.rarefied_4maaslin_sampledata$Lot, sep="_")

#make genus level otu table
HDI.rarefied_4maaslin.counts <- count_seqs(HDI.rarefied_4maaslin@otu_table)

abund.HDI.rarefied_4maaslin <- HDI.rarefied_4maaslin@otu_table %>% as.matrix() %>% as.data.frame() %>%
  rownames_to_column("ASV_") %>%
  pivot_longer(-ASV_, names_to = "sample") %>%
  merge(HDI.rarefied_4maaslin.counts$counts, by="sample") %>%
  dplyr::select(-n) %>%
  mutate(newcol=sample)%>%
  separate(newcol, c("X16S","Month", "Lot","Day","Rep","Loc","Actual_Day"), "_",remove = TRUE)%>%
  group_by(sample) %>%
  mutate(rel_abund = value/counts) %>% #relative abundance by treatment/day/rep
  ungroup() 

#check that relative abundance is calculated correctly
abund.HDI.rarefied_4maaslin %>% group_by(sample) %>% summarise(sum = sum(rel_abund))

#get taxonomy
taxHDI.rarefied_4maaslin <- HDI.rarefied_4maaslin@tax_table %>%
  as.data.frame() %>%
  rownames_to_column("ASV_") 
head(taxHDI.rarefied_4maaslin)

#merge taxonomy with rarefied relative abundance
count_tax.HDI.rarefied_4maaslin <- merge(abund.HDI.rarefied_4maaslin, taxH, by="ASV_")
head(count_tax.HDI.rarefied_4maaslin)

# Get summary abundance of each taxa
taxa.meta.sumHDI.rarefied_4maaslin <- count_tax.HDI.rarefied_4maaslin%>%
  dplyr::select(-c(value, counts)) %>%
  pivot_longer(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "ASV"),
               names_to = "level", 
               values_to = "taxon") 

#format phylum names
genus_rel_abundHDI.rarefied_4maaslin <- taxa.meta.sumHDI.rarefied_4maaslin %>%
  filter(level == "Genus") %>%
  group_by(taxon, sample,Month, Lot,Day,Rep,Loc,Actual_Day) %>% 
  summarise(rel_abund = 100*sum(rel_abund), .groups="drop") %>%
  mutate(taxon = str_replace(taxon,
                             "^unclassified (.*)", "Unclassified *\\1*"),
         taxon = str_replace(taxon, "^(\\S*)$", "\\1")) %>%
  ungroup()
#REPLACE UNCLASSIFIED TAXA 'NA' WITH 'UNCLASSSIFIED'
genus_rel_abundHDI.rarefied_4maaslin$taxon <-genus_rel_abundHDI.rarefied_4maaslin$taxon%>%replace_na('Unclassified')
genus_rel_abundHDI.rarefied_4maaslin_nodcast<-genus_rel_abundHDI.rarefied_4maaslin
genus_rel_abundHDI.rarefied_4maaslin<-dcast(genus_rel_abundHDI.rarefied_4maaslin,sample~taxon,value.var="rel_abund")
genus_rel_abundHDI.rarefied_4maaslin<-genus_rel_abundHDI.rarefied_4maaslin%>% column_to_rownames(var="sample")
class(HDI.rarefied_4maaslin_sampledata)

#run maaslin2
HDI_byDAY_genus_rarefit = Maaslin2(input_data = genus_rel_abundHDI.rarefied_4maaslin, 
                                   input_metadata = HDI.rarefied_4maaslin_sampledata,
                                   analysis_method = "CPLM",
                                   transform= "NONE",
                                   min_prevalence = 0,
                                   min_abundance = 0,
                                   normalization  = "NONE",
                                   random_effects = c("sampling_code"),
                                   output         = "/local1/workdir1/tw488/CIDA/16S_Maaslin2_Output_11172024/HDI_byDAY_genus_rarefit_output", 
                                   fixed_effects  = c("Day"),
                                   reference      = c("Day,H")) 

#IMPORT RESULTS AND make a legible data tbl
HDI_byDAY<-read.table("/local1/workdir1/tw488/CIDA/16S_Maaslin2_Output_11172024/HDI_byDAY_genus_rarefit_output/all_results.tsv",header=TRUE)
HDI_byDAY0.05<-HDI_byDAY%>%filter(qval<0.05)
#total table of all enriched genera
HDIvector<-c(HDI_byDAY0.05$feature)

#make table
genus_rel_abundHDI.rarefied_4maaslin_no0_taxaMEANbyDAY<-genus_rel_abundHDI.rarefied_4maaslin_nodcast %>% group_by(taxon,Day) %>% summarise(AvgbyDay= mean(rel_abund))
genus_rel_abundHDI.rarefied_4maaslin_no0_taxaMEANbyDAY1<-genus_rel_abundHDI.rarefied_4maaslin_no0_taxaMEANbyDAY%>%filter(taxon %in% HDIvector)
genus_rel_abundHDI.rarefied_4maaslin_no0_taxaMEANbyDAY1_allothers<-genus_rel_abundHDI.rarefied_4maaslin_no0_taxaMEANbyDAY%>%filter(taxon %in% c('Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium', 'Clostridium sensu stricto 8', 'Clostridium sensu stricto 1', 'Clostridium sensu stricto 10', 'Clostridium sensu stricto 13','Clostridium sensu stricto 7','Candidatus Entotheonella','Candidatus Methylopumilus','hgcI clade','OM27 clade','KCM-B-112','Clostridium sensu stricto 7'))
genus_rel_abundHDI.rarefied_4maaslin_no0_taxaMEANbyDAY1<-dcast(genus_rel_abundHDI.rarefied_4maaslin_no0_taxaMEANbyDAY1,taxon~Day,value.var="AvgbyDay")
genus_rel_abundHDI.rarefied_4maaslin_no0_taxaMEANbyDAY1_allothers<-dcast(genus_rel_abundHDI.rarefied_4maaslin_no0_taxaMEANbyDAY1_allothers,taxon~Day,value.var="AvgbyDay")
genus_rel_abundHDI.rarefied_4maaslin_no0_taxaMEANbyDAY4export<-rbind(genus_rel_abundHDI.rarefied_4maaslin_no0_taxaMEANbyDAY1,genus_rel_abundHDI.rarefied_4maaslin_no0_taxaMEANbyDAY1_allothers)
colnames(HDI_byDAY0.05)[1] <- "taxon"
tabl4_supp<-merge(HDI_byDAY0.05,genus_rel_abundHDI.rarefied_4maaslin_no0_taxaMEANbyDAY4export,by="taxon", all = TRUE)
tabl4_supp<-tabl4_supp %>%
  select("taxon", "coef","qval","DI","H")%>%mutate(log2_fold_change = log2(exp(coef)))
colnames(tabl4_supp)[1] <- "Genus"
colnames(tabl4_supp)[4] <- "Relative Abundance (%) at Day Initial"
colnames(tabl4_supp)[5] <- "Harvest"
colnames(tabl4_supp)[3] <- "Significance Level"
tabl4_supp<-select(tabl4_supp,1,5,4,6,3)
write.csv(tabl4_supp,"/local1/workdir1/tw488/CIDA/16S_Maaslin2_Output_11172024/HDI_byDAY_genus_rarefit_output/tabl4_supp_05232025", row.names = TRUE)

#___Shelf life differential abundance ____________________________________________________________________________________________

#make sample data and taxa frames
DI_D22_D28_noSEUSA <- subset_samples(DI_D22_D28,Loc!= "FL")
DI_D22_D28_noSEUSA <- subset_samples(DI_D22_D28_noSEUSA,Loc!= "GA")
DI_D22_D28_noSEUSA<- subset_samples(DI_D22_D28_noSEUSA,SampleName!="ITS_1222_A_D7_3_AZ_contig.gz")

#check to verify there are no SE USA samples
View(DI_D22_D28_noSEUSA@sam_data)

allShelfLifeandLoc_rarefied_4maaslin_TAX<-data.frame(DI_D22_D28_noSEUSA@tax_table)
allShelfLifeandLoc_rarefied_4maaslin_sampledata<-data.frame(DI_D22_D28_noSEUSA@sam_data)
allShelfLifeandLoc_rarefied_4maaslin_sampledata$sampling_code = paste(allShelfLifeandLoc_rarefied_4maaslin_sampledata$Month, allShelfLifeandLoc_rarefied_4maaslin_sampledata$Lot, sep="_")
#make a separate column based off of Actual_Day to make Day a continuous variable
allShelfLifeandLoc_rarefied_4maaslin_sampledata<-allShelfLifeandLoc_rarefied_4maaslin_sampledata%>%mutate(Actual_Day_noD=Actual_Day)%>%filter(Rep!="NC")
allShelfLifeandLoc_rarefied_4maaslin_sampledata$Actual_Day_noD<-gsub("DI", "1", allShelfLifeandLoc_rarefied_4maaslin_sampledata$Actual_Day_noD)
allShelfLifeandLoc_rarefied_4maaslin_sampledata$Actual_Day_noD<-gsub("D", "", allShelfLifeandLoc_rarefied_4maaslin_sampledata$Actual_Day_noD)
allShelfLifeandLoc_rarefied_4maaslin_sampledata$Actual_Day_noD<-as.integer(allShelfLifeandLoc_rarefied_4maaslin_sampledata$Actual_Day_noD)

#check to make sure 'Actual-Day_noD' is an integer now
str(allShelfLifeandLoc_rarefied_4maaslin_sampledata)

#how many ASVs in shelf life samples?
allShelfLifeandLoc_rarefied_4maaslin.otu_df <- data.frame(DI_D22_D28_noSEUSA@otu_table)
allShelfLifeandLoc_rowsums<-as.data.frame(rowSums(allShelfLifeandLoc_rarefied_4maaslin.otu_df))
colnames(allShelfLifeandLoc_rowsums)[1] <- "rowsums"
allShelfLifeandLoc_rowsums<-allShelfLifeandLoc_rowsums%>%filter(rowsums!=0)

#make genus level otu table
allShelfLifeandLoc_rarefied_4maaslin.counts <- count_seqs(DI_D22_D28_noSEUSA@otu_table)

abund.allShelfLifeandLoc_rarefied_4maaslin <- DI_D22_D28_noSEUSA@otu_table %>% as.matrix() %>% as.data.frame() %>%
  rownames_to_column("ASV_") %>%
  pivot_longer(-ASV_, names_to = "sample") %>%
  merge(allShelfLifeandLoc_rarefied_4maaslin.counts$counts, by="sample") %>%
  dplyr::select(-n) %>%
  mutate(newcol=sample)%>%
  separate(newcol, c("X16S","Month", "Lot","Day","Rep","Loc","Actual_Day"), "_",remove = TRUE)%>%
  group_by(sample) %>%
  mutate(rel_abund = value/counts) %>% #relative abundance by treatment/day/rep
  ungroup() 

#check that relative abundance is calculated correctly
abund.allShelfLifeandLoc_rarefied_4maaslin %>% group_by(sample) %>% summarise(sum = sum(rel_abund))

#get taxonomy
allShelfLifeandLoc_rarefied_4maaslin_tax <- DI_D22_D28_noSEUSA@tax_table %>%
  as.data.frame() %>%
  rownames_to_column("ASV_") 
head(allShelfLifeandLoc_rarefied_4maaslin_tax)

#merge taxonomy with rarefied relative abundance
count_tax_allShelfLifeandLoc_rarefied_4maaslin <- merge(abund.allShelfLifeandLoc_rarefied_4maaslin, allShelfLifeandLoc_rarefied_4maaslin_tax, by="ASV_")
head(count_tax_allShelfLifeandLoc_rarefied_4maaslin)

#### Relative abundance plot ####

# Get summary abundance of each taxa
taxa.meta.sum.allShelfLifeandLoc_rarefied_4maaslin <- count_tax_allShelfLifeandLoc_rarefied_4maaslin%>%
  dplyr::select(-c(value, counts)) %>%
  pivot_longer(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "ASV"),
               names_to = "level", 
               values_to = "taxon") 

#format phylum names
genus_rel_abund_allShelfLifeandLoc_rarefied_4maaslin <- taxa.meta.sum.allShelfLifeandLoc_rarefied_4maaslin %>%
  filter(level == "Genus") %>%
  group_by(taxon, sample,Month, Lot,Day,Rep,Loc,Actual_Day) %>% 
  summarise(rel_abund = 100*sum(rel_abund), .groups="drop") %>%
  mutate(taxon = str_replace(taxon,
                             "^unclassified (.*)", "Unclassified *\\1*"),
         taxon = str_replace(taxon, "^(\\S*)$", "\\1")) %>%
  ungroup()
genus_rel_abund_allShelfLifeandLoc_rarefied_4maaslin$taxon <-genus_rel_abund_allShelfLifeandLoc_rarefied_4maaslin$taxon%>%replace_na('Unclassified')
genus_rel_abund_allShelfLifeandLoc_rarefied_4maaslin_nodcast<-genus_rel_abund_allShelfLifeandLoc_rarefied_4maaslin
genus_rel_abund_allShelfLifeandLoc_rarefied_4maaslin<-dcast(genus_rel_abund_allShelfLifeandLoc_rarefied_4maaslin,sample~taxon,value.var="rel_abund")
genus_rel_abund_allShelfLifeandLoc_rarefied_4maaslin<-genus_rel_abund_allShelfLifeandLoc_rarefied_4maaslin%>% column_to_rownames(var="sample")

#run maaslin

AllShelfLife_Loc_genus_rarefitv2 = Maaslin2(input_data = genus_rel_abund_allShelfLifeandLoc_rarefied_4maaslin, 
                                            input_metadata = allShelfLifeandLoc_rarefied_4maaslin_sampledata,
                                            analysis_method = "CPLM",
                                            transform= "NONE",
                                            min_prevalence = 0,
                                            min_abundance = 0,
                                            normalization  = "NONE",
                                            max_pngs = 200,
                                            save_models="TRUE",
                                            save_scatter="TRUE",
                                            cores=2,
                                            random_effects = c("sampling_code"),
                                            output         = "/local1/workdir1/tw488/CIDA/16S_Maaslin2_Output_11172024/AllShelfLife_Loc_genus_rarefitv2_output", 
                                            fixed_effects  = c("Actual_Day_noD","Loc"),
                                            reference      = c("Loc,CA")) 

#Import results
AllShelfLife_Loc_genusv2<-read.table("/local1/workdir1/tw488/CIDA/16S_Maaslin2_Output_11172024/AllShelfLife_Loc_genus_rarefitv2_output/all_results.tsv",header=TRUE)
AllShelfLife_Loc_genus0.05v2<-AllShelfLife_Loc_genusv2%>%filter(qval<0.05)%>%filter(value=="Actual_Day_noD")%>%
  select("feature", "value", "coef","qval")%>% mutate(log2_fold_change = log2(exp(coef)))
colnames(AllShelfLife_Loc_genus0.05v2)[1] <- "Genus"

#look at top avg genera over shelf life
genus_rel_abund_allShelfLifeandLoc_rarefied_4maaslin_nodcast1<-genus_rel_abund_allShelfLifeandLoc_rarefied_4maaslin_nodcast%>%filter(Rep!="NC")%>%group_by(taxon)%>%mutate(avgrelab=mean(rel_abund))%>%distinct(taxon,avgrelab)
top20_genovershelflife<-head(arrange(genus_rel_abund_allShelfLifeandLoc_rarefied_4maaslin_nodcast1,desc(avgrelab)), n = 21)
top20_genovershelflife
#get avg DI and D22 values for each of the top 20 genera
genus_rel_abund_allShelfLifeandLoc_rarefied_4maaslin_nodcast2<-genus_rel_abund_allShelfLifeandLoc_rarefied_4maaslin_nodcast%>%
  filter(Day!="D7")%>%filter(Day!="D17")%>%
  filter(Day!="D12")%>%group_by(taxon,Day) %>%
  mutate(avgbtwnreps=mean(rel_abund))%>%distinct(taxon,Day,avgbtwnreps)

top20aslist<-top20_genovershelflife$Genus
genus_rel_abund_allShelfLifeandLoc_rarefied_4maaslin_nodcast2_b<-genus_rel_abund_allShelfLifeandLoc_rarefied_4maaslin_nodcast2%>%filter(taxon %in% top20aslist)
genus_rel_abund_allShelfLifeandLoc_rarefied_4maaslin_nodcast3<-dcast(genus_rel_abund_allShelfLifeandLoc_rarefied_4maaslin_nodcast2_b,taxon~Day,value.var="avgbtwnreps")
colnames(genus_rel_abund_allShelfLifeandLoc_rarefied_4maaslin_nodcast3)[1] <- "Genus"

#merge maaslin output and relative abundance
AllShelfLife_df<-merge(x=genus_rel_abund_allShelfLifeandLoc_rarefied_4maaslin_nodcast3, y=AllShelfLife_Loc_genus0.05v2, by='Genus',all.x = TRUE)
colnames(top20_genovershelflife)[1] <- "Genus"
AllShelfLife_df<-merge(x=top20_genovershelflife, y=AllShelfLife_df, by='Genus',all.x = TRUE)
colnames(AllShelfLife_df)[2] <- "Average Relative Abundance Over Shelf Life"
AllShelfLife_df1 <- AllShelfLife_df[,c("Genus","Average Relative Abundance Over Shelf Life", "DI", "D22","log2_fold_change", "qval","coef")]
write.csv(AllShelfLife_df1,"/local1/workdir1/tw488/CIDA/16S_Maaslin2_Output_11172024/AllShelfLife_Loc_genus_rarefitv2_output/AllShelfLife_df1", row.names = TRUE)

#find log fold change of genera that were not sig diff: Erwinia,Paenibacillus,Pseudomonas
top20genera_notsigdiff<-c("Erwinia", "Paenibacillus", "Pseudomonas")

#filter for top 20 genera
AllShelfLife_Loc_genusv2_notsigdiff<-AllShelfLife_Loc_genusv2%>%filter(feature %in% top20genera_notsigdiff)

#filter for 
AllShelfLife_Loc_genusv2_notsigdiff<-AllShelfLife_Loc_genusv2_notsigdiff%>%filter(value=="Actual_Day_noD")%>%
  select("feature", "value", "coef","qval")%>% mutate(log2_fold_change = log2(exp(coef)))

colnames(AllShelfLife_Loc_genusv2_notsigdiff)[1] <- "Genus"
colnames(AllShelfLife_Loc_genusv2_notsigdiff)[4] <- "False Discovery Rate"

write.csv(AllShelfLife_Loc_genusv2_notsigdiff,"/local1/workdir1/tw488/CIDA/16S_Maaslin2_Output_11172024/AllShelfLife_Loc_genus_rarefitv2_output/Not_Significant_Genera_Log2foldchange", row.names = TRUE)

#genera diff abundant by area: dont need to make a table, include in text only
AllShelfLife_Loc_genus0.05v2_AZ<-AllShelfLife_Loc_genusv2%>%filter(value=="AZ")%>%
  filter(qval<"0.05")%>%
  select("feature", "value", "coef","qval")%>% mutate(log2_fold_change = log2(exp(coef)))

#change clostridia names so they will match when filtering
AllShelfLife_Loc_genus0.05v2_AZ$feature<-gsub("Clostridium.sensu.stricto.8", "Clostridium sensu stricto 8", AllShelfLife_Loc_genus0.05v2_AZ$feature)
AllShelfLife_Loc_genus0.05v2_AZ$feature<-gsub("Clostridium.sensu.stricto.1", "Clostridium sensu stricto 1", AllShelfLife_Loc_genus0.05v2_AZ$feature)
AllShelfLife_Loc_genus0.05v2_AZ$feature<-gsub("Clostridium.sensu.stricto.10", "Clostridium sensu stricto 10", AllShelfLife_Loc_genus0.05v2_AZ$feature)

#find rel ab and match it with the area data
genus_rel_abund_allShelfLifeandLoc_rarefied_4maaslin_nodcast_loc<-genus_rel_abund_allShelfLifeandLoc_rarefied_4maaslin_nodcast%>%filter(Rep!="NC")%>%group_by(taxon, Loc)%>%mutate(avgrelab=mean(rel_abund))%>%distinct(taxon,Loc,avgrelab)
genus_rel_abund_allShelfLifeandLoc_rarefied_4maaslin_nodcast_loc<-dcast(genus_rel_abund_allShelfLifeandLoc_rarefied_4maaslin_nodcast_loc,taxon~Loc,value.var="avgrelab")

#filter for genera in area data
shelflife_areagenera<-c(AllShelfLife_Loc_genus0.05v2_AZ$feature)
print(shelflife_areagenera)
#filter
genus_rel_abund_allShelfLifeandLoc_rarefied_4maaslin_nodcast_loc1<-genus_rel_abund_allShelfLifeandLoc_rarefied_4maaslin_nodcast_loc%>%filter(taxon %in% shelflife_areagenera)
colnames(genus_rel_abund_allShelfLifeandLoc_rarefied_4maaslin_nodcast_loc1)[1] <- "Genus"
#merge
colnames(AllShelfLife_Loc_genus0.05v2_AZ)[1] <- "Genus"

AllShelfLife_Loc_ASV_location_tbl <- genus_rel_abund_allShelfLifeandLoc_rarefied_4maaslin_nodcast_loc1%>% 
  full_join(AllShelfLife_Loc_genus0.05v2_AZ, by = "Genus")%>% select("Genus", "qval", "log2_fold_change","AZ","CA")

#reorganize and export 
AllShelfLife_Loc_ASV_location_tbl<-select(AllShelfLife_Loc_ASV_location_tbl,1,4,5,3,2)
colnames(AllShelfLife_Loc_ASV_location_tbl)[5] <- "Significance Level"

write.csv(AllShelfLife_Loc_ASV_location_tbl,"/local1/workdir1/tw488/CIDA/16S_Maaslin2_Output_11172024/AllShelfLife_Loc_genus_rarefitv2_output/ShelfLife_Loc_ASV_location_tbl", row.names = TRUE)

#__CALCULATE TOP 5 GENERA FOR A TIMEPOINT BY AVG RELATIVE ABUNDANCE___________________________________________________________
#most common genera (TOP 5) in D7 samples
D7toprelabund<-genus_rel_abund_allShelfLifeandLoc_rarefied_4maaslin_nodcast%>%filter(Day=="D7")%>%filter(Rep!="NC")%>%group_by(taxon)%>%mutate(avgrelab=mean(rel_abund))%>%distinct(taxon,avgrelab)
head(arrange(D7toprelabund,desc(avgrelab)), n = 5)

#most common genera (TOP 5) in D21/22 samples
topD22_DI<-dcast(genus_rel_abund_allShelfLifeandLoc_rarefied_4maaslin_nodcast2,taxon~Day,value.var="avgbtwnreps")
topD22<-topD22_DI%>%subset(select=-c(DI))

#total genera in shelf life samples
totalshelflifegenera<-genus_rel_abund_allShelfLifeandLoc_rarefied_4maaslin_nodcast%>%filter(rel_abund!="0")%>%distinct(taxon)
totalshelflifegenera
  
#most common genera (TOP 5) in DI samples
DItop5<-genus_rel_abundHDI.rarefied_4maaslin_no0_taxaMEANbyDAY%>%filter(Day=="DI")
head(dplyr::arrange(DItop5,desc(AvgbyDay)), n = 5)

#how many genera in DI samples total?
totalDIgenera_before5

#MAKE A LIST OF ALL GENERA IN DI SAMPLES >5% RELATIVE ABUNDANCE
genus_rel_abundDI_AVG

#most common genera (TOP 5) in H samples
Htop5<-genus_rel_abundHDI.rarefied_4maaslin_no0_taxaMEANbyDAY%>%filter(Day=="H")
head(dplyr::arrange(Htop5,desc(AvgbyDay)), n = 5)

#how many genera in H samples total?
genus_rel_abundH_total

#MAKE A LIST OF ALL GENERA IN H SAMPLES >5% RELATIVE ABUNDANCE
totalHgenera_over5


#______ASV level differential abundance in Shelf life Samples____________________________________________________________________________
allShelfLifeandLoc_rarefied_4maaslin_ASV<-data.frame(DI_D22_D28_noSEUSA@otu_table)

AllShelfLife_Loc_genus_rarefitv2 = Maaslin2(input_data = allShelfLifeandLoc_rarefied_4maaslin_ASV, 
                                            input_metadata = allShelfLifeandLoc_rarefied_4maaslin_sampledata,
                                            analysis_method = "CPLM",
                                            transform= "NONE",
                                            min_prevalence = 0.1,
                                            min_abundance = 0.01,
                                            normalization  = "NONE",
                                            max_pngs = 200,
                                            save_models="TRUE",
                                            save_scatter="TRUE",
                                            cores=2,
                                            random_effects = c("sampling_code"),
                                            output         = "/local1/workdir1/tw488/CIDA/16S_Maaslin2_Output_11172024/AllShelfLife_Loc_ASV_output", 
                                            fixed_effects  = c("Actual_Day_noD","Loc"),
                                            reference      = c("Loc,CA")) 

AllShelfLife_Loc_ASV<-read.table("/local1/workdir1/tw488/CIDA/16S_Maaslin2_Output_11172024/AllShelfLife_Loc_ASV_output/all_results.tsv",header=TRUE)
AllShelfLife_Loc_ASV0.05<-AllShelfLife_Loc_ASV%>%filter(qval<0.05)%>%filter(value=="Actual_Day_noD")%>%
  select("feature", "value", "coef","qval")%>% mutate(log2_fold_change = log2(exp(coef)))

allShelfLifeandLoc_rarefied_4maaslin_TAX1<-allShelfLifeandLoc_rarefied_4maaslin_TAX%>%rownames_to_column("feature") 

AllShelfLife_Loc_ASV0.05_1 <- allShelfLifeandLoc_rarefied_4maaslin_TAX1 %>% 
  full_join(AllShelfLife_Loc_ASV0.05, by = "feature")%>%na.omit%>% select("feature", "value", "coef","qval", "log2_fold_change","Family","Genus")

#TOTAL ENRICHED ASVs
#total + enriched ASVs 

library(dplyr)
AllShelfLife_Loc_ASV0.05_1_plus<-AllShelfLife_Loc_ASV0.05_1%>%filter(coef>0)%>% group_by(Genus)%>% dplyr::count(Genus)
colnames(AllShelfLife_Loc_ASV0.05_1_plus)[2] <- "Increasing"

#total - enriched ASVs 
AllShelfLife_Loc_ASV0.05_1_neg<-AllShelfLife_Loc_ASV0.05_1%>%filter(coef<0)%>% group_by(Genus) %>% dplyr::count(Genus)
colnames(AllShelfLife_Loc_ASV0.05_1_neg)[2] <- "Decreasing"

#TOTAL ASVS PER GENUS
allShelfLifeandLoc_rarefied_4maaslin_TAX_taxasum<-allShelfLifeandLoc_rarefied_4maaslin_TAX%>% group_by(Genus) %>% dplyr::count(Genus)
colnames(allShelfLifeandLoc_rarefied_4maaslin_TAX_taxasum)[2] <- "Total Number of ASVs assigned to Genus"

AllShelfLife_Loc_ASV0.05_1_alltogether <- AllShelfLife_Loc_ASV0.05_1_neg%>% 
  full_join(AllShelfLife_Loc_ASV0.05_1_plus, by = "Genus")%>%replace(is.na(.), 0)

AllShelfLife_Loc_ASV0.05_1_alltogether1 <- AllShelfLife_Loc_ASV0.05_1_alltogether %>% 
  full_join(allShelfLifeandLoc_rarefied_4maaslin_TAX_taxasum, by = "Genus")%>%na.omit()
colnames(AllShelfLife_Loc_ASV0.05_1_alltogether1)[1] <- "Genus of Enriched ASVs"
AllShelfLife_Loc_ASV0.05_1_alltogether1 <- AllShelfLife_Loc_ASV0.05_1_alltogether1[, c(1, 4, 2, 3)]

#total ASVs by Area
AllShelfLife_Loc_ASV0.05_Loc<-AllShelfLife_Loc_ASV%>%filter(qval<0.05)%>%filter(value=="AZ")%>%
  select("feature", "value", "coef","qval")%>% mutate(log2_fold_change = log2(exp(coef)))
#find taxa
AllShelfLife_Loc_ASV0.05_Loc_1 <- allShelfLifeandLoc_rarefied_4maaslin_TAX1 %>% 
  full_join(AllShelfLife_Loc_ASV0.05_Loc, by = "feature")%>%na.omit%>% select("feature", "value", "coef","qval", "log2_fold_change","Family","Genus")

#total + enriched ASVs 
AllShelfLife_Loc_ASV0.05_Loc_plus<-AllShelfLife_Loc_ASV0.05_Loc_1%>%filter(coef>0)%>% group_by(Genus) %>% dplyr::count(Genus)
colnames(AllShelfLife_Loc_ASV0.05_Loc_plus)[2] <- "Yuma"
#total - enriched ASVs 
AllShelfLife_Loc_ASV0.05_Loc_neg<-AllShelfLife_Loc_ASV0.05_Loc_1%>%filter(coef<0)%>% group_by(Genus) %>% dplyr::count(Genus)
colnames(AllShelfLife_Loc_ASV0.05_Loc_neg)[2] <- "Salinas"

#combine
AllShelfLife_Loc_ASV0.05_Loc_1_alltogether <- AllShelfLife_Loc_ASV0.05_Loc_neg%>% 
  full_join(AllShelfLife_Loc_ASV0.05_Loc_plus, by = "Genus")%>%replace(is.na(.), 0)
colnames(AllShelfLife_Loc_ASV0.05_Loc_1_alltogether)[1] <- "Genus of Enriched ASVs"

AllShelfLife_ASV_loc_tfsos <- AllShelfLife_Loc_ASV0.05_1_alltogether1  %>% 
  full_join(AllShelfLife_Loc_ASV0.05_Loc_1_alltogether, by = "Genus of Enriched ASVs")#%>%na.omit()


#TOTAL ASVs PER GENUS for Lelliottia
Lelliottia<-allShelfLifeandLoc_rarefied_4maaslin_TAX_taxasum%>% filter(Genus=="Lelliottia")
colnames(Lelliottia)[1] <- "Genus of Enriched ASVs"

AllShelfLife_ASV_loc_tfsos <- AllShelfLife_ASV_loc_tfsos %>% 
  full_join(Lelliottia, by = "Genus of Enriched ASVs")

write.csv(AllShelfLife_ASV_loc_tfsos ,"/local1/workdir1/tw488/CIDA/16S_Maaslin2_Output_11172024/AllShelfLife_Loc_ASV_output/ShelfLifeASVs.csv", row.names = TRUE)

#Find the # of reads per ASV
allShelfLifeandLoc_rarefied_4maaslin_ASV_rowsums<-as.data.frame(rowSums(allShelfLifeandLoc_rarefied_4maaslin_ASV))

#________PERMANOVA FOR select FIVE DAY SAMPLES, NMDS FOR FIVE DAY SAMPLES______________________________________________-

#WITH BRAY CURTIS
SHELFLIFEotu<-pstoveg_otu(fiveday)
#vegansample1<-pstoveg_sample(ps.noncontam_noNC.rarefied)
SHELFLIFEsd<-pstoveg_sd(fiveday)
SHELFLIFEsd$Monthlot = paste(SHELFLIFEsd$Month, SHELFLIFEsd$Lot, sep="_")
SHELFLIFEvegdist<-vegdist(SHELFLIFEotu, method="bray")

#NMDS
SHELFLIFE_NMDS <- metaMDS(SHELFLIFEvegdist)
stressplot(SHELFLIFE_NMDS)

SHELFLIFEsdNMDS<-SHELFLIFEsd%>%rownames_to_column("sample_name")
#SHELFLIFEsdNMDS<-SHELFLIFEsd
plot_SHELFLIFE_NMDS <- scores(SHELFLIFE_NMDS, display = "sites") %>% 
  as.data.frame() %>% 
  rownames_to_column("sample_name") %>% 
  full_join(SHELFLIFEsdNMDS, by = "sample_name")

fivedaypalNMDS<- c("chocolate", "#B2DF8A","red","#FF7F00", "#CAB2D6","black","darkblue","grey39","blue")

plot_SHELFLIFE_NMDS$Day<-factor(plot_SHELFLIFE_NMDS$Day,
                           levels=c("DI","D7","D12","D17","D22"))
str(plot_SHELFLIFE_NMDS$Day)

plot_SHELFLIFE_NMDS<-plot_SHELFLIFE_NMDS%>%filter(SampleName!="ITS_1222_A_D7_3_AZ")


plot_SHELFLIFE_NMDS_nmds <- ggplot(plot_SHELFLIFE_NMDS, aes(x = NMDS1, y = NMDS2, color = Day, shape=Loc)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_manual(values = fivedaypalNMDS) +
  stat_ellipse(aes(color=Loc), linetype=2, level=0.5)+
  theme_classic()
plot_SHELFLIFE_NMDS_nmds
ggsave("fiveday_Bact_NMDS_04162025.tiff", units="in", width=8, height=6, dpi=300, compression = 'lzw')

#PERMANOVA code

SHELFLIFEsd_otu<-cbind(SHELFLIFEsd,SHELFLIFEotu)
SHELFLIFEsd_otu<-SHELFLIFEsd_otu%>%mutate(MonthlotDay=paste(Monthlot,Day,sep="_"),.after=Monthlot)

SHELFLIFEsd_otu%>%group_by(MonthlotDay)%>%dplyr::count()%>%View()
SHELFLIFEsd_otu<-SHELFLIFEsd_otu%>%filter(SampleName!="16S_1222_A_D7_3_AZ")

#how function from permute package
hsubSHELFLIFE <- with(SHELFLIFEsd_otu, how(within=Within(type="free"),blocks=Monthlot,plots=Plots(strata=MonthlotDay, type="free")))
hsubSHELFLIFE

SHELFLIFEsd_permanova<-SHELFLIFEsd_otu%>%select(SampleType,quant_reading,is.neg,SampleName,Month,Lot,Day,Rep,Loc,Actual_Day,Monthlot,MonthlotDay)
table(SHELFLIFEsd_permanova$Monthlot)
table(table(SHELFLIFEsd_permanova$MonthlotDay))
table(SHELFLIFEsd_permanova$Actual_Day)

SHELFLIFEotu_permanova<-SHELFLIFEsd_otu%>%ungroup()%>%select(-SampleType,-quant_reading,-is.neg,-SampleName,-Month,-Lot,-Day,-Rep,-Loc,-Actual_Day,-Monthlot,-MonthlotDay)

#ensure replicates are permuted together
testsd<-SHELFLIFEsd_permanova
testsd$permute1<-shuffle(178,control=hsubSHELFLIFE)
testsd$origorder<-1:178

#visualize permutations
testsd%>%pivot_longer(permute1:origorder, names_to="column", values_to="position", cols_vary = "slowest") %>%
  ggplot(aes(x=column, y=position, group=SampleName, color=Day))+geom_line()

#make in to distance object
SHELFLIFEvegdist4<-vegdist(SHELFLIFEotu_permanova, method="bray")

#plug in here
SHELFLIFE_perm_trial <- adonis2(SHELFLIFEvegdist4~Actual_Day, data=SHELFLIFEsd_permanova, permutations = hsubSHELFLIFE, method = "bray", 
                          sqrt.dist = FALSE, add = FALSE, by = "terms",
                          parallel = getOption("mc.cores"), na.action = na.fail)
SHELFLIFE_perm_trial

fiveday_16s_dayloc <- adonis2(SHELFLIFEvegdist4~Actual_Day*Loc, data=SHELFLIFEsd_permanova, permutations = hsubSHELFLIFE, method = "bray", 
                                sqrt.dist = FALSE, add = FALSE, by = "terms",
                                parallel = getOption("mc.cores"), na.action = na.fail)
fiveday_16s_dayloc


# pairwise adonis
#cant do pairwise adonis w/o pesudoreplication of the two reps, so doing it one by one
#DID7--significant
SHELFLIFE_perm_DID7 <- adonis2(vegdist(SHELFLIFEotu_permanova[SHELFLIFEsd_permanova$Actual_Day%in%c("DI","D7"),], method="bray")~Actual_Day, 
                                data=SHELFLIFEsd_permanova[SHELFLIFEsd_permanova$Actual_Day%in%c("DI","D7"),], 
                                permutations = with(SHELFLIFEsd_otu[SHELFLIFEsd_permanova$Actual_Day%in%c("DI","D7"),], how(within=Within(type="free"),blocks=Monthlot,plots=Plots(strata=MonthlotDay, type="free"))), 
                               method = "bray", 
                                sqrt.dist = FALSE, add = FALSE, by = "terms",
                                parallel = getOption("mc.cores"), na.action = na.fail)
SHELFLIFE_perm_DID7

#DID12--significant
SHELFLIFE_perm_DID12 <- adonis2(vegdist(SHELFLIFEotu_permanova[SHELFLIFEsd_permanova$Actual_Day%in%c("DI","D12"),], method="bray")~Actual_Day, 
                               data=SHELFLIFEsd_permanova[SHELFLIFEsd_permanova$Actual_Day%in%c("DI","D12"),], 
                               permutations = with(SHELFLIFEsd_otu[SHELFLIFEsd_permanova$Actual_Day%in%c("DI","D12"),], how(within=Within(type="free"),blocks=Monthlot,plots=Plots(strata=MonthlotDay, type="free"))), 
                               method = "bray", 
                               sqrt.dist = FALSE, add = FALSE, by = "terms",
                               parallel = getOption("mc.cores"), na.action = na.fail)
SHELFLIFE_perm_DID12

#DID17--significant
SHELFLIFE_perm_DID17 <- adonis2(vegdist(SHELFLIFEotu_permanova[SHELFLIFEsd_permanova$Actual_Day%in%c("DI","D17"),], method="bray")~Actual_Day, 
                               data=SHELFLIFEsd_permanova[SHELFLIFEsd_permanova$Actual_Day%in%c("DI","D17"),], 
                               permutations = with(SHELFLIFEsd_otu[SHELFLIFEsd_permanova$Actual_Day%in%c("DI","D17"),], how(within=Within(type="free"),blocks=Monthlot,plots=Plots(strata=MonthlotDay, type="free"))), 
                               method = "bray", 
                               sqrt.dist = FALSE, add = FALSE, by = "terms",
                               parallel = getOption("mc.cores"), na.action = na.fail)
SHELFLIFE_perm_DID17

#DID22-significant
SHELFLIFE_perm_DID22 <- adonis2(vegdist(SHELFLIFEotu_permanova[SHELFLIFEsd_permanova$Actual_Day%in%c("DI","D22"),], method="bray")~Actual_Day, 
                              data=SHELFLIFEsd_permanova[SHELFLIFEsd_permanova$Actual_Day%in%c("DI","D22"),], 
                               permutations = with(SHELFLIFEsd_otu[SHELFLIFEsd_permanova$Actual_Day%in%c("DI","D22"),], how(within=Within(type="free"),blocks=Monthlot,plots=Plots(strata=MonthlotDay, type="free"))), 
                               method = "bray", 
                               sqrt.dist = FALSE, add = FALSE, by = "terms",
                               parallel = getOption("mc.cores"), na.action = na.fail)
SHELFLIFE_perm_DID22

#D7D12-significant

SHELFLIFE_perm_D7D12 <- adonis2(vegdist(SHELFLIFEotu_permanova[SHELFLIFEsd_permanova$Actual_Day%in%c("D7","D12"),], method="bray")~Actual_Day, 
                               data=SHELFLIFEsd_permanova[SHELFLIFEsd_permanova$Actual_Day%in%c("D7","D12"),], 
                               permutations = with(SHELFLIFEsd_otu[SHELFLIFEsd_permanova$Actual_Day%in%c("D7","D12"),], how(within=Within(type="free"),blocks=Monthlot,plots=Plots(strata=MonthlotDay, type="free"))), 
                               method = "bray", 
                               sqrt.dist = FALSE, add = FALSE, by = "terms",
                               parallel = getOption("mc.cores"), na.action = na.fail)
SHELFLIFE_perm_D7D12

#D7D17-significant
SHELFLIFE_perm_D7D17 <- adonis2(vegdist(SHELFLIFEotu_permanova[SHELFLIFEsd_permanova$Actual_Day%in%c("D7","D17"),], method="bray")~Actual_Day, 
                                data=SHELFLIFEsd_permanova[SHELFLIFEsd_permanova$Actual_Day%in%c("D7","D17"),], 
                                permutations = with(SHELFLIFEsd_otu[SHELFLIFEsd_permanova$Actual_Day%in%c("D7","D17"),], how(within=Within(type="free"),blocks=Monthlot,plots=Plots(strata=MonthlotDay, type="free"))), 
                                method = "bray", 
                                sqrt.dist = FALSE, add = FALSE, by = "terms",
                                parallel = getOption("mc.cores"), na.action = na.fail)
SHELFLIFE_perm_D7D17

#D7D22- significant
SHELFLIFE_perm_D7D22 <- adonis2(vegdist(SHELFLIFEotu_permanova[SHELFLIFEsd_permanova$Actual_Day%in%c("D7","D22"),], method="bray")~Actual_Day, 
                                data=SHELFLIFEsd_permanova[SHELFLIFEsd_permanova$Actual_Day%in%c("D7","D22"),], 
                                permutations = with(SHELFLIFEsd_otu[SHELFLIFEsd_permanova$Actual_Day%in%c("D7","D22"),], how(within=Within(type="free"),blocks=Monthlot,plots=Plots(strata=MonthlotDay, type="free"))), 
                                method = "bray", 
                                sqrt.dist = FALSE, add = FALSE, by = "terms",
                                parallel = getOption("mc.cores"), na.action = na.fail)
SHELFLIFE_perm_D7D22

#D12D17- significant

SHELFLIFE_perm_D12D17 <- adonis2(vegdist(SHELFLIFEotu_permanova[SHELFLIFEsd_permanova$Actual_Day%in%c("D12","D17"),], method="bray")~Actual_Day, 
                                 data=SHELFLIFEsd_permanova[SHELFLIFEsd_permanova$Actual_Day%in%c("D12","D17"),], 
                                 permutations = with(SHELFLIFEsd_otu[SHELFLIFEsd_permanova$Actual_Day%in%c("D12","D17"),], how(within=Within(type="free"),blocks=Monthlot,plots=Plots(strata=MonthlotDay, type="free"))), 
                                 method = "bray", 
                                 sqrt.dist = FALSE, add = FALSE, by = "terms",
                                 parallel = getOption("mc.cores"), na.action = na.fail)
SHELFLIFE_perm_D12D17

#D12D22- significant
SHELFLIFE_perm_D12D22 <- adonis2(vegdist(SHELFLIFEotu_permanova[SHELFLIFEsd_permanova$Actual_Day%in%c("D12","D22"),], method="bray")~Actual_Day, 
                                data=SHELFLIFEsd_permanova[SHELFLIFEsd_permanova$Actual_Day%in%c("D12","D22"),], 
                                permutations = with(SHELFLIFEsd_otu[SHELFLIFEsd_permanova$Actual_Day%in%c("D12","D22"),], how(within=Within(type="free"),blocks=Monthlot,plots=Plots(strata=MonthlotDay, type="free"))), 
                                method = "bray", 
                                sqrt.dist = FALSE, add = FALSE, by = "terms",
                                parallel = getOption("mc.cores"), na.action = na.fail)
SHELFLIFE_perm_D12D22

#d17D22- not sig

SHELFLIFE_perm_D17D22 <- adonis2(vegdist(SHELFLIFEotu_permanova[SHELFLIFEsd_permanova$Actual_Day%in%c("D17","D22"),], method="bray")~Actual_Day, 
                                 data=SHELFLIFEsd_permanova[SHELFLIFEsd_permanova$Actual_Day%in%c("D17","D22"),], 
                                 permutations = with(SHELFLIFEsd_otu[SHELFLIFEsd_permanova$Actual_Day%in%c("D17","D22"),], how(within=Within(type="free"),blocks=Monthlot,plots=Plots(strata=MonthlotDay, type="free"))), 
                                 method = "bray", 
                                 sqrt.dist = FALSE, add = FALSE, by = "terms",
                                 parallel = getOption("mc.cores"), na.action = na.fail)
SHELFLIFE_perm_D17D22
#______SEVEN DAY INTERVAL PERMANOVA AND NMDS_____________________________________________________________________

#WITH BRAY CURTIS
sevendayotu<-pstoveg_otu(sevenday)
sevendaysd<-pstoveg_sd(sevenday)
sevendaysd$Monthlot = paste(sevendaysd$Month, sevendaysd$Lot, sep="_")
sevendayvegdist<-vegdist(sevendayotu, method="bray") 

#NMDS
sevenday_NMDS <- metaMDS(sevendayvegdist)
stressplot(sevenday_NMDS)

sevendaysdNMDS<-sevendaysd%>%rownames_to_column("sample")
plot_sevenday_NMDS <- scores(sevenday_NMDS, display = "sites") %>% 
  as.data.frame() %>% 
  rownames_to_column("sample") %>% 
  full_join(sevendaysdNMDS, by = "sample")%>% 
  filter(Actual_Day!="H")

plot_sevenday_NMDS_nmds <- ggplot(plot_sevenday_NMDS, aes(x = NMDS1, y = NMDS2, color = Actual_Day, shape = Loc)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(title = "Seven Day Shelf Life Samples NMDS by Month and Location _Bacteria")+
  theme_classic()
plot_sevenday_NMDS_nmds

#PERMANOVA
sevendaysd_otu<-cbind(sevendaysd,sevendayotu)
sevendaysd_otu<-sevendaysd_otu%>%mutate(MonthlotDay=paste(Monthlot,Day,sep="_"),.after=Monthlot)

#check to make sure there ar eveent # samples
table(SHELFLIFEsd_otu$Monthlot)
# Must subsample plots to be balanced, but do this 100 times because of random sampling
hsubsevenday <- with(sevendaysd_otu, how(within=Within(type="free"),blocks=Monthlot,plots=Plots(strata=MonthlotDay, type="free")))
hsubsevenday

sevendaysd_permanova<-sevendaysd_otu%>%select(SampleType,quant_reading,is.neg,SampleName,Month,Lot,Day,Rep,Loc,Actual_Day,Monthlot,MonthlotDay)
table(sevendaysd_permanova$Monthlot)
table(sevendaysd_permanova$Actual_Day)

sevendayotu_permanova<-sevendaysd_otu%>%ungroup()%>%select(-SampleType,-quant_reading,-is.neg,-SampleName,-Month,-Lot,-Day,-Rep,-Loc,-Actual_Day,-Monthlot,-MonthlotDay)

#ensure replicates are permuted together
sevendaytestsd<-sevendaysd_permanova
sevendaytestsd$permute1<-shuffle(48,control=hsubsevenday)
sevendaytestsd$origorder<-1:48

#visualize permutations
sevendaytestsd%>%pivot_longer(permute1:origorder, names_to="column", values_to="position", cols_vary = "slowest") %>%
  ggplot(aes(x=column, y=position, group=SampleName, color=Day))+geom_line()

#make in to distance object
sevendayvegdist4<-vegdist(sevendayotu_permanova, method="bray")

#plug in here
sevenday_perm_trial <- adonis2(sevendayvegdist4~Actual_Day, data=sevendaysd_permanova, permutations = hsubsevenday, method = "bray", 
                                sqrt.dist = FALSE, add = FALSE, by = "terms",
                                parallel = getOption("mc.cores"), na.action = na.fail)
sevenday_perm_trial

#with day and loc
sevenday_perm_trial_dayloc <- adonis2(sevendayvegdist4~Actual_Day*Loc, data=sevendaysd_permanova, permutations = hsubsevenday, method = "bray", 
                               sqrt.dist = FALSE, add = FALSE, by = "terms",
                               parallel = getOption("mc.cores"), na.action = na.fail)
sevenday_perm_trial_dayloc

# pairwise adonis OF 7 DAY INTERVAL SAMPLES
#cant do pairwise adonis w/o pesudoreplication of the two reps, so doing it one by one

#DID7--significant
sevenday_perm_DID7 <- adonis2(vegdist(sevendayotu_permanova[sevendaysd_permanova$Actual_Day%in%c("DI","D7"),], method="bray")~Actual_Day, 
                               data=sevendaysd_permanova[sevendaysd_permanova$Actual_Day%in%c("DI","D7"),], 
                               permutations = with(sevendaysd_otu[sevendaysd_permanova$Actual_Day%in%c("DI","D7"),], how(within=Within(type="free"),blocks=Monthlot,plots=Plots(strata=MonthlotDay, type="free"))), 
                               method = "bray", 
                               sqrt.dist = FALSE, add = FALSE, by = "terms",
                               parallel = getOption("mc.cores"), na.action = na.fail)
sevenday_perm_DID7

#DID14--not significant
sevenday_perm_DID14 <- adonis2(vegdist(sevendayotu_permanova[sevendaysd_permanova$Actual_Day%in%c("DI","D14"),], method="bray")~Actual_Day, 
                                data=sevendaysd_permanova[sevendaysd_permanova$Actual_Day%in%c("DI","D14"),], 
                                permutations = with(sevendaysd_otu[sevendaysd_permanova$Actual_Day%in%c("DI","D14"),], how(within=Within(type="free"),blocks=Monthlot,plots=Plots(strata=MonthlotDay, type="free"))), 
                                method = "bray", 
                                sqrt.dist = FALSE, add = FALSE, by = "terms",
                                parallel = getOption("mc.cores"), na.action = na.fail)
sevenday_perm_DID14

#DID21--significant
sevenday_perm_DID21 <- adonis2(vegdist(sevendayotu_permanova[sevendaysd_permanova$Actual_Day%in%c("DI","D21"),], method="bray")~Actual_Day, 
                                data=sevendaysd_permanova[sevendaysd_permanova$Actual_Day%in%c("DI","D21"),], 
                                permutations = with(sevendaysd_otu[sevendaysd_permanova$Actual_Day%in%c("DI","D21"),], how(within=Within(type="free"),blocks=Monthlot,plots=Plots(strata=MonthlotDay, type="free"))), 
                                method = "bray", 
                                sqrt.dist = FALSE, add = FALSE, by = "terms",
                                parallel = getOption("mc.cores"), na.action = na.fail)
sevenday_perm_DID21

#DID22-not significant
sevenday_perm_DID28 <- adonis2(vegdist(sevendayotu_permanova[sevendaysd_permanova$Actual_Day%in%c("DI","D28"),], method="bray")~Actual_Day, 
                                data=sevendaysd_permanova[sevendaysd_permanova$Actual_Day%in%c("DI","D28"),], 
                                permutations = with(sevendaysd_otu[sevendaysd_permanova$Actual_Day%in%c("DI","D28"),], how(within=Within(type="free"),blocks=Monthlot,plots=Plots(strata=MonthlotDay, type="free"))), 
                                method = "bray", 
                                sqrt.dist = FALSE, add = FALSE, by = "terms",
                                parallel = getOption("mc.cores"), na.action = na.fail)
sevenday_perm_DID28

#D7D14-significant

sevenday_perm_D7D14 <- adonis2(vegdist(sevendayotu_permanova[sevendaysd_permanova$Actual_Day%in%c("D7","D14"),], method="bray")~Actual_Day, 
                                data=sevendaysd_permanova[sevendaysd_permanova$Actual_Day%in%c("D7","D14"),], 
                                permutations = with(sevendaysd_otu[sevendaysd_permanova$Actual_Day%in%c("D7","D14"),], how(within=Within(type="free"),blocks=Monthlot,plots=Plots(strata=MonthlotDay, type="free"))), 
                                method = "bray", 
                                sqrt.dist = FALSE, add = FALSE, by = "terms",
                                parallel = getOption("mc.cores"), na.action = na.fail)
sevenday_perm_D7D14

#D7D21-not significant
sevenday_perm_D7D21 <- adonis2(vegdist(sevendayotu_permanova[sevendaysd_permanova$Actual_Day%in%c("D7","D21"),], method="bray")~Actual_Day, 
                                data=sevendaysd_permanova[sevendaysd_permanova$Actual_Day%in%c("D7","D21"),], 
                                permutations = with(sevendaysd_otu[sevendaysd_permanova$Actual_Day%in%c("D7","D21"),], how(within=Within(type="free"),blocks=Monthlot,plots=Plots(strata=MonthlotDay, type="free"))), 
                                method = "bray", 
                                sqrt.dist = FALSE, add = FALSE, by = "terms",
                                parallel = getOption("mc.cores"), na.action = na.fail)
sevenday_perm_D7D21

#D7D28- not significant
sevenday_perm_D7D28 <- adonis2(vegdist(sevendayotu_permanova[sevendaysd_permanova$Actual_Day%in%c("D7","D28"),], method="bray")~Actual_Day, 
                                data=sevendaysd_permanova[sevendaysd_permanova$Actual_Day%in%c("D7","D28"),], 
                                permutations = with(sevendaysd_otu[sevendaysd_permanova$Actual_Day%in%c("D7","D28"),], how(within=Within(type="free"),blocks=Monthlot,plots=Plots(strata=MonthlotDay, type="free"))), 
                                method = "bray", 
                                sqrt.dist = FALSE, add = FALSE, by = "terms",
                                parallel = getOption("mc.cores"), na.action = na.fail)
sevenday_perm_D7D28

#D14D21- not significant

sevenday_perm_D14D21 <- adonis2(vegdist(sevendayotu_permanova[sevendaysd_permanova$Actual_Day%in%c("D14","D21"),], method="bray")~Actual_Day, 
                                 data=sevendaysd_permanova[sevendaysd_permanova$Actual_Day%in%c("D14","D21"),], 
                                 permutations = with(sevendaysd_otu[sevendaysd_permanova$Actual_Day%in%c("D14","D21"),], how(within=Within(type="free"),blocks=Monthlot,plots=Plots(strata=MonthlotDay, type="free"))), 
                                 method = "bray", 
                                 sqrt.dist = FALSE, add = FALSE, by = "terms",
                                 parallel = getOption("mc.cores"), na.action = na.fail)
sevenday_perm_D14D21

#D14D28- not significant
sevenday_perm_D14D28 <- adonis2(vegdist(sevendayotu_permanova[sevendaysd_permanova$Actual_Day%in%c("D14","D28"),], method="bray")~Actual_Day, 
                                 data=sevendaysd_permanova[sevendaysd_permanova$Actual_Day%in%c("D14","D28"),], 
                                 permutations = with(sevendaysd_otu[sevendaysd_permanova$Actual_Day%in%c("D14","D28"),], how(within=Within(type="free"),blocks=Monthlot,plots=Plots(strata=MonthlotDay, type="free"))), 
                                 method = "bray", 
                                 sqrt.dist = FALSE, add = FALSE, by = "terms",
                                 parallel = getOption("mc.cores"), na.action = na.fail)
sevenday_perm_D14D28

#d21D28- not sig

sevenday_perm_D21D28 <- adonis2(vegdist(sevendayotu_permanova[sevendaysd_permanova$Actual_Day%in%c("D21","D28"),], method="bray")~Actual_Day, 
                                 data=sevendaysd_permanova[sevendaysd_permanova$Actual_Day%in%c("D21","D28"),], 
                                 permutations = with(sevendaysd_otu[sevendaysd_permanova$Actual_Day%in%c("D21","D28"),], how(within=Within(type="free"),blocks=Monthlot,plots=Plots(strata=MonthlotDay, type="free"))), 
                                 method = "bray", 
                                 sqrt.dist = FALSE, add = FALSE, by = "terms",
                                 parallel = getOption("mc.cores"), na.action = na.fail)
sevenday_perm_D21D28

#Alpha diversity OF 7 DAY INTERVAL SAMPLES

sevendayalphadiv<-estimate_richness(sevenday, split = TRUE, measures = NULL)
sevendayalphadiv$sample <- row.names(sevendayalphadiv)  
sevendayalphadiv <- sevendayalphadiv %>% separate(sample, c("x16s", "month", "lot", "day", "rep", "loc","r1"), "_")
sevendayalphadiv$sample <- row.names(sevendayalphadiv) 
sevendayalphadiv$sampling_code = paste(sevendayalphadiv$month, sevendayalphadiv$lot, sep="_")
sevendayalphadiv <- sevendayalphadiv %>% 
  as.data.frame() %>% 
  full_join(TFSOSv2, by = "sampling_code")%>%na.omit

#calculate Pielou's Evenness
sevendayalphadiv <-sevendayalphadiv %>% mutate (Pielou=(sevendayalphadiv$Shannon)/log(sevendayalphadiv$Observed))

#Visualize
#Observed
sevendayalphadiv  %>% ggplot(aes(x=day, y=Observed)) +
  geom_jitter(size=0.25) +
  labs(x="Day",
       y="Number of Observed ASVs") +
  guides(color = guide_legend(override.aes = list(size=1))) +
  theme_classic()

#spearman test to correlate time and alpha div

#Pielou
sevendayalphadiv%>% ggplot(aes(x=day, y=Pielou)) +
  geom_jitter(size=0.25) +
  labs(x="Day",
       y="Pielou's Evenness") +
  guides(color = guide_legend(override.aes = list(size=1))) +
  theme_classic()

#SHANNON
sevendayalphadiv%>% ggplot(aes(x=day, y=Shannon)) +
  geom_jitter(size=0.25) +
  labs(x="Day",
       y="Shannon Index") +
  guides(color = guide_legend(override.aes = list(size=1))) +
  theme_classic()


#ANOVA for seven day interval alpha diversity values as categorical variable
sevendayalphadiv_noH<-sevendayalphadiv%>% filter(day!="H") 

#order Days as factors
sevendayalphadiv_noH$day<-factor(sevendayalphadiv_noH$day,
                            levels=c("DI","D7","D12","D17","D22"))

#ANOVA to compare days
#Shannon
ShnANOVA<- aov(Shannon ~ day, data =sevendayalphadiv_noH)
em_ShnANOVA<-emmeans(ShnANOVA,~day)

pairs.em_ShnANOVA<-as.data.frame(pairs(em_ShnANOVA))
print(pairs.em_ShnANOVA)
cldList(p.value~contrast,data=pairs.em_ShnANOVA, remove.space="TRUE")

#Richness
ObsANOVA<- aov(Observed~ day, data =sevendayalphadiv_noH)
em_ObsANOVA<-emmeans(ObsANOVA,~day)

pairs.em_ObsANOVA<-as.data.frame(pairs(em_ObsANOVA))
print(pairs.em_ObsANOVA)
cldList(p.value~contrast,data=pairs.em_ObsANOVA, remove.space="TRUE")

#Evenness
PieANOVA<- aov(Pielou ~ day, data =sevendayalphadiv_noH)
em_PieANOVA<-emmeans(PieANOVA,~day)

pairs.em_PieANOVA<-as.data.frame(pairs(em_PieANOVA))
print(pairs.em_PieANOVA)
cldList(p.value~contrast,data=pairs.em_PieANOVA, remove.space="TRUE")

#calculate average values for sevendayalphadiv for the data table
sevendayalphadiv_noH <-sevendayalphadiv_noH %>% group_by(day)%>%mutate(avgPielou=mean(Pielou))%>%mutate(avgShn=mean(Shannon))%>%mutate(avgObs=mean(Observed))
sevendayalphadiv_avgs <-sevendayalphadiv_noH %>%distinct(day,avgPielou,avgShn,avgObs)
write.csv(sevendayalphadiv_avgs,"/local1/workdir1/tw488/CIDA/16S_Maaslin2_Output_11172024/Allshelflife_D7_alphadivTable_avgsforeachalphadivindex", row.names = TRUE)

#Five day interval samples alpha diversity__________________________________________

fivedayalphadiv<-estimate_richness(fiveday, split = TRUE, measures = NULL)
fivedayalphadiv$sample <- row.names(fivedayalphadiv)  
fivedayalphadiv <- fivedayalphadiv %>% 
  separate(sample, c("x16s", "month", "lot", "day", "rep", "loc","r1"), "_")%>%
  ungroup()%>%
  filter (SampleName!="X16S_1222_A_D7_3_AZ_R1")
fivedayalphadiv$sample <- row.names(fivedayalphadiv) 
fivedayalphadiv$sampling_code = paste(fivedayalphadiv$month, fivedayalphadiv$lot, sep="_")

#calculate Pielou'e Evenness
fivedayalphadiv <-fivedayalphadiv %>% mutate (Pielou=(fivedayalphadiv$Shannon)/log(fivedayalphadiv$Observed))

#ANOVA for d7 alpha div values as categorical variable

#correct days
fivedayalphadiv["day"][fivedayalphadiv["day"] == "D7"] <- "Day 7"
fivedayalphadiv["day"][fivedayalphadiv["day"] == "DI"] <- "Day Initial"
fivedayalphadiv["day"][fivedayalphadiv["day"] == "D12"] <- "Day 12"
fivedayalphadiv["day"][fivedayalphadiv["day"] == "D22"] <- "Day 22"
fivedayalphadiv["day"][fivedayalphadiv["day"] == "D17"] <- "Day 17"
fivedayalphadiv$day<-factor(fivedayalphadiv$day,
                            levels=c("Day Initial","Day 7","Day 12","Day 17","Day 22"))

#Pielou
fivedayalphadiv%>% ggplot(aes(x=day, y=Pielou)) +
  geom_boxplot(size=0.25) +
  labs(x="Day",
       y="Pielou's Evenness") +
  guides(color = guide_legend(override.aes = list(size=1))) +
  theme_classic()
#Richenss
fivedayalphadiv%>% ggplot(aes(x=day, y=Observed)) +
  geom_boxplot(size=0.25) +
  labs(x="Day",
       y="Richnness") +
  guides(color = guide_legend(override.aes = list(size=1))) +
  theme_classic()

fivedayalphadiv%>% ggplot(aes(x=day, y=Shannon)) +
  geom_boxplot(size=0.25) +
  labs(x="Day",
       y="Shannon") +
  guides(color = guide_legend(override.aes = list(size=1))) +
  theme_classic()

#ANOVA for 5 day interval alpha diversity values as categorical variable
fivedayalphadiv_noH<-fivedayalphadiv%>% filter(day!="H") 

#Shannon
fiveShnANOVA<- aov(Shannon ~ day, data =fivedayalphadiv_noH)
em_fiveShnANOVA<-emmeans(fiveShnANOVA,~day)

pairs.em_fiveShnANOVA<-as.data.frame(pairs(em_fiveShnANOVA))
print(pairs.em_fiveShnANOVA)
cldList(p.value~contrast,data=pairs.em_fiveShnANOVA, remove.space="TRUE")

#Richness
fiveObsANOVA<- aov(Observed~ day, data =fivedayalphadiv_noH)
em_fiveObsANOVA<-emmeans(fiveObsANOVA,~day)

pairs.em_fiveObsANOVA<-as.data.frame(pairs(em_fiveObsANOVA))
print(pairs.em_fiveObsANOVA)
cldList(p.value~contrast,data=pairs.em_fiveObsANOVA, remove.space="TRUE")

#Evenness
fivePieANOVA<- aov(Pielou ~ day, data =fivedayalphadiv_noH)
em_fivePieANOVA<-emmeans(fivePieANOVA,~day)

pairs.em_fivePieANOVA<-as.data.frame(pairs(em_fivePieANOVA))
print(pairs.em_fivePieANOVA)
cldList(p.value~contrast,data=pairs.em_fivePieANOVA, remove.space="TRUE")

#calculate average values for fivedayalphadiv
fivedayalphadiv_noH <-fivedayalphadiv_noH %>% group_by(day)%>%mutate(avgPielou=mean(Pielou))%>%mutate(avgShn=mean(Shannon))%>%mutate(avgObs=mean(Observed))
fivedayalphadiv_avgs <-fivedayalphadiv_noH %>%distinct(day,avgPielou,avgShn,avgObs)
write.csv(fivedayalphadiv_avgs,"/local1/workdir1/tw488/CIDA/16S_Maaslin2_Output_11172024/Allshelflife_5day_alphadivTable_avgsforeachalphadivindex", row.names = TRUE)


