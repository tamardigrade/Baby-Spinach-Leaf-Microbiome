#ITS analyses for CIDA sequencing paper

#Loading the environment
load('/local1/workdir1/tw488/RData_CIDA/REnvr_from_Github_Code.RData')
#VIEW THE SAVED FILE IN A THE DIRECTORY
dir()
#SAVE ENVIRONMENT
save.image(file='/local1/workdir1/tw488/RData_CIDA/REnvr_from_Github_Code.RData')

#PACKAGES NEEDED____________________________________________________________________
library(ggplot2)
library(tidyr)
library(dplyr)
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
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
#if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
#devtools::install_github("jbisanz/qiime2R") # current version is 0.99.20
library(devtools)
library("qiime2R")
BiocManager::install("decontam")
library(decontam); packageVersion("decontam")
library(pbkrtest)
library(lme4)
library(emmeans)
library(multcompView)
library(rcompanion)
library(forcats)
#Just for reference:for testing normality
library(ggpubr)

#CODE IMPORTED FROM ITS_Runs1_7.R_______________________________________
ITSps<-qza_to_phyloseq(
  features="/local1/workdir1/tw488/CIDA_Sequencing_Runs/ITS_RUNS1_7_allinone/Trimmomaticoutput/Pair_usethis/uchime-dn-o/table-nonchimeric-wo-borderline.qza",
  taxonomy="/local1/workdir1/tw488/CIDA_Sequencing_Runs/ITS_RUNS1_7_allinone/Trimmomaticoutput/Pair_usethis/UNITEref-taxonomy.qza",
  metadata = "/local1/workdir1/tw488/CIDA_Sequencing_Runs/ITS_RUNS1_7_allinone/Trimmomaticoutput/Pair_usethis/ITSmetadata_05292024.txt")

#check out the ps object
ITSsamdf<-as.data.frame(ITSps@sam_data)
ITSotutbl<-as.data.frame(ITSps@otu_table)
ITStaxtbl<-as.data.frame(ITSps@tax_table)
sample_variables(ITSps)

#prune mock samples from the ps object for decontam

ITSps1 <- subset_samples(ITSps, Month != "mock")

#check out new ps object
ITSsamdf1<-as.data.frame(ITSps1@sam_data)
ITSotutbl1<-as.data.frame(ITSps1@otu_table)
ITStaxtbl1<-as.data.frame(ITSps1@tax_table)

#run decontam wth combined method
contITS<-isContaminant(ITSps1, conc = "quant.reading", neg = "is.neg", method = "combined",
                       batch = NULL, threshold = 0.1,normalize = TRUE,detailed = TRUE)
contITS

#check out score dist on ggplot histogram
p<-ggplot(contITS, aes(x=p.freq)) + geom_histogram()
p

o<-ggplot(contITS, aes(x=p.prev)) + geom_histogram()
o

i<-ggplot(contITS, aes(x=p)) + geom_histogram()
i

contITS1<-contITS%>%
  filter(contaminant=="TRUE")
#check to see if contaminants make sense by looking at the taxonomic info
ITScontamTaxtbl <- merge(contITS1, ITStaxtbl, by = 'row.names', all = TRUE)

ITScontamTaxtbl <-ITScontamTaxtbl %>%filter(contaminant=="TRUE")
#write.csv(ITScontamTaxtbl,"/local1/workdir1/tw488/CIDA_Sequencing_Runs/ITS_RUNS1_6_allinoneCOPY/Trimmomaticoutput/Pair_usethis/ITS_contam_Tax_tbl.csv")

toBremovITS<-ITScontamTaxtbl$Row.names
toBremovITS

#Prune Contaminant Taxa from ps object
ITSps.noncontam <- prune_taxa(!(taxa_names(ITSps1) %in% toBremovITS), ITSps1)
ITSps.noncontam

#check out the ps object
pssamnonc<-as.data.frame(ITSps.noncontam@sam_data)
ITSotutblnonc<-as.data.frame(ITSps.noncontam@otu_table)
ITStaxtblnonc<-as.data.frame(ITSps.noncontam@tax_table)

#____________________________________________________________________
#Retain ps object with spinach reads for future reference
ITSps.noncontam_withspinach <-ITSps.noncontam

#take a look at what plants were amplified
#PlantsAmplified <- subset_taxa(ITSps.noncontam_withspinach, Kingdom== "Viridiplantae")
#PlantsAmplifiedtaxtblnonc<-as.data.frame(PlantsAmplified@tax_table)

#code to remove ASVs<2 and chloroplast/mitochondria
ITSps.noncontam <- subset_taxa(ITSps.noncontam, Family!= "mitochondria")
ITSps.noncontam <- subset_taxa(ITSps.noncontam, Class!="Chloroplast" )
ITSps.noncontam <- subset_taxa(ITSps.noncontam, Kingdom!= "Viridiplantae" &Kingdom!= "Unassigned")

# If taxa with 1 counts are removed  from https://deneflab.github.io/Diversity_Productivity/analysis/OTU_Removal_Analysis.html
ITSps.noncontam<-prune_taxa(taxa_sums(ITSps.noncontam) > 1, ITSps.noncontam) 

#LOOK AT READ COUNTS
pssamnonc<-as.data.frame(ITSps.noncontam@sam_data)
ITSotutblnonc<-as.data.frame(ITSps.noncontam@otu_table)
ITStaxtblnonc<-as.data.frame(ITSps.noncontam@tax_table)
ITSotu_tbl_rowsums<-as.data.frame(rowSums(ITSotutblnonc))
ITSotu_tbl_colsums<-as.data.frame(colSums(ITSotutblnonc))

#remove ITS_0722_A_D7_3_CA_r_contig.gz,ITS_1222_A_DI_1_AZ_r_contig.gz, ITS_0522_A_H_3_CA_R2_contig.gz
#NOTE: these samples are not included in the NCBI submission
ITSps.noncontam = subset_samples(ITSps.noncontam, sample != "ITS_0522_A_H_3_CA_R2" )
ITSps.noncontam = subset_samples(ITSps.noncontam, sample != "ITS_1222_A_DI_1_AZ_r")
ITSps.noncontam = subset_samples(ITSps.noncontam, sample!= "ITS_0722_A_D7_3_CA_r")

#LOOK AT READ COUNTS
pssamnonc<-as.data.frame(ITSps.noncontam@sam_data)
ITSotutblnonc<-as.data.frame(ITSps.noncontam@otu_table)
ITStaxtblnonc<-as.data.frame(ITSps.noncontam@tax_table)
ITSotu_tbl_rowsums<-as.data.frame(rowSums(ITSotutblnonc))
ITSotu_tbl_colsums<-as.data.frame(colSums(ITSotutblnonc))
ITSps.noncontam<-prune_taxa(taxa_sums(ITSps.noncontam) > 1, ITSps.noncontam) 

#Look AGAIN- all good?
ITSotutblnonc<-as.data.frame(ITSps.noncontam@otu_table)
ITStaxtblnonc<-as.data.frame(ITSps.noncontam@tax_table)
ITSotu_tbl_rowsums<-as.data.frame(rowSums(ITSotutblnonc))
ITSotu_tbl_colsums<-as.data.frame(colSums(ITSotutblnonc))

#subset mock samples
MOCK <- subset_samples(ITSps, Month == "mock")

#____ITS MOCK COMMUNITY ANALYSIS FOR ALL SEVEN SEQUENCING RUNS____________________

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
set.seed(1993)
MOCK.counts <- count_seqs(MOCK@otu_table)

abund.MOCK <- MOCK@otu_table %>% as.matrix() %>% as.data.frame() %>%
  rownames_to_column("OTU") %>%
  pivot_longer(-OTU, names_to = "sample") %>%
  merge(MOCK.counts$counts, by="sample") %>%
  dplyr::select(-n) %>%
  mutate(newcol=sample)%>%
  separate(newcol, c("ITS","mock","dil1", "dil2","dil3","dil4"), "_",remove = TRUE)%>%
  group_by(sample) %>%
  mutate(rel_abund = value/counts) %>% #relative abundance by treatment/day/rep
  ungroup() 

#check that relative abundance is calculated correctly
abund.MOCK %>% group_by(sample) %>% summarise(sum = sum(rel_abund))

#get taxonomy
tax <- MOCK@tax_table %>%
  as.data.frame() %>%
  rownames_to_column("OTU") 
head(tax)

#merge taxonomy with relative abundance
count_tax.MOCK <- merge(abund.MOCK, tax, by="OTU")
head(count_tax.MOCK)

#### Relative abundance plot ####

# Get summary abundance of each taxa
taxa.meta.sumMOCK <- count_tax.MOCK%>%
  dplyr::select(-c(value, counts)) %>%
  pivot_longer(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "OTU"),
               names_to = "level", 
               values_to = "taxon")

#format phylum names
genus_rel_abundMOCK <- taxa.meta.sumMOCK %>%
  filter(level == "Genus") %>%
  group_by(taxon, sample,mock,dil1, dil2,dil3,dil4) %>% 
  summarise(rel_abund = 100*sum(rel_abund), .groups="drop") %>%
  #group_by(land, day, taxon) %>%
  #summarize(mean_rel_abund = mean(rel_abund) *100, .groups="drop") %>%
  mutate(taxon = str_replace(taxon,
                             "^unclassified (.*)", "Unclassified *\\1*"),
         taxon = str_replace(taxon, "^(\\S*)$", "\\1")) %>%
  ungroup()
genus_rel_abundMOCK$taxon <-genus_rel_abundMOCK$taxon%>%replace_na('Unclassified')

#pool taxon <1% relative abundance into "other" category
taxon_poolMOCK <- genus_rel_abundMOCK %>%
  group_by(sample, taxon) %>%
  summarise(mean = mean(rel_abund), .groups="drop") %>%
  group_by(taxon) %>%
  summarise(pool = max(mean) <2,
            mean=mean(mean), 
            .groups="drop") 

#arrange data frame to group samples by relative abundance within each taxon; 
#this will aid visual comparison between the barplots
join_taxon_rel_abundMOCK <- inner_join(genus_rel_abundMOCK, taxon_poolMOCK, by="taxon") %>%
  mutate(taxon = ifelse(pool, "< 2 %", taxon)) %>%
  group_by(sample, taxon,mock,dil1, dil2,dil3,dil4) %>%
  summarise(rel_abund = sum(rel_abund),
            mean = mean(mean), .groups="drop") %>%
  mutate(taxon = factor(taxon),
         taxon = fct_reorder(taxon, rel_abund, .desc=TRUE),
         taxon = fct_shift(taxon, n=1)) %>%
  ungroup() %>%
  unique() %>%
  group_by(taxon) %>%
  arrange(desc(rel_abund), .by_group=TRUE) 

levels(join_taxon_rel_abundMOCK $taxon)

color_map <- c( "< 2 %"="black","Alternaria" = "red", "Fusarium" = "blue","Aureobasidium"="tan1", 
                "Cladosporium"= "plum1", "Cryptococcus"="darkgoldenrod1", "Neurospora"="orangered2", 
                "Curvularia"= "darkolivegreen2", "Plectosphaerella"="chocolate", "Pleiochaeta"="magenta",
                "Stemphylium"= "gray28", "Vishniacozyma"= "skyblue"," Penicillium"="darkorchid4", 
                "Sporobolomyces"= "gray64", "Botrytis"="green4", "Gibellulopsis"="papayawhip",
                "Neocamarosporium"= "rosybrown1", "Phoma"= "thistle3", "Stereum"= "yellow",
                "Cercospora"="navy", "Saccharomyces"="purple", "Cutaneotrichosporon"="pink", 
                "Pleosporales_gen_Incertae_sedis"= "palevioletred", "Aspergillus"= "sienna",  
                "Nakaseomyces"= "deeppink", "Candida"="aquamarine", "Fungi_gen_Incertae_sedis"= "grey99")
join_taxon_rel_abundMOCK  %>%
  ggplot()+
  geom_col(mapping = aes(x = sample, y = rel_abund, fill = taxon), color = "black", position = "fill", show.legend = TRUE, width=.4)+
  #facet_grid(cols = vars(WMU),scales="free")+
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


#make a chart w/ rel ab for fungi mock samples
genus_rel_abundMOCK <- genus_rel_abundMOCK %>% filter(rel_abund!="0") 
genus_rel_abundMOCK<-dcast(genus_rel_abundMOCK,taxon~sample,value.var="rel_abund")

# Create variables for the Defined Composition of the ITS mock community
Genus <- c("Aspergillus", "Cryptococcus","Trichophyton" , "Penicillium" ,"Fusarium","Malassezia","Saccharomyces","Cutaneotrichosporon","Candida","Nakaseomyces")
Original_Composition <- c(10.0, 10.0, 10.0, 10.0,10.0, 10.0,10.0, 10.0,10.0,10.0)

# Join the variables to create a data frame
original <- data.frame(Genus,Original_Composition)
colnames(genus_rel_abundMOCK)[1] <-"Genus"
totalmock_orig <- merge(original,genus_rel_abundMOCK,by=c("Genus"), all=T)
totalmock <-totalmock_orig%>%arrange(desc(Original_Composition))%>%replace(is.na(.), 0)
totalmock<-totalmock%>%select(1,2,9,10,11,12,13,14,15)
totalmock<-totalmock[, c(1,2,9,3,4,5,6,7,8)]
colnames(totalmock)[3] <-"Sequencing Run 1 1:100 Dilution Relative Abundance"
colnames(totalmock)[4] <-"Sequencing Run 2 1:100 Dilution Relative Abundance"
colnames(totalmock)[5] <-"Sequencing Run 3 1:100 Dilution Relative Abundance"
colnames(totalmock)[6] <-"Sequencing Run 4 1:100 Dilution Relative Abundance"
colnames(totalmock)[7] <-"Sequencing Run 5 1:100 Dilution Relative Abundance"
colnames(totalmock)[8] <-"Sequencing Run 6 1:100 Dilution Relative Abundance"
colnames(totalmock)[9] <-"Sequencing Run 7 1:100 Dilution Relative Abundance"
colnames(totalmock)[2] <-"Defined_Composition"
totalmock<-totalmock%>%filter(Defined_Composition!="0")
write.csv(totalmock, "/local1/workdir1/tw488/CIDA/For_Publication/ITS/mockITSsampletable_forsupp_03072025.csv", row.names=FALSE)

#_____MAKE Phyloseq OBJECTS FOR EACH ANALYSIS___________________________________________________________

#REMOVE NC SAMPLES, samples being excluded AND rarefy ITSps.noncontam
ITSps.noncontam_noNC <- subset_samples(ITSps.noncontam, Rep!= "NC")
#REMOVE 0722_B
ITSps.noncontam_noNC_exclsamples <- subset_samples(ITSps.noncontam_noNC, sample!= "ITS_0722_B_H_3_CA")
#REMOVE 0922_B shelf life samples (Harvest is ok)
ITSps.noncontam_noNC_exclsamples <- subset_samples(ITSps.noncontam_noNC_exclsamples, sample!= "ITS_0922_B_DI_2_CA")
ITSps.noncontam_noNC_exclsamples <- subset_samples(ITSps.noncontam_noNC_exclsamples, sample!= "ITS_0922_B_DI_3_CA")
ITSps.noncontam_noNC_exclsamples <- subset_samples(ITSps.noncontam_noNC_exclsamples, sample!= "ITS_0922_B_D7_3_CA")
ITSps.noncontam_noNC_exclsamples <- subset_samples(ITSps.noncontam_noNC_exclsamples, sample!= "ITS_0922_B_D7_2_CA")
ITSps.noncontam_noNC_exclsamples <- subset_samples(ITSps.noncontam_noNC_exclsamples, sample!= "ITS_0922_B_D22_2_CA")
ITSps.noncontam_noNC_exclsamples <- subset_samples(ITSps.noncontam_noNC_exclsamples, sample!= "ITS_0922_B_D22_1_CA")
ITSps.noncontam_noNC_exclsamples <- subset_samples(ITSps.noncontam_noNC_exclsamples, sample!= "ITS_0922_B_D17_3_CA")
ITSps.noncontam_noNC_exclsamples <- subset_samples(ITSps.noncontam_noNC_exclsamples, sample!= "ITS_0922_B_D17_1_CA")
ITSps.noncontam_noNC_exclsamples <- subset_samples(ITSps.noncontam_noNC_exclsamples, sample!= "ITS_0922_B_D12_2_CA")
ITSps.noncontam_noNC_exclsamples <- subset_samples(ITSps.noncontam_noNC_exclsamples, sample!= "ITS_0922_B_D12_1_CA")

#remove the third 0422_B sample
ITSps.noncontam_noNC_exclsamples <- subset_samples(ITSps.noncontam_noNC_exclsamples, sample!= "ITS_0422_B_H_1_CA")
ITSps.noncontam_noNC.rarefied = rarefy_even_depth(ITSps.noncontam_noNC_exclsamples, rngseed=1, sample.size=min(sample_sums(ITSps.noncontam_noNC)), replace=F)
ITSps.noncontam_noNC.rarefiedotu<-ITSps.noncontam_noNC.rarefied@otu_table
ITSps.noncontam_noNC.rarefiedotu_colsums<-as.data.frame(colSums(ITSps.noncontam_noNC.rarefiedotu))
ITSps.noncontam_noNC.rarefiedotu_rowsums<-as.data.frame(rowSums(ITSps.noncontam_noNC.rarefiedotu))

#MAKE PS OBJECTS FOR EACH ANALYSIS
#HARVEST AND DAY INITIAL ONLY
HDI <- subset_samples(ITSps.noncontam_noNC.rarefied,Day!= "D22")
HDI <- subset_samples(HDI,Day!= "D17")
HDI <- subset_samples(HDI,Day!= "D12")
HDI <- subset_samples(HDI,Day!= "D7")
#HARVEST ONLY
H <- subset_samples(HDI,Day!= "DI")
#SHELF LIFE SAMPLES
DI_D22_D28 <- subset_samples(ITSps.noncontam_noNC.rarefied,Day!= "H")

# SEVEN DAY INTERVAL SHELF LIFE SAMPLES
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
sevenday <- subset_samples(sevenday,sample!= "ITS_0622_B_D7_2_CA")
sevenday <- subset_samples(sevenday,sample!= "ITS_0622_B_D7_3_CA")
sevenday <- subset_samples(sevenday,sample!= "ITS_0622_B_DI_1_CA")
sevenday <- subset_samples(sevenday,sample!= "ITS_0622_B_DI_3_CA")
sevenday <- subset_samples(sevenday,sample!= "ITS_0622_B_H_2_CA")
sevenday <- subset_samples(sevenday,sample!= "ITS_0622_B_H_3_CA")

#FIVE DAY INTERVAL SHELF LIFE SAMPLES
fiveday <- subset_samples(DI_D22_D28,Month!= "0522")
fiveday <- subset_samples(fiveday,Month!= "0422")
fiveday <- subset_samples(fiveday,sample!= "ITS_0622_A_D12_1_CA")
fiveday <- subset_samples(fiveday,sample!= "ITS_0622_A_D12_2_CA")
fiveday <- subset_samples(fiveday,sample!= "ITS_0622_A_D17_1_CA")
fiveday <- subset_samples(fiveday,sample!= "ITS_0622_A_D17_2_CA")
fiveday <- subset_samples(fiveday,sample!= "ITS_0622_A_D22_2_CA")
fiveday <- subset_samples(fiveday,sample!= "ITS_0622_A_D22_3_CA")
fiveday <- subset_samples(fiveday,sample!= "ITS_0622_A_D7_3_CA")
fiveday <- subset_samples(fiveday,sample!= "ITS_0622_A_D7_1_CA")
fiveday <- subset_samples(fiveday,sample!= "ITS_0622_A_DI_1_CA")
fiveday <- subset_samples(fiveday,sample!= "ITS_0622_A_DI_3_CA")

#___Calculate  relative abundance for H samples_______________________________________

H.counts <- count_seqs(H@otu_table)

abund.H <- H@otu_table %>% as.matrix() %>% as.data.frame() %>%
  rownames_to_column("OTU") %>%
  pivot_longer(-OTU, names_to = "sample") %>%
  merge(H.counts$counts, by="sample") %>%
  dplyr::select(-n) %>%
  mutate(newcol=sample)%>%
  separate(newcol, c("X16S","Month", "Lot","Day","Rep","Loc","Actual_Day"), "_",remove = TRUE)%>%
  group_by(sample) %>%
  mutate(rel_abund = value/counts) %>% #relative abundance by treatment/day/rep
  ungroup() 

#check that relative abundance is calculated correctly
abund.H %>% group_by(sample) %>% summarise(sum = sum(rel_abund))

#get taxonomy
taxH <- H@tax_table %>%
  as.data.frame() %>%
  rownames_to_column("OTU") 
head(taxH)

#how many otus are there total?
nrow(taxH)

#merge taxonomy with rarefied relative abundance
count_tax.H <- merge(abund.H, taxH, by="OTU")
head(count_tax.H)

#### Relative abundance plot ####

# Get summary abundance of each taxa
taxa.meta.sumH <- count_tax.H%>%
  dplyr::select(-c(value, counts)) %>%
  pivot_longer(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "OTU"),
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
genus_rel_abundH$taxon <-gsub("Fungi_gen_Incertae_sedis", "Unclassified", genus_rel_abundH$taxon)
genus_rel_abundH$taxon <-gsub("Pleosporaceae_gen_Incertae_sedis", "Unclassified Pleosporaceae", genus_rel_abundH$taxon)
genus_rel_abundH$taxon <-gsub("Sporormiaceae_gen_Incertae_sedis", "Unclassified Sporormiaceae", genus_rel_abundH$taxon)
genus_rel_abundH$taxon <-gsub("Pleosporales_gen_Incertae_sedis", "Unclassified Pleosporales", genus_rel_abundH$taxon)

#how many genera are there total?
totalHgenera_before5<-genus_rel_abundH%>%group_by(taxon)%>%filter(rel_abund!=0)%>%distinct(taxon)

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

#Look at unclassified taxa specifically for the mean, max and min
uncl<-join_taxon_rel_abundH%>%filter(taxon=="Unclassified")

#MAKE A LIST OF ALL GENERA IN HARVEST SAMPLES >5% RELATIVE ABUNDANCE
totalHgenera<-join_taxon_rel_abundH%>%group_by(taxon)%>%distinct(taxon)


#make different columns: monthlot, month-lot-day-rep (MLDR) and day-rep (DR)
join_taxon_rel_abundH$Monthlot = paste(join_taxon_rel_abundH$Month, join_taxon_rel_abundH$Lot, sep="_")
join_taxon_rel_abundH$MLDR = paste(join_taxon_rel_abundH$Month,join_taxon_rel_abundH$Lot,join_taxon_rel_abundH$Day,join_taxon_rel_abundH$Rep,sep="_")
join_taxon_rel_abundH$DR = paste(join_taxon_rel_abundH$Day,join_taxon_rel_abundH$Rep,sep="_")
join_taxon_rel_abundH$Month<-factor(join_taxon_rel_abundH$Month,
                                    levels=c("1221","0122","0222","0322","0422","0522","0622","0722","0822","0922","1022","1122","1222"))

join_taxon_rel_abundH$DR<-factor(join_taxon_rel_abundH$DR,
                                 levels=c("H_1","H_2","H_3", "DI_1","DI_2","DI_3", "D7_1","D7_2","D7_3", "D12_1","D12_2","D12_3", "D17_1","D17_2","D17_3","D22_1","D22_2","D22_3"))

color_map<- c( "< 5 %"="black","Alternaria" = "grey86", "Fungi_gen_Incertae_sedis" = "blue",
               "Aureobasidium"="chartreuse", "Cladosporium"= "darkorchid4", "Cryptococcus"="darkgoldenrod1",
               "Neurospora"="orange", "Curvularia"= "darkolivegreen", 
               "Unclassified Pleosporales"="sienna", "Pleiochaeta"="green4", 
               "Stemphylium"= "yellow", "Vishniacozyma"= "skyblue", "Sporobolomyces"= "gray64", 
               "Unclassified Pleosporaceae"="magenta", "Comoclathris"="palevioletred", 
               "Neocamarosporium"= "rosybrown3", "Unclassified Sporormiaceae"= "thistle3",
               "Celosporium"="aquamarine","Comoclathris"= "yellow3",
               "Fusarium"= "navy","Epicoccum"= "darkcyan","Botrytis"= "deeppink",
               "Unclassified Sporormiaceae"="firebrick","Unclassified"="red")


#MAKE THE BAR PLOT
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

ggsave("Fungal_RelAB_Barplot_05302025.tiff", units="in", width=20, height=15, dpi=500, compression = 'lzw')


#HARVEST SAMPLE DATA ANALYSIS_____________________________________________________________________________________

#PERMANOVA
Hotu<-pstoveg_otu(H)
Hsd<-pstoveg_sd(H)
Hsd$Monthlot = paste(Hsd$Month, Hsd$Lot, sep="_")
Hvegdist<-vegdist(Hotu, method="bray") 

#make Monthlot in to a plot
h <- with(Hsd, how(within=Within(type="free"),plots=Plots(strata=Monthlot, type="free"),nperm=999))
h
#check to make sure all months have 2 samples
table(Hsd$Monthlot,Hsd$Loc)
H_perm <- adonis2(Hvegdist ~Location, data=Hsd, permutations = h, method = "bray", 
                  sqrt.dist = FALSE, add = FALSE, by = "terms",
                  parallel = getOption("mc.cores"), na.action = na.fail)
H_perm

#Non-metric Multidimensional Scaling (NMDS) for H samples_________________________________________________________

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

plot_H_nmds <- ggplot(plot_dfH, aes(x = NMDS1, y = NMDS2, color = Month, shape = Location)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = palNMDS) +
  stat_ellipse(aes(color=Location), linetype=2, level=0.5)+
  theme_classic()
plot_H_nmds

ggsave("H_FUNGAL_NMDS_04292025.tiff", units="in", width=8, height=6, dpi=300, compression = 'lzw')

#_________remove FL from H PS object_________________________________________________________
H_noFL <- subset_samples(H,Location!= "FL")
HOTU<-as.data.frame(H_noFL@otu_table)
HOTU_colsums<-as.data.frame(colSums(HOTU))
HTAX<-data.frame(H_noFL@tax_table)%>%tibble::rownames_to_column("OTU")
Hsampledata<-data.frame(H_noFL@sam_data)
Hsampledata$Monthlot = paste(Hsampledata$Month, Hsampledata$Lot, sep="_")

#maaslin2 for ASV level in H samples
Hrarefit_ITS_CA = Maaslin2(input_data = HOTU, 
                             input_metadata = Hsampledata,
                             analysis_method = "CPLM",
                             transform= "NONE",
                             min_prevalence = 0.1,
                             min_abundance = 0.01,
                             normalization  = "NONE",
                             output         = "/local1/workdir1/tw488/CIDA/ITS_Maaslin2_Output/H_ASV_ITS_output_08022024", 
                             fixed_effects  = "Location",
                             reference      = c("Location,CA"),
                             max_pngs = 200,
                             save_models="TRUE",
                             save_scatter="TRUE",
                             cores=2,
                             random_effects = c("Monthlot"))

Hsig_results<-read.table("/local1/workdir1/tw488/CIDA/ITS_Maaslin2_Output/H_ASV_ITS_output_08022024/significant_results.tsv",header=TRUE)
names(Hsig_results)[1]<-paste("OTU") 
Hsig_results1<-Hsig_results%>%filter(qval<0.05)%>%merge(HTAX, by='OTU', all = TRUE)%>%filter(qval!="NA") 
Hsiggenera<-c(Hsig_results1$Genus)

#TOTAL ENRICHED ASVS PER LOCATION

#total enrched ASVs in AZ
H_ASV_sig_results_enrichAZ<-Hsig_results1%>%filter(coef>0)%>% group_by(Genus) %>% count(Genus)
names(H_ASV_sig_results_enrichAZ)[2]<-paste("Number of OTUs Differentially Abundant in AZ")
#total enrched ASVs in CA
H_ASV_sig_results_enrichCA<-Hsig_results1%>%filter(coef<0)%>% group_by(Genus) %>% count(Genus)
names(H_ASV_sig_results_enrichCA)[2]<-paste("Number of OTUs Differentially Abundant in CA")
#combine the  enriched ASVs in to one
H_OTUcount_byArea <- merge(H_ASV_sig_results_enrichCA,H_ASV_sig_results_enrichAZ, by="Genus", all=TRUE)

#TOTAL ASVS PER GENUS
HTAX$Genus<-HTAX$Genus %>% replace_na('Unclassified')
HTAX_taxasum<-HTAX%>% group_by(Genus) %>% count(Genus)
HTAX_taxasum1<-HTAX_taxasum%>%filter(Genus %in% Hsiggenera)
H_OTUcount_byArea <- merge(HTAX_taxasum1,H_OTUcount_byArea, by="Genus", all=TRUE)
names(H_OTUcount_byArea)[2]<-paste("Total Number of OTUs assigned to Genus")
#write.csv(H_OTUcount_byArea,"/local1/workdir1/tw488/CIDA/ITS_Maaslin2_Output/H_ASV_ITS_output_08022024/ITS_OTU_H_Table", row.names = TRUE,na = "0")

#Calculate relative abundance values from H phyloseq object w/o FL______________________

HnoFL.counts <- count_seqs(H_noFL@otu_table)

abund.H_noFL <- H_noFL@otu_table %>% as.matrix() %>% as.data.frame() %>%
  rownames_to_column("OTU_") %>%
  pivot_longer(-OTU_, names_to = "sample") %>%
  merge(HnoFL.counts$counts, by="sample") %>%
  dplyr::select(-n) %>%
  mutate(newcol=sample)%>%
  separate(newcol, c("X16S","Month", "Lot","Day","Rep","Loc","Actual_Day"), "_",remove = TRUE)%>%
  group_by(sample) %>%
  mutate(rel_abund = value/counts) %>% #relative abundance by treatment/day/rep
  ungroup() 

#check that relative abundance is calculated correctly
abund.H_noFL %>% group_by(sample) %>% summarise(sum = sum(rel_abund))

#get taxonomy
taxH.noFL <- H_noFL@tax_table %>%
  as.data.frame() %>%
  rownames_to_column("OTU_") 
head(taxH.noFL)

#merge taxonomy with rarefied relative abundance
count_tax.H_noFL <- merge(abund.H_noFL, taxH.noFL, by="OTU_")
head(count_tax.H_noFL)

# Get summary abundance of each taxa
taxa.meta.sumH_noFL <- count_tax.H_noFL%>%
  dplyr::select(-c(value, counts)) %>%
  pivot_longer(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "OTU_"),
               names_to = "level", 
               values_to = "taxon") 
#format phylum names
genus_rel_abundH_noFL <- taxa.meta.sumH_noFL %>%
  filter(level == "Genus") %>%
  group_by(taxon, sample,Month, Lot,Day,Rep,Loc,Actual_Day) %>% 
  summarise(rel_abund = 100*sum(rel_abund), .groups="drop") %>%
  mutate(taxon = str_replace(taxon,
                             "^unclassified (.*)", "Unclassified *\\1*"),
         taxon = str_replace(taxon, "^(\\S*)$", "\\1")) %>%
  ungroup()
genus_rel_abundH_noFL$taxon <-genus_rel_abundH_noFL$taxon%>%replace_na('Unclassified')
genus_rel_abundH_noFL_4maaslin<-dcast(genus_rel_abundH_noFL,sample~taxon,value.var="rel_abund")
genus_rel_abundH_noFL_4maaslin<-genus_rel_abundH_noFL_4maaslin%>% column_to_rownames(var="sample")
class(genus_rel_abundH_noFL_4maaslin)

#ITS diversity metrics by time into planting season_____________________________________

#DFSOS= day from start of season
#DFSOS v2
DFSOSv2<-read.csv("/local1/workdir1/tw488/CIDA_Sequencing_Runs/DFSOSv2.csv", header=TRUE)
DFSOSv2$first_planting_date_from_TF<-gsub("_", "-", DFSOSv2$first_planting_date_from_TF)
DFSOSv2$date_of_harvest_yyyy_mm_dd<-gsub("_", "-", DFSOSv2$date_of_harvest_yyyy_mm_dd)
DFSOSv2$first_planting_date_from_TF<-as.Date(DFSOSv2$first_planting_date_from_TF)
DFSOSv2$date_of_harvest_yyyy_mm_dd<-as.Date(DFSOSv2$date_of_harvest_yyyy_mm_dd)
DFSOSv2$time_from_first_planting<-difftime(DFSOSv2$date_of_harvest_yyyy_mm_dd,DFSOSv2$first_planting_date_from_TF,units = "days")
DFSOSv2$time_from_first_planting<-as.integer(DFSOSv2$time_from_first_planting)
DFSOSv2<-DFSOSv2%>%filter(location!="FL"&location!="GA")%>%filter(sample_type=="r")
DFSOSv2<-DFSOSv2[c("sampling_code","time_from_first_planting")]

#ALPHA DIVERSITY OF H SAMPLES
Halphadiv<-estimate_richness(H, split = TRUE, measures = NULL)
Halphadiv$sample <- row.names(Halphadiv)  
Halphadiv <- Halphadiv %>% separate(sample, c("x16s", "month", "lot", "day", "rep", "loc","r1"), "_")
Halphadiv$sample <- row.names(Halphadiv) 
Halphadiv$sampling_code = paste(Halphadiv$month, Halphadiv$lot, sep="_")
Halphadiv <- Halphadiv %>% 
  as.data.frame() %>% 
  full_join(DFSOSv2, by = "sampling_code")

#calculate Pielou's Eveness
Halphadiv <- Halphadiv %>% mutate (Pielou=(Halphadiv$Shannon)/log(Halphadiv$Observed))
#filter out FL, 0722B and 1222B
Halphadiv1 <- Halphadiv %>% filter(sampling_code!="1222_B")%>% filter(sampling_code!="0722_B")%>% filter(sampling_code!="1221_A")

#visualize
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

#COMPARING ALPHA DIVERSITY REULTS BY LOCATION_____________________________________________________

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

#____Which genera are over-represented in a specific time of planting season (maybe specific by region)?___________________________________________________________________________________________

#add DFSOS to data frames
#FILTER DFSOS WITH SAMPLES THAT WERE NOT SEQUENCED OR BEING EXCLUDED

Hsampledata$sampling_code = paste(Hsampledata$Month,Hsampledata$Lot, sep="_")
Hsampledata_dfsos<- Hsampledata %>%
  tibble::rownames_to_column("rowsample") %>% 
  inner_join(DFSOSv2, by = "sampling_code")%>%column_to_rownames( var = "rowsample")

#run maaslin2 w/ updated data frames with DFSOS for H data at genus level
#for the code below, 'sampling_code' is equivalent to the earlier 'monthlot'

H_genus_dfsos_CAref_08162024 = Maaslin2(input_data =genus_rel_abundH_noFL_4maaslin, 
                                        input_metadata = Hsampledata_dfsos,
                                        analysis_method = "CPLM",
                                        transform= "NONE",
                                        min_prevalence = 0,
                                        min_abundance = 0,
                                        normalization  = "NONE",
                                        random_effects = c("Monthlot"),
                                        output         = ("/local1/workdir1/tw488/CIDA/ITS_Maaslin2_Output/H_GENUS_ITS_dfsos_output_08162024"), 
                                        fixed_effects  = c("Location","time_from_first_planting"),
                                        reference      = c("Location,CA"),
                                        max_pngs = 200,
                                        save_models="TRUE",
                                        save_scatter="TRUE",
                                        cores=2)

#Import results
H_LocDFSOS_results<-read.table("/local1/workdir1/tw488/CIDA/ITS_Maaslin2_Output/H_GENUS_ITS_dfsos_output_08162024/significant_results.tsv",header=TRUE)
H_LocDFSOS_results0.05<-H_LocDFSOS_results%>%filter(qval<0.05)
#select columns from diff abundance data, calc log 2 fold change
H_LocDFSOS_results0.05<-H_LocDFSOS_results0.05%>%select(feature,value,coef,qval)%>% mutate(log2_fold_change = log2(exp(coef)))%>%select(feature,value,qval,log2_fold_change)

#merge w/ relative abundance data to make a table 
#take the time/dfsos data, split time in to 1st half and second half of season by dividing month associated with a location
H_LocDFSOS_results0.05_time<-H_LocDFSOS_results0.05%>%filter(value!="AZ")
colnames(H_LocDFSOS_results0.05_time)[1] <- "Genus"
colnames(H_LocDFSOS_results0.05_time)[3] <- "False Discovery Rate"

#filter out unneeded taxa
H_LocDFSOS_results0.05_time_taxon<-c(H_LocDFSOS_results0.05_time$Genus)
genus_rel_abundH_noFL_4datatbl<-genus_rel_abundH_noFL%>%filter(taxon %in% H_LocDFSOS_results0.05_time_taxon)

#separate in to CA and AZ
genus_rel_abundH_noFL_4datatbl_sampNameCA<-genus_rel_abundH_noFL_4datatbl%>%filter(Loc=="CA")%>%mutate(whichhalfofszn="First")
genus_rel_abundH_noFL_4datatbl_sampNameAZ<-genus_rel_abundH_noFL_4datatbl%>%filter(Loc=="AZ")%>%mutate(whichhalfofszn="First")

#add in first or second in whichhalfofszn column
genus_rel_abundH_noFL_4datatbl_sampNameCA$whichhalfofszn[genus_rel_abundH_noFL_4datatbl_sampNameCA$Month == "0722"] <- "Second"
genus_rel_abundH_noFL_4datatbl_sampNameCA$whichhalfofszn[genus_rel_abundH_noFL_4datatbl_sampNameCA$Month == "0822"] <- "Second"
genus_rel_abundH_noFL_4datatbl_sampNameCA$whichhalfofszn[genus_rel_abundH_noFL_4datatbl_sampNameCA$Month == "0922"] <- "Second"
genus_rel_abundH_noFL_4datatbl_sampNameCA$whichhalfofszn[genus_rel_abundH_noFL_4datatbl_sampNameCA$Month == "1022"] <- "Second"

#now AZ
genus_rel_abundH_noFL_4datatbl_sampNameAZ$whichhalfofszn[genus_rel_abundH_noFL_4datatbl_sampNameAZ$Month == "0222"] <- "Second"
genus_rel_abundH_noFL_4datatbl_sampNameAZ$whichhalfofszn[genus_rel_abundH_noFL_4datatbl_sampNameAZ$Month == "0322"] <- "Second"
genus_rel_abundH_noFL_4datatbl_sampNameAZ$whichhalfofszn[genus_rel_abundH_noFL_4datatbl_sampNameAZ$Month == "0422"] <- "Second"

#combine the two areas and get average relative abundance bywhichhalfofszn
genus_rel_abundH_noFL_4datatbl_combo<-rbind(genus_rel_abundH_noFL_4datatbl_sampNameCA,genus_rel_abundH_noFL_4datatbl_sampNameAZ)
genus_rel_abundH_noFL_4datatbl_combo1<-genus_rel_abundH_noFL_4datatbl_combo%>%group_by(taxon,whichhalfofszn)%>%mutate(AvybyWhos=mean(rel_abund))%>%distinct(taxon,whichhalfofszn,AvybyWhos)
genus_rel_abundH_noFL_4datatbl_combo1<-dcast(genus_rel_abundH_noFL_4datatbl_combo1,taxon~whichhalfofszn,value.var="AvybyWhos")

#reformat relative abundance data and merge with maaslin2 output
colnames(genus_rel_abundH_noFL_4datatbl_combo1)[1] <- "Genus"
colnames(genus_rel_abundH_noFL_4datatbl_combo1)[2] <- "Relative Abundance for First Half of Planting Season (%)"
colnames(genus_rel_abundH_noFL_4datatbl_combo1)[3] <- "Relative Abundance for Second Half of Planting Season (%)"
genus_rel_abundH_noFL_4datatbl_combo2<-merge(genus_rel_abundH_noFL_4datatbl_combo1,H_LocDFSOS_results0.05_time, by="Genus")
genus_rel_abundH_noFL_4datatbl_combo2<-genus_rel_abundH_noFL_4datatbl_combo2[,-c(4)]
#This data was inlcuded in the results text, but not as a table
#write.csv(genus_rel_abundH_noFL_4datatbl_combo2,"/local1/workdir1/tw488/CIDA/ITS_Maaslin2_Output/H_GENUS_ITS_dfsos_output_08162024/Time_from_start_of_season_Maaslin_results_with_RelAb", row.names = TRUE)

#Area/Location data
H_LocDFSOS_results0.05_AZ<-H_LocDFSOS_results0.05%>%filter(value=="AZ")
colnames(H_LocDFSOS_results0.05_AZ)[1] <- "Genus"
colnames(H_LocDFSOS_results0.05_AZ)[3] <- "False Discovery Rate"

#filter out unneeded taxa
H_LocDFSOS_results0.05_AZ_taxon<-c(H_LocDFSOS_results0.05_AZ$Genus)
genus_rel_abundH_noFL__loc_4datatbl<-genus_rel_abundH_noFL%>%filter(taxon %in% H_LocDFSOS_results0.05_AZ_taxon)

#separate in to CA and AZ
genus_rel_abundH_noFL__loc_4datatbl_sampNameCA<-genus_rel_abundH_noFL__loc_4datatbl%>%filter(Loc=="CA")
genus_rel_abundH_noFL__loc_4datatbl_sampNameAZ<-genus_rel_abundH_noFL__loc_4datatbl%>%filter(Loc=="AZ")

#relative abundance for CA
genus_rel_abundH_noFL__loc_4datatbl_sampNameCA<-genus_rel_abundH_noFL__loc_4datatbl_sampNameCA%>%group_by(taxon)%>%mutate(CA_Avg=mean(rel_abund))%>%distinct(taxon,CA_Avg)

#relative abundance for AZ
genus_rel_abundH_noFL__loc_4datatbl_sampNameAZ<-genus_rel_abundH_noFL__loc_4datatbl_sampNameAZ%>%group_by(taxon)%>%mutate(AZ_Avg=mean(rel_abund))%>%distinct(taxon,AZ_Avg)

#combine the two areas and get average relative abundance bywhichhalfofszn
genus_rel_abundH_noFL__loc_4datatbl_combo<-merge(genus_rel_abundH_noFL__loc_4datatbl_sampNameAZ,genus_rel_abundH_noFL__loc_4datatbl_sampNameCA, by="taxon")

# merge with maaslin2 output
colnames(genus_rel_abundH_noFL__loc_4datatbl_combo)[1] <- "Genus"
genus_rel_abundH_noFL__loc_4datatbl_combo1<-merge(genus_rel_abundH_noFL__loc_4datatbl_combo,H_LocDFSOS_results0.05_AZ, by="Genus")
colnames(genus_rel_abundH_noFL__loc_4datatbl_combo1)[2] <- "AZ Relative Abundance (%)"
colnames(genus_rel_abundH_noFL__loc_4datatbl_combo1)[3] <- "CA Relative Abundance (%)"
genus_rel_abundH_noFL__loc_4datatbl_combo1 <- genus_rel_abundH_noFL__loc_4datatbl_combo1[, -4] 
#write.csv(genus_rel_abundH_noFL__loc_4datatbl_combo1,"/local1/workdir1/tw488/CIDA/ITS_Maaslin2_Output/H_GENUS_ITS_dfsos_output_08162024/Location_Maaslin_results_with_RelAb", row.names = TRUE)

#_______________H VS DI______________________________________________________________________

#Find most common genera (TOP 5) in DI samples
DI <- subset_samples(HDI,Day!= "H")

#### Calculate relative abundance ####

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


#H_DI PERMANOVA_________________________________________________________

#sort out HDI data, for permanova 
HDIotu<-pstoveg_otu(HDI)
HDIsd<-pstoveg_sd(HDI)
HDIsd$Monthlot = paste(HDIsd$Month, HDIsd$Lot, sep="_")
HDIvegdist<-vegdist(HDIotu, method="bray") 

#remove 0122A, 0822B,0922B,1122A, 1222B, 1222C b/c they only have two samples (either H or DI)
HDI_trim<-subset_samples(HDI,sample!= "ITS_0122_A_DI_1_AZ")
HDI_trim<-subset_samples(HDI_trim,sample!= "ITS_0122_A_DI_2_AZ")

HDI_trim<-subset_samples(HDI_trim,sample!= "ITS_0822_B_DI_2_CA")
HDI_trim<-subset_samples(HDI_trim,sample!= "ITS_0822_B_DI_3_CA")

HDI_trim<-subset_samples(HDI_trim,sample!= "ITS_0922_B_H_1_CA")
HDI_trim<-subset_samples(HDI_trim,sample!= "ITS_0922_B_H_2_CA")

HDI_trim<-subset_samples(HDI_trim,sample!= "ITS_1122_A_DI_2_CA")
HDI_trim<-subset_samples(HDI_trim,sample!= "ITS_1122_A_DI_3_CA")

HDI_trim<-subset_samples(HDI_trim,sample!= "ITS_1222_B_DI_1_AZ")
HDI_trim<-subset_samples(HDI_trim,sample!= "ITS_1222_B_DI_3_AZ")

HDI_trim<-subset_samples(HDI_trim,sample!= "ITS_1222_C_H_2_AZ")
HDI_trim<-subset_samples(HDI_trim,sample!= "ITS_1222_C_H_3_AZ")

HDI_trimotu<-pstoveg_otu(HDI_trim)
HDI_trimsd<-pstoveg_sd(HDI_trim)
HDI_trimsd$Monthlot = paste(HDI_trimsd$Month, HDI_trimsd$Lot, sep="_")
HDI_trimvegdist<-vegdist(HDI_trimotu, method="bray") 

#check to make sure all months have 4 samples
table(HDI_trimsd$Monthlot,HDI_trimsd$Loc)

#make Monthlot in to a plot
hdi <- with(HDI_trimsd, how(within=Within(type="free"),plots=Plots(strata=Monthlot, type="free"),nperm=999))
hdi

#Run PERMANOVA
HDI_perm <- adonis2(HDI_trimvegdist ~Location*Day, data=HDI_trimsd, permutations = hdi, method = "bray", 
                    sqrt.dist = FALSE, add = FALSE, by = "terms",
                    parallel = getOption("mc.cores"), na.action = na.fail)
HDI_perm


#pairwise PERMANOVA

HDI_perm_pairwise<-pairwise.adonis2(HDI_trimvegdist~Location, data=HDI_trimsd)#strata='Monthlot'
HDI_perm_pairwise

##Non-metric Multidimensional Scaling (NMDS)_________________________________________________________

HDI_NMDS <- metaMDS(HDIvegdist)
HDI_NMDS
stressplot(HDI_NMDS)

HDIsd<-HDIsd %>% rename(SampleName = sample)

HDIplot_df <- scores(HDI_NMDS, display = "sites") %>% 
  as.data.frame() %>% 
  rownames_to_column("sample") %>% 
  full_join(HDIsdNMDS, by = "sample")

HDIplot_df$Day<-factor(HDIplot_df$Day,
                          levels=c("H","DI"))

palNMDS_HDI <- c( "#1F78B4", "#B2DF8A","black","grey","grey39")

plot_HDI_nmds2 <- ggplot(HDIplot_df, aes(x = NMDS1, y = NMDS2, color = Day, shape = Location)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = palNMDS_HDI) +
  stat_ellipse(linetype=2, level=0.5)+
  theme_classic()
plot_HDI_nmds2
ggsave("HDI_FUNGAL_NMDS_04292025v2.tiff", units="in", width=8, height=6, dpi=300, compression = 'lzw')

#HDI ALPHA DIVERSITY

#calculate alpha diversity stats
HDIrichnessstats<-estimate_richness(HDI, split = TRUE, measures = NULL)
HDIrichnessstats$sample <- row.names(HDIrichnessstats)  
HDIrichnessstats <- HDIrichnessstats %>% separate(sample, c("x16s", "month", "lot", "day", "rep", "loc","R2","r1"), "_")

#calculate Pielou
HDIrichnessstats <- HDIrichnessstats %>% mutate (Pielou=(HDIrichnessstats$Shannon)/log(HDIrichnessstats$Observed))

ShTTEST<-t.test(Shannon ~ day, data = HDIrichnessstats)
ShTTEST #p=7.255e-06

ObsTTEST<-t.test(Observed ~ day, data = HDIrichnessstats)
ObsTTEST #p=1.336e-07

PielOUTTEST<-t.test(Pielou ~ day, data = HDIrichnessstats)
PielOUTTEST #p=0.002245

#Visualize
#shannon
HDIrichnessstats%>% 
  ggplot(aes(x=day, y=Shannon)) +
  geom_boxplot(size=0.25) +
  labs(x="Day",
       y="Shannon Index") +
  theme_classic()
#observed
HDIrichnessstats%>% 
  ggplot(aes(x=day, y=Observed)) +
  geom_boxplot(size=0.25) +
  labs(x="Day",
       y="Richness") +
  theme_classic()
#Pielou
HDIrichnessstats%>% 
  ggplot(aes(x=day, y=Pielou)) +
  geom_boxplot(size=0.25) +
  labs(x="Day",
       y="Pielou's Evenness") +
  theme_classic()

#Genera that show the greatest shift in relative abundance between H and DI
#sort out HDI data, remove FL for maaslin2

HDI_noFL<-subset_samples(HDI,Location!= "FL")
HDI_noFLotu<-pstoveg_otu(HDI_noFL)
HDI_noFLsd<-pstoveg_sd(HDI_noFL)
HDI_noFLsd$Monthlot = paste(HDI_noFLsd$Month, HDI_noFLsd$Lot, sep="_")
HDI_noFLvegdist<-vegdist(HDI_noFLotu, method="bray") 

#### Calculate relative abundance w/o FL ####

HDI_noFL.counts <- count_seqs(HDI_noFL@otu_table)

abund.HDI_noFL <- HDI_noFL@otu_table %>% as.matrix() %>% as.data.frame() %>%
  rownames_to_column("OTU_") %>%
  pivot_longer(-OTU_, names_to = "sample") %>%
  merge(HDI_noFL.counts$counts, by="sample") %>%
  dplyr::select(-n) %>%
  mutate(newcol=sample)%>%
  separate(newcol, c("X16S","Month", "Lot","Day","Rep","Loc","Actual_Day"), "_",remove = TRUE)%>%
  group_by(sample) %>%
  mutate(rel_abund = value/counts) %>% #relative abundance by treatment/day/rep
  ungroup() 

#check that relative abundance is calculated correctly
abund.HDI_noFL %>% group_by(sample) %>% summarise(sum = sum(rel_abund))

#get taxonomy
taxHDI_noFL <- HDI_noFL@tax_table %>%
  as.data.frame() %>%
  rownames_to_column("OTU_") 
head(taxHDI_noFL)

#merge taxonomy with rarefied relative abundance
count_tax.HDI_noFL <- merge(abund.HDI_noFL, taxHDI_noFL, by="OTU_")
head(count_tax.HDI_noFL)

# Get summary abundance of each taxa
taxa.meta.sumHDI_noFL <- count_tax.HDI_noFL%>%
  dplyr::select(-c(value, counts)) %>%
  pivot_longer(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "OTU_"),
               names_to = "level", 
               values_to = "taxon") 
#format phylum names
genus_rel_abundHDI_noFL<- taxa.meta.sumHDI_noFL %>%
  filter(level == "Genus") %>%
  group_by(taxon, sample,Month, Lot,Day,Rep,Loc,Actual_Day) %>% 
  summarise(rel_abund = 100*sum(rel_abund), .groups="drop") %>%
  mutate(taxon = str_replace(taxon,
                             "^unclassified (.*)", "Unclassified *\\1*"),
         taxon = str_replace(taxon, "^(\\S*)$", "\\1")) %>%
  ungroup()
genus_rel_abundHDI_noFL$taxon <-genus_rel_abundHDI_noFL$taxon%>%replace_na('Unclassified')
genus_rel_abundHDI_noFL_4maaslin<-dcast(genus_rel_abundHDI_noFL,sample~taxon,value.var="rel_abund")
genus_rel_abundHDI_noFL_4maaslin<-genus_rel_abundHDI_noFL_4maaslin%>% column_to_rownames(var="sample")

HDI_noFLL_samdata <- data.frame(HDI_noFL@sam_data)
HDI_noFLL_samdata$Monthlot = paste(HDI_noFLL_samdata$Month, HDI_noFLL_samdata$Lot, sep="_")
HDI_noFLL_samdata <-HDI_noFLL_samdata%>%select(-R2,-SampleType,-is.neg)
class(HDI_noFLL_samdata$Day)
HDI_noFLL_samdata$Day<-as.factor(HDI_noFLL_samdata$Day)

#run maaslin2 for HDI at genus level
HDI_genus_ITS_output_11162024 = Maaslin2(input_data = genus_rel_abundHDI_noFL_4maaslin, 
                                         input_metadata = HDI_noFLL_samdata,
                                         analysis_method = "CPLM",
                                         transform= "NONE",
                                         min_prevalence = 0,
                                         min_abundance = 0,
                                         max_pngs = 200,
                                         save_models="TRUE",
                                         save_scatter="TRUE",
                                         cores=2,
                                         normalization  = "NONE",
                                         random_effects = c("Monthlot"),
                                         output         = "/local1/workdir1/tw488/CIDA/ITS_Maaslin2_Output/HDI_genus_ITS_output_11162024", 
                                         fixed_effects  = c("Day"),
                                         reference      = c("Day,H")) 


#Import results
ITS_genus_H_DI_results<-read.table("/local1/workdir1/tw488/CIDA/ITS_Maaslin2_Output/HDI_genus_ITS_output_11162024/all_results.tsv",header=TRUE)

#select columns from diff abundance data, calc log 2 fold change
ITS_genus_H_DI_results0.05<-ITS_genus_H_DI_results%>%filter(qval<0.05)%>%select(feature,value,coef,qval)%>% mutate(log2_fold_change = log2(exp(coef)))%>%select(feature,value,qval,log2_fold_change)

#find average relative abundance by day and merge with maaslin2 output
HDI_byDAY0.05vector<-c(ITS_genus_H_DI_results0.05$feature)
genus_rel_abundHDI_noFL_taxaMEANbyDAY<-genus_rel_abundHDI_noFL%>%group_by(taxon,Day) %>% summarise(AvgbyDay= mean(rel_abund))
genus_rel_abundHDI_noFL_taxaMEANbyDAY1<-genus_rel_abundHDI_noFL_taxaMEANbyDAY%>%filter(taxon %in% HDI_byDAY0.05vector)
genus_rel_abundHDI_noFL_taxaMEANbyDAY2<-dcast(genus_rel_abundHDI_noFL_taxaMEANbyDAY1,taxon~Day,value.var="AvgbyDay")
colnames(genus_rel_abundHDI_noFL_taxaMEANbyDAY2)[1] <- "Genus"
colnames(genus_rel_abundHDI_noFL_taxaMEANbyDAY2)[2] <- "Day Initial Relative Abundance (%)"
colnames(genus_rel_abundHDI_noFL_taxaMEANbyDAY2)[3] <- "Harvest Relative Abundance (%)"
colnames(ITS_genus_H_DI_results0.05)[1] <- "Genus"
colnames(ITS_genus_H_DI_results0.05)[3] <- "FDR"
ITS_genus_H_DI_results0.05_4export<-merge(genus_rel_abundHDI_noFL_taxaMEANbyDAY2,ITS_genus_H_DI_results0.05, by="Genus")
ITS_genus_H_DI_results0.05_4export<-select(ITS_genus_H_DI_results0.05_4export,1,3,2,6,5)
#write.csv(ITS_genus_H_DI_results0.05_4export,"/local1/workdir1/tw488/CIDA/ITS_Maaslin2_Output/HDI_genus_ITS_output_11162024/ITS_Genus_H_DI", row.names = TRUE)

#____SHELF LIFE SAMPLES_________________________________________________________________________

#NMDS for FIVEDAY shelf life samples

fivedayotu<-pstoveg_otu(fiveday)
fivedaysd<-pstoveg_sd(fiveday)
fivedaysd$Monthlot = paste(fivedaysd$Month,fivedaysd$Lot, sep="_")
fivedayvegdist<-vegdist(fivedayotu, method="bray")

fiveday_NMDS <- metaMDS(fivedayvegdist)
fiveday_NMDS
stressplot(fiveday_NMDS)

fivedaysd<-fivedaysd %>% rename(SampleName = sample)
fivedaysdNMDS<-fivedaysd%>%rownames_to_column("sample")

fivedayplot_df <- scores(fiveday_NMDS, display = "sites") %>% 
  as.data.frame() %>% 
  rownames_to_column("sample") %>% 
  full_join(fivedaysdNMDS, by = "sample")
fivedaypalNMDS<- c("chocolate", "#B2DF8A","red","#FF7F00", "#CAB2D6","black","darkblue","grey39","blue")

fivedayplot_df$Day<-factor(fivedayplot_df$Day,
                           levels=c("DI","D7","D12","D17","D22"))

fivedayplot_df<-fivedayplot_df%>%filter(SampleName!="ITS_1222_A_D7_3_AZ")

plot_fiveday_nmds <- ggplot(fivedayplot_df, aes(x = NMDS1, y = NMDS2, color = Day, shape = Location)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = fivedaypalNMDS) +
  stat_ellipse(aes(color=Location), linetype=2, level=0.5)+
  theme_classic()
plot_fiveday_nmds
ggsave("FIVEDAY_FUNGAL_NMDS_04292025.tiff", units="in", width=8, height=6, dpi=300, compression = 'lzw')

#Maaslin2  to find diff abundance over shelf life____________________________________________
#make correct data tables, filter out FL and GA
DI_D22_D28_noFLGA<- subset_samples(DI_D22_D28,Location!= "FL")
DI_D22_D28_noFLGA<- subset_samples(DI_D22_D28_noFLGA,Location!= "GA")
DI_D22_D28_noFLGA<- subset_samples(DI_D22_D28_noFLGA,sample!="ITS_1222_A_D7_3_AZ_contig.gz")

#make sample data and taxa frames
allShelfLife_4maaslin_TAX<-data.frame(DI_D22_D28_noFLGA@tax_table)
allShelfLife_4maaslin_sampledata<-data.frame(DI_D22_D28_noFLGA@sam_data)
allShelfLife_4maaslin_sampledata$Monthlot = paste(allShelfLife_4maaslin_sampledata$Month, allShelfLife_4maaslin_sampledata$Lot, sep="_")
#make a separate column based off of Actual_Day to make Day a continuous variable
allShelfLife_4maaslin_sampledata<-allShelfLife_4maaslin_sampledata%>%mutate(Actual_Day_noD=Actual_Day)%>%filter(Rep!="NC")
allShelfLife_4maaslin_sampledata$Actual_Day_noD<-gsub("DI", "1", allShelfLife_4maaslin_sampledata$Actual_Day_noD)
allShelfLife_4maaslin_sampledata$Actual_Day_noD<-gsub("D", "", allShelfLife_4maaslin_sampledata$Actual_Day_noD)
allShelfLife_4maaslin_sampledata$Actual_Day_noD<-as.integer(allShelfLife_4maaslin_sampledata$Actual_Day_noD)

#check to make sure 'Actual-Day_noD' is an integer now
str(allShelfLife_4maaslin_sampledata)

#make genus level otu table
allShelfLife.counts <- count_seqs(DI_D22_D28_noFLGA@otu_table)

abund.allShelfLife <- DI_D22_D28_noFLGA@otu_table %>% as.matrix() %>% as.data.frame() %>%
  rownames_to_column("OTU_") %>%
  pivot_longer(-OTU_, names_to = "sample") %>%
  merge(allShelfLife.counts$counts, by="sample") %>%
  dplyr::select(-n) %>%
  mutate(newcol=sample)%>%
  separate(newcol, c("X16S","Month", "Lot","Day","Rep","Loc","Actual_Day"), "_",remove = TRUE)%>%
  group_by(sample) %>%
  mutate(rel_abund = value/counts) %>% #relative abundance by treatment/day/rep
  ungroup() 

#check that relative abundance is calculated correctly
abund.allShelfLife %>% group_by(sample) %>% summarise(sum = sum(rel_abund))

#get taxonomy
allShelfLife_tax <- DI_D22_D28@tax_table %>%
  as.data.frame() %>%
  rownames_to_column("OTU_") 
head(allShelfLife_tax)

#merge taxonomy with rarefied relative abundance
count_tax_allShelfLife<- merge(abund.allShelfLife, allShelfLife_tax, by="OTU_")
head(count_tax_allShelfLife)

#### Relative abundance plot ####

# Get summary abundance of each taxa
taxa.meta.sum.allShelfLife <- count_tax_allShelfLife%>%
  dplyr::select(-c(value, counts)) %>%
  pivot_longer(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
               names_to = "level", 
               values_to = "taxon") 

#format phylum names
genus_rel_abund_allShelfLife <- taxa.meta.sum.allShelfLife %>%
  filter(level == "Genus") %>%
  group_by(taxon, sample,Month, Lot,Day,Rep,Loc,Actual_Day) %>% 
  summarise(rel_abund = 100*sum(rel_abund), .groups="drop") %>%
  mutate(taxon = str_replace(taxon,
                             "^unclassified (.*)", "Unclassified *\\1*"),
         taxon = str_replace(taxon, "^(\\S*)$", "\\1")) %>%
  ungroup()
genus_rel_abund_allShelfLife$taxon <-genus_rel_abund_allShelfLife$taxon%>%replace_na('Unclassified')
genus_rel_abund_allShelfLife_nodcast<-genus_rel_abund_allShelfLife
genus_rel_abund_allShelfLife_4maaslin<-dcast(genus_rel_abund_allShelfLife,sample~taxon,value.var="rel_abund")
genus_rel_abund_allShelfLife_4maaslin<-genus_rel_abund_allShelfLife_4maaslin%>% column_to_rownames(var="sample")

#how many genera are there total? 999999999999999999
totalHgenera_before5<-genus_rel_abund_allShelfLife_nodcast%>%group_by(taxon)%>%distinct(taxon)

#run maaslin2 for shelf life samples at genus level

AllShelfLife_genus_ITS_08162024 = Maaslin2(input_data =genus_rel_abund_allShelfLife_4maaslin, 
                                           input_metadata = allShelfLife_4maaslin_sampledata,
                                           analysis_method = "CPLM",
                                           transform= "NONE",
                                           min_prevalence = 0,
                                           min_abundance = 0,
                                           normalization  = "NONE",
                                           max_pngs = 200,
                                           save_models="TRUE",
                                           save_scatter="TRUE",
                                           cores=2,
                                           random_effects = c("Monthlot"),
                                           output         = "/local1/workdir1/tw488/CIDA/ITS_Maaslin2_Output/AllShelfLife_genus_ITS_08162024", 
                                           fixed_effects  = c("Actual_Day_noD","Location"),
                                           reference      = c("Location,CA")) 

#import results
AllShelfLife_genus_ITS_results<-read.table("/local1/workdir1/tw488/CIDA/ITS_Maaslin2_Output/AllShelfLife_genus_ITS_08162024/significant_results.tsv",header=TRUE)
AllShelfLife_genus_ITS_results0.05<-AllShelfLife_genus_ITS_results%>%filter(qval<0.05)%>%
  select("feature", "value", "coef","qval")%>% mutate(log2_fold_change = log2(exp(coef)))
colnames(AllShelfLife_genus_ITS_results0.05)[1] <- "Genus"
colnames(AllShelfLife_genus_ITS_results0.05)[4] <- "False Discovery Rate"

# Divy up in to day and location
#Take the time data, split time in to 1st half and second half of season by dividing month associated with a location
AllShelfLife_genus_ITS_results0.05_time<-AllShelfLife_genus_ITS_results0.05%>%filter(value!="AZ")
AllShelfLife_genus_ITS_results0.05_loc<-AllShelfLife_genus_ITS_results0.05%>%filter(value=="AZ")

#look at top average genera over shelf life
genus_rel_abund_allShelfLife_nodcast1<-genus_rel_abund_allShelfLife_nodcast%>%group_by(taxon)%>%mutate(avgrelab=mean(rel_abund))%>%distinct(taxon,avgrelab)
top20_genovershelflife<-head(arrange(genus_rel_abund_allShelfLife_nodcast1,desc(avgrelab)), n = 20)
top20_genovershelflife
top20_genovershelflifetaxon<-c(top20_genovershelflife$taxon)

#how many genera in all shelf life samples?
shelflifeonly<-genus_rel_abund_allShelfLife_nodcast1%>%filter(avgrelab!=0)

#get average DI and D22 values for each of the top 20 genera
#need to distinguish btwn 7 day tested samples and 5 day tested samples
genus_rel_abund_allShelfLife_nodcast2<-genus_rel_abund_allShelfLife_nodcast%>%filter(Day!="D7")%>%filter(Day!="D12")
genus_rel_abund_allShelfLife_nodcast2$Monthlot = paste(genus_rel_abund_allShelfLife_nodcast2$Month,genus_rel_abund_allShelfLife_nodcast2$Lot, sep="_")
#get average relative abundance over shelf life
genus_rel_abund_allShelfLife_nodcast2_overshelflife<-genus_rel_abund_allShelfLife_nodcast2%>%group_by(taxon) %>%
  mutate(avgbtwnreps=mean(rel_abund))%>%distinct(taxon,avgbtwnreps)
colnames(genus_rel_abund_allShelfLife_nodcast2_overshelflife)[1] <- "Genus"
#7 day; filter Day 22 b/c thats actually Day 28
genus_rel_abund_allShelfLife_nodcast2_7day<-genus_rel_abund_allShelfLife_nodcast2%>%
  filter(Month!="0122")%>%filter(Month!="0222")%>%filter(Month!="0322")%>%filter(Month!="0722")%>%
  filter(Monthlot!="0622_B")%>%filter(Month!="0822")%>%filter(Month!="0922")%>%
  filter(Month!="1022")%>%filter(Month!="1122")%>%filter(Month!="1222")%>%filter(Day!="D22")
genus_rel_abund_allShelfLife_nodcast2_7day$Day<-gsub("D17", "D22", genus_rel_abund_allShelfLife_nodcast2_7day$Day)
#5 day interval
genus_rel_abund_allShelfLife_nodcast2_5day<-genus_rel_abund_allShelfLife_nodcast2%>%
  filter(Month!="0522")%>%filter(Month!="0422")%>%filter(Monthlot!="0622_A")%>%filter(Day!="D17")
#merge the 5 day interval and 7 day interval back together and calculate relative abundance by day
genus_rel_abund_allShelfLife_nodcast2v2<-rbind(genus_rel_abund_allShelfLife_nodcast2_5day,genus_rel_abund_allShelfLife_nodcast2_7day)
genus_rel_abund_allShelfLife_nodcast2v2<-genus_rel_abund_allShelfLife_nodcast2v2%>%group_by(taxon,Day) %>%
  mutate(avgbtwnreps=mean(rel_abund))%>%distinct(taxon,Day,avgbtwnreps)
#dcast for the data frame,add in average relative abundance over shelf life
genus_rel_abund_allShelfLife_2v2<-dcast(genus_rel_abund_allShelfLife_nodcast2v2,taxon~Day,value.var="avgbtwnreps")
colnames(genus_rel_abund_allShelfLife_2v2)[1] <- "Genus"
genus_rel_abund_allShelfLife_2v2<-merge(genus_rel_abund_allShelfLife_2v2,genus_rel_abund_allShelfLife_nodcast2_overshelflife, by="Genus")
#filter for top 20 genera
genus_rel_abund_allShelfLife_2v2_filter<-genus_rel_abund_allShelfLife_2v2%>%filter(Genus %in% top20_genovershelflifetaxon)
AllShelfLife_genus_ITS_results0.05_time_filter<-AllShelfLife_genus_ITS_results0.05_time%>%filter(Genus %in% top20_genovershelflifetaxon)

#merge with maaslin2 output
allShelfLife_top10_time_df<-left_join(genus_rel_abund_allShelfLife_2v2_filter, AllShelfLife_genus_ITS_results0.05_time_filter, by="Genus")
allShelfLife_top10_time_df<-allShelfLife_top10_time_df%>%select(-5,-6)
colnames(allShelfLife_top10_time_df)[3] <- "Average Relative Abundance at Day Initial"
colnames(allShelfLife_top10_time_df)[2] <- "Average Relative Abundance at Day 22"
colnames(allShelfLife_top10_time_df)[4] <- "Average Relative Abundance Over Shelf Life"
#write.csv(allShelfLife_top10_time_df,"/local1/workdir1/tw488/CIDA/ITS_Maaslin2_Output/AllShelfLife_genus_ITS_08162024/allShelfLife_top20_time_df", row.names = TRUE)

#Results significant by Area
AllShelfLife_genus_ITS_results0.05_loc
AllShelfLife_genus_ITS_results0.05_loc_df<-AllShelfLife_genus_ITS_results0.05_loc%>%select(-2,-3)

#calculate shelf life relative abundance by area
#look at top avg genera over shelf life
genus_rel_abund_allShelfLife_nodcastA<-genus_rel_abund_allShelfLife_nodcast%>%group_by(taxon,Loc)%>%mutate(avgrelab=mean(rel_abund))%>%distinct(taxon,Loc,avgrelab)

#dcast for the data frame,add in average relative abundance over shelf life
genus_rel_abund_allShelfLife_dcastA<-dcast(genus_rel_abund_allShelfLife_nodcastA,taxon~Loc,value.var="avgrelab")
colnames(genus_rel_abund_allShelfLife_dcastA)[1] <- "Genus"

#make exportable list
genus_rel_abund_allShelfLife_dcastA<-merge(genus_rel_abund_allShelfLife_dcastA,AllShelfLife_genus_ITS_results0.05_loc_df, by="Genus", all.y=TRUE)
colnames(genus_rel_abund_allShelfLife_dcastA)[2] <- "Yuma, AZ area"
colnames(genus_rel_abund_allShelfLife_dcastA)[3] <- "Salinas, CA area"
#write.csv(genus_rel_abund_allShelfLife_dcastA,"/local1/workdir1/tw488/CIDA/ITS_Maaslin2_Output/AllShelfLife_genus_ITS_08162024/allShelfLife_byAREA", row.names = TRUE)

#find log fold change of genera that were not significantly different: Cladosporium,Stemphylium,Papiliotrema,Cystofilobasidium
AllShelfLife_genus_ITS_allresults<-read.table("/local1/workdir1/tw488/CIDA/ITS_Maaslin2_Output/AllShelfLife_genus_ITS_08162024/all_results.tsv",header=TRUE)
top20genera_notsigdiff<-c("Cladosporium", "Stemphylium", "Papiliotrema", "Cystofilobasidium")

#filter for top 20 genera
AllShelfLife_genus_ITS_results_genera_notsigdiff<-AllShelfLife_genus_ITS_allresults%>%filter(feature %in% top20genera_notsigdiff)

#filter for 
AllShelfLife_genus_ITS_results_genera_notsigdiff<-AllShelfLife_genus_ITS_results_genera_notsigdiff%>%filter(value=="Actual_Day_noD")%>%
  select("feature", "value", "coef","qval")%>% mutate(log2_fold_change = log2(exp(coef)))

colnames(AllShelfLife_genus_ITS_results_genera_notsigdiff)[1] <- "Genus"
colnames(AllShelfLife_genus_ITS_results_genera_notsigdiff)[4] <- "False Discovery Rate"


#Maaslin2 with OTUs_____________________________________________________________________

allShelfLife_4maaslin_OTUtbl<-as.data.frame(DI_D22_D28_noFLGA@otu_table)

AllShelfLife_OTU_ITS_08192024 = Maaslin2(input_data =DI_D22_D28_noFLGA@otu_table, 
                                         input_metadata = allShelfLife_4maaslin_sampledata,
                                         analysis_method = "CPLM",
                                         transform= "NONE",
                                         min_prevalence = 0,
                                         min_abundance = 0,
                                         normalization  = "NONE",
                                         max_pngs = 200,
                                         save_models="TRUE",
                                         save_scatter="TRUE",
                                         cores=2,
                                         random_effects = c("Monthlot"),
                                         output         = "/local1/workdir1/tw488/CIDA/ITS_Maaslin2_Output/AllShelfLife_OTU_ITS_08192024", 
                                         fixed_effects  = c("Actual_Day_noD","Location"),
                                         reference      = c("Location,CA")) 

#import results
AllShelfLife_OTU_ITS_results<-read.table("/local1/workdir1/tw488/CIDA/ITS_Maaslin2_Output/AllShelfLife_OTU_ITS_08192024/significant_results.tsv",header=TRUE)
AllShelfLife_OTU_ITS_results0.05<-AllShelfLife_OTU_ITS_results%>%filter(qval<0.05)%>%
  select("feature", "value", "coef","qval")%>% mutate(log2_fold_change = log2(exp(coef)))%>%select("feature", "value","qval","log2_fold_change")
colnames(AllShelfLife_OTU_ITS_results0.05)[1] <- "OTU"
colnames(AllShelfLife_OTU_ITS_results0.05)[3] <- "False Discovery Rate"

#separate Day and dfsos results
AllShelfLife_OTU_ITS_results0.05_time<-AllShelfLife_OTU_ITS_results0.05%>%filter(value!="AZ")
AllShelfLife_OTU_ITS_results0.05_area<-AllShelfLife_OTU_ITS_results0.05%>%filter(value=="AZ")

# match OTUs to taxonomy
allShelfLife_4maaslin_TAX_edit<-allShelfLife_4maaslin_TAX%>%select(-1,-2,-3,-4,-5,-7)%>%rownames_to_column()
colnames(allShelfLife_4maaslin_TAX_edit)[1] <- "OTU"

AllShelfLife_OTU_ITS_results0.05_time<-right_join(allShelfLife_4maaslin_TAX_edit,AllShelfLife_OTU_ITS_results0.05_time, by="OTU")
AllShelfLife_OTU_ITS_results0.05_area<-right_join(allShelfLife_4maaslin_TAX_edit,AllShelfLife_OTU_ITS_results0.05_area, by="OTU")

#get total OTUs for TIME based 
AllShelfLife_OTU_ITS_results0.05_time_counted_incr<-AllShelfLife_OTU_ITS_results0.05_time%>%filter(log2_fold_change>0)%>%group_by(Genus)%>% dplyr::count(Genus)
colnames(AllShelfLife_OTU_ITS_results0.05_time_counted_incr)[2] <- "Number of OTUs of Increasing Differential Abundance"
AllShelfLife_OTU_ITS_results0.05_time_counted_dec<-AllShelfLife_OTU_ITS_results0.05_time%>%filter(log2_fold_change<0)%>%group_by(Genus)%>%dplyr::count(Genus)
names(AllShelfLife_OTU_ITS_results0.05_time_counted_dec)[2]<-paste("Number of OTUs of Decreasing Differential Abundance")
AllShelfLife_OTU_ITS_results0.05_time_counted_incr_decr<-merge(AllShelfLife_OTU_ITS_results0.05_time_counted_incr,AllShelfLife_OTU_ITS_results0.05_time_counted_dec,by="Genus",all=TRUE)

#get total OTUs for AREA based
AllShelfLife_OTU_ITS_results0.05_area_counted_AZ<-AllShelfLife_OTU_ITS_results0.05_area%>%filter(log2_fold_change>0)%>%group_by(Genus)%>% dplyr::count(Genus)
colnames(AllShelfLife_OTU_ITS_results0.05_area_counted_AZ)[2] <- "Number of OTUs Differentially Abundant in AZ"
AllShelfLife_OTU_ITS_results0.05_area_counted_CA<-AllShelfLife_OTU_ITS_results0.05_area%>%filter(log2_fold_change<0)%>%group_by(Genus)%>% dplyr::count(Genus)
names(AllShelfLife_OTU_ITS_results0.05_area_counted_CA)[2]<-paste("Number of OTUs Differentially Abundant in CA")
AllShelfLife_OTU_ITS_results0.05_area_counted_CA_AZ<-merge(AllShelfLife_OTU_ITS_results0.05_area_counted_CA,AllShelfLife_OTU_ITS_results0.05_area_counted_AZ,by="Genus",all=TRUE)

#TOTAL ASVS PER GENUS
allShelfLife_4maaslin_TAX_edit$Genus<-allShelfLife_4maaslin_TAX_edit$Genus %>% replace_na('Unclassified')
allShelfLife_4maaslin_TAX_edit_taxasum<-allShelfLife_4maaslin_TAX_edit%>% group_by(Genus) %>% dplyr::count(Genus)
names(AllShelfLife_OTU_ITS_results0.05_area_counted_CA)[2]<-paste("Number of OTUs Differentially Abundant in CA")

#AREA based genera
area<-c(AllShelfLife_OTU_ITS_results0.05_area$Genus)
allShelfLife_4maaslin_TAX_edit_taxasum_area<-allShelfLife_4maaslin_TAX_edit_taxasum%>%filter(Genus %in% area)
names(allShelfLife_4maaslin_TAX_edit_taxasum_area)[2]<-paste("Total Number of OTUs in Area")

#time based genera
time<-c(AllShelfLife_OTU_ITS_results0.05_time$Genus)
allShelfLife_4maaslin_TAX_edit_taxasum_time<-allShelfLife_4maaslin_TAX_edit_taxasum%>%filter(Genus %in% time)
names(allShelfLife_4maaslin_TAX_edit_taxasum_time)[2]<-paste("Total Number of OTUs in Time ")

#merge in to final df
#time
allShelfLife_4maaslin_TAX_edit_taxasum_time_final <- merge(allShelfLife_4maaslin_TAX_edit_taxasum_time,AllShelfLife_OTU_ITS_results0.05_time_counted_incr_decr, by="Genus", all=TRUE)
allShelfLife_4maaslin_TAX_edit_taxasum_time_final[is.na(allShelfLife_4maaslin_TAX_edit_taxasum_time_final)] <- 0
#write.csv(allShelfLife_4maaslin_TAX_edit_taxasum_time_final,"/local1/workdir1/tw488/CIDA/ITS_Maaslin2_Output/AllShelfLife_OTU_ITS_08192024/ITS_OTU_shelflife_time_Table", row.names = TRUE,na = "0")

#area
allShelfLife_4maaslin_TAX_edit_taxasum_area_final <- merge(allShelfLife_4maaslin_TAX_edit_taxasum_area,AllShelfLife_OTU_ITS_results0.05_area_counted_CA_AZ, by="Genus", all=TRUE)
allShelfLife_4maaslin_TAX_edit_taxasum_area_final[is.na(allShelfLife_4maaslin_TAX_edit_taxasum_area_final)] <- 0
#write.csv(allShelfLife_4maaslin_TAX_edit_taxasum_area_final,"/local1/workdir1/tw488/CIDA/ITS_Maaslin2_Output/AllShelfLife_OTU_ITS_08192024/ITS_OTU_shelflife_area_Table", row.names = TRUE,na = "0")

#combo
combo<-merge(allShelfLife_4maaslin_TAX_edit_taxasum_area_final,allShelfLife_4maaslin_TAX_edit_taxasum_time_final,by="Genus",all=TRUE)
#write.csv(combo,"/local1/workdir1/tw488/CIDA/ITS_Maaslin2_Output/AllShelfLife_OTU_ITS_08192024/ITS_OTU_shelflife_combo_Table", row.names = TRUE)

#TOP OTUs________________________________________________________________________________

allShelfLife_4maaslin_OTUtbl_rowsums<-as.data.frame(rowSums(allShelfLife_4maaslin_OTUtbl))
allShelfLife_4maaslin_OTUtbl_rowsums_top<-head(arrange(allShelfLife_4maaslin_OTUtbl_rowsums,desc(rowSums(allShelfLife_4maaslin_OTUtbl))), n = 3)
allShelfLife_4maaslin_OTUtbl_rowsums_top<-merge(allShelfLife_4maaslin_OTUtbl_rowsums_top,allShelfLife_4maaslin_TAX, by="row.names")

#________#PERMANOVA of 7 Day INTERVAL tested samples__________________________________________________________________

#WITH BRAY CURTIS
sevendaySHELFLIFEotu<-pstoveg_otu(sevenday)
sevendaySHELFLIFEsd<-pstoveg_sd(sevenday)
sevendaySHELFLIFEsd$Monthlot = paste(sevendaySHELFLIFEsd$Month, sevendaySHELFLIFEsd$Lot, sep="_")
sevendaySHELFLIFEvegdist<-vegdist(sevendaySHELFLIFEotu, method="bray")

#permanova
sevendaysd_otu<-cbind(sevendaySHELFLIFEsd,sevendaySHELFLIFEotu)
cbind(rownames(sevendaySHELFLIFEsd),  rownames (sevendaySHELFLIFEotu))
sevendaysd_otu<-sevendaysd_otu%>%mutate(MonthlotDay=paste(Monthlot,Day,sep="_"),.after=Monthlot)
sevendaysd_otu<-sevendaysd_otu%>%arrange(MonthlotDay)

#check to make sure there are even # samples
table(sevendaysd_otu$Monthlot)
table(sevendaysd_otu$MonthlotDay)
# Must subsample plots to be balanced, but do this 100 times because of random sampling
hsubsevenday <- with(sevendaysd_otu, how(within=Within(type="free"),blocks=Monthlot,plots=Plots(strata=MonthlotDay, type="free")))
hsubsevenday

sevendaysd_permanova<-sevendaysd_otu%>%select(SampleType,quant.reading,is.neg,sample,Month,Lot,Day,Rep,Location,Actual_Day,Monthlot,MonthlotDay)
table(sevendaysd_permanova$Monthlot)
table(sevendaysd_permanova$Actual_Day)
table(sevendaysd_permanova$Monthlot,sevendaysd_permanova$Actual_Day)

sevendayotu_permanova<-sevendaysd_otu%>%ungroup()%>%select(-SampleType,-quant.reading,-is.neg,-sample,-Month,-Lot,-Day,-Rep,-Location,-Actual_Day,-Monthlot,-MonthlotDay, -R2)

#ensure replicates are permuted together
sevendaytestsd<-sevendaysd_permanova
sevendaytestsd$permute1<-shuffle(48,control=hsubsevenday)
sevendaytestsd$origorder<-1:48

#visualize permutations
sevendaytestsd%>%pivot_longer(permute1:origorder, names_to="column", values_to="position", cols_vary = "slowest") %>%
  ggplot(aes(x=column, y=position, group=sample, color=Day))+geom_line()

#make in to distance object
sevendayvegdist4<-vegdist(sevendayotu_permanova, method="bray")

#plug in here
sevenday_perm_trial <- adonis2(sevendayvegdist4~Actual_Day, data=sevendaysd_permanova, permutations = hsubsevenday, method = "bray", 
                               sqrt.dist = FALSE, add = FALSE, by = "terms",
                               parallel = getOption("mc.cores"), na.action = na.fail)
sevenday_perm_trial
#not significant!

#Day, location, and the interaction of
sevenday_perm_trial_dayloc <- adonis2(sevendayvegdist4~Actual_Day*Location, data=sevendaysd_permanova, permutations = hsubsevenday, method = "bray", 
                               sqrt.dist = FALSE, add = FALSE, by = "terms",
                               parallel = getOption("mc.cores"), na.action = na.fail)
sevenday_perm_trial_dayloc

#alpha diversity OF 7 DAY INTERVAL TESTED SAMPLES_____________________________________________________-
sevendayalphadiv<-estimate_richness(sevenday, split = TRUE, measures = NULL)
sevendayalphadiv$sample <- row.names(sevendayalphadiv)  
sevendayalphadiv <- sevendayalphadiv %>% separate(sample, c("x16s", "month", "lot", "day", "rep", "loc","r1"), "_")
sevendayalphadiv$sample <- row.names(sevendayalphadiv) 
sevendayalphadiv$sampling_code = paste(sevendayalphadiv$month, sevendayalphadiv$lot, sep="_")
sevendayalphadiv <- sevendayalphadiv %>% 
  as.data.frame() %>% 
  full_join(DFSOSv2, by = "sampling_code")%>%na.omit

#calculate Pielou
sevendayalphadiv <-sevendayalphadiv %>% mutate (Pielou=(sevendayalphadiv$Shannon)/log(sevendayalphadiv$Observed))

#Observed
sevendayalphadiv  %>% ggplot(aes(x=day, y=Observed)) +
  geom_jitter(size=0.25) +
  labs(x="Day",
       y="Number of Observed ASVs") +
  guides(color = guide_legend(override.aes = list(size=1))) +
  theme_classic()

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

sevendayalphadiv<-sevendayalphadiv%>%mutate(Actual_Day_noD=day)
sevendayalphadiv$Actual_Day_noD<-gsub("DI", "1", sevendayalphadiv$Actual_Day_noD)
sevendayalphadiv$Actual_Day_noD<-gsub("D", "", sevendayalphadiv$Actual_Day_noD)
sevendayalphadiv$Actual_Day_noD<-as.integer(sevendayalphadiv$Actual_Day_noD)
sevendayalphadiv$day<-as.factor(sevendayalphadiv$day)
str(sevendayalphadiv$day)
str(sevendayalphadiv$Shannon)

#ANOVA to compare days
#Shannon
sevendayalphadiv$day1<-as.character(sevendayalphadiv$day)
str(sevendayalphadiv$day1)

#order Days as factors
sevendayalphadiv$day1<-factor(sevendayalphadiv$day1,
                            levels=c("DI","D7","D12","D17","D22"))

ShnANOVA_lm<- aov(Shannon ~ day1, data =sevendayalphadiv)
em_ShnANOVA_lm<-emmeans(ShnANOVA_lm,~day1)
em_ShnANOVA_lm

pairs.em_ShnANOVA_lm<-as.data.frame(pairs(em_ShnANOVA_lm))
print(pairs.em_ShnANOVA_lm)
cldList(p.value~contrast,data=pairs.em_ShnANOVA_lm, remove.space="TRUE")

#Richness
ObservedANOVA_lm<- aov(Observed ~ day1, data =sevendayalphadiv)
em_ObservedANOVA_lm<-emmeans(ObservedANOVA_lm,~day1)
em_ObservedANOVA_lm

pairs.em_ObservedANOVA_lm<-as.data.frame(pairs(em_ObservedANOVA_lm))
print(pairs.em_ObservedANOVA_lm)
cldList(p.value~contrast,data=pairs.em_ObservedANOVA_lm, remove.space="TRUE")

#Evenness
PieANOVA<- aov(Pielou ~ day1, data =sevendayalphadiv)
em_PieANOVA<-emmeans(PieANOVA,~day1)
em_PieANOVA

pairs.em_PieANOVA<-as.data.frame(pairs(em_PieANOVA))
print(pairs.em_PieANOVA)
cldList(p.value~contrast,data=pairs.em_PieANOVA, remove.space="TRUE")

#calc average values for sevendayalphadiv
sevendayalphadiv <-sevendayalphadiv %>% group_by(day)%>%mutate(avgPielou=mean(Pielou))%>%mutate(avgShn=mean(Shannon))%>%mutate(avgObs=mean(Observed))
sevendayalphadiv_avgs <-sevendayalphadiv %>%distinct(day,avgPielou,avgShn,avgObs)
#write.csv(sevendayalphadiv_avgs,"/local1/workdir1/tw488/CIDA/For_Publication/ITS/ShelfLifeANOVA/Allshelflife_D7_alphadivTable_avgsforeachalphadivindex", row.names = TRUE)

#________#PERMANOVA of 5 Day Interval tested samples__________________________________________________________________-

#WITH BRAY CURTIS
fivedaySHELFLIFEotu<-pstoveg_otu(fiveday)
fivedaySHELFLIFEsd<-pstoveg_sd(fiveday)
fivedaySHELFLIFEsd$Monthlot = paste(fivedaySHELFLIFEsd$Month, fivedaySHELFLIFEsd$Lot, sep="_")
fivedaySHELFLIFEvegdist<-vegdist(fivedaySHELFLIFEotu, method="bray")

table(fivedaySHELFLIFEsd$Monthlot)
table(fivedaySHELFLIFEsd$Actual_Day)

#________PERMANOVA FOR select FIVE DAY SAMPLES, ______________________________________________-

fivedaySHELFLIFEsd_otu<-cbind(fivedaySHELFLIFEsd,fivedaySHELFLIFEotu)
fivedaySHELFLIFEsd_otu<-fivedaySHELFLIFEsd_otu%>%mutate(MonthlotDay=paste(Monthlot,Day,sep="_"),.after=Monthlot)

fivedaySHELFLIFEsd_otu%>%group_by(MonthlotDay)%>%dplyr::count()%>%View()
fivedaySHELFLIFEsd_otu<-fivedaySHELFLIFEsd_otu%>%filter(sample!="ITS_1222_A_D7_3_AZ")
fivedaySHELFLIFEsd_otu<-fivedaySHELFLIFEsd_otu%>%arrange(MonthlotDay)

#how function from permute package

fivedayhsubSHELFLIFE <- with(fivedaySHELFLIFEsd_otu, how(within=Within(type="free"),blocks=Monthlot,plots=Plots(strata=MonthlotDay, type="free")))
fivedayhsubSHELFLIFE

fivedaySHELFLIFEsd_permanova<-fivedaySHELFLIFEsd_otu%>%select(SampleType,quant.reading,is.neg,sample,Month,Lot,Day,Rep,Location,Actual_Day,Monthlot,MonthlotDay)
table(fivedaySHELFLIFEsd_permanova$Monthlot)
table(table(fivedaySHELFLIFEsd_permanova$MonthlotDay))
table(fivedaySHELFLIFEsd_permanova$Actual_Day)

fivedaySHELFLIFEotu_permanova<-fivedaySHELFLIFEsd_otu%>%ungroup()%>%select(-SampleType,-quant.reading,-is.neg,-sample,-Month,-Lot,-Day,-Rep,-Location,-Actual_Day,-Monthlot,-MonthlotDay,-R2)


#ensure replicates are permuted together
fivedaytestsd<-fivedaySHELFLIFEsd_permanova
fivedaytestsd$permute1<-shuffle(178,control=fivedayhsubSHELFLIFE)
fivedaytestsd$origorder<-1:178

#visualize permutations
fivedaytestsd%>%pivot_longer(permute1:origorder, names_to="column", values_to="position", cols_vary = "slowest") %>%
  ggplot(aes(x=column, y=position, group=sample, color=Day))+geom_line()

#make in to distance object
fivedaySHELFLIFEvegdist4<-vegdist(fivedaySHELFLIFEotu_permanova, method="bray")

#plug in here
fivedaySHELFLIFE_perm_trial <- adonis2(fivedaySHELFLIFEvegdist4~Actual_Day, data=fivedaySHELFLIFEsd_permanova, permutations = fivedayhsubSHELFLIFE, method = "bray", 
                                sqrt.dist = FALSE, add = FALSE, by = "terms",
                                parallel = getOption("mc.cores"), na.action = na.fail)
fivedaySHELFLIFE_perm_trial

#plug in here
fivedaySHELFLIFE_perm_dayloc <- adonis2(fivedaySHELFLIFEvegdist4~Actual_Day*Location, data=fivedaySHELFLIFEsd_permanova, permutations = fivedayhsubSHELFLIFE, method = "bray", 
                                       sqrt.dist = FALSE, add = FALSE, by = "terms",
                                       parallel = getOption("mc.cores"), na.action = na.fail)
fivedaySHELFLIFE_perm_dayloc

# pairwise adonis
#cant do pairwise adonis w/o pesudoreplication of the two reps, so doing it one by one
#DID7--significant
fivedaySHELFLIFE_perm_DID7 <- adonis2(vegdist(fivedaySHELFLIFEotu_permanova[fivedaySHELFLIFEsd_permanova$Actual_Day%in%c("DI","D7"),], method="bray")~Actual_Day, 
                               data=fivedaySHELFLIFEsd_permanova[fivedaySHELFLIFEsd_permanova$Actual_Day%in%c("DI","D7"),], 
                               permutations = with(fivedaySHELFLIFEsd_otu[fivedaySHELFLIFEsd_permanova$Actual_Day%in%c("DI","D7"),], how(within=Within(type="free"),blocks=Monthlot,plots=Plots(strata=MonthlotDay, type="free"))), 
                               method = "bray", 
                               sqrt.dist = FALSE, add = FALSE, by = "terms",
                               parallel = getOption("mc.cores"), na.action = na.fail)
fivedaySHELFLIFE_perm_DID7

#DID12--significant
fivedaySHELFLIFE_perm_DID12 <- adonis2(vegdist(fivedaySHELFLIFEotu_permanova[fivedaySHELFLIFEsd_permanova$Actual_Day%in%c("DI","D12"),], method="bray")~Actual_Day, 
                                data=fivedaySHELFLIFEsd_permanova[fivedaySHELFLIFEsd_permanova$Actual_Day%in%c("DI","D12"),], 
                                permutations = with(fivedaySHELFLIFEsd_otu[fivedaySHELFLIFEsd_permanova$Actual_Day%in%c("DI","D12"),], how(within=Within(type="free"),blocks=Monthlot,plots=Plots(strata=MonthlotDay, type="free"))), 
                                method = "bray", 
                                sqrt.dist = FALSE, add = FALSE, by = "terms",
                                parallel = getOption("mc.cores"), na.action = na.fail)
fivedaySHELFLIFE_perm_DID12

#DID17--significant
fivedaySHELFLIFE_perm_DID17 <- adonis2(vegdist(fivedaySHELFLIFEotu_permanova[fivedaySHELFLIFEsd_permanova$Actual_Day%in%c("DI","D17"),], method="bray")~Actual_Day, 
                                data=fivedaySHELFLIFEsd_permanova[fivedaySHELFLIFEsd_permanova$Actual_Day%in%c("DI","D17"),], 
                                permutations = with(fivedaySHELFLIFEsd_otu[fivedaySHELFLIFEsd_permanova$Actual_Day%in%c("DI","D17"),], how(within=Within(type="free"),blocks=Monthlot,plots=Plots(strata=MonthlotDay, type="free"))), 
                                method = "bray", 
                                sqrt.dist = FALSE, add = FALSE, by = "terms",
                                parallel = getOption("mc.cores"), na.action = na.fail)
fivedaySHELFLIFE_perm_DID17

#DID22-significant
fivedaySHELFLIFE_perm_DID22 <- adonis2(vegdist(fivedaySHELFLIFEotu_permanova[fivedaySHELFLIFEsd_permanova$Actual_Day%in%c("DI","D22"),], method="bray")~Actual_Day, 
                                data=fivedaySHELFLIFEsd_permanova[fivedaySHELFLIFEsd_permanova$Actual_Day%in%c("DI","D22"),], 
                                permutations = with(fivedaySHELFLIFEsd_otu[fivedaySHELFLIFEsd_permanova$Actual_Day%in%c("DI","D22"),], how(within=Within(type="free"),blocks=Monthlot,plots=Plots(strata=MonthlotDay, type="free"))), 
                                method = "bray", 
                                sqrt.dist = FALSE, add = FALSE, by = "terms",
                                parallel = getOption("mc.cores"), na.action = na.fail)
fivedaySHELFLIFE_perm_DID22

#D7D12-not significant

fivedaySHELFLIFE_perm_D7D12 <- adonis2(vegdist(fivedaySHELFLIFEotu_permanova[fivedaySHELFLIFEsd_permanova$Actual_Day%in%c("D7","D12"),], method="bray")~Actual_Day, 
                                data=fivedaySHELFLIFEsd_permanova[fivedaySHELFLIFEsd_permanova$Actual_Day%in%c("D7","D12"),], 
                                permutations = with(fivedaySHELFLIFEsd_otu[fivedaySHELFLIFEsd_permanova$Actual_Day%in%c("D7","D12"),], how(within=Within(type="free"),blocks=Monthlot,plots=Plots(strata=MonthlotDay, type="free"))), 
                                method = "bray", 
                                sqrt.dist = FALSE, add = FALSE, by = "terms",
                                parallel = getOption("mc.cores"), na.action = na.fail)
fivedaySHELFLIFE_perm_D7D12

#D7D17-not significant
fivedaySHELFLIFE_perm_D7D17 <- adonis2(vegdist(fivedaySHELFLIFEotu_permanova[fivedaySHELFLIFEsd_permanova$Actual_Day%in%c("D7","D17"),], method="bray")~Actual_Day, 
                                data=fivedaySHELFLIFEsd_permanova[fivedaySHELFLIFEsd_permanova$Actual_Day%in%c("D7","D17"),], 
                                permutations = with(fivedaySHELFLIFEsd_otu[fivedaySHELFLIFEsd_permanova$Actual_Day%in%c("D7","D17"),], how(within=Within(type="free"),blocks=Monthlot,plots=Plots(strata=MonthlotDay, type="free"))), 
                                method = "bray", 
                                sqrt.dist = FALSE, add = FALSE, by = "terms",
                                parallel = getOption("mc.cores"), na.action = na.fail)
fivedaySHELFLIFE_perm_D7D17

#D7D22- not significant
fivedaySHELFLIFE_perm_D7D22 <- adonis2(vegdist(fivedaySHELFLIFEotu_permanova[fivedaySHELFLIFEsd_permanova$Actual_Day%in%c("D7","D22"),], method="bray")~Actual_Day, 
                                data=fivedaySHELFLIFEsd_permanova[fivedaySHELFLIFEsd_permanova$Actual_Day%in%c("D7","D22"),], 
                                permutations = with(fivedaySHELFLIFEsd_otu[fivedaySHELFLIFEsd_permanova$Actual_Day%in%c("D7","D22"),], how(within=Within(type="free"),blocks=Monthlot,plots=Plots(strata=MonthlotDay, type="free"))), 
                                method = "bray", 
                                sqrt.dist = FALSE, add = FALSE, by = "terms",
                                parallel = getOption("mc.cores"), na.action = na.fail)
fivedaySHELFLIFE_perm_D7D22

#D12D17- not significant

fivedaySHELFLIFE_perm_D12D17 <- adonis2(vegdist(fivedaySHELFLIFEotu_permanova[fivedaySHELFLIFEsd_permanova$Actual_Day%in%c("D12","D17"),], method="bray")~Actual_Day, 
                                 data=fivedaySHELFLIFEsd_permanova[fivedaySHELFLIFEsd_permanova$Actual_Day%in%c("D12","D17"),], 
                                 permutations = with(fivedaySHELFLIFEsd_otu[fivedaySHELFLIFEsd_permanova$Actual_Day%in%c("D12","D17"),], how(within=Within(type="free"),blocks=Monthlot,plots=Plots(strata=MonthlotDay, type="free"))), 
                                 method = "bray", 
                                 sqrt.dist = FALSE, add = FALSE, by = "terms",
                                 parallel = getOption("mc.cores"), na.action = na.fail)
fivedaySHELFLIFE_perm_D12D17

#D12D22- not significant
fivedaySHELFLIFE_perm_D12D22 <- adonis2(vegdist(fivedaySHELFLIFEotu_permanova[fivedaySHELFLIFEsd_permanova$Actual_Day%in%c("D12","D22"),], method="bray")~Actual_Day, 
                                 data=fivedaySHELFLIFEsd_permanova[fivedaySHELFLIFEsd_permanova$Actual_Day%in%c("D12","D22"),], 
                                 permutations = with(fivedaySHELFLIFEsd_otu[fivedaySHELFLIFEsd_permanova$Actual_Day%in%c("D12","D22"),], how(within=Within(type="free"),blocks=Monthlot,plots=Plots(strata=MonthlotDay, type="free"))), 
                                 method = "bray", 
                                 sqrt.dist = FALSE, add = FALSE, by = "terms",
                                 parallel = getOption("mc.cores"), na.action = na.fail)
fivedaySHELFLIFE_perm_D12D22

#d17D22- not sig

fivedaySHELFLIFE_perm_D17D22 <- adonis2(vegdist(fivedaySHELFLIFEotu_permanova[fivedaySHELFLIFEsd_permanova$Actual_Day%in%c("D17","D22"),], method="bray")~Actual_Day, 
                                 data=fivedaySHELFLIFEsd_permanova[fivedaySHELFLIFEsd_permanova$Actual_Day%in%c("D17","D22"),], 
                                 permutations = with(fivedaySHELFLIFEsd_otu[fivedaySHELFLIFEsd_permanova$Actual_Day%in%c("D17","D22"),], how(within=Within(type="free"),blocks=Monthlot,plots=Plots(strata=MonthlotDay, type="free"))), 
                                 method = "bray", 
                                 sqrt.dist = FALSE, add = FALSE, by = "terms",
                                 parallel = getOption("mc.cores"), na.action = na.fail)
fivedaySHELFLIFE_perm_D17D22

#alpha diversity OF 5 DAY INTERVAL TESTED SAMPLES_____________________________________________________-

fivedayalphadiv<-estimate_richness(fiveday, split = TRUE, measures = NULL)
fivedayalphadiv$sample <- row.names(fivedayalphadiv)  
fivedayalphadiv <- fivedayalphadiv %>% separate(sample, c("x16s", "month", "lot", "day", "rep", "loc","r1"), "_")
fivedayalphadiv$sample <- row.names(fivedayalphadiv) 
fivedayalphadiv$sampling_code = paste(fivedayalphadiv$month,fivedayalphadiv$lot, sep="_")
fivedayalphadiv <- fivedayalphadiv %>% 
  as.data.frame() %>% 
  full_join(DFSOSv2, by = "sampling_code")%>%na.omit

fivedayalphadiv<-fivedayalphadiv%>%filter(sample!="ITS_1222_A_D7_3_AZ_contig.gz")

#calculate Pielou
fivedayalphadiv <-fivedayalphadiv%>%ungroup() %>% mutate (Pielou=(fivedayalphadiv$Shannon)/log(fivedayalphadiv$Observed))

#Observed
fivedayalphadiv  %>% ggplot(aes(x=day, y=Observed)) +
  geom_jitter(size=0.25) +
  #geom_smooth(aes(group=monthsample), se=FALSE, size=4) +
  labs(x="Day",
       y="Number of Observed ASVs") +
  #scale_x_continuous() +
  guides(color = guide_legend(override.aes = list(size=1))) +
  theme_classic()

#Pielou
fivedayalphadiv%>% ggplot(aes(x=day, y=Pielou)) +
  geom_jitter(size=0.25) +
  #geom_smooth(aes(group=monthsample), se=FALSE, size=4) +
  labs(x="Day",
       y="Pielou's Evenness") +
  #scale_x_continuous() +
  guides(color = guide_legend(override.aes = list(size=1))) +
  theme_classic()

#SHANNON
fivedayalphadiv%>% ggplot(aes(x=day, y=Shannon)) +
  geom_jitter(size=0.25) +
  #geom_smooth(aes(group=monthsample), se=FALSE, size=4) +
  labs(x="Day",
       y="Shannon Index") +
  #scale_x_continuous() +
  guides(color = guide_legend(override.aes = list(size=1))) +
  theme_classic()

#order Days as factors
fivedayalphadiv$day<-factor(fivedayalphadiv$day,
                                    levels=c("DI","D7","D12","D17","D22"))

#ANOVA to compare days
#Shannon
ShnANOVA5<- aov(Shannon ~ day, data =fivedayalphadiv)
em_ShnANOVA5<-emmeans(ShnANOVA5,~day)

pairs.em_ShnANOVA5<-as.data.frame(pairs(em_ShnANOVA5))
print(pairs.em_ShnANOVA5)
cldList(p.value~contrast,data=pairs.em_ShnANOVA5, remove.space="TRUE")

#Richness
ObsANOVA5<- aov(Observed~ day, data =fivedayalphadiv)
em_ObsANOVA5<-emmeans(ObsANOVA5,~day)

pairs.em_ObsANOVA5<-as.data.frame(pairs(em_ObsANOVA5))
print(pairs.em_ObsANOVA5)
cldList(p.value~contrast,data=pairs.em_ObsANOVA5, remove.space="TRUE")

#Evenness
PieANOVA5<- aov(Pielou ~ day, data =fivedayalphadiv)
em_PieANOVA5<-emmeans(PieANOVA5,~day)

pairs.em_PieANOVA5<-as.data.frame(pairs(em_PieANOVA5))
print(pairs.em_PieANOVA5)
cldList(p.value~contrast,data=pairs.em_PieANOVA5, remove.space="TRUE")

#calculate average values for 5dayalphadiv table in the manuscript
fivedayalphadiv <-fivedayalphadiv %>% group_by(day)%>%mutate(avgPielou=mean(Pielou))%>%mutate(avgShn=mean(Shannon))%>%mutate(avgObs=mean(Observed))
fivedayalphadiv_avgs <-fivedayalphadiv %>%distinct(day,avgPielou,avgShn,avgObs)
#write.csv(fivedayalphadiv_avgs,"/local1/workdir1/tw488/CIDA/For_Publication/ITS/ShelfLifeANOVA/Allshelflife_D5_alphadivTable_avgsforeachalphadivindex", row.names = TRUE)

#__TOP 5 GENERA FOR A TIMEPOINT BY AVG RELATIVE ABUNDANCE___________________________________________________________

#H and DI average relative abundance
genus_rel_abundHDI_noFL_taxaMEANbyDAY_b<-dcast(genus_rel_abundHDI_noFL_taxaMEANbyDAY,taxon~Day,value.var="AvgbyDay")

#H only
Honly<-genus_rel_abundHDI_noFL_taxaMEANbyDAY_b%>%select(-2)%>%arrange(desc(H))
Honly<-head(arrange(Honly), n = 5)
colnames(Honly)[2] <- "Average Relative Abundance at Harvest"

#DI only
DIonly<-genus_rel_abundHDI_noFL_taxaMEANbyDAY_b%>%select(-3)%>% arrange(desc(DI))
DIonly<-head(arrange(DIonly), n = 5)
colnames(DIonly)[2] <- "Average Relative Abundance at Day Initial"

#D7 average relative abundance
genus_rel_abund_allShelfLife_nodcast
D7toprelabund<-genus_rel_abund_allShelfLife_nodcast%>%filter(Day=="D7")%>%filter(sample!="ITS_1222_A_D7_3_AZ_contig.gz")%>%filter(Rep!="NC")%>%group_by(taxon)%>%mutate(avgrelab=mean(rel_abund))%>%distinct(taxon,avgrelab)
D7toprelabund_top5<-head(arrange(D7toprelabund,desc(avgrelab)), n = 5)
colnames(D7toprelabund_top5)[2] <- "Average Relative Abundance at Day 7"

#D21/22 average relative abundance
D21_22toprelabund<-genus_rel_abund_allShelfLife_2v2_filter%>%select(-3,-4)
D21_22toprelabund<-head(arrange(D21_22toprelabund,desc(D22)), n = 5)
colnames(D21_22toprelabund)[2] <- "Average Relative Abundance at Day 21/22"

#export all together now
write.csv(Honly,"/local1/workdir1/tw488/CIDA/For_Publication/ITS/Table10_top5avggenera_08192024/H", row.names = TRUE)
write.csv(DIonly,"/local1/workdir1/tw488/CIDA/For_Publication/ITS/Table10_top5avggenera_08192024/DayI", row.names = TRUE)
write.csv(D7toprelabund_top5,"/local1/workdir1/tw488/CIDA/For_Publication/ITS/Table10_top5avggenera_08192024/Day7", row.names = TRUE)
write.csv(D21_22toprelabund,"/local1/workdir1/tw488/CIDA/For_Publication/ITS/Table10_top5avggenera_08192024/Day22", row.names = TRUE)

