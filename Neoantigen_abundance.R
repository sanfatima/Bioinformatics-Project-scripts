#install packages 
#ggplot2 , vcfR , ggpubr , reshape2 , RColorBrewer
library(reshape)
library(ggplot2)
library(vcfR)
library (ggpubr)
library (reshape2)
library (RColorBrewer)
library(pgirmess)
library(plyr)
library(dplyr)
library(data.table)
# Install/load the tidyverse
if (require("tidyverse") == FALSE) {
  install.packages("tidyverse")
  library("tidyverse")
} else {
  library("tidyverse")
}


# Install/load ggthemr
if (require("ggthemr") == FALSE) {
  devtools::install_github('cttobin/ggthemr')
  library("ggthemr")
} else {
  library("ggthemr")
}

# Install/load pals
if (require("pals") == FALSE) {
  install.packages("pals")
  library("pals")
} else {
  library("pals")
}

# Install/load gridExtra
if (require("gridExtra") == FALSE) {
  install.packages("gridExtra")
  library("gridExtra")
} else {
  library("gridExtra")
}

# Install/load gridExtra
if (require("grid") == FALSE) {
  install.packages("grid")
  library("grid")
} else {
  library("grid")
}

# Install/load scales
if (require("scales") == FALSE) {
  install.packages("scales")
  library("scales")
} else {
  library("scales")
}

# Install/load extrafont
if (require("extrafont") == FALSE) {
  install.packages("extrafont")
  library("extrafont")
} else {
  library("extrafont")
}

library("Biostrings")

library("ggpubr")

library("cocor")



######################## SPCRC Filtration & Normalisation ##############################

###Neopred antigen table
#snv table
dir = '/Users/syedaanamfatima/Downloads/Neopredpipe/summary_table/spCRCs'
SPCRCneo_snv <- read.table(paste0(dir, '/Total.neoantigensj.txt'), header=F,
                           sep = '\t',stringsAsFactors = F, fill=T)

names(SPCRCneo_snv) <- c('Sample','R1','R2','R3','4','5','6','7','8','9','10','11','12','13','14','15','16', 'LineID', 'Chrom', 'allelepos',
                         'RefAll', 'AltAll', 'Genename:refseqID', 'pos', 'hla', 'peptide', 'core', 'Of', 'Gp',
                         'Gl', 'Ip', 'Il', 'Icore', 'Identity', 'Score', 'Affinity','Rank', 'Cand', 'BindLevel','Novelty')

#create antigen ID for each table based on Sample and Identity
SPCRCneo_snv$AntigenID<-apply(SPCRCneo_snv, 1,function(x) paste0(x['Sample'], ':',x['Identity']) )



#SPCRC recopo table for snvs
#snv table
SPCRCreco_snv <- read.table('/Users/syedaanamfatima/Downloads/Neopredpipe/summary_table/spCRCs/PredictedRecognitionPotentialsj.txt',
                            stringsAsFactors = F, header=T)


#create antigen ID for each table based on Sample and Mutation
SPCRCreco_snv$AntigenID<-apply(SPCRCreco_snv, 1,function(x) paste0(x['Sample'], ':',x['Mutation']))


#merge neoantigen and recopo tables for snvs together to allow filtration from both tables

SPCRCneoreco_snv<-full_join(SPCRCneo_snv,SPCRCreco_snv, by="AntigenID")


#Filter neoantigens based on nove, novel & SB, novel and IRP, novel, SB and IRP

#filter 1: Novel neoantigens----
SPCRCnovel_snv<- filter(SPCRCneoreco_snv, Novelty!=0)
#deduplicate neoantigens
SPCRCnovel_snv.dedup<- SPCRCnovel_snv[!duplicated(SPCRCnovel_snv$AntigenID),]
SPCRCnovel_snv.dedup<- data.frame(table(SPCRCnovel_snv.dedup$AntigenID))

#column names
colnames(SPCRCnovel_snv.dedup)=c("Sample","SPCRC_novel_snv")
#sub sample names
SPCRCnovel_snv.dedup$Sample <- gsub(":.*$", "", SPCRCnovel_snv.dedup$Sample)

#sum the values for similar samples
SPCRCnovel_snv.dedup<- ddply(SPCRCnovel_snv.dedup,"Sample",numcolwise(sum))


#Filter 2: novel and SBs (snvs)----
SPCRCnovel_sb_snv<-filter(SPCRCneoreco_snv,(Novelty!=0 & BindLevel=="SB"))
#dedulpicate
SPCRCnovel_sb_snv.dedup<-SPCRCnovel_sb_snv[!duplicated(SPCRCnovel_sb_snv$AntigenID),]
SPCRCnovel_sb_snv.dedup<-data.frame(table(SPCRCnovel_sb_snv.dedup$AntigenID))
#column names
colnames(SPCRCnovel_sb_snv.dedup)=c("Sample","SPCRC_novel_SB_snv")
#sub sample names
SPCRCnovel_sb_snv.dedup$Sample <- gsub(":.*$", "", SPCRCnovel_sb_snv.dedup$Sample)
#sum the values for similar samples
SPCRCnovel_sb_snv.dedup<- ddply(SPCRCnovel_sb_snv.dedup,"Sample",numcolwise(sum))



#filter 3: novel, IRP (snvs)----
SPCRCnovel_IRP_snv<-filter(SPCRCneoreco_snv, (Novelty!=0 & NeoantigenRecognitionPotential>1e-1))
#deduplicate
SPCRCnovel_IRP_snv.dedup<-SPCRCnovel_IRP_snv[!duplicated(SPCRCnovel_IRP_snv$AntigenID),]
SPCRCnovel_IRP_snv.dedup<-data.frame(table(SPCRCnovel_IRP_snv.dedup$AntigenID))
#column names
colnames(SPCRCnovel_IRP_snv.dedup)=c("Sample","SPCRC_novel_IRP_snv")
#sub sample names
SPCRCnovel_IRP_snv.dedup$Sample <- gsub(":.*$", "", SPCRCnovel_IRP_snv.dedup$Sample)
#sum the values for similar samples
SPCRCnovel_IRP_snv.dedup<- ddply(SPCRCnovel_IRP_snv.dedup,"Sample",numcolwise(sum))
SPCRCnovel_IRP_snv.dedup


# Plot recognition potential values for snvs
ggplot(SPCRCneoreco_snv, aes(x=NeoantigenRecognitionPotential)) + geom_density(fill='#9b8049', alpha=0.7) +
  scale_x_log10(limits=c(1e-6,1e3)) +
  theme_bw() +
  labs(x='Recognition Potential', y='Density') + geom_vline(xintercept = 1e-1, colour='firebrick', size=1.2, linetype='dashed')



#filter 4: novel,sb and IRP (snvs)----
SPCRCnovel_sb_IRP_snv<-filter(SPCRCneoreco_snv, (Novelty!=0 & BindLevel=="SB" & NeoantigenRecognitionPotential>1e-1))
#deduplicate
SPCRCnovel_sb_IRP_snv.dedup<-SPCRCnovel_sb_IRP_snv[!duplicated(SPCRCnovel_sb_IRP_snv$AntigenID),]
SPCRCnovel_sb_IRP_snv.dedup<-data.frame(table(SPCRCnovel_sb_IRP_snv.dedup$AntigenID))
#column names
colnames(SPCRCnovel_sb_IRP_snv.dedup)=c("Sample","SPCRC_novel_sb_IRP_snv")
#sub sample names
SPCRCnovel_sb_IRP_snv.dedup$Sample <- gsub(":.*$", "", SPCRCnovel_sb_IRP_snv.dedup$Sample)
#sum the values for similar samples
SPCRCnovel_sb_IRP_snv.dedup<- ddply(SPCRCnovel_sb_IRP_snv.dedup,"Sample",numcolwise(sum))



#merge all snvs filtered values together----
SPCRC_novel_novel_sb<-full_join (SPCRCnovel_snv.dedup, SPCRCnovel_sb_snv.dedup, by= "Sample") 
SPCRC_novelsb_novelIRP<-full_join(SPCRCnovel_IRP_snv.dedup, SPCRCnovel_sb_IRP_snv.dedup, by="Sample")
SPCRC_fil_snvs<-full_join(SPCRC_novel_novel_sb,SPCRC_novelsb_novelIRP, by="Sample")

#Determine total mutations input into the NeoPredPipe for normalisation

dir= '/Users/syedaanamfatima/Downloads/Neopredpipe/summary_table/spCRCs/fastaFilessp/'
fastafiles<- list.files(dir, pattern="*.tmp.10.fasta")
#read in all the snv mutations fo rfast amutations
fasta<-
  for (i in 1:length(fastafiles)){
    temp <- readLines(fastafiles[i])
    ff<-sum(grepl('>', temp))
    b<-print(paste0(gsub('.tmp.10.fasta', ' ',fastafiles[i])," ",ff))
    c<-lapply(b, write, "/Users/syedaanamfatima/Downloads/Neopredpipe/summary_table/spCRCs/fastaFilessp/SPCRC_snv_fasta_count.txt", append=TRUE, ncolumns=1000) 
    
  }


#add column names to the table
SPCRCsnv<- read.table('/Users/syedaanamfatima/Downloads/Neopredpipe/summary_table/spCRCs/fastaFilessp/SPCRC_snv_fasta_count.txt',
                      stringsAsFactors = F, header=T)
colnames(SPCRCsnv)=c("Sample","SPCRC_snv_fasta_mutations")
write.table(SPCRCsnv, file= "/Users/syedaanamfatima/Downloads/Neopredpipe/summary_table/spCRCs/fastaFilessp/SPCRC_snv_fasta_count.txt", sep="\t", quote = FALSE, row.names = FALSE)
SPCRC<-full_join(SPCRC_fil_snvs, SPCRCsnv, by="Sample")



#SPCRC NORMALISATION for each filteration level by dividing the neoantigen count by the total number of mutations----
#1. Novel(snvs)
SPCRC<- SPCRC %>%
  mutate(SPCRC_nml_novel_snv = SPCRC_novel_snv / SPCRC_snv_fasta_mutations)



#2. Novel, SBs
SPCRC<- SPCRC %>%
  mutate(SPCRC_nml_novel_sb_snv = SPCRC_novel_SB_snv / SPCRC_snv_fasta_mutations)



#3. NOVEL, IRP
SPCRC<- SPCRC %>%
  mutate(SPCRC_nml_novel_IRP_snv = SPCRC_novel_IRP_snv / SPCRC_snv_fasta_mutations)



#4. novel,sb,irp
SPCRC<- SPCRC %>%
  mutate(SPCRC_nml_novel_sb_IRP_snv= SPCRC_novel_sb_IRP_snv/SPCRC_snv_fasta_mutations)


SPCRC<- add_column(SPCRC, cancerType = "SPCRC")
write.table(SPCRC, file= "/Users/syedaanamfatima/Downloads/Neopredpipe/summary_table/SPCRC.txt", sep="\t", quote = FALSE, row.names = FALSE)


#calculate means
SPCRC_mean_count <- SPCRC %>%
  select(1, 7 :10) %>% #the last 8 columns #7 :10 #12:19
  summarise_at(vars(2:5), tibble::lst(base::mean), na.rm = TRUE) #2:5 #2:9


###Neopred antigen table
#snv table
dir = '/Users/syedaanamfatima/Downloads/Neopredpipe/summary_table/TCGA'

TCGAneo_snv<-read.table(paste0(dir, '/Total.neoantigens.txt'), header=T,
                        sep = '\t',stringsAsFactors = F, fill=T)
names(TCGAneo_snv) <- c('Sample', 'LineID', 'Chrom', 'allelepos',
                        'RefAll', 'AltAll', 'Genename:refseqID', 'pos', 'hla', 'peptide', 'core', 'Of', 'Gp',
                        'Gl', 'Ip', 'Il', 'Icore', 'Identity', 'Score', 'Affinity','Rank', 'Cand', 'BindLevel','Novelty')

#create antigen ID for each table based on Sample and Identity
TCGAneo_snv$AntigenID<-apply(TCGAneo_snv, 1,function(x) paste0(x['Sample'], ':',x['Identity']) )



#TCGA recopo table to filter neoantigens above 1 x 10^-1
#snv table
TCGAreco_snv <- read.table('/Users/syedaanamfatima/Downloads/Neopredpipe/summary_table/TCGA/PredictedRecognitionPotentials.txt',
                           stringsAsFactors = F, header=T)


#create antigen ID for each table based on Sample and Mutation
TCGAreco_snv$AntigenID<-apply(TCGAreco_snv, 1,function(x) paste0(x['Sample'], ':',x['Mutation']))


#merge neoantigen and recopo tables for snvs together to allow filtration from both tables

TCGAneoreco_snv<-full_join(TCGAneo_snv,TCGAreco_snv, by="AntigenID")
TCGAneoreco_ind<-full_join(TCGAneo_indel,TCGAreco_ind, by="AntigenID")

#SNVs
#filter 1: Novel neoantigens----
novel_snv<- filter(TCGAneoreco_snv, Novelty!=0)
#deduplicate neoantigens
novel_snv.dedup<- novel_snv[!duplicated(novel_snv$AntigenID),]
novel_snv.dedup<- data.frame(table(novel_snv.dedup$AntigenID))
#column names
colnames(novel_snv.dedup)=c("Sample","TCGA_novel_snv")
#sub sample names
novel_snv.dedup$Sample <- gsub(":.*$", "", novel_snv.dedup$Sample)
#sum the values for similar samples
novel_snv.dedup<- ddply(novel_snv.dedup,"Sample",numcolwise(sum))
novel_snv.dedup



#Filter 2: novel and SBs (snvs)----
novel_sb_snv<-filter(TCGAneoreco_snv,(Novelty!=0 & BindLevel=="SB"))
#dedulpicate
novel_sb_snv.dedup<-novel_sb_snv[!duplicated(novel_sb_snv$AntigenID),]
novel_sb_snv.dedup<-data.frame(table(novel_sb_snv.dedup$AntigenID))
#column names
colnames(novel_sb_snv.dedup)=c("Sample","TCGA_novel_SB_snv")
#sub sample names
novel_sb_snv.dedup$Sample <- gsub(":.*$", "", novel_sb_snv.dedup$Sample)
#sum the values for similar samples
novel_sb_snv.dedup<- ddply(novel_sb_snv.dedup,"Sample",numcolwise(sum))



#filter 3: novel, IRP (snvs)----
novel_IRP_snv<-filter(TCGAneoreco_snv, (Novelty!=0 & NeoantigenRecognitionPotential>1e-1))
#deduplicate
novel_IRP_snv.dedup<-novel_IRP_snv[!duplicated(novel_IRP_snv$AntigenID),]
novel_IRP_snv.dedup<-data.frame(table(novel_IRP_snv.dedup$AntigenID))
#column names
colnames(novel_IRP_snv.dedup)=c("Sample","TCGA_novel_IRP_snv")
#sub sample names
novel_IRP_snv.dedup$Sample <- gsub(":.*$", "", novel_IRP_snv.dedup$Sample)
#sum the values for similar samples
novel_IRP_snv.dedup<- ddply(novel_IRP_snv.dedup,"Sample",numcolwise(sum))
novel_IRP_snv.dedup


# Plot recognition potential values for snvs
ggplot(TCGAneoreco_snv, aes(x=NeoantigenRecognitionPotential)) + geom_density(fill='#9b8049', alpha=0.7) +
  scale_x_log10(limits=c(1e-6,1e3)) +
  theme_bw() +
  labs(x='Recognition Potential', y='Density') + geom_vline(xintercept = 1e-1, colour='firebrick', size=1.2, linetype='dashed')



#filter 4: novel,sb and IRP (snvs)----
novel_sb_IRP_snv<-filter(TCGAneoreco_snv, (Novelty!=0 & BindLevel=="SB" & NeoantigenRecognitionPotential>1e-1))
#deduplicate
novel_sb_IRP_snv.dedup<-novel_sb_IRP_snv[!duplicated(novel_sb_IRP_snv$AntigenID),]
novel_sb_IRP_snv.dedup<-data.frame(table(novel_sb_IRP_snv.dedup$AntigenID))
#column names
colnames(novel_sb_IRP_snv.dedup)=c("Sample","TCGA_novel_sb_IRP_snv")
#sub sample names
novel_sb_IRP_snv.dedup$Sample <- gsub(":.*$", "", novel_sb_IRP_snv.dedup$Sample)
#sum the values for similar samples
novel_sb_IRP_snv.dedup<- ddply(novel_sb_IRP_snv.dedup,"Sample",numcolwise(sum))


#merge all snvs filtered values together----
TCGA_novel_novel_sb<-full_join (novel_snv.dedup, novel_sb_snv.dedup, by= "Sample") 
TCGA_novelsb_novelIRP<-full_join(novel_IRP_snv.dedup, novel_sb_IRP_snv.dedup, by="Sample")
TCGA_fil_snvs<-full_join(TCGA_novel_novel_sb,TCGA_novelsb_novelIRP, by="Sample")



#antigen normalisation TCGA
### number of fasta headers in each .tmp.10.fasta

dir= '/Users/syedaanamfatima/Downloads/Neopredpipe/summary_table/TCGA/TCGA_fasta'
fastafiles<- list.files(dir, pattern="*.tmp.10.fasta")
fastafiles

tcfasta<-
  for (i in 1:length(fastafiles)){
    temp <- readLines(fastafiles[i])
    ff<-sum(grepl('>', temp))
    b<-print(paste0(gsub('.tmp.10.fasta', ' ',fastafiles[i])," ",ff))
    c<-lapply(b, write, "/Users/syedaanamfatima/Downloads/Neopredpipe/summary_table/TCGA_fasta_count.txt", append=TRUE, ncolumns=1000) 
  }   


fastafilesindel<- list.files(dir, pattern="*.tmp.10.Indels.fasta")
fastafilesindel
tcfastaindel<-
  for (i in 1:length(fastafilesindel)){
    temp <- readLines(fastafilesindel[i])
    ff_indel<-sum(grepl('>', temp))
    b_indel<-print(paste0(gsub('.tmp.10.Indels.fasta', ' ',fastafilesindel[i])," ",ff_indel))
    c<-lapply(b_indel,write,"/Users/syedaanamfatima/Downloads/Neopredpipe/summary_table/TCGA_ind_fasta_count.txt", append=TRUE, ncolumns=1000) 
    
  }

#merge snv and indel fasta header counts files together
snv<- read.table('/Users/syedaanamfatima/Downloads/Neopredpipe/summary_table/TCGA_fasta_count.txt',
                 stringsAsFactors = F, header=T)
colnames(snv)=c("Sample","TCGA_snv_fasta_mutations")
write.table(snv, file= "/Users/syedaanamfatima/Downloads/Neopredpipe/summary_table/TCGA_fasta_count.txt", sep="\t", quote = FALSE, row.names = FALSE)

indel<-read.table('/Users/syedaanamfatima/Downloads/Neopredpipe/summary_table/TCGA_ind_fasta_count.txt',
                  stringsAsFactors = F, header=T)
colnames(indel)=c("Sample","TCGA_indel_fasta_mutations")
write.table(indel, file= "/Users/syedaanamfatima/Downloads/Neopredpipe/summary_table/TCGA_ind_fasta_count.tx", sep="\t", quote = FALSE, row.names = FALSE)


snv_indel_fasta<-full_join(snv,indel, by="Sample")
snv_indel_fasta

#join neoantigen score for eac level with total mutation count for snvs and indels
TCGA<-full_join(TCGA_final, snv_indel_fasta, by="Sample")
TCGA


#TCGA NORMALISATION for each filteration level----
#1. Novel(snvs)

TCGA<- TCGA %>%
  mutate(TCGA_nml_novel_snv = TCGA_novel_snv / TCGA_snv_fasta_mutations)


#2. Novel, SBs
TCGA<- TCGA %>%
  mutate(TCGA_nml_novel_sb_snv = TCGA_novel_SB_snv / TCGA_snv_fasta_mutations)


#3. NOVEL, IRP
TCGA<- TCGA %>%
  mutate(TCGA_nml_novel_IRP_snv = TCGA_novel_IRP_snv / TCGA_snv_fasta_mutations)


#4. novel,sb,irp
TCGA<- TCGA %>%
  mutate(TCGA_nml_novel_sb_IRP_snv= TCGA_novel_sb_IRP_snv/TCGA_snv_fasta_mutations)

TCGA<- add_column(TCGA, cancerType = "TCGA")
TCGA$TCGA_nml_novel_IRP_ind
value<-subset(TCGA, (TCGA$TCGA_nml_novel_IRP_ind > 0.9))
write.table(TCGA, file= "/Users/syedaanamfatima/Downloads/Neopredpipe/summary_table/TCGA.txt", sep="\t", quote = FALSE, row.names = FALSE)


#calculate means
TCGA_mean_count <- TCGA %>%
  select(1, 12:19) %>% #the last 8 columns
  summarise_at(vars(2:9), tibble::lst(base::mean), na.rm = TRUE) 







#####################  CACRC_s2s Filtration & Normalisation----######################################

#specify directory for the summary table
dir = '/Users/syedaanamfatima/Downloads/Neopredpipe/summary_table/jay'

#Neoantigen tables
#snv
CACRC_s2_neo_snv<-read.table(paste0(dir, '/Total.neoantigens.txt'), header=F,
                          sep = '\t',stringsAsFactors = F, fill=T)

names(CACRC_s2_neo_snv) <- c('Sample','R1','R2','R3','4','5','6','7','8','9','10','11','LineID', 'Chrom', 'allelepos',
                          'RefAll', 'AltAll', 'Genename:refseqID', 'pos', 'hla', 'peptide', 'core', 'Of', 'Gp',
                          'Gl', 'Ip', 'Il', 'Icore', 'Identity', 'Score', 'Affinity','Rank', 'Cand', 'BindLevel','Novelty')


#create antigen ID for each table based on Sample and Identity
CACRC_s2_neo_snv$AntigenID<-apply(CACRC_s2_neo_snv, 1,function(x) paste0(x['Sample'], ':',x['Identity']) )


#recopo table (snvs)
CACRC_s2reco_snv  <- read.table('/Users/syedaanamfatima/Downloads/Neopredpipe/summary_table/jay/PredictedRecognitionPotentials.txt',
                             stringsAsFactors = F, header=T)

#create antigen ID for each table based on Sample and Mutation
CACRC_s2reco_snv$AntigenID<-apply(CACRC_s2reco_snv, 1,function(x) paste0(x['Sample'], ':',x['Mutation']))

#merge neoantigen and recopo tables for snvs together to allow filtration from both tables
CACRC_s2neoreco_snv<-full_join(CACRC_s2_neo_snv,CACRC_s2reco_snv, by="AntigenID")




#SNVs
#filter 1: Novel neoantigens----
CACRC_s2novel_snv<- filter(CACRC_s2neoreco_snv, Novelty!=0)
#deduplicate neoantigens
CACRC_s2novel_snv.dedup<- CACRC_s2novel_snv[!duplicated(CACRC_s2novel_snv$AntigenID),]
CACRC_s2novel_snv.dedup<- data.frame(table(CACRC_s2novel_snv.dedup$AntigenID))
#column names
colnames(CACRC_s2novel_snv.dedup)=c("Sample","CACRC_s2_novel_snv")
#sub sample names
CACRC_s2novel_snv.dedup$Sample <- gsub(":.*$", "", CACRC_s2novel_snv.dedup$Sample)
#sum the values for similar samples
CACRC_s2novel_snv.dedup<- ddply(CACRC_s2novel_snv.dedup,"Sample",numcolwise(sum))





#Filter 2: novel and SBs (snvs)----
CACRC_s2novel_sb_snv<-filter(CACRC_s2neoreco_snv,(Novelty!=0 & BindLevel=="SB"))
#dedulpicate
CACRC_s2novel_sb_snv.dedup<-CACRC_s2novel_sb_snv[!duplicated(CACRC_s2novel_sb_snv$AntigenID),]
CACRC_s2novel_sb_snv.dedup<-data.frame(table(CACRC_s2novel_sb_snv.dedup$AntigenID))
#column names
colnames(CACRC_s2novel_sb_snv.dedup)=c("Sample","CACRC_s2_novel_SB_snv")
#sub sample names
CACRC_s2novel_sb_snv.dedup$Sample <- gsub(":.*$", "", CACRC_s2novel_sb_snv.dedup$Sample)
#sum the values for similar samples
CACRC_s2novel_sb_snv.dedup<- ddply(CACRC_s2novel_sb_snv.dedup,"Sample",numcolwise(sum))




#filter 3: novel, IRP (snvs)----
CACRC_s2novel_IRP_snv<-filter(CACRC_s2neoreco_snv, (Novelty!=0 & NeoantigenRecognitionPotential>1e-1))

#deduplicate
CACRC_s2novel_IRP_snv.dedup<-CACRC_s2novel_IRP_snv[!duplicated(CACRC_s2novel_IRP_snv$AntigenID),]
CACRC_s2novel_IRP_snv.dedup<-data.frame(table(CACRC_s2novel_IRP_snv.dedup$AntigenID))
#column names
colnames(CACRC_s2novel_IRP_snv.dedup)=c("Sample","CACRC_s2_novel_IRP_snv")
#sub sample names
CACRC_s2novel_IRP_snv.dedup$Sample <- gsub(":.*$", "", CACRC_s2novel_IRP_snv.dedup$Sample)
#sum the values for similar samples
CACRC_s2novel_IRP_snv.dedup<- ddply(CACRC_s2novel_IRP_snv.dedup,"Sample",numcolwise(sum))



# Plot recognition potential values for snvs
ggplot(CACRC_s2neoreco_snv, aes(x=NeoantigenRecognitionPotential)) + geom_density(fill='#9b8049', alpha=0.7) +
  scale_x_log10(limits=c(1e-6,1e3)) +
  theme_bw() +
  labs(x='Recognition Potential', y='Density') + geom_vline(xintercept = 1e-1, colour='firebrick', size=1.2, linetype='dashed')




#filter 4: novel,sb and IRP (snvs)----
CACRC_s2novel_sb_IRP_snv<-filter(CACRC_s2neoreco_snv, (Novelty!=0 & BindLevel=="SB" & NeoantigenRecognitionPotential>1e-1))
#deduplicate
CACRC_s2novel_sb_IRP_snv.dedup<-CACRC_s2novel_sb_IRP_snv[!duplicated(CACRC_s2novel_sb_IRP_snv$AntigenID),]
CACRC_s2novel_sb_IRP_snv.dedup<-data.frame(table(CACRC_s2novel_sb_IRP_snv.dedup$AntigenID))
#column names
colnames(CACRC_s2novel_sb_IRP_snv.dedup)=c("Sample","CACRC_s2_novel_sb_IRP_snv")
#sub sample names
CACRC_s2novel_sb_IRP_snv.dedup$Sample <- gsub(":.*$", "", CACRC_s2novel_sb_IRP_snv.dedup$Sample)
#sum the values for similar samples
CACRC_s2novel_sb_IRP_snv.dedup<- ddply(CACRC_s2novel_sb_IRP_snv.dedup,"Sample",numcolwise(sum))


#merge all snvs filtered values together----
CACRC_s2_novel_novel_sb<-full_join (CACRC_s2novel_snv.dedup, CACRC_s2novel_sb_snv.dedup, by= "Sample") 
CACRC_s2_novelsb_novelIRP<-full_join(CACRC_s2novel_IRP_snv.dedup, CACRC_s2novel_sb_IRP_snv.dedup, by="Sample")
CACRC_s2_fil_snvs<-full_join(CACRC_s2_novel_novel_sb,CACRC_s2_novelsb_novelIRP, by="Sample")



### number of fasta headers in each .tmp.10.fasta

dir= '/Users/syedaanamfatima/Downloads/Neopredpipe/summary_table/jay/fastaFiles/'
fastafiles<- list.files(dir, pattern="*.tmp.10.fasta")
# read in all the snv fasta muattaions
fasta_snv<-
  for (i in 1:length(fastafiles)){
    temp <- readLines(fastafiles[i])
    ff<-sum(grepl('>', temp))
    b<-print(paste0(gsub('.tmp.10.fasta', ' ',fastafiles[i])," ",ff))
    c<-lapply(b, write, "/Users/syedaanamfatima/Downloads/Neopredpipe/summary_table/CACRs_snv_mut.txt", append=TRUE, ncolumns=1000) 
    
  }


#Add column names to total mutation counts
ca_snv<- read.table('/Users/syedaanamfatima/Downloads/Neopredpipe/summary_table/CACRs_snv_mut.txt',
                    stringsAsFactors = F, header=T)
colnames(ca_snv)=c("Sample","CACRC_s2_snv_mutations")

write.table(ca_snv, file= "/Users/syedaanamfatima/Downloads/Neopredpipe/summary_table/CACRs_sp_snv_mut.txt", sep="\t", quote = FALSE, row.names = FALSE)





#join neoantigen score for each filtration level with total mutation count
CACRC_s2<-full_join(CACRC_s2_fil_snvs, ca_snv, by="Sample")
#CACRC_s2<-CACRC_s2_fil_snvs
CACRC_s2
#CACRC_s2s NORMALISATION for each filteration level----
#1. Novel(snvs)
CACRC_s2<- CACRC_s2 %>%
  mutate(CACRC_s2_nml_novel_snv = CACRC_s2_novel_snv / CACRC_s2_snv_mutations)


#2. Novel, SBs
CACRC_s2<- CACRC_s2 %>%
  mutate(CACRC_s2_nml_novel_sb_snv = CACRC_s2_novel_SB_snv / CACRC_s2_snv_mutations)
CACRC_s2$CACRC_s2_snv_mutations

#3. NOVEL, IRP
CACRC_s2<- CACRC_s2 %>%
  mutate(CACRC_s2_nml_novel_IRP_snv = CACRC_s2_novel_IRP_snv / CACRC_s2_snv_mutations)



#4. novel,sb,irp
CACRC_s2<- CACRC_s2 %>%
  mutate(CACRC_s2_nml_novel_sb_IRP_snv= CACRC_s2_novel_sb_IRP_snv/CACRC_s2_snv_mutations)

CACRC_s2<- add_column(CACRC_s2, cancerType = "CA-CRCs")
CACRC_s2
write.table(CACRC_s2, file= "/Users/syedaanamfatima/Downloads/Neopredpipe/summary_table/CACRC_s2.txt", sep="\t", quote = FALSE, row.names = FALSE)


#calculate means nml for wilcoxon test-------
CACRC_s2_mean_count <- CACRC_s2 %>%
  select(1, 7:10) %>% #the last 8 columns
  summarise_at(vars(2:5), tibble::lst(base::mean), na.rm = TRUE) 


###################### Box plot for neoantigen number at each filtration level ################

# set plot theme
ggthemr('flat', layout = 'clean', spacing = 1)

#1. novel, SB & immunogenic----

#for TCGAs--
TCGA
TCGA_neoAg_normCount_SBIRP <- TCGA %>%
  gather(key = key, value = value,
         TCGA_nml_novel_sb_IRP_snv) %>%#, TCGA_nml_novel_sb_ind) %>%
  arrange(Sample) %>%
  select(-TCGA_novel_snv, -TCGA_novel_SB_snv, -TCGA_novel_IRP_snv, -TCGA_novel_sb_IRP_snv,
         -TCGA_snv_fasta_mutations,  -TCGA_nml_novel_snv,
         -TCGA_nml_novel_sb_snv, -TCGA_nml_novel_IRP_snv )%>%
  mutate(varType = sub(".*nml_", "", key))



#for SPCRC--
SPCRC
SPCRC_neoAg_normCount_SBIRP <- SPCRC %>%
  gather(key = key, value = value,
         SPCRC_nml_novel_sb_IRP_snv)%>%
  arrange(Sample) %>%  
  select(-SPCRC_novel_snv, -SPCRC_novel_SB_snv, -SPCRC_novel_IRP_snv, -SPCRC_novel_sb_IRP_snv,
         -SPCRC_snv_fasta_mutations,  -SPCRC_nml_novel_snv, 
         -SPCRC_nml_novel_sb_snv,  -SPCRC_nml_novel_IRP_snv)%>%
  mutate(varType = sub(".*nml_", "", key))



#for IBD- ass. CRCs
CACRC_s2

CACRC_s2_neoAg_normCount_SBIRP <- CACRC_s2 %>%
  gather(key = key, value = value,
         CACRC_s2_nml_novel_sb_IRP_snv)%>%
  arrange(Sample) %>%
  select(-CACRC_s2_novel_snv, -CACRC_s2_novel_SB_snv, -CACRC_s2_novel_IRP_snv, -CACRC_s2_novel_sb_IRP_snv,
         -CACRC_s2_snv_mutations,  -CACRC_s2_nml_novel_snv,
         -CACRC_s2_nml_novel_sb_snv, -CACRC_s2_nml_novel_IRP_snv)%>%
  mutate(varType = sub(".*nml_", "", key))
CACRC_s2_neoAg_normCount_SBIRP
neoAg_normCounts_SBIRP <- bind_rows(CACRC_s2_neoAg_normCount_SBIRP,SPCRC_neoAg_normCount_SBIRP,TCGA_neoAg_normCount_SBIRP)
neoAg_normCounts_SBIRP


##plotting box plot to compare distribution between FF-CACRCs, spCRCs and sp-TCGA
my_comparisons <- list( c("SPCRC_nml_novel_sb_IRP_snv", "CACRC_s2_nml_novel_sb_IRP_snv"), c("CACRC_s2_nml_novel_sb_IRP_snv","TCGA_nml_novel_sb_IRP_snv"))#,
give.n <- function(x){
  return(c(y = mean(x), label = length(x)))
}
p5<-ggplot() +
  geom_boxplot() +
  aes(x = neoAg_normCounts_SBIRP$key, y = neoAg_normCounts_SBIRP$value, fill = neoAg_normCounts_SBIRP$cancerType) +
  geom_jitter(color="purple", size=0.4, alpha=0.9) +
  scale_fill_manual(breaks = neoAg_normCounts_SB$value, values = c("maroon", "green","gold"))+
  theme(text=element_text(size=14)) +
  stat_compare_means(comparisons = my_comparisons, method = 'wilcox.test', label.y=c(0.15,0.18),check_overlap = T)+
  xlab('Sample Type') +
  ylab('Normalised neoantigen count') +
  ggtitle('', subtitle = 'Novel, strong-binder & immunogenic neoantigens') +
  theme(plot.title = element_text(face = 'italic')) +
  labs(fill = 'Cancer type') +
  scale_x_discrete(labels =c("FF-CACRCs","spCRCs","TCGA"))
#save plot
ggsave(plot = p5, width = 5, height = 4, dpi = 300, filename = "/Users/syedaanamfatima/Downloads/Neopredpipe/pdfs/s2nsirp.pdf")

# 2. novel & SB-------------
TCGA_neoAg_normCount_SB <- TCGA %>%
  gather(key = key, value = value,
         TCGA_nml_novel_sb_snv) %>%
  arrange(Sample) %>%
  select(-TCGA_novel_snv, -TCGA_novel_SB_snv, -TCGA_novel_IRP_snv, -TCGA_novel_sb_IRP_snv,
         -TCGA_snv_fasta_mutations, -TCGA_nml_novel_snv,
         -TCGA_nml_novel_IRP_snv, -TCGA_nml_novel_sb_IRP_snv, )%>%
  mutate(varType = sub(".*nml_", "", key))


#for SPCRCs--
SPCRC
SPCRC_neoAg_normCount_SB <- SPCRC %>%
  gather(key = key, value = value,
         SPCRC_nml_novel_sb_snv)%>%
  arrange(Sample) %>%
  select(-SPCRC_novel_snv, -SPCRC_novel_SB_snv, -SPCRC_novel_IRP_snv, -SPCRC_novel_sb_IRP_snv,
         -SPCRC_snv_fasta_mutations, -SPCRC_nml_novel_snv,
         -SPCRC_nml_novel_IRP_snv, -SPCRC_nml_novel_sb_IRP_snv)%>%
  mutate(varType = sub(".*nml_", "", key))

#for CARCs
CACRC_s2
CACRC_s2_neoAg_normCount_SB <- CACRC_s2 %>%
  gather(key = key, value = value,
         CACRC_s2_nml_novel_sb_snv)%>%
  arrange(Sample) %>%
  select(-CACRC_s2_novel_snv, -CACRC_s2_novel_SB_snv, -CACRC_s2_novel_IRP_snv, -CACRC_s2_novel_sb_IRP_snv,
        
         -CACRC_s2_snv_mutations, -CACRC_s2_nml_novel_snv, 
         -CACRC_s2_nml_novel_IRP_snv, -CACRC_s2_nml_novel_sb_IRP_snv)%>%
  mutate(varType = sub(".*nml_", "", key))

neoAg_normCounts_SB <- bind_rows(SPCRC_neoAg_normCount_SB, 
                                 CACRC_s2_neoAg_normCount_SB,TCGA_neoAg_normCount_SB)



#plottingSB
my_comparisons <- list( c("SPCRC_nml_novel_sb_snv", "CACRC_s2_nml_novel_sb_snv"),c("CACRC_s2_nml_novel_sb_snv",'TCGA_nml_novel_sb_snv'))#,
give.n <- function(x){
  return(c(y = mean(x), label = length(x)))
}
p6<-ggplot() +
  geom_boxplot() +
  aes(x = neoAg_normCounts_SB$key, y = neoAg_normCounts_SB$value, fill = neoAg_normCounts_SB$cancerType) +
  scale_fill_manual(breaks = neoAg_normCounts_SB$value, values = c("maroon", "green","gold"))+
  geom_jitter(color="purple", size=0.4, alpha=0.9) +
  theme(text=element_text(size=14)) +

  stat_compare_means(comparisons = my_comparisons,method = 'wilcox.test', label.y = c(0.70,0.75),check_overlap = T)+
  xlab('Sample type') +
  ylab('Normalised neoantigen count') +
  ggtitle('Novel & strong-binder neoantigens')+
  theme(plot.title = element_text(face = 'italic')) +
  labs(fill = 'Cancer type') +
  scale_x_discrete(labels =c("FF-CACRCs","spCRCs","TCGA"))
ggsave(plot = p6, width = 5, height = 4, dpi = 300, filename = "/Users/syedaanamfatima/Downloads/Neopredpipe/pdfs/s2nS.pdf")


# 3. novel & IRP----
#TCGA
TCGA_neoAg_normCount_IRP <- TCGA %>%
  gather(key = key, value = value,
         TCGA_nml_novel_IRP_snv) %>%
  arrange(Sample) %>%
  select(-TCGA_novel_snv, -TCGA_novel_SB_snv, -TCGA_novel_IRP_snv, -TCGA_novel_sb_IRP_snv,
         -TCGA_snv_fasta_mutations, -TCGA_indel_fasta_mutations, -TCGA_nml_novel_snv,
         -TCGA_nml_novel_sb_snv, -TCGA_nml_novel_sb_IRP_snv)%>%
  mutate(varType = sub(".*nml_", "", key))

#for spCRCs
SPCRC
SPCRC_neoAg_normCount_IRP <- SPCRC %>%
  gather(key = key, value = value,
         SPCRC_nml_novel_IRP_snv) %>%# ,SPCRC_nml_novel_IRP_ind) %>%
  arrange(Sample) %>%
  select(-SPCRC_novel_snv, -SPCRC_novel_SB_snv, -SPCRC_novel_IRP_snv, -SPCRC_novel_sb_IRP_snv,
         -SPCRC_snv_fasta_mutations, -SPCRC_nml_novel_snv, 
         -SPCRC_nml_novel_sb_snv, -SPCRC_nml_novel_sb_IRP_snv)%>%
  mutate(varType = sub(".*nml_", "", key))


#for CRCs
CACRC_s2
CACRC_s2_neoAg_normCount_IRP <- CACRC_s2%>%
  gather(key = key, value = value,
         CACRC_s2_nml_novel_IRP_snv) %>%
  arrange(Sample) %>%
  select(-CACRC_s2_novel_snv, -CACRC_s2_novel_SB_snv, -CACRC_s2_novel_IRP_snv, -CACRC_s2_novel_sb_IRP_snv,
         
         -CACRC_s2_snv_mutations, -CACRC_s2_nml_novel_snv, 
         -CACRC_s2_nml_novel_sb_snv, -CACRC_s2_nml_novel_sb_IRP_snv)%>%
  mutate(varType = sub(".*nml_", "", key))
neoAg_normCounts_IRP <- bind_rows(SPCRC_neoAg_normCount_IRP, CACRC_s2_neoAg_normCount_IRP,TCGA_neoAg_normCount_IRP)

#plotting it  
my_comparisons <- list( c("SPCRC_nml_novel_IRP_snv", "CACRC_s2_nml_novel_IRP_snv"), c("CACRC_s2_nml_novel_IRP_snv","TCGA_nml_novel_IRP_snv"))#,
give.n <- function(x){
  return(c(y = mean(x), label = length(x)))
}
p7<-ggplot() +
  geom_boxplot() +
  aes(x = neoAg_normCounts_IRP$key, y = neoAg_normCounts_IRP$value, fill = neoAg_normCounts_IRP$cancerType) +
  scale_fill_manual(breaks = neoAg_normCounts_IRP$value, values = c("maroon", "green","gold"))+
  geom_jitter(color="purple", size=0.4, alpha=0.9) +
  theme(text=element_text(size=14)) +
  
  stat_compare_means(comparisons = my_comparisons,method = 'wilcox.test', label.y = c(0.15,0.18),check_overlap = T)+
  xlab('Sample type') +
  ylab('Normalised neoantigen count') +
  ggtitle('', subtitle = 'Novel & immunogenic neoantigens') +
  theme(plot.title = element_text(face = 'italic')) +
  labs(fill = 'Cancer type') +
  scale_x_discrete(labels =c("FF-CACRCs","spCRCs","TCGA"))

ggsave(plot = p7, width = 5, height = 4, dpi = 300, filename = "/Users/syedaanamfatima/Downloads/Neopredpipe/pdfs/s2nirp.pdf")


# 4. novel----
TCGA_neoAg_normCount_novel <- TCGA %>%
  gather(key = key, value = value,
         TCGA_nml_novel_snv) %>%
  arrange(Sample) %>%
  select(-TCGA_novel_snv, -TCGA_novel_SB_snv, -TCGA_novel_IRP_snv, -TCGA_novel_sb_IRP_snv,
         -TCGA_novel_ind, -TCGA_novel_SB_ind, -TCGA_novel_IRP_ind, -TCGA_novel_sb_IRP_ind,
         -TCGA_snv_fasta_mutations, -TCGA_indel_fasta_mutations,
         -TCGA_nml_novel_sb_snv, -TCGA_nml_novel_sb_ind, -TCGA_nml_novel_IRP_snv, -TCGA_nml_novel_IRP_ind,
         -TCGA_nml_novel_sb_IRP_snv, -TCGA_nml_novel_sb_IRP_ind)%>%
  mutate(varType = sub(".*nml_", "", key))
#for SPCRCs
SPCRC
SPCRC_neoAg_normCount_novel <- SPCRC %>%
  gather(key = key, value = value,
         SPCRC_nml_novel_snv) %>%#, SPCRC_nml_novel_ind) %>%
  arrange(Sample) %>%
  select(-SPCRC_novel_snv, -SPCRC_novel_SB_snv, -SPCRC_novel_IRP_snv, -SPCRC_novel_sb_IRP_snv,
         -SPCRC_snv_fasta_mutations, 
         -SPCRC_nml_novel_sb_snv,  -SPCRC_nml_novel_IRP_snv,
         -SPCRC_nml_novel_sb_IRP_snv)%>%
  mutate(varType = sub(".*nml_", "", key))


#for CRCs
CACRC_s2
CACRC_s2_neoAg_normCount_novel <- CACRC_s2 %>%
  gather(key = key, value = value,
         CACRC_s2_nml_novel_snv) %>%
  arrange(Sample) %>%
  select(-CACRC_s2_novel_snv, -CACRC_s2_novel_SB_snv, -CACRC_s2_novel_IRP_snv, -CACRC_s2_novel_sb_IRP_snv,
         
         -CACRC_s2_snv_mutations, 
         -CACRC_s2_nml_novel_sb_snv,  -CACRC_s2_nml_novel_IRP_snv,
         -CACRC_s2_nml_novel_sb_IRP_snv)%>%
  mutate(varType = sub(".*nml_", "", key))

neoAg_normCounts_novel <- bind_rows(SPCRC_neoAg_normCount_novel, 
                                    CACRC_s2_neoAg_normCount_novel,TCGA_neoAg_normCount_novel)

#plotting comparisons between spCRCs, FF-CACRCs and TCGAs  
my_comparisons <- list( c("SPCRC_nml_novel_snv", "CACRC_s2_nml_novel_snv"),c("CACRC_s2_nml_novel_snv","TCGA_nml_novel_snv"))
give.n <- function(x){
  return(c(y = mean(x), label = length(x)))
}
p8<-ggplot() +
  geom_boxplot() +
  aes(x = neoAg_normCounts_novel$key, y = neoAg_normCounts_novel$value, fill = neoAg_normCounts_novel$cancerType) +
  scale_fill_manual(breaks = neoAg_normCounts_novel$value, values = c("maroon", "green","gold"))+
  geom_jitter(color="purple", size=0.4, alpha=0.9) +
  theme(text=element_text(size=14)) +
  stat_compare_means(comparisons = my_comparisons,method = 'wilcox.test', label.y = c(0.9, 0.95),check_overlap = T)+
  xlab('Sample type') +
  ylab('Normalised neoantigen count') +
  ggtitle('', subtitle = 'Novel neoantigens') +
  theme(plot.title = element_text(face = 'italic')) +
  labs(fill = 'Cancer type') +
  scale_x_discrete(labels =c("FF-CACRCs","spCRCs","TCGA"))

#save as a pdf
ggsave(plot = p8, width = 5, height = 4, dpi = 300, filename = "/Users/syedaanamfatima/Downloads/Neopredpipe/pdfs/s2n.pdf") 







