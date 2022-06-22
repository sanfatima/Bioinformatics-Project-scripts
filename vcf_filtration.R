# Install/load vcfR package
if (require("vcfR") == FALSE) {
  install.packages("vcfR")
  library("vcfR")
} else {
  library("vcfR")
}
library(gridExtra)
library(ggplot2)
library(dplyr)
library ("tidyverse")

#locate files 
vcf_file <- "/Users/syedaanamfatima/Downloads/Edinburgh_IBD_vcfs_2/9N.vcf"
#store vcf files in a list vector
#vcf_file <- list.files(vcfFileDir, pattern = "[.]vcf")
vcf_file 
proc_postfix <- 'filtered' # new name to show vcf file was processed
norm<-for (i in 1:length(vcf_file)) { 
  vcf<- read.vcfR(vcf_file, verbose=FALSE)
  # Filtering step 1: keep only mutations that are absent (or present in very low read-count) in the normal sample
  filterNorm <- function(normGT){
    normGT <- sapply(normGT, function(x) unlist(strsplit(x,':'))[2])
    refC <- as.numeric(map(strsplit(normGT, ','),1))
    altC <- as.numeric(map(strsplit(normGT, ','),2))
  # These are the numbers used in e.g. Will's papers and for EPICC
    somatic <- (refC>5 & refC<30 & altC < 1 | refC>29 & (altC+refC)<100 & altC < 2 | (altC+refC)>99 & altC<3)
    return(somatic)
  }
}


# Filtering step 2: keep only mutations that are high-confidence: have enough evidence suggesting they are true variants, not false positives
# can also plot diagnostic plots, showing how much of the mutations are kept in at the defined thresholds
# parameters are:
# - dpLim = >X number of reads (depth) at the variant position
# - acLim = >X variant reads (alternative count) at the position
# - vafLim = >X variant/total reads
tum<-for (i in 1:length(vcf_file)) { 
  vcf<- read.vcfR(vcf_file, verbose=FALSE)
  filterTumour <- function(tumGT, sampleName, dpLim, acLim, vafLim, plotting=F){
    # Filter out mutations with depth < dpLim, altCount < acLim and VAF < vafLim
      tumGT <- sapply(tumGT, function(x) unlist(strsplit(x,':'))[2])
      refC <- as.numeric(map(strsplit(tumGT, ','),1))
      altC <- as.numeric(map(strsplit(tumGT, ','),2))
      #somatic_tum <- (refC>5 & refC<30 & altC < 1 | refC>29 & (altC+refC)<100 & altC < 2 | (altC+refC)>99 & altC<3)
      ggplot(data.frame(x=(refC+altC)), aes(x=x)) +
        geom_histogram() + theme_bw() + scale_x_log10() +
        labs(x='Sequencing depth at pos', y='Number of mutations', title=sampleName) +
        geom_vline(xintercept = dpLim, colour='darkred')
      ggplot(data.frame(x=altC/(refC+altC)), aes(x=x)) +
        geom_histogram() + theme_bw() +
        labs(x='Variant Allele Frequency', y='Number of mutations', title=sampleName) +
        geom_vline(xintercept = vafLim, colour='darkred')
      ggplot(data.frame(r=refC,a=altC), aes(x=(r+a),y=a/(r+a), colour=a>acLim)) +
        labs(x='Sequencing depth at pos', y='Variant Allele Frequency', title=sampleName) +
        geom_point() + scale_x_log10() + theme_bw() +
        scale_colour_manual(values=c('grey20','darkred'))
      if(plotting){grid.arrange(p0,p1,p2)}
    
      filtered <- ((refC+altC) > dpLim & altC > acLim & altC/(refC+altC) > vafLim)
      return(filtered)
  }
}

# Keep only variants which PASSed filters in the vcf
vcf <- vcf[vcf@fix[,'FILTER']=='PASS',]
# Keep only variants with single nucleotide changes and biallelic sites
# This is because the filtering after this cannot process where there is more than 1 variant allele
# You could probably modify this to be able to process that too, and it is okay with NeoPredPipe
# You can check how many variants you throw away by doing nrow(vcf) and nrow(vcf.singlealt)
vcf.singlealt <- vcf[!grepl(',',vcf@fix[,'ALT']),]

# Get which column is the normal and use filter step 1

normName <- '9N_normal'#change this to name of the normal column in the vcf
vcf.somatic <- vcf.singlealt[filterNorm(vcf.singlealt@gt[,normName]),]

# Get which column is the tumor and use filter step 1

tumName<- "9N_tumour"
# Change the numbers here to change what cut-off is used for 1) coverage depth 2) # of reads supporting variant 3) VAF of variant
# E.g. setting any of them to 0 means no variants are thrown away according to that criterion
# You can also skip this step and use: vcf.tumfilt <- vcf.somatic
vcf.tumfilt <- vcf.somatic[filterTumour(vcf.somatic@gt[,tumName], tumName, 9, 2,0.05),]
#vcf.tumfilt <- vcf.somatic[filterTumour(vcf.somatic@gt[,tumName]),]
# Save the processed/filtered file


write.vcf(vcf.tumfilt, file = paste0("/Users/syedaanamfatima/Downloads/Edinburgh_IBD_vcfs_2/9N_filtered", colnames(vcf@gt)[i], ".vcf.gz"))


#To split multialleles into separate rows ----------------------------


vcfFileDir <- "/Users/syedaanamfatima/Downloads/Edinburgh_IBD_vcfs_2"
setwd(vcfFileDir)
vcf_file <-"9N_filteredFORMAT.vcf"
vcf_file
vcf<- read.vcfR(vcf_file, verbose=FALSE)


#check number of variants with >1 alternative variant
for (i in 1:length(vcf_file)) {
  print(paste("Currently working on the following file:", vcf_file))
  vcf <- read.vcfR(vcf_file)
  print(paste0("The number of variant rows for ", vcf_file, " is currently ", nrow(vcf@fix)))
  print(paste0("The number of variant rows with >1 alternative variant is ", length(grep(x = vcf@fix[, "ALT"], pattern = ".{2}"))))
  
  # check variant row numbers are consistent between fix and gt regions of vcf
  if (nrow(vcf@fix) == nrow(vcf@gt)) {
    print(paste("The variant row number is the same for fix and gt regions of the file"))
  } else {
    print(paste("The variant row number is NOT the same for fix and gt regions of the file"))
  }
  
  origRowNo <- nrow(vcf@fix)
  altVarRowIndex <- grep(x = vcf@fix[, "ALT"], pattern = ".{2}")
  
  for (j in 1:origRowNo) { # for every variant row
    if (grepl(x = vcf@fix[j, "ALT"], pattern = ".{2}")) { # if it has >1 alt allele
      comNum <- length(grep(x = vcf@fix[j, "ALT"], pattern = ".{2}")) 
      comNum <- as.numeric(comNum)                                   
      splitRowNum <- comNum + 1
      splitRows_fix <- setNames(data.frame(matrix(ncol = length(colnames(vcf@fix)), nrow = splitRowNum)), 
                                colnames(vcf@fix)) # create empty rows for no of alt alleles
      splitRows_gt <- setNames(data.frame(matrix(ncol = length(colnames(vcf@gt)), nrow = splitRowNum)), 
                               colnames(vcf@gt)) # create empty rows for no of alt alleles
      
      # populate empty rows with 1st, 2nd... values of alt col and NRNV data from sample cols
      for (k in 1:splitRowNum) { # for every split row that needs to be populated
        altValues <- unlist(sapply(vcf@fix[j, "ALT"], function(x) strsplit(x, '')))
        splitRows_fix[k, "ALT"] <- altValues[k]
        # populate cols that don't need modification
        splitRows_fix[k, "CHROM"] <- vcf@fix[j, "CHROM"]
        splitRows_fix[k, "POS"] <- vcf@fix[j, "POS"]
        splitRows_fix[k, "ID"] <- vcf@fix[j, "ID"]
        splitRows_fix[k, "REF"] <- vcf@fix[j, "REF"]
        splitRows_fix[k, "QUAL"] <- vcf@fix[j, "QUAL"]
        splitRows_fix[k, "FILTER"] <- vcf@fix[j, "FILTER"]
        splitRows_fix[k, "INFO"] <- vcf@fix[j, "INFO"]
        
        for (l in 2:length(colnames(vcf@gt))) { # for every sample col in every split row that needs to be populated
          GTvalues <- sapply(vcf@gt[j, l], 
                             function(x) unlist(strsplit(x, ':')))
          ADvalues <- tail(GTvalues, 9)
          someGTvalues <- GTvalues[1:4]
          # separate out the ADR, ADVV values
          ADRvalues <- ADvalues[1]
          ADVvalues <- ADvalues[2]
          
          if (grepl(",{1,}", ADRvalues)) {
            # split the strings by comma
            ADRvalues <- sapply(ADRvalues, 
                                function(x) unlist(strsplit(x, ','))) #finish her substitute ADRvalues with AD
            ADVvalues <- sapply(ADVvalues, 
                                function(x) unlist(strsplit(x, ',')))
            # for loop to create GT strings for all alt alleles in row m/k, col l
            for (m in 1:splitRowNum) {
              tempGTstring <- paste(c(someGTvalues, ADRvalues[m], ADVvalues[m]), collapse = ':')
              splitRows_gt[m, l] <- tempGTstring
              # populate FORMAT col
              splitRows_gt[m, "FORMAT"] <- vcf@gt[j, "FORMAT"]
            }
          } else {
            # for loop to create GT strings for all alt alleles in row m/k, col l
            for (m in 1:splitRowNum) {
              tempGTstring <- paste(c(someGTvalues, ADRvalues, ADVvalues), collapse = ':')
              splitRows_gt[m, l] <- tempGTstring
              # populate FORMAT col
              splitRows_gt[m, "FORMAT"] <- vcf@gt[j, "FORMAT"]
            }
          }
          
          
        }
        
      }
      
      # append using rbind the splitRows_fix and splitRows_gt objects to the vcf
      splitRows_fix <- as.matrix(splitRows_fix)
      vcf@fix <- rbind(vcf@fix[, colnames(splitRows_fix)], splitRows_fix)
      splitRows_gt <- as.matrix(splitRows_gt)
      vcf@gt <- rbind(vcf@gt[, colnames(splitRows_gt)], splitRows_gt)
    } else {
      next
    }
  }
  # print number of rows now that alt alleles have been separated into individual rows and added
  print(paste0("After adding split rows, the number of variant rows for ", vcf_file, " is now ", nrow(vcf@fix)))
  print(paste0("After adding split rows, the number of variant rows with >1 alternative variant is now ", length(grep(x = vcf@fix[, "ALT"], pattern = ".{2}"))))
  # delete orig variant rows with multiple alt alleles from fix and gt regions of vcf
  vcf@fix <- vcf@fix[-altVarRowIndex, ]
  vcf@gt <- vcf@gt[-altVarRowIndex, ]
  # print number of rows now that original alt allele rows have been removed
  print(paste0("After removing original >1 alt variant rows, the number of variant rows for ", vcf_file, " is now ", nrow(vcf@fix)))
  print(paste0("After removing original >1 alt variant rows, the number of variant rows with >1 alternative variant is now ", length(grep(x = vcf@fix[, "ALT"], pattern = ".{2}"))))
  # check variant row numbers are consistent between fix and gt regions of vcf
  if (nrow(vcf@fix) == nrow(vcf@gt)) {
    print(paste("After processing, the variant row number is the same for fix and gt regions of the file"))
  } else {
    print(paste("After processing, the variant row number is NOT the same for fix and gt regions of the file"))
  }
  
  # save split row vcf so that original PASSfiltered VCFs are overwritten
  #save<- '/Users/syedaanamfatima/Downloads/Edinburgh_IBD_vcfs_2/deep_filteres_vcfs'
  #setwd(save)
  write.vcf(vcf, file = paste0("/Users/syedaanamfatima/Downloads/Edinburgh_IBD_vcfs_2/deep_filtered_vcfs/9N_final",
                               vcf_file,
                               ".gz"))
}

 
