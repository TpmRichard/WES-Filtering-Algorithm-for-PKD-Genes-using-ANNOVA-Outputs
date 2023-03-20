# The following script filters the WES data for variants of interest within specific genes. The filter comes in seven parts:
# Part 1 = Notice on Startup
# Part 2 = Start up Code
# Part 3 = Extend and Clean
# Part 4 = Basic Filters
# Part 5 = User Area Cleanup
# Part 6 = Advanced Filters
# Part 7 = Dumps

#### PART 1: A Note ####

### Assigning Inputs ###

# This top section requires a little user input to ensure the code can open the data frames and knows were to save the files. 

## Inputs ##

S <- ("Sample Name as String")

{# Start Automated Code by running this line. 

#### PART 2: Setting Up ####

print("Setting up.", quote = FALSE)

## Load Libraries ##

library(tidyverse)

## Load Documents ##

NameString1INDEL <- ".GATK.indel.annovar.hg19_multianno.csv"
NameString1SNP <- ".GATK.snp.annovar.hg19_multianno.csv"
NameString2 <- S
NameINDEL <- paste0(S, NameString1INDEL)
NameSNP <- paste0(S, NameString1SNP)
 
# The block above makes loading the sample data easier. Instead of typing a full name string everytime, the sample number can be assigned and the rest of the file name appended by script. This will differ between projects.
 
FPsnp <- NameSNP        # Assigns file Path for SNP data frame.
FPindel <- NameINDEL      # Assigns file Path for INDEL data frame.
DFsnp <- read.csv(FPsnp) # Loads data frame
DFindel <- read.csv(FPindel) # Loads data frame

# The block above loads the data sets. 

## Outputs ##

# Determine file name and output paths to be used in the rest of the script. FN is the file name, FO is the final output path and save location.

DIR1 <- getwd()
x <- "Processed_Folder"
dir.create(x, recursive = TRUE)
DIR2 <- "/Processed_Folder"
FOutput <- paste0(DIR1, DIR2)      # Assigns file path for output files. 

Fn_cleanSNP <- '/Clean_SNP.csv'         
Fn_cleanINDEL <- '/Clean_INDEL.csv'
Fn_varSNP <- '/ExonicVariants_SNP.csv'
Fn_varINDEL <-'/ExonicVariants_INDEL.csv'
Fn_rareSNP <- '/RareExonicVariants_SNP.csv'
Fn_rareINDEl <- '/RareExonicVariants_INDEL.csv'

FOcleanSNP <- paste(FOutput,Fn_cleanSNP,sep = "") 
FOcleanINDEL <- paste(FOutput,Fn_cleanINDEL,sep = "") 
FOvarSNP <- paste(FOutput,Fn_varSNP,sep = "") 
FOvarINDEL <- paste(FOutput,Fn_varINDEL,sep = "") 
FOrareSNP <- paste(FOutput,Fn_rareSNP,sep = "") 
FOrareINDEl <- paste(FOutput,Fn_rareINDEl,sep = "") 

## The settings below are altered to allow some functions to work better ##

options(scipen = 999) # Disable scientific notation (e^x) which may break a few columns while the script is running. 

#### PART 3: Extend and Clean ####

### Clean Document and Assign Filters ###

# This part will create new filters [Such as Average Population and Predicted Pathogenic] which will be needed in parts 4 and 6 to filter down the information interested in. 

## Extend and Clean SNP ##


print("Extending and cleaning SNP data.", quote = FALSE)

options(warn = - 1) 
 
# A warning does occur when converting columns from text to number due to the presence of ".", However it still works. To stop the console flooding with errors, warnings are disabled while this fragment runs. 

DFsnpProcessing <- DFsnp %>% 
  rename(MiscSampleInfo = 54) %>%
  mutate(Homozygous = ifelse(grepl("0/1", MiscSampleInfo), "No", "Yes"))
# Renames a Column header that is named after the sample ID to a more uniform name, allowing easier script handling. 

DFsnpProcessing <- DFsnpProcessing %>%
  mutate_at(vars(c(30:32)), funs(as.numeric))
# Converts data in columns 30 - 32 from text to numbers. These columns correspond to population and needs to happen to allow for the next line of code to take place. 

DFsnpProcessing$MeanPopAllData <- rowMeans(DFsnpProcessing[,c('X1000g2015aug_all', 'esp6500siv2_all', 'ExAC_ALL')], na.rm=TRUE)
# Calculates an average of 3 databases which is used as a cutoff.

DFsnpProcessing <- DFsnpProcessing %>%
  mutate(GwasCKDFilter = ifelse(grepl("Chronic kidney disease", gwasCatalog), "Yes", "No"))
# This line creates the column used for the genome wide association study filter. 

DFsnpProcessing <- DFsnpProcessing %>%
  mutate(ClinvarPathogenicFilter = ifelse(grepl('pathogenic', clinvar_20150330) & !grepl('non', clinvar_20150330), "Yes", "No"))
# This line creates the column used in the ClinVar filter. 

DFsnpProcessing <- DFsnpProcessing %>%
  mutate(HGMDPathogenicFilter = ifelse(HGMD_ID_Diseasename != ".", "Yes", "No"))
# This line creates the column used in the HGMD filter. 


DFsnpProcessing <- DFsnpProcessing %>% 
  mutate(SIFT_Used = ifelse(SIFT != ".", 1, 0)) %>%
  mutate(Polyphen2_HVAR_Used = ifelse(Polyphen2_HVAR != ".", 1, 0)) %>%
  mutate(MutationTaster_Used = ifelse(MutationTaster != ".", 1, 0)) %>%
  mutate(LRT_Used = ifelse(LRT != ".", 1, 0)) %>%
  mutate(MutationAssessor_Used = ifelse(MutationAssessor != ".", 1, 0)) %>%
  mutate(FATHMM_Used = ifelse(FATHMM != ".", 1, 0)) %>%
  mutate(SIFT_Passed = ifelse(grepl('D', SIFT), 1, 0)) %>%
  mutate(Polyphen2_HVAR_Passed = ifelse(grepl('D', Polyphen2_HVAR) | grepl('P', Polyphen2_HVAR), 1, 0)) %>%
  mutate(MutationTaster_Passed = ifelse(grepl('A', MutationTaster) | grepl('D', MutationTaster), 1, 0)) %>%
  mutate(LRT_Passed = ifelse(grepl('D', LRT), 1, 0)) %>%
  mutate(MutationAssessor_Passed = ifelse(grepl('H', MutationAssessor) | grepl('M', MutationAssessor), 1, 0)) %>%
  mutate(FATHMM_Passed = ifelse(grepl('D', FATHMM), 1, 0))
# This fragment generates a series of Yes (1) and No (0) which are sumed and then divided in later code lines to create the percentage filter. This filter can be used as a cutoff.

DFsnpProcessing$Patho_Programs_Avaliable <- rowSums(DFsnpProcessing[, c('SIFT_Used', 'Polyphen2_HVAR_Used', 'MutationTaster_Used','LRT_Used', 'MutationAssessor_Used', 'FATHMM_Used')])
DFsnpProcessing$Patho_Programs_Passed <- rowSums(DFsnpProcessing[, c('SIFT_Passed', 'Polyphen2_HVAR_Passed', 'MutationTaster_Passed','LRT_Passed', 'MutationAssessor_Passed', 'FATHMM_Passed')])
DFsnpProcessing$Patho_Programs_PercentagePassed <- (DFsnpProcessing$Patho_Programs_Passed/DFsnpProcessing$Patho_Programs_Avaliable)*100
# These 3 lines calculate the percentage values. 

DFsnpProcessing <- DFsnpProcessing %>%
  mutate(GerpFilterPassed = ifelse(gerp..gt2 == ".", "No", "Yes")) %>%
  mutate_at(vars(c(47:49)), funs(as.numeric)) %>%
  mutate(SiPhy_FilterPassed = ifelse(SiPhy_29way_logOdds >= 2, "Yes", "No")) %>%
  mutate(PhyloP_ConservedPassed = ifelse(phyloP20way_mammalian >= 0.5, "Yes", "No")) %>%
  mutate_at(vars(SiPhy_29way_logOdds), ~replace(., is.na(.), "No Data")) %>%
  mutate_at(vars(phyloP20way_mammalian), ~replace(., is.na(.), "No Data")) %>%
  mutate_at(vars(SiPhy_FilterPassed:PhyloP_ConservedPassed), ~replace(., is.na(.), "No")) %>%
  mutate_at(vars(MeanPopAllData), ~replace(., is.nan(.), "No Data")) %>%
  mutate_at(vars(Patho_Programs_PercentagePassed), ~replace(., is.na(.), "No Data")) %>%
  mutate(CADDFilterPassed = ifelse(CADD == ".", "No", "Yes"))
# This block of code makes a series of Pass filters for the conservation programs. Not used to filter data in this code but can determine if a variant is conserved. 

DFsnpProcessed <- DFsnpProcessing %>%
  select(GeneName,Gene,CHROM,POS,cytoBand,ID,Gencode,Func,GeneDetail,ExonicFunc,AAChange,Homozygous,avsnp147,X1000g2015aug_all,esp6500siv2_all,ExAC_ALL,MeanPopAllData,gwasCatalog,
         GwasCKDFilter,clinvar_20150330,ClinvarPathogenicFilter,HGMD_ID_Diseasename,HGMD_mutation,HGMDPathogenicFilter,SIFT,Polyphen2_HVAR,MutationTaster,LRT,MutationAssessor,FATHMM,
         Patho_Programs_Avaliable, Patho_Programs_PercentagePassed,phyloP20way_mammalian,PhyloP_ConservedPassed,SiPhy_29way_logOdds,SiPhy_FilterPassed,gerp..gt2,GerpFilterPassed,CADD,
		 CADDFilterPassed,wgRna,targetScanS,tfbsConsSites,genomicSuperDups,Repeat,OMIM,GWAS_Pubmed_pValue,GO_BP,GO_CC,GO_MF,KEGG_PATHWAY,REACTOME_PATHWAY,REF,ALT,QUAL,FILTER,MiscSampleInfo)
# This line deletes unneeded (such as columns used in calculations only) and unused (such as specific ethnicity population data) columns to make a cleaner spreadsheet that is easier to follow. 

write.csv(DFsnpProcessed, FOcleanSNP)

# Saves the SNP document. 

print("SNP document has been cleaned and extended, ready for filtering.", quote = FALSE)
 
options(warn = 0) 
 
# Returns warnings back to normal
 
# This runs for SNP only. The code is the same for INDEL data but runs in a separate block. 

## Extend and Clean INDEL ## 

print("Extending and cleaning INDEL data.", quote = FALSE)
   
options(warn = - 1) 
  
# A warning does occur when converting columns from text to number due to the presence of ".", However it still works. To stop the console flooding with errors, warnings are disabled while this fragment runs. 
  
DFindelProcessing <- DFindel %>% 
   rename(MiscSampleInfo = 54) %>%
   mutate(Homozygous = ifelse(grepl("0/1", MiscSampleInfo), "No", "Yes"))
# Renames a Column header that is named after the sample ID to a more uniform name, allowing easier script handling. 

DFindelProcessing <- DFindelProcessing %>%
  mutate_at(vars(c(30:32)), funs(as.numeric))
# Converts data in columns 30 - 32 from text to numbers. These columns correspond to population and needs to happen to allow for the next line of code to take place. 

DFindelProcessing$MeanPopAllData <- rowMeans(DFindelProcessing[,c('X1000g2015aug_all', 'esp6500siv2_all', 'ExAC_ALL')], na.rm=TRUE)
# Calculates an average of 3 databases which is used as a cutoff.

DFindelProcessing <- DFindelProcessing %>%
  mutate(GwasCKDFilter = ifelse(grepl("Chronic kidney disease", gwasCatalog), "Yes", "No"))
# This line creates the column used for the genome wide association study filter. 

DFindelProcessing <- DFindelProcessing %>%
  mutate(ClinvarPathogenicFilter = ifelse(grepl('pathogenic', clinvar_20150330) & !grepl('non', clinvar_20150330), "Yes", "No"))
# This line creates the column used in the ClinVar filter. 

DFindelProcessing <- DFindelProcessing %>%
  mutate(HGMDPathogenicFilter = ifelse(HGMD_ID_Diseasename != ".", "Yes", "No"))
# This line creates the column used in the HGMD filter.


DFindelProcessing <- DFindelProcessing %>% 
  mutate(SIFT_Used = ifelse(SIFT != ".", 1, 0)) %>%
  mutate(Polyphen2_HVAR_Used = ifelse(Polyphen2_HVAR != ".", 1, 0)) %>%
  mutate(MutationTaster_Used = ifelse(MutationTaster != ".", 1, 0)) %>%
  mutate(LRT_Used = ifelse(LRT != ".", 1, 0)) %>%
  mutate(MutationAssessor_Used = ifelse(MutationAssessor != ".", 1, 0)) %>%
  mutate(FATHMM_Used = ifelse(FATHMM != ".", 1, 0)) %>%
  mutate(SIFT_Passed = ifelse(grepl('D', SIFT), 1, 0)) %>%
  mutate(Polyphen2_HVAR_Passed = ifelse(grepl('D', Polyphen2_HVAR) | grepl('P', Polyphen2_HVAR), 1, 0)) %>%
  mutate(MutationTaster_Passed = ifelse(grepl('A', MutationTaster) | grepl('D', MutationTaster), 1, 0)) %>%
  mutate(LRT_Passed = ifelse(grepl('D', LRT), 1, 0)) %>%
  mutate(MutationAssessor_Passed = ifelse(grepl('H', MutationAssessor) | grepl('M', MutationAssessor), 1, 0)) %>%
  mutate(FATHMM_Passed = ifelse(grepl('D', FATHMM), 1, 0))
# This fragment generates a series of Yes (1) and No (0) which are sumed and then divided in later code lines to create the percentage filter. This filter can be used as a cutoff.

DFindelProcessing$Patho_Programs_Avaliable <- rowSums(DFindelProcessing[, c('SIFT_Used', 'Polyphen2_HVAR_Used', 'MutationTaster_Used','LRT_Used', 'MutationAssessor_Used', 'FATHMM_Used')])
DFindelProcessing$Patho_Programs_Passed <- rowSums(DFindelProcessing[, c('SIFT_Passed', 'Polyphen2_HVAR_Passed', 'MutationTaster_Passed','LRT_Passed', 'MutationAssessor_Passed', 'FATHMM_Passed')])
DFindelProcessing$Patho_Programs_PercentagePassed <- (DFindelProcessing$Patho_Programs_Passed/DFindelProcessing$Patho_Programs_Avaliable)*100
# These 3 lines calculate the percentage values.  

DFindelProcessing <- DFindelProcessing %>%
  mutate(GerpFilterPassed = ifelse(gerp..gt2 == ".", "No", "Yes")) %>%
  mutate_at(vars(c(47:49)), funs(as.numeric)) %>%
  mutate(SiPhy_FilterPassed = ifelse(SiPhy_29way_logOdds >= 2, "Yes", "No")) %>%
  mutate(PhyloP_ConservedPassed = ifelse(phyloP20way_mammalian >= 0.5, "Yes", "No")) %>%
  mutate_at(vars(SiPhy_29way_logOdds), ~replace(., is.na(.), "No Data")) %>%
  mutate_at(vars(phyloP20way_mammalian), ~replace(., is.na(.), "No Data")) %>%
  mutate_at(vars(SiPhy_FilterPassed:PhyloP_ConservedPassed), ~replace(., is.na(.), "No")) %>%
  mutate_at(vars(MeanPopAllData), ~replace(., is.nan(.), "No Data")) %>%
  mutate_at(vars(Patho_Programs_PercentagePassed), ~replace(., is.na(.), "No Data")) %>%
  mutate(CADDFilterPassed = ifelse(CADD == ".", "No", "Yes"))
# This block of code makes a series of Pass filters for the conservation programs. Not used to filter data in this code but can determine if a variant is conserved. 

DFindelProcessed <- DFindelProcessing %>%
  select(GeneName,Gene,CHROM,POS,cytoBand,ID,Gencode,Func,GeneDetail,ExonicFunc,AAChange,Homozygous,avsnp147,X1000g2015aug_all,esp6500siv2_all,ExAC_ALL,MeanPopAllData,gwasCatalog,
         GwasCKDFilter,clinvar_20150330,ClinvarPathogenicFilter,HGMD_ID_Diseasename,HGMD_mutation,HGMDPathogenicFilter,SIFT,Polyphen2_HVAR,MutationTaster,LRT,MutationAssessor,FATHMM,
         Patho_Programs_Avaliable,Patho_Programs_PercentagePassed,phyloP20way_mammalian,PhyloP_ConservedPassed,SiPhy_29way_logOdds,SiPhy_FilterPassed,gerp..gt2,GerpFilterPassed,CADD,
		 CADDFilterPassed,wgRna,targetScanS,tfbsConsSites,genomicSuperDups,Repeat,OMIM,GWAS_Pubmed_pValue,GO_BP,GO_CC,GO_MF,KEGG_PATHWAY,REACTOME_PATHWAY,REF,ALT,QUAL,FILTER,MiscSampleInfo)
# This line deletes unneeded (such as columns used in calculations only) and unused (such as specific ethnicity population data) columns to make a cleaner spreadsheet that is easier to follow. 

write.csv(DFindelProcessed, FOcleanINDEL)

# Saves the INDEL document. 

print("INDEL document has been cleaned and extended, ready for filtering.", quote = FALSE)
 
options(warn = 0) 
 
# Returns warnings back to normal
 
print("Cleaning up user area.", quote = FALSE)
 
#### PART 4: Basic Filters ####
 
print("Basic filters started.", quote = FALSE)

SNPev <- filter(DFsnpProcessed, Func == "splicing" | Func == "exonic;splicing" | Func == "exonic")                                                                                     
INDELev <- filter(DFindelProcessed, Func == "splicing" | Func == "exonic;splicing" | Func == "exonic")

# The filter above focuses on the genome annotations of interest. 

print("Basic filters finished.", quote = FALSE)

## Population Filter ##

print("Rare variant filters started.", quote = FALSE)
 
# This filter filters out all variants with a greater than 5% population occurrence. 

options(warn = - 1) 

# A warning does occur when converting columns from text to number due to the presence of ".", However it still works. To stop the console flooding with errors, warnings are disabled while this fragment runs. 

SNPrp <- SNPev %>%   # First block sets the column to number so the number filter can work. 
  mutate_at(vars(c(17)), funs(as.numeric)) %>%
  mutate_at(vars(MeanPopAllData), ~replace(., is.na(.), "0")) %>%
  mutate_at(vars(c(31)), funs(as.numeric))
SNPr <- filter(SNPrp, MeanPopAllData <= 0.05) # Second line filters all data with a less than 5% population occurrence when averaging accross the three databases used (Average is calculated in part 3).  
INDELrp <- INDELev %>%
  mutate_at(vars(c(17)), funs(as.numeric)) %>%
  mutate_at(vars(MeanPopAllData), ~replace(., is.na(.), "0")) %>%
  mutate_at(vars(c(31)), funs(as.numeric))
INDELr <- filter(INDELrp, MeanPopAllData <= 0.05)

print("Rare variant filters finished.", quote = FALSE)

options(warn = 0) # Returns warnings back to normal

write.csv(SNPr, FOrareSNP)
write.csv(INDELr, FOrareINDEl)
write.csv(SNPev, FOvarSNP)
write.csv(INDELev, FOvarINDEL)

#### PART 5: User Area Cleanup ####

print("Cleaning up.", quote = FALSE)

# Just renames documents the user will need

SNPr <- SNPr
INDELr <- INDELr
SNPev <- SNPev
INDELev <- INDELev
SNPc <- DFsnpProcessed
INDELc <- DFindelProcessed
SNPo <- DFsnp
INDELo <- DFindel

### Shutdown Unneeded Parts ###

# This is the part that actually cleans up the user area

rm(NameString1SNP, NameString1INDEL, NameString2, NameSNP, NameINDEL, FPsnp, FPindel)
rm(DIR1, DIR2, x, FOutput, Fn_cleanSNP, Fn_cleanINDEL, FOcleanSNP, FOcleanINDEL, DFindel, DFsnp)
rm(DFindelProcessing, DFindelProcessed, DFsnpProcessing, DFsnpProcessed)
rm(SNPrp, INDELrp, Fn_varSNP, Fn_varINDEL, Fn_rareSNP, Fn_rareINDEl, FOvarSNP, FOvarINDEL, FOrareSNP, FOrareINDEl)

print("Clean up has completed. Ready for variant filtering.", quote = FALSE)
 
#### Part 6: Variant Filters ####

# This Filter will filter the INDEL and SNP data sets for a given Sample (Defined as "S(Number)"). The output will consist of 1 new folder and 6 Sheets. 1 Variant sheet, 1 Rare variant sheet and 1 potential mutation sheet for both SNP and INDEL data. 

## ARPKD Variant Filter ##

# This filter setup is duplicated to allow for the investigation of different gene sets. The main setup is the same but the filter == "gene name" parts are replaced with different sets of genes, such as ADPKD genes. 

### Setup ###
	
print("Setting up ARPKD gene filter.", quote = FALSE)

## Setting File Output Path ## 

NFolder <- "Variants_ARPKD"
dir.create(NFolder, recursive = TRUE)
DIR1 <- getwd()
DIR2 <- "/Variants_ARPKD"
FOutput <- paste0(DIR1, DIR2)
	
# The block above creates the new folder and defines where to save the outputs. 

## Setting Output File Names ## 

FnS1 <- '/ARPKDGenes_Var_SNP.csv'  
FnS2 <- '/ARPKDGenes_RareVar_SNP.csv' 
FnS3 <- '/ARPKDGenes_Mutations_SNP.csv' 
FnI1 <- '/ARPKDGenes_Var_INDEL.csv' 
FnI2 <- '/ARPKDGenes_RareVar_INDEL.csv' 
FnI3 <- '/ARPKDGenes_Mutations_INDEL.csv'

# Defines the names of the outputs. 

FOs1  <- paste(FOutput,FnS1,sep = "")
FOs2  <- paste(FOutput,FnS2,sep = "")
FOs3  <- paste(FOutput,FnS3,sep = "")
FOi1  <- paste(FOutput,FnI1,sep = "")
FOi2  <- paste(FOutput,FnI2,sep = "")
FOi3  <- paste(FOutput,FnI3,sep = "")

# Defines the output paths for each file. 

### Gene Filter ###

print("The ARPKD filters have started.", quote = FALSE)
	
## Variants Filter ##

print("Filtering variants within ARPKD genes.", quote = FALSE)

SNPvar <- SNPev %>% 
    filter(GeneName == "PKHD1" | GeneName == "DZIP1L") 
INDELvar <- INDELev %>%
	filter(GeneName == "PKHD1" | GeneName == "DZIP1L") 

# The above functions filter down the data frame to just PKHD1 and DZIP1L variants.
	
## Rare variants filter ##

print("Filtering rare variants within ARPKD genes.", quote = FALSE)

SNPrare <- SNPr %>% 
    filter(GeneName == "PKHD1" | GeneName == "DZIP1L") 
INDELrare <- INDELr %>% 
    filter(GeneName == "PKHD1" | GeneName == "DZIP1L") 

# The above functions filter down the data frame to just PKHD1 and DZIP1L rare variants.
	
## Mutation Filter ##

print("Filtering potential mutations within ARPKD genes.", quote = FALSE)

options(warn = - 1) 
S0 <- SNPr %>%                                                         
    mutate_at(vars(c(32)), funs(as.numeric)) %>%
    mutate_at(vars(Patho_Programs_PercentagePassed), ~replace(., is.na(.), -1))
options(warn = 0) 

# Converts Patho_Programs_PercentagePassed into numbers as the phrase no data makes this a text column. 
# All references to na that forms after converting to numbers are relisted as -1, to stop them appearing in the filter. 
# Options(Warn) disable warning messages (all it says is it has converted the no data into na), these are renabled after the conversion has finished.

options(warn = - 1) 
I0 <- INDELr %>%                                                         
    mutate_at(vars(c(32)), funs(as.numeric)) %>%
    mutate_at(vars(Patho_Programs_PercentagePassed), ~replace(., is.na(.), -1))
options(warn = 0) 

# Repeat of before for INDEL dataset.

S0 <- S0 %>%
    filter(GeneName == "PKHD1" | GeneName == "DZIP1L") 
S1 <- S0 %>% 
    filter(ClinvarPathogenicFilter == "Yes" | HGMDPathogenicFilter == "Yes")
S2 <- S0 %>%
    filter(Func == "splicing" | Func == "exonic;splicing" | ExonicFunc == "stoploss" | ExonicFunc == "stopgain" | ExonicFunc == "frameshift deletion" | 
           ExonicFunc == "frameshift insertion")
S3 <- S0 %>%
    dplyr::filter(ExonicFunc == "missense SNV" & Patho_Programs_PercentagePassed >= 30)
S4 <- S0 %>% 
    dplyr::filter(ExonicFunc == "missense SNV" & Patho_Programs_Avaliable == 0)

# These are the filters. 
# First the rare variant table is filtered down to the genes of interest. 
# This is then followed by a quick scan of ClinVar and HGMD to see if any of the rare variants are known mutations. 
# This is then followed by simple observations such as mutation types highly likely to cause a disease phenotype (Frameshifts and Truncations).
# Next a simple filter is applied to locate missense mutations that fit our disease causing criteria (Prediction of greater than 30%)
# Finally, missense mutations without any information (high chance of being novel) are also scanned for. 
# This only hits the exonic regions due to earlier filters removing Intronic regions.

I0 <- I0 %>%
    filter(GeneName == "PKHD1" | GeneName == "DZIP1L") 
I1 <- I0 %>% 
    filter(ClinvarPathogenicFilter == "Yes" | HGMDPathogenicFilter == "Yes")
I2 <- I0 %>%
    filter(Func == "splicing" | Func == "exonic;splicing" | ExonicFunc == "stoploss" | ExonicFunc == "stopgain" | ExonicFunc == "frameshift deletion" | 
           ExonicFunc == "frameshift insertion")
I3 <- I0 %>%
    dplyr::filter(ExonicFunc == "missense SNV" & Patho_Programs_PercentagePassed >= 30)
I4 <- I0 %>% 
    dplyr::filter(ExonicFunc == "missense SNV" & Patho_Programs_Avaliable == 0)

# Repeat of before for the INDEL dataset.

SMASTER <- full_join(S1, S2)
SMASTER <- full_join(SMASTER, S3)
SMASTER <- full_join(SMASTER, S4)
SMASTER <- filter(SMASTER, ExonicFunc != "synonymous SNV")
IMASTER <- full_join(I1, I2)
IMASTER <- full_join(IMASTER, I3)
IMASTER <- full_join(IMASTER, I4)
IMASTER <- filter(IMASTER, ExonicFunc != "synonymous SNV")
SMASTER <- SMASTER %>%
	filter(!grepl("non", clinvar_20150330)) 
IMASTER <- IMASTER %>% 
	filter(!grepl("non", clinvar_20150330)) 

# The last block of code merges all of the tables together and removes stray synonymous SNV. Also removes confirmed non pathogenic variants from the list. 

## Saving ##

print("Saving.", quote = FALSE)
	
write.csv(SNPvar, FOs1)
write.csv(INDELvar, FOi1)
write.csv(SNPrare, FOs2)
write.csv(INDELrare, FOi2)	
write.csv(SMASTER, FOs3)
write.csv(IMASTER, FOi3)

# These lines just save the tables using the names defined earlier in the code. 

print("Saving completed.", quote = FALSE)
	
print("Cleaning up.", quote = FALSE)

ARPKDsnpVAR <- SNPvar
ARPKDsnpRARE <- SNPrare
ARPKDsnpMutations <- SMASTER
ARPKDindelVAR <- INDELvar
ARPKDindelRARE <- INDELrare
ARPKDindelMutations <- IMASTER

# Gives better to understand names in global view

rm(DIR1, DIR2, FOutput, FnS1, FnS2, FnS3, FnI1, FnI2, FnI3, FOs1, FOs2, FOs3, FOi1, FOi2, FOi3, NFolder)
rm(SNPvar, INDELvar, SNPrare, INDELrare)
rm(S0, S1, S2, S3, S4, SMASTER, I0, I1, I2, I3, I4, IMASTER)

# rm() shuts down data frames in the global environment no longer being used. 
	
print("Clean up finished.", quote = FALSE)
print("ARPKD filtering task completed.", quote = FALSE)

#### PART 7: Mutation Dump ####

### Setup ###

print("Setting up filters for potential mutation dump.", quote = FALSE)

## Setting File Output Path ## 

NFolder <- "Dump_Mutations"
dir.create(NFolder, recursive = TRUE)
DIR1 <- getwd()
DIR2 <- "/Dump_Mutations"
FOutput <- paste0(DIR1, DIR2)
	
# Creates the new folder and defines were to save the outputs. 

## Setting Output File Names ## 

FnS1 <- '/Dump_Mutations_SNP.csv'
FnI1 <- '/Dump_Mutations_INDEL.csv'
FnS3 <- '/Dump_Predicted Mutations_SNP.csv' 
FnI3 <- '/Dump_Predicted Mutations_INDEL.csv'

# Defines the names of the outputs. 

FOs1  <- paste(FOutput,FnS1,sep = "")
FOi1  <- paste(FOutput,FnI1,sep = "")	
FOs3  <- paste(FOutput,FnS3,sep = "")
FOi3  <- paste(FOutput,FnI3,sep = "")

# Defines the Output paths for each file. 

### Mutation Filter ###

print("Creating mutations dump.", quote = FALSE)

options(warn = - 1) 
S0 <- SNPr %>%                                                         
    mutate_at(vars(c(32)), funs(as.numeric)) %>%
    mutate_at(vars(Patho_Programs_PercentagePassed), ~replace(., is.na(.), -1))
options(warn = 0) 

# Converts Patho_Programs_PercentagePassed into numbers as the phrase no data makes this a text column. 
# All references to na that forms after converting to numbers are relisted as -1, to stop them appearing in the filter. 
# Options(Warn) disable warning messages (all it says is it has converted the no data into na), these are renabled after the conversion has finished.

options(warn = - 1) 
I0 <- INDELr %>%                                                         
    mutate_at(vars(c(32)), funs(as.numeric)) %>%
    mutate_at(vars(Patho_Programs_PercentagePassed), ~replace(., is.na(.), -1))
options(warn = 0) 

# Repeat of before for INDEL dataset.

S0 <- S0 %>%
    filter(GeneName != "PKHD1") %>%
    filter(GeneName != "DZIP1L")
	
S1 <- S0 %>% 
    filter(ClinvarPathogenicFilter == "Yes" | HGMDPathogenicFilter == "Yes")
S2 <- S0 %>%
    filter(Func == "splicing" | Func == "exonic;splicing" | ExonicFunc == "stoploss" | ExonicFunc == "stopgain" | ExonicFunc == "frameshift deletion" | 
           ExonicFunc == "frameshift insertion")
S3 <- S0 %>%
    dplyr::filter(ExonicFunc == "missense SNV" & Patho_Programs_PercentagePassed >= 30)
S4 <- S0 %>% 
    dplyr::filter(ExonicFunc == "missense SNV" & Patho_Programs_Avaliable == 0)

# These are the filters. 
# First the rare variant table is filtered down to the genes of interest. 
# This is then followed by a quick scan of ClinVar and HGMD to see if any of the rare variants are known mutations. 
# This is then followed by simple observations such as mutations types highly likely to cause a disease phenotype (Out of Frameshits and Truncations).
# Next a simple filter is applied to locate missense mutations that fit our disease causing criteria (Prediction of greater than 30%)
# Finally, missense mutations without any information (high chance of being novel) are also scanned for. 
# This only hits the exonic regions due to earlier filters removing Intronic regions.

I0 <- I0 %>%
    filter(GeneName != "PKHD1") %>%
	filter(GeneName != "DZIP1L")
	
I1 <- I0 %>% 
    filter(ClinvarPathogenicFilter == "Yes" | HGMDPathogenicFilter == "Yes")
I2 <- I0 %>%
    filter(Func == "splicing" | Func == "exonic;splicing" | ExonicFunc == "stoploss" | ExonicFunc == "stopgain" | ExonicFunc == "frameshift deletion" | 
           ExonicFunc == "frameshift insertion")
I3 <- I0 %>%
    dplyr::filter(ExonicFunc == "missense SNV" & Patho_Programs_PercentagePassed >= 30)
I4 <- I0 %>% 
    dplyr::filter(ExonicFunc == "missense SNV" & Patho_Programs_Avaliable == 0)

# Repeat of before for the INDEL dataset.

SMASTER <- full_join(S2, S3)
SMASTER <- full_join(SMASTER, S4)
SMASTER <- filter(SMASTER, ExonicFunc != "synonymous SNV")
SMASTER <- SMASTER %>% 
	filter(ClinvarPathogenicFilter != "Yes") %>%
    filter(HGMDPathogenicFilter != "Yes") %>%
    filter(!grepl("non", clinvar_20150330)) 
IMASTER <- full_join(I2, I3)
IMASTER <- full_join(IMASTER, I4)
IMASTER <- filter(IMASTER, ExonicFunc != "synonymous SNV")
IMASTER <- IMASTER %>% 
	filter(ClinvarPathogenicFilter != "Yes") %>%
    filter(HGMDPathogenicFilter != "Yes") %>%
    filter(!grepl("non", clinvar_20150330)) 
		
		
# The last block of code merges all of the tables togeather and removes stray synonymous SNV. Also removes confirmed non pathogenic variants from the list. 
# Removes confirmed from the Master files as they are getting saved into there own file. 

## Saving ##

print("Saving.", quote = FALSE)

write.csv(S1,FOs1)
write.csv(I1,FOi1)
write.csv(SMASTER, FOs3)
write.csv(IMASTER, FOi3)

# These lines just save the tables using the names defined earlier in the code. 

print("Saving completed.", quote = FALSE)
	
print("Cleaning up.", quote = FALSE)
    
MutationsDumpSNP <- S1
MutationsDumpINDEL <- I1
PMutationsDumpSNP <- SMASTER
PMutationsDumpINDEL <- IMASTER

# Gives better to understand names in global view

rm(DIR1, DIR2, FOutput, FnS3, FnI3, FOs3, FOi3, NFolder, FnS1, FnI1, FOs1, FOi1)
rm(S0, S1, S2, S3, S4, SMASTER, I0, I1, I2, I3, I4, IMASTER)

SNPrarevariants <- SNPr
INDELrarevariants <- INDELr
SNPvariants <- SNPev
INDELvariants <- INDELev
SNPprocessed <- SNPc
INDELprocessed <- INDELc
SNPoriginal <- SNPo
INDELoriginal <- INDELo

rm(SNPr, INDELr, SNPev, INDELev, SNPc, INDELc, SNPo, INDELo)

# rm() shuts down data frames in the global environment no longer being used. 
	
print("Clean up finished.", quote = FALSE)
print("Mutations dump completed.", quote = FALSE)
}
