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
