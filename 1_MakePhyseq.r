#Create Physeq

#dependencies
library(tidyverse)
library(phyloseq)

#read data
dat<-read.csv("Data/Qry_Ab_Tha.csv") %>% 
    as_tibble %>%
    filter(Ab!=0) %>%
    mutate(MajorSurvey=paste(Method, Site, Year, Seas, sep="_")) %>%
    mutate(RegionalSurvey=paste(Method, AWNm_Uniq, Year, Seas, sep="_"))



#make phyloseq
    #occurance table
        Occ_table<-dat %>%
            select(MajorSurvey, Latin_NoSubSpp, Ab) %>%
            group_by(MajorSurvey ,Latin_NoSubSpp) %>%
            summarise(Ab=sum(Ab)) %>%
            pivot_wider(names_from = Latin_NoSubSpp, values_from = Ab, values_fill=0) %>%
            ungroup

        otuTable<-Occ_table %>%
            select(!MajorSurvey) %>%
            otu_table(taxa_are_rows=F)

        rownames(otuTable)<-Occ_table$MajorSurvey

    #tax table
        taxonomyTable<-dat %>%
            group_by(Latin_NoSubSpp, Order, Family, EUFG, FMFG) %>%
            summarise %>%
            ungroup

        taxaTable <- taxonomyTable %>%  
            #select(!Latin_NoSubSpp) %>%
            mutate_all(as.factor) %>%
            tax_table

        colnames(taxaTable)<-colnames(taxonomyTable)
        rownames(taxaTable)<-taxonomyTable$Latin_NoSubSpp

    #sampledata
        surveyData<-dat %>%
            group_by(AWNm_Uniq, Site, Method, RegionalSurvey, MajorSurvey, Seas, Year) %>%
            summarise(TotalCatch=sum(Ab)) %>%
            mutate(Decade=paste0(floor(Year/10), "0s"))

        sampleData<-surveyData %>%sample_data
        rownames(sampleData)<-surveyData$MajorSurvey

    #combine
        ps<-phyloseq(otuTable, taxaTable, sampleData)

# save output
    saveRDS(file="Outputs/ps.RDS", ps)