#Some analysis

#dependencies
library(tidyverse)
library(phyloseq)
#ggthemr::ggthemr("fresh")

#load data
ps<-readRDS(file="Outputs/ps.RDS")


#abundance analysis
    #total
        sampleSumData<-sample_sums(ps) %>%
            cbind(sample_data(ps)) 

        TotalAbundance<-sampleSumData %>% 
            ggplot(aes(x=Year, y=log10(.), color=Method)) +
            geom_point()+
            geom_smooth(method=lm)+
            coord_cartesian(ylim=c(0,5))+
            facet_grid(AWNm_Uniq~1)+
            ylab("log10(Abundance)")

    #split by feeding functional group
        psFMFG<-ps
        tax_table(psFMFG)<-cbind(tax_table(ps)[,5], NA)

        psFMFG<-tax_glom(psFMFG, "FMFG")
        tax_table(psFMFG)<-cbind(tax_table(psFMFG)[,1], c("Benthivore", "Omnivore", "Pescivore", "Detritivore", "Zooplanktivore"))

        FMFGAbundance<-psmelt(psFMFG) %>%  
            ggplot(aes(y=log10(Abundance),  color=Method, x=Year)) +
                geom_point()+
                geom_smooth(method=lm)+
                coord_cartesian(ylim=c(0,5))+
                facet_grid(AWNm_Uniq~V2)

    #split by estuary use functional group
        psEUFG<-ps
        tax_table(psEUFG)<-cbind(tax_table(ps)[,4], NA)

        psEUFG<-tax_glom(psEUFG, "EUFG")
        tax_table(psEUFG)<-cbind(tax_table(psEUFG)[,1], c("Coastal", "Marine migrant", "Estuarine", "Freshwater", "Anadromous", "Freshwater migrant", "Marine"))

        EUFGAbundance<-psmelt(psEUFG) %>%  
            ggplot(aes(y=log10(Abundance),  color=Method, x=Year)) +
                geom_point()+
                geom_smooth(method=lm)+
                coord_cartesian(ylim=c(0,5))+
                facet_grid(AWNm_Uniq~V2)

    #split by taxa 

        regionalTaxaAbundances<-list()
        for (RiverRegion in unique(sample_data(ps)$AWNm_Uniq)) {
            psRegional<-ps %>% prune_samples(sample_data(.)$AWNm_Uniq==RiverRegion, .) %>% prune_taxa(taxa_sums(.)>0, .)

            topTax<-sort(taxa_sums(psRegional), TRUE)[which(as.numeric(sort(taxa_sums(psRegional), TRUE))>100)] %>% names


            regionalTaxaAbundances[[RiverRegion]]<-psRegional %>%
                prune_taxa( taxa_names(.) %in% topTax, .) %>%
                psmelt(.) %>%  
                arrange(Latin_NoSubSpp) %>%   # First sort by val. This sort the dataframe but NOT the factor levels
                mutate(Latin_NoSubSpp=factor(Latin_NoSubSpp, levels=unique(Latin_NoSubSpp))) %>% 
                    ggplot(aes(y=log10(Abundance),  color=Method, x=Year)) +
                        geom_point()+
                        geom_smooth(method=lm)+
                        coord_cartesian(ylim=c(0,log10(max(otu_table(psRegional)))))+
                        facet_wrap(~Latin_NoSubSpp)+
                        ggtitle(RiverRegion)
        }

        regionalTaxaAbundances[["THA.L"]]
        regionalTaxaAbundances[["THA.M"]]
        regionalTaxaAbundances[["THA.U"]]




#richness
    #basic data visualisations
    sampleSumData<-sample_sums(ps) %>%
        cbind(sample_data(ps)) 

    #proportion of samples removed if 10 sample depth required for richness estimation
    SmallSamples<-sampleSumData%>% 
        ggplot(aes(x=Method, y=log10(.), color=Method)) +
        geom_boxplot()+
        geom_jitter(width=0.3)+
        geom_hline(yintercept=log10(10))+
        facet_wrap(~AWNm_Uniq)+
            ylab("log10(Abundance)")

    #distribution of those samples across time
    SmallSamplesTemporalDistribution<-sampleSumData %>% 
        filter(.<10) %>%
        ggplot(aes(x=Year, y=log10(.), color=Method)) +
        geom_point(data=sampleSumData, color="grey", alpha=0.5)+
        geom_point()+
        geom_hline(yintercept=log10(10))+
        facet_grid(Method~AWNm_Uniq)+
            ylab("log10(Abundance)")



#inext
    psBigSamples<-ps %>% prune_samples(sample_sums(.)>=10, .) %>% prune_taxa(taxa_sums(.)>0, .) # drops 1 taxa (gobius auratus. mediterranean endemic) and 132 samples (mostly BT15 and KCK, along with some SEN)


    metadata<-sample_data(psBigSamples)
    metadata$Sample<-rownames(metadata)

    otuTable<-otu_table(psBigSamples)
    otuTable<-unclass(otuTable) %>% t %>%as.data.frame
        richness<-iNEXT::estimateD(otuTable,
                        q=c(0,1,2),
                        datatype="abundance",
                        base="coverage",
                        level=NULL,
                        nboot=0
                        )

        head(richness)

    richness<-richness %>%
        left_join(metadata, by=c("Assemblage"="MajorSurvey")) %>% 
        filter(Method.y!="KCK") %>%
        filter(!(Method.y=="SEN" & AWNm_Uniq=="THA.L")) %>%
        filter(!(Method.y=="BT15" & AWNm_Uniq=="THA.M")) 


# New facet label names for supp variable
supp.labs <- c("Richness Estimator", "Shannon Estimator", "Simpson Estimator")
names(supp.labs) <- c("0", "1", "2")


    AlphaTemporalChange<-richness%>%
    ggplot(aes(x=Year, y=qD, color=Method.y))+
            geom_point()+
            geom_smooth(method=lm)+
            facet_grid(Order.q~AWNm_Uniq, , 
            labeller = labeller( Order.q = supp.labs))+
            coord_cartesian(ylim=c(0,10))









# ordination

    psRA<-rarefy_even_depth(ps, replace=F, rngseed=1, sample.size=100) %>%
            prune_samples(!sample_data(.)$Method %in% c("KCK", "BT2"), .) %>% prune_taxa(taxa_sums(.)>0, .) 


        dist<-distance(psRA, "bray")
        ord <- ordinate(psRA, "NMDS", distance=dist , trymax=500)

        CompositionalTemporalChange<-plot_ordination(psRA, ord, type="samples", shape="Seas", color="Decade", title=paste0("Bray-Curtis Dissimilarity, Stress=", signif(ord$stress, 2)))+
            facet_grid(Method~AWNm_Uniq)+
            viridis::scale_color_viridis(discrete=TRUE)


        CompositionalMethodDifferences<-plot_ordination(psRA, ord, type="samples",  color="Method", title=paste0("Bray-Curtis Dissimilarity, Stress=", signif(ord$stress, 2)))+
          stat_ellipse(type = "t")



#printing
    pdf(file = file.path("Figures",paste0("ThamesFishSurveysDataExploration",Sys.Date(),".pdf")), width=20, height =20 ) # The height of the plot in inches        

        plot(0:10, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")
        text(5, 8, "Initial Exploration of Thames Fish Survey Data")
        text(5, 7, "Survey methods, abundance, alpha diversity and compositional turnover.")
        text(5, 6, paste0("Jake Williams   ", Sys.Date()))

            plot(0:10, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")
            text(5, 6, "Difference in Survey Methods")

        TotalAbundance+ labs(
                subtitle="First we visualise the total number of surveys of each type, and their catch rates, through time \n
this reveals systematically different catch rates and hints at different compositions of those catches \n
as the change through time differs by survey method. \n
It also suggests the possibility of abundance trends, which we will come back to."
        )

        CompositionalMethodDifferences+ labs(
                subtitle="To visualise how composition differs by method we create an NMDS plot. First we subsample every \n
sample to the same catch rate (100 in this case), and discard samples with lower catch rates ('rarefying'). This allows \n
the effect of survey method on composition to be isolated from catch rate. Though it does result in losing \n
all samples from the KCK and BT2 survey method. We then calculate the Bray-Curtis dissimlarity between the \n
compositions of every pair of samples. Finally, we reduce this to 2 dimensional space for plotting. \n
The result shows us that the survey methods systematically catch different parts of the fish community. \n
particularly seine netting"
        )

            plot(0:10, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")
            text(5, 6, "Abundance")

        FMFGAbundance+ labs(
                subtitle="Abundance trends for each feeding method functional group within each river region."
        )
        EUFGAbundance+ labs(
                subtitle="Abundance trends for each estuarine use functional group within each river region."
        )
        regionalTaxaAbundances[["THA.L"]]+ labs(
                title="Lower Thames",
                subtitle="Abundance trends for all species for which more than 100 individuals were caught in this river region."
        )
        regionalTaxaAbundances[["THA.M"]]+ labs(
                title="Middle Thames",
                subtitle="Abundance trends for all species for which more than 100 individuals were caught in this river region."
        )
        regionalTaxaAbundances[["THA.U"]]+ labs(
                title="Upper Thames",
                subtitle="Abundance trends for all species for which more than 100 individuals were caught in this river region."
        )
        
            plot(0:10, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")
            text(5, 6, "Alpha Diversity")

            SmallSamples+ labs(
                subtitle="To estimate metrics of alpha diversity (broadly - how many taxa are in a community), we will normalise \n
all samples to equal sample coverage (i.e. a depth that represents the same proportion of all taxa predicted to be in\n
the community). To perform this normalisations samples must include enough catches to make an estimate of coverage, we \n
 therefore need to discard samples with too few catches. We discard all samples with fewer than 10 catches, those falling below the \n
 indicated line"
        )

            SmallSamplesTemporalDistribution+ labs(
                subtitle="To confirm we haven't inadvertently discarded all samples from a particulalr period where catch rates were low \n
we check the distribution of these samples across time (they are the colored point in this plot). We see they are \n
randomly distributed."
        )



    AlphaTemporalChange+ labs(
                subtitle="From the remaining samples we normalsie to even sample coverage and then calculate hill numbers. Hill numbers \n
are a robust metric of alpha diversity that unites richness and eveness in a single framework. For our purposes we \n
calculate hill numbers corresponding to arobust measure of species richness, Shannon diversity (evenness matters more) \n
and Simpson diversity (evenness matters even more). This reveals some trends in the middle and lower thames, but any \n
effect in the upper thames is very small. Note the differet trends in otter trawling and seine netting; this likely \n
results from the different parts of the community they pick up, as seen in the above NMDS plot and taxa abundance \n
breakdown for the middle thames. "
        )

            plot(0:10, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")
            text(5, 6, "Compositional Turnover")

    CompositionalTemporalChange+ labs(
                subtitle="I also include a quick NMDS plot (again based on Bray-Curtis Dissimilarity of rarefied samples) to look at whether \n
composition is systamtically shifting through time. Any effect here must be small given what we can see."
        )


                dev.off()
