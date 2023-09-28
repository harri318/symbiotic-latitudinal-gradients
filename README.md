# README 

## Analysis for "Is there a latitudinal diversity gradient in symbiotic microbes? A case study with sensitive partridge peas"
Tia L. Harrison, Zoe A. Parshuram, Megan E. Frederickson, John R. Stinchcombe 

This repository stores code and data for analyses to be submitted for publication at Molecular Ecology 

## Abstract 

Mutualism is thought to be more prevalent in the tropics than temperate zones and may therefore play an important role in generating and maintaining high species richness found at lower latitudes. However, results on the impact of mutualism on latitudinal diversity gradients are mixed, and few empirical studies sample both temperate and tropical regions. We investigated whether a latitudinal diversity gradient exists in the symbiotic microbial community associated with the legume Chamaecrista nictitans. We sampled bacteria DNA from nodules and the surrounding soil of plant roots across a latitudinal gradient (38.64 °N to 8.68 °N). Using 16S rRNA sequence data, we identified many non-rhizobial species within C. nictitans nodules that cannot form nodules or fix nitrogen. Species richness increased towards lower latitudes in the non-rhizobial portion of the nodule community but not in the rhizobial community. The microbe community in the soil did not effectively predict the non-rhizobia community inside nodules, indicating that host selection is important for structuring non-rhizobia communities in nodules. We next factorially manipulated the presence of three non-rhizobia strains in greenhouse experiments and found that co-inoculations of non-rhizobia strains with rhizobia had a marginal effect on nodule number and no effect on plant growth. Our results suggest that these non-rhizobia bacteria are likely commensals – species that benefit from associating with a host but are neutral for host fitness. Overall, our study suggests that temperate C. nictitans plants are more selective in their associations with the non-rhizobia community, potentially due to differences in soil nitrogen across latitude.

## Description of folders 
1. Coinoculation Experiments

  Data from greenhouse and growth chamber experiments testing the effects of rhizobia and non-rhizobia strains on host plant fitness. R code available for analysis of greenhouse data (nodule number, nodule weight, biomass)

2. Latitude Analyses

  Data includes ASV feature tables, ASV phylogenies, and taxonomy tables created from 16S sequence data using QIIME2. Metadata on samples and populations are also available. R scripts are included for analyzing ASV data and metadata for culture, nodule, and soil samples.

3. QIIME2 Analyses

  Code available for analyzing raw 16S sequence data in Linux. ASV feature tables, phylogenies, and taxonomy files created using this code was used for analysis in R (Latitude Analyses folder) 


