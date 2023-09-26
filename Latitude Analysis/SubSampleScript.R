# AUTHOR:   Tia Lene Harrison 
# DATE:     25 November 2022

# PURPOSE:  Create a data frame that includes the sub-sampled groups of tropical and subtropical populations 

# INPUT: 
# data = dataframe from the overall phyloseq object 
# phylo_obj = the phyloseq object with all the information 

# OUTPUT: 
# dist_data = one dataframe that is all the values for the different geographic regions with subsampling 

# FUNCTION: sub.sample
sub.sample<- function(data, phylo_obj){
  
  # Code the categories and re attach the sample dataset again 
  data2 <- data %>%
    mutate(region=case_when(latitude < 11 ~ "Tropical", 
                            latitude > 35 ~ "Temperate", 
                            latitude > 25 | latitude <33 ~ "Subtropical"))

sub_data<- data2 %>%
  filter(region=="Subtropical") %>%
  dplyr::select(sample.id)

## Do subsampling and pairwise distances in the same loop for the subtropical group 
## Make the empty vector to store the means 
sub_means<- vector(mode="numeric", length=100)

for(i in 1:100) {
  # Sample the 17 plants 
  sub_data_random<-sub_data %>%
    sample_n(17)
  
  # Prune the phyloseq object 
  sub_obj<-prune_samples(sub_data_random$sample.id, phylo_obj)
  
  # Prune out asvs that are zero across all samples 
  sub_obj<-prune_taxa(taxa_sums(sub_obj)>0, sub_obj)
  
  # Get the tree 
  sub_tree<-phy_tree(sub_obj)
  
  # Calculate distances and save 
  sub_dist<-cophenetic(sub_tree)
  sub_dist2<-as.data.frame(sub_dist)
  
  # Transform into a list 
  sub_list<-data.frame(asv1=rownames(sub_dist)[row(sub_dist)],   asv2=colnames(sub_dist)[col(sub_dist)], values=c(sub_dist))
  
  # Take out NAs and give region label 
  sub_list <- sub_list %>%
    filter(!is.na(values)) %>%
    filter(values>0) %>%
    mutate(region="subtropical")
  
  # Calculate the mean 
  sub_means[i]<-mean(sub_list$values)
}

sub_means_data<-as.data.frame(sub_means)
sub_means_data


# Repeat the process for tropical data  
trop_data<- data2 %>%
  filter(region=="Tropical") %>%
  dplyr::select(sample.id)

# Run the loop 
trop_means<- vector(mode="numeric", length=100)

for(i in 1:100) {
  # Sample the 17 plants 
  trop_data_random<-trop_data %>%
    sample_n(17)
  
  # Prune the phyloseq object 
  trop_obj<-prune_samples(trop_data_random$sample.id, phylo_obj)
  
  # Prune out asvs that are zero across all samples 
  trop_obj<-prune_taxa(taxa_sums(trop_obj)>0, trop_obj)
  
  # Get the tree 
  trop_tree<-phy_tree(trop_obj)
  
  # Calculate distances and save 
  trop_dist<-cophenetic(trop_tree)
  trop_dist2<-as.data.frame(trop_dist)
  
  # Transform into a list 
  trop_list<-data.frame(asv1=rownames(trop_dist)[row(trop_dist)],   asv2=colnames(trop_dist)[col(trop_dist)], values=c(trop_dist))
  
  # Take out NAs and give region label 
  trop_list <- trop_list %>%
    filter(!is.na(values)) %>%
    filter(values>0) %>%
    mutate(region="tropical")
  
  # Calculate the mean 
  trop_means[i]<-mean(trop_list$values)
}

trop_means_data<-as.data.frame(trop_means)
trop_means_data

# Add a column for the region and merge the datasets 
trop_means_data <- trop_means_data %>%
  mutate(region="tropical") %>%
  rename(dist_mean=trop_means)

sub_means_data <- sub_means_data %>%
  mutate(region="subtropical") %>%
  rename(dist_mean=sub_means)

# Get the temperate data and calculate the mean - there will only be one though??
temp_data<- data2 %>%
  filter(region=="Temperate") %>%
  dplyr::select(sample.id)

# Prune the phyloseq object 
temp_obj<-prune_samples(temp_data$sample.id, phylo_obj)

# Prune out asvs that are zero across all samples 
temp_obj<-prune_taxa(taxa_sums(temp_obj)>0, temp_obj)

# Get the tree 
temp_tree<-phy_tree(temp_obj)

# Calculate distances and save 
temp_dist<-cophenetic(temp_tree)
temp_dist2<-as.data.frame(temp_dist)

# Transform into a list 
temp_list<-data.frame(asv1=rownames(temp_dist)[row(temp_dist)],   asv2=colnames(temp_dist)[col(temp_dist)], values=c(temp_dist))

# Take out NAs and give region label 
temp_list <- temp_list %>%
  filter(!is.na(values)) %>%
  filter(values>0) %>%
  mutate(region="temperate")

# Calculate the mean 
temp_means<-mean(trop_list$values)

# Save the data 
temp_means_data<-as.data.frame(temp_means)
temp_means_data <- temp_means_data %>%
  mutate(region="temperate") %>%
  rename(dist_mean=temp_means)

# Combine the data sets for plotting 
dist_data<- do.call("rbind", list(sub_means_data, trop_means_data, temp_means_data))

# Final dataset to return 
return(dist_data)

}