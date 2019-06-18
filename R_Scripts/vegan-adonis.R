# Running a permanova test using the adonis function in vegan. 
# Tests the hypothesis that the samples have different centroids

table(sample_data(acalymma.peponapis)$State)
table(sample_data(acalymma.peponapis)$State2)

################# 
### Peponapis ###
#################

set.seed(1)

TopNOTUs = names(sort(taxa_sums(peponapis), TRUE)[1:20])
peponapis <- prune_taxa(TopNOTUs, peponapis)
peponapis

# Calculate bray curtis distance matrix
peponapis_bray <- phyloseq::distance(peponapis, method = "bray")

# Calculate UniFrac distance matrix
peponapis_unifrac <- phyloseq::distance(peponapis, method = "unifrac")

# make a data frame from the sample_data
peponapis.sampledf <- data.frame(sample_data(peponapis))

# Adonis test
peponapis.adonis <- adonis(peponapis_unifrac ~ State2, 
                           permutations = 99,
                           data = peponapis.sampledf)
peponapis.adonis

# Homogeneity of dispersion test
beta <- betadisper(peponapis_unifrac, peponapis.sampledf$State2)
permutest(beta)

# This output tells us that our adonis test is 
# significant so we can reject the null hypothesis that 
# our three sites have the same centroid.

# Additionally, our betadisper results are not significant, 
# meaning we cannot reject the null hypothesis that our groups 
# have the same dispersions. This means we can be more confident 
# that our adonis result is a real result, and not due to differences 
# in group dispersions
# There is a lot more analysis that can be done here. 
# We could use a distance metric other than Bray-curtis, 
# we could test different grouping variables, or we could 
# create a more complex permanova by testing a model that 
# combines multiple variables. 
# Unfortunately, there are currently no post-hoc tests 
# developed for adonis.

################################ 
## Acalymma - non Time Course ##
################################

## Subset top 20 otus in Acalymma

TopNOTUs = names(sort(taxa_sums(acalymma.non.timecourse), TRUE)[1:20])
acalymma <- prune_taxa(TopNOTUs, acalymma.non.timecourse) 

# Calculate bray curtis distance matrix
acalymma.non.timecourse_bray <- phyloseq::distance(acalymma.non.timecourse, method = "bray")

# Calculate UniFrac distance matrix
acalymma.non.timecourse.unifrac <- phyloseq::distance(acalymma.non.timecourse, method = "unifrac")

# make a data frame from the sample_data
sampledf.acalymma <- data.frame(sample_data(acalymma.non.timecourse))

# Adonis test
acalymma.adonis <- adonis(acalymma.non.timecourse.unifrac ~ State2, 
                   #       by = State2,
                          permutations = 99,
                          data = sampledf.acalymma)
acalymma.adonis

# Homogeneity of dispersion test
beta <- betadisper(acalymma_bray, sampledf.acalymma$State2)
permutest(beta)

############################ 
## Acalymma - Time Course ##
############################

## Subset top 20 otus in Acalymma

TopNOTUs = names(sort(taxa_sums(acalymma.timecourse), TRUE)[1:20])
acalymma.timecourse <- prune_taxa(TopNOTUs, acalymma.timecourse) 

# Calculate bray curtis distance matrix
acalymma.timeCourse_bray <- phyloseq::distance(acalymma.timecourse, method = "bray")

# Calculate UniFrac distance matrix
acalymma.timeCourse.unifrac <- phyloseq::distance(acalymma.timecourse, method = "unifrac")

# make a data frame from the sample_data
acalymma.timeCourse.sampledf <- data.frame(sample_data(acalymma.timecourse))

# Adonis test
acalymma.timeCourse.adonis <- adonis(acalymma.timeCourse.unifrac ~ State, 
                                     permutations = 99,
                                     data = acalymma.timeCourse.sampledf)
acalymma.timeCourse.adonis

# Homogeneity of dispersion test
beta <- betadisper(acalymma.timeCourse_bray, acalymma.timeCourse.sampledf$State)
permutest(beta)

############################################

# https://rstudio-pubs-static.s3.amazonaws.com/268156_d3ea37937f4f4469839ab6fa2c483842.html#alpha-diversity
## http://microbiome.github.io/microbiome/Core.html
# https://bioconductor.org/packages/release/bioc/html/microbiome.html