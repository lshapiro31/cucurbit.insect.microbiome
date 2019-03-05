# Running a permanova test using the adonis function in vegan. 
# Tests the hypothesis that the samples have different centroids

################# 
### Peponapis ###
#################

set.seed(1)

insect.peponapis <- subset_samples(insect2, Genus=="Peponapis")
insect.peponapis <- prune_samples(sample_sums(insect.peponapis) > 0, insect.peponapis)
TopNOTUs = names(sort(taxa_sums(insect.peponapis), TRUE)[1:20])
insect.peponapis <- prune_taxa(TopNOTUs, insect.peponapis)
insect.peponapis

# Calculate bray curtis distance matrix
insect_bray <- phyloseq::distance(insect.peponapis, method = "bray")

# make a data frame from the sample_data
sampledf <- data.frame(sample_data(insect.peponapis))

# Adonis test
peponapis.adonis <- adonis(insect_bray ~ State2, data = sampledf)
peponapis.adonis

# Homogeneity of dispersion test
beta <- betadisper(insect_bray, sampledf$State2)
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

insect.acalymma <- subset_samples(insect2, Genus=="Acalymma" & Time_Course=="No")
TopNOTUs = names(sort(taxa_sums(insect.acalymma), TRUE)[1:20])
insect.acalymma <- prune_taxa(TopNOTUs, insect.acalymma) 

# Calculate bray curtis distance matrix
insect_bray <- phyloseq::distance(insect.acalymma, method = "bray")

# make a data frame from the sample_data
sampledf <- data.frame(sample_data(insect.acalymma))

# Adonis test
acalymma.adonis <- adonis(insect_bray ~ State2, data = sampledf)
acalymma.adonis

# Homogeneity of dispersion test
beta <- betadisper(insect_bray, sampledf$State2)
permutest(beta)


############################ 
## Acalymma - Time Course ##
############################


## Subset top 20 otus in Acalymma

insect.timeCourse.acalymma <- subset_samples(insect2, Genus=="Acalymma" & Time_Course=="Yes")
TopNOTUs = names(sort(taxa_sums(insect.timeCourse.acalymma), TRUE)[1:20])
insect.timeCourse.acalymma <- prune_taxa(TopNOTUs, insect.timeCourse.acalymma) 
  

# Calculate bray curtis distance matrix
insect_bray <- phyloseq::distance(insect.timeCourse.acalymma, method = "bray")

# make a data frame from the sample_data
sampledf <- data.frame(sample_data(insect.timeCourse.acalymma))

# Adonis test
acalymma.timeCourse.adonis <- adonis(insect_bray ~ State2, data = sampledf)
acalymma.timeCourse.adonis

# Homogeneity of dispersion test
beta <- betadisper(insect_bray, sampledf$State2)
permutest(beta)

############################################




insect_b <- insectR
insect_b <- subset_samples(insect_b, Time_Course=="Yes")
insect_b <- prune_taxa(taxa_sums(insect_b) > 2, insect_b)

set.seed(1)

# Calculate bray curtis distance matrix
insect_bray <- phyloseq::distance(insect_b, method = "bray")

# make a data frame from the sample_data
sampledf <- data.frame(sample_data(insect_b))

# Adonis test
time.course.results <- adonis(insect_bray ~ State2, data = sampledf)
time.course.results
summary(time.course.results)

# Homogeneity of dispersion test
beta <- betadisper(insect_bray, sampledf$State2)
permutest(beta)

## http://microbiome.github.io/microbiome/Core.html
# https://bioconductor.org/packages/release/bioc/html/microbiome.html