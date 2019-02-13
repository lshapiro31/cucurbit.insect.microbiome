# Here is an example of how to run a permanova test using the 
# adonis function in vegan. 
# In this example we are testing the hypothesis that the 
# three stations we collected samples from have different centroids

set.seed(1)

# Calculate bray curtis distance matrix
insect_bray <- phyloseq::distance(insect_b, method = "bray")

# make a data frame from the sample_data
sampledf <- data.frame(sample_data(insect_b))

# Adonis test
adonis(insect_bray ~ insect_Type, data = sampledf)

# Homogeneity of dispersion test
beta <- betadisper(insect_bray, sampledf$insect_Type)
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