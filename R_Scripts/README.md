#### Pipeline for 16S gut microbiome analysis of beetles (Acalymma vittatum) and bees (Peponapis pruinosa) that are adapted to Cucurbita plants

## Parameters for pre-processing of raw data in qiime
qiime.txt

## Import pre-processed Qiime biome file and metadata into R
## Assign global variables and remove rare OTUs
importBiomMapTre.R

## Create map to show were samples were collected
Sample.map.R

## Calculate alpha diversity
alphaRichness.R

## Create barchart plots to show community composition
barplots.R

## Create PCoA plots to cluster samples based on community similarities
beta-diversity.R

## Use the adonis function in vegan to test whether the centroids of the different sample groups in the PCoA are statistically different
vegan-adonis.R

## erwiniaInAcalymma.R
Subset just the Erwinia tracheiphila OTU to compare abundance and frequency in all Acalymma samples

## Create tree to compare which of the most common OTUs occur in bees and/or beetles
plotTrees.R


# insect.microbiome useful stuff for me

git init
git add README.md
git commit -m "first commit"
git remote add origin https://github.com/lshapiro31/insect.microbiome.git
git push -u origin master
