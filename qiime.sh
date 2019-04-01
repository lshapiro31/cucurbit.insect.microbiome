#!/bin/bash 

#SBATCH -n 5 #Number of cores 
#SBATCH -N 1 #run all on the same node
#SBATCH -t 3000  #Runtime in minutes
#SBATCH --mail-user=lori.r.shapiro@gmail.com
#SBATCH --mail-type=END
#SBATCH -o qiime.out
#SBATCH -e qiime.err
#SBATCH -p serial_requeue,general #Partition to submit to 
#SBATCH --mem-per-cpu=30000 #Memory per cpu in MB (see also --mem)

module load python/2.7.14-fasrc01
module load uchime/4.2.40-fasrc01
#source activate ody
source activate qiime1

# validate_mapping_file.py -m New_mapping_file.txt -o validate_mapping_file_output

# split_libraries_fastq.py -o slout/ -i Undetermined_S0_L001_R1_001.fastq.gz,Undetermined_S0_L001_R2_001.fastq.gz -b Undetermined_S0_L001_I1_001.fastq.gz -m map.txt

# join_paired_ends.py -f Undetermined_S0_L001_R1_001.fastq.gz -r Undetermined_S0_L001_R2_001.fastq.gz -b Undetermined_S0_L001_I1_001.fastq.gz -o joined/

# split_libraries_fastq.py -m New_mapping_file.txt -i joined/fastqjoin.join.fastq -b joined/fastqjoin.join_barcodes.fastq --barcode_type 12 -o demultiplexed/

#################
## OTU picking ##
#################

## Open OTUs
# pick_open_reference_otus.py -o otus_open/ -i demultiplexed/seqs.fna

## Closed OTUs (this takes a long time)
# pick_de_novo_otus.py -i demultiplexed/seqs.fna -o otus

# assign_taxonomy.py -i /otus/rep_set/repr_set_seqs.fasta -m rdp

# make_otu_table.py -i otus/uclust_picked_otus/seqs_otus.txt -t otus/uclust_assigned_taxonomy/seqs_rep_set_tax_assignments.txt -o otus/otu_table.biom -m New_mapping_file.txt

### Filter low abundance OTUs from table
# filter_otus_from_otu_table.py -i otus/otu_table.biom -n 5 -o otus/otu_table_filter.biom

### Add sample metadata to the biome file
# biom add-metadata -i otu_table_filter.biom -o otu_table_filter_meta.biom --sample-metadata-fp New_mapping_file.txt --observation-metadata-fp otus/uclust_assigned_taxonomy/seqs_rep_set_tax_assignments.txt --sc-separated taxonomy --observation-header OTUID,taxonomy,confidence

#####################
## Remove Chimeras ##
#####################

## Remove Chimeras
identify_chimeric_seqs.py -i demultiplexed/seqs.fna -m usearch61 -o usearch_checked_chimeras/ -r otus/rep_set/seqs_rep_set.fasta
#identify_chimeric_seqs.py -i repr_set_seqs.fasta -t taxonomy_assignment.txt -r ref_seq_set.fna -m blast_fragments -o chimeric_seqs_blast.txt


#########################################################
## Add both mapping file and taxonomy to the biom file ##
#########################################################

# biom add-metadata -i otu_table_noChloroMito.biom -o otu_filter_noChloroMito_meta.biom --sample-metadata-fp New_mapping_file.txt --observation-metadata-fp otus/uclust_assigned_taxonomy/seqs_rep_set_tax_assignments.txt --sc-separated taxonomy --observation-header OTUID,taxonomy,confidence

# biom add-metadata -i otus/otu_table_filter.biom -o otu_filter_meta.biom --sample-metadata-fp New_mapping_file.txt --observation-metadata-fp otus/uclust_assigned_taxonomy/seqs_rep_set_tax_assignments.txt --sc-separated taxonomy --observation-header OTUID,taxonomy,confidence

# biom convert -i otus/otu_filter_meta.biom -o table.from_biom.txt --to-tsv --header-key taxonomy







####################
## Make phylogeny ##
####################
# https://groups.google.com/forum/#!topic/qiime-forum/CVkmERqGnUw

# parallel_align_seqs_pynast.py -i otus/rep_set/seqs_rep_set.fasta -o otus/rep_set/rep_set_aligned_fasta
# filter_alignment.py -i otus/rep_set/rep_set_aligned.fasta -o otus/rep_set/rep_set_aligned_pfiltered.fasta 
# make_phylogeny.py -i otus/rep_set/rep_set_aligned_pfiltered.fasta -o otus/rep_set/rep_phylo.tre

# Filter tree for tips to keep or remove
# filter_tree.py -i

# alpha_rarefaction.py -i otus/otu_table_filter.biom -o alpha_rarefaction -m New_mapping_file.txt 

# alpha_diversity.py -f -i otu_filter_noChloroMito_noLow_meta.biom -o alpha_diversity.txt -m chao1,goods_coverage,shannon,PD_whole_tree -t otus/rep_se$

# beta_diversity.py -f -i otu_filter_noChloroMito_noLow_meta.biom -m New_mapping_file.txt -o beta_diversity -t otus/rep_set.tre

# beta_diversity_through_plots.py -f -i otu_filter_noChloroMito_noLow_meta.biom -m New_mapping_file.txt -o bdiv/

# summarize_taxa_through_plots.py -f -i otus/otu_filter_meta.biom -o taxa_summary -m New_mapping_file.txt

# make_otu_table.py -i otus/uclust_picked_otus/seqs_otus.txt -t otus/uclust_assigned_taxonomy/seqs_rep_set_tax_assignments.txt -o otu_table.txt

# assign_taxonomy.py -i repr_set_seqs.fasta -r ref_seq_set.fna -t id_to_taxonomy.txt

# make_2d_plots.py -i beta_diversity/unweighted_unifrac_pc.txt -m New_mapping_file.txt -o beta_diversity/scree_unweighted 

# make_2d_plots.py -i beta_diversity/weighted_unifrac_pc.txt -m New_mapping_file.txt -o beta_diversity/scree_weighted


### Remove chimeras
identify_chimeric_seqs.py -m ChimeraSlayer -i otus/rep_set/rep_set_aligned_fasta/seqs_rep_set_aligned.fasta -a otus/pynast_aligned_seqs/seqs_rep_set_aligned.fasta -o chimeric_seqs.txt

