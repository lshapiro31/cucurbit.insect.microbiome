

## To see which otus are in the dataset

cat Top_bacteria_in_Striped_Cucumber_Beetles.txt | awk '{print $2}' | sed 's/"//g' | less
grep -f otus.bugs -F seqs_rep_set_tax_assignments.txt  # Need a version of unix that is not mine...

grep "Cyanobacteria" seqs_rep_set_tax_assignments-bugs.txt > Cyanobacteria.list
grep "Wolbachia" seqs_rep_set_tax_assignments-bugs.txt > Wolbachia.list


#### To make table of abundances

for lin in `cat list.txt`
do
echo $lin
cat $lin | awk '{print $2}' | sed 's/"//g' > lin.otus
grep -f lin.otus -F seqs_rep_set_tax_assignments.txt > $lin.prev.txt
done
