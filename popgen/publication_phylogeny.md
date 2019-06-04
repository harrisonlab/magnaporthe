# Busco Phylogeny

Commands used to build a consensus maximum likihood phylogeny from conserved
BUSCO genes identified in sequenced isolates.

# 0 prepare data for BUSCO analysis

## Collect genome for analysis

Assembly files were provided on a memory stick and uploaded to the cluster:

```bash
cd /Volumes/TAIWO/Genomes\ file
mkdir ../cluster_upload
for File in $(ls *.fasta); do
NewName=$(echo $File | sed 's/_Nuclear//g' | sed 's/Nuclear//g');
cp $File ../cluster_upload/$NewName;
done
scp -r ../cluster_upload cluster:/home/groups/harrisonlab/project_files/magnaporthe/assembly/publication_busco_phylogeny
```


## Run BUSCO

Quast and busco were run to assess assembly quality:

```bash
for Assembly in $(ls assembly/publication_busco_phylogeny/*.fasta); do
# Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
Strain=$(basename ${Assembly%.fasta})
# Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism="M.grisea"
Jobs=$(qstat | grep 'sub_busco' | grep 'qw' | wc -l)
while [ $Jobs -gt 2 ]; do
sleep 30s
# printf "."
Jobs=$(qstat | grep 'sub_busco' | grep 'qw' | wc -l)
done		
printf "\n"
# Quast
OutDir=analysis_AA/popgen/publication_busco_phylogeny/$Organism/$Strain/assembly
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
qsub $ProgDir/sub_quast.sh $Assembly $OutDir
# Busco
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```

```bash
printf "Strain\tComplete\tDuplicated\tFragmented\tMissing\tTotal\n"
for File in $(ls analysis_AA/popgen/publication_busco_phylogeny/*/*/assembly/run_*/short_summary_*.txt); do  
FileName=$(basename $File)
Strain=$(echo $File | cut -f5 -d '/')
Complete=$(cat $File | grep "(C)" | cut -f2)
Duplicated=$(cat $File | grep "(D)" | cut -f2)
Fragmented=$(cat $File | grep "(F)" | cut -f2)
Missing=$(cat $File | grep "(M)" | cut -f2)
Total=$(cat $File | grep "Total" | cut -f2)
printf "$Strain\t$Complete\t$Duplicated\t$Fragmented\t$Missing\t$Total\n"
done
```
```
98-06   3623    10      69      33      3725
B51     3608    2       84      33      3725
B71     3628    1       69      28      3725
BdJes16-1       3472    1       194     59      3725
BdMeh16-1       3612    2       81      32      3725
D15s47  3612    0       77      36      3725
D1s11   3616    0       76      33      3725
D8s15   3614    0       76      35      3725
E17     3609    1       84      32      3725
E1      3614    0       81      30      3725
E34     3568    9       109     48      3725
E7      3613    3       84      28      3725
EI9064  3622    0       69      34      3725
EI9411  3619    3       74      32      3725
FJ81278 3635    2       42      48      3725
G22     3611    2       82      32      3725
Guy11   3620    1       75      30      3725
K1367   3608    0       86      31      3725
K23123  3641    3       30      54      3725
K38     3603    0       89      33      3725
K48115n 3605    4       85      35      3725
K5      3604    2       85      36      3725
KJ201   3498    9       157     70      3725
MG12    3558    2       128     39      3725
MZ5-1-6 3628    2       67      30      3725
NI907   3633    2       64      28      3725
PH42    3550    1       123     52      3725
Py22.1  3628    1       65      32      3725
T11     3607    0       83      35      3725
T12     3607    0       85      33      3725
T17     3610    0       83      32      3725
T5      3613    0       76      36      3725
U13     3618    2       74      33      3725
U47     3616    0       75      34      3725
WBKY11  3625    1       72      28      3725
Y34     3602    3       86      37      3725
```


# 1 Find single copy busco genes in P.cactorum assemblies



Create a list of all BUSCO IDs

```bash
OutDir=analysis_AA/popgen/publication_busco_phylogeny/tree
mkdir -p $OutDir
BuscoDb="sordariomyceta_odb9"
ls -1 /home/groups/harrisonlab/dbBusco/$BuscoDb/hmms/*hmm | rev | cut -f1 -d '/' | rev | sed -e 's/.hmm//' > $OutDir/all_buscos_"$BuscoDb".txt
```

For each busco gene create a folder and move all single copy busco hits from
each assembly to the folder.
Then create a fasta file containing all the aligned reads for each busco gene for
alignment later.

```bash

printf "" > analysis_AA/popgen/publication_busco_phylogeny/tree/single_hits.txt
for Busco in $(cat analysis_AA/popgen/busco_phylogeny2/all_buscos_*.txt); do
echo $Busco
OutDir=analysis_AA/popgen/publication_busco_phylogeny/tree/$Busco
mkdir -p $OutDir
for Fasta in $(ls analysis_AA/popgen/publication_busco_phylogeny/*/*/assembly/*/single_copy_busco_sequences/$Busco*.fna | grep -v -w  -e 'NC'); do
# for Fasta in $(ls gene_pred/busco/*/*/assembly/*/single_copy_busco_sequences/$Busco*.fna | grep -v -e 'A.gaisen' -e 'dauci' -e 'Alternaria_destruens' -e 'Alternaria_porri' -e 'BMP0308' -e 'BMP2338'); do
Strain=$(echo $Fasta | rev | cut -f5 -d '/' | rev)
Organism=$(echo $Fasta | rev | cut -f6 -d '/' | rev)
FileName=$(basename $Fasta)
cat $Fasta | sed "s/:.*.fasta:/:"$Strain":/g" | sed "s/:.*.fa:/:"$Strain":/g" > $OutDir/"$Organism"_"$Strain"_"$Busco".fasta
done
cat $OutDir/*_*_"$Busco".fasta > $OutDir/"$Busco"_appended.fasta
SingleBuscoNum=$(cat $OutDir/"$Busco"_appended.fasta | grep '>' | cut -f2 -d ':' | sort | uniq | wc -l)
printf "$Busco\t$SingleBuscoNum\n" >> analysis_AA/popgen/publication_busco_phylogeny/tree/single_hits.txt
done
```

If all isolates have a single copy of a busco gene, move the appended fasta to
a new folder

```bash
  OutDir=analysis_AA/popgen/publication_busco_phylogeny/alignments
  mkdir -p $OutDir
  OrganismNum=$(cat analysis_AA/popgen/publication_busco_phylogeny/tree/single_hits.txt | cut -f2 | sort -nr | head -n1)
  for Busco in $(cat analysis_AA/popgen/busco_phylogeny2/all_buscos_*.txt); do
  echo $Busco
  HitNum=$(cat analysis_AA/popgen/publication_busco_phylogeny/tree/single_hits.txt | grep "$Busco" | cut -f2)
  printf "$HitNum\t$OrganismNum\n"
  if [ $HitNum == $OrganismNum ]; then
    # cp analysis_AA/popgen/busco_phylogeny2/$Busco/"$Busco"_appended.fasta $OutDir/.
    cat analysis_AA/popgen/publication_busco_phylogeny/tree/$Busco/"$Busco"_appended.fasta \
    | sed "s/$Busco://g" \
    | sed "s/genome.ctg.fa://g" \
    | sed "s/_contigs_unmasked.fa//g" \
    | sed -E "s/:.*//g" \
    | tr '.,:' '_' \
    > $OutDir/"$Busco"_appended.fasta
  fi
  done
```

Submit alignment for single copy busco genes with a hit in each organism


```bash
  AlignDir=analysis_AA/popgen/publication_busco_phylogeny/alignments
  CurDir=$PWD
  cd $AlignDir
  ProgDir=/home/armita/git_repos/emr_repos/scripts/popgen/phylogenetics
  qsub $ProgDir/sub_mafft_alignment.sh
  cd $CurDir
```

Trimming sequence alignments using Trim-Al
* Note - automated1 mode is optimised for ML tree reconstruction

```bash
  OutDir=analysis_AA/popgen/publication_busco_phylogeny/trimmed_alignments
  mkdir -p $OutDir
  for Alignment in $(ls analysis_AA/popgen/publication_busco_phylogeny/alignments/*_appended_aligned.fasta); do
    TrimmedName=$(basename $Alignment .fasta)"_trimmed.fasta"
    echo $Alignment
    trimal -in $Alignment -out $OutDir/$TrimmedName -automated1
  done
```

```bash
for Alignment in $(ls analysis_AA/popgen/publication_busco_phylogeny/trimmed_alignments/*aligned_trimmed.fasta); do
Jobs=$(qstat | grep 'sub_RAxML' | grep 'qw' | wc -l)
while [ $Jobs -gt 2 ]; do
sleep 2s
# printf "."
Jobs=$(qstat | grep 'sub_RAxML' | grep 'qw' | wc -l)
done		
printf "\n"
echo $Prefix
Prefix=$(basename $Alignment | cut -f1 -d '_')
OutDir=analysis_AA/popgen/publication_busco_phylogeny/RAxML/$Prefix
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/phylogenetics
qsub $ProgDir/sub_RAxML.sh $Alignment $Prefix $OutDir
done
```


Run Astral to build a consensus phylogeny from a collective set of
"best phylogenies" from each BUSCO locus

* Note - "Recent versions of ASTRAL output a branch support value even without bootstrapping. Our analyses have revealed that this form of support is more reliable than bootstrapping (under the conditions we explored). Nevertheless, you may want to run bootstrapping as well."

Tutorial tips:
https://github.com/smirarab/ASTRAL/blob/master/astral-tutorial.md#running-with-unresolved-gene-trees


```bash
OutDir=analysis_AA/popgen/publication_busco_phylogeny/ASTRAL
mkdir -p $OutDir
# cat analysis_AA/popgen/busco_phylogeny2/RAxML/*/RAxML_bestTree.* > $OutDir/Pcac_phylogeny.appended.tre
# Taxa names were noted to be incorrect at this point and were corrected
cat analysis_AA/popgen/publication_busco_phylogeny/RAxML/*/RAxML_bestTree.*  | sed -r "s/CTG.\w+:/:/g" | sed 's/__/_/g' | sed -r "s/M_\w*_//g" > $OutDir/Mg_phylogeny.appended.tre
# InTree=$(ls /home/armita/prog/Astral/Astral/test_data/song_primates.424.gene.tre)
# -
# Trimm back branches that have less than 10% bootstrap support for each tree
# in the given file
# -
/home/armita/prog/newick_utilities/newick_utils/src/nw_ed $OutDir/Mg_phylogeny.appended.tre 'i & b<=10' o > $OutDir/Mg_phylogeny.appended.trimmed.tre
# -
# Calculate combined tree
# -
ProgDir=/home/armita/prog/Astral/Astral
# java -Xmx1000M -jar $ProgDir/astral.5.6.1.jar -i $OutDir/Pcac_phylogeny.appended.trimmed.tre -o $OutDir/Pcac_phylogeny.consensus.tre | tee 2> $OutDir/Pcac_phylogeny.consensus.log
java -Xmx1000M -jar $ProgDir/astral.5.6.1.jar -i $OutDir/Mg_phylogeny.appended.tre -o $OutDir/Mg_phylogeny.consensus.tre | tee 2> $OutDir/Mg_phylogeny.consensus.log
java -Xmx1000M -jar $ProgDir/astral.5.6.1.jar -q $OutDir/Mg_phylogeny.consensus.tre -i $OutDir/Mg_phylogeny.appended.tre -o $OutDir/Mg_phylogeny.consensus.scored.tre 2> $OutDir/Mg_phylogeny.consensus.scored.log
```

```
This is ASTRAL version 5.6.1
Gene trees are treated as unrooted
944 trees read from analysis_AA/popgen/busco_phylogeny2/ASTRAL/Mg_phylogeny.appended.tre
All output trees will be *arbitrarily* rooted at 67
```

GGtree was used to make a plot.

* Note- Tips can be found here: https://bioconnector.org/r-ggtree.html

* Note- Tips can be found here: https://bioconnector.org/r-ggtree.html

The consensus tree was downloaded to my local machine

* Note - I had to import into geneious and export again in newick format to get around polytomy branches having no branch length.
* Terminal branch lengths are meanlingless from ASTRAL and should all be set to an arbitrary value. This will be done by geneious (set to 1), but it also introduces a branch length of 2 for one isolate that needs to be corrected with sed

```bash
cat /Users/armita/Downloads/Mg/Mg_phylogeny.consensus.scored.geneious.tre | sed 's/:2.0/:1.0/g' | sed 's/:1.0/:0.1/g' > Mg_phylogeny.consensus.scored.geneious2.tre
# cat /Users/armita/Downloads/Mg/Mg_phylogeny.consensus.scored.geneious.tre | sed 's/:2.0/:1.0/g' | sed 's/:1.0/:0.1/g' | sed 's/FJ81278:0.9999999999999996/FJ81278:0.1/g'> Mg_phylogeny.consensus.scored.geneious2.tre
```


```r
setwd("/Users/armita/Downloads/Mg")
#===============================================================================
#       Load libraries
#===============================================================================

library(ape)
library(ggplot2)

library(ggtree)
library(phangorn)
library(treeio)

tree <- read.tree("/Users/armita/Downloads/Mg/Mg_phylogeny.consensus.scored.geneious2.tre")

mydata <- read.csv("/Users/armita/Downloads/Mg/Andrew_traits2.csv", stringsAsFactors=FALSE,  na.strings=c("","NA"))
rownames(mydata) <- mydata$"Names.of.isolates.as.they.appear.in.assembly.files"
mydata <- mydata[match(tree$tip.label,rownames(mydata)),]

t <- ggtree(tree, aes(linetype=nodes$support)) # Core tree
# Adjust terminal branch lengths:
branches <- t$data
tree$edge.length[branches$isTip] <- 0.1

#Tree <- branches$branch.length
# rescale_tree(t, branches$branch.length)

t <- t + geom_treescale(offset=-1.0, fontsize = 3) # Add scalebar
# t <- t + xlim(0, 0.025) # Add more space for labels



# Colouring labels by values in another df
t <- t %<+% mydata # Allow colouring of nodes by another df
#t <- t + geom_tiplab(aes(color=Source), size=3, hjust=0) +
scale_color_manual(values=c("gray39","black")) # colours as defined by col2rgb

tips <- data.frame(t$data)
tips$label <- paste(tips$X.The.names.to.be.displayed.on.the.phylogeny, tips$Country, tips$Host,  sep = '   ')
t <- t + geom_tiplab(data=tips, aes(color="black"), size=3, hjust=0, offset = +0.05) +
scale_color_manual(values=c("gray39","black")) # colours as defined by col2rgb

# Add in a further set of labels that can also be used to extend the right hand plot area
#tree_mod <- data.frame(t$data)
#tree_mod$label <- paste(tips$Country, tips$Host, sep = '   ')
t <- t + geom_tiplab(data=tips, aes(label="", color="black"), size=3, hjust=0, offset = +0.75) +
scale_color_manual(values=c("gray39","black"))

t <- t + geom_tippoint(data=tips, aes(shape=Grasshopper..Grh.), size=2)

# Format nodes by values
nodes <- data.frame(t$data)
#nodes <- nodes[!nodes$isTip,]
nodes$label <- as.numeric(nodes[nodes$label,])
as.numeric(nodes$label)
#nodes$label[nodes$label < 0.80] <- ''
nodes$support[nodes$isTip] <- 'supported'
nodes$support[(!nodes$isTip) & (nodes$label > 0.80)] <- 'supported'
nodes$support[(!nodes$isTip) & (nodes$label < 0.80)] <- 'unsupported'
nodes$support[(!nodes$isTip) & (nodes$label == '')] <- 'supported'
t <- t + aes(linetype=nodes$support)
nodes$label[nodes$label > 0.80] <- ''
t <- t + geom_nodelab(data=nodes, size=2, hjust=-0.05) # colours as defined by col2rgb


# Annotate a clade with a bar line
# t <- t + geom_cladelabel(node=43, label='sect. Alternaria', align=T, colour='black', offset=-0.0)
# t <- t + geom_cladelabel(node=70, label='gaisen clade', align=T, colour='black', offset=-2.0)
# t <- t + geom_cladelabel(node=46, label='tenuissima clade', align=T, colour='black', offset=-2.0)
# t <- t + geom_cladelabel(node=65, label='arborescens clade', align=T, colour='black', offset=-2.0)

# Save as PDF and force a 'huge' size plot
ggsave("Mg.pdf", t, width =20, height = 30, units = "cm", limitsize = FALSE)

```` -->
