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
scp -r ../cluster_upload cluster:/home/groups/harrisonlab/project_files/magnaporthe/assembly/publication_busco_phylogeny_final
```

Add in the missing MG03 isolate

```bash
cp assembly/external/M.sp./MG03/genome.ctg.fa assembly/publication_busco_phylogeny_final/MG03.fasta
```

```bash
scp ../cluster_upload/* cluster:/home/groups/harrisonlab/project_files/magnaporthe/assembly/NCBI_new_genomes_files/.
scp /Users/armita/Downloads/DsLIZ-digitaria.fasta cluster:/home/groups/harrisonlab/project_files/magnaporthe/assembly/NCBI_new_genomes_files/.
scp /Users/armita/Downloads/VO107-digitaria.fasta cluster:/home/groups/harrisonlab/project_files/magnaporthe/assembly/NCBI_new_genomes_files/.
```

```bash
cd /home/groups/harrisonlab/project_files/magnaporthe
mkdir assembly/final_analysis
for File in $(ls $PWD/assembly/NCBI_new_genomes_files/*.fasta); do
  NewName=$(echo $File | rev | cut -f1 -d '/' | rev | sed 's/_Nuclear//g' | sed 's/-digitaria//g')
  cp -s $File assembly/final_analysis/$NewName
done
for File in $(ls $PWD/assembly/publication_busco_phylogeny_final/*.fasta | grep -v -w -e 'D15s47' -e 'D1s11' -e 'D8s15' -e 'E17' -e 'E1' -e 'E34' -e 'E7' -e 'K1367' -e 'K23123' -e 'K38' -e 'K48115n' -e 'K5' -e 'T11' -e 'T12' -e 'T17' -e 'T5' -e 'U13' -e 'U47'); do
  cp -s $File assembly/final_analysis/.
done

```

## Run BUSCO

Quast and busco were run to assess assembly quality:

```bash
# for Assembly in $(ls assembly/publication_busco_phylogeny_final/*.fasta | grep 'MG03'); do
for Assembly in $(ls assembly/final_analysis/*.fasta | grep -v -e 'DsLIZ' -e 'VO107'); do
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
OutDir=analysis_AA/popgen/publication_busco_phylogeny_final_final/$Organism/$Strain/assembly
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
qsub $ProgDir/sub_quast.sh $Assembly $OutDir
# Busco
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/sordariomyceta_odb9)
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done

for Assembly in $(ls assembly/final_analysis/*.fasta | grep -e 'DsLIZ' -e 'VO107'); do
# Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
Strain=$(basename ${Assembly%.fasta})
# Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism="M.digitaria"
Jobs=$(qstat | grep 'sub_busco' | grep 'qw' | wc -l)
while [ $Jobs -gt 2 ]; do
sleep 30s
# printf "."
Jobs=$(qstat | grep 'sub_busco' | grep 'qw' | wc -l)
done		
printf "\n"
# Quast
OutDir=analysis_AA/popgen/publication_busco_phylogeny_final_final/$Organism/$Strain/assembly
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
for File in $(ls analysis_AA/popgen/publication_busco_phylogeny_final/*/*/assembly/run_*/short_summary_*.txt); do  
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
DsLIZ   3577    51      103     45      3725
VO107   3013    0       449     263     3725
98-06   3623    10      69      33      3725
B51     3608    2       84      33      3725
B71     3628    1       69      28      3725
BdJes16-1       3472    1       194     59      3725
BdMeh16-1       3612    2       81      32      3725
D15s47  3609    0       80      36      3725
D1s11   3609    0       82      34      3725
D8s15   3617    0       74      34      3725
E17     3615    3       75      35      3725
E1      3613    0       80      32      3725
E34     3568    9       109     48      3725
E7      3618    4       77      30      3725
EI9064  3622    0       69      34      3725
EI9411  3619    3       74      32      3725
FJ81278 3635    2       42      48      3725
G22     3611    2       82      32      3725
Guy11   3620    1       75      30      3725
K1367   3615    1       77      33      3725
K23123  3626    3       70      29      3725
K38     3611    0       80      34      3725
K48115n 3620    4       73      32      3725
K5      3613    3       80      32      3725
KJ201   3498    9       157     70      3725
MG03    3457    2       166     102     3725
MG12    3558    2       128     39      3725
MZ5-1-6 3628    2       67      30      3725
NI907   3633    2       64      28      3725
PH42    3550    1       123     52      3725
Py22.1  3628    1       65      32      3725
T11     3598    1       89      38      3725
T12     3612    1       83      30      3725
T17     3609    0       83      33      3725
T5      3611    0       80      34      3725
U13     3628    2       66      31      3725
U47     3624    2       67      34      3725
WBKY11  3625    1       72      28      3725
Y34     3602    3       86      37      3725
```


# 1 Find single copy busco genes in Magnaporthe assemblies



Create a list of all BUSCO IDs

```bash
OutDir=analysis_AA/popgen/publication_busco_phylogeny_final_final/tree
mkdir -p $OutDir
BuscoDb="sordariomyceta_odb9"
ls -1 /home/groups/harrisonlab/dbBusco/$BuscoDb/hmms/*hmm | rev | cut -f1 -d '/' | rev | sed -e 's/.hmm//' > $OutDir/all_buscos_"$BuscoDb".txt
```

For each busco gene create a folder and move all single copy busco hits from
each assembly to the folder.
Then create a fasta file containing all the aligned reads for each busco gene for
alignment later.

```bash
printf "" > analysis_AA/popgen/publication_busco_phylogeny_final_final/tree/single_hits.txt
for Busco in $(cat analysis_AA/popgen/publication_busco_phylogeny_final_final/tree/all_buscos_*.txt); do
echo $Busco
OutDir=analysis_AA/popgen/publication_busco_phylogeny_final_final/tree/$Busco
mkdir -p $OutDir
for Fasta in $(ls analysis_AA/popgen/publication_busco_phylogeny_final_final/*/*/assembly/*/single_copy_busco_sequences/$Busco*.fna | grep -v -w  -e 'NC'); do
Strain=$(echo $Fasta | rev | cut -f5 -d '/' | rev)
Organism=$(echo $Fasta | rev | cut -f6 -d '/' | rev)
FileName=$(basename $Fasta)
cat $Fasta | sed "s/:.*.fasta:/:"$Strain":/g" | sed "s/:.*.fa:/:"$Strain":/g" > $OutDir/"$Organism"_"$Strain"_"$Busco".fasta
done
cat $OutDir/*_*_"$Busco".fasta > $OutDir/"$Busco"_appended.fasta
SingleBuscoNum=$(cat $OutDir/"$Busco"_appended.fasta | grep '>' | cut -f2 -d ':' | sort | uniq | wc -l)
printf "$Busco\t$SingleBuscoNum\n" >> analysis_AA/popgen/publication_busco_phylogeny_final_final/tree/single_hits.txt
done
```

If all isolates have a single copy of a busco gene, move the appended fasta to
a new folder

```bash
  OutDir=analysis_AA/popgen/publication_busco_phylogeny_final_final/alignments
  mkdir -p $OutDir
  OrganismNum=$(cat analysis_AA/popgen/publication_busco_phylogeny_final_final/tree/single_hits.txt | cut -f2 | sort -nr | head -n1)
  for Busco in $(cat analysis_AA/popgen/publication_busco_phylogeny_final_final/tree/all_buscos_*.txt); do
  echo $Busco
  HitNum=$(cat analysis_AA/popgen/publication_busco_phylogeny_final_final/tree/single_hits.txt | grep "$Busco" | cut -f2)
  printf "$HitNum\t$OrganismNum\n"
  if [ $HitNum == $OrganismNum ]; then
    # cp analysis_AA/popgen/busco_phylogeny2/$Busco/"$Busco"_appended.fasta $OutDir/.
    cat analysis_AA/popgen/publication_busco_phylogeny_final_final/tree/$Busco/"$Busco"_appended.fasta \
    | sed "s/$Busco://g" \
    | sed "s/genome.ctg.fa://g" \
    | sed "s/_contigs_unmasked.fa//g" \
    | sed -E "s/:.*//g" \
    | tr '.,:' '_' \
    > $OutDir/"$Busco"_appended.fasta
  fi
  done

  ls $OutDir | wc -l
```

```
  2491
```

Submit alignment for single copy busco genes with a hit in each organism


```bash
  AlignDir=analysis_AA/popgen/publication_busco_phylogeny_final_final/alignments
  CurDir=$PWD
  cd $AlignDir
  ProgDir=/home/armita/git_repos/emr_repos/scripts/popgen/phylogenetics
  qsub $ProgDir/sub_mafft_alignment.sh
  cd $CurDir
```

Trimming sequence alignments using Trim-Al
* Note - automated1 mode is optimised for ML tree reconstruction

```bash
  OutDir=analysis_AA/popgen/publication_busco_phylogeny_final/trimmed_alignments
  mkdir -p $OutDir
  for Alignment in $(ls analysis_AA/popgen/publication_busco_phylogeny_final/alignments/*_appended_aligned.fasta); do
    TrimmedName=$(basename $Alignment .fasta)"_trimmed.fasta"
    echo $Alignment
    trimal -in $Alignment -out $OutDir/$TrimmedName -automated1
  done
```

```bash
for Alignment in $(ls analysis_AA/popgen/publication_busco_phylogeny_final/trimmed_alignments/*aligned_trimmed.fasta); do
Jobs=$(qstat | grep 'sub_RAxML' | grep 'qw' | wc -l)
while [ $Jobs -gt 2 ]; do
sleep 2s
# printf "."
Jobs=$(qstat | grep 'sub_RAxML' | grep 'qw' | wc -l)
done		
printf "\n"
echo $Prefix
Prefix=$(basename $Alignment | cut -f1 -d '_')
OutDir=analysis_AA/popgen/publication_busco_phylogeny_final/RAxML/$Prefix
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/phylogenetics
qsub $ProgDir/sub_RAxML.sh $Alignment $Prefix $OutDir
done
```


Run Astral to build a consensus phylogeny from a collective set of
"best phylogenies" from each BUSCO locus

* Note - "Recent versions of ASTRAL output a branch support value even without bootstrapping. Our analyses have revealed that this form of support is more reliable than bootstrapping (under the conditions we explored). Nevertheless, you may want to run bootstrapping as well."

Tutorial tips:
https://github.com/smirarab/ASTRAL/blob/master/astral-tutorial.md#running-with-unresolved-gene-trees

```
Branch length and support

ASTRAL measures branch length in coalescent units and also has a fast way of measuring support without a need for bootstrapping. The algorithms to compute branch lengths and support and the meaning of support outputted is further described in this paper. We will return to these in later sections. Some points have to be emphasized:

ASTRAL only estimates branch lengths for internal branches and those terminal branches that correspond to species with more than one individuals sampled.
Branch lengths are in coalescent units and are a direct measure of the amount of discordance in the gene trees. As such, they are prone to underestimation because of statistical noise in gene tree estimation.
Branch support values measure the support for a quadripartition (the four clusters around a branch) and not the bipartition, as is commonly done.
```

```bash
OutDir=analysis_AA/popgen/publication_busco_phylogeny_final/ASTRAL
mkdir -p $OutDir
# cat analysis_AA/popgen/busco_phylogeny2/RAxML/*/RAxML_bestTree.* > $OutDir/Pcac_phylogeny.appended.tre
# Taxa names were noted to be incorrect at this point and were corrected
cat analysis_AA/popgen/publication_busco_phylogeny_final/RAxML/*/RAxML_bestTree.*  | sed -r "s/CTG.\w+:/:/g" | sed 's/__/_/g' | sed -r "s/M_\w*_//g" > $OutDir/Mg_phylogeny.appended.tre
# InTree=$(ls /home/armita/prog/Astral/Astral/test_data/song_primates.424.gene.tre)
# -
# Trimm back branches that have less than 10% bootstrap support for each tree
# in the given file
# -
/home/armita/prog/newick_utilities/newick_utils/src/nw_ed $OutDir/Mg_phylogeny.appended.tre 'i & b<=10' o > $OutDir/Mg_phylogeny.appended.trimmed.tre
# -
# Calculate combined tree
# -
ProgDir=/home/armita/prog/Astral/ASTRAL
# java -Xmx1000M -jar $ProgDir/astral.5.6.1.jar -i $OutDir/Pcac_phylogeny.appended.trimmed.tre -o $OutDir/Pcac_phylogeny.consensus.tre | tee 2> $OutDir/Pcac_phylogeny.consensus.log
java -Xmx1000M -jar $ProgDir/astral.5.6.3.jar -i $OutDir/Mg_phylogeny.appended.tre -o $OutDir/Mg_phylogeny.consensus.tre | tee 2> $OutDir/Mg_phylogeny.consensus.log
```

```
Total Number of elements weighted: 243493
Normalized score (portion of input quartet trees satisfied before correcting for multiple individuals): 0.49759591612864623
Optimization score: 75225226
Optimal tree inferred in 550.1 secs.
(BdMeh16-1,(BdJes16-1,(B71,((Py22_1,WBKY11),(((NI907,(VO107,DsLIZ)),(Guy11,((Y34,KJ201),(FJ81278,98-06)))),(((MG12,MG03),(EI9411,(MZ5-1-6,G22))),(((K48115n,T12),(E17,(E7,E34))),((PH42,EI9064),(((T17,(T5,T11)),(U47,(U13,K5))),((B51,K1367),(E1,(K38,(D8s15,(K23123,(D1s11,D15s47)))))))))))))));
Final quartet score is: 75225226
Final normalized quartet score is: 0.49759591612864623
(NI907,((VO107,DsLIZ)1:0.16357504788118793,((Guy11,((Y34,KJ201)1:0.06708039574542396,(FJ81278,98-06)0.98:0.04246437410410593)1:0.11123798299094884)1:1.3802870766228756,((((MG12,MG03)1:0.10969520866115125,(EI9411,(MZ5-1-6,G22)1:0.08889163484194768)0.96:0.03801186728590763)1:0.20339670700860193,(((K48115n,T12)0.93:0.033195336373058874,(E17,(E7,E34)1:0.07621906524484112)0.66:0.021962380397008686)1:0.13729362897751055,((PH42,EI9064)1:0.11113802322248813,(((T17,(T5,T11)1:0.10318275859015538)0.96:0.038524924818800994,(U47,(U13,K5)0.75:0.023080151296769033)0.73:0.019450719278929735)0.42:0.004628694811233028,((B51,K1367)0.4:0.004060223353179192,(E1,(K38,(D8s15,(K23123,(D1s11,D15s47)0.5:0.008738449925658108)0.65:0.015772771599607254)0.59:0.013160785444352929)0.36:0.001162373846389581)0.48:0.008046812956183981)0.54:0.010739444913158976)0.44:0.0059545058293864256)1:0.1377811056111765)0.59:0.03409910386834665)1:0.18127549492548786,((Py22_1,WBKY11)0.82:0.04706317532157305,(B71,(BdJes16-1,BdMeh16-1)0.68:0.016954593268949113)1:0.18408397642257737)1:0.29671769580317814)1:0.38792070570654474)1:6.502448053146836));
Weight calculation took 467.285112287 secs
ASTRAL finished in 551.844 secs
```

GGtree was used to make a plot.

* Note- Tips can be found here: https://bioconnector.org/r-ggtree.html

* Note- Tips can be found here: https://bioconnector.org/r-ggtree.html

The consensus tree was downloaded to my local machine

* Note - I had to import into geneious and export again in newick format to get around polytomy branches having no branch length.
* Terminal branch lengths are meanlingless from ASTRAL and should all be set to an arbitrary value. This will be done by geneious (set to 1), but it also introduces a branch length of 2 for one isolate that needs to be corrected with sed

## Rooted style:

As a distance of 1.0 was added to the NI907 branch, this was removed in 2x0.05 edits

<!-- ```bash
cat /Users/armita/Downloads/Mg_publication_phylogeny/Mg_phylogeny.consensus.geneious3.tre | sed 's/:2.0/:1.0/g' | sed 's/:1.0/:0.0/g' | sed -E "s/0.99+8/0.00/g" | sed -E "s/0.99+6/0.00/g" > /Users/armita/Downloads/Mg_publication_phylogeny/Mg_phylogeny.consensus.geneious3_edited.tre
``` -->

On local computer

```bash
cd /Users/armita/Downloads/Mg
codsworth:Mg armita$ scp -r cluster:/home/groups/harrisonlab/project_files/magnaporthe/analysis_AA/popgen/publication_busco_phylogeny_final/ASTRAL .
mkdir /Users/armita/Downloads/Mg/ASTRAL/Mg_publication_phylogeny
```

```bash
cat /Users/armita/Downloads/Mg/ASTRAL/Mg_phylogeny.consensus.geneious.tre | sed 's/:2.0/:1.0/g' | sed 's/:1.0/:0.0/g' | sed -E "s/0.99+8/0.00/g" | sed -E "s/0.99+6/0.00/g" | sed 's/0.0000000000000009/0.0/g' | sed 's/0.9999999999999991/0.0/g'  > /Users/armita/Downloads/Mg/ASTRAL/Mg_phylogeny.consensus.geneious_edits.tre
```

<!--
 ```r
# setwd("/Users/armita/Downloads/Mg_publication_phylogeny")
setwd("/Users/armita/Downloads/Mg/ASTRAL/Mg_publication_phylogeny")
#===============================================================================
#       Load libraries
#===============================================================================

library(ape)
library(ggplot2)

library(ggtree)
library(phangorn)
library(treeio)

t <- read.tree("/Users/armita/Downloads/Mg/ASTRAL/Mg_phylogeny.consensus.tre")
layout(matrix(1:4, 2, 2))
ape::plot.phylo(t, type = "unrooted")
ape::plot.phylo(t, type = "unrooted", use.edge.length = FALSE)
ape::plot.phylo(t, type = "fan")
ape::plot.phylo(t, type = "fan", use.edge.length = FALSE)
layout(matrix(1))

# p1 <- ggtree(t, layout="unrooted")
p2 <- ggtree(t, layout="unrooted", branch.length="none")
# p3 <- ggtree(t, layout="circular")
p4 <- ggtree(t, layout="circular", branch.length="none")
# multiplot(p1, p2, p3, p4, ncol = 2)
multiplot(p2, p4, ncol = 2)

ggtree(t, layout="circular") + geom_tiplab2(aes(angle=angle), color='blue')
```
-->


```R
# setwd("/Users/armita/Downloads/Mg_publication_phylogeny")
setwd("/Users/armita/Downloads/Mg/ASTRAL/Mg_publication_phylogeny")
#===============================================================================
#       Load libraries
#===============================================================================

library(ape)
library(ggplot2)

library(ggtree)
library(phangorn)
library(treeio)


tree <- read.tree("/Users/armita/Downloads/Mg/ASTRAL/Mg_phylogeny.consensus.geneious_edits.tre")
# tree <- read.tree("/Users/armita/Downloads/Mg_publication_phylogeny/Mg_phylogeny.consensus_rooted_edited.tre")

mydata <- read.csv("/Users/armita/Downloads/Mg/ASTRAL/Andrew_traits_for_the_MO_project.csv", stringsAsFactors=FALSE,  na.strings=c("","NA"))
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
#t <- t + geom_tiplab(data=tips, aes(color="black"), size=3, hjust=0, offset = +0.05)
t <- t + geom_tiplab(data=tips, aes(label=tips$X.The.names.to.be.displayed.on.the.phylogeny), align=T, linetype = 'dotted', size=3, hjust=0, offset = +0.05)
t <- t + geom_tiplab(data=tips, aes(label=tips$Country), align=T, linetype = NULL, size=3, hjust=0, offset = +0.9)
t <- t + geom_tiplab(data=tips, aes(label=tips$Host), align=T, linetype = NULL, size=3, hjust=0, offset = +1.8)

#+
#scale_color_manual(values=c("gray39","black")) # colours as defined by col2rgb

# Add in a further set of labels that can also be used to extend the right hand plot area
#tree_mod <- data.frame(t$data)
#tree_mod$label <- paste(tips$Country, tips$Host, sep = '   ')
t <- t + geom_tiplab(data=tips, aes(label="", color="black"), size=3, hjust=0, offset = +0.75) +
scale_color_manual(values=c("gray39","black"))

# t <- t + geom_tippoint(data=tips, aes(shape=Grasshopper..Grh.), size=2)

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

# Use this to add an invisable line on to the image to help extend the page if text doesnt fit
t <- t + geom_cladelabel(node=38, label='', colour='NA', offset=+4.25)


# Annotate a clade with a bar line
# t <- t + geom_cladelabel(node=43, label='sect. Alternaria', align=T, colour='black', offset=-0.0)
# t <- t + geom_cladelabel(node=70, label='gaisen clade', align=T, colour='black', offset=-2.0)
# t <- t + geom_cladelabel(node=46, label='tenuissima clade', align=T, colour='black', offset=-2.0)
# t <- t + geom_cladelabel(node=65, label='arborescens clade', align=T, colour='black', offset=-2.0)

# Save as PDF and force a 'huge' size plot
ggsave("Mg_rooted_10-10-19.pdf", t, width =20, height = 30, units = "cm", limitsize = FALSE)
````

<!--
## Unrooted sstyle:

```bash
cat /Users/armita/Downloads/Mg_publication_phylogeny/Mg_phylogeny.consensus_unrooted.tre | sed 's/:1.0/:0.1/g' | sed -E "s/0.99+6/0.00/g" > /Users/armita/Downloads/Mg_publication_phylogeny/Mg_phylogeny.consensus_unrooted_edited.tre
```


```r
setwd("/Users/armita/Downloads/Mg_publication_phylogeny")
#===============================================================================
#       Load libraries
#===============================================================================

library(ape)
library(ggplot2)

library(ggtree)
library(phangorn)
library(treeio)

tree <- read.tree("/Users/armita/Downloads/Mg_publication_phylogeny/Mg_phylogeny.consensus_unrooted_edited.tre")
# tree <- read.tree("/Users/armita/Downloads/Mg_publication_phylogeny/Mg_phylogeny.consensus_unrooted.tre")

mydata <- read.csv("/Users/armita/Downloads/Mg_publication_phylogeny/traits.csv", stringsAsFactors=FALSE,  na.strings=c("","NA"))
rownames(mydata) <- mydata$"Names.of.isolates.as.they.appear.in.assembly.files"
mydata <- mydata[match(tree$tip.label,rownames(mydata)),]

t <- ggtree(tree, aes(linetype=nodes$support), layout="unrooted") # Core tree
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
#t <- t + geom_tiplab(data=tips, aes(color="black"), size=3, hjust=0, offset = +0.05)
t <- t + geom_tiplab(data=tips, aes(label=tips$X.The.names.to.be.displayed.on.the.phylogeny), align=F, linetype = 'dotted', size=3, hjust=0)
#t <- t + geom_tiplab(data=tips, aes(label=tips$Country), align=F, linetype = NULL, size=3, hjust=0)
#t <- t + geom_tiplab(data=tips, aes(label=tips$Host), align=F, linetype = NULL, size=3, hjust=0)

#+
#scale_color_manual(values=c("gray39","black")) # colours as defined by col2rgb

# Add in a further set of labels that can also be used to extend the right hand plot area
#tree_mod <- data.frame(t$data)
#tree_mod$label <- paste(tips$Country, tips$Host, sep = '   ')
t <- t + geom_tiplab2(data=tips, aes(label="", color="black"), size=3, hjust=0, offset = +0.75) +
scale_color_manual(values=c("gray39","black"))

# t <- t + geom_tippoint(data=tips, aes(shape=Grasshopper..Grh.), size=2)

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
# t <- t + geom_cladelabel(node=38, label='', colour='NA', offset=+1.55)


# Annotate a clade with a bar line
# t <- t + geom_cladelabel(node=43, label='sect. Alternaria', align=T, colour='black', offset=-0.0)
# t <- t + geom_cladelabel(node=70, label='gaisen clade', align=T, colour='black', offset=-2.0)
# t <- t + geom_cladelabel(node=46, label='tenuissima clade', align=T, colour='black', offset=-2.0)
# t <- t + geom_cladelabel(node=65, label='arborescens clade', align=T, colour='black', offset=-2.0)

# Save as PDF and force a 'huge' size plot
ggsave("Mg_unrooted.pdf", t, width =20, height = 20, units = "cm", limitsize = FALSE)

````
-->
