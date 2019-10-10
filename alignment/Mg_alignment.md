
# 1. Alignment of Mg raw reads vs the each genome


Read alignment was performed by Taiwo and sent on Hard Drive to NIAB EMR.

```bash
# From local machine
scp -r /Volumes/TAIWO3/Files\ for\ Andrew cluster:/data/scratch/armita/magnaporthe/alignment/.
```

<!-- Alignment of reads from a single run:

```bash
for Reference in $(ls repeat_masked/*/*/*/*_contigs_unmasked.fa | grep -e '1166' -e '650' -e '97.0013'); do
for StrainPath in $(ls -d qc_dna/paired/*/*); do
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/spades
Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
F_Read=$(ls $StrainPath/F/*_trim.fq.gz)
R_Read=$(ls $StrainPath/R/*_trim.fq.gz)
for FileF in $(ls qc_rna/paired/*/*/F/*.fastq.gz); do
Jobs=$(qstat | grep 'sub_bwa' | grep 'qw'| wc -l)
while [ $Jobs -gt 1 ]; do
sleep 1m
printf "."
Jobs=$(qstat | grep 'sub_bwa' | grep 'qw'| wc -l)
done
echo $F_Read
echo $R_Read
Prefix="${Organism}_${Strain}"
OutDir=/home/scratch/groups/harrisonlab/alternaria/analysis/genome_alignment/bwa/$Organism/$Strain/vs_${Reference}
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment
qsub $ProgDir/bwa/sub_bwa.sh $Prefix $Reference $F_Read $R_Read $OutDir
done
done
done
``` -->

Identify read coverage over each bp

```bash
  cd /data/scratch/armita/magnaporthe
  for Bam in $(ls alignment/*/*_sorted.bam | tail -n+2); do
    Ref=$(echo $Bam | cut -f2 -d '/' | sed 's/Align//g')
    Strain=$(basename ${Bam%_sorted.bam})
    Organism="M.grisea"
    echo "$Organism - $Strain"
    # OutDir=$(dirname $Bam)
    OutDir="alignment/genome_alignment/vs_${Ref}/$Organism/$Strain/"
    mkdir -p $OutDir
    samtools depth -aa $Bam > $OutDir/${Organism}_${Strain}_vs_${Ref}_depth.tsv
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment/coverage_analysis
    $ProgDir/cov_by_window.py --cov $OutDir/${Organism}_${Strain}_vs_${Ref}_depth.tsv > $OutDir/${Organism}_${Strain}_vs_${Ref}_depth_10kb.tsv
    sed -i "s/$/\t$Strain/g" $OutDir/${Organism}_${Strain}_vs_${Ref}_depth_10kb.tsv
  done
for RefDir in $(ls -d alignment/genome_alignment/vs_*); do
Ref=$(echo $RefDir | rev | cut -f1 -d '/'| rev | sed 's/vs_//g')
echo $Ref
cat $RefDir/*/*/*_vs_${Ref}_depth_10kb.tsv > $RefDir/vs_${Ref}_grouped_depth.tsv
done
```

From my local computer open an R-studio session.

```r
library(readr)
# install.packages("ggplot2")
library(ggplot2)
require(scales)
setwd("~/Downloads/Mg/coverage")

myFun <- function(x) {
c(min = min(x), max = max(x),
mean = mean(x), median = median(x),
std = sd(x))
}

# df_WE34 <- read_delim("~/Downloads/Mg/vs_WE34_grouped_depth.tsv", "\t", escape_double = FALSE, col_names = FALSE, col_types = cols(X4 = col_factor(levels = c("675", "97.0013", "97.0016", "650", "648", "24350", "1082", "1164", "635", "743", "1166", "1177"))), trim_ws = TRUE)
# df_WE34 <- read_delim("~/Downloads/Mg/vs_WE34_grouped_depth.tsv", "\t", escape_double = FALSE, col_names = FALSE, col_types = cols(X4 = col_factor(levels = c("D15s47", "D1s11", "D8s15", "E1", "E17", "E34", "E7", "K1367", "K23123", "K38", "K48115n", "K5", "T11", "T12", "T17", "T5", "U13", "U47"))), trim_ws = TRUE)
df_WE34 <- read.delim("~/Downloads/Mg/vs_WE34_grouped_depth.tsv", header=FALSE)

colnames(df_WE34) <- c("contig","position", "depth", "strain")
# install.packages("readr")
library(readr)
df_WE34$treatment <- paste(df_WE34$strain , df_WE34$contig)
tapply(df_WE34$depth, df_WE34$treatment, myFun)
df2 <- cbind(do.call(rbind, tapply(df_WE34$depth, df_WE34$treatment, myFun)))
write.csv(df2, 'WE34_contig_coverage.csv')
df_WE34$depth <- ifelse(df_WE34$depth > 200, 200, df_WE34$depth)



for (i in 1:9){
contig = paste("Scaffold0", i, sep = "")
p0 <- ggplot(data=df_WE34[df_WE34$contig == contig, ], aes(x=`position`, y=`depth`, group=1)) +
geom_line() +
labs(x = "Position (bp)", y = "Coverage") +
scale_y_continuous(breaks=seq(0,200,50), limits=c(0,200)) +
facet_wrap(~strain, nrow = 18, ncol = 1, strip.position = "left")
outfile = paste("WE34_contig", i, "cov.jpg", sep = "_")
ggsave(outfile , plot = p0, device = 'jpg', path = NULL,
scale = 1, width = 500, height = 500, units = 'mm',
dpi = 150, limitsize = TRUE)
}
for (i in 10:31){
contig = paste("Scaffold", i, sep = "")
p0 <- ggplot(data=df_WE34[df_WE34$contig == contig, ], aes(x=`position`, y=`depth`, group=1)) +
geom_line() +
labs(x = "Position (bp)", y = "Coverage") +
scale_y_continuous(breaks=seq(0,200,50), limits=c(0,200)) +
facet_wrap(~strain, nrow = 18, ncol = 1, strip.position = "left")
outfile = paste("WE34_contig", i, "cov.jpg", sep = "_")
ggsave(outfile , plot = p0, device = 'jpg', path = NULL,
scale = 1, width = 500, height = 500, units = 'mm',
dpi = 150, limitsize = TRUE)
}




df_WK23123 <- read.delim("~/Downloads/Mg/vs_WK23123_grouped_depth.tsv", header=FALSE)

colnames(df_WK23123) <- c("contig","position", "depth", "strain")
# install.packages("readr")
library(readr)
df_WK23123$treatment <- paste(df_WK23123$strain , df_WK23123$contig)
tapply(df_WK23123$depth, df_WK23123$treatment, myFun)
df2 <- cbind(do.call(rbind, tapply(df_WK23123$depth, df_WK23123$treatment, myFun)))
write.csv(df2, 'WK23123_contig_coverage.csv')
df_WK23123$depth <- ifelse(df_WK23123$depth > 200, 200, df_WK23123$depth)



for (contig in levels(df_WK23123$contig)){
# contig = paste("tig0000000", i, "_pilon_pilon_pilon", sep = "")
p0 <- ggplot(data=df_WK23123[df_WK23123$contig == contig, ], aes(x=`position`, y=`depth`, group=1)) +
geom_line() +
labs(x = "Position (bp)", y = "Coverage") +
scale_y_continuous(breaks=seq(0,200,50), limits=c(0,200)) +
facet_wrap(~strain, nrow = 18, ncol = 1, strip.position = "left")
outfile = paste(contig, "cov.jpg", sep = "_")
ggsave(outfile , plot = p0, device = 'jpg', path = NULL,
scale = 1, width = 500, height = 500, units = 'mm',
dpi = 150, limitsize = TRUE)
}



df_WMZ5 <- read.delim("~/Downloads/Mg/vs_WMZ5_grouped_depth.tsv", header=FALSE)

colnames(df_WMZ5) <- c("contig","position", "depth", "strain")
# install.packages("readr")
library(readr)
df_WMZ5$treatment <- paste(df_WMZ5$strain , df_WMZ5$contig)
tapply(df_WMZ5$depth, df_WMZ5$treatment, myFun)
df2 <- cbind(do.call(rbind, tapply(df_WMZ5$depth, df_WMZ5$treatment, myFun)))
write.csv(df2, 'WMZ5_contig_coverage.csv')
df_WMZ5$depth <- ifelse(df_WMZ5$depth > 200, 200, df_WMZ5$depth)



for (contig in levels(df_WMZ5$contig)){
# contig = paste("tig0000000", i, "_pilon_pilon_pilon", sep = "")
p0 <- ggplot(data=df_WMZ5[df_WMZ5$contig == contig, ], aes(x=`position`, y=`depth`, group=1)) +
geom_line() +
labs(x = "Position (bp)", y = "Coverage") +
scale_y_continuous(breaks=seq(0,200,50), limits=c(0,200)) +
facet_wrap(~strain, nrow = 19, ncol = 1, strip.position = "left")
outfile = paste(contig, "cov.jpg", sep = "_")
ggsave(outfile , plot = p0, device = 'jpg', path = NULL,
scale = 1, width = 500, height = 500, units = 'mm',
dpi = 150, limitsize = TRUE)
}



```
