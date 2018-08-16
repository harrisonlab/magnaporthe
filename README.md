# magnaporthe
Commands used in the analysis of Magnaporthe genomes


## Make project directory

```bash
mkdir -p /home/groups/harrisonlab/project_files/magnaporthe
cd /home/groups/harrisonlab/project_files/magnaporthe
mkdir -p assembly/taiwo/M.grisea
```

from local machine
```bash
scp * cluster:/home/groups/harrisonlab/project_files/magnaporthe/assembly/taiwo/M.grisea/.
```

```bash
for File in $(ls assembly/taiwo/M.grisea/*.fasta); do
  Strain=$(basename ${File%.fasta})
  mkdir assembly/taiwo/M.grisea/$Strain
  mv $File assembly/taiwo/M.grisea/$Strain/.
done
```

## Download Reference genomes

https://www.ncbi.nlm.nih.gov/Traces/wgs/?val=LNTK01  -MG04
ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/fungi/Magnaporthe_sp._MG03/latest_assembly_versions/    - MG03
ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/fungi/Magnaporthe_sp._MG12/latest_assembly_versions/     -MG12
https://www.ncbi.nlm.nih.gov/Traces/wgs/?val=AACU03  -Mo70-15
https://www.ncbi.nlm.nih.gov/Traces/wgs/?val=MQOP01  -Guy11
https://www.ncbi.nlm.nih.gov/Traces/wgs/?val=MQOO01  -FJ81278
https://www.ncbi.nlm.nih.gov/Traces/wgs/?val=UELY01     -EI
https://www.ncbi.nlm.nih.gov/Traces/wgs/?val=UEMA01    -WT
https://www.ncbi.nlm.nih.gov/Traces/wgs/?val=UCOF01    -SI
https://www.ncbi.nlm.nih.gov/Traces/wgs/?val=AABX03    - NC (Outgroup)

```bash
CurDir=$PWD

# MG04
Orgnaism=M.oryzae
Strain=MG04
OutDir=assembly/external/$Organism/$Strain
mkdir -p $OutDir
cd $OutDir
wget ftp://ftp.ncbi.nlm.nih.gov/sra/wgs_aux/LN/TK/LNTK01/LNTK01.1.fsa_nt.gz
cat *.fsa_nt.gz | gunzip -cf > genome.ctg.fa
cd $CurDir

# MG03
Organism=M.sp.
Strain=MG03
OutDir=assembly/external/$Organism/$Strain
mkdir -p $OutDir
cd $OutDir
wget ftp://ftp.ncbi.nlm.nih.gov/sra/wgs_aux/LN/TJ/LNTJ01/LNTJ01.1.fsa_nt.gz
cat *.fsa_nt.gz | gunzip -cf > genome.ctg.fa
cd $CurDir

# MG12
Orgnaism=M.oryzae
Strain=MG12
OutDir=assembly/external/$Organism/$Strain
mkdir -p $OutDir
cd $OutDir
wget ftp://ftp.ncbi.nlm.nih.gov/sra/wgs_aux/LN/TO/LNTO01/LNTO01.1.fsa_nt.gz
cat *.fsa_nt.gz | gunzip -cf > genome.ctg.fa
cd $CurDir

# Mo70-15
Orgnaism=M.oryzae
Strain=M070-15
OutDir=assembly/external/$Organism/$Strain
mkdir -p $OutDir
cd $OutDir
wget ftp://ftp.ncbi.nlm.nih.gov/sra/wgs_aux/AA/CU/AACU03/AACU03.1.fsa_nt.gz
cat *.fsa_nt.gz | gunzip -cf > genome.ctg.fa
cd $CurDir

# Guy11
Orgnaism=M.oryzae
Strain=Guy11
OutDir=assembly/external/$Organism/$Strain
mkdir -p $OutDir
cd $OutDir
wget ftp://ftp.ncbi.nlm.nih.gov/sra/wgs_aux/MQ/OP/MQOP01/MQOP01.1.fsa_nt.gz
cat *.fsa_nt.gz | gunzip -cf > genome.ctg.fa
cd $CurDir

# FJ81278
Orgnaism=M.oryzae
Strain=FJ81278
OutDir=assembly/external/$Organism/$Strain
mkdir -p $OutDir
cd $OutDir
wget ftp://ftp.ncbi.nlm.nih.gov/sra/wgs_aux/MQ/OO/MQOO01/MQOO01.1.fsa_nt.gz
cat *.fsa_nt.gz | gunzip -cf > genome.ctg.fa
cd $CurDir

# EI
Orgnaism=M.oryzae
Strain=EI
OutDir=assembly/external/$Organism/$Strain
mkdir -p $OutDir
cd $OutDir
wget ftp://ftp.ncbi.nlm.nih.gov/sra/wgs_aux/UE/LY/UELY01/UELY01.1.fsa_nt.gz
cat *.fsa_nt.gz | gunzip -cf > genome.ctg.fa
cd $CurDir

# WT
Orgnaism=M.oryzae
Strain=WT
OutDir=assembly/external/$Organism/$Strain
mkdir -p $OutDir
cd $OutDir
wget ftp://ftp.ncbi.nlm.nih.gov/sra/wgs_aux/UE/MA/UEMA01/UEMA01.1.fsa_nt.gz
cat *.fsa_nt.gz | gunzip -cf > genome.ctg.fa
cd $CurDir

# SI
Orgnaism=M.oryzae
Strain=SI
OutDir=assembly/external/$Organism/$Strain
mkdir -p $OutDir
cd $OutDir
wget ftp://ftp.ncbi.nlm.nih.gov/sra/wgs_aux/UC/OF/UCOF01/UCOF01.1.fsa_nt.gz
cat *.fsa_nt.gz | gunzip -cf > genome.ctg.fa
cd $CurDir

# NC (Outgroup)
Orgnaism=N.crassa
Strain=NC
OutDir=assembly/external/$Organism/$Strain
mkdir -p $OutDir
cd $OutDir
wget ftp://ftp.ncbi.nlm.nih.gov/sra/wgs_aux/AA/BX/AABX03/AABX03.1.fsa_nt.gz
cat *.fsa_nt.gz | gunzip -cf > genome.ctg.fa
cd $CurDir
```
