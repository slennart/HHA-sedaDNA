# HHA-sedaDNA

Analysis of shotgun sequencing &amp; hybridization capture data of Holocene High Arctic sedaDNA

### Overview 

Here, I provide the code that was used for the processing, filtering, and analysis of the shotgun sequencing and hybridization capture data of sedimentary ancient DNA from four marine sediment core from around Northern Greenland.

The processing of the sequencing data and taxonomic classification were carried out on a HPC. The code to process the data is part of this README. Further data processing, filtering and the analysis of the data was carried out in R. The R-scripts are part of this repository and the README guides through the order of running them.

Overview of used programs, versions used for the analysis, and the source publication


| Program  | Version | Publication |
| ------------- | ------------- | ------------- |
| leeHom  | 1.2.15  | Renaud et al. 2014 |
| sga  | 0.10.15  | Simpson &amp; Durbin 2012 |
| bowtie2  | 2.4.2  | Langmead &amp; Salzberg 2012 |
| samtools  | 1.20  | Li et al. 2009 |
| metaDMG  | 0.38.0-2 | Michelsen et al. 2022 |

### processing of raw fastq files 
```
leeHom -fq1 $i -fq2 $b -fqo $bname2 --ancientdna -t 4")
```

### filtering of merged reads using sga (duplicate removal, min. length 30 bp, homopolymer removal)
```
sga preprocess --dust-threshold=1 -m 30 $input/$bname  --no-primer-check -o $output/$bname2.dust.fq
sga index --algorithm=ropebwt --threads=4 $bname2.dust.fq
sga filter --threads=4  --no-kmer-check --homopolymer-check $bname2.dust.fq -o $bname2.dust.sga.fq
```

### bowtie2 alignment of filtered reads 
```
bowtie2 -k 1000 --threads 2 -x $DB -U <(zcat $basefolder/$bname) --no-unal | samtools sort -n -o  $output/$bname2.nsorted.bam
```

### additional duplicate removal step using samtools 
```
#sort by coordinates (required for samtools rmdup)
samtools sort $output/$bname2.nsorted.bam -o $output/$bname2.sorted.bam
#remove duplicates identical by start/stop position
samtools rmdup -s $output/$bname2.sorted.bam $output/$bname2.rmdup.sorted.bam
#sort by name (required for ngsLCA)
samtools sort -n $output/$bname2.rmdup.sorted.bam -o $output/$bname2.rmdup.nsorted.bam
```

### metaDMG (includes lca filtering step)
```
metaDMG config 03Sep24/*.filt.sorted.bam -s 0.98 -S 1 --max-position 25 -c config_PoolCap01 --names ncbi_docs/names.dmp --nodes ncbi_docs/nodes.dmp --acc2tax ncbi_docs/acc2taxid_mtgenomes_v4_171023.gz --custom-database --metaDMG-cpp metaDMG-cpp --output-dir metadmg/PoolCap01_11/03Sep24
metaDMG compute config_PoolCap01.yaml
```

### metaDMG converting damage files
```
metaDMG mismatch-to-mapDamage $i -o $input/mapDamage_files/$bname2.tsv;
```

### R scripts for processing the species-sample matrix, plotting &amp; correlating DNA detections with palaeoenvironmental proxies

1. `palaeoenvironment.R` for processing and plotting palaeoenvironmental proxies
2. `PoolScreen01_Analysis.R` for processing and plotting shotgun sequencing data
3. `PoolCap01_13_Analysis.R` for processing and plotting hybridization capture sequencing data
4. `HHA_sedaDNA_statistics.R` for running correlation analysis and redundancy analysis on the DNA detections and palaeoenvironmental proxies










