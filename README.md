# 16S_Pipeline
# 16S sequencing data pipeline

```bash
for R1 in Raw_data/Bacteria*_1.fq.gz; do
    R2=$(echo $R1|sed 's/_1/_2/')
    name=$(echo $R1|sed 's/_raw_1.fq.gz//'|sed 's/Raw_data\///')
    echo "fastp -i $R1 \
                 -I $R2 \
                 -o Raw_data/Trimmed_Reads/$name.trim.R1.fq.gz \
                 -O Raw_data/Trimmed_Reads/$name.trim.R2.fq.gz \
                 -w 6 \
                 --detect_adapter_for_pe \
                 -j Raw_data/logfile/$name.json \
                 -h Raw_data/logfile/$name.html \
                 &> Raw_data/logfile/$name.log.txt"
done > Raw_data/fastp.commands

```

```bash
source Raw_data/fastp.commands

## refer to https://wikis.utexas.edu/display/bioiteam/fastp+-+GVA2023
  
  ls Downloads/Raw_data/Trimmed_Reads/*.fq.gz > Downloads/Raw_data/file_list.txt
```

```r
**## Build manifest file

mydata<-read.csv("/Users/yunliangli/Downloads/Raw_data/file_list.txt",header=F)
manifest <- data.frame(mydata,stringsAsFactors = F)colnames(manifest)[1]<-c("absolute-filepath")
manifest$`absolute-filepath`<- paste("/Users/yunliangli/",manifest$`absolute-filepath`,sep='')  # absolute path to the target files is necessary
sample_id <- sub(".*Trimmed_Reads/(.*).trim.*", "\\1", manifest$`absolute-filepath`)
sample_id <- gsub('Bacteria_','Bacteria.',sample_id) #deblur cannot handle underscore in sample id
direction <- rep(c("forward","reverse"), dim(mydata)[1]/2)
manifest1 <- data.frame('sample-id'=sample_id,'absolute-filepath'=manifest$`absolute-filepath`,'direction'=direction,check.names = F) # check.names = F keep '-' from 
changing into '.'
write_csv(manifest1, file="/Users/yunliangli/Downloads/Raw_data/manifest.csv",quote='none') # NOT write.csv, write_csv function is from readr package, which can remove quote 
"" of names which is required for the data import into Qiime using 'PairedEndFastqManifestPhred33'**

                       
```

```bash
**## import data

conda activate qiime2-amplicon-2024.2
mkdir Downloads/Raw_data/reads_qza

#Use the import command with the corresponding parameters in qiime2
      
qiime tools import \
      --type "SampleData[PairedEndSequencesWithQuality]" \
      --input-path Downloads/Raw_data/manifest.csv \
      --output-path Downloads/Raw_data/reads_qza/reads_trimmed.qza \
      --input-format PairedEndFastqManifestPhred33
      
#############################
# What is the difference between PairedEndFastqManifestPhred33 and CasavaOneEightSingleLanePerSampleDirFmt
##########################
 
## Reference:  Qiime importing data

##PairedEndFastqManifestPhred33.      In this variant of the fastq manifest format, there must be forward and reverse read fastq.gz / fastq files for each sample id. As a 
result, each sample id is represented twice in this file: once for its forward reads, and once for its reverse reads. This format assumes that the PHRED offset used for the 
positional quality scores in all of the fastq.gz / fastq files is 33.

##Casava 1.8 paired-end demultiplexed fastq 
##Format description

## In this format, there are two fastq.gz file for each sample in the study, and the file name includes the sample identifier. The forward and reverse read file names for a 
single sample might look like L2S357_15_L001_R1_001.fastq.gz and L2S357_15_L001_R2_001.fastq.gz, respectively. The underscore-separated fields in this file name are the 
sample identifier, the barcode sequence or a barcode identifier, the lane number, the read number, and the set number.**
                                         
qiime tools validate **Downloads/Raw_data/**reads_qza/reads_trimmed.qza  

qiime vsearch merge-pairs \
  --i-demultiplexed-seqs Downloads/Raw_data/reads_qza/reads_trimmed.qza \
  --p-allowmergestagger \
  --o-merged-sequences Downloads/Raw_data/reads_qza/reads_trimmed_merged.qza \
  --o-unmerged-sequences Downloads/Raw_data/reads_qza/reads_trimmed_unmerged.qza                                     
  
  #--p-maxdiffs 4 \ 
   
## Make a qiime2 visualization object to look at using qiime2 view. You'll use this to set your --p-trim-length in the deblur step that follows.
qiime demux summarize \
  --i-data Downloads/Raw_data/reads_qza/reads_trimmed_merged.qza \
  --o-visualization Downloads/Raw_data/reads_qza/reads_trimmed_merged_summary.qzv

qiime demux summarize \
  --i-data Downloads/Raw_data/reads_qza/reads_trimmed_unmerged.qza \
  --o-visualization Downloads/Raw_data/reads_qza/reads_trimmed_unmerged_summary.qzv 
  
## unmerged reads: 49444.    merged reads:  4210481
#Demultiplexed sequence length summary

#  Total Sequences Sampled	10000.0
#  2%	444 nts
#  9%	444 nts
#  25%	446 nts
#  50% (Median)	467 nts
#  75%	469 nts
#  91%	470 nts
#  98%	471 nts

```

```bash

# 
qiime deblur denoise-16S \
  --i-demultiplexed-seqs Downloads/Raw_data/reads_qza/reads_trimmed_merged.qza \
  --p-trim-length 444 \
  --o-representative-sequences Downloads/Raw_data/reads_qza/rep-seqs-deblur.qza \
  --o-table Downloads/Raw_data/reads_qza/table-deblur.qza \
  --p-sample-stats \
  --o-stats Downloads/Raw_data/reads_qza/deblur-stats.qza

#Plugin error from deblur:

 # Deblur cannot operate on sample IDs that contain underscores. The following ID(s) contain one or more underscores: Bacteria_NRF1, Bacteria_NRF2, Bacteria_NRF3, 
Bacteria_NRK1, Bacteria_NRK2, Bacteria_NRK3, Bacteria_NRN1, Bacteria_NRN2, Bacteria_NRN3.  
 
 Replace '_' with '.' for the sample-ID column in manifest file
 

```

 

```bash
qiime deblur visualize-stats \
  --i-deblur-stats Downloads/Raw_data/reads_qza/deblur-stats.qza \
  --o-visualization Downloads/Raw_data/reads_qza/deblur-stats.qzv
  
  

mv Downloads/Raw_data/reads_qza/rep-seqs-deblur.qza Downloads/Raw_data/reads_qza/rep-seqs.qza
mv Downloads/Raw_data/reads_qza/table-deblur.qza Downloads/Raw_data/reads_qza/table.qza
```

```bash
qiime feature-table tabulate-seqs \
  --i-data Downloads/Raw_data/reads_qza/rep-seqs.qza \
  --o-visualization Downloads/Raw_data/reads_qza/rep-seqs.qzv
```

```bash
qiime feature-classifier extract-reads \
  --i-sequences Downloads/silva-138-99-seqs.qza \
  --p-f-primer ACTCCTACGGGAGGCAGCA \
  --p-r-primer GGACTACHVGGGTWTCTAAT \
  --p-n-jobs 6 \
  --o-reads Downloads/ref_seqs.qza
  
  # Process silva data by yourself, refer to RESCRIPt: Reproducible sequence taxonomy reference database management
  # This command can handle degenerate primers.  Refer to https://forum.qiime2.org/t/allowed-primer-read-mismatch-in-extract-reads/12650
```

```bash
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads Downloads/ref_seqs.qza \
  --i-reference-taxonomy Downloads/silva-138-99-tax.qza \
  --o-classifier Downloads/SILVA_V3_V4_338F_806R_qiime2_2024_02_classifier.qza
  
 # RESCRIPt: Reproducible sequence taxonomy reference database management
```

```bash
qiime feature-classifier classify-sklearn \
  --i-classifier Downloads/SILVA_V3_V4_338F_806R_qiime2_2024_02_classifier.qza \
  --p-n-jobs 6 \
  --i-reads Downloads/Raw_data/reads_qza/rep-seqs.qza \
  --o-classification Downloads/Raw_data/reads_qza/taxonomy.qza

```

```
qiime metadata tabulate \
  --m-input-file Downloads/Raw_data/reads_qza/taxonomy.qza \
  --o-visualization Downloads/Raw_data/reads_qza/taxonomy.qzv

qiime tools export --input-path Downloads/Raw_data/reads_qza/table.qza --output-path Downloads/Raw_data/reads_qza
qiime tools export --input-path Downloads/Raw_data/reads_qza/taxonomy.qza --output-path Downloads/Raw_data/reads_qza

biom convert -i Downloads/Raw_data/reads_qza/feature-table.biom -o Downloads/Raw_data/reads_qza/feature-table.tsv --to-tsv
```
