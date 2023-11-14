# Scripts to filter and merge VCFs as well as make some basic diagnostic plots.

1. Use mamba to build environment from environment.yml
2. Place scripts and bed file in the top directory.
    - For nextflow's rnavar output, use scripts in the rnavar folder to accommodate rnavar's folder structure this would be the directory called variant_calling. (Target files will be in a subfolder with the same name as the file prefix.)
4. Run merge.py to perform some basic filtering and merge resulting vcfs into all.vcf.gz. Will automatically run filt_plot.R to perform further filtering and output plots.
5. filt_plot.R can be run on its own if you have a merged vcf file already.

## merge.py command:

python3 ./merge.py DEPTH UNZIPPED FINAL_VAR_COUNT SAMPLE_NAMES

DEPTH = Minimum sequencing depth to filter by (default = 2000)
UNZIPPED = Whether or not input vcfs are zipped (0 = zipped (.vcf.gz)) or (1 = unzipped (.vcf)). (default = 0)
FINAL_VAR_COUNT = How many variants per sample to be included in the vcf. Randomly sampled. (default = 100)
SAMPLE_NAMES = Names of specific samples to include. Default is to include all samples.

i.e.
python3 ./merge.py 500 1 100 
    - Sequencing depth of 500 reads, target files have .vcf extension, will output 100 variants per sample, will target all files in current folder
python3 ./merge_rnavar.py 300 0 75 foo bar 
    - Sequencing depth of 300 reads, target files have .vcf.gz extension, will output 75 unique variants per sample, will target files foo.haplotypecaller.filtered.vcf.gz and bar.haplotypecaller.filtered.vcf.gz


## filt_plot.R command:
Rscript ./filt_plot.R FILENAME.vcf.gz FINAL_VAR_COUNT SAMPLE_NAMES


## folder structure based on rnavar output. Target files should have same sample name as their parent folder. 
```
./
+-- environment.yml
+-- filt_merge.sh
+-- filt_plot.R
+-- merge.py
+-- bedfilename.bed
+-- sample_name_one/
|   +-- sample_name_one.haplotypecaller.filtered.vcf*
|   +-- sample_name_one.haplotypecaller.filtered.vcf.gz**
|   +-- sample_name_one.haplotypecaller.filtered.vcf.gz.tbi
|   +-- sample_name_one.haplotypecaller.vcf.gz
|   +-- sample_name_one.haplotypecaller.vcf.gz.tbi
+-- sample_name_two/
.
.
.
+-- sample_name_N
```
*This file is made by the script. Will only be present if you have run it before
**This is the minimally required file
