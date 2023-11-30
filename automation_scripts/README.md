# Scripts to filter and merge VCFs as well as make some basic diagnostic plots
# Instructions
1. Use mamba to build environment from environment.yml
2. Place scripts and bed file in the same folder as vcfs
4. Run merge.py to perform some basic filtering and merge resulting vcfs into all.vcf.gz. Will automatically run filt_plot.R to perform further filtering and output plots.
   ```
   python3 ./merge.py DEPTH UNZIPPED FINAL_VAR_COUNT SAMPLE_NAMES
   ```
6. filt_plot.R can be run on its own if you have a merged vcf file already.
   ```
   Rscript ./filt_plot.R FILENAME.vcf.gz FINAL_VAR_COUNT SAMPLE_NAMES
   ```

# Parameters
* DEPTH = Minimum sequencing depth to filter by (default = 2000)
* UNZIPPED = Whether or not input vcfs are zipped (0 = zipped (.vcf.gz)) or (1 = unzipped (.vcf)). (default = 0)
* FINAL_VAR_COUNT = How many variants per sample to be included in the vcf. Randomly sampled. (default = 100)
* SAMPLE_NAMES = Names of specific samples to include. Default is to include all samples.


# Examples
```
python3 ./merge.py 500 1 100 
```
Sequencing depth of 500 reads, target files have .vcf extension, will output 100 variants per sample, will target all files in current folder
```
python3 ./merge_rnavar.py 300 0 75 foo bar 
```
Sequencing depth of 300 reads, target files have .vcf.gz extension, will output 75 unique variants per sample, will target files foo.vcf.gz and bar.vcf.gz
