import sys
import subprocess
import glob
import os


if __name__ == "__main__":
    #Set default values
    depth = 2000
    if len(sys.argv)>1:
        #check if input passed by user is safe and sensical
        if sys.argv[1].is_integer() & 0 <= sys.argv[1]:
            depth = sys.argv[1]
        else:
            raise ScriptExit("Depth must be an integer larger than 0")
    
    #get list of files with this suffix in current directory's subfolders
    FILES = [os.path.basename(x) for x in glob.glob('./**/*.haplotypecaller.filtered.vcf')]
    MERGE_ARGS = ['bcftools', 'merge']

    for F in FILES:
        folder = F.removesuffix('.haplotypecaller.filtered.vcf')
        path = folder + '/' + F
        MERGE_ARGS.append(folder + '/' + folder + '.DP2000.fordemux.filt.recode.vcf.gz')
        runinfo  = ' INFO/DP>' + str(depth)
        subprocess.run(['./filt_merge.sh', path, runinfo])

    MERGE_ARGS.append('-o')
    MERGE_ARGS.append('all.vcf.gz')

    subprocess.run(MERGE_ARGS)
    subprocess.run(['Rscript', './filt_plot.R'])
