import sys
import subprocess
import glob
import os


if __name__ == "__main__":
    #Set default values
    depth = 2000
    zipped = 0
    allfolders = 1
    if len(sys.argv)>1:
        #check if input passed by user is safe and sensical
        if 0 <= int(sys.argv[1]):
            depth = int(sys.argv[1])
        else:
            sys.exit("Invalid Args: Depth must be an integer larger than 0")
        #check arg determining if vcfs are already unzipped (i.e. if this script has already been run before)
        if len(sys.argv)>2:
            if (int(sys.argv[2]) == 1):
                zipped = 1
            if (int(sys.argv[2]) not in [0, 1]):
                sys.exit("Invalid Args: please use 0 (zipped (.vcf.gz)) or 1 (unzipped (.vcf)) to indicate source file type")
        #check if specific samples are specified
        if len(sys.argv)>3:
            allfolders = 0
            FILES  = sys.argv[3:]
            if (int(sys.argv[2]) == 1):
                for i in range(0,len(FILES)):
                    FILES[i] = (FILES[i] + '.haplotypecaller.filtered.vcf')
            else:
                for i in range(0,len(FILES)):
                    FILES[i] = (FILES[i] + '.haplotypecaller.filtered.vcf.gz')
           
    #get list of files with this suffix in current directory's subfolders
    if allfolders == 1:
        if (zipped == 1):
            FILES = [os.path.basename(x) for x in glob.glob('./**/*.haplotypecaller.filtered.vcf')]
        else:
            FILES = [os.path.basename(x) for x in glob.glob('./**/*.haplotypecaller.filtered.vcf.gz')]
    
    MERGE_ARGS = ['bcftools', 'merge']

    for F in FILES:
        if (zipped == 0):
            folder = F.removesuffix('.haplotypecaller.filtered.vcf.gz')
            path = folder + '/' + F
            subprocess.run(['gunzip', '-k', path])
            path = folder + '/' + (F.removesuffix('.gz'))
        else:
            folder = F.removesuffix('.haplotypecaller.filtered.vcf')
            path = folder + '/' + F

        MERGE_ARGS.append(folder + '/' + folder + ('.DP%s.fordemux.filt.recode.vcf.gz' % depth))
        runinfo  = ' INFO/DP>' + str(depth)
        subprocess.run(['./filt_merge.sh', path, runinfo, str(depth)])

    MERGE_ARGS.append('-o')
    MERGE_ARGS.append('all.vcf.gz')

    subprocess.run(MERGE_ARGS)
    subprocess.run(['Rscript', './filt_plot.R', 'all.vcf.gz'])
