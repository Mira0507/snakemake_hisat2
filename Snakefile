
#################################### Defined by users #################################
configfile:"config/config_single.yaml"    # Sets path to the config file

#######################################################################################

THREADS=config["THREADS"]

shell.prefix('set -euo pipefail; ')
shell.executable('/bin/bash')

rule all:
    input:
        expand("hisat2_output/{sample}.bam", sample=config['FASTQ_PREFIX'])



rule align_hisat2: 
    """
    This rule aligns the reads using HISAT2
    """
    output:
        expand("hisat2_output/{sample}.bam", sample=config['FASTQ_PREFIX']),
        temp(expand("hisat2_output/{sample}.sam", sample=config['FASTQ_PREFIX']))
    params:
        indir=config['FASTQ_DIR'],
        files=config["FASTQ_PREFIX"], 
        read_ends=config['END'],
        ext=config['FASTQ_EXT'],
        index=config['INDEX_HISAT']
    run:
        for i in range(len(params.files)):
            p=params.files[i]
            r1= params.indir + params.files[i] + "_1" + params.ext 
            r2=""
            read="-U " + r1
            if len(params.read_ends) == 2: 
                r2= params.indir + params.files[i] + "_2" + params.ext  
                read="-1 " + r1 + " -2 " + r2
            shell("hisat2 -q -p {THREADS} "
                  "--seed 23 "
                  "-x {params.index} "
                  "{read} "
                  "-S hisat2_output/{p}.sam && "
                  "samtools view -bS "
                  "-@ {THREADS} "
                  "hisat2_output/{p}.sam > hisat2_output/{p}.bam")
