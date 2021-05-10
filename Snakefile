
#################################### Defined by users #################################
configfile:"config/config_single.yaml"    # Sets path to the config file
#######################################################################################


# This workflow requires fastq.gz files in fastq directory 
# e.g. paired-end: DMSO_rep1_1.fastq.gz, DMSO_rep1_2.fastq.gz, Drug_rep1_1.fastq.gz, Drug_rep1_2.fastq.gz 
#      single-end: DMSO_rep1_1.fastq.gz, Drug_rep1_1.fastq.gz


shell.prefix('set -euo pipefail; ')
shell.executable('/bin/bash')

rule all:
    input:
        expand("hisat2_output/{sample}.bam", sample=config['FASTQ_PREFIX'])


rule get_reference:    
    """
    This rule downloads reference files
    """
    params:
        gen_link=config['REFERENCE_LINK']['GENOME'][0],   # Gencode reference genome file link 
        gen_name=config['REFERENCE_LINK']['GENOME'][1],   # Output reference genome location & name 
        anno_link=config['REFERENCE_LINK']['ANNOTATION'][0],  # Gencode GTF (annotation) file link
        anno_name=config['REFERENCE_LINK']['ANNOTATION'][1]   # Output GTF file location & name
    output:
        gen=expand("reference/{gen}", gen=config['REFERENCE_LINK']['GENOME'][2]),  # Decompressed reference genome file 
        anno=expand("reference/{anno}", anno=config['REFERENCE_LINK']['ANNOTATION'][2])  # Decompressed GTF file
    shell:
        "set +o pipefail; "
        "wget -c {params.gen_link} -O reference/{params.gen_name} && "
        "wget -c {params.anno_link} -O reference/{params.anno_name} && "
        "gzip -d reference/*.gz"

rule index_hisat2:
    """
    This rule constructs HISAT2 index files
    """
    input: 
        expand("reference/{gen}", gen=config['REFERENCE_LINK']['GENOME'][2])    # Decompressed reference genome file
    output:
        expand("reference/hisat2_index/hisat2_index.{number}.ht2", number=[x+1 for x in range(8)])   # HISAT2 indexing files
    threads: 8
    shell:
        "set +o pipefail; "
        "hisat2-build -f -o 4 "
        "-p {threads} "
        "--seed 67 "
        "{input} "
        "hisat2_index && "
        "mv *.ht2 reference/hisat2_index"

rule align_hisat2:    # Creates bam files in hisat2_output directory"
    """
    This rule aligns the reads using HISAT2    
    """
    input:
        fastq=expand("fastq/{sample}_{end}.fastq.gz", sample=config['FASTQ_PREFIX'], end=config['END']),  # Gzipped FASTQ files
        index=expand("reference/hisat2_index/hisat2_index.{number}.ht2", number=[x+1 for x in range(8)])  # HISAT2 indexing files
    output:
        expand("hisat2_output/{sample}.bam", sample=config['FASTQ_PREFIX']),   # Bam files
        temp(expand("hisat2_output/{sample}.sam", sample=config['FASTQ_PREFIX']))
    params:
        files=config["FASTQ_PREFIX"],  # e.g. Ctrl, Treatment
        ext=config['FASTQ_EXT'],       # extension of the FASTQ files (e.g. .fastq.gz)
        indexing=config['INDEX_HISAT'] # HISAT2 indexing location and file name prefix
    threads: 8
    run:
        for i in range(len(params.files)):
            p=params.files[i]
            r1= "fastq/" + params.files[i] + "_1" + params.ext 
            r2=""
            read="-U " + r1
            if len(input.fastq) == 2 * len(params.files): 
                r2= "fastq/" + params.files[i] + "_2" + params.ext  
                read="-1 " + r1 + " -2 " + r2
            shell("hisat2 -q -p {threads} "
                  "--seed 23 "
                  "--summary-file hisat2_output/summary_{p}.txt "
                  "-x {params.indexing} "
                  "{read} "
                  "-S hisat2_output/{p}.sam && "
                  "samtools view -bS "
                  "-@ {threads} "
                  "hisat2_output/{p}.sam > hisat2_output/{p}.bam")

