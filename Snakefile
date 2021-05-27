
#################################### Defined by users #################################
configfile:"config/config_paired1.yaml"    # Sets path to the config file
#######################################################################################


# This workflow requires fastq.gz files in fastq directory 
# e.g. paired-end: DMSO_rep1_1.fastq.gz, DMSO_rep1_2.fastq.gz, Drug_rep1_1.fastq.gz, Drug_rep1_2.fastq.gz 
#      single-end: DMSO_rep1_1.fastq.gz, Drug_rep1_1.fastq.gz


rule all: 
    input: 
        expand("reference/{ref}", ref=config['REFERENCE']['FILE']['GENOME']),  # Reference genome and annotation (GTF) files
        expand("hisat2_output/{sample}.bam", sample=list(config['SAMPLE'].keys())),   # HISAT2 output BAM files


rule get_reference:    
    """
    This rule downloads and decompresses reference files
    """
    params:
        link=config['REFERENCE']['LINK'] + config['REFERENCE']['FILE']['GENOME']
    output:
        "reference/{ref}"   # Decompressed reference genome
    shell:
        "wget -c {params.link}.gz -O {output}.gz && "
        "gzip -d {output}.gz"


rule index_hisat2:
    """
    This rule constructs HISAT2 index files
    """
    input: 
        expand("reference/{ref}", ref=config['REFERENCE']['FILE']['GENOME'])   # Decompressed reference genome
    output:
        expand("{path}.{number}.ht2", path=config['INDEX_HISAT'], number=[x for x in range(1, 9)])   # HISAT2 indexing files
    threads: 16
    shell:
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
        fastq=expand("fastq/{{sample}}_{end}.fastq.gz", end=config["END"]),   # Gzipped FASTQ files
        index=expand("{path}.{number}.ht2", path=config['INDEX_HISAT'], number=[x for x in range(1, 9)])  # HISAT2 indexing files
    output:
        temp("hisat2_output/{sample}.sam"),
        "hisat2_output/{sample}.bam"    # Bam files
    params:
        ext=config['FASTQ_EXT'],       # extension of the FASTQ files (e.g. .fastq.gz)
        indexing=config['INDEX_HISAT'] # HISAT2 indexing location and file name prefix
    threads: 16
    run:
        r1="fastq/" + wildcards.sample + "_1" + params.ext
        r2=""
        read="-U " + r1
        if len(input.fastq) == 2:    # if paired-end
            r2= "fastq/" + wildcards.sample + "_2" + params.ext  
            read="-1 " + r1 + " -2 " + r2
        shell("hisat2 -q -p {threads} "
              "--seed 23 "
              "--summary-file hisat2_output/summary_{wildcards.sample}.txt "
              "-x {params.indexing} "
              "{read} "
              "-S {output[0]} && "
              "samtools view -bS "
              "-@ {threads} "
              "{output[0]} > {output[1]}")

