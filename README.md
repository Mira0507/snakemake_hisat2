
### HISAT2 alignment using Snakemake

#### 1. Conda environment

- References: [Conda doc](https://docs.conda.io/projects/conda/en/latest/index.html), [sra-tools](https://github.com/ncbi/sra-tools), [Snakemake Installation Guide](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
- config: [config/conda_env.yaml](https://github.com/Mira0507/snakemake_hisat2/blob/master/config/conda_env.yaml)


#### 2. Snakemake 

- Reference: [Snakemake doc](https://snakemake.readthedocs.io/en/stable), [HISAT2 manual](http://daehwankimlab.github.io/hisat2/manual), [samtools manual](http://www.htslib.org/doc/samtools.html)

- [Snakefile](https://github.com/Mira0507/snakemake_hisat2/blob/master/Snakefile)

```



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
```

- [config/config_single.yaml (single-end testing)](https://github.com/Mira0507/snakemake_hisat2/blob/master/config/config_single.yaml)


```yaml

###################### Sample info ######################


SAMPLE:
  DMSO_rep1: SRR13190144
  DMSO_rep2: SRR13190145
  DMSO_rep3: SRR13190146
  SR0813_rep1: SRR13190150
  SR0813_rep2: SRR13190151
  SR0813_rep3: SRR13190152


END: 
  - 1

###################### Reference info ######################

REFERENCE:
  LINK: 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/'
  FILE: 
    GENOME: 'GRCh38.primary_assembly.genome.fa'
    ANNOTATION: 'gencode.v37.primary_assembly.annotation.gtf'
    


  
# e.g. 
# Reference genome link: ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/GRCh38.primary_assembly.genome.fa.gz
# Reference annotation (GTF) link: ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.primary_assembly.annotation.gtf.gz 


###################### Extra-setting info ######################

INDEX_HISAT: "reference/hisat2_index/hisat2_index" # Assigns hisat2 index files (e.g.reference/hisat2_index/hisat2_index.1.ht2, reference/hisat2_index/hisat2_index.2.ht2, ...)


FASTQ_EXT: '.fastq.gz'
```


- [config/config_paired1.yaml (paired-end testing)](https://github.com/Mira0507/snakemake_hisat2/blob/master/config/config_paired1.yaml)


```yaml

###################### Sample info ######################


SAMPLE:
  Treated_rep1: SRR6461133
  Treated_rep2: SRR6461134
  Treated_rep3: SRR6461135
  Control_rep1: SRR6461139 
  Control_rep2: SRR6461140
  Control_rep3: SRR6461141



END: 
  - 1
  - 2




###################### Reference info ######################

REFERENCE:
  LINK: 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/'
  FILE: 
    GENOME: 'GRCh38.primary_assembly.genome.fa'
    ANNOTATION: 'gencode.v37.primary_assembly.annotation.gtf'
  
# e.g. 
# Reference genome link: ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/GRCh38.primary_assembly.genome.fa.gz
# Reference annotation (GTF) link: ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.primary_assembly.annotation.gtf.gz 

###################### Extra-setting info ######################

INDEX_HISAT: "reference/hisat2_index/hisat2_index" # Assigns hisat2 index files (e.g.reference/hisat2_index/hisat2_index.1.ht2, reference/hisat2_index/hisat2_index.2.ht2, ...)

FASTQ_EXT: '.fastq.gz'
```

#### 3. Running the Snakemake workflow

- Reference: [Snakemake Command Line Arguments](https://snakemake.readthedocs.io/en/stable/executing/cli.html)

- **Dry run**


```bash
#!/bin/bash

snakemake -n

```


- **DAG visualization** 


```bash
#!/bin/bash


# The dot commend requires graphviz (downloadable via conda)
snakemake --dag | dot -Tpdf > dag.pdf


```


- **Run**


```bash
#!/bin/bash

# Either -j or --cores assignes the number of cores
snakemake -j 16

```
