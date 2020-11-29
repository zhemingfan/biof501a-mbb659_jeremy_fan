rule all:
    input:
        "results/lineage_report.csv",
        "results/covid_histogram.png"

rule retrieve_reads:
    output:
        reads="data/SRR11801823_combined.fastq"
    shell:
        """
        fasterq-dump --split-files SRR11801823 -O data 
        cat data/SRR11801823_1.fastq data/SRR11801823_2.fastq > data/SRR11801823_combined.fastq && rm data/SRR11801823_1.fastq data/SRR11801823_2.fastq
        """

rule map_reads:
    input:
        ref="data/reference_sequence.fasta",
        reads=rules.retrieve_reads.output.reads
        #reads="data/SRR11801823_combined.fastq"
    output:
        bam="SRR11801823.alignments.sorted.bam"
    threads: 1
    shell:
        """
        minimap2 -t {threads} -ax sr {input.ref} {input.reads} | samtools sort -o {output.bam} 
        """

rule create_fasta:
    input:
        ref="data/reference_sequence.fasta",
        bam=rules.map_reads.output.bam
        #bam="SRR11801823.alignments.sorted.bam"
    output:
        fasta="SRR11801823.fasta"
    conda:
        "envs/create_fasta.yml"
    threads:
        1
    shell:
        """
        samtools mpileup -uf {input.ref} {input.bam}| bcftools call -c | vcfutils.pl vcf2fq > covid.fastq
        seqtk seq -aQ64 -q20 -n N covid.fastq > SRR11801823.fasta
        """ 

rule create_csv:
    input:
        visualize_fasta=rules.create_fasta.output.fasta
	    #visualize_fasta="SRR11801823.fasta"
    output:
	    csv="results/lineage_report.csv"
    shell:
        """
        pangolin {input.visualize_fasta} -o results
        """

rule produce_histogram:
    input:
        #input_csv=rules.create_csv.output.csv
        input_csv="results/lineage_report.csv"
    output:
        var_plot="results/covid_histogram.png",
    shell:
        """
        python3 plot_histogram.py {input.input_csv}
        """


### finall plot the pictures
### TODO: 

# 1) MAKE SURE PIPELINE CAN RE-RUN
# 3) add in more read files.. 
# 4) get rid of warning messages 

# rule dag:
#     output:
#         dag="snakemake_workflow_dag.png"
#     shell:
#         """
#         snakemake parse_variants --dag | dot -Tpng > {output.dag}
#         """