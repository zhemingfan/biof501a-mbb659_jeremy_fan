rule all:
    input:
        "results/lineage_report.csv",
        "results/covid_histogram.png"

rule retrieve_reads:
    output:
        reads_1="data/SRR12960723_combined.fastq",
        reads_2="data/SRR12960724_combined.fastq",
        reads_3="data/SRR12960725_combined.fastq",
        reads_4="data/SRR12960726_combined.fastq",
        reads_5="data/SRR12960727_combined.fastq"
    shell:
        """
        fasterq-dump --split-files SRR12960723 -O data 
        cat data/SRR12960723_1.fastq data/SRR12960723_2.fastq > data/SRR12960723_combined.fastq && rm data/SRR12960723_1.fastq data/SRR12960723_2.fastq

        fasterq-dump --split-files SRR12960724 -O data 
        cat data/SRR12960724_1.fastq data/SRR12960724_2.fastq > data/SRR12960724_combined.fastq && rm data/SRR12960724_1.fastq data/SRR12960724_2.fastq

        fasterq-dump --split-files SRR12960725 -O data 
        cat data/SRR12960725_1.fastq data/SRR12960725_2.fastq > data/SRR12960725_combined.fastq && rm data/SRR12960725_1.fastq data/SRR12960725_2.fastq

        fasterq-dump --split-files SRR12960726 -O data 
        cat data/SRR12960726_1.fastq data/SRR12960726_2.fastq > data/SRR12960726_combined.fastq && rm data/SRR12960726_1.fastq data/SRR12960726_2.fastq

        fasterq-dump --split-files SRR12960727 -O data 
        cat data/SRR12960727_1.fastq data/SRR12960727_2.fastq > data/SRR12960727_combined.fastq && rm data/SRR12960727_1.fastq data/SRR12960727_2.fastq
        """


rule map_reads:
    input:
        ref="data/reference_sequence.fasta",
        reads_1=rules.retrieve_reads.output.reads_1,
        reads_2=rules.retrieve_reads.output.reads_2,
        reads_3=rules.retrieve_reads.output.reads_3,
        reads_4=rules.retrieve_reads.output.reads_4,
        reads_5=rules.retrieve_reads.output.reads_5

    output:
        bam_1="SRR12960723.alignments.sorted.bam",
        bam_2="SRR12960724.alignments.sorted.bam",
        bam_3="SRR12960725.alignments.sorted.bam",
        bam_4="SRR12960726.alignments.sorted.bam",
        bam_5="SRR12960727.alignments.sorted.bam"
    threads: 1
    shell:
        """
        minimap2 -t {threads} -ax sr {input.ref} {input.reads_1} | samtools sort -o {output.bam_1} 
        minimap2 -t {threads} -ax sr {input.ref} {input.reads_2} | samtools sort -o {output.bam_2} 
        minimap2 -t {threads} -ax sr {input.ref} {input.reads_3} | samtools sort -o {output.bam_3} 
        minimap2 -t {threads} -ax sr {input.ref} {input.reads_4} | samtools sort -o {output.bam_4} 
        minimap2 -t {threads} -ax sr {input.ref} {input.reads_5} | samtools sort -o {output.bam_5} 
        """

rule create_fasta:
    input:
        ref="data/reference_sequence.fasta",
        bam_1=rules.map_reads.output.bam_1,
        bam_2=rules.map_reads.output.bam_2,
        bam_3=rules.map_reads.output.bam_3,
        bam_4=rules.map_reads.output.bam_4,
        bam_5=rules.map_reads.output.bam_5

    output:
        fasta="combined.fasta"
    threads:
        1
    shell:
        """
        samtools mpileup -uf {input.ref} {input.bam_1}| bcftools call -c | vcfutils.pl vcf2fq > covid_1.fastq
        seqtk seq -aQ64 -q20 -n N covid_1.fastq > SRR11801823.fasta && rm covid_1.fastq

        samtools mpileup -uf {input.ref} {input.bam_2}| bcftools call -c | vcfutils.pl vcf2fq > covid_2.fastq
        seqtk seq -aQ64 -q20 -n N covid_2.fastq > SRR11801824.fasta && rm covid_2.fastq

        samtools mpileup -uf {input.ref} {input.bam_3}| bcftools call -c | vcfutils.pl vcf2fq > covid_3.fastq
        seqtk seq -aQ64 -q20 -n N covid_3.fastq > SRR11801825.fasta && rm covid_3.fastq

        samtools mpileup -uf {input.ref} {input.bam_4}| bcftools call -c | vcfutils.pl vcf2fq > covid_4.fastq
        seqtk seq -aQ64 -q20 -n N covid_4.fastq > SRR11801826.fasta && rm covid_4.fastq

        samtools mpileup -uf {input.ref} {input.bam_5}| bcftools call -c | vcfutils.pl vcf2fq > covid_5.fastq
        seqtk seq -aQ64 -q20 -n N covid_5.fastq > SRR11801827.fasta && rm covid_5.fastq

        cat SRR11801823.fasta SRR11801824.fasta SRR11801825.fasta SRR11801826.fasta SRR11801827.fasta > combined.fasta && rm SRR11801823.fasta SRR11801824.fasta SRR11801825.fasta SRR11801826.fasta SRR11801827.fasta
        """ 

rule create_csv:
    input:
        visualize_fasta=rules.create_fasta.output.fasta
    output:
	    csv="results/lineage_report.csv"
    shell:
        """
        python3 replace_header.py
        pangolin new_combined.fasta  -o results
        """

rule produce_histogram:
    input:
        input_csv=rules.create_csv.output.csv
    output:
        var_plot="results/covid_histogram.png",
    shell:
        """
        python3 plot_histogram.py {input.input_csv}
        """