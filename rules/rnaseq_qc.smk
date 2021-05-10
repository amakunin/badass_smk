
rule cutadapt:
    input:
        R1="{species}/transcriptomic_data/{sample}/rna-seq/{sample}.R1.fastq.gz",
        R2="{species}/transcriptomic_data/{sample}/rna-seq/{sample}.R2.fastq.gz"
    output:
        R1="{species}/transcriptomic_data/{sample}/rna-seq/trimmed/{sample}.R1.fastq.gz",
        R2="{species}/transcriptomic_data/{sample}/rna-seq/trimmed/{sample}.R2.fastq.gz"
    log: "{species}/transcriptomic_data/{sample}/rna-seq/trimmed/cutadapt.log"
    params:
        min_len=20
    conda: "rnaseq_qc.yml"
    shell:
        "cutadapt -m {params.min_len} -a AGATCGGAAGAGC -A AGATCGGAAGAGC -o {output.R1} -p {output.R2} "
        "{input.R1} {input.R2} > {log}"

rule star_index:
    input:
        "{species}/working/{dna_sample}.{asm}.{date}/{purge_dir}/purged.fa"
    output:
        touch("{species}/working/{dna_sample}.{asm}.{date}/{purge_dir}/star_index/index.done")
    params:
        genome_dir="{species}/working/{dna_sample}.{asm}.{date}/{purge_dir}/star_index/"
    conda: "rnaseq_qc.yml"
    threads: 8
    resources: mem_mb=10000
    shell:
        "STAR --runThreadN {threads} --runMode genomeGenerate --genomeSAindexNbases 13 "
        "--genomeDir {params.genome_dir} --genomeFastaFiles {input}"

rule star_align:
    input:
        R1="{species}/transcriptomic_data/{rna_sample}/rna-seq/trimmed/{rna_sample}.R1.fastq.gz",
        R2="{species}/transcriptomic_data/{rna_sample}/rna-seq/trimmed/{rna_sample}.R2.fastq.gz",
        ref="{species}/working/{dna_sample}.{asm}.{date}/{purge_dir}/star_index/index.done"
    output:
        touch("{species}/working/{rna_sample}.star.{dna_sample}.{asm}.{date}.{purge_dir}/{rna_sample}.star.done"),
        "{species}/working/{rna_sample}.star.{dna_sample}.{asm}.{date}.{purge_dir}/{rna_sample}.Aligned.sortedByCoord.out.bam"
    params:
        genome_dir="{species}/working/{dna_sample}.{asm}.{date}/{purge_dir}/star_index/",
        out_prefix="{species}/working/{rna_sample}.star.{dna_sample}.{asm}.{date}.{purge_dir}/{rna_sample}."
    conda: "rnaseq_qc.yml"
    threads: 16
    resources: mem_mb=24000
    shell:
        "STAR --runThreadN {threads} --genomeDir {params.genome_dir} "
        "--readFilesIn {input.R1} {input.R2} --readFilesCommand zcat "
        "--outFileNamePrefix {params.out_prefix} --outSAMtype BAM SortedByCoordinate "
        "--outSAMunmapped Within --outSAMattributes NH HI NM MD AS"

rule mark_duplicates:
    input:
        "{species}/working/{rna_sample}.star.{dna_sample}.{asm}.{date}.{purge_dir}/{rna_sample}.Aligned.sortedByCoord.out.bam"
    output:
        bam="{species}/working/{rna_sample}.star.{dna_sample}.{asm}.{date}.{purge_dir}/{rna_sample}.markdup.bam",
        metrics="{species}/working/{rna_sample}.star.{dna_sample}.{asm}.{date}.{purge_dir}/{rna_sample}.markdup.txt"
    conda: "rnaseq_qc.yml"
    resources: mem_mb=2048
    shell:
        "picard MarkDuplicates I={input} O={output.bam} M={output.metrics}"

rule samtools_stats:
    input:
        bam="{species}/working/{rna_sample}.star.{dna_sample}.{asm}.{date}.{purge_dir}/{rna_sample}.markdup.bam",
        ref="{species}/working/{dna_sample}.{asm}.{date}/{purge_dir}/purged.fa.gz"
    output:
        "{species}/working/{rna_sample}.star.{dna_sample}.{asm}.{date}.{purge_dir}/{rna_sample}.markdup.stats"
    conda: "rnaseq_qc.yml"
    resources: mem_mb=2048
    shell:
        "samtools stats -r {input.ref} {input.bam} > {output}"
