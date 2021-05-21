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

## Reference-based

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
        ref="{species}/working/{dna_sample}.{asm}.{date}/{purge_dir}/purged.fa"
    output:
        "{species}/working/{rna_sample}.star.{dna_sample}.{asm}.{date}.{purge_dir}/{rna_sample}.markdup.stats"
    conda: "rnaseq_qc.yml"
    resources: mem_mb=2048
    shell:
        "samtools stats -r {input.ref} {input.bam} > {output}"

## assembly-based

rule reads_for_trinity:
    input:
        "{species}/transcriptomic_data/{sample}/rna-seq/trimmed/{sample}.R{read}.fastq.gz"
    output:
        temp("{species}/working/{sample}.trinity/{sample}.R{read}.fasta")
    params: read="{read}"
    conda: "fastool.yml"
    shell:
        "zcat {input} | fastool --to-fasta --append /{params.read} > {output}"

rule trinity:
    input:
        left="{species}/working/{sample}.trinity/{sample}.R1.fasta",
        right="{species}/working/{sample}.trinity/{sample}.R2.fasta"
    output:
        "{species}/working/{sample}.trinity/Trinity.fasta"
    log:
        "{species}/working/{sample}.trinity/Trinity.log"
    params:
        queue="long",
        max_mem="20G", # link to resource
        out_dir="{species}/working/{sample}.trinity/",
        extra=""
    threads: 8
    resources:
        mem_mb=20 * 1024
    # wrapper:
    #     "0.74.0/bio/trinity"
    singularity: "scripts/trinityrnaseq-v2.12.0.sif"
    shell:
        "Trinity --left {input.left} --right {input.right} --CPU {threads} "
        " --max_memory {params.max_mem} --seqType fa --SS_lib_type RF "
        " --output {params.out_dir} {params.extra} 2> {log}"

rule trinity_archive:
    input:
        "{species}/working/{sample}.trinity/Trinity.fasta"
    output:
        "{species}/working/{sample}.trinity/read_partitions.tar.gz"
    params:
        arc_dir="{species}/working/{sample}.trinity/read_partitions"
    shell:
        "tar -czf {output} {params.arc_dir} && rm -rf {params.arc_dir}"
