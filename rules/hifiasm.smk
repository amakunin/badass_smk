localrules: hifiasm_with_hic_init

# workaround - singularity does not mkdir
rule hifiasm_with_hic_init:
    input: "{species}/working/{sample}.hicanu.20210327/{sample}.trimmedReads.fasta.gz"
    output: "{species}/working/{sample}.hifiasm.{date}/init.done"
    shell: "touch {output}"

rule hifiasm_with_hic:
    input:
        pb="{species}/working/{sample}.hicanu.20210327/{sample}.trimmedReads.fasta.gz",
        f="{species}/genomic_data/{sample_hic}/hic-arima2/{sample_hic}.R1.fastq.gz",
        r="{species}/genomic_data/{sample_hic}/hic-arima2/{sample_hic}.R2.fastq.gz",
        init="{species}/working/{sample}.hifiasm.{date}/init.done"
    output:
        touch("{species}/working/{sample}.hifiasm.{date}/{sample}.hifiasm.{sample_hic}.hic.done"),
    params:
        outfolder="{species}/working/{sample}.hifiasm.{date}",
        pb="../{sample}.hicanu.20210327/{sample}.trimmedReads.fasta.gz",
        f="../../genomic_data/{sample_hic}/hic-arima2/{sample_hic}.R1.fastq.gz",
        r="../../genomic_data/{sample_hic}/hic-arima2/{sample_hic}.R2.fastq.gz",
        sample="{sample}",
        log="{sample}.hifiasm.log"
    threads: 24
    resources: mem_mb=100000
    singularity: "scripts/hifiasm-0.15.sif"
    shell:
        "cd {params.outfolder} && "
        "hifiasm --primary -t{threads} -o {params.sample} --h1 {params.f} --h2 {params.r} {params.pb} 2> {params.log}"

rule hifiasm_with_hic_gfa_to_fa:
    # do not use haplotype-resolved assembly - incompatible with polishing
    input: 
        primary='{species}/working/{sample}.hifiasm.{date}/{sample}.hic.p_ctg.gfa',
        alternate='{species}/working/{sample}.hifiasm.{date}/{sample}.hic.a_ctg.gfa',
        # hap1='{species}/working/{sample}.hifiasm.{date}/{sample}.hic.hap1.p_ctg.gfa',
        # hap2='{species}/working/{sample}.hifiasm.{date}/{sample}.hic.hap2.p_ctg.gfa'
    output:
        primary='{species}/working/{sample}.hifiasm.{date}/{sample}.p_ctg.fa.gz',
        alternate='{species}/working/{sample}.hifiasm.{date}/{sample}.a_ctg.fa.gz',
        # hap1='{species}/working/{sample}.hifiasm.{date}/{sample}.p_ctg.fa.gz',
        # hap2='{species}/working/{sample}.hifiasm.{date}/{sample}.a_ctg.fa.gz'
    conda: 'tabix.yml'
    shell:
        'awk \'/^S/{{print ">"$2"\\n"$3}}\' {input.primary} | fold | bgzip -c > {output.primary} && '
        'awk \'/^S/{{print ">"$2"\\n"$3}}\' {input.alternate} | fold | bgzip -c > {output.alternate} &&'
        # 'awk \'/^S/{{print ">"$2"\\n"$3}}\' {input.hap1} | fold | bgzip -c > {output.hap1} && '
        # 'awk \'/^S/{{print ">"$2"\\n"$3}}\' {input.hap2} | fold | bgzip -c > {output.hap2}'