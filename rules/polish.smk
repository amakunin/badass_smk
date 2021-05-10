import os

rule index_ref:
    input:
        "{species}/working/{sample}.{assembler}.{date}/mito-{purge_dir}/purged_and_htigs_and_mito.fa"
    output:
        "{species}/working/{sample}.{assembler}.{date}/mito-{purge_dir}/purged_and_htigs_and_mito.fa.fai"
    conda: "samtools.yml"
    shell:
        "samtools faidx {input}"

rule polish:
    input:
        bam="{species}/working/{sample}.{assembler}.{date}/wdl-{purge_dir}/{sample}/outs/possorted_bam.bam",
        summary="{species}/working/{sample}.{assembler}.{date}/wdl-{purge_dir}/{sample}/outs/summary.csv",
        ref="{species}/working/{sample}.{assembler}.{date}/mito-{purge_dir}/purged_and_htigs_and_mito.fa.fai",
    output:
        consensus="{species}/working/{sample}.{assembler}.{date}/polished-{purge_dir}/consensus.fasta"
    params:
        queue="long",
        working_dir="{species}/working/{sample}.{assembler}.{date}/",
        freebayes="/software/tola/installs/tol-workflows/vr-runner/run-freebayes-consensus",
        user=os.environ["USER"],
        bam="wdl-{purge_dir}/{sample}/outs/possorted_bam.bam",
        summary="wdl-{purge_dir}/{sample}/outs/summary.csv",
        ref="mito-{purge_dir}/purged_and_htigs_and_mito.fa",
        out_dir="polished-{purge_dir}"
    # conda: "purge_dups.yml"
    shell:
        """
        set +e
        cd {params.working_dir} && 
        {params.freebayes} +maxjobs 0 +retries 7 +loop 60 +mail {params.user} \
            -b {params.bam} -s {params.summary} -f {params.ref} -o {params.out_dir}
        exitcode=$? 
        >&2 echo $exitcode
        if [ $exitcode -eq 1 ]
        then
            exit 1
        else
            exit 0
        fi
        """

rule split_polished:
    input:
        consensus="{species}/working/{sample}.{assembler}.{date}/polished-{purge_dir}/consensus.fasta",
        primary_list="{species}/working/{sample}.{assembler}.{date}/wdl-{purge_dir}/primary.txt"
    output:
        primary="{species}/working/{sample}.{assembler}.{date}/polished-{purge_dir}/primary.fasta",
        htigs="{species}/working/{sample}.{assembler}.{date}/polished-{purge_dir}/haplotigs.fasta",
    params:
        filterfasta="/software/team311/mu2/filterfasta.py"
    conda: "wdl_polish_salsa.yml"
    shell:
        "python {params.filterfasta} -i {input.primary_list} {input.consensus} {output.primary} && "
        "python {params.filterfasta} -n -i {input.primary_list} {input.consensus} {output.htigs}"
    