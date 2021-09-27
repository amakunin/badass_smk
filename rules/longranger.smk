
rule longranger_make_ref:
    input:
        ref="{species}/working/{sample}.{assembler}.{date}/mito-{purge_dir}/purged_and_htigs_and_mito.fa"
    output:
        touch("{species}/working/{sample}.{assembler}.{date}/wdl-{purge_dir}/make_ref.done"),
        directory("{species}/working/{sample}.{assembler}.{date}/wdl-{purge_dir}/refdata-purged_and_htigs_and_mito")
    params:
        # longranger="../../../../scripts/longranger-2.2.2/longranger-cs/2.2.2/bin/longranger",
        longranger="/software/tola/bin/longranger-2.2.2/longranger-cs/2.2.2/bin/longranger",
        working_dir="{species}/working/{sample}.{assembler}.{date}/wdl-{purge_dir}/",
        ref="../mito-{purge_dir}/purged_and_htigs_and_mito.fa"
    threads: 2
    shell:
        "cd {params.working_dir} && "
        "{params.longranger} mkref {params.ref}"

rule longranger_align:
    input:
        # fofn_10x="{species}/working/{sample}.{assembler}.{date}/wdl-{purge_dir}/{sample}.10x.fofn",
        ref="{species}/working/{sample}.{assembler}.{date}/wdl-{purge_dir}/make_ref.done"
    output:
        bam="{species}/working/{sample}.{assembler}.{date}/wdl-{purge_dir}/{sample}/outs/possorted_bam.bam",
        summary="{species}/working/{sample}.{assembler}.{date}/wdl-{purge_dir}/{sample}/outs/summary.csv"
    params:
        queue="long",
        # longranger="../../../../scripts/longranger-2.2.2/longranger-cs/2.2.2/bin/longranger",
        longranger="/software/tola/bin/longranger-2.2.2/longranger-cs/2.2.2/bin/longranger",
        working_dir="{species}/working/{sample}.{assembler}.{date}/wdl-{purge_dir}/",
        refdata="refdata-purged_and_htigs_and_mito",
        fastqs_10x_dir="../../../genomic_data/{sample}/10x/",
        # fofn_10x="{sample}.10x.fofn",
        sample="{sample}",
        override="/software/tola/bin/longranger-2.2.2/override.json"
    # alternative - run local with 32 cores and 60 gb (test did not use over 4 cores)
    threads: 2
    # resources: mem_mb=60000 # also change in cmd
    shell:
        "cd {params.working_dir} && "
        "rm -rf {params.sample} && "
        "{params.longranger} align --id={params.sample} --fastqs={params.fastqs_10x_dir} "
        "--sample={params.sample} --reference={params.refdata} "
        "--disable-ui --nopreflight --override={params.override} "
        "--jobmode=lsf --maxjobs=1000 --jobinterval=100 "
        "--localcores=2 --localmem=1 "