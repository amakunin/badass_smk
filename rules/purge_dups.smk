import os

envvars:
    "USER"

rule purge_dups:
    input:
        contigs="{species}/working/{sample}.{assembler}.{date}/{sample}.p_ctg.fa.gz",
        htigs="{species}/working/{sample}.{assembler}.{date}/{sample}.a_ctg.fa.gz",
        reads="{species}/working/{sample}.hicanu.20210327/{sample}.trimmedReads.fasta.gz"
    output:
        contigs="{species}/working/{sample}.{assembler}.{date}/{purge_dir}/purged.fa.gz",
        htigs="{species}/working/{sample}.{assembler}.{date}/{purge_dir}/purged.htigs.fa.gz"
    params:
        runner = ("/lustre/scratch116/tol/teams/lawniczak/data/badass/scripts/run-purge_dups_e" 
                    if config["purge_no_e"] else
                    "/software/tola/installs/tol-workflows/vr-runner/run-purge_dups"),
        purge_dir="{species}/working/{sample}.{assembler}.{date}/{purge_dir}/",
        cutoffs=config["purge_cutoffs"],
        user=os.environ["USER"],
        contigs="../{sample}.p_ctg.fa.gz",
        htigs="../{sample}.a_ctg.fa.gz",
        reads="../../{sample}.hicanu.20210327/{sample}.trimmedReads.fasta.gz"
    conda: 'purge_dups.yml'
    # need to source env manually - does not work inside smk shell 
    # "source /software/tola/installs/vr-runner/etc/bashrc && "
    # also
    # successful exit code >1  - https://snakemake.readthedocs.io/en/stable/project_info/faq.html#id12
    shell:
        """
        set +e
        mkdir -p {params.purge_dir} 
        cd {params.purge_dir}
        {params.runner} +mail {params.user} +loop 60 +retries 5 +maxjobs 0 \
            -a {params.contigs} -h {params.htigs} -F {params.reads} \
            -x ccs -o . -c {params.cutoffs} 
        exitcode=$? 
        >&2 echo $exitcode
        #cp -v seqs/purged.* .
        if [ $exitcode -eq 1 ]
        then
            exit 1
        else
            exit 0
        fi
        """