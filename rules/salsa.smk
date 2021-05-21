import os

# key output files are protected in case of incomplete run

rule salsa_scaffolding:
    input:
        purged="{species}/working/{sample}.{assembler}.{date}/{purge_dir}/purged.fa",
        hic_fofn="{species}/working/{sample}.{assembler}.{date}/wdl-{purge_dir}/{hic_sample}.hic-arima2.fofn"
    output:
        touch("{species}/working/{sample}.{assembler}.{date}/scaff.{purge_dir}.hic.{hic_sample}.done"),
        protected("{species}/working/{sample}.{assembler}.{date}/scaff.{purge_dir}.hic.{hic_sample}/out.break.salsa/scaffolds_FINAL.fasta"),
        protected("{species}/working/{sample}.{assembler}.{date}/scaff.{purge_dir}.hic.{hic_sample}/out.break.salsa/salsa_scaffolds.hic")
    params:
        queue="long",
        salsa="/software/tola/installs/tol-workflows/vr-runner/run-salsa",
        arima_motif="GATC,GANTC,CTNAG,TTAA",
        working_dir="{species}/working/{sample}.{assembler}.{date}",
        out_dir="scaff.{purge_dir}.hic.{hic_sample}",
        purged="{purge_dir}/purged.fa", 
        hic_fofn="wdl-{purge_dir}/{hic_sample}.hic-arima2.fofn",
        user=os.environ["USER"]
        # filterfasta="/software/team311/mu2/filterfasta.py",
    # conda: 'purge_dups.yml'
    # need to source env manually - does not work inside smk shell 
    # "source /software/tola/installs/vr-runner/etc/bashrc && "
    # also
    # successful exit code >1  - https://snakemake.readthedocs.io/en/stable/project_info/faq.html#id12
    shell:
        """
        set +e
        cd {params.working_dir}
        {params.salsa} +mail {params.user} +loop 60 +retries 5 +maxjobs 0 \
            --break --motif "{params.arima_motif}" --fofn {params.hic_fofn} \
            --ref-fa {params.purged} --outdir {params.out_dir}
        exitcode=$? 
        >&2 echo $exitcode
        if [ $exitcode -eq 111 ] || [ $exitcode -eq 0 ]
        then
            exit 0
        else
            exit 1
        fi
        """ 

rule salsa_scaffolding_polished:
    input:
        polished_primary="{species}/working/{sample}.{assembler}.{date}/polished-{purge_dir}/primary.fasta",
        hic_fofn="{species}/working/{sample}.{assembler}.{date}/wdl-{purge_dir}/{hic_sample}.hic-arima2.fofn"
    output:
        touch("{species}/working/{sample}.{assembler}.{date}/scaff_polished.{purge_dir}.hic.{hic_sample}.done"),
        # hic failed several times for idAnoFuneDA-408_05, removing
        # protected("{species}/working/{sample}.{assembler}.{date}/scaff_polished.{purge_dir}.hic.{hic_sample}/out.break.salsa/salsa_scaffolds.hic")
    # output should not be removed, putting in as log
    log: "{species}/working/{sample}.{assembler}.{date}/scaff_polished.{purge_dir}.hic.{hic_sample}/out.break.salsa/scaffolds_FINAL.fasta",
    params:
        queue="long",
        salsa="/software/tola/installs/tol-workflows/vr-runner/run-salsa",
        arima_motif="GATC,GANTC,CTNAG,TTAA",
        working_dir="{species}/working/{sample}.{assembler}.{date}",
        out_dir="scaff_polished.{purge_dir}.hic.{hic_sample}",
        polished_primary="polished-{purge_dir}/primary.fasta", 
        hic_fofn="wdl-{purge_dir}/{hic_sample}.hic-arima2.fofn",
        user=os.environ["USER"]
        # filterfasta="/software/team311/mu2/filterfasta.py",
    # conda: 'purge_dups.yml'
    # need to source env manually - does not work inside smk shell 
    # "source /software/tola/installs/vr-runner/etc/bashrc && "
    # also
    # successful exit code >1  - https://snakemake.readthedocs.io/en/stable/project_info/faq.html#id12
    shell:
        """
        set +e
        cd {params.working_dir}
        {params.salsa} +mail {params.user} +loop 60 +retries 5 +maxjobs 0 \
            --break --motif "{params.arima_motif}" --fofn {params.hic_fofn} \
            --ref-fa {params.polished_primary} --outdir {params.out_dir}
        exitcode=$? 
        >&2 echo $exitcode
        if [ $exitcode -eq 1 ]
        then
            exit 1
        else
            exit 0
        fi
        """

rule salsa_refine_prep:
    input:
        "{species}/working/{sample}.{assembler}.{date}/{scaff_dir}.{purge_dir}.hic.{hic_sample}/out.break.salsa/scaffolds_FINAL.fasta"
    output:
        touch("{species}/working/{sample}.{assembler}.{date}/{scaff_dir}.{purge_dir}.hic.{hic_sample}/out.break.salsa2/salsa2.prep.done")
    params:
        ref="{species}/working/{sample}.{assembler}.{date}/{scaff_dir}.{purge_dir}.hic.{hic_sample}/out.break.salsa/assembly.cleaned.fasta"
    conda: "samtools.yml"
    shell: "samtools faidx {params.ref}"

rule salsa_refine:
    input:
        "{species}/working/{sample}.{assembler}.{date}/{scaff_dir}.{purge_dir}.hic.{hic_sample}/out.break.salsa/scaffolds_FINAL.fasta",
        "{species}/working/{sample}.{assembler}.{date}/{scaff_dir}.{purge_dir}.hic.{hic_sample}/out.break.salsa2/salsa2.prep.done"
    output:
        touch("{species}/working/{sample}.{assembler}.{date}/{scaff_dir}.{purge_dir}.hic.{hic_sample}/out.break.salsa2/salsa2.done"),
    log: "{species}/working/{sample}.{assembler}.{date}/{scaff_dir}.{purge_dir}.hic.{hic_sample}/out.break.salsa2/scaffolds_FINAL.fasta"
    params:
        # previous implementation unifying salsa_refine and postsalsa_refined
        # exits upon postsalsa submission
        # salsa2="/software/team311/cz3/install/SALSA/run-salsa2",
        arima_motif="GATC,GANTC,CTNAG,TTAA",
        out_dir="{species}/working/{sample}.{assembler}.{date}/{scaff_dir}.{purge_dir}.hic.{hic_sample}/out.break.salsa2/",
        ref="{species}/working/{sample}.{assembler}.{date}/{scaff_dir}.{purge_dir}.hic.{hic_sample}/out.break.salsa/assembly.cleaned.fasta",
        bed="{species}/working/{sample}.{assembler}.{date}/{scaff_dir}.{purge_dir}.hic.{hic_sample}/out.break.salsa/alignment_iteration_1.bed"
    singularity: "/software/team311/cz3/images/salsa2-2.2-4c80ac1.sif"
    shell:
        # "{params.salsa2} -m {params.arima_motif} -o {params.out_dir} {params.ref} {params.bed}"
        "run_pipeline.py -a {params.ref} -l {params.ref}.fai -b {params.bed} "
        "-e {params.arima_motif} -m yes -p yes -o {params.out_dir}"

rule postsalsa_refined:
    input:
        "{species}/working/{sample}.{assembler}.{date}/{scaff_dir}.{purge_dir}.hic.{hic_sample}/out.break.salsa2/scaffolds_FINAL.fasta"
    output:
        touch("{species}/working/{sample}.{assembler}.{date}/{scaff_dir}.{purge_dir}.hic.{hic_sample}/out.break.salsa2/postsalsa.done")
    params:
        queue="long",
        salsa_post="/software/team311/cz3/install/SALSA/run-salsa2-post",
        working_dir="{species}/working/{sample}.{assembler}.{date}/{scaff_dir}.{purge_dir}.hic.{hic_sample}/out.break.salsa2",
        out_dir=".",
        user=os.environ["USER"]
        # filterfasta="/software/team311/mu2/filterfasta.py",
    # conda: 'purge_dups.yml'
    # need to source env manually - does not work inside smk shell 
    # "source /software/tola/installs/vr-runner/etc/bashrc && "
    # also
    # successful exit code >1  - https://snakemake.readthedocs.io/en/stable/project_info/faq.html#id12
    shell:
        """
        set +e
        cd {params.working_dir}
        {params.salsa_post} +mail {params.user} +loop 60 +retries 5 +maxjobs 0 \
            -d {params.out_dir}
        exitcode=$? 
        >&2 echo $exitcode
        if [ $exitcode -eq 111 ] || [ $exitcode -eq 0 ]
        then
            exit 0
        else
            exit 1
        fi
        """