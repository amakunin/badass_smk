import json
import glob
import os

pdir = "/lustre/scratch116/tol/teams/lawniczak/data/badass/"

localrules: primary_contig_names, hic_fofn, fofn_10x, create_wdl_config

rule primary_contig_names:
    input: "{species}/working/{sample}.{assembler}.{date}/{purge_dir}/purged.fa.gz"
    output: "{species}/working/{sample}.{assembler}.{date}/wdl-{purge_dir}/primary.txt"
    shell:
        "zcat {input} | grep '>' | sed 's/>//g' > {output}"

# sort hic crams by name - applicable to gam/fun
# remove secondary alignments
rule collate_hic:
    input: "{species}/genomic_data/{hic_sample}/hic-arima2/IRODS.{hic_sample}.hic-arima2.fofn"
    output: touch("{species}/genomic_data/{hic_sample}/hic-arima2/coord_sorted/collate.done")
    params:
        queue="long",
        working_dir="{species}/genomic_data/{hic_sample}/hic-arima2"
    conda: "samtools.yml"
    threads: 8
    resources: mem_mb=8000
    shell: 
        """
        cd {params.working_dir}

        mkdir -pv coord_sorted

        for f in *.cram
        do
            mv -v $f* coord_sorted/.
            samtools view -h -F 0x800 coord_sorted/$f | samtools collate --threads {threads} -o $f - temp_collate
        done
        """

rule hic_fofn:
    input: "{species}/genomic_data/{hic_sample}/hic-arima2/IRODS.{hic_sample}.hic-arima2.fofn"
    output: "{species}/working/{sample}.{assembler}.{date}/wdl-{purge_dir}/{hic_sample}.hic-arima2.fofn"
    params:
        search_pattern="{species}/genomic_data/{hic_sample}/hic-arima2/*.cram"
    run:
        crams = glob.glob(params.search_pattern)
        with open(output[0], 'w') as outfile:
            for fn in crams:
                outfile.write(os.path.abspath(fn))
                outfile.write('\n')

rule fofn_10x:
    input: "{species}/working/{sample}.{assembler}.{date}/wdl-{purge_dir}/primary.txt"
    output: "{species}/working/{sample}.{assembler}.{date}/wdl-{purge_dir}/{sample}.10x.fofn"
    params:
        reads_dir="{species}/genomic_data/{sample}/10x/"
    shell:
        "python badass_smk/scripts/10x_fofn.py {params.reads_dir} > {output}"


rule create_wdl_config:
    input:
        primary="{species}/working/{sample}.{assembler}.{date}/wdl-{purge_dir}/primary.txt",
        hic_fofn="{species}/working/{sample}.{assembler}.{date}/wdl-{purge_dir}/{hic_sample}.hic-arima2.fofn",
        irods_10x="{species}/genomic_data/{sample}/10x/IRODS.{sample}.10x.fofn",
        comb_ref="{species}/working/{sample}.{assembler}.{date}/mito-{purge_dir}/purged_and_htigs_and_mito.fa",
        fofn_10x="{species}/working/{sample}.{assembler}.{date}/wdl-{purge_dir}/{sample}.10x.fofn"
    output:
        "{species}/working/{sample}.{assembler}.{date}/aps.{purge_dir}.hic.{hic_sample}.json"
    params:
        user = "am60",
        working_dir=pdir+"{species}/working/{sample}.{assembler}.{date}/wdl-{purge_dir}",
        primary=pdir+"{species}/working/{sample}.{assembler}.{date}/wdl-{purge_dir}/primary.txt",
        hic_fofn=pdir+"{species}/working/{sample}.{assembler}.{date}/wdl-{purge_dir}/{hic_sample}.hic-arima2.fofn",
        sample="{sample}",
        fastqs_10x_dir=pdir+"{species}/genomic_data/{sample}/10x/",
        comb_ref="../mito-{purge_dir}/purged_and_htigs_and_mito.fa",
        fofn_10x=pdir+"{species}/working/{sample}.{assembler}.{date}/wdl-{purge_dir}/{sample}.10x.fofn",
        meryl_db=pdir+"{species}/genomic_data/{sample}/10x/preqc/{sample}.10x.k21.meryl/"
    run:
        data = {
            "aps.user":params.user,
            "aps.working_dir":params.working_dir,
            "aps.primary": params.primary,
            "aps.fofn": params.hic_fofn,
            "aps.filterfasta": "/software/team311/mu2/filterfasta.py",
            "aps.sample": params.sample,
            "aps.fastqs": params.fastqs_10x_dir,
            "aps.reference": params.comb_ref,
            "aps.salsa": "/software/tola/installs/tol-workflows/vr-runner/run-salsa",
            "aps.freebayes": "/lustre/scratch116/tol/teams/lawniczak/data/badass/scripts/run_freebayes.py",
            "aps.longranger": "/lustre/scratch116/tol/teams/lawniczak/data/badass/scripts/longranger-2.2.2/longranger-cs/2.2.2/bin/longranger",
            "aps.longranger_align": "/lustre/scratch116/tol/teams/lawniczak/data/badass/scripts/run_longranger_edit.py",
            "aps.arima_motif": "GATC,GANTC,CTNAG,TTAA",
            "aps.busco": "/software/tola/images/busco-5.0.0_cv1.sif",
            "aps.fofn_10x": params.fofn_10x,
            "aps.kmc": "/lustre/scratch116/vr/projects/vgp/user/kk16/bin/kmc.py",
            "aps.busco_lineages": "/lustre/scratch116/tol/resources/busco/v5/lineages/",
            "aps.busco_lib": "diptera_odb10",
            "aps.samtools": "/software/tola/images/htstools-1.11.sif",
            "aps.meryl_db": params.meryl_db,
            "aps.merqury": "/software/tola/images/merqury-1.1.sif"  
        }
        with open(output[0], 'w') as outfile:
            json.dump(data, outfile, indent=4)


# perl environment is not set correctly, switch to run-bsub
rule run_polish_salsa:
    input: "{species}/working/{sample}.{assembler}.{date}/aps.{purge_dir}.hic.{hic_sample}.json"
    output: touch("{species}/working/{sample}.{assembler}.{date}/aps.{purge_dir}.hic.{hic_sample}.done")
    threads: 1
    resources: mem_mb=6000
    # hacky old way of setting queue, see https://github.com/Snakemake-Profiles/lsf/issues/27
    params: queue='long'
    conda: 'wdl_polish_salsa.yml'
    shell:
        "java -Dconfig.file=scripts/Workflow.LSF.aps.conf -jar scripts/cromwell-47.jar "
        "run scripts/Workflow.aps.busco.merqury.farm5.wdl -i {input}"
