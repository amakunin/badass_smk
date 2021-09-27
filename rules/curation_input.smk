import glob
import yaml
import os

localrules: draft_assembly_polished_scaffolds, generate_yaml_polished_scaffolds, draft_assembly_scaffolds, generate_yaml_scaffolds, pretext_mq0_to_asm_polished, pretext_mq0_to_asm

asm_dir=os.environ["WD"]

def get_orig_path(path, wildcards):
    """
    Helper renaming function matching sample names requested by JT to original ones
    e.g., idAnoMarsDA-429_01 (in working) to idAnoMarsDA_429_01 (in assembly/draft)
    """
    split_sample = wildcards.sample.split('_', maxsplit=1)
    orig_sample = split_sample[0] + '-' + split_sample[1]
    # expand with all wildcards plus orig_sample
    orig_path = expand(path, orig_sample=orig_sample, **wildcards)

    return orig_path


localrules: test_replace
rule test_replace:
    input:
        lambda wildcards: get_orig_path("{species}/working/{orig_sample}.{assembler}.{date}/scaff_polished.{purge_dir}.hic.{hic_sample}/out.break.salsa2/scaffolds_FINAL.fasta", wildcards),
    output:
        "{species}/assembly/draft/{sample}.{today}/{sample}.{assembler}.{date}.scaff_polished.{purge_dir}.hic.{hic_sample}.test"
    shell:
        "cp {input} {output}"


rule draft_assembly_polished_scaffolds:
    input:
        kmc=lambda wildcards: get_orig_path("{species}/working/{orig_sample}.{assembler}.{date}/"
            "scaff_polished.{purge_dir}.hic.{hic_sample}/out.break.salsa2/kmc/{orig_sample}.scaff_polished.ccs.k21.kat.png", wildcards),
        busco=lambda wildcards: get_orig_path("{species}/working/{orig_sample}.{assembler}.{date}/"
            "scaff_polished.{purge_dir}.hic.{hic_sample}/out.break.salsa2/busco5/busco.done", wildcards),
        scaff=lambda wildcards: get_orig_path("{species}/working/{orig_sample}.{assembler}.{date}/"
            "scaff_polished.{purge_dir}.hic.{hic_sample}/out.break.salsa2/scaffolds_FINAL.fasta", wildcards),
        postsalsa=lambda wildcards: get_orig_path("{species}/working/{orig_sample}.{assembler}.{date}/"
            "scaff_polished.{purge_dir}.hic.{hic_sample}/out.break.salsa2/postsalsa.done", wildcards),
        # pretext_mq0=lambda wildcards: get_orig_path("{species}/working/{orig_sample}.{assembler}.{date}/"
        #     "scaff_polished.{purge_dir}.hic.{hic_sample}/out.break.salsa2/scaffolds_FINAL.mq0.pretext", wildcards)    
    output:
        touch("{species}/assembly/draft/{sample}.{today}/{sample}.{assembler}.{date}.scaff_polished.{purge_dir}.hic.{hic_sample}.done")
    params:
        in_htigs=lambda wildcards: get_orig_path("{species}/working/{orig_sample}.{assembler}.{date}/"
            "polished-{purge_dir}/haplotigs.fasta", wildcards),
        in_mito=lambda wildcards: get_orig_path("{species}/working/{orig_sample}.{assembler}.{date}/"
            "mito-{purge_dir}/final_mitogenome.fasta", wildcards),
        in_agp=lambda wildcards: get_orig_path("{species}/working/{orig_sample}.{assembler}.{date}/"
            "scaff_polished.{purge_dir}.hic.{hic_sample}/out.break.salsa2/scaffolds_FINAL.agp", wildcards),
        in_hic=lambda wildcards: get_orig_path("{species}/working/{orig_sample}.{assembler}.{date}/"
            "scaff_polished.{purge_dir}.hic.{hic_sample}/out.break.salsa2/salsa_scaffolds.hic", wildcards),
        in_pretext=lambda wildcards: get_orig_path("{species}/working/{orig_sample}.{assembler}.{date}/"
            "scaff_polished.{purge_dir}.hic.{hic_sample}/out.break.salsa2/scaffolds_FINAL.pretext", wildcards),
        in_kat=lambda wildcards: get_orig_path("{species}/working/{orig_sample}.{assembler}.{date}/"
            "scaff_polished.{purge_dir}.hic.{hic_sample}/out.break.salsa2/kmc/{orig_sample}.scaff_polished.ccs.k21.stacked.kat.png", wildcards),
        in_stats=lambda wildcards: get_orig_path("{species}/working/{orig_sample}.{assembler}.{date}/"
            "scaff_polished.{purge_dir}.hic.{hic_sample}/out.break.salsa2/scaffolds_FINAL.fasta.stats", wildcards),
        out_dir="{species}/assembly/draft/{sample}.{today}",
        out_primary="{species}/assembly/draft/{sample}.{today}/{sample}.{today}.scaff_polished.fa.gz",
        out_htigs="{species}/assembly/draft/{sample}.{today}/{sample}.{today}.haplotigs.fa.gz",
        out_mito="{species}/assembly/draft/{sample}.{today}/{sample}.{today}.mito.fa.gz",
        out_agp="{species}/assembly/draft/{sample}.{today}/{sample}.{today}.scaff_polished.agp",
        out_hic="{species}/assembly/draft/{sample}.{today}/{sample}.{today}.scaff_polished.hic",
        out_pretext="{species}/assembly/draft/{sample}.{today}/{sample}.{today}.scaff_polished.pretext",
        # out_pretext_mq0="{species}/assembly/draft/{sample}.{today}/{sample}.{today}.scaff_polished.mq0.pretext",
        out_kat="{species}/assembly/draft/{sample}.{today}/{sample}.{today}.scaff_polished.ccs.k21.stacked.kat.png",
        out_stats="{species}/assembly/draft/{sample}.{today}/{sample}.{today}.scaff_polished.stats"
    conda: 'tabix.yml'
    # copy mt only if exists
    shell:
        "mkdir -p {params.out_dir} && "
        "chmod g+w {params.out_dir} && "
        "bgzip -c {input.scaff} > {params.out_primary} && "
        "bgzip -c {params.in_htigs} > {params.out_htigs} && "
        "cp {params.in_agp} {params.out_agp} && "
        "cp {params.in_hic} {params.out_hic} && "
        "cp {params.in_pretext} {params.out_pretext} && "
        # "cp {input.pretext_mq0} {params.out_pretext_mq0} && "
        "cp {params.in_kat} {params.out_kat} && "
        "cp {params.in_stats} {params.out_stats} && "
        "[ -f {params.in_mito} ] && bgzip -c {params.in_mito} > {params.out_mito} || true"

rule generate_yaml_polished_scaffolds:
    input:
        "{species}/assembly/draft/{sample}.{today}/{sample}.{assembler}.{date}.scaff_polished.{purge_dir}.hic.{hic_sample}.done",
        busco=lambda wildcards: get_orig_path("{species}/working/{orig_sample}.{assembler}.{date}/"
            "scaff_polished.{purge_dir}.hic.{hic_sample}/out.break.salsa2/busco5/busco.done", wildcards),
    output:
        "{species}/assembly/draft/{sample}.{today}/{sample}.{assembler}.{date}.scaff_polished.{purge_dir}.hic.{hic_sample}.draft.yaml"
    params:
        species="{species}",
        specimen="{sample}",
        projects=["badass"],
        primary=asm_dir+"{species}/assembly/draft/{sample}.{today}/{sample}.{today}.scaff_polished.fa.gz",
        haplotigs=asm_dir+"{species}/assembly/draft/{sample}.{today}/{sample}.{today}.haplotigs.fa.gz",
        mito=asm_dir+"{species}/assembly/draft/{sample}.{today}/{sample}.{today}.mito.fa.gz",
        agp=asm_dir+"{species}/assembly/draft/{sample}.{today}/{sample}.{today}.scaff_polished.agp",
        hic=asm_dir+"{species}/assembly/draft/{sample}.{today}/{sample}.{today}.scaff_polished.hic",
        pretext=asm_dir+"{species}/assembly/draft/{sample}.{today}/{sample}.{today}.scaff_polished.pretext",
        kmer_spectra_img=asm_dir+"{species}/assembly/draft/{sample}.{today}/{sample}.{today}.scaff_polished.ccs.k21.stacked.kat.png",
        stats=asm_dir+"{species}/assembly/draft/{sample}.{today}/{sample}.{today}.scaff_polished.stats",
        hic_map_img="",
        jira_queue="GRIT",
        reads_pacbio=lambda wildcards: get_orig_path(asm_dir+"{species}/genomic_data/{orig_sample}/pacbio/fasta/*.filtered.fasta.gz", wildcards),
        reads_10x=lambda wildcards: get_orig_path(asm_dir+"{species}/genomic_data/{orig_sample}/10x/*.fastq.gz", wildcards),
        reads_hic=asm_dir+"{species}/genomic_data/{hic_sample}/hic-arima2/*cram",
        pipeline=['hifiasm (version 0.14)','purge_dups (version 1.2.3)','longranger (version 2.2.2)',
                  'freebayes (v1.3.1)','MitoHiFi (v2)','salsa (v2.2-4c80ac1)'],
        notes=""
    run:
        assert os.path.isfile(params.primary), params.primary
        assert os.path.isfile(params.haplotigs), params.haplotigs
        # skipping mito check to be validated manually
        if not os.path.isfile(params.mito):
            params.mito=''
        assert os.path.isfile(params.agp), params.agp
        assert os.path.isfile(params.hic), params.hic
        assert os.path.isfile(params.pretext), params.pretext
        assert os.path.isfile(params.kmer_spectra_img), params.kmer_spectra_img
        assert os.path.isfile(params.stats), params.stats

        params.species = params.species.replace("_"," ")

        # expand within get_orig_path creates a list
        params.reads_pacbio = list(glob.glob(params.reads_pacbio[0]))
        params.reads_10x = list(glob.glob(params.reads_10x[0]))
        # wildcard matching creates a string
        params.reads_hic = list(glob.glob(params.reads_hic))

        with open(output[0], 'w') as outfile:
            yaml.dump(dict(params), outfile, sort_keys=False, indent=2)
            # manually add stats
            outfile.write("stats: |\n")
            # expand within get_orig_path creates a list
            # need to replace as workflow tracks flag, not output file
            with open(input.busco[0][:-10] + "short_summary.specific.diptera_odb10...txt") as infile:
                for line in infile:
                    line=line.strip()
                    if line.startswith('C'):
                        # add indent
                        outfile.write('  ' + line + '\n')
            with open(params.stats) as infile:
                for line in infile:
                    outfile.write('  ' + line)

rule draft_assembly_scaffolds:
    input:
        kmc=lambda wildcards: get_orig_path("{species}/working/{orig_sample}.{assembler}.{date}/"
            "scaff.{purge_dir}.hic.{hic_sample}/out.break.salsa2/kmc/{orig_sample}.scaff.ccs.k21.kat.png", wildcards),
        busco=lambda wildcards: get_orig_path("{species}/working/{orig_sample}.{assembler}.{date}/"
            "scaff.{purge_dir}.hic.{hic_sample}/out.break.salsa2/busco5/busco.done", wildcards),
        scaff=lambda wildcards: get_orig_path("{species}/working/{orig_sample}.{assembler}.{date}/"
            "scaff.{purge_dir}.hic.{hic_sample}/out.break.salsa2/scaffolds_FINAL.fasta", wildcards),
        postsalsa=lambda wildcards: get_orig_path("{species}/working/{orig_sample}.{assembler}.{date}/"
            "scaff.{purge_dir}.hic.{hic_sample}/out.break.salsa2/postsalsa.done", wildcards),
        # pretext_mq0=lambda wildcards: get_orig_path("{species}/working/{orig_sample}.{assembler}.{date}/"
        #   "scaff.{purge_dir}.hic.{hic_sample}/out.break.salsa2/scaffolds_FINAL.mq0.pretext", wildcards),        
    output:
        touch("{species}/assembly/draft/{sample}.{today}/{sample}.{assembler}.{date}.scaff.{purge_dir}.hic.{hic_sample}.done")
    params:
        in_htigs=lambda wildcards: get_orig_path("{species}/working/{orig_sample}.{assembler}.{date}/"
            "{purge_dir}/purged.htigs.fa.gz", wildcards),
        in_mito=lambda wildcards: get_orig_path("{species}/working/{orig_sample}.{assembler}.{date}/"
            "mito-{purge_dir}/final_mitogenome.fasta", wildcards),
        in_agp=lambda wildcards: get_orig_path("{species}/working/{orig_sample}.{assembler}.{date}/"
            "scaff.{purge_dir}.hic.{hic_sample}/out.break.salsa2/scaffolds_FINAL.agp", wildcards),
        in_hic=lambda wildcards: get_orig_path("{species}/working/{orig_sample}.{assembler}.{date}/"
            "scaff.{purge_dir}.hic.{hic_sample}/out.break.salsa2/salsa_scaffolds.hic", wildcards),
        in_pretext=lambda wildcards: get_orig_path("{species}/working/{orig_sample}.{assembler}.{date}/"
            "scaff.{purge_dir}.hic.{hic_sample}/out.break.salsa2/scaffolds_FINAL.pretext", wildcards),
        in_kat=lambda wildcards: get_orig_path("{species}/working/{orig_sample}.{assembler}.{date}/"
            "scaff.{purge_dir}.hic.{hic_sample}/out.break.salsa2/kmc/{orig_sample}.scaff.ccs.k21.stacked.kat.png", wildcards),
        in_stats=lambda wildcards: get_orig_path("{species}/working/{orig_sample}.{assembler}.{date}/"
            "scaff.{purge_dir}.hic.{hic_sample}/out.break.salsa2/scaffolds_FINAL.fasta.stats", wildcards),
        out_dir="{species}/assembly/draft/{sample}.{today}",
        out_primary="{species}/assembly/draft/{sample}.{today}/{sample}.{today}.scaff.fa.gz",
        out_htigs="{species}/assembly/draft/{sample}.{today}/{sample}.{today}.haplotigs.fa.gz",
        out_mito="{species}/assembly/draft/{sample}.{today}/{sample}.{today}.mito.fa.gz",
        out_agp="{species}/assembly/draft/{sample}.{today}/{sample}.{today}.scaff.agp",
        out_hic="{species}/assembly/draft/{sample}.{today}/{sample}.{today}.scaff.hic",
        out_pretext="{species}/assembly/draft/{sample}.{today}/{sample}.{today}.scaff.pretext",
        # out_pretext_mq0="{species}/assembly/draft/{sample}.{today}/{sample}.{today}.scaff.mq0.pretext",
        out_kat="{species}/assembly/draft/{sample}.{today}/{sample}.{today}.scaff.ccs.k21.stacked.kat.png",
        out_stats="{species}/assembly/draft/{sample}.{today}/{sample}.{today}.scaff.stats"
    conda: 'tabix.yml'
    # copy mt only if exists
    shell:
        "mkdir -p {params.out_dir} && "
        "bgzip -c {input.scaff} > {params.out_primary} && "
        "cp {params.in_htigs} {params.out_htigs} && "
        "cp {params.in_agp} {params.out_agp} && "
        "cp {params.in_hic} {params.out_hic} && "
        "cp {params.in_pretext} {params.out_pretext} && "
        # "cp {input.pretext_mq0} {params.out_pretext_mq0} && "
        "cp {params.in_kat} {params.out_kat} && "
        "cp {params.in_stats} {params.out_stats} && "
        "[ -f {params.in_mito} ] && bgzip -c {params.in_mito} > {params.out_mito} || true"

rule generate_yaml_scaffolds:
    input:
        "{species}/assembly/draft/{sample}.{today}/{sample}.{assembler}.{date}.scaff.{purge_dir}.hic.{hic_sample}.done",
        busco=lambda wildcards: get_orig_path("{species}/working/{orig_sample}.{assembler}.{date}/"
            "scaff.{purge_dir}.hic.{hic_sample}/out.break.salsa2/busco5/busco.done", wildcards),
    output:
        "{species}/assembly/draft/{sample}.{today}/{sample}.{assembler}.{date}.scaff.{purge_dir}.hic.{hic_sample}.draft.yaml"
    params:
        species="{species}",
        specimen="{sample}",
        projects=["badass"],
        primary=asm_dir+"{species}/assembly/draft/{sample}.{today}/{sample}.{today}.scaff.fa.gz",
        haplotigs=asm_dir+"{species}/assembly/draft/{sample}.{today}/{sample}.{today}.haplotigs.fa.gz",
        mito=asm_dir+"{species}/assembly/draft/{sample}.{today}/{sample}.{today}.mito.fa.gz",
        agp=asm_dir+"{species}/assembly/draft/{sample}.{today}/{sample}.{today}.scaff.agp",
        hic=asm_dir+"{species}/assembly/draft/{sample}.{today}/{sample}.{today}.scaff.hic",
        pretext=asm_dir+"{species}/assembly/draft/{sample}.{today}/{sample}.{today}.scaff.pretext",
        kmer_spectra_img=asm_dir+"{species}/assembly/draft/{sample}.{today}/{sample}.{today}.scaff.ccs.k21.stacked.kat.png",
        stats=asm_dir+"{species}/assembly/draft/{sample}.{today}/{sample}.{today}.scaff.stats",
        hic_map_img="",
        jira_queue="RC",
        reads_pacbio=lambda wildcards: get_orig_path(asm_dir+"{species}/genomic_data/{orig_sample}/pacbio/fasta/*.filtered.fasta.gz", wildcards),
        reads_hic=asm_dir+"{species}/genomic_data/{hic_sample}/hic-arima2/*cram",
        pipeline=['hifiasm (version 0.14)','purge_dups (version 1.2.3)','MitoHiFi (v2)','salsa (v2.2-4c80ac1)'],
        notes=""
    run:
        assert os.path.isfile(params.primary), params.primary
        assert os.path.isfile(params.haplotigs), params.haplotigs
        # skipping mito check to be validated manually
        if not os.path.isfile(params.mito):
            params.mito=''
        assert os.path.isfile(params.agp), params.agp
        assert os.path.isfile(params.hic), params.hic
        assert os.path.isfile(params.pretext), params.pretext
        assert os.path.isfile(params.kmer_spectra_img), params.kmer_spectra_img
        assert os.path.isfile(params.stats), params.stats

        params.species = params.species.replace("_"," ")

        params.reads_pacbio = list(glob.glob(params.reads_pacbio[0]))
        params.reads_hic = list(glob.glob(params.reads_hic[0]))

        with open(output[0], 'w') as outfile:
            yaml.dump(dict(params), outfile, sort_keys=False, indent=2)
            # manually add stats
            outfile.write("stats: |\n")
            with open(input.busco[0][:-10] + "short_summary.specific.diptera_odb10...txt") as infile:
                for line in infile:
                    line=line.strip()
                    if line.startswith('C'):
                        outfile.write(line + '\n')
            with open(params.stats) as infile:
                for line in infile:
                    outfile.write('  ' + line)

rule bwa_index:
    input: "{prefix}.fasta"
    output: "{prefix}.fasta.bwt"
    conda: "pretext.yml"
    shell: "bwa index {input}"

# all crams to a single interleaved fastq file (temp)
rule crams_to_fastq:
    input:
        hic_fofn="{species}/working/{sample}.{assembler}.{date}/wdl-{purge_dir}/{hic_sample}.hic-arima2.fofn"
    output:
        hic_fq=temp("{species}/working/{sample}.{assembler}.{date}/wdl-{purge_dir}/{hic_sample}.fastq")
    params:
        queue="long",
        cram_pattern="{species}/genomic_data/{hic_sample}/hic-arima2/*cram"
    threads: 8
    # resources: mem_mb=8000
    shell:
        "rm -f {output.hic_fq} && "
        "for file in {params.cram_pattern}; do samtools fastq -@{threads} -s /dev/null -0 /dev/null $file >> {output.hic_fq}; done"

rule pretext_mq0:
    input:
        scaff="{species}/working/{sample}.{assembler}.{date}/{scaff_dir}.{purge_dir}.hic.{hic_sample}/out.break.salsa2/scaffolds_FINAL.fasta",
        scaff_bwa_index="{species}/working/{sample}.{assembler}.{date}/{scaff_dir}.{purge_dir}.hic.{hic_sample}/out.break.salsa2/scaffolds_FINAL.fasta.bwt",
        hic_fofn="{species}/working/{sample}.{assembler}.{date}/wdl-{purge_dir}/{hic_sample}.hic-arima2.fofn",
        hic_fq="{species}/working/{sample}.{assembler}.{date}/wdl-{purge_dir}/{hic_sample}.fastq"
    output:
        "{species}/working/{sample}.{assembler}.{date}/{scaff_dir}.{purge_dir}.hic.{hic_sample}/out.break.salsa2/scaffolds_FINAL.mq0.pretext"
    params:
        queue="long",
        cram_pattern="{species}/genomic_data/{hic_sample}/hic-arima2/*cram",
    conda: "pretext.yml"
    threads: 32
    resources: mem_mb=24000
    shell:
        "bwa mem -t {threads} -p {input.scaff} {input.hic_fq} | "
        "PretextMap --sortby length --mapq 0 -o {output}"

rule pretext_mq0_to_asm_polished:
    input:
        pretext_mq0=lambda wildcards: get_orig_path("{species}/working/{orig_sample}.{assembler}.{date}/"
             "scaff_polished.{purge_dir}.hic.{hic_sample}/out.break.salsa2/scaffolds_FINAL.mq0.pretext", wildcards)
    output:
        touch("{species}/assembly/draft/{sample}.{today}/{sample}.{assembler}.{date}.scaff_polished.{purge_dir}.hic.{hic_sample}.mq0.done")
    params:
        out_pretext_mq0="{species}/assembly/draft/{sample}.{today}/{sample}.{today}.scaff_polished.mq0.pretext"
    shell:
        "cp {input.pretext_mq0} {params.out_pretext_mq0}"

rule pretext_mq0_to_asm:
    input:
        pretext_mq0=lambda wildcards: get_orig_path("{species}/working/{orig_sample}.{assembler}.{date}/"
             "scaff.{purge_dir}.hic.{hic_sample}/out.break.salsa2/scaffolds_FINAL.mq0.pretext", wildcards)
    output:
        touch("{species}/assembly/draft/{sample}.{today}/{sample}.{assembler}.{date}.scaff.{purge_dir}.hic.{hic_sample}.mq0.done")
    params:
        out_pretext_mq0="{species}/assembly/draft/{sample}.{today}/{sample}.{today}.scaff.mq0.pretext"
    shell:
        "cp {input.pretext_mq0} {params.out_pretext_mq0}"

