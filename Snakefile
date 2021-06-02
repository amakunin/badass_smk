# from datetime import datetime

configfile: 'badass_smk/config.yaml'

localrules: gunzip_fa, purge_and_qc, mitohifi_asm_and_reads, salsa_scaffold_and_qc, salsa_scaffold_and_qc_polished

# commonly used wildcards should not include dots or slashes
common_constraint = '[^./]+'
wildcard_constraints:
    sample=common_constraint,
    species=common_constraint,
    assembler=common_constraint,
    purge_dir=common_constraint,
    data_dir=common_constraint,
    data_type=common_constraint,
    hic_sample=common_constraint,
    scaff_dir=common_constraint,
    # date - digits only
    date='\d+',
    # today is today
    # today=datetime.today().strftime('%Y%m%d')
    fa_suffix='fa|fasta'

# rule per tool
include: 'rules/mitohifi.smk'
include: 'rules/busco.smk'
include: 'rules/dfguan_kat.smk'
include: 'rules/merqury.smk'
include: 'rules/purge_dups.smk'
include: 'rules/hifiasm.smk'
include: 'rules/wdl_polish_salsa.smk'
include: 'rules/rnaseq_qc.smk'
include: 'rules/salsa.smk'
include: 'rules/longranger.smk'
include: 'rules/polish.smk'
include: 'rules/nucmer_dotplot.smk'
include: 'rules/curation_input.smk'


# rules shared between tools

# multiple tools rely on uncompressed fasta files
# which can be deleted upon execution
rule gunzip_fa:
    input:
        "{fa_prefix}.{fa_suffix}.gz"
    output:
        temp("{fa_prefix}.{fa_suffix}")
    shell:
        "gunzip -c {input} > {output}"

rule cram_to_fastq:
    input: "{species}/{data_dir}/{sample}/{data_type}/IRODS.{sample}.{data_type}.fofn"
    output: 
        R1="{species}/{data_dir}/{sample}/{data_type}/{sample}.R1.fastq.gz",
        R2="{species}/{data_dir}/{sample}/{data_type}/{sample}.R2.fastq.gz"
    params:
        # lazily expect single cram per sample
        # in case of multiple crams will convert only the first one
        cram_pattern="{species}/{data_dir}/{sample}/{data_type}/*.cram",
        temp="{species}/{data_dir}/{sample}/{data_type}/temp.{sample}"
    conda: "rules/samtools.yml"
    shell:
        # preliminary sort in case cram was aligned to reference
        "samtools collate -O -u {params.cram_pattern} {params.temp} | "
        "samtools fastq -1 {output.R1} -2 {output.R2} -n -0 /dev/null -s /dev/null -"

rule purge_and_qc:
    input:
        "{species}/working/{sample}.{assembler}.{date}/{purge_dir}/purged.fa.gz",
        "{species}/working/{sample}.{assembler}.{date}/{purge_dir}/purged.diptera_odb10.busco/busco.done",
        "{species}/working/{sample}.{assembler}.{date}/{purge_dir}/purged.10x.kat/dfguan_kat.done",
        "{species}/working/{sample}.{assembler}.{date}/{purge_dir}/purged.10x.meryl/merqury.done"
    output:
        "{species}/working/{sample}.{assembler}.{date}/{purge_dir}/purge_and_qc.done"
    shell:
        "touch {output}"

rule mitohifi_asm_and_reads:
    input:
        "{species}/working/{sample}.{assembler}.{date}/mito-{purge_dir}/final_mitogenome.fasta",
        "{species}/working/{sample}.{assembler}.{date}/mito-{purge_dir}/purged_and_htigs_and_mito.fasta",
        # current limitation - same date in mitohifi on reads as in original assembly
        "{species}/working/{sample}.mitohifi-v2.{date}/final_mitogenome.fasta"
    output:
        "{species}/working/{sample}.{assembler}.{date}/mito-{purge_dir}/mitohifi_asm_and_reads.done"
    shell:
        "touch {output}"

rule salsa_scaffold_and_qc:
    input:
        # salsa scaffold
        "{species}/working/{sample}.{assembler}.{date}/scaff.{purge_dir}.hic.{hic_sample}.done",
        # # qc
        # "{species}/working/{sample}.{assembler}.{date}/scaff.{purge_dir}.hic.{hic_sample}/out.break.salsa/busco5/busco.done",
        # "{species}/working/{sample}.{assembler}.{date}/scaff.{purge_dir}.hic.{hic_sample}/out.break.salsa/kmc/{sample}.scaff.ccs.k21.kat.png",
        # refine scaffolds
        "{species}/working/{sample}.{assembler}.{date}/scaff.{purge_dir}.hic.{hic_sample}/out.break.salsa2/postsalsa.done",
        # refined qc
        "{species}/working/{sample}.{assembler}.{date}/scaff.{purge_dir}.hic.{hic_sample}/out.break.salsa2/busco5/busco.done",
        "{species}/working/{sample}.{assembler}.{date}/scaff.{purge_dir}.hic.{hic_sample}/out.break.salsa2/kmc/{sample}.scaff.ccs.k21.kat.png"
    output:
        "{species}/working/{sample}.{assembler}.{date}/scaff.{purge_dir}.hic.{hic_sample}.qc.done"
    shell:
        "touch {output}"

rule salsa_scaffold_and_qc_polished:
    input:
        # longranger align 10x
        "{species}/working/{sample}.{assembler}.{date}/wdl-{purge_dir}/make_ref.done",
        "{species}/working/{sample}.{assembler}.{date}/wdl-{purge_dir}/{sample}/outs/possorted_bam.bam",
        # freebayes polish and split primary from htigs
        "{species}/working/{sample}.{assembler}.{date}/polished-{purge_dir}/primary.fasta",
        # salsa scaffold primary
        "{species}/working/{sample}.{assembler}.{date}/scaff_polished.{purge_dir}.hic.{hic_sample}.done",
        # refine scaffolds
        "{species}/working/{sample}.{assembler}.{date}/scaff_polished.{purge_dir}.hic.{hic_sample}/out.break.salsa2/postsalsa.done",
        # qc
        "{species}/working/{sample}.{assembler}.{date}/scaff_polished.{purge_dir}.hic.{hic_sample}/out.break.salsa2/busco5/busco.done",
        "{species}/working/{sample}.{assembler}.{date}/scaff_polished.{purge_dir}.hic.{hic_sample}/out.break.salsa2/kmc/{sample}.scaff_polished.ccs.k21.kat.png",
        "{species}/working/{sample}.{assembler}.{date}/scaff_polished.{purge_dir}.hic.{hic_sample}/out.break.salsa2/merqury/merqury.done",
        # draft curation yaml
        "{species}/working/{sample}.{assembler}.{date}/scaff_polished.{purge_dir}.hic.{hic_sample}/{sample}.draft.yaml"
    output:
        "{species}/working/{sample}.{assembler}.{date}/scaff_polished.{purge_dir}.hic.{hic_sample}.qc.done"
    shell:
        "touch {output}"

