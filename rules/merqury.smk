localrules: init_merqury_purged, init_merqury_scaffolds

# singularity cannot mkdir, workaround
rule init_merqury_purged:
    input: "{species}/working/{sample}.{assembler}.{date}/{purge_dir}/purged.fa"
    output: "{species}/working/{sample}.{assembler}.{date}/{purge_dir}/purged.10x.meryl/init.done"
    shell: "touch {output}"

rule merqury_purged:
    input: 
        purged="{species}/working/{sample}.{assembler}.{date}/{purge_dir}/purged.fa",
        htigs="{species}/working/{sample}.{assembler}.{date}/{purge_dir}/purged.htigs.fa",
        flag="{species}/working/{sample}.{assembler}.{date}/{purge_dir}/purged.10x.meryl/init.done"
    output: touch("{species}/working/{sample}.{assembler}.{date}/{purge_dir}/purged.10x.meryl/merqury.done")
    params: 
        outfolder="{species}/working/{sample}.{assembler}.{date}/{purge_dir}/purged.10x.meryl",
        # not checking that these data exist
        meryl_10x="../../../genomic_data/{sample}/10x/preqc/{sample}.10x.k21.meryl/",
        purged="../purged.fa",
        htigs="../purged.htigs.fa",
        prefix="meryl"
    singularity: "scripts/merqury-1.1.sif" 
    resources:
        mem_mb=4000
    threads: 4
    shell:
        "cd {params.outfolder} && "
        "merqury.sh {params.meryl_10x} {params.purged} {params.htigs} {params.prefix}"


# singularity cannot mkdir, workaround
rule init_merqury_scaffolds:
    input: "{species}/working/{sample}.{assembler}.{date}/scaff_polished.{purge_dir}.hic.{hic_sample}/{salsa_dir}/scaffolds_FINAL.fasta"
    output: "{species}/working/{sample}.{assembler}.{date}/scaff_polished.{purge_dir}.hic.{hic_sample}/{salsa_dir}/merqury/init.done"
    shell: "touch {output}"

rule merqury_scaffolds:
    input:
        scaff="{species}/working/{sample}.{assembler}.{date}/scaff_polished.{purge_dir}.hic.{hic_sample}/{salsa_dir}/scaffolds_FINAL.fasta",
        htigs="{species}/working/{sample}.{assembler}.{date}/polished-{purge_dir}/haplotigs.fasta",
        flag="{species}/working/{sample}.{assembler}.{date}/scaff_polished.{purge_dir}.hic.{hic_sample}/{salsa_dir}/merqury/init.done"
    output:
        touch("{species}/working/{sample}.{assembler}.{date}/scaff_polished.{purge_dir}.hic.{hic_sample}/{salsa_dir}/merqury/merqury.done")
    params: 
        outfolder="{species}/working/{sample}.{assembler}.{date}/scaff_polished.{purge_dir}.hic.{hic_sample}/{salsa_dir}/merqury",
        # not checking that these data exist
        meryl_10x="../../../../../genomic_data/{sample}/10x/preqc/{sample}.10x.k21.meryl/",
        purged="../scaffolds_FINAL.fasta",
        htigs="../../../polished-{purge_dir}/haplotigs.fasta",
        prefix="meryl"
    singularity: "scripts/merqury-1.1.sif" 
    resources:
        mem_mb=4000
    threads: 4
    shell:
        "cd {params.outfolder} && "
        "merqury.sh {params.meryl_10x} {params.purged} {params.htigs} {params.prefix}"
 