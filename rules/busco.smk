localrules: init_busco_purged, init_busco_scaff, init_busco_transcriptome

# singularity cannot mkdir, workaround
rule init_busco_purged:
    input: "{species}/working/{sample}.{assembler}.{date}/{purge_dir}/purged.fa"
    output: "{species}/working/{sample}.{assembler}.{date}/{purge_dir}/purged.diptera_odb10.busco/init.done"
    shell: "touch {output}"

rule busco_purged:
    input: 
        fa="{species}/working/{sample}.{assembler}.{date}/{purge_dir}/purged.fa",
        flag="{species}/working/{sample}.{assembler}.{date}/{purge_dir}/purged.diptera_odb10.busco/init.done"
    output: touch("{species}/working/{sample}.{assembler}.{date}/{purge_dir}/purged.diptera_odb10.busco/busco.done")
    params: 
        outfolder="{species}/working/{sample}.{assembler}.{date}/{purge_dir}/purged.diptera_odb10.busco",
        fa="../purged.fa",
        lineage="diptera_odb10"
    singularity: "scripts/busco5.sif" 
    resources:
        mem_mb=16000
    threads: 8
    # adding -r to continue analysis (each run is considered existing for an unclear reason)
    shell:
        "cd {params.outfolder} && "
        "busco -i {params.fa} -c {threads} -o . -m genome -l /lineages/diptera_odb10 --offline -r"

# singularity cannot mkdir, workaround
rule init_busco_scaff:
    input: 
        "{species}/working/{sample}.{assembler}.{date}/{scaff_dir}.{purge_dir}.hic.{hic_sample}/out.break.salsa/scaffolds_FINAL.fasta"
    output: 
        "{species}/working/{sample}.{assembler}.{date}/{scaff_dir}.{purge_dir}.hic.{hic_sample}/busco5/init.done"
    shell: "touch {output}"


rule busco_scaff:
    input: 
        fa="{species}/working/{sample}.{assembler}.{date}/{scaff_dir}.{purge_dir}.hic.{hic_sample}/out.break.salsa/scaffolds_FINAL.fasta",
        flag="{species}/working/{sample}.{assembler}.{date}/{scaff_dir}.{purge_dir}.hic.{hic_sample}/busco5/init.done"
    output: touch("{species}/working/{sample}.{assembler}.{date}/{scaff_dir}.{purge_dir}.hic.{hic_sample}/busco5/busco.done")
    params: 
        outfolder="{species}/working/{sample}.{assembler}.{date}/{scaff_dir}.{purge_dir}.hic.{hic_sample}/busco5",
        fa="../out.break.salsa/scaffolds_FINAL.fasta",
        lineage="diptera_odb10"
    singularity: "scripts/busco5.sif" 
    resources:
        mem_mb=16000
    threads: 8
    # adding -r to continue analysis (each run is considered existing for an unclear reason)
    shell:
        "cd {params.outfolder} && "
        "busco -i {params.fa} -c {threads} -o . -m genome -l /lineages/diptera_odb10 --offline -r"

rule init_busco_transcriptome:
    input: 
        "{species}/working/{sample}.trinity/Trinity.fasta"
    output: 
        "{species}/working/{sample}.trinity/busco5/init.done"
    shell: "touch {output}"


rule busco_transcriptome:
    input: 
        fa="{species}/working/{sample}.trinity/Trinity.fasta",
        flag="{species}/working/{sample}.trinity/busco5/init.done"
    output: touch("{species}/working/{sample}.trinity/busco5/busco.done")
    params: 
        outfolder="{species}/working/{sample}.trinity/busco5/",
        fa="../Trinity.fasta",
        lineage="diptera_odb10"
    singularity: "scripts/busco5.sif" 
    resources:
        mem_mb=16000
    threads: 8
    # adding -r to continue analysis (each run is considered existing for an unclear reason)
    shell:
        "cd {params.outfolder} && "
        "busco -i {params.fa} -c {threads} -o . -m transcriptome -l /lineages/diptera_odb10 --offline -r"
