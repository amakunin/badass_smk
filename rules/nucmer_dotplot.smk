rule nucmer_scaff_to_ref:
    input: 
        q="{species}/working/{sample}.{assembler}.{date}/{scaff_dir}.{purge_dir}.hic.{hic_sample}/out.break.salsa2/scaffolds_FINAL.fasta",
        r="{species}/working/{ref}/{ref}.fa"
    output: "{species}/working/dot/{sample}.{assembler}.{date}.{scaff_dir}.{purge_dir}.hic.{hic_sample}.vs.{ref}.delta"
    params: 
        prefix="{species}/working/dot/{sample}.{assembler}.{date}.{scaff_dir}.{purge_dir}.hic.{hic_sample}.vs.{ref}",
        queue="long"
    log: "{species}/working/dot/{sample}.{assembler}.{date}.{scaff_dir}.{purge_dir}.hic.{hic_sample}.vs.{ref}.log"
    conda: "nucmer_dotplot.yml"
    resources:
        mem_mb=12000
    shell:
        "nucmer -p {params.prefix} {input.r} {input.q} 2> {log}"

rule nucmer_scaff_to_scaff:
    input: 
        q="{species}/working/{sample_q}.{assembler}.{date}/{scaff_dir_q}.{purge_dir_q}.hic.{hic_sample_q}/out.break.salsa2/scaffolds_FINAL.fasta",
        r="{species}/working/{sample_r}.{assembler}.{date}/{scaff_dir_r}.{purge_dir_r}.hic.{hic_sample_r}/out.break.salsa2/scaffolds_FINAL.fasta"
    output: "{species}/working/dot/{sample_r}.{assembler}.{date}.{scaff_dir_r}.{purge_dir_r}.hic.{hic_sample_r}/"
            "vs.{sample_q}.{scaff_dir_q}.{purge_dir_q}.hic.{hic_sample_q}.delta"
    params: 
        prefix="{species}/working/dot/{sample_r}.{assembler}.{date}.{scaff_dir_r}.{purge_dir_r}.hic.{hic_sample_r}/"
               "vs.{sample_q}.{scaff_dir_q}.{purge_dir_q}.hic.{hic_sample_q}",
        queue="long"
    log: "{species}/working/dot/{sample_r}.{assembler}.{date}.{scaff_dir_r}.{purge_dir_r}.hic.{hic_sample_r}/"
            "vs.{sample_q}.{scaff_dir_q}.{purge_dir_q}.hic.{hic_sample_q}.log"
    conda: "nucmer_dotplot.yml"
    resources:
        mem_mb=12000
    shell:
        "nucmer -p {params.prefix} {input.r} {input.q} 2> {log}"

rule nucmer_dotprep:
    input: "{prefix}.delta"
    output: touch("{prefix}.dotprep.done")
    conda: "nucmer_dotplot.yml"
    params: 
        dotprep="/nfs/users/nfs_a/am60/install/DotPrep.py",
        prefix="{prefix}"
    shell: 'python {params.dotprep} --delta {input} --out {params.prefix}'