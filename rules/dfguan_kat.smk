localrules: generate_fofn

rule generate_fofn:
    input:
        "{species}/working/{sample}.{assembler}.{date}/{purge_dir}/purged.fa"
    output:
        "{species}/working/{sample}.{assembler}.{date}/{purge_dir}/purged.10x.kat/10x.fofn"
    params:
        reads_dir="{species}/genomic_data/{sample}/10x/"
    shell:
        "python badass_smk/scripts/10x_fofn.py {params.reads_dir} > {output}"

# for 10x, use kk16 wrapper
rule dfguan_kat_10x:
    input:
        genome="{species}/working/{sample}.{assembler}.{date}/{purge_dir}/purged.fa",
        fofn="{species}/working/{sample}.{assembler}.{date}/{purge_dir}/purged.10x.kat/10x.fofn"
    output:
        touch("{species}/working/{sample}.{assembler}.{date}/{purge_dir}/purged.10x.kat/dfguan_kat.done")
    params:
        outfolder="{species}/working/{sample}.{assembler}.{date}/{purge_dir}/purged.10x.kat/"
    conda: "dfguan_kat.yml"
    resources:
        mem_mb=20000
    threads: 4
    shell:
        "scripts/kmc.py -p {input.genome} -r {input.fofn} -o {params.outfolder}"

# for ccs and scaffolds, use am60 drafts
rule dfguan_kmc_count_scaffolds:
    input: "{species}/working/{sample}.{assembler}.{date}/{scaff_dir}.{purge_dir}.hic.{hic_sample}/{salsa_dir}/scaffolds_FINAL.fasta"
    output:
        "{species}/working/{sample}.{assembler}.{date}/{scaff_dir}.{purge_dir}.hic.{hic_sample}/{salsa_dir}/kmc/{sample}.{scaff_dir}.k21.kmc_pre",
        "{species}/working/{sample}.{assembler}.{date}/{scaff_dir}.{purge_dir}.hic.{hic_sample}/{salsa_dir}/kmc/{sample}.{scaff_dir}.k21.kmc_suf",
    params:
        kmc_bin="/nfs/users/nfs_d/dg30/luster_dg30/pub/kmc_bins/bin/kmc",
        kmer_size=21,
        output_prefix="{species}/working/{sample}.{assembler}.{date}/{scaff_dir}.{purge_dir}.hic.{hic_sample}/{salsa_dir}/kmc/{sample}.{scaff_dir}.k21",
        work_dir="{species}/working/{sample}.{assembler}.{date}/{scaff_dir}.{purge_dir}.hic.{hic_sample}/{salsa_dir}/kmc/tmp.{sample}.{scaff_dir}.k21",
        mem_gb=16 # don"t forget to link to resources.mem_mb
    resources:
        mem_mb=16 * 1024
    conda: "dfguan_kat.yml" 
    threads: 4
    shell:
        "mkdir -p {params.work_dir} && "
        "{params.kmc_bin} -k{params.kmer_size} -t{threads} -m{params.mem_gb} -ci0 -cs10000 -fm "
        "-sm {input} {params.output_prefix} {params.work_dir} && "
        "rm -r {params.work_dir}"

rule dfguan_kmc_count_ccs_reads:
    input: "{species}/genomic_data/{sample}/pacbio/IRODS.{sample}.fofn"
    output:
        "{species}/genomic_data/{sample}/pacbio/kmc/{sample}.ccs.k21.kmc_pre",
        "{species}/genomic_data/{sample}/pacbio/kmc/{sample}.ccs.k21.kmc_suf"
    params:
        kmc_bin="/nfs/users/nfs_d/dg30/luster_dg30/pub/kmc_bins/bin/kmc",
        kmer_size=21,
        out_prefix="{species}/genomic_data/{sample}/pacbio/kmc/{sample}.ccs.k21",
        work_dir="{species}/genomic_data/{sample}/pacbio/kmc/tmp_k21", # unique for each task to avoid same-file writing SegFault
        fa_mask="{species}/genomic_data/{sample}/pacbio/fasta/*.fasta.gz",
        fofn="{species}/genomic_data/{sample}/pacbio/kmc/tmp_k21/fasta.fofn",
        mem_gb=12 # don"t forget to link to resources.mem_mb
    resources:
        mem_mb=12 * 1024
    conda: "dfguan_kat.yml" 
    threads: 4
    shell:
        # The lower (-ci) and upper (-cs) bounds exclude k-mers with counts outside these boundaries
        "mkdir -p {params.work_dir} && "
        "find {params.fa_mask} -type f > {params.fofn} && "
        "{params.kmc_bin} -k{params.kmer_size} -t{threads} -m{params.mem_gb} -ci0 -cs10000 "
        "-sm -fa @{params.fofn} {params.out_prefix} {params.work_dir} && "
        "rm -r {params.work_dir}"


rule dfguan_kmc_analyse_scaffolds_ccs:
    input: 
        assembly="{species}/working/{sample}.{assembler}.{date}/{scaff_dir}.{purge_dir}.hic.{hic_sample}/{salsa_dir}/kmc/{sample}.{scaff_dir}.k21.kmc_suf",
        reads="{species}/genomic_data/{sample}/pacbio/kmc/{sample}.ccs.k21.kmc_suf"
    output:
        matrix="{species}/working/{sample}.{assembler}.{date}/{scaff_dir}.{purge_dir}.hic.{hic_sample}/{salsa_dir}/kmc/{sample}.{scaff_dir}.ccs.k21.mx"
    params:
        kmc_bin="/nfs/users/nfs_d/dg30/luster_dg30/pub/kmc_bins/bin/kmc_tools",
        assembly_prefix="{species}/working/{sample}.{assembler}.{date}/{scaff_dir}.{purge_dir}.hic.{hic_sample}/{salsa_dir}/kmc/{sample}.{scaff_dir}.k21",
        reads_prefix="{species}/genomic_data/{sample}/pacbio/kmc/{sample}.ccs.k21"
    resources:
        mem_mb=4096
    conda: "dfguan_kat.yml"
    shell:
        "{params.kmc_bin} analyze {params.reads_prefix} {params.assembly_prefix} {output.matrix}"

rule dfguan_kat_plot:
    input:
        matrix="{species}/working/{sample}.{assembler}.{date}/{scaff_dir}.{purge_dir}.hic.{hic_sample}/{salsa_dir}/kmc/{sample}.{scaff_dir}.ccs.k21.mx"
    output: 
        plot="{species}/working/{sample}.{assembler}.{date}/{scaff_dir}.{purge_dir}.hic.{hic_sample}/{salsa_dir}/kmc/{sample}.{scaff_dir}.ccs.k21.kat.png",
        stacked_plot="{species}/working/{sample}.{assembler}.{date}/{scaff_dir}.{purge_dir}.hic.{hic_sample}/{salsa_dir}/kmc/{sample}.{scaff_dir}.ccs.k21.stacked.kat.png"
    conda: "dfguan_kat.yml"
    params:
        spectra="/software/grit/tools/KMC/bin/spectra.py"
    shell:
        "python3 {params.spectra} {input.matrix} {output.plot} && "
        "python3 {params.spectra} -s {input.matrix} {output.stacked_plot} "

