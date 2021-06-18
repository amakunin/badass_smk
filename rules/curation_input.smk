import glob
import yaml
import os

localrules: draft_assembly_polished_scaffolds, generate_yaml_polished_scaffolds, draft_assembly_scaffolds, generate_yaml_scaffolds


asm_dir="/lustre/scratch116/tol/teams/lawniczak/data/badass/"

rule draft_assembly_polished_scaffolds:
    input:
        kmc="{species}/working/{sample}.{assembler}.{date}/scaff_polished.{purge_dir}.hic.{hic_sample}/out.break.salsa2/kmc/{sample}.scaff_polished.ccs.k21.kat.png",
        busco="{species}/working/{sample}.{assembler}.{date}/scaff_polished.{purge_dir}.hic.{hic_sample}/out.break.salsa2/busco5/busco.done",
        scaff="{species}/working/{sample}.{assembler}.{date}/scaff_polished.{purge_dir}.hic.{hic_sample}/out.break.salsa2/scaffolds_FINAL.fasta",
        postsalsa="{species}/working/{sample}.{assembler}.{date}/scaff_polished.{purge_dir}.hic.{hic_sample}/out.break.salsa2/postsalsa.done",
        pretext_mq0="{species}/working/{sample}.{assembler}.{date}/scaff_polished.{purge_dir}.hic.{hic_sample}/out.break.salsa2/scaffolds_FINAL.mq0.pretext"        
    output:
        touch("{species}/assembly/draft/{sample}.PB.asm1.purge1.polish1.scaff1/{sample}.{assembler}.{date}.{purge_dir}.polished.hic.{hic_sample}.done")
    params:
        in_htigs="{species}/working/{sample}.{assembler}.{date}/polished-{purge_dir}/haplotigs.fasta",
        in_mito="{species}/working/{sample}.{assembler}.{date}/mito-{purge_dir}/final_mitogenome.fasta",
        in_agp="{species}/working/{sample}.{assembler}.{date}/scaff_polished.{purge_dir}.hic.{hic_sample}/out.break.salsa2/scaffolds_FINAL.agp",
        in_hic="{species}/working/{sample}.{assembler}.{date}/scaff_polished.{purge_dir}.hic.{hic_sample}/out.break.salsa2/salsa_scaffolds.hic",
        in_pretext="{species}/working/{sample}.{assembler}.{date}/scaff_polished.{purge_dir}.hic.{hic_sample}/out.break.salsa2/scaffolds_FINAL.pretext",
        in_kat="{species}/working/{sample}.{assembler}.{date}/scaff_polished.{purge_dir}.hic.{hic_sample}/out.break.salsa2/kmc/{sample}.scaff_polished.ccs.k21.stacked.kat.png",
        in_stats="{species}/working/{sample}.{assembler}.{date}/scaff_polished.{purge_dir}.hic.{hic_sample}/out.break.salsa2/scaffolds_FINAL.fasta.stats",
        out_dir="{species}/assembly/draft/{sample}.PB.asm1.purge1.polish1.scaff1",
        out_primary="{species}/assembly/draft/{sample}.PB.asm1.purge1.polish1.scaff1/{sample}.PB.asm1.purge1.polish1.scaff1.scaff_polished.fa.gz",
        out_htigs="{species}/assembly/draft/{sample}.PB.asm1.purge1.polish1.scaff1/{sample}.PB.asm1.purge1.polish1.scaff1.haplotigs.fa.gz",
        out_mito="{species}/assembly/draft/{sample}.PB.asm1.purge1.polish1.scaff1/{sample}.PB.asm1.purge1.polish1.scaff1.mito.fa.gz",
        out_agp="{species}/assembly/draft/{sample}.PB.asm1.purge1.polish1.scaff1/{sample}.PB.asm1.purge1.polish1.scaff1.scaff_polished.agp",
        out_hic="{species}/assembly/draft/{sample}.PB.asm1.purge1.polish1.scaff1/{sample}.PB.asm1.purge1.polish1.scaff1.scaff_polished.hic",
        out_pretext="{species}/assembly/draft/{sample}.PB.asm1.purge1.polish1.scaff1/{sample}.PB.asm1.purge1.polish1.scaff1.scaff_polished.pretext",
        out_pretext_mq0="{species}/assembly/draft/{sample}.PB.asm1.purge1.polish1.scaff1/{sample}.PB.asm1.purge1.polish1.scaff1.scaff_polished.mq0.pretext",
        out_kat="{species}/assembly/draft/{sample}.PB.asm1.purge1.polish1.scaff1/{sample}.PB.asm1.purge1.polish1.scaff1.scaff_polished.ccs.k21.stacked.kat.png",
        out_stats="{species}/assembly/draft/{sample}.PB.asm1.purge1.polish1.scaff1/{sample}.PB.asm1.purge1.polish1.scaff1.scaff_polished.stats"
    conda: 'tabix.yml'
    # copy mt only if exists
    shell:
        "mkdir -p {params.out_dir} && "
        "bgzip -c {input.scaff} > {params.out_primary} && "
        "bgzip -c {params.in_htigs} > {params.out_htigs} && "
        "cp {params.in_agp} {params.out_agp} && "
        "cp {params.in_hic} {params.out_hic} && "
        "cp {params.in_pretext} {params.out_pretext} && "
        "cp {input.pretext_mq0} {params.out_pretext_mq0} && "
        "cp {params.in_kat} {params.out_kat} && "
        "cp {params.in_stats} {params.out_stats} && "
        "[ -f {params.in_mito} ] && bgzip -c {params.in_mito} > {params.out_mito} || true"

rule generate_yaml_polished_scaffolds:
    input:
        "{species}/assembly/draft/{sample}.PB.asm1.purge1.polish1.scaff1/{sample}.{assembler}.{date}.{purge_dir}.polished.hic.{hic_sample}.done",
        busco="{species}/working/{sample}.{assembler}.{date}/scaff_polished.{purge_dir}.hic.{hic_sample}/out.break.salsa2/busco5/busco.done"
    output:
        "{species}/assembly/draft/{sample}.PB.asm1.purge1.polish1.scaff1/{sample}.{assembler}.{date}.{purge_dir}.polished.hic.{hic_sample}.draft.yaml"
    params:
        species="{species}",
        specimen="{sample}",
        projects=["badass"],
        primary=asm_dir+"{species}/assembly/draft/{sample}.PB.asm1.purge1.polish1.scaff1/{sample}.PB.asm1.purge1.polish1.scaff1.scaff_polished.fa.gz",
        haplotigs=asm_dir+"{species}/assembly/draft/{sample}.PB.asm1.purge1.polish1.scaff1/{sample}.PB.asm1.purge1.polish1.scaff1.haplotigs.fa.gz",
        mito=asm_dir+"{species}/assembly/draft/{sample}.PB.asm1.purge1.polish1.scaff1/{sample}.PB.asm1.purge1.polish1.scaff1.mito.fa.gz",
        agp=asm_dir+"{species}/assembly/draft/{sample}.PB.asm1.purge1.polish1.scaff1/{sample}.PB.asm1.purge1.polish1.scaff1.scaff_polished.agp",
        hic=asm_dir+"{species}/assembly/draft/{sample}.PB.asm1.purge1.polish1.scaff1/{sample}.PB.asm1.purge1.polish1.scaff1.scaff_polished.hic",
        pretext=asm_dir+"{species}/assembly/draft/{sample}.PB.asm1.purge1.polish1.scaff1/{sample}.PB.asm1.purge1.polish1.scaff1.scaff_polished.pretext",
        kmer_spectra_img=asm_dir+"{species}/assembly/draft/{sample}.PB.asm1.purge1.polish1.scaff1/{sample}.PB.asm1.purge1.polish1.scaff1.scaff_polished.ccs.k21.stacked.kat.png",
        stats=asm_dir+"{species}/assembly/draft/{sample}.PB.asm1.purge1.polish1.scaff1/{sample}.PB.asm1.purge1.polish1.scaff1.scaff_polished.stats",
        hic_map_img="",
        # reads_pacbio_hifi=asm_dir+"{species}/genomic_data/{sample}/pacbio/fasta/*.filtered.fasta.gz",
        # reads_10x=asm_dir+"{species}/genomic_data/{sample}/10x/*.fastq.gz",
        # reads_hic_arima2=asm_dir+"{species}/genomic_data/{hic_sample}/hic-arima2/*cram",
        pipeline=['hifiasm (version 0.14)','purge_dups (version 1.2.3)','longranger (version 2.2.2)',
                  'freebayes (v1.3.1)','MitoHiFi(v2)','salsa (v2.2-4c80ac1)'],
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

        # params.reads_pacbio_hifi = list(glob.glob(params.reads_pacbio_hifi))
        # params.reads_10x = list(glob.glob(params.reads_10x))
        # params.reads_hic_arima2 = list(glob.glob(params.reads_hic_arima2))

        with open(output[0], 'w') as outfile:
            yaml.dump(dict(params), outfile, sort_keys=False, indent=2)
            # manually add stats
            outfile.write("stats: |\n")
            with open(input.busco[:-10] + "short_summary.specific.diptera_odb10...txt") as infile:
                for line in infile:
                    line=line.strip()
                    if line.startswith('C'):
                        outfile.write(line + '\n')
            with open(params.stats) as infile:
                for line in infile:
                    outfile.write(line)

rule draft_assembly_scaffolds:
    input:
        kmc="{species}/working/{sample}.{assembler}.{date}/scaff.{purge_dir}.hic.{hic_sample}/out.break.salsa2/kmc/{sample}.scaff_polished.ccs.k21.kat.png",
        busco="{species}/working/{sample}.{assembler}.{date}/scaff.{purge_dir}.hic.{hic_sample}/out.break.salsa2/busco5/busco.done",
        scaff="{species}/working/{sample}.{assembler}.{date}/scaff.{purge_dir}.hic.{hic_sample}/out.break.salsa2/scaffolds_FINAL.fasta",
        postsalsa="{species}/working/{sample}.{assembler}.{date}/scaff.{purge_dir}.hic.{hic_sample}/out.break.salsa2/postsalsa.done",
        pretext_mq0="{species}/working/{sample}.{assembler}.{date}/scaff.{purge_dir}.hic.{hic_sample}/out.break.salsa2/scaffolds_FINAL.mq0.pretext"        
    output:
        touch("{species}/assembly/draft/{sample}.PB.asm1.purge1.scaff1/{sample}.{assembler}.{date}.{purge_dir}.polished.hic.{hic_sample}.done")
    params:
        # in_htigs="{species}/working/{sample}.{assembler}.{date}/polished-{purge_dir}/haplotigs.fasta",
        in_mito="{species}/working/{sample}.{assembler}.{date}/mito-{purge_dir}/final_mitogenome.fasta",
        in_agp="{species}/working/{sample}.{assembler}.{date}/scaff.{purge_dir}.hic.{hic_sample}/out.break.salsa2/scaffolds_FINAL.agp",
        in_hic="{species}/working/{sample}.{assembler}.{date}/scaff.{purge_dir}.hic.{hic_sample}/out.break.salsa2/salsa_scaffolds.hic",
        in_pretext="{species}/working/{sample}.{assembler}.{date}/scaff.{purge_dir}.hic.{hic_sample}/out.break.salsa2/scaffolds_FINAL.pretext",
        in_kat="{species}/working/{sample}.{assembler}.{date}/scaff.{purge_dir}.hic.{hic_sample}/out.break.salsa2/kmc/{sample}.scaff_polished.ccs.k21.stacked.kat.png",
        in_stats="{species}/working/{sample}.{assembler}.{date}/scaff.{purge_dir}.hic.{hic_sample}/out.break.salsa2/scaffolds_FINAL.fasta.stats",
        out_dir="{species}/assembly/draft/{sample}.PB.asm1.purge1.scaff1",
        out_primary="{species}/assembly/draft/{sample}.PB.asm1.purge1.scaff1/{sample}.PB.asm1.purge1.scaff1.scaff_polished.fa.gz",
        out_htigs="{species}/assembly/draft/{sample}.PB.asm1.purge1.scaff1/{sample}.PB.asm1.purge1.scaff1.haplotigs.fa.gz",
        out_mito="{species}/assembly/draft/{sample}.PB.asm1.purge1.scaff1/{sample}.PB.asm1.purge1.scaff1.mito.fa.gz",
        out_agp="{species}/assembly/draft/{sample}.PB.asm1.purge1.scaff1/{sample}.PB.asm1.purge1.scaff1.scaff_polished.agp",
        out_hic="{species}/assembly/draft/{sample}.PB.asm1.purge1.scaff1/{sample}.PB.asm1.purge1.scaff1.scaff_polished.hic",
        out_pretext="{species}/assembly/draft/{sample}.PB.asm1.purge1.scaff1/{sample}.PB.asm1.purge1.scaff1.scaff_polished.pretext",
        out_pretext_mq0="{species}/assembly/draft/{sample}.PB.asm1.purge1.scaff1/{sample}.PB.asm1.purge1.scaff1.scaff_polished.mq0.pretext",
        out_kat="{species}/assembly/draft/{sample}.PB.asm1.purge1.scaff1/{sample}.PB.asm1.purge1.scaff1.scaff_polished.ccs.k21.stacked.kat.png",
        out_stats="{species}/assembly/draft/{sample}.PB.asm1.purge1.scaff1/{sample}.PB.asm1.purge1.scaff1.scaff_polished.stats"
    conda: 'tabix.yml'
    # copy mt only if exists
    shell:
        "mkdir -p {params.out_dir} && "
        "bgzip -c {input.scaff} > {params.out_primary} && "
        "bgzip -c {params.in_htigs} > {params.out_htigs} && "
        "cp {params.in_agp} {params.out_agp} && "
        "cp {params.in_hic} {params.out_hic} && "
        "cp {params.in_pretext} {params.out_pretext} && "
        "cp {input.pretext_mq0} {params.out_pretext_mq0} && "
        "cp {params.in_kat} {params.out_kat} && "
        "cp {params.in_stats} {params.out_stats} && "
        "[ -f {params.in_mito} ] && bgzip -c {params.in_mito} > {params.out_mito} || true"

rule generate_yaml_scaffolds:
    input:
        "{species}/assembly/draft/{sample}.PB.asm1.purge1.scaff1/{sample}.{assembler}.{date}.{purge_dir}.polished.hic.{hic_sample}.done",
        busco="{species}/working/{sample}.{assembler}.{date}/scaff.{purge_dir}.hic.{hic_sample}/out.break.salsa2/busco5/busco.done"
    output:
        "{species}/assembly/draft/{sample}.PB.asm1.purge1.scaff1/{sample}.{assembler}.{date}.{purge_dir}.polished.hic.{hic_sample}.draft.yaml"
    params:
        species="{species}",
        specimen="{sample}",
        projects=["badass"],
        primary=asm_dir+"{species}/assembly/draft/{sample}.PB.asm1.purge1.scaff1/{sample}.PB.asm1.purge1.scaff1.scaff_polished.fa.gz",
        haplotigs=asm_dir+"{species}/assembly/draft/{sample}.PB.asm1.purge1.scaff1/{sample}.PB.asm1.purge1.scaff1.haplotigs.fa.gz",
        mito=asm_dir+"{species}/assembly/draft/{sample}.PB.asm1.purge1.scaff1/{sample}.PB.asm1.purge1.scaff1.mito.fa.gz",
        agp=asm_dir+"{species}/assembly/draft/{sample}.PB.asm1.purge1.scaff1/{sample}.PB.asm1.purge1.scaff1.scaff_polished.agp",
        hic=asm_dir+"{species}/assembly/draft/{sample}.PB.asm1.purge1.scaff1/{sample}.PB.asm1.purge1.scaff1.scaff_polished.hic",
        pretext=asm_dir+"{species}/assembly/draft/{sample}.PB.asm1.purge1.scaff1/{sample}.PB.asm1.purge1.scaff1.scaff_polished.pretext",
        kmer_spectra_img=asm_dir+"{species}/assembly/draft/{sample}.PB.asm1.purge1.scaff1/{sample}.PB.asm1.purge1.scaff1.scaff_polished.ccs.k21.stacked.kat.png",
        stats=asm_dir+"{species}/assembly/draft/{sample}.PB.asm1.purge1.scaff1/{sample}.PB.asm1.purge1.scaff1.scaff_polished.stats",
        hic_map_img="",
        # reads_pacbio_hifi=asm_dir+"{species}/genomic_data/{sample}/pacbio/fasta/*.filtered.fasta.gz",
        # reads_hic_arima2=asm_dir+"{species}/genomic_data/{hic_sample}/hic-arima2/*cram",
        pipeline=['hifiasm (version 0.14)','purge_dups (version 1.2.3)','longranger (version 2.2.2)',
                  'freebayes (v1.3.1)','MitoHiFi(v2)','salsa (v2.2-4c80ac1)'],
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

        # params.reads_pacbio_hifi = list(glob.glob(params.reads_pacbio_hifi))
        # params.reads_hic_arima2 = list(glob.glob(params.reads_hic_arima2))

        with open(output[0], 'w') as outfile:
            yaml.dump(dict(params), outfile, sort_keys=False, indent=2)
            # manually add stats
            outfile.write("stats: |\n")
            with open(input.busco[:-10] + "short_summary.specific.diptera_odb10...txt") as infile:
                for line in infile:
                    line=line.strip()
                    if line.startswith('C'):
                        outfile.write(line + '\n')
            with open(params.stats) as infile:
                for line in infile:
                    outfile.write(line)

rule bwa_index:
    input: "{prefix}.fasta"
    output: "{prefix}.fasta.bwt"
    conda: "pretext.yml"
    shell: "bwa index {input}"

rule pretext_mq0:
    input:
        scaff="{species}/working/{sample}.{assembler}.{date}/{scaff_dir}.{purge_dir}.hic.{hic_sample}/out.break.salsa2/scaffolds_FINAL.fasta",
        scaff_index="{species}/working/{sample}.{assembler}.{date}/{scaff_dir}.{purge_dir}.hic.{hic_sample}/out.break.salsa2/scaffolds_FINAL.fasta.bwt",
        hic_fofn="{species}/working/{sample}.{assembler}.{date}/wdl-{purge_dir}/{hic_sample}.hic-arima2.fofn"
    output:
        "{species}/working/{sample}.{assembler}.{date}/{scaff_dir}.{purge_dir}.hic.{hic_sample}/out.break.salsa2/scaffolds_FINAL.mq0.pretext"
    params:
        queue="long"
    conda: "pretext.yml"
    # only align first cram
    threads: 16
    resources: mem_mb=16000
    shell:
        "cram=$(head -n 1 {input.hic_fofn}) && "
        "samtools bam2fq -@ {threads} -0 /dev/null $cram | "
        "bwa mem -t {threads} -p {input.scaff} - | "
        "PretextMap --sortby length --mapq 0 -o {output}"
