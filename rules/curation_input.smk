import yaml
import os

localrules: generate_yaml_polished

asm_dir="/lustre/scratch116/tol/teams/lawniczak/data/badass/"
        

rule generate_yaml_polished:
    input:
        kmc="{species}/working/{sample}.{assembler}.{date}/scaff_polished.{purge_dir}.hic.{hic_sample}/out.break.salsa2/kmc/{sample}.scaff_polished.ccs.k21.kat.png",
        busco="{species}/working/{sample}.{assembler}.{date}/scaff_polished.{purge_dir}.hic.{hic_sample}/out.break.salsa2/busco5/busco.done",
        scaff="{species}/working/{sample}.{assembler}.{date}/scaff_polished.{purge_dir}.hic.{hic_sample}/out.break.salsa2/scaffolds_FINAL.fasta"
    output:
        "{species}/working/{sample}.{assembler}.{date}/scaff_polished.{purge_dir}.hic.{hic_sample}/{sample}.draft.yaml"
    params:
        species="{species}",
        specimen="{sample}",
        projects=["badass"],
        primary=asm_dir+"{species}/working/{sample}.{assembler}.{date}/scaff_polished.{purge_dir}.hic.{hic_sample}/out.break.salsa2/scaffolds_FINAL.fasta",
        haplotigs=asm_dir+"{species}/working/{sample}.{assembler}.{date}/polished-{purge_dir}/haplotigs.fasta",
        mito=asm_dir+"{species}/working/{sample}.{assembler}.{date}/mito-{purge_dir}/final_mitogenome.fasta",
        agp=asm_dir+"{species}/working/{sample}.{assembler}.{date}/scaff_polished.{purge_dir}.hic.{hic_sample}/out.break.salsa2/scaffolds_FINAL.agp",
        hic=asm_dir+"{species}/working/{sample}.{assembler}.{date}/scaff_polished.{purge_dir}.hic.{hic_sample}/out.break.salsa2/salsa_scaffolds.hic",
        pretext=asm_dir+"{species}/working/{sample}.{assembler}.{date}/scaff_polished.{purge_dir}.hic.{hic_sample}/out.break.salsa2/scaffolds_FINAL.pretext",
        hic_map_img="",
        kmer_spectra_img=asm_dir+"{species}/working/{sample}.{assembler}.{date}/scaff_polished.{purge_dir}.hic.{hic_sample}/out.break.salsa2/kmc/{sample}.scaff_polished.ccs.k21.stacked.kat.png",
        pipeline=['hifiasm (version 0.14)','purge_dups (version 1.2.3)','longranger (version 2.2.2)',
                  'freebayes (v1.3.1)','MitoHiFi(v2)','salsa (v2.2-4c80ac1)'],
        notes=""
    run:
        assert os.path.isfile(params.primary)
        assert os.path.isfile(params.haplotigs)
        # skipping mito check to be validated manually
        if not os.path.isfile(params.mito):
            params.mito=''
        assert os.path.isfile(params.agp)
        assert os.path.isfile(params.hic)
        assert os.path.isfile(params.pretext)
        assert os.path.isfile(params.kmer_spectra_img)

        params.species = params.species.replace("_"," ")

        with open(output[0], 'w') as outfile:
            yaml.dump(dict(params), outfile, sort_keys=False, indent=2)
            # manually add stats
            outfile.write("stats: |\n")
            with open(input.busco[:-10] + "short_summary.specific.diptera_odb10...txt") as infile:
                for line in infile:
                    line=line.strip()
                    if line.startswith('C'):
                        outfile.write(line + '\n')
            with open(input.scaff + ".stats") as infile:
                for line in infile:
                    outfile.write(line)


