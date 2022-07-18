# target files
# for purged assembly, "{species}/working/{sample}.{assembler}.{date}/mito-{purge_dir}/mitohifi-v2.done"
# for reads, "{species}/working/{sample}.mitohifi-v2.{date}/mitohifi-v2.done"

localrules: combine_purged_and_htigs, init_mito_reads, \
            replace_mt_in_asm, find_mito_reference_reads, \
            find_mito_reference_asm, init_mito_asm_nopb, check_mito, rotfix_final_mt, replace_mt_in_curated

rule combine_purged_and_htigs:
    input:
        "{species}/working/{sample}.{assembler}.{date}/{purge_dir}/purged.fa.gz",
        "{species}/working/{sample}.{assembler}.{date}/{purge_dir}/purged.htigs.fa.gz"
    output:
        "{species}/working/{sample}.{assembler}.{date}/mito-{purge_dir}/purged_and_htigs.fa"
    shell:
        "zcat {input} > {output}"

rule find_mito_reference_asm:
    input: 
        "{species}/working/{sample}.{assembler}.{date}/mito-{purge_dir}/purged_and_htigs.fa"
    output:
        touch("{species}/working/{sample}.{assembler}.{date}/mito-{purge_dir}/find_mito_ref.done")
    params:
        species="{species}",
        email=config["email"],
        outfolder="{species}/working/{sample}.{assembler}.{date}/mito-{purge_dir}/"
    singularity:
        "scripts/mitohifi-v2.sif"
    shell:
        "cd {params.outfolder} && "
        "findMitoReference.py --species {params.species} --email {params.email} "
        "--outfolder ./ --min_length 14500"

# ideally, we should use dynamic output names from find_mito_reference_asm
rule mitohifi_asm:
    input:
        asm="{species}/working/{sample}.{assembler}.{date}/mito-{purge_dir}/purged_and_htigs.fa",
        ref="{species}/working/{sample}.{assembler}.{date}/mito-{purge_dir}/find_mito_ref.done"
    output:
        "{species}/working/{sample}.{assembler}.{date}/mito-{purge_dir}/final_mitogenome.fasta",
        "{species}/working/{sample}.{assembler}.{date}/mito-{purge_dir}/final_mitogenome.gb",
    params:
        outfolder="{species}/working/{sample}.{assembler}.{date}/mito-{purge_dir}/",
        asm="purged_and_htigs.fa"
    singularity:
        "scripts/mitohifi-v2.sif"
    resources:
        mem_mb=8000
    threads: 4
    shell:
        "cd {params.outfolder} && "
        "mitohifi_v2.py -c {params.asm} "
        "-f *.fasta -g *.gb -t {threads} -o 5"

rule replace_mt_in_asm:
    input:
        asm="{species}/working/{sample}.{assembler}.{date}/mito-{purge_dir}/purged_and_htigs.fa",
        mt="{species}/working/{sample}.{assembler}.{date}/mito-{purge_dir}/final_mitogenome.fasta",
        check="{species}/working/{sample}.{assembler}.{date}/mito-{purge_dir}/final_mitogenome.check.txt"
    output:
        "{species}/working/{sample}.{assembler}.{date}/mito-{purge_dir}/purged_and_htigs_and_mito.fasta"
    conda: 
        "mitohifi_post.yml"
    shell:
        "ctg_name=$(grep '>' {input.mt} | sed -e 's/^>//' -e 's/_rotated$//' -e 's/_rc$//') && " 
        "echo $ctg_name && "
        "faFilter -v name=$ctg_name {input.asm} /dev/stdout | "
        "cat /dev/stdin {input.mt} > {output}"

# cut purge_e assemblies by Ns
# rule replace_mt_in_asm_cut_by_n:
#     input:
#         primary="{species}/working/{sample}.{assembler}.{date}/{purge_dir}/purged.fa.gz",
#         htigs="{species}/working/{sample}.{assembler}.{date}/{purge_dir}/purged.htigs.fa.gz",
#         mt="{species}/working/{sample}.{assembler}.{date}/mito_{purge_dir}/final_mitogenome.fasta"


rule replace_mt_with_reads:
    input:
         asm="{species}/working/{sample}.{assembler}.{date}/mito-{purge_dir}/purged_and_htigs.fa",
         asm_mt="{species}/working/{sample}.{assembler}.{date}/mito-{purge_dir}/final_mitogenome.fasta",
         read_mt="{species}/working/{sample}.mitohifi-v2.{date}/final_mitogenome.fasta",
         asm_check="{species}/working/{sample}.{assembler}.{date}/mito-{purge_dir}/final_mitogenome.check.txt",
         read_check="{species}/working/{sample}.mitohifi-v2.{date}/final_mitogenome.check.txt"
    output:
        "{species}/working/{sample}.{assembler}.{date}/mito-{purge_dir}/purged_and_htigs_and_mito_from_reads.fasta"
    conda: 
        "mitohifi_post.yml"
    shell:
        "ctg_name=$(grep '>' {input.asm_mt} | sed -e 's/^>//' -e 's/_rotated$//' -e 's/_rc$//') && " 
        "echo $ctg_name && "
        "faFilter -v name=$ctg_name {input.asm} /dev/stdout | "
        "cat /dev/stdin {input.read_mt} > {output}"

rule replace_mt_with_hicanu:
    input:
        asm="{species}/working/{sample}.{assembler}.{date}/mito-{purge_dir}/purged_and_htigs.fa",
        asm_mt="{species}/working/{sample}.{assembler}.{date}/mito-{purge_dir}/final_mitogenome.fasta",
        hicanu_mt="{species}/working/{sample}.hicanu.20210327/mito-purging/final_mitogenome.fasta",
        asm_check="{species}/working/{sample}.{assembler}.{date}/mito-{purge_dir}/final_mitogenome.check.txt",
        hicanu_check="{species}/working/{sample}.hicanu.20210327/mito-purging/final_mitogenome.check.txt"
    output:
        "{species}/working/{sample}.{assembler}.{date}/mito-{purge_dir}/purged_and_htigs_and_mito_from_hicanu.fasta"
    conda: 
        "mitohifi_post.yml"
    shell:
        "ctg_name=$(grep '>' {input.asm_mt} | sed -e 's/^>//' -e 's/_rotated$//' -e 's/_rc$//') && " 
        "echo $ctg_name && "
        "faFilter -v name=$ctg_name {input.asm} /dev/stdout | "
        "cat /dev/stdin {input.hicanu_mt} > {output}"

# singularity cannot mkdir, workaround
rule init_mito_reads:
    input:
        "{species}/working/{sample}.hicanu.20210327/{sample}.trimmedReads.fasta.gz"
    output:
        "{species}/working/{sample}.mitohifi-v2.{date}/init.done"
    shell:
        "touch {output}"

rule find_mito_reference_reads:
    input: 
        "{species}/working/{sample}.mitohifi-v2.{date}/init.done"
    output:
        touch("{species}/working/{sample}.mitohifi-v2.{date}/find_mito_ref.done")
    params:
        species="{species}",
        email=config["email"],
        outfolder="{species}/working/{sample}.mitohifi-v2.{date}/"
    singularity:
        "scripts/mitohifi-v2.sif"
    shell:
        "cd {params.outfolder} && "
        "findMitoReference.py --species {params.species} --email {params.email} "
        "--outfolder ./ --min_length 14500"

# ideally, we should use dynamic output names from find_mito_reference_reads
# and non-hard-coded name for trimmed reads directory
# note date mismatch between hicanu and reads
rule mitohifi_reads:
    input:
        reads="{species}/working/{sample}.hicanu.20210327/{sample}.trimmedReads.fasta",
        ref="{species}/working/{sample}.mitohifi-v2.{date}/find_mito_ref.done"
    output:
        "{species}/working/{sample}.mitohifi-v2.{date}/final_mitogenome.fasta"
    params:
        outfolder="{species}/working/{sample}.mitohifi-v2.{date}/",
        reads="../{sample}.hicanu.20210327/{sample}.trimmedReads.fasta"
    singularity:
        "scripts/mitohifi-v2.sif"
    resources:
        mem_mb=8000
    threads: 4
    shell:
        "cd {params.outfolder} && "
        "mitohifi_v2.py -r {params.reads} "
        "-f *.fasta -g *.gb -t {threads} -o 5"

rule check_mito:
    input: "{directory}/final_mitogenome.fasta"
    output: "{directory}/final_mitogenome.check.txt"
    params:
        query="{directory}/final_mitogenome.gb",
        # bash regex matching expecting one reference
        ref="{directory}/[A-Z]*.[0-9].gb"
    conda: "../scripts/check_mito.yml"
    shell:
        "python3 badass_smk/scripts/check_mito.py {params.query} {params.ref} > {output}"



# singularity cannot mkdir, workaround
rule init_mito_asm_nopb:
    input:
        "{species}/working/{sample}.tshea.20210511/assembly.scaffolds.fasta"
    output:
        "{species}/working/{sample}.tshea.20210511/mitohifi/init.done"
    shell:
        "touch {output}"

rule find_mito_reference_asm_nopb:
    input: 
        "{species}/working/{sample}.tshea.20210511/mitohifi/init.done"
    output:
        touch("{species}/working/{sample}.tshea.20210511/mitohifi/find_mito_ref.done")
    params:
        species="{species}",
        email=config["email"],
        outfolder="{species}/working/{sample}.tshea.20210511/mitohifi/"
    singularity:
        "scripts/mitohifi-v2.sif"
    shell:
        "cd {params.outfolder} && "
        "findMitoReference.py --species {params.species} --email {params.email} "
        "--outfolder ./ --min_length 14500"

# ideally, we should use dynamic output names from find_mito_reference_asm
rule mitohifi_asm_nopb:
    input:
        asm="{species}/working/{sample}.tshea.20210511/assembly.scaffolds.fasta",
        ref="{species}/working/{sample}.tshea.20210511/mitohifi/find_mito_ref.done"
    output:
        "{species}/working/{sample}.tshea.20210511/mitohifi/final_mitogenome.fasta",
        "{species}/working/{sample}.tshea.20210511/mitohifi/final_mitogenome.gb",
    params:
        outfolder="{species}/working/{sample}.tshea.20210511/mitohifi/",
        asm="../assembly.scaffolds.fasta"
    singularity:
        "scripts/mitohifi-v2.sif"
    resources:
        mem_mb=8000
    threads: 4
    shell:
        "cd {params.outfolder} && "
        "mitohifi_v2.py -c {params.asm} "
        "-f *.fasta -g *.gb -t {threads} -o 5"

# rotate mt to start with tRNA-Ile
# pick up start point and complement status from .gb
# NB does not work for alternate mt contigs
rule rotfix_final_mt:
    input:
        fa="{prefix}/final_mitogenome.fasta",
        gb="{prefix}/final_mitogenome.gb"
    output:
        "{prefix}/final_mitogenome_rotfix.fa"
    conda:
        "rotate.yml"
    shell:
        """
        ile_coord=$(grep -B1 "Ile" {input.gb} | grep "gene")
        echo $ile_coord
        if [[ $ile_coord == *"complement"* ]]; then
            start=$(echo $ile_coord | cut -d '.' -f3 | cut -d ')' -f1)
            echo Rotating with rc from $start
            python scripts/rotate.py -c -r $start -i {input.fa} > {output}
        elif [[ $ile_coord == *".."* ]]; then
            start=$(echo $ile_coord | rev | cut -d ' ' -f1 | rev | cut -d '.' -f1)
            echo Rotating from $start
            python scripts/rotate.py -r $start -i {input.fa} > {output}
        else
            exit 1
        fi
        """

# given final mitogenome rotfix and a curated assembly
# replace scaff_MT in curated assembly with rotfix final mitogenome
# NB introduced check for bak file to avoid overwriting
# NB lots of directories created in output - workaround for differences in curated dir naming
rule replace_mt_in_curated:
    input:
        asm="{species}/assembly/curated/{asm_prefix}.curated_primary.fa",
        mt="{species}/working/{sample}.{assembler}.{date}/mito-{purge_dir}/final_mitogenome_rotfix.fa"
    output:
        touch("{species}/working/{sample}.{assembler}.{date}/mito-{purge_dir}/{asm_prefix}.replace_curated.done")
    params:
        asm_bak="{species}/assembly/curated/{asm_prefix}.curated_primary.fa.bak",
        reheader_mt="{species}/working/{sample}.{assembler}.{date}/mito-{purge_dir}/final_mitogenome_rotfix.fa.tmp"
    conda:
        "mitohifi_post.yml"
    shell:
        """
        if [[ -f {params.asm_bak} ]]; then
            echo "assembly backup exists, exiting"
            exit 1
        fi

        echo '>scaffold_MT' > {params.reheader_mt}
        tail -n +2 {input.mt} >> {params.reheader_mt}

        mv {input.asm} {params.asm_bak}
        faFilter -v name='scaffold_MT' {params.asm_bak} /dev/stdout | cat /dev/stdin {params.reheader_mt} > {input.asm}

        rm {params.reheader_mt}
        """

