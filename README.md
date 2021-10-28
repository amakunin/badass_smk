# badass_smk

Ad-hoc snakemake recipes for manual assembly steps in BAdAss project.

To use, install snakemake (tested on 5.27.4) and configure an [LSF profile](https://github.com/snakemake-profiles/doc)

## Examples

All commands should be executed from top level directory of the project with snakemake in environment
```
cd /lustre/scratch116/tol/teams/lawniczak/data/badass
conda activate smk
```

### Contig adjustments

Purge dups with changed thresholds, purging within contigs, and subsequent QC (BUSCO, Merqury, KAT)

```
source /software/tola/installs/vr-runner/etc/bashrc
bash badass_smk/submit.sh -n --config purge_cutoffs=\"5,40,140\" purge_no_e=True -p Anopheles_coustani/working/idAnoZiCoDA-A2_x.hifiasm.20210327/purging_40_e/purge_and_qc.done
```

Mitohifi on assembly and reads, create an assembly with fixed mitochondrial genome, and check the mitogenome annotation for obvious problems
```
bash badass_smk/submit.sh -n Anopheles_coustani/working/idAnoCousDA-361_x.hifiasm.20210327/mito-purging/mitohifi_asm_and_reads.done
```

In case mitogenome needs to be substituted from hicanu or reads:
```
bash badass_smk/submit.sh -n Sabethes_cyaneus/working/idSabCyanKW18_F2.hifiasm.20210516/mito-purging/purged_and_htigs_and_mito_from_{hicanu,reads}.fasta
```

Output filenames of mitohifi and input filenames for downstream processing are currently incompatible, so one needs to rename manually before next step of processing
```
cd Anopheles_coustani/working/idAnoCousDA-361_x.hifiasm.20210327/mito-purging/
mv purged_and_htigs_and_mito.fasta purged_and_htigs_and_mito.fa
``` 

Fix mitogenome rotation and orientation to match Anopheles mt in GenBank which start from tRNA-Ile on + strand
```
bash badass_smk/submit.sh Anopheles_funestus/working/idAnoFuneDA-386_01.hicanu.20210327/mito-purging/final_mitogenome_rotfix.fa
```

Rotation issue was discovered after bulk of curation requests were submitted, so we need to manually replace mt in draft directory
```
gzip -c Anopheles_aquasalis/working/idAnoAquaMG-Q_14.hicanu.20210327/mito-purging/final_mitogenome_rotfix.fa > Anopheles_aquasalis/assembly/draft/idAnoAquaMG_Q_14.20210712/idAnoAquaMG_Q_14.20210712.mito.fa.gz
```

TODO guidelines on replacing mt in curated genomes

### Scaffolding

Prior to scaffolding in funestus and gambiae, need to collate crams and remove secondary alignments as those were aligned to reference and sorted by coordinate
```
bash badass_smk/submit.sh -n Anopheles_gambiae/genomic_data/idAnoGambDA-150_06/hic-arima2/coord_sorted/collate.done
```

Prior to scaffolding in purge_e contigs, need to cut by Ns
```
cd Anopheles_funestus/working/idAnoFuneDA-402_03.hifiasm.20210327/mito-purging_e
seqtk cutN -n1 purged_and_htigs_and_mito.fasta | sed 's/:/_/g' | sed 's/-/_/g' > purged_and_htigs_and_mito_cutN.fasta
grep ">ptg" purged_and_htigs_and_mito_cutN.fasta | sed 's/>//g' > ../wdl-purging_e/primary.txt 
cp purged_and_htigs_and_mito_cutN.fasta purged_and_htigs_and_mito.fa
```

Polishing and scaffolding with snakemake (translated from wdl due to difficulties with debugging). Now also includes tweaked salsa scaffolds.
```
source /software/tola/installs/vr-runner/etc/bashrc
bash badass_smk/submit.sh -n Anopheles_funestus/working/idAnoFuneDA-408_06.hifiasm.20210327/scaff_polished.purging.hic.idAnoFuneDA-408_05.qc.done
```

Scaffolding without polishing where 10x data is not available
```
source /software/tola/installs/vr-runner/etc/bashrc
bash badass_smk/submit.sh -n Anopheles_funestus/working/idAnoFuneDA-402_05.hifiasm.20210327/scaff.purging.hic.idAnoFuneDA-402_09.qc.done
```

Re-scaffolding tshea assemblies with fixed Salsa - similar to scaffolding without polishing
```
source /software/tola/installs/vr-runner/etc/bashrc
bash badass_smk/submit.sh -n Anopheles_bellator/working/idAnoBellAS-SP24_06.tshea.20210511/scaff.scaffolds.hic.idAnoBellAS-SP24_05.qc.done
```

### Curation

Curation inputs (including MAPQ0 pretext)
```
bash badass_smk/submit.sh -n Anopheles_funestus/assembly/draft/idAnoFuneDA_408_06.20210712/idAnoFuneDA_408_06.hifiasm.20210327.scaff_polished.purging.hic.idAnoFuneDA-408_05.draft.yaml
```
where `idAnoFuneDA_408_06` is derived from `idAnoFuneDA-408_06` by replacement of dash with underscore, `20210712` is today's date


In cases when mitochondrial genome in hifiasm assembly is not good enough, replace it manually with mitogenome from hicanu assembly or reads
```
gzip -c Anopheles_aquasalis/working/idAnoAquaMG-Q_14.hicanu.20210327/mito-purging/final_mitogenome.fasta > Anopheles_aquasalis/assembly/draft/idAnoAquaMG_Q_14.20210712/idAnoAquaMG_Q_14.20210712.mito.fa.gz
```

For assemblies without PacBio data, mitochondrial contigs should be identified separately
```
badass_smk/submit.sh -n Anopheles_bellator/working/idAnoBellAS-SP24_06.tshea.20210511/mitohifi/final_mitogenome.check.txt
```

### Comparative

Nucmer dotplots vs reference
```
bash badass_smk/submit.sh -n Anopheles_funestus/working/dot/idAnoFuneDA-402_05.hifiasm.20210327.scaff.purging.hic.idAnoFuneDA-402_09.vs.afunf3.dotprep.done
```

### RNAseq data QC

Adapter trimming and alignment to ref
```
bash badass_smk/submit.sh -n Anopheles_aquasalis/working/idAnoAquaMG-Q_16.star.idAnoAquaMG-Q_14.hifiasm.20210327.purging/idAnoAquaMG-Q_16.markdup.stats
```

Transcriptome assembly and busco
```
bash badass_smk/submit.sh -n Anopheles_darlingi/working/idAnoDarlMG-G_01.trinity/busco5/busco.done
```

Compress assembly intermediate files (up to 100k per sample)
```
bash badass_smk/submit.sh -n Anopheles_darlingi/working/idAnoDarlMG-G_01.trinity/read_partitions.tar.gz
```


## Results

see [this table](https://docs.google.com/spreadsheets/d/1q0uxXEVn_n8Z8_hkwaqNoP9qYQbC0cooIZHnpUmZNjM/edit?pli=1#gid=1776822398)


### Re-purge 

Significant proportion of two haplotypes retained in 9 samples, 3 of those have low coverage (below 15x).
Two solutions tested:
- Increase coverage threshold for separating haplotype and diploid peaks
- Enable purging within contigs (remove -e option) - might introduce mis-assemblies, need to split 

Summary
- idAnoZiCoDA-A2_x - 3.5% busco dups in all runs (internal 12Mbp dup, use default)
- idAnoFuneDA-402_03 - 2.0->0.5% busco dups in no-e (use no-e, need N-split in downstream processing)
- idAnoFuneDA-408_07 - 1.5% busco dups in all runs (use default)
- idAnoGambNW-F1_2 - 2.2->1.8% busco dups in no-e (use no-e)
- idAnoMoucSN-F20_07 - 1.3% busco dups in all runs (use default)
- idAnoNiliSN-F31_02 - 1.1% busco dups in all runs (use default)
- Low coverage samples idAnoFuneDA-402_05, idAnoFuneDA-416_05, idAnoGambDA-150_04 - use default purge


### Mitohifi 

Anopheles_aquasalis/working/idAnoAquaMG-Q_14.hifiasm.20210327/mito-purging/final_mitogenome.check.txt
- frameshift indel in ND5, no assembly in reads, tried finding in hicanu assembly 
- TODO check annotation and read alignments 
- hicanu assembly looks good, replacing 

Anopheles_aquasalis/working/idAnoAquaMG-Q_19.hifiasm.20210327/mito-purging/final_mitogenome.check.txt
- OK

Anopheles_coluzzii/working/idAnoColuKW18-c001_1.hifiasm.20210327/mito-purging/final_mitogenome.check.txt
- OK

Anopheles_coluzzii/working/idAnoColuKW18-c001_2.hifiasm.20210327/mito-purging/final_mitogenome.check.txt
- frameshift in COX3, two mt contigs hap_ptg000450l_1_1 (with frameshift) and ptg000420l_1 (without frameshift) - remove both, add only one - purged_and_htigs_and_mito_alt.fasta 

Anopheles_coustani/working/idAnoCousDA-361_x.hifiasm.20210327/mito-purging/final_mitogenome.check.txt
- COX3 22bp shorter in query, lack 3' - re-evaluated using coustani ref, cheked annotation in IGV - no stop - use as is 

Anopheles_coustani/working/idAnoZiCoDA-A2_x.hifiasm.20210327/mito-purging/final_mitogenome.check.txt
- multiple frameshifts in both potential mt contigs hap_ptg000459l_1_1 and ptg000097l_1. Three contigs in reads, among those, ptg000001l has the least problems - two genes shortened, one of those with stop codon 
- reads have similar problems
- run mitohifi on hicanu with coustani reference - looks good, replace two contigs with hicanu

Anopheles_coustani/working/idAnoZiCoDA-A2_x.hifiasm.20210327/mito-purging_40/final_mitogenome.check.txt
- skip, errors as in idAnoZiCoDA-A2_x purging

Anopheles_coustani/working/idAnoZiCoDA-A2_x.hifiasm.20210327/mito-purging_40_e/final_mitogenome.check.txt
- skip, errors as in idAnoZiCoDA-A2_x purging

Anopheles_coustani/working/idAnoZiCoDA-A2_x.hifiasm.20210327/mito-purging_e/final_mitogenome.check.txt
- skip, errors as in idAnoZiCoDA-A2_x purging

Anopheles_darlingi/working/idAnoDarlJC-H15_27.hifiasm.20210327/mito-purging/final_mitogenome.check.txt
- OK

Anopheles_darlingi/working/idAnoDarlMG-H_01.hifiasm.20210327/mito-purging/final_mitogenome.check.txt
- OK

Anopheles_darlingi/working/idAnoDarlMG-H_07.hifiasm.20210327/mito-purging/final_mitogenome.check.txt
- OK

Anopheles_funestus/working/idAnoFuneDA-386_01.hifiasm.20210327/mito-purging/final_mitogenome.check.txt
- no mt in reads or assembly, but there are hits in individual reads. Thus, problem with assembly in mitohifi - mt was assembled by hicanu - add to hifiasm assembly

Anopheles_funestus/working/idAnoFuneDA-386_06.hifiasm.20210327/mito-purging/final_mitogenome.check.txt
- OK

Anopheles_funestus/working/idAnoFuneDA-402_03.hifiasm.20210327/mito-purging/final_mitogenome.check.txt
- skip

Anopheles_funestus/working/idAnoFuneDA-402_03.hifiasm.20210327/mito-purging_30/final_mitogenome.check.txt
- skip

Anopheles_funestus/working/idAnoFuneDA-402_03.hifiasm.20210327/mito-purging_e/final_mitogenome.check.txt
- OK

Anopheles_funestus/working/idAnoFuneDA-402_05.hifiasm.20210327/mito-purging/final_mitogenome.check.txt
- no mt in assembly, mt in reads have multiple shorter CDSs 
- annotation indicates stop codons, 
- tried hicanu mt - both potential contigs look bad
- read+assembly alignments - reads have long deletion and minor issues; hicanu has small indels and SNPs supported by a single read (in both potental contigs)
- use alt mt from hicanu as it seemingly has less errors 

Anopheles_funestus/working/idAnoFuneDA-408_06.hifiasm.20210327/mito-purging/final_mitogenome.check.txt
- assembly mt has 2 CDSs shorter than expected, reads have no mt
- hicanu mt has multiple shortened CDSs, presumably frameshift at 5' ND5 and 3' ND4 
- read+assembly alignments suggest untrue deletions in hifiasm, hicanu looks good, but circularisaion is off 
- replace hifiasm mt with hicanu

Anopheles_funestus/working/idAnoFuneDA-408_07.hifiasm.20210327/mito-purging/final_mitogenome.check.txt
- OK

Anopheles_funestus/working/idAnoFuneDA-408_07.hifiasm.20210327/mito-purging_40/final_mitogenome.check.txt
- skip

Anopheles_funestus/working/idAnoFuneDA-408_07.hifiasm.20210327/mito-purging_e/final_mitogenome.check.txt
- skip

Anopheles_funestus/working/idAnoFuneDA-414_04.hifiasm.20210327/mito-purging/final_mitogenome.check.txt
- OK

Anopheles_funestus/working/idAnoFuneDA-416_04.hifiasm.20210327/mito-purging/final_mitogenome.check.txt
- OK

Anopheles_funestus/working/idAnoFuneDA-416_05.hifiasm.20210327/mito-purging/final_mitogenome.check.txt
- no mt in assembly, mt asm in reads ok - added mt from reads to asm without replacement

Anopheles_gambiae/working/idAnoGambDA-150_04.hifiasm.20210327/mito-purging/
- no mt in assembly, mt asm in reads ok - added mt from reads to asm without replacement

Anopheles_gambiae/working/idAnoGambDA-407_04.hifiasm.20210327/mito-purging/final_mitogenome.check.txt
- many frameshifts in assembly, good assembly in reads - replaced with mito from reads

Anopheles_gambiae/working/idAnoGambDA-407_05.hifiasm.20210327/mito-purging/final_mitogenome.check.txt
- no mt in assembly, good assembly in reads - added mt from reads to asm without replacement

Anopheles_gambiae/working/idAnoGambDA-407_15.hifiasm.20210327/mito-purging/final_mitogenome.check.txt
- no mt in assembly, good assembly in reads - added mt from reads to asm without replacement

Anopheles_gambiae/working/idAnoGambNW-F1_1.hifiasm.20210327/mito-purging/final_mitogenome.check.txt
- two contigs in assembly hap_ptg000342l_1_1  ptg000384l_1, both with frameshifts, good assembly in reads - remove two contigs and add mito from reads

Anopheles_gambiae/working/idAnoGambNW-F1_2.hifiasm.20210327/mito-purging/final_mitogenome.check.txt
- skip

Anopheles_gambiae/working/idAnoGambNW-F1_2.hifiasm.20210327/mito-purging_33/final_mitogenome.check.txt
- skip

Anopheles_gambiae/working/idAnoGambNW-F1_2.hifiasm.20210327/mito-purging_e/final_mitogenome.check.txt
- two contigs in assembly hap_ptg000367l  ptg000329l, both with frameshifts, relatively good top assembly (atg000001l only has shorter ATP8) in reads - remove two contigs and add mito from reads

Anopheles_maculipalpis/working/idAnoMacuDA-375_x.hifiasm.20210327/mito-purging/final_mitogenome.check.txt
- OK

Anopheles_marshallii/working/idAnoMarsDA-429_01.hifiasm.20210327/mito-purging/final_mitogenome.check.txt
- OK

Anopheles_marshallii/working/idAnoMarsDA-429_02.hifiasm.20210327/mito-purging/final_mitogenome.check.txt
- OK

Anopheles_moucheti/working/idAnoMoucSN-F20_07.hifiasm.20210327/mito-purging/final_mitogenome.check.txt
- OK

Anopheles_moucheti/working/idAnoMoucSN-F20_07.hifiasm.20210327/mito-purging_35/final_mitogenome.check.txt
- skip

Anopheles_moucheti/working/idAnoMoucSN-F20_07.hifiasm.20210327/mito-purging_e/final_mitogenome.check.txt
- skip

Anopheles_nili/working/idAnoNiliSN-F31_02.hifiasm.20210327/mito-purging/final_mitogenome.check.txt
- ok

Anopheles_nili/working/idAnoNiliSN-F31_02.hifiasm.20210327/mito-purging_36/final_mitogenome.check.txt
- skip

Anopheles_nili/working/idAnoNiliSN-F31_02.hifiasm.20210327/mito-purging_e/final_mitogenome.check.txt
- skip 

Anopheles_nili/working/idAnoNiliSN-F5_01.hifiasm.20210327/mito-purging/final_mitogenome.check.txt
- OK

Anopheles_ziemanni/working/idAnoZiemDA-56_x.hifiasm.20210327/mito-purging/final_mitogenome.check.txt
- two contigs in assembly hap_ptg000505l_1_1  ptg000383l_1, both with frameshifts, reads annotation has truncated ND4 (\~500 bp deletion due to frameshift), hicanu mt assembly looks good - replace with hicanu


### Polishing and scaffolding

Polish only samples with 10x, others only scaffold

### Nucmer dotplots

Scaffolds vs public ref


### RNAseq QC

Done vs purged contigs where available. Otherwise vs existing ref.

