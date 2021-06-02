from collections import defaultdict
import pandas as pd
import sys
import os
import glob

wd = '/lustre/scratch116/tol/teams/lawniczak/data/badass/'

cutadapt_pat = wd + '*/transcriptomic_data/*/rna-seq/trimmed/cutadapt.log'
samstats_pat = wd + '*/working/*.star.*/*.markdup.stats'
busco_pat = wd + '*/working/*.trinity/busco5/short_summary.*.txt'
asm_stats_pat = wd + '*/working/*.trinity/Trinity.stats'

stats = defaultdict(dict)
# do not include full samtools stats report
samtools_stats = defaultdict(dict)

for fn in glob.glob(cutadapt_pat):
    sample = fn.split('/')[-4]
    with open(fn) as f:
        for l in f:
            if l.startswith('Total read pairs processed'):
                stats[sample]['total_read_pairs'] = l.split(':')[1].strip().replace(',','')
            elif l.startswith('Pairs written (passing filters)'):
                stats[sample]['trimmed_read_pairs'] = l.split(':')[1].strip().split(' ')[0].replace(',','')
            elif l.startswith('Total basepairs processed'):
                stats[sample]['total_bases'] = l.split(':')[1].strip().split(' ')[0].replace(',','')
            elif l.startswith('Total written (filtered)'):
                stats[sample]['trimmed_bases'] = l.split(':')[1].strip().split(' ')[0].replace(',','')
                break

for fn in glob.glob(busco_pat):
    sample = fn.split('/')[-3].split('.')[0]
    with open(fn) as f:
        for l in f:
            l = l.strip()
            if l.startswith('C:'):
                stats[sample]['busco'] = l
                break

for fn in glob.glob(asm_stats_pat):
    sample = fn.split('/')[-2].split('.')[0]
    with open(fn) as f:
        for l in f:
            l = l.strip()
            if l.startswith('sum'):
                for s in l.split(','):
                    k,v=s.split('=')
                    k = 'asm_' + k.strip()
                    stats[sample][k] = v.strip()
            elif l.startswith('N50'):
                for s in l.split(','):
                    k,v=s.split(' = ')
                    k = 'asm_N50' if k.startswith('N50') else 'asm_L50'
                    stats[sample][k] = v.strip()
                break

for fn in glob.glob(samstats_pat):
    sample = fn.split('/')[-2].split('.')[0]
    with open(fn) as f:
        for l in f:
            # parse Summary Numbers only
            if l.startswith('SN'):
                k,v=l.split(':')
                k = 'sn_'+k[3:].replace(' ','_')
                v = v.split('#')[0].strip()
                samtools_stats[sample][k] = v
            elif l.startswith('FFQ'):
                break



stats_df = pd.DataFrame.from_dict(stats).T.sort_index()
ss_df = pd.DataFrame.from_dict(samtools_stats, dtype='float').T.sort_index()

stats_df['sn_percent_mapped'] = ss_df.sn_reads_mapped / ss_df.sn_sequences * 100
stats_df['sn_percent_duplicated'] = ss_df.sn_reads_duplicated / ss_df.sn_reads_mapped * 100
stats_df['sn_percent_MQ0'] = ss_df.sn_reads_MQ0 / ss_df.sn_reads_mapped * 100
stats_df['sn_percent_mismatch'] = ss_df.sn_error_rate * 100
stats_df['sn_insert_size_average'] = ss_df.sn_insert_size_average
stats_df['sn_insert_size_standard_deviation'] = ss_df.sn_insert_size_standard_deviation
stats_df['sn_average_length'] = ss_df.sn_average_length

sys.stdout.write(stats_df.to_csv(float_format='%.2f'))

#sys.stderr.write(ss_df.to_csv())