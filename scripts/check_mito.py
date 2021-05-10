# Usage example:
# python check_mito.py final_mitogenome.gb MT917182.1.gb

from Bio import SeqIO
import sys
import os

def check_mito(query_fn, ref_fn):


    q = SeqIO.read(query_fn, "genbank")
    r = SeqIO.read(ref_fn, "genbank")

    sn = os.path.abspath(query_fn).split('/')[10].split('.')[0]
    rsp = r.features[0].qualifiers['organism'][0]
    print('INFO: query sample {}, ref species {}'.format(sn, rsp))
    print('INFO: query length {}, ref length {}'.format(len(q.seq),len(r.seq)))

    # remove excessive features
    q.features = [feature for feature in q.features if feature.type != 'gene']
    r.features = [feature for feature in r.features if feature.type != 'gene']
    r.features = [feature for feature in r.features if feature.type != 'source']

    # features to match mitofinder names
    frename = {
        'COI':'COX1',
        'COII':'COX2',
        'COIII':'COX3',
        '12S ribosomal RNA':'rrnS',
        '16S ribosomal RNA':'rrnL',
        'COB':'CYTB',
        'NAD1':'ND1',
        'NAD2':'ND2',
        'NAD3':'ND3',
        'NAD4':'ND4',
        'NAD4L':'ND4L',
        'NAD5':'ND5',
        'NAD6':'ND6'
    }

    for s in q,r:
        for feature in s.features:
            if feature.type == 'CDS':
                fname = feature.qualifiers['gene'][0].upper()
                if fname in frename.keys():
                    fname = frename[fname]
                feature.qualifiers['id'] = fname
            elif feature.type == 'misc_feature':
                feature.qualifiers['id'] = feature.qualifiers['note'][0]
            else:
                fname = feature.qualifiers['product'][0]
                if fname in frename.keys():
                    fname = frename[fname]            
                # duplicate trna
                if fname.startswith('tRNA') and fname.endswith('2'):
                    fname = fname[:-1]
                feature.qualifiers['id'] = fname

    q_feature_ids = set([feature.qualifiers['id'] for feature in q.features])
    r_feature_ids = set([feature.qualifiers['id'] for feature in r.features])


    print('INFO: {} features in query ({} unique), {} features in ref ({} unique)'.format(
                len(q.features), len(q_feature_ids), len(r.features), len(r_feature_ids)))
    q_extra_features = q_feature_ids-r_feature_ids
    if len(q_extra_features) > 0:
        print('ERROR: Extra features in query: {}'.format(q_extra_features))
    r_extra_features = r_feature_ids-q_feature_ids
    if len(r_extra_features) > 0:
        print('ERROR: Missing features in query: {}'.format(r_extra_features))

    for q_feature in q.features:
        qid = q_feature.qualifiers['id']
        if q_feature.type == 'CDS':
            qprot = q_feature.qualifiers['translation'][0]
            # stop acceptable at last position
            if '*' in qprot[:-1]:
                print('ERROR: Found internal stop codon in {}'.format(qid))

        ref_features = [r_feature for r_feature in r.features 
                            if r_feature.qualifiers['id'] == qid]
        if len(ref_features) == 0:
            print('WARN: Could not find reference feature for {}, skipping length check'.format(qid))
            continue
        if len(ref_features) > 1:
            print('INFO: Multiple ref features for {}'.format(qid))
        r_feature = ref_features[0]

        qlen = q_feature.location.end - q_feature.location.start
        rlen = r_feature.location.end - r_feature.location.start

        if qlen - rlen != 0:
            if abs(qlen - rlen) == 1:
                diff = 'INFO: 1bp'
            elif abs(qlen - rlen) > 9:
                diff = 'ERROR: >9bp'
            else:
                diff = 'WARN: 2-9bp'
            print('{} length mismatch for {} between query '
                  '({}, {}) and ref ({},{})'.format(diff,qid,qlen,q_feature.location,
                            rlen,r_feature.location))

if __name__ == "__main__":
    check_mito(sys.argv[1], sys.argv[2])