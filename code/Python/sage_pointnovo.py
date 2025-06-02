from sys import argv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from utils import parse_sage, aa2mass
from nltk import edit_distance
from collections import defaultdict
import os

plt.rcParams.update({'axes.titlesize': 14,
                     'axes.labelsize': 14,
                     'xtick.labelsize': 13,
                     'ytick.labelsize': 13,
                     'legend.fontsize': 13})

BEST_PSM = True
DIST_CUTOFF = 0
Q_CUTOFF = 0.01
TAG_NAME = 'PSMtag'

psm_dfs = []
ion_dfs = []
for d in argv[1:]:
    psms = parse_sage(os.path.join(d, 'results.sage.tsv'))
    lfq = pd.read_csv(os.path.join(d, 'lfq.tsv'), sep='\t')
    ions = pd.read_csv(os.path.join(d, 'matched_fragments.sage.tsv'), sep='\t')
    lfq['intensity'] = lfq.iloc[:,-1:]
    psms = psms.merge(lfq, how='left', on=['peptide', 'charge'])
    [filename] = psms['filename'].unique()
    ions['filename'] = filename
    psms = psms[psms['spectrum_q'] < Q_CUTOFF]
    pointnovo = pd.read_csv(os.path.join(d, 'pointnovo.csv'), sep='\t')
    pointnovo = pointnovo[~pointnovo['predicted_position_score'].isna()]
    
    scores = []
    scans = []
    bare_sequences = []
    for i, row in pointnovo.iterrows():
        aa_scores = [float(s) for s in row['predicted_position_score'].split(',')]
        scores.append(row['predicted_score'])
        scans.append(int(row['feature_id'].split(':')[-1]))
        bare_sequences.append(''.join(aa[0] for aa in row['predicted_sequence'].split(',')))
    pointnovo = pointnovo.assign(scan=scans, score=scores, bare_sequence=bare_sequences)

    merged = pointnovo.merge(psms, on='scan')
    psm_dfs.append(merged)
    ion_dfs.append(ions)

merged = pd.concat(psm_dfs)
ions = pd.concat(ion_dfs)

distances = []
aa_matched = []
for i, row in merged.iterrows():
    distance = edit_distance(row['bare_sequence_x'].replace('I', 'L'), row['bare_sequence_y'].replace('I', 'L'), transpositions=True)
    distances.append(distance)
    
    sage_aas = []
    novo_aas = []
    for i, aa in enumerate(row['bare_sequence_y']):
        sage_aas.append({'aa':aa2mass[aa], 'prefix':sum(aa2mass[aap] for aap in row['bare_sequence_y'][:i])})
    for i, aa in enumerate(row['bare_sequence_x']):
        novo_aas.append({'aa':aa2mass[aa], 'prefix':sum(aa2mass[aap] for aap in row['bare_sequence_x'][:i])})
        
    aa_match = [False] * len(row['bare_sequence_x'])
    for i, aa in enumerate(novo_aas):
        if any(abs(aa['aa']-aa2['aa'])<0.02 and abs(aa['prefix']-aa2['prefix'])<0.1 for aa2 in sage_aas):
            aa_match[i] = True
    aa_matched.append(aa_match)
    
merged = merged.assign(distance=distances, aa_matched=aa_matched)

tagged_psms = merged[merged['tag'] > 0]
lf_psms = merged[merged['tag'] == 0]
if BEST_PSM:
    tagged_psms = tagged_psms.sort_values('score_x', ascending=False).drop_duplicates(['peptide', 'charge'])
    lf_psms = lf_psms.sort_values('score_x', ascending=False).drop_duplicates(['peptide', 'charge'])

intersect = lf_psms.merge(tagged_psms, on=['untagged_sequence','charge'])
intersect = intersect[intersect['filename_x'] != intersect['filename_y']]

aa2matched_lf = defaultdict(list)
aa2matched_tagged = defaultdict(list)
pep_len_filter = 0
for _, row in merged[merged['peptide_len']>=pep_len_filter].iterrows():
    for i in [0,1,2,3,4]:
        if i >= len(row['aa_matched']):
            continue
        if row['tag'] == 0:
            aa2matched_lf[f'{i}'].append(row['aa_matched'][i])
            aa2matched_lf[f'-{i}'].append(row['aa_matched'][-i-1])
        else:
            aa2matched_tagged[f'{i}'].append(row['aa_matched'][i])
            aa2matched_tagged[f'-{i}'].append(row['aa_matched'][-i-1])

fig = plt.figure()
aa_is = ['0','1','2','3','4','...','-4','-3','-2','-1','-0']
hatches = [None]*6 + ['//']*5
probs_lf = np.asarray([np.mean(aa2matched_lf[i])*100 for i in aa_is])
probs_tagged = np.asarray([np.mean(aa2matched_tagged[i])*100 for i in aa_is])
plt.bar(aa_is, probs_lf, width=-0.4, align='edge', hatch=hatches)
plt.bar(aa_is, probs_tagged, width=0.4, align='edge', hatch=hatches)
plt.bar(['...'], [0], color='gray', label='N-terminal AA')
plt.bar(['...'], [0], color='gray', hatch='//', label='C-terminal AA')
plt.xlabel('AA distance from terminus')
plt.ylabel('% AAs correctly matched at site')
plt.ylim(0,100)
plt.legend()
fig.set_size_inches(5,5)

fig = plt.figure()
fig.set_size_inches(5,5)
pep_lens = list(range(4, 26))
tagged_complete = [sum(intersect[intersect['peptide_len_y']==l]['distance_y']==0) for l in pep_lens]
lf_complete = [sum(intersect[intersect['peptide_len_x']==l]['distance_x']==0) for l in pep_lens]
plt.plot(pep_lens, lf_complete, label='LF (intersected)')
plt.plot(pep_lens, tagged_complete, label=f'{TAG_NAME} (intersected)')
plt.xticks(pep_lens, [str(i) if i%2==0 else '' for i in pep_lens])
plt.ylabel('Sequences matching database search')
tagged_complete = [sum(tagged_psms[tagged_psms['peptide_len']==l]['distance']==0) for l in pep_lens]
lf_complete = [sum(lf_psms[lf_psms['peptide_len']==l]['distance']==0) for l in pep_lens]
plt.plot(pep_lens, lf_complete, label='LF (all)', linestyle='--', color='tab:blue')
plt.plot(pep_lens, tagged_complete, label=f'{TAG_NAME} (all)', linestyle='--', color='tab:orange')
plt.xlabel('Peptide length')
plt.legend()

precision_lf = []
recall_lf = []
precision_t6 = []
recall_t6 = []
precision_lf_all = []
recall_lf_all = []
precision_t6_all = []
recall_t6_all = []
fdr_lf = []
fdr_t6 = []

ap_lf = len(intersect)
ap_t6 = len(intersect)

ap_lf_all = len(lf_psms)
ap_t6_all = len(tagged_psms)

tps_lf = []
fps_lf = []
tps_t6 = []
fps_t6 = []
tps_lf_all = []
fps_lf_all = []
tps_t6_all = []
fps_t6_all = []

for score in np.arange(0, -2, -0.01):
    filtered_lf = lf_psms[lf_psms['score_x'] >= score]
    filtered_t6 = tagged_psms[tagged_psms['score_x'] >= score]
    tp_lf_all = sum(filtered_lf['distance'] <= DIST_CUTOFF)
    fp_lf_all = sum(filtered_lf['distance'] > DIST_CUTOFF)
    tp_t6_all = sum(filtered_t6['distance'] <= DIST_CUTOFF)
    fp_t6_all = sum(filtered_t6['distance'] > DIST_CUTOFF)

    filtered_lf = intersect[intersect['score_x_x'] >= score]
    filtered_t6 = intersect[intersect['score_x_y'] >= score]
    tp_lf = sum(filtered_lf['distance_x'] <= DIST_CUTOFF)
    fp_lf = sum(filtered_lf['distance_x'] > DIST_CUTOFF)
    tp_t6 = sum(filtered_t6['distance_y'] <= DIST_CUTOFF)
    fp_t6 = sum(filtered_t6['distance_y'] > DIST_CUTOFF)

    tps_lf.append(tp_lf)
    fps_lf.append(fp_lf)
    tps_t6.append(tp_t6)
    fps_t6.append(fp_t6)
    tps_lf_all.append(tp_lf_all)
    fps_lf_all.append(fp_lf_all)
    tps_t6_all.append(tp_t6_all)
    fps_t6_all.append(fp_t6_all)

    try: precision_lf.append(tp_lf / (tp_lf + fp_lf))
    except: precision_lf.append(1)
    try: precision_t6.append(tp_t6 / (tp_t6 + fp_t6))
    except: precision_t6.append(1)
    recall_lf.append(tp_lf / ap_lf)
    recall_t6.append(tp_t6 / ap_t6)

    try: precision_lf_all.append(tp_lf_all / (tp_lf_all + fp_lf_all))
    except: precision_lf_all.append(1)
    try: precision_t6_all.append(tp_t6_all / (tp_t6_all + fp_t6_all))
    except: precision_t6_all.append(1)
    recall_lf_all.append(tp_lf_all / ap_lf_all)
    recall_t6_all.append(tp_t6_all / ap_t6_all)

fig = plt.figure(constrained_layout=True)
fig.set_size_inches(5,5)
axs = fig.subplot_mosaic([['a','a'],['b','c']], gridspec_kw={'height_ratios':[2,1], 'hspace':0.05})
axs['a'].set_xlabel('Peptide recall')
axs['a'].set_ylabel('Peptide precision')
axs['a'].plot(recall_lf, precision_lf, label='LF (intersected)')
axs['a'].plot(recall_t6, precision_t6, label=f'{TAG_NAME} (intersected)')
axs['a'].plot(recall_lf_all, precision_lf_all, label='LF (all)', linestyle='--', color='tab:blue')
axs['a'].plot(recall_t6_all, precision_t6_all, label=f'{TAG_NAME} (all)', linestyle='--', color='tab:orange')
axs['a'].set_xlim(0,0.8)
axs['a'].set_ylim(0.2,1.1)

axs['b'].set_title('Total recall')
axs['b'].set_ylabel('# of peptides')
axs['b'].bar([0], [max(tps_lf)])
axs['b'].bar([1], [max(tps_t6)])
axs['b'].bar([2], [max(tps_lf_all)], color='white', edgecolor='tab:blue', linestyle='--')
axs['b'].bar([3], [max(tps_t6_all)], color='white', edgecolor='tab:orange', linestyle='--')
axs['b'].get_xaxis().set_visible(False)

axs['c'].set_title('Area under curve')
axs['c'].bar([0], [np.trapezoid(precision_lf, recall_lf)])
axs['c'].bar([1], [np.trapezoid(precision_t6, recall_t6)])
axs['c'].bar([2], [np.trapezoid(precision_lf_all, recall_lf_all)], color='white', edgecolor='tab:blue', linestyle='--')
axs['c'].bar([3], [np.trapezoid(precision_t6_all, recall_t6_all)], color='white', edgecolor='tab:orange', linestyle='--')
axs['c'].get_xaxis().set_visible(False)

plt.show()

