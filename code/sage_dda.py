import pandas
from sys import argv
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import seaborn
import numpy
import os
import re
from collections import defaultdict

plt.rcParams.update({'axes.labelsize': 12,
                     'xtick.labelsize': 11,
                     'ytick.labelsize': 11,
                     'legend.fontsize': 11})

BEST_PSM = True
CTERM_FILTER = None
Q_CUTOFF = 0.01

def parse_sage(csv_path, tag_mass_min=80):
    psms = pandas.read_csv(csv_path, delimiter='\t')
    bare_sequences = []
    varmod_masses = []
    tags = []
    num_tags = []
    matched_peak_pcts = []
    modmass_re = re.compile(r'\[(.*?)\]')
    scans = []
    for i, row in psms.iterrows():
        scans.append(int(row['scannr'].split('scan=')[-1]))
        bare_sequence = ''.join(c for c in row['peptide'] if c.isupper() and c.isalpha())
        mod_masses = [float(s) for s in modmass_re.findall(row['peptide'])]
        tag_masses = set([m for m in mod_masses if m>tag_mass_min])
        num_tags.append(len([m for m in mod_masses if m>tag_mass_min]))
        if len(tag_masses) == 0:
            tags.append(0)
        elif len(tag_masses) == 1:
            [tag_mass] = tag_masses
            tags.append(tag_mass)
        else:
            tags.append(-1)
        bare_sequences.append(bare_sequence)
        varmod_masses.append(mod_masses)
        possible_ions = ((row['peptide_len']-1) * 2 * (row['charge']-1)) or 1
        matched_peak_pcts.append((row['matched_peaks']/possible_ions) * 100)
    psms = psms.assign(scan=scans, bare_sequence=bare_sequences, varmod_masses=varmod_masses, tag=tags, num_tags=num_tags, matched_peak_pct=matched_peak_pcts)
    return psms

def ion_key(ion):
    ion, z = ion.split('+')
    by = ion[0]
    ind = int(ion[1:])
    return (by, ind, z)


#read Sage output results into concatenated DF
psm_dfs = []
ion_dfs = []
for d in argv[1:]:
    psms = parse_sage(os.path.join(d, 'results.sage.tsv'))
    lfq = pandas.read_csv(os.path.join(d, 'lfq.tsv'), sep='\t')
    ions = pandas.read_csv(os.path.join(d, 'matched_fragments.sage.tsv'), sep='\t')
    lfq['intensity'] = lfq.iloc[:,-1:]
    psms = psms.merge(lfq, how='left', on=['peptide', 'charge'])
    [filename] = psms['filename'].unique()
    ions['filename'] = filename
    psm_dfs.append(psms)
    ion_dfs.append(ions)
    
psms = pandas.concat(psm_dfs)
psms = psms[psms['spectrum_q'] < Q_CUTOFF]
if CTERM_FILTER is not None:
    psms = psms[psms['bare_sequence'].str.endswith(CTERM_FILTER)]
ions = pandas.concat(ion_dfs)


#iterate through PSMs to log matched fragments and fragmentation sites
tag2frag_total = defaultdict(lambda: defaultdict(int))
tag2frag_matched = defaultdict(lambda: defaultdict(int))
tag2frag_int = defaultdict(lambda: defaultdict(list))
unique_frags = set()
psm_id2sites = defaultdict(set)
psm_id2frags = defaultdict(set)
for i, row in psms.merge(ions, on=['filename', 'psm_id']).iterrows():
    if row['fragment_type'] == 'y':
        site = len(row['bare_sequence']) - row['fragment_ordinals']
    elif row['fragment_type'] == 'b':
        site = row['fragment_ordinals']
    psm_id2sites[(row['filename'], row['psm_id'])].add(site)
    psm_id2frags[(row['filename'], row['psm_id'])].add((row['fragment_type'], row['fragment_ordinals']))
    ion = f"{row['fragment_type']}{row['fragment_ordinals']}+{row['fragment_charge']}"
    tag2frag_matched[row['tag']][ion] += 1
    tic = row['ms2_intensity'] / (row['matched_intensity_pct']/100)
    tag2frag_int[row['tag']][ion].append(row['fragment_intensity']/tic)
    unique_frags.add(ion)

missing_sites = []
matched_pct = []
for i, row in psms.iterrows():
    missing_sites.append(len(row['bare_sequence'])-1 - len(psm_id2sites[(row['filename'], row['psm_id'])]))
    matched_pct.append(len(psm_id2frags[(row['filename'], row['psm_id'])])/((len(row['bare_sequence'])-1)*2)*100)
    for z in range(1, row['charge'] + 1):
        for series in ('b','y'):
            for i in range(1, row['peptide_len']):
                tag2frag_total[row['tag']][f'{series}{i}+{z}'] += 1
psms = psms.assign(missing_sites=missing_sites,
                   matched_pct=matched_pct)

tagged_psms = psms[psms['tag'] > 0]
lf_psms = psms[psms['tag'] == 0]

#take the best PSM per peptide if this flag is enabled
if BEST_PSM:
    tagged_psms = tagged_psms.sort_values('hyperscore', ascending=False).drop_duplicates(['peptide', 'charge'])
    lf_psms = lf_psms.sort_values('hyperscore', ascending=False).drop_duplicates(['peptide', 'charge'])

#Venn diagram of unique backbone sequences    
fig = plt.figure()
plt.title('Unique peptide sequences')
venn2([set(lf_psms['bare_sequence'].unique()), set(tagged_psms['bare_sequence'].unique())],
      set_labels=['Label free', 'Tag6 labeled'],
      set_colors=['steelblue', 'darkorange'])
fig.set_size_inches(5,5)

#perform intersection on (backbone sequence, charge)
intersect = lf_psms.merge(tagged_psms, on=['bare_sequence','charge'])
intersect = intersect[intersect['filename_x'] != intersect['filename_y']]
print('score delta=', numpy.median(intersect['hyperscore_y'] - intersect['hyperscore_x']))
int_ratios = numpy.log10(intersect['intensity_y']/intersect['intensity_x'])
score_ratios = numpy.log2(intersect['hyperscore_y']/intersect['hyperscore_x'])

#for hexbin: compute number of datapoints in each quadrant
q1 = sum(intensity<0 and score>0 for intensity, score in zip(int_ratios, score_ratios))
q2 = sum(intensity>0 and score>0 for intensity, score in zip(int_ratios, score_ratios))
q3 = sum(intensity<0 and score<0 for intensity, score in zip(int_ratios, score_ratios))
q4 = sum(intensity>0 and score<0 for intensity, score in zip(int_ratios, score_ratios))
n_all = len(int_ratios)

#hexbin of intensity vs. score ratios
fig = plt.figure()
plt.hexbin(int_ratios, score_ratios, bins=50, mincnt=1, cmap='autumn')
plt.title('Tag6 / unlabeled Human digest: intersected peptides')
plt.xlabel('Ratio of precursor intensity, log₁₀')
plt.ylabel('Ratio of score, log₂')
plt.xlim(-4, 4)
plt.ylim(-2.5, 2.5)
plt.axvline(0, color='gray', linestyle='--')
plt.axhline(0, color='gray', linestyle='--')
plt.gca().set_aspect(4/2.5)
fig.set_size_inches(5,5)
plt.text(-3.75, 2.25, f'n={q1}\n({(q1/n_all*100):.2f}%)', ha='left', va='top', fontsize=11, weight='bold')
plt.text(3.75, 2.25, f'n={q2}\n({(q2/n_all*100):.2f}%)', ha='right', va='top', fontsize=11, weight='bold')
plt.text(-3.75, -2.25, f'n={q3}\n({(q3/n_all*100):.2f}%)', ha='left', va='bottom', fontsize=11, weight='bold')
plt.text(3.75, -2.25, f'n={q4}\n({(q4/n_all*100):.2f}%)', ha='right', va='bottom', fontsize=11, weight='bold')

#fragment frequency plots for +1 and +2 fragment ions
for z in (1, 2):
    fig, axs = plt.subplots(2)
    fig.set_size_inches(10,5)
    plt.subplots_adjust(hspace=0.4)
    for i, series in enumerate(['y','b']):
        [tag] = [tag for tag in tag2frag_total.keys() if tag > 0] 
        fig.suptitle(f'Frequency of fragment production (all peptides)')
        frag_names = []
        frag_ratios = []
        frag_ratios_lf = []
        for frag in sorted(unique_frags, key=ion_key):
            if not frag.startswith(series):
                continue
            if frag.endswith(f'+{z}'):
                frag_names.append(frag.split('+')[0][1:])
                frag_ratios.append(tag2frag_matched[tag][frag] / (tag2frag_total[tag][frag] or 1)*100)
                frag_ratios_lf.append(tag2frag_matched[0][frag] / (tag2frag_total[0][frag] or 1)*100)
        ax = axs[i-1]
        ax.set_ylim([0, 100/z])
        fig.text(0.04, 0.5, '% of theoretical fragments observed', va='center', rotation='vertical', fontsize=12)
        ax.set_title(f'{series} ion series (+{z} charge state)')
        ax.bar(frag_names, frag_ratios_lf, label='Label free', width=-0.4, align='edge')
        ax.bar(frag_names, frag_ratios, label='Tag6 labeled', width=0.4, align='edge')
        ax.legend()

#fragment intensity plots for +1 and +2 fragment ions 
for z in (1, 2):
    fig, axs = plt.subplots(2)
    fig.set_size_inches(10,5)
    plt.subplots_adjust(hspace=0.4)
    for i, series in enumerate(['y','b']):
        [tag] = [tag for tag in tag2frag_total.keys() if tag > 0] 
        fig.suptitle(f'Fragment relative intensity (all peptides)')
        frag_names = []
        frag_ratios = []
        frag_ratios_lf = []
        for frag in sorted(unique_frags, key=ion_key):
            if not frag.startswith(series):
                continue
            if frag.endswith(f'+{z}'):
                frag_names.append(frag.split('+')[0][1:])
                frag_ratios.append(numpy.mean(tag2frag_int[tag][frag])*100)
                frag_ratios_lf.append(numpy.mean(tag2frag_int[0][frag])*100)
        ax = axs[i-1]
        ax.set_ylim([0, 6.5/z])
        fig.text(0.04, 0.5, 'Mean % of TIC', va='center', rotation='vertical', fontsize=12)
        ax.set_title(f'{series} ion series (+{z} charge state)')
        ax.bar(frag_names, frag_ratios_lf, label='Label free', width=-0.4, align='edge')
        ax.bar(frag_names, frag_ratios, label='Tag6 labeled', width=0.4, align='edge')
        ax.legend()

#fragmentation site coverage (intersected)
fig = plt.figure()
plt.title('Fragmentation site coverage (intersected peptides)')
bins = numpy.arange(0,17,1)
tagged_weights = [100/len(tagged_psms)]*len(tagged_psms)
lf_weights = [100/len(lf_psms)]*len(lf_psms)
plt.hist([intersect['missing_sites_x'], intersect['missing_sites_y']], label=['Label free', 'Tag6 labeled'], weights=[[100/len(intersect)]*len(intersect), [100/len(intersect)]*len(intersect)], bins=bins-0.5)
fig.set_size_inches(5,5)
plt.xticks(bins[:-1])
plt.xlabel('Number of missing fragmentation sites')
plt.ylabel('% of total peptides')
plt.legend()

#complete fragmentation site coverage vs. length (non-intersected)
fig = plt.figure()
plt.title('Complete fragmentation site coverage')
pep_lens = list(range(6, 26))
tagged_complete = [sum(tagged_psms[tagged_psms['peptide_len']==l]['missing_sites']==0)/(sum(tagged_psms['peptide_len']==l) or 1)*100 for l in pep_lens]
lf_complete = [sum(lf_psms[lf_psms['peptide_len']==l]['missing_sites']==0)/(sum(lf_psms['peptide_len']==l) or 1)*100 for l in pep_lens]
plt.plot(pep_lens, lf_complete, label='Label free')
plt.plot(pep_lens, tagged_complete, label='Tag6 labeled')
fig.set_size_inches(5,5)
plt.xticks(pep_lens, [str(i) if i%2==0 else '' for i in pep_lens])
plt.xlabel('Peptide length')
plt.ylabel('% complete fragmentation site coverage')
plt.legend()

#gain in matched peaks histogram
fig = plt.figure()
plt.title('Tag6 gain in matched peaks (intersected peptides)')
bins = numpy.arange(-25,25,1)
plt.hist(intersect['matched_peaks_y'] - intersect['matched_peaks_x'], weights=[100/len(intersect)]*len(intersect), bins=bins-0.5)
fig.set_size_inches(5,5)
mean = numpy.mean(intersect['matched_peaks_y'] - intersect['matched_peaks_x'])
plt.axvline(mean, color='red', label=f'mean = {mean:.2f}')
plt.xlabel('Matched peak delta vs. label free')
plt.ylabel('% of total peptides')
plt.legend()

#density plot of % fragments matched
fig = plt.figure()
plt.title('Fragments matched (intersected peptides)')
bins = numpy.arange(0,107.5,5)
tagged_weights = [100/len(tagged_psms)]*len(tagged_psms)
lf_weights = [100/len(lf_psms)]*len(lf_psms)
seaborn.kdeplot(intersect['matched_pct_x'], fill=True, label='Label free', bw=0.2, clip=(0,100))
seaborn.kdeplot(intersect['matched_pct_y'], fill=True, label='Tag6 labeled', bw=0.2, clip=(0,100))
fig.set_size_inches(5,5)
plt.xticks(range(0,110,10))
plt.xlabel('% of theoretical fragments matched')
plt.ylabel('Density')
plt.legend()

#density plot of % MS2 intensity matched
fig = plt.figure()
plt.title('MS2 intensity matched (intersected peptides)')
bins = numpy.arange(0,107.5,5)
tagged_weights = [100/len(tagged_psms)]*len(tagged_psms)
lf_weights = [100/len(lf_psms)]*len(lf_psms)
seaborn.kdeplot(intersect['matched_intensity_pct_x'], fill=True, label='Label free', bw=0.2, clip=(0,100))
seaborn.kdeplot(intersect['matched_intensity_pct_y'], fill=True, label='Tag6 labeled', bw=0.2, clip=(0,100))
fig.set_size_inches(5,5)
plt.xticks(range(0,110,10))
plt.xlabel('% of MS2 intensity matched')
plt.ylabel('Density')
plt.legend()
    
plt.show()
