import pandas
import re

aa2mass = {
    'A':71.037113805,
    'C':103.009184505 + 57.021464,
    'D':115.026943065,
    'E':129.042593135,
    'F':147.068413945,
    'G':57.021463735,
    'H':137.058911875,
    'I':113.084064015,
    'K':128.094963050,
    'L':113.084064015,
    'M':131.040484645,
    'N':114.042927470,
    'O':237.147726925,
    'P':97.052763875,
    'Q':128.058577540,
    'R':156.101111050,
    'S':87.032028435,
    'T':101.047678505,
    'U':150.953633405,
    'V':99.068413945,
    'W':186.079312980,
    'Y':163.063328575,
}

def parse_sage(csv_path, tag_mass_min=80):
    psms = pandas.read_csv(csv_path, delimiter='\t')
    bare_sequences = []
    untagged_sequences = []
    varmod_masses = []
    tags = []
    num_tags = []
    matched_peak_pcts = []
    modmass_re = re.compile(r'\[(.*?)\]')
    scans = []
    for i, row in psms.iterrows():
        scans.append(int(row['scannr'].split('scan=')[-1]))
        bare_sequence = ''.join(c for c in row['peptide'] if c.isupper() and c.isalpha())
        modmass_strs = modmass_re.findall(row['peptide'])
        mod_masses = [float(s) for s in modmass_strs]
        tag_masses = set([m for m in mod_masses if m>tag_mass_min])
        untagged_sequence = row['peptide']
        for s in modmass_strs:
            if float(s) in tag_masses:
                untagged_sequence = untagged_sequence.replace(f'[{s}]','')
        if untagged_sequence.startswith('-'):
            untagged_sequence = untagged_sequence[1:]
        num_tags.append(len([m for m in mod_masses if m>tag_mass_min]))
        if len(tag_masses) == 0:
            tags.append(0)
        elif len(tag_masses) == 1:
            [tag_mass] = tag_masses
            tags.append(tag_mass)
        else:
            tags.append(-1)
        bare_sequences.append(bare_sequence)
        untagged_sequences.append(untagged_sequence)
        varmod_masses.append(mod_masses)
        possible_ions = ((row['peptide_len']-1) * 2 * (row['charge']-1)) or 1
        matched_peak_pcts.append((row['matched_peaks']/possible_ions) * 100)
    psms = psms.assign(scan=scans, bare_sequence=bare_sequences, untagged_sequence=untagged_sequences, varmod_masses=varmod_masses, tag=tags, num_tags=num_tags, matched_peak_pct=matched_peak_pcts)
    return psms

