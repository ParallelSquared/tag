import numpy as np

def xics(mzml, xic_mzs, tolerance_min=-6, tolerance_max=6, units='ppm'):
    """
    Given an iterable of pyteomics scan objects, compute an XIC for each provided M/z.

    Parameters
    ----------
    mzml : iterable
        An iterable of pyteomics scan objects (e.g. a pyteomics.mzml.MzML).
        This need only contain MS1 scans.
    xic_mzs : list
        A list of M/z values for which to produce XICs.
    tolerance_min/max : float
        Bounds for XIC extraction tolerance (min value should usually be negative).
    units : string
        tolerance_min/max is interpreted as ppm if this value is 'ppm', else Da.

    Returns a dict containing the following:
        xic_mzs : sorted array of the provided M/z values.
        xic_ints : 2D array with dimensions len(xic_mzs) * len(ms1_scans).
                   An array of XICs, each an array of extracted intensity values per scan.
        scan_inds : 1D array of len(ms1_scans).
                    MS1 scan indexes (compatible with pyteomics mzml get_by_index).
        scan_rts : 1D array of len(ms1_scans).
                   MS1 scan start times in seconds.
    
    """
    xic_mzs = sorted(xic_mzs)
    tolerance_min = float(tolerance_min)
    tolerance_max = float(tolerance_max)

    #calculate mz max/min values based on tolerance
    mz_windows = []
    for xic_mz in xic_mzs:
        if units == 'ppm':
            mz_windows.append((xic_mz+(xic_mz*tolerance_min)/1000000, xic_mz+(xic_mz*tolerance_max)/1000000))
        else:
            mz_windows.append((xic_mz+tolerance_min, xic_mz+tolerance_max))

    #extract MS1 scans from mzml
    ms1_scans = [scan for scan in mzml if scan['ms level'] == 1]

    #initialize XIC intensity arrays for each target mz
    xic_ints = np.zeros([len(xic_mzs), len(ms1_scans)])
    scan_inds = []
    scan_rts = []
    for si, scan in enumerate(ms1_scans):
        scan_inds.append(scan['index'])
        scan_rts.append(scan['scanList']['scan'][0]['scan start time'] * 60)
        mzs = scan['m/z array']
        intensities = scan['intensity array']

        #for each mz, extract and sum intensity from all ions in its mz tolerance window 
        l_ind = 0
        for (xic, (l_mz, h_mz)) in zip(xic_ints, mz_windows):
            int_sum = 0
            try:
                while mzs[l_ind] < l_mz:
                    l_ind += 1
                h_ind = l_ind
                while mzs[h_ind] <= h_mz:
                    int_sum += intensities[h_ind]
                    h_ind += 1
            except IndexError:
                pass
            xic[si] = int_sum
            
    return {'xic_mzs' : xic_mzs,
            'xic_ints': xic_ints,
            'scan_inds': np.asarray(scan_inds, dtype=int),
            'scan_rts': np.asarray(scan_rts, dtype=float)
            }

def quantify_xic(xic_ints, xic_rts, rt_min, rt_max, threshold=0, gate=2):
    """
    Simple XIC peak boundary detection and trapezoidal integration.
    
    Given an XIC (defined by an intensity array and RT array), seek to the highest peak
        within the given RT range. Then, walk outward left and right until finding [gate]
        consecutive scans at or below [threshold] intensity. Place peak boundaries and
        perform trapezoidal XIC area integration.

    Parameters
    ----------
    xic_ints : 1D array
        XIC as a 1D array of scanwise intensity values.
    xic_rts : 1D array
        Scan start times in seconds for each scanwise intensity value in xic_ints.
    rt_min/max : float
        RT boundaries in seconds for locating the XIC peak to quantify.
        Iterative peak-boundary-setting algorithm begins at the max intensity in this
            RT range. Area integration is NOT restricted to these boundaries. Reasonable
            boundaries for a peptide would be the max and min RTs at which the peptide
            was quantified. rt_min and rt_max can be the same if only one RT value is
            available.
    threshold : float
        Intensity threshold for setting peak boundaries.
    gate : int
        Number of consecutive scans required to be at or below [threshold] in order to
        place peak boundary.

    Returns a dict containing the following:
        area : total integrated area
        start_ind : array index of starting peak boundary
        end_ind : array index of end peak boundary
    """
    
    i_low, i_high = np.searchsorted(xic_rts, [rt_min, rt_max])
    if i_low == len(xic_ints):
        #handle case where RT exceeds RT range of scan
        return {'area':0, 'start_ind':i_low, 'end_ind':i_high}
    i_start = i_low + np.argmax(xic_ints[i_low:i_high+1])
    low_ind = i_start
    high_ind = i_start
    while low_ind > 0 and any(xic_ints[low_ind-i] > threshold for i in range(1,gate+1) if low_ind-i >= 0):
        low_ind -= 1
    while high_ind < len(xic_ints)-1 and any(xic_ints[high_ind+i] > threshold for i in range(1,gate+1) if high_ind+i < len(xic_ints)):
        high_ind += 1
    return {
        'area':np.trapz(xic_ints[low_ind:high_ind+1], xic_rts[low_ind:high_ind+1]),
        'start_ind':low_ind,
        'end_ind':high_ind
    }
