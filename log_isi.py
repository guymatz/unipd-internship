"""
File: log_isi.py
Author: Guy Matz
Email: gmatz@matz.org
Github: https://www.github.com/guymatz/
Description: LogISI Burst Detection Algorithm?  This Python implementation
    was translated from:
    https://github.com/ellesec/burstanalysis/blob/master/Burst_detection_methods/logisi_pasq_method.R
    The original may or may not be correct!
"""
import math
import numpy as np

def get_peaks(histogram:list[float], pd:int =2, threshold:float =0, np=NULL):
    
    """
    Finds peaks in logISI histogram

    Args:
        histogram: 
        pd: 
        threshold: 
        np:

    Returns:
        return description

    Raises:
    """

    m:int = 0
    L:int = length(h$density)
    j:int = 0
    Np:int = ifelse(is.null(Np), L, Np)
    peaks:list[float] = []
    locs:list[float] = []
    void_threshold:float = 0.7
    while((j<L) and (m<Np)):
        j = j+1
        endL = max(1,j-Pd)
        if (m>0  and  j<min(c(locs[m]+Pd, L-1))):
            j = min(c(locs[m]+Pd, L-1))
            endL = j-Pd
        endR = min(L, j+Pd)
        temp = h$density[endL:endR]
        aa = which(j==endL:endR)
        temp[aa] = -Inf
        if (Pd>1):
            idx1 = max(1, aa-2)
            idx2 = min(aa+2, length(temp))
            idx3 = max(1, aa-1)
            idx4 = min(aa+1, length(temp))
            if (sum((h$density[j]>(temp[c(1:idx1, idx2:length(temp))]+Th))==FALSE)==0  and  sum((h$density[j]>(temp[idx3:idx4]))==FALSE)==0  and  j!=1  and  j!=L):
                m = m+1
                peaks[m] = h$density[j]
                locs[m] = j
        elif (sum((h$density[j]>(temp+Th))==FALSE)==0 ) :
            m = m+1
            peaks[m] = h$density[j]
            locs[m] = j

    return data.frame(peaks=peaks, locs=locs)

def _find_bursts(spike_times: list[float], min_ibi: int,
                 min_duration: int, min_spikes: int, isi_low: int) -> list[list[float]]:
    nspikes: int = len(spike_times)

    ## Create a temp array for the storage of the bursts.  Assume that
    ## it will not be longer than Nspikes/2 since we need at least two
    ## spikes to be in a burst.

    max_bursts: int = math.floor(nspikes / 2)
    bursts: list[list[float]] = [None] * max_bursts  # type: ignore[list-item]
    burst_num: int = 0 #current burst number

    ## Phase 1 -- burst detection. Each interspike interval of the data
    ## is compared with the threshold THRESHOLD. If the interval is greater
    ## than the threshold value, it can not be part of a burst; if the
    ## interval is smaller or equal to the threshold, the interval may be
    ## part of a burst.

    ## last_end is the time of the last spike in the previous burst.
    ## This is used to calculate the IBI.
    ## For the first burst, this is no previous IBI
    last_end: float = None #for first burst, there is no IBI.

    eps: float = 10**(-10)
    in_burst: bool = False

    next_isi: float = 0.0
    end: float = 0.0
    end: float = 0.0
    res: list[float] = []
    ibi: float = 0
    for idx, spike in  enumerate(spikes_times[1:]):
        next_isi = spike_times[idx] - spike_times[idx - 1]
        if in_burst:
            if (next_isi - isi_low > eps):
                ## end of burst
                end = idx - 1
                in_burst = False

                ibi = spike_times[beg] - last_end
                last_end = spike_times[end]
                res = [beg, end, ibi]
                burst = burst + 1
                if burst > max_bursts:
                    raise ValueError("too many bursts!!!")
                bursts.append(res)
        else:
            ## not yet in burst.
            if (next_isi - isi_low <= eps):
                ## Found the start of a new burst.
                beg = n - 1
                in_burst = TRUE

    ## At the end of the burst, check if we were in a burst when the
    ## train finished.
    if in_burst:
        end = nspikes
        ibi = spike_times[beg] - last_end
        res = [beg, end, ibi]
        burst_num = burst_num + 1
        if (burst > max_bursts):
            raise ValueError("too many bursts!!!")
        bursts.append(res)

    ## Check if any bursts were found.
    if burst_num > 0:
        ## truncate to right length, as bursts will typically be very long.
        bursts = bursts[0:burst_num]
    else:
        ## no bursts were found, so return an empty structure.
        return []

    if debug:
        print("End of phase1\n")
        print(bursts)

    ## Phase 2 -- merging of bursts.    Here we see if any pair of bursts
    ## have an IBI less than MIN.IBI; if so, we then merge the bursts.
    ## We specifically need to check when say three bursts are merged
    ## into one.

    #ibis = bursts[, "IBI"]
    ibis:list[float] = [ bursts ]
    merged_bursts:list[float] = [ ibi for ibi in ibis if ibi < min_ibi ]

    if (any(merge.bursts)):
        ## Merge bursts efficiently.    Work backwards through the list, and
        ## then delete the merged lines afterwards.    This works when we
        ## have say 3+ consecutive bursts that merge into one.

        for (burst in rev(merge.bursts)):
            bursts[burst - 1, "end"] = bursts[burst, "end"]
            bursts[burst, "end"] = NA #not needed, but helpful.
        bursts = bursts[-merge.bursts, , drop = FALSE] #delete the unwanted info.

    if (debug):
        print("End of phase 2\n")
        print(bursts)

    ## Phase 3 -- remove small bursts: less than min duration (MIN.DURN), or
    ## having too few spike_times (less than MIN.spike_times).
    ## In this phase we have the possibility of deleting all spike_times.

    ## LEN = number of spike_times in a burst.
    ## DURN = duration of burst.
    len = bursts[, "end"] - bursts[, "beg"] + 1
    durn = spike_times[bursts[, "end"]] - spike_times[bursts[, "beg"]]
    bursts = cbind(bursts, len, durn)

    rejects = which((durn < min.durn) | (len < min.spike_times))

    if (any(rejects)):
        bursts = bursts[-rejects, , drop = FALSE]

    if (nrow(bursts) == 0):
        ## All the bursts were removed during phase 3.
        bursts = no.bursts
    else:
        ## Compute mean ISIS
        len = bursts[, "end"] - bursts[, "beg"] + 1
        durn = spike_times[bursts[, "end"]] - spike_times[bursts[, "beg"]]
        mean.isis = durn / (len - 1)

        ## Recompute IBI (only needed if phase 3 deleted some cells).
        if (nrow(bursts) > 1):
            ibi2 = c(NA, calc.ibi(spike_times, bursts))
        else:
            ibi2 = NA
        bursts[, "IBI"] = ibi2

        SI = rep(1, length(mean.isis))
        bursts = cbind(bursts, mean.isis, SI)

    ## End -- return burst structure.
    bursts
