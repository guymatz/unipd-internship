"""
File: max_interval.py
Author: Guy Matz
Email: gmatz@matz.org
Github: https://www.github.com/guymatz/
Description: MaxInterval Burst Detection Algorithm?  This Python implementation
    was translated from:
    https://github.com/igm-team/meaRtools/blob/master/meaRtools/R/maxinterval.R
    The original may or may not be correct!
"""

# import typing as t
import json
import math
import numbers
import numpy as np

from utils import _calc_ibi

# import numpy.typing as npt


###  WHY DOESN"T mypy ACCEPT THIS?!?!?
# def max_interval(spike_times: list[float] | npt.NDArray[np.float64],
# pylint: disable=too-many-positional-arguments,too-many-arguments
def max_interval(
    spike_times: list[float],
    max_begin_isi: float,
    max_end_isi: float,
    min_ibi: float,
    min_burst_duration: float,
    min_spikes_in_burst: int,
    sampling_rate: int,
) -> list[list[float]]:
    """
    Convenience function wrapping the methods below

    Args:
    ----------
    - spike_times: 1D list or NumpPy array of spike timestamps (sorted)
    - max_begin_isi: Max interval between first two spikes to start a burst
    - max_end_isi: Max interval between spikes to continue a burst
    - min_ibi: Smallest allowable Inter-Burst Interval
    - min_burst_duration: Minimum duration (ms/sec) for a cluster to be a burst
    - min_spikes_in_burst: Minimum number of spikes required in a burst
    - sampling_rate: Samples per second

    Returns:
    -------
    - bursts: Dicts of [str:List] containing:
          'start_time', 'end_time', 'duration', and 'num_spikes'
    """
    if isinstance(spike_times, np.ndarray):
        spike_times = spike_times.tolist()

    mi_bursts: list[list[float]] = _burst_detection(
        spike_times, max_begin_isi, max_end_isi, sampling_rate
    )
    mi_merged_bursts: list[list[float]] = _merge_bursts(mi_bursts, min_ibi)
    mi_burst_data: list[list[float]] = _quality_control(
        spike_times,
        sampling_rate,
        mi_merged_bursts,
        min_spikes_in_burst,
        min_burst_duration,
    )
    return mi_burst_data


# pylint: disable=too-many-locals
def _burst_detection(
    spike_times: list[float],
    beg_isi: float,
    end_isi: float,
    sampling_rate: float,
) -> list[list[float]]:
    """

    Args:
    - spike_times: 1D list or NumpPy array of spike timestamps (sorted)
    - max_begin_isi: Max interval between first two spikes to start a burst
    - max_end_isi: Max interval between spikes to continue a burst
    - min_burst_duration: Minimum duration (ms/sec) for a cluster to be a burst
    - min_spikes_in_burst: Minimum number of spikes required in a burst

    Returns:
        - bursts: Numpy array of burst times
    """

    spikes: list[float] = [t / sampling_rate for t in spike_times]
    nspikes = len(spikes)
    max_bursts = math.floor(nspikes / 2)

    bursts: list[list[float]] = [None] * max_bursts  # type: ignore[list-item]
    in_burst: bool = False
    next_isi: float = 0.0
    last_end: float = 0.0

    burst: int = 0
    beg: int = 0
    end: int = 0

    for n in range(1, len(spikes)):
        next_isi = spikes[n] - spikes[n - 1]
        if in_burst:
            # Store spike in list for this burst
            if next_isi > end_isi:
                # Burst has ended . . .  add the list of spikes in the burst
                # to all_burst_data
                # ... so we start anew
                end = n - 1
                in_burst = False
                if last_end == 0:
                    ibi = None
                else:
                    ibi = spikes[beg] - last_end
                last_end = spikes[end]
                res = [beg, end, ibi]
                if burst > max_bursts:
                    raise ValueError("Too Many Bursts")
                # print(f"Burst Ended: {n}")
                # Must ignore due to needing a None in first burst
                bursts[burst] = res  # type: ignore[assignment]
                burst = burst + 1
        else:
            if next_isi < beg_isi:
                beg = n - 1
                in_burst = True
                # print(f"Burst Started: {n}")
    if in_burst:
        end = nspikes
        ibi = spikes[beg] - last_end
        res = [beg, end, ibi]
        burst = burst + 1
        if burst > max_bursts:
            raise ValueError("Too Many Bursts")
        # Must ignore due to needing a None in first burst
        bursts[burst] = res  # type: ignore[assignment]

    if burst > 0:
        bursts = bursts[0:burst]

    return bursts


def _merge_bursts(mb_bursts: list[list[float]], min_ibi: float) -> list[list[float]]:
    """
    Phase 2 - Merging of Bursts

    Here we see if any pair of bursts have an IBI less than min_IBI; if so, we
    then merge the bursts. We specifically need to check when say three bursts
    are merged into one.

    Args:
        bursts: Numpy array of array of bursts
        ibis: Numpy array of ISIs for each bursts in bursts
        min_ibi (float): Minimum

    Returns:
        bursts: Numpy array of spike durations
    """
    mb_merged_bursts: list[int] = []
    burst: int = -1

    # Get indexes of bursts to merge
    mb_merged_bursts = [n for n, x in enumerate(mb_bursts) if x[2] and x[2] < min_ibi]

    if len(mb_merged_bursts) > 0:
        for burst in mb_merged_bursts[::-1]:
            # combine the two bursts, setting end of burst n-1 to end of burst n
            mb_bursts[burst - 1][1] = mb_bursts[burst][1]
            mb_bursts[burst] = [0, 0, 0]
        # now we delete the burst merged burst
        for burst in mb_merged_bursts[::-1]:
            b = mb_bursts.pop(burst)

    # breakpoint()
    for b in mb_bursts:
        if isinstance(b[0], numbers.Number):
            b[0] = b[0] + 1
        if isinstance(b[1], numbers.Number):
            b[1] = b[1] + 1

    return mb_bursts


def _quality_control(
    spikes: list[float],
    sampling_rate: int,
    bursts: list[list[float]],
    min_spikes_in_burst: float,
    min_burst_duration: float,
) -> list[list[float]]:
    """
    Phase 3 - Quality Control

    Remove small bursts less than min_bursts_duration or having too few
    spikes less than min_spikes_in_bursts. In this phase we have the
    possibility of deleting all spikes.
    """

    # Normalize spikes
    spikes = [s / sampling_rate for s in spikes]

    burst_lens: list[float] = [
        bursts[n][1] - bursts[n][0] + 1 for n in range(len(bursts))
    ]
    burst_durns: list[float] = [
        spikes[bursts[n][1] - 1]  # type: ignore[call-overload]
        - spikes[bursts[n][0] - 1]  # type: ignore[call-overload]
        for n in range(len(bursts))
    ]
    burst_len_durn: list[list[float]] = bursts
    for idx, b_l_d in enumerate(burst_len_durn):
        b_l_d.extend([burst_lens[idx], burst_durns[idx]])
    # Only burst that qualify
    selected_bursts: list[list[float]] = [
        b
        for b in burst_len_durn
        # because IBI for first entry is None
        if (b[3] >= min_spikes_in_burst and b[4] >= min_burst_duration)
    ]

    if len(selected_bursts) == 0:
        return selected_bursts

    # Compute mean ISIs
    # recompute len & durn
    burst_lens = [
        selected_bursts[n][1] - selected_bursts[n][0] + 1
        for n in range(len(selected_bursts))
    ]
    burst_durns = [
        spikes[selected_bursts[n][1] - 1]  # type: ignore[call-overload]
        - spikes[selected_bursts[n][0] - 1]  # type: ignore[call-overload]
        for n in range(len(selected_bursts))
    ]
    mean_isis = [
        burst_durns[n] / (burst_lens[n] - 1) for n in range(len(selected_bursts))
    ]

    # Recompute IBI (since some bursts may have just been deleted)
    ibis: list[float] = []
    if len(selected_bursts) > 1:
        ibis = _calc_ibi(spikes, selected_bursts)  # type: ignore[assignment]
    # Replace current IBI with new IBI
    for ibi_num, _ in enumerate(ibis):
        selected_bursts[ibi_num][2] = ibis[ibi_num]

    # `si` is Spike Interval?  No idea what this is for
    # See https://github.com/igm-team/meaRtools/blob/master/meaRtools/R/maxinterval.R#L154
    si: int = 1
    for idx, sbs in enumerate(selected_bursts):
        sbs.extend([mean_isis[idx], si])

    return selected_bursts


if __name__ == "__main__":

    # pass
    SAMPLING_RATE: int = 1000
    MAX_BEGIN_ISI: float = 0.17
    MAX_END_ISI: float = 0.3
    MIN_BURST_DURATION: float = 0.01
    MIN_IBI: float = 0.4
    MIN_SPIKES_IN_BURST: int = 3

    SPIKE_TRAIN: list[float] = []
    with open("tests/data/cult.json", "r", encoding="utf-8") as f:
        SPIKE_TRAIN = json.load(f)

    # bursts: list[list[float]] = max_interval( SPIKE_TRAIN, MAX_BEGIN_ISI, MAX_END_ISI,
    #           MIN_IBI, MIN_BURST_DURATION, MIN_SPIKES_IN_BURST, SAMPLING_RATE)
    detected_bursts: list[list[float]] = _burst_detection(
        SPIKE_TRAIN, MAX_BEGIN_ISI, MAX_END_ISI, SAMPLING_RATE
    )
    merged_bursts: list[list[float]] = _merge_bursts(detected_bursts, MIN_IBI)
    burst_data: list[list[float]] = _quality_control(
        SPIKE_TRAIN,
        SAMPLING_RATE,
        merged_bursts,
        MIN_SPIKES_IN_BURST,
        MIN_BURST_DURATION,
    )
    # breakpoint()
    print(len(detected_bursts))
    print(detected_bursts[0])
    print(detected_bursts[-1])
