"""
File: isin.py
Author: Guy Matz
Email: gmatz@matz.org
Github: https://www.github.com/guymatz/
Description: ISIn Burst Detection Algorithm?  This Python implementation
    was translated from code found in:
        https://pmc.ncbi.nlm.nih.gov/articles/instance/3915237/bin/Presentation1.PDF
"""

# import typing as t
import json
import math
import numbers
import numpy as np

# import numpy.typing as npt
from utils import _calc_ibi

"""
The ISIn Burst detection algorithm described in "Parameters for burst detection":
    https://pmc.ncbi.nlm.nih.gov/articles/PMC3915237/

Args:
    spike_times (list[float]): of spike times
    min_spikes_in_burst (int): Minimum number of spikes in a burst
    min_burst_duration (int): Minimum duration of burst in milliseconds

Returns:
    bursts: A list of lists continaing
        start
        end
        spikes (number of spikes in the burst)
        channels (number of channels in the burst)

    spike_burst_number: Burst # for each spike (-1 if spike is not in a burst)
Raises:
"""

def burst_detect(spike_train: list[float], N:int, ISI_n: int, spike_channels: list[float] = []) -> list[list[float]]:
    pass


def _get_lagged_diffs(spike_train: list[float], n: int) -> list[list[float]]:
    # Create a LoLoF with inf everywhere
    delta_times: list[list[float]] = [[math.inf] * n for _ in range(len(spike_train))]
    delta: float
    spike_a: int
    spike_b: int
    # A sliding window of spike_train[n + N] - spike_train[n]
    for j in range(n):
        for spike_num in range(n, len(spike_train)-1):
            spike_b = spike_num + j - 1
            spike_a = spike_num + j - n
            delta = spike_train[spike_b] - spike_train[spike_a]
            delta_times[spike_num-1][j] = delta
    return delta_times

def _filter_within_isin(spikes: list[float], delta_times: list[float], n: int, isin: int) -> list[list[float]]:
    """Get a List of spikes that fit ISIn criteria"""

    spike_burst_number: list[float] = [-1 for _ in range(len(spikes))]

    in_burst: bool = False
    burst_num_iterator: int = 0
    burst_num_assigned: int = -1
    burst_length: int = 0

    criteria: list[bool] =  [ True if min(x) < isin else False for x in delta_times]

    for idx, _ in enumerate(spikes):
        if in_burst:
            if criteria[idx]:
                in_burst = True
                burst_num_iterator += 1
                burst_num_assigned = burst_num_iterator
                burst_length = 1
            else:
                continue
        if not criteria[idx]:
            in_burst = False
            if burst_length < n:
                [burst_num] = -1
                spike_burst_number[:] = [-1 if x == burst_num_assigned else x for x in spike_burst_number]
                burst_num -= 1
            elif (spike_train[idx] - spike_train[idx - (n-1)]) > isin and burst_length >= n:

                burst_num_iterator += 1
                burst_num_assigned = burst_num_iterator
                burst_length = 1
            else:
                burst_length += 1
        spike_burst_number[idx] = burst_num_assigned

    # Assign Burst Information
    # Number of bursts:
    max_burst_number: int = max(spike_burst_number)

    # Burst.T_start: Burst Start Time
    # Burst.T_end: Burst End Time
    # Burst.S: Size (total Spikes)
    # Burst.C: Size (total channels)
    bursts: list[list[float]] = [[0.0] * 4] * max_burst_number
    for i in range(max_burst_number):
        spike_indexes_in_burst: list[int] = [i for i, x in enumerate(spike_burst_number) if x == 6]
        bursts[i][0] = spike_indexes_in_burst[0]    # start spike
        bursts[i][1] = spike_indexes_in_burst[-1]   # end spike
        bursts[i][2] = len(spike_indexes_in_burst)  # total spikes
        bursts[i][3] = len(spike_indexes_in_burst)  # total channels TODO

    return bursts
