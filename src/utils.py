"""Utility functions

Translated from https://github.com/igm-team/meaRtools/blob/master/meaRtools/R/burst_stats.R
"""


def _calc_ibi(spikes: list[float], bursts: list[list[float]]) -> list[float | None]:
    """
    Utiliy function for computing Inter-Burst Intervals

    Args:
        spikes: list of spike times
        bursts: list of burst each containing, start index, stop index, length & duration

    Returns:
        ibis: List of Inter-Burst Intervals

    Raises:
    """
    if len(bursts) < 2:
        # non IBIs
        return []

    ibis: list[float | None] = []
    burst_ends: list[float] = [b[1] for b in bursts]

    start_spikes: list[int] = [int(b[0]) - 1 for b in bursts[1:]]
    end_spikes: list[int] = [int(b) - 1 for b in burst_ends[:-1]]

    assert len(start_spikes) == len(
        end_spikes
    ), "start_spikes and end_spikes need to be the same length"

    ibis = [spikes[start_spikes[n]] - spikes[end_spikes[n]] for n in range(len(start_spikes))]
    # Prepend a None to the ibis list, since the first burst has no ibi
    ibis.insert(0, None)
    return ibis


if __name__ == "__main__":
    pass
