import numpy as np
from numpy.random import rand


def sample(alphabet, dist):
    """ This method produce a new discrete sample list by alphabet with probability
    distribution given by dist.
    The length of alphabet and dist must be same."""
    res = None
    cum_dist = np.cumsum(dist)
    r = rand()
    for i in range(len(dist)):
        if r < cum_dist[i]:
            res = alphabet[i]
            break

    return res