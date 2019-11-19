
import numpy as np

#remove data jumps (steps/spikes)
def fixdata(x, jump):
    delta = np.concatenate([np.zeros(1), np.abs(np.diff(x, axis=1))])
    delta[delta>jump] = 0.
    y = x[0]+np.cumsum(delta)
    return y


# test problem - not sure if fixdata actually works correctly
#x = np.arange(1, 21)
#x[10] = 30
#print x
#print fixdata(x, 10)
