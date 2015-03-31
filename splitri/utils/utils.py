import numpy as np



def nearest_int(u):
    iu = int(u)
    if u-iu < 0.5:
        return iu
    else:
        return iu+1

def iround(x):
    """
    iround(number) -> integer
    Round a number to the nearest integer.
    """
    return np.array([nearest_int(_x) for _x in x])
#    return np.array([int(round(_x) - .5) + (_x > 0) for _x in x])
