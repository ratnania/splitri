import numpy as np
from caid.cad_geometry import cad_geometry, cad_nurbs

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


def cad_geometry_levels(geo_i, geo_e, d):
    """
    constructs a set of (d-2) curves that lie between
    geo_i (interior) and geo_e (exterior)
    TODO: treate the NURBS case. Only the spline (weights=1) is treated here
    Args:
        geo_i (cad_geometry) : describing a spline curve
        geo_e (cad_geometry) : describing a spline curve
    Output:
        cad_geometry object that contains a list of curves, ordered from the
        interior to the exterior
    """
    geo = cad_geometry()
    nrb_i = geo_i[0]
    nrb_e = geo_e[0]

    U = np.linspace(0.,1.,d+1)
    for u in U[1:-1]:
        P_i = nrb_i.points
        P_e = nrb_e.points

        P   = u * P_i + (1.-u) * P_e

        cad_nrb = cad_nurbs(nrb_i.knots, P)
        geo.append(cad_nrb)

    return geo

def reduce_degree_spline(geo):
    """
    this routine creates a new spline of degree -1 than geo (a B-spline curve)
    """
