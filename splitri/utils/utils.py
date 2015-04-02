import numpy as np
from caid.cad_geometry import cad_geometry, cad_nurbs
from pigasus.utils.manager import context
from pigasus.fit.curfit import curfit
from pigasus.fit.curfit import compute_uk
from caid.cad_geometry import line

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


def construct_curve_from_points(x,y, px, \
                                knots=None, method='chord', alpha=0.1, nx=None):
    with context():
        #-----------------------------------
        def MakeConstraint(cond, face=None, value=None):
            if cond.lower() == "closed":
                constraint = {}
                constraint['patch_id_m'] = 0
                constraint['face_m']     = 0
                constraint['patch_id_s'] = 0
                constraint['face_s']     = 1
            if cond.lower() == "c0":
                constraint = {}
                constraint['patch_id_m'] = 0
                constraint['face_m']     = face
                constraint['type']       = "C0"
                constraint['values']     = [value]
            if cond.lower() == "c1":
                constraint = {}
                constraint['patch_id_m'] = 0
                constraint['face_m']     = face
                constraint['type']       = "C1"
                constraint['values']     = [value]
            return constraint
        #-----------------------------------

#        #-----------------------------------
#        print(("Approximation using " + method + " knots"))
#        print(("n " + str(nx) + " p ", str(px)))
#        #-----------------------------------

        #-----------------------------------
        # ...
        list_x = list(x) ; list_y = list(y)
        # ...

        # ...
        list_Q = list(zip(list_x, list_y))
        uk = compute_uk(list_Q, method=method)
        U1       = []
        U1      += list(uk)
        list_xk  = []     ; list_yk  = []
        list_xk += list_x ; list_yk += list_y

        lists_uk = [U1]
        lists_xk = [list_xk]
        lists_yk = [list_yk]
        # ...
        #-----------------------------------

        #-----------------------------------
        constraints = []

        # ... C0_COND
        constraint = MakeConstraint("C0", 0, [x[0], y[0]])
        constraints.append(constraint)

        constraint = MakeConstraint("C0", 1, [x[-1], y[-1]])
        constraints.append(constraint)
        #-----------------------------------

        #-----------------------------------
        geo = line(p=[px])
        if knots is not None:
            u = knots
            geo.refine(list_t=[u])
        else:
            ub = patch.knots[0][0]
            ue = patch.knots[0][-1]
            u = np.linspace(ub,ue,nx+1)[1:-1]
            geo.refine(list_t=[u])
        #-----------------------------------

        #-----------------------------------
        fit = curfit(geometry=geo, constraints=constraints, alpha=alpha)
        #-----------------------------------

        #-----------------------------------
        patch_id = 0
        xk = lists_xk[patch_id]
        yk = lists_yk[patch_id]

        geo = fit.construct([xk, yk], uk=lists_uk)
        #-----------------------------------

        fit.PDE.free()

        return geo

