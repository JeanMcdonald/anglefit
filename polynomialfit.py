import zfit
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcsetup
from zfit import z
import zfit.z.numpy as znp
import math




def integral(limits, params, model):
    lower = limits.limits[0]
    upper = limits.limits[1]

    c0 = params['c0']
    c1 = params['c1']
    c2 = params['c2']

    def f(x):
        return c0*x + (c1*x**2)/2 + (2*c2*x**3)/3 - c2*x

    return f(upper) - f(lower)


class Polyfit(zfit.pdf.ZPDF):
    _PARAMS = ['c0', 'c1', 'c2']

    @zfit.supports()
    def _unnormalized_pdf(self, x, params):
        var = x[0]

        c0 = params['c0']
        c1 = params['c1']
        c2 = params['c2']

        return c0 + c1 * var + c2 * (2*var**2 - 1)


def create_distribution(obs):

    l0 = zfit.Parameter('l0', 0.3, -3, 3)
    l1 = zfit.Parameter('l1', 0.2, -3, 3)
    l2 = zfit.Parameter('l2', 0.1, -3, 3)

    k0 = zfit.Parameter('k0', 0.3, -3, 3)
    k1 = zfit.Parameter('k1', 0.2, -3, 3)
    k2 = zfit.Parameter('k2', 0.1, -3, 3)

    p0 = zfit.Parameter('p0', 0.3, -3, 3)
    p1 = zfit.Parameter('p1', 0.2, -3, 3)
    p2 = zfit.Parameter('p2', 0.1, -3, 3)

    limit0 = zfit.Space(axes=0, lower=zfit.Space.ANY_LOWER, upper=zfit.Space.ANY_UPPER)
    limit1 = zfit.Space(axes=1, lower=zfit.Space.ANY_LOWER, upper=zfit.Space.ANY_UPPER)
    limit2 = zfit.Space(axes=2, lower=zfit.Space.ANY_LOWER, upper=zfit.Space.ANY_UPPER)
    limits = limit0 * limit1 * limit2
    Polyfit.register_analytic_integral(integral, limits=limits)

    l = Polyfit(obs=obs, c0=l0, c1=l1, c2=l2)
    k = Polyfit(obs=obs, c0=k0, c1=k1, c2=k2)
    p = Polyfit(obs=obs, c0=p0, c1=p1, c2=p2)

    return l * k * p

def fitpdf(angledata):
    xobs = zfit.Space('thetal', (-1, 1))
    yobs = zfit.Space('thetak', (-1, 1))
    zobs = zfit.Space('phi', (-math.pi, math.pi))
    obs = xobs * yobs * zobs


    angular_distribution = create_distribution(angledata, obs)
    #angular_distribution.update_integration_options(tol=10e-5)
    minimizer = zfit.minimize.Minuit(tol=0.000010)
    data = zfit.data.Data.from_numpy(array=np.transpose(np.array(angledata)), obs=obs)

    nll = zfit.loss.UnbinnedNLL(model=angular_distribution, data=data)
    result = minimizer.minimize(nll)
    return (result, angular_distribution.create_projection_pdf(obs=(yobs * zobs)),
            angular_distribution.create_projection_pdf(obs=(xobs * zobs)),
            angular_distribution.create_projection_pdf(obs=(xobs * yobs)))





