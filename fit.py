import zfit
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcsetup
from zfit import z
import zfit.z.numpy as znp
import math

import polynomialfit


def integral(limits, params, model):
    lower = limits.limits[0][0]
    upper = limits.limits[1][0]

    print(lower, upper, "AAAA")
    K1ss = params['K1ss']
    K1cc = params['K1cc']
    K1c = params['K1c']
    K2ss = params['K2ss']
    K2cc = params['K2cc']
    K2c = params['K2c']
    K3sc = params['K3sc']
    K3s = params['K3s']
    K4sc = params['K4sc']
    K4s = params['K4s']



    def f(x, y, z):
        return ((K1ss*(x - (x**3)/3) + K1cc*(x**3)/3 + K1c*(x**2)/2)*y*z
                 + (K2ss*(x - (x**3)/3) + K2cc*(x**3)/3 + K2c*(x**2)/2)*(y**2)*z/2
                 - (1/12)*((2*K3sc*(x**2 - 1) + 3*K3s*x) * (1 - x**2)**0.5 + 3*math.asin(x)) * (y*(1-y**2)**0.5 + math.asin(y)) * math.cos(z)
                 + (1/12)*((2*K4sc*(x**2 - 1) + 3*K4s*x) * (1 - x**2)**0.5 + 3*math.asin(x)) * (y*(1-y**2)**0.5 + math.asin(y)) * math.sin(z)
                + 0.5*(-K4sc*((1-x**2)**3/2)/3 + K4s*(x*(1-x**2)**0.5 + math.asin(x))) * (y*(1-y**2)**0.5 + math.asin(y)) * math.sin(z)
        )

    return f(*upper) - f(*lower)
class AngularDistribution(zfit.pdf.ZPDF):
    _PARAMS = ['K1ss', 'K1cc', 'K1c', 'K2ss', 'K2cc', 'K2c', 'K3sc', 'K3s', 'K4sc', 'K4s']

    @zfit.supports()
    def _unnormalized_pdf(self, x, params):
        costhetal = x[0]
        costhetak = x[1]
        cosphi = znp.cos(x[2])
        sinthetal = znp.sqrt(1 - znp.square(costhetal))
        sinthetak = znp.sqrt(1 - znp.square(costhetak))
        sinphi = znp.sin(x[2])

        K1ss = params['K1ss']
        K1cc = params['K1cc']
        K1c = params['K1c']
        K2ss = params['K2ss']
        K2cc = params['K2cc']
        K2c = params['K2c']
        K3sc = params['K3sc']
        K3s = params['K3s']
        K4sc = params['K4sc']
        K4s = params['K4s']

        return ((K1ss * sinthetal**2 + K1cc * costhetal**2 + K1c * costhetal)
                + (K2ss * sinthetal**2 + K2cc * costhetal**2 + K2c * costhetal) * costhetak
                + (K3sc * costhetal + K3s) * sinthetal * sinthetak * sinphi
                + (K4sc * costhetal + K4s) * sinthetal * sinthetak * cosphi)

def create_distribution(angledata, obs):
    K1ss = zfit.Parameter("K1ss", 0.1, -10, 10)
    K1cc = zfit.Parameter("K1cc", 0, -10, 10)
    K1c = zfit.Parameter("K1c", 0, -10, 10)
    K2ss = zfit.Parameter("K2ss", 0, -10, 10)
    K2cc = zfit.Parameter("K2cc", 0, -10, 10)
    K2c = zfit.Parameter("K2c", 0, -10, 10)
    K3sc = zfit.Parameter("K3sc", 0, -10, 10)
    K3s = zfit.Parameter("K3s", 0, -10, 10)
    K4sc = zfit.Parameter("K4sc", 0, -10, 10)
    K4s = zfit.Parameter("K4s", 0, -10, 10)

    return AngularDistribution(obs=obs, K1ss=K1ss, K1cc=K1cc, K1c=K1c, K2ss=K2ss, K2cc=K2cc, K2c=K2c, K3sc=K3sc, K3s=K3s, K4sc=K4sc, K4s=K4s)
    #AngularDistribution.register_analytic_integral(integral, limits=limits)

def fitpdf(angledata):
    xobs = zfit.Space('thetal', (-1, 1))
    yobs = zfit.Space('thetak', (-1, 1))
    zobs = zfit.Space('phi', (-math.pi, math.pi))
    obs = xobs * yobs * zobs
    data = zfit.data.Data.from_numpy(array=np.transpose(np.array(angledata)), obs=obs)
    angular_distribution = create_distribution(angledata, obs)
    angular_distribution.update_integration_options(tol=10e-5)
    minimizer = zfit.minimize.Minuit(tol=0.000010)
    nll = zfit.loss.UnbinnedNLL(model=angular_distribution, data=data)
    result = minimizer.minimize(nll)
    return result


def fit_with_bkg(angledata):
    xobs = zfit.Space('thetal', (-1, 1))
    yobs = zfit.Space('thetak', (-1, 1))
    zobs = zfit.Space('phi', (-math.pi, math.pi))
    obs = xobs * yobs * zobs
    data = zfit.data.Data.from_numpy(array=np.transpose(np.array(angledata)), obs=obs)

    signal = create_distribution(angledata, obs)
    background = polynomialfit.create_distribution(angledata, obs)

    frac = zfit.Parameter('frac', 0.9, 0, 1)

    total_pdf = zfit.pdf.SumPDF([signal, background], [frac])
    minimizer = zfit.minimize.Minuit(tol=0.0001)
    nll = zfit.loss.UnbinnedNLL(model=total_pdf, data=data)
    result = minimizer.minimize(nll)
    return result
