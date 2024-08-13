import zfit
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcsetup
from zfit import z
import zfit.z.numpy as znp
import math



def integral(limits, params, model):
    lower = limits.limits[0][0]
    upper = limits.limits[1][0]
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
                - 0.5*(-K3sc*((1-x**2)**3/2)/3 + K3s*(x*(1-x**2)**0.5 + math.asin(x))) * (y*(1-y**2)**0.5 + math.asin(y)) * math.cos(z)
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

        return ((K1ss * sinthetal**2 + K1cc * costhetal**2 +K1c * costhetal)
                + (K2ss * sinthetal**2 + K2cc * costhetal**2 +K2c * costhetal) * costhetak
                + (K3sc * costhetal + K3s) * sinthetal * sinthetak * sinphi
                + (K4sc * costhetal + K4s) * sinthetal * sinthetak * cosphi)


def fitpdf(angledata):
    xobs = zfit.Space('thetal', (-1, 1))
    yobs = zfit.Space('thetak', (-1, 1))
    zobs = zfit.Space('phi', (-math.pi, math.pi))
    obs = xobs * yobs * zobs
    data = zfit.data.Data.from_numpy(array=np.transpose(np.array(angledata)), obs=obs)  # obs is automatically used as limits here.

    K1ss = zfit.Parameter("K1ss", 3, -20, 20)
    K1cc = zfit.Parameter("K1cc", 1, -20, 20)
    K1c = zfit.Parameter("K1c", 1, -20, 20)
    K2ss = zfit.Parameter("K2ss", 1, -20, 20)
    K2cc = zfit.Parameter("K2cc", 1, -20, 20)
    K2c = zfit.Parameter("K2c", 1, -20, 20)
    K3sc = zfit.Parameter("K3sc", 1, -20, 20)
    K3s = zfit.Parameter("K3s", 1, -20, 20)
    K4sc = zfit.Parameter("K4sc", 1, -20, 20)
    K4s = zfit.Parameter("K4s", 1, -20, 20)

    limit0 = zfit.Space(axes=0, lower=zfit.Space.ANY_LOWER, upper=zfit.Space.ANY_UPPER)
    limit1 = zfit.Space(axes=1, lower=zfit.Space.ANY_LOWER, upper=zfit.Space.ANY_UPPER)
    limit2 = zfit.Space(axes=2, lower=zfit.Space.ANY_LOWER, upper=zfit.Space.ANY_UPPER)
    limits = limit0 * limit1 * limit2

    angular_distribution = AngularDistribution(obs=obs, K1ss=K1ss, K1cc=K1cc, K1c=K1c, K2ss=K2ss, K2cc=K2cc, K2c=K2c, K3sc=K3sc, K3s=K3s, K4sc=K4sc, K4s=K4s)
    AngularDistribution.register_analytic_integral(integral, limits=limits)
    # angular_distribution.update_integration_options(tol=10e-5)
    # pdf = angular_distribution.pdf(data)
    minimizer = zfit.minimize.Minuit(tol=0.0001)
    nll = zfit.loss.UnbinnedNLL(model=angular_distribution, data=data)
    result = minimizer.minimize(nll)
    param_errors = result.hesse()
    print(result)
    # print(pdf[:20])



#
# nparr = np.array(thetal_array)
# obs1 = zfit.Space("thetal", limits=(-10, 10))
#
# data = zfit.Data.from_numpy(obs=obs1, array=nparr)
#
# data_plot = zfit.run(z.unstack_x(data))
#
# plt.hist(data_plot, bins=50)
#
# plt.savefig("testl.pdf")
# plt.show()
#
# nparr = np.array(thetak_array)
# obs1 = zfit.Space("thetak", limits=(-10, 10))
#
# data = zfit.Data.from_numpy(obs=obs1, array=nparr)
#
# data_plot = zfit.run(z.unstack_x(data))
#
# plt.hist(data_plot, bins=50)
#
# plt.savefig("testk.pdf")
# plt.show()
#
# nparr = np.array(phi_array)
# obs1 = zfit.Space("phi", limits=(-10, 10))
#
# data = zfit.Data.from_numpy(obs=obs1, array=nparr)
#
# data_plot = zfit.run(z.unstack_x(data))
#
# plt.hist(data_plot, bins=50)
#
# plt.savefig("testphi.pdf")
# plt.show()