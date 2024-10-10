import zfit
import numpy as np
import zfit.z.numpy as znp
import math

import plot
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



    # def f(x, y, z):
        # return ((K1ss*(x - (x**3)/3) + K1cc*(x**3)/3 + K1c*(x**2)/2)*y*z
        #          + (K2ss*(x - (x**3)/3) + K2cc*(x**3)/3 + K2c*(x**2)/2)*(y**2)*z/2
        #          - (1/12)*((2*K3sc*(x**2 - 1) + 3*K3s*x) * (1 - x**2)**0.5 + 3*math.asin(x)) * (y*(1-y**2)**0.5 + math.asin(y)) * math.cos(z)
        #          + (1/12)*((2*K4sc*(x**2 - 1) + 3*K4s*x) * (1 - x**2)**0.5 + 3*math.asin(x)) * (y*(1-y**2)**0.5 + math.asin(y)) * math.sin(z)
        #         + 0.5*(-K4sc*((1-x**2)**3/2)/3 + K4s*(x*(1-x**2)**0.5 + math.asin(x))) * (y*(1-y**2)**0.5 + math.asin(y)) * math.sin(z)
        # )
    return (4*K1ss/3 + 2*K1cc/3)*4*math.pi

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

def create_distribution(obs):
    K1ss = zfit.Parameter("K1ss", 5, -10, 10)
    K1cc = zfit.Parameter("K1cc", 0, -10, 10)
    K1c = zfit.Parameter("K1c", 0, -10, 10)
    K2ss = zfit.Parameter("K2ss", 0, -10, 10)
    K2cc = zfit.Parameter("K2cc", 0, -10, 10)
    K2c = zfit.Parameter("K2c", 0, -10, 10)
    K3sc = zfit.Parameter("K3sc", 0, -10, 10)
    K3s = zfit.Parameter("K3s", 0, -10, 10)
    K4sc = zfit.Parameter("K4sc", 0, -10, 10)
    K4s = zfit.Parameter("K4s", 0, -10, 10)
    limit0 = zfit.Space(axes=0, lower=zfit.Space.ANY_LOWER, upper=zfit.Space.ANY_UPPER)
    limit1 = zfit.Space(axes=1, lower=zfit.Space.ANY_LOWER, upper=zfit.Space.ANY_UPPER)
    limit2 = zfit.Space(axes=2, lower=zfit.Space.ANY_LOWER, upper=zfit.Space.ANY_UPPER)
    limits = limit0 * limit1 * limit2
    AngularDistribution.register_analytic_integral(integral, limits=limits)
    return AngularDistribution(obs=obs, K1ss=K1ss, K1cc=K1cc, K1c=K1c, K2ss=K2ss, K2cc=K2cc, K2c=K2c, K3sc=K3sc, K3s=K3s, K4sc=K4sc, K4s=K4s)

def fitpdf(angledata):
    xobs = zfit.Space('thetal', (-1, 1))
    yobs = zfit.Space('thetak', (-1, 1))
    zobs = zfit.Space('phi', (-math.pi, math.pi))
    obs = xobs * yobs * zobs
    data = zfit.data.Data.from_numpy(array=np.transpose(np.array(angledata)), obs=obs)
    angular_distribution = create_distribution(angledata, obs)
    #angular_distribution.update_integration_options(tol=10e-5)
    minimizer = zfit.minimize.Minuit(tol=0.000010)
    nll = zfit.loss.UnbinnedNLL(model=angular_distribution, data=data)
    result = minimizer.minimize(nll)
    return result

def crystalball(mass_obs):
    mu = zfit.Parameter("mu", 5.5, 5.3, 5.7)
    sigma = zfit.Parameter("sigma", 0.05, 0, 0.5)
    alpha = zfit.Parameter("alpha", 5, 1, 7)
    n = zfit.Parameter("n", 4, -5, 10)
    sigma2 = zfit.Parameter("sigma2", 0.05, 0, 0.5)
    alpha2 = zfit.Parameter("alpha2", 5, 1, 7)
    n2 = zfit.Parameter("n2", 4, -5, 10)
    f = zfit.Parameter("f", 0.6, 0, 1)
    crystalball1 = zfit.pdf.CrystalBall(mu=mu, sigma=sigma, obs=mass_obs, alpha=alpha, n=n)
    crystalball2 = zfit.pdf.CrystalBall(mu=mu, sigma=sigma2, obs=mass_obs, alpha=alpha2, n=n2)
    crystalball = zfit.pdf.SumPDF([crystalball1, crystalball2], [f])
    return crystalball

def mass_background(mass_obs):
    lam = zfit.Parameter('lambda', -0.1, -5, 0)
    return zfit.pdf.Exponential(lam=lam, obs=mass_obs)

def fit_with_bkg(angledata, bkg_frac):
    xobs = zfit.Space('thetal', (-1, 1))
    yobs = zfit.Space('thetak', (-1, 1))
    zobs = zfit.Space('phi', (-math.pi, math.pi))
    mass_obs = zfit.Space('mass', (5.5, 6))
    obs = xobs * yobs * zobs
    num_events = int(len(angledata[0])*bkg_frac)
    thetal_bkg = np.random.uniform(size=num_events, low=-1, high=1)
    thetak_bkg = np.random.uniform(size=num_events, low=-1, high=1)
    phi_bkg = np.random.uniform(size=num_events, low=-math.pi, high=math.pi)
    mass_bkg = np.random.exponential(scale=1/8, size=num_events) + 5.5
    bkg_array = np.transpose(np.array([thetal_bkg, thetak_bkg, phi_bkg, mass_bkg]))

    data = zfit.data.Data.from_numpy(array=np.concatenate((np.transpose(np.array(angledata)), bkg_array)), obs=obs*mass_obs)#array=np.transpose(np.array(angledata)), obs=obs*mass_obs)

        #array=np.concatenate((np.transpose(np.array(angledata)), bkg_array)), obs=obs*mass_obs)

    # plot.plotdata(data[0], '1')
    # plot.plotdata(data[1], '2')
    # plot.plotdata(data[2], '3')
    # plot.plotdata(data[3], '4')

    signal = create_distribution(obs)*crystalball(mass_obs)
    background = polynomialfit.create_distribution(obs)*mass_background(mass_obs)

    frac = zfit.Parameter('frac', 0.5, 0, 1)
    total_pdf = zfit.pdf.SumPDF([signal, background], [frac], obs=obs*mass_obs)
    minimizer = zfit.minimize.Minuit(tol=0.0001)
    nll = zfit.loss.UnbinnedNLL(model=total_pdf, data=data)
    #nll = zfit.loss.UnbinnedNLL(model=signal, data=data)

    result = minimizer.minimize(nll)
    return result, data