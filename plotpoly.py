import numpy as np
import zfit
import matplotlib.pyplot as plt
from zfit import z
from polynomialfit import Polyfit
import zfit.z.numpy as znp
import math

# def fit(phi_array):
#     npphi = np.transpose(np.array(phi_array))
#     obs = zfit.Space("phi", limits=(-math.pi, math.pi))
#     data = zfit.Data.from_numpy(obs=obs, array=npphi)
#     # K1ss = zfit.Parameter("K1ss", 0.1, -10, 10)
#     # K1cc = K1ss#zfit.Parameter("K1cc", 0, -10, 10)
#     # K3s = zfit.Parameter("K3s", 0, -10, 10)
#     # K4s = zfit.Parameter("K4s", 0, -10, 10)
#     # phi = PhiPDF(obs=obs, K1ss=K1ss, K1cc=K1cc, K3s=K3s, K4s=K4s)
#     # minimizer = zfit.minimize.Minuit(tol=0.000010)
#     # nll = zfit.loss.UnbinnedNLL(model=phi, data=data)
#     # result = minimizer.minimize(nll)
#     # print(result)
#     data_plot = zfit.run(z.unstack_x(data))
#     plt.hist(data_plot, bins=500)
#     plt.savefig('phi')
#     plt.show()
#     def plot_model(model, data, name):
#         size_normal = len(phi_array)
#         nbins = 500
#         lower, upper = data.data_range.limit1d
#         x = np.linspace(lower, upper, num=nbins)
#         y = model.pdf(x) * size_normal / nbins * data.data_range.area()
#         plt.plot(x, y)
#
#         data_plot = zfit.run(z.unstack_x(data))
#         plt.hist(data_plot, bins=nbins)
#         plt.savefig(name)
#         plt.show()
#
#     #plot_model(phi, data, "phi")

def plot(thetal_array, thetak_array, phi_array, result):

    def plot_model(model, data, name):
        size_normal = len(thetal_array)
        nbins = 500
        lower, upper = data.data_range.limit1d
        x = np.linspace(lower, upper, num=nbins)
        y = model.pdf(x) * size_normal / nbins * data.data_range.area()
        plt.plot(x, y)

        data_plot = zfit.run(z.unstack_x(data))
        plt.hist(data_plot, bins=nbins)
        plt.savefig(name)
        plt.show()

    npthetal = np.transpose(np.array(thetal_array))
    obs = zfit.Space("costhetal", limits=(-1, 1))
    data = zfit.Data.from_numpy(obs=obs, array=npthetal)
    thetal = Polyfit(obs=obs, c0 = result.params['l0']['value'], c1=result.params['l1']['value'], c2=result.params['l2']['value'])
    plot_model(thetal, data, "costhetal")

    npthetak = np.transpose(np.array(thetak_array))
    obs = zfit.Space("costhetak", limits=(-1, 1))
    data = zfit.Data.from_numpy(obs=obs, array=npthetak)
    thetak = Polyfit(obs=obs, c0 = result.params['k0']['value'], c1=result.params['k1']['value'], c2=result.params['k2']['value'])
    plot_model(thetak, data, "costhetak")

    npphi = np.transpose(np.array(phi_array))
    obs = zfit.Space("phi", limits=(-math.pi, math.pi))
    data = zfit.Data.from_numpy(obs=obs, array=npphi)
    phi = Polyfit(obs=obs, c0 = result.params['p0']['value'], c1=result.params['p1']['value'], c2=result.params['p2']['value'])
    plot_model(phi, data, "phi")
