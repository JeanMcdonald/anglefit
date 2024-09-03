import numpy as np
import zfit
import matplotlib.pyplot as plt
from zfit import z
import zfit.z.numpy as znp
import math
class ThetaLPDF(zfit.pdf.ZPDF):
    """1-dimensional PDF implementing the exp(alpha * x) shape."""
    _PARAMS = ("K1ss", "K1cc", "K1c")  # specify which parameters to take
    @zfit.supports()
    def _unnormalized_pdf(self, x, params):  # implement function
        data = x[0]  # axis 0
        K1ss = params["K1ss"]
        K1cc = params["K1cc"]
        K1c = params["K1c"]
        return 4*math.pi*(K1ss*(1 - znp.square(data)) + K1cc*znp.square(data) + K1c*data)

class ThetaLPDF(zfit.pdf.ZPDF):
    """1-dimensional PDF implementing the exp(alpha * x) shape."""
    _PARAMS = ("K1ss", "K1cc", "K1c")  # specify which parameters to take
    @zfit.supports()
    def _unnormalized_pdf(self, x, params):  # implement function
        data = x[0]  # axis 0
        K1ss = params["K1ss"]
        K1cc = params["K1cc"]
        K1c = params["K1c"]
        return 4*math.pi*(K1ss*(1 - znp.square(data)) + K1cc*znp.square(data) + K1c*data)

class ThetaKPDF(zfit.pdf.ZPDF):
    _PARAMS = ("K1ss", "K1cc", "K2ss", "K2cc")  # specify which parameters to take
    @zfit.supports()
    def _unnormalized_pdf(self, x, params):  # implement function
        data = x[0]  # axis 0
        K1ss = params["K1ss"]
        K1cc = params["K1cc"]
        K2ss = params["K2ss"]
        K2cc = params["K2cc"]
        return 2 * ((K1ss*4/3 + K1cc*2/3) + (K2ss*4/3 + K2cc*2/3)*data)

class PhiPDF(zfit.pdf.ZPDF):
    _PARAMS = ("K1ss", "K1cc", "K3s", "K4s")  # specify which parameters to take
    @zfit.supports()
    def _unnormalized_pdf(self, x, params):  # implement function
        data = x[0]  # axis 0
        K1ss = params["K1ss"]
        K1cc = params["K1cc"]
        K3s = params["K3s"]
        K4s = params["K4s"]
        return 2*(K1ss*4/3 + K1cc*2/3) + (math.pi**2) * (K3s*znp.sin(data) + K4s*znp.cos(data)) / 4
        #return (math.pi**2)/4 * (K3s*znp.sin(data) + K4s*znp.cos(data))

#
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


def plotdata(data, name):
    data_plot = zfit.run(z.unstack_x(data))
    plt.hist(data_plot, bins=100)
    plt.savefig(name)
    plt.show()


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
    thetal = ThetaLPDF(obs=obs, K1ss=result.params['K1ss']['value'], K1c=result.params['K1c']['value'], K1cc=result.params['K1cc']['value'])
    plot_model(thetal, data, "costhetal")

    npthetak = np.transpose(np.array(thetak_array))
    obs = zfit.Space("costhetak", limits=(-1, 1))
    data = zfit.Data.from_numpy(obs=obs, array=npthetak)
    thetak = ThetaKPDF(obs=obs, K1ss=result.params['K1ss']['value'], K1cc=result.params['K1cc']['value'], K2cc=result.params['K2cc']['value'], K2ss=result.params['K2ss']['value'])
    plot_model(thetak, data, "costhetak")

    npphi = np.transpose(np.array(phi_array))
    obs = zfit.Space("phi", limits=(-math.pi, math.pi))
    data = zfit.Data.from_numpy(obs=obs, array=npphi)
    phi = PhiPDF(obs=obs, K1ss=result.params['K1ss']['value'], K1cc=result.params['K1cc']['value'], K3s=result.params['K3s']['value'], K4s=result.params['K4s']['value'])
    plot_model(phi, data, "phi")

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


