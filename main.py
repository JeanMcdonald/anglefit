import data
import fit
import plot
import plotpoly
import polynomialfit


d = data.collect("clean.root")

res = fit.fit_with_bkg(d[:4], 1)

print(res[0])
plot.plot(*res[1], res[0])
print(res[0])
#res = polynomialfit.fitpdf(d)

#res = fit.fit_with_bkg(d)
# print(res[0])
# plot.plot(*d, res)
#plotpoly.plot(*d, res)