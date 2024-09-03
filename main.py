import data
import fit
import plot
import plotpoly
import polynomialfit

d = data.collect("clean.root")

res = fit.fitpdf(d[:3])

plot.plot(*d[:3], res)
print(res)
#res = polynomialfit.fitpdf(d)

#res = fit.fit_with_bkg(d)
# print(res[0])
# plot.plot(*d, res)
#plotpoly.plot(*d, res)