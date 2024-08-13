import data
import fit
import plot

d = data.collect("clean.root")

res = fit.fitpdf(d)
plot.plot(*d, res)