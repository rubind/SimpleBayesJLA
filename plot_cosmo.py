import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import pickle
from IPython import embed
plt.rcParams["font.family"] = "serif"

def reflect(samps, othersamps = None, reflect_cut = 0.2):
    the_min = min(samps)
    the_max = max(samps)

    inds = np.where((samps < the_min*(1. - reflect_cut) + the_max*reflect_cut) & (samps > the_min))
    pad_samples = np.concatenate((samps, the_min - (samps[inds] - the_min)))
    if othersamps != None:
        pad_other = np.concatenate((othersamps, othersamps[inds]))

    inds = np.where((samps > the_min*reflect_cut + the_max*(1. - reflect_cut)) & (samps < the_max))
    pad_samples = np.concatenate((pad_samples, the_max + (the_max - samps[inds])))
    if othersamps != None:
        pad_other = np.concatenate((pad_other, othersamps[inds]))
        return pad_samples, pad_other

    return pad_samples

def reflect_2D(samps1, samps2, reflect_cut = 0.2):
    pad_samps1, pad_samps2 = reflect(samps1, samps2, reflect_cut = reflect_cut)
    pad_samps2, pad_samps1 = reflect(pad_samps2, pad_samps1, reflect_cut = reflect_cut)

    return pad_samps1, pad_samps2

def every_other_tick(ticks):
    """Matplotlib loves tick labels!"""

    labels = []

    for i in range(len(ticks) - 1):
        if i % 2 == len(ticks) % 2:
            labels.append(ticks[i])
        else:
            labels.append("")
    labels.append("")
    return labels

contours = [0.317311, 0.0455003]

grayscales = np.linspace(0.8, 0.4, len(contours))
colors = [[item]*3 for item in grayscales]

samples = pickle.load(open('./results.pickle', 'rb'))
om = samples[1]['Om']
ol = samples[1]['OL']

pad_om, pad_ol = reflect_2D(om, ol)
kernel = gaussian_kde(np.array([pad_om, pad_ol]), bw_method=0.1)
xvals, yvals = np.meshgrid(np.linspace(min(om), max(om), 100), np.linspace(min(ol), max(ol), 100))
eval_points = np.array([xvals.reshape(10000), yvals.reshape(10000)])
kernel_eval = kernel(eval_points)
kernel_eval /= kernel_eval.sum()
kernel_sort = np.sort(kernel_eval)
kernel_eval = np.reshape(kernel_eval, (100, 100))

kernel_cum = np.cumsum(kernel_sort)
levels = [kernel_sort[np.argmin(abs(kernel_cum - item))] for item in contours[::-1]]

ax = plt.axes()
ax.contourf(xvals, yvals, kernel_eval, levels = levels + [1], colors = colors)
ax.contour(xvals, yvals, kernel_eval, levels = levels, colors = 'k')
plt.show()
embed()