import matplotlib.pyplot as plt
import numpy as np


def draw_hist(ax, arr, bins, range):
    h, edges = np.histogram(arr, bins=bins, range=range)
    hden, _ = np.histogram(arr, bins=bins, range=range, density=True)
    e = np.sqrt(h)
    eden = e * hden.sum() / h.sum()
    centers = 0.5 * (edges[1:] + edges[:-1])
    ax.errorbar(centers, hden, yerr=eden, color="black", fmt=".", label="data")
    return h, hden, e, eden, edges


arr_data = np.loadtxt("tddp4_data.txt")
arr_mc = np.loadtxt("tddp4_mc.txt")

# dtime
plt.clf()
ax = plt.gca()
range = (0, 3.25)
draw_hist(ax, arr_data[:, -1], 100, range)
ax.hist(
    arr_mc[:, -3],
    bins=100,
    range=range,
    density=True,
    color="tab:blue",
    histtype="step",
    weights=arr_mc[:, -2] * arr_mc[:, -1],
)
ax.set_yscale("log")
ax.set_xlim(range[0], range[1])
ax.set_xlabel(r"$t(D^0)$ [ps]")
plt.savefig("tddp4_dtime.png")

# m12
plt.clf()
ax = plt.gca()
range = (0.25, 1.25)
draw_hist(ax, arr_data[:, 0], 100, range)
ax.hist(
    arr_mc[:, 0],
    bins=100,
    range=range,
    density=True,
    color="tab:blue",
    histtype="step",
    weights=arr_mc[:, -1] * arr_mc[:, -2],
)
ax.set_xlim(range[0], range[1])
ax.set_xlabel(r"$m(\pi\pi)$ [GeV]")
plt.savefig("tddp4_m12.png")

# m34
plt.clf()
ax = plt.gca()
range = (0.6, 1.6)
draw_hist(ax, arr_data[:, 1], 100, (0.6, 1.6))
ax.hist(
    arr_mc[:, 1],
    bins=100,
    range=range,
    density=True,
    color="tab:blue",
    histtype="step",
    weights=arr_mc[:, -1] * arr_mc[:, -2],
)
ax.set_xlim(range[0], range[1])
ax.set_xlabel(r"$m(K\pi)$ [GeV]")
plt.savefig("tddp4_m34.png")

# cos12
plt.clf()
ax = plt.gca()
range = (-1, 1)
draw_hist(ax, arr_data[:, 2], 100, range)
ax.hist(
    arr_mc[:, 2],
    bins=100,
    range=range,
    density=True,
    color="tab:blue",
    histtype="step",
    weights=arr_mc[:, -1] * arr_mc[:, -2],
)
ax.set_xlim(range[0], range[1])
ax.set_xlabel(r"$\cos\theta_{\pi\pi}$")
plt.savefig("tddp4_cos12.png")

# cos34
plt.clf()
ax = plt.gca()
range = (-1, 1)
draw_hist(ax, arr_data[:, 3], 100, range)
ax.hist(
    arr_mc[:, 3],
    bins=100,
    range=range,
    density=True,
    color="tab:blue",
    histtype="step",
    weights=arr_mc[:, -1] * arr_mc[:, -2],
)
ax.set_xlim(range[0], range[1])
ax.set_xlabel(r"$\cos\theta_{K\pi}$")
plt.savefig("tddp4_cos34.png")

# phi
plt.clf()
ax = plt.gca()
range = (-np.pi, np.pi)
draw_hist(ax, arr_data[:, 4], 100, range)
ax.hist(
    arr_mc[:, 4],
    bins=100,
    range=range,
    density=True,
    color="tab:blue",
    histtype="step",
    weights=arr_mc[:, -1] * arr_mc[:, -2],
)
ax.set_xlim(range[0], range[1])
ax.set_xlabel(r"$\phi$")
plt.savefig("tddp4_phi.png")
