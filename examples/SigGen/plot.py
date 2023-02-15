#!/usr/bin/python
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

m = np.genfromtxt("../DP4/ToyMC.txt", skip_header=1)
sig = np.genfromtxt("SigGen.txt", skip_header=2)

m12 = m[:, 0]
m34 = m[:, 1]
c12 = m[:, 2]
c34 = m[:, 3]
phi = m[:, 4]

m12_1 = sig[:, 0]
m34_1 = sig[:, 1]
c12_1 = sig[:, 2]
c34_1 = sig[:, 3]
phi_1 = sig[:, 4]
w_1 = sig[:, 5]

weights1 = 1.0 / len(m12_1) * np.ones_like(m12_1)
weights2 = 1.0 / len(m12) * np.ones_like(m12)

bins = np.histogram(np.hstack((m12, m12_1)), bins=100)[1]
width = round((bins[-1] - bins[0]) / 100, 3)
plt.figure(figsize=(12, 6))
plt.hist(m12_1, bins, weights=weights1, alpha=1, label="GooFit", color="red")
plt.hist(m12, bins, weights=weights2, alpha=0.45, label="MINT", color="blue")
plt.ylabel(f"Entries/{width}", ha="right", y=1)
plt.xlabel(r"$m(\pi^+\pi^-)/GeV$", ha="right", x=1)
plt.minorticks_on()
plt.ticklabel_format(style="sci", axis="y", scilimits=(0, 0))
plt.legend(loc="best")
plt.savefig("SigGen_m12.png")

bins = np.histogram(np.hstack((m34, m34_1)), bins=100)[1]
width = round((bins[-1] - bins[0]) / 100, 3)
plt.figure(figsize=(12, 6))
plt.hist(m34_1, bins, weights=weights1, alpha=1, label="GooFit", color="red")
plt.hist(m34, bins, weights=weights2, alpha=0.45, label="MINT", color="blue")
plt.ylabel(f"Entries/{width}", ha="right", y=1)
plt.xlabel(r"$m(K^-\pi^-)/GeV$", ha="right", x=1)
plt.minorticks_on()
plt.ticklabel_format(style="sci", axis="y", scilimits=(0, 0))
plt.legend(loc="best")
plt.savefig("SigGen_m34.png")

bins = np.histogram(np.hstack((phi, phi_1)), bins=100)[1]
width = round((bins[-1] - bins[0]) / 100, 3)
plt.figure(figsize=(12, 6))
plt.hist(phi_1, bins, weights=weights1, alpha=1, label="GooFit", color="red")
plt.hist(phi, bins, weights=weights2, alpha=0.45, label="MINT", color="blue")
plt.ylabel(f"Entries/{width}", ha="right", y=1)
plt.xlabel(r"$\phi/ \degree $", ha="right", x=1)
plt.minorticks_on()
plt.ticklabel_format(style="sci", axis="y", scilimits=(0, 0))
plt.xlim(-np.pi, np.pi)
plt.legend(loc="best")
plt.savefig("SigGen_phi.png")

bins = np.histogram(np.hstack((c12, c12_1)), bins=100)[1]
width = round((bins[-1] - bins[0]) / 100, 3)
plt.figure(figsize=(12, 6))
plt.hist(c12_1, bins, weights=weights1, alpha=1, label="GooFit", color="red")
plt.hist(c12, bins, weights=weights2, alpha=0.45, label="MINT", color="blue")
plt.ylabel(f"Entries/{width}", ha="right", y=1)
plt.xlabel(r"$\cos(\theta_{12})/ \degree $", ha="right", x=1)
plt.minorticks_on()
plt.ticklabel_format(style="sci", axis="y", scilimits=(0, 0))
plt.legend(loc="best")
plt.savefig("SigGen_c12.png")

bins = np.histogram(np.hstack((c34, c34_1)), bins=100)[1]
width = round((bins[-1] - bins[0]) / 100, 3)
plt.figure(figsize=(12, 6))
plt.hist(c34_1, bins, weights=weights1, alpha=1, label="GooFit", color="red")
plt.hist(c34, bins, weights=weights2, alpha=0.45, label="MINT", color="blue")
plt.ylabel(f"Entries/{width}", ha="right", y=1)
plt.xlabel(r"$\cos(\theta_{34})/ \degree $", ha="right", x=1)
plt.minorticks_on()
plt.ticklabel_format(style="sci", axis="y", scilimits=(0, 0))
plt.legend(loc="best")
plt.savefig("SigGen_c34.png")
