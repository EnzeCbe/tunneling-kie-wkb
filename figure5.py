"""Figure 5: KIE at 300 K over 117 barriers, tunneling only vs. ZPE in the barrier."""
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from rate_sweep import rates

Vs = np.linspace(0.5, 1.5, 9)
Ls = np.linspace(0.2, 0.8, 13)
res = {}
for zpe in (False, True):
    G = np.zeros((len(Vs), len(Ls), 2))
    for i, V in enumerate(Vs):
        for j, LA in enumerate(Ls):
            kH, kD = rates(V, LA, [300, 220], zpe=zpe)
            G[i, j] = kH / kD
    res[zpe] = G
np.save("kie_grid_tunneling.npy", res[False]); np.save("kie_grid_zpe.npy", res[True])

plt.rcParams.update({"font.size": 12})
fig, axs = plt.subplots(1, 2, figsize=(13, 5), dpi=300)
lo = np.log10(min(res[False][..., 0].min(), res[True][..., 0].min()))
hi = np.log10(max(res[False][..., 0].max(), res[True][..., 0].max()))
for ax, zpe, title in zip(axs, (False, True), ("Tunneling only", "Tunneling with ZPE in the barrier")):
    Z = np.log10(res[zpe][..., 0])
    im = ax.imshow(Z, origin="lower", aspect="auto", extent=[Ls[0], Ls[-1], Vs[0], Vs[-1]],
                   cmap="viridis", vmin=lo, vmax=hi)
    cs = ax.contour(Ls, Vs, Z, levels=[np.log10(7), np.log10(100)], colors="white", linewidths=1.4)
    ax.clabel(cs, fmt={np.log10(7): "7", np.log10(100): "100"})
    ax.plot(0.25, 1.5, "r*", ms=14)
    ax.set_xlabel("Barrier width parameter, $L$ (Å)")
    ax.set_ylabel("Classical barrier height, $V_{max}$ (eV)")
    ax.set_title(f"{title}, 300 K")
fig.colorbar(im, ax=axs, label=r"$\log_{10}$ KIE ($k_H/k_D$)")
fig.savefig("figure5_parameter_map.png", bbox_inches="tight")

t, z = res[False], res[True]
print("Tunneling only, 300 K: KIE>7 in %d of %d, KIE>100 in %d" % ((t[..., 0] > 7).sum(), t[..., 0].size, (t[..., 0] > 100).sum()))
print("Tunneling only, 0.5-0.875 eV with L>=0.4: KIE300 %.2f to %.2f" % (t[:4, 4:, 0].min(), t[:5, 4:, 0].max()))
print("ZPE in barrier, 300 K: 7-100 in %d, >100 in %d, <7 in %d" % (((z[..., 0] >= 7) & (z[..., 0] <= 100)).sum(), (z[..., 0] > 100).sum(), (z[..., 0] < 7).sum()))
print("ZPE in barrier, 300 K range: %.2f to %.3g" % (z[..., 0].min(), z[..., 0].max()))
