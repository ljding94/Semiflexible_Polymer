import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches
from scipy.optimize import curve_fit
from matplotlib import rc
from config_plot import *


def get_obs_data(folder, param):
    R2s, Rg2s, R2_errs, Rg2_errs = [], [], [], []
    Rxxs, Rxzs, Rxx_errs, Rxz_errs = [], [], [], []
    for L, kappa, f, gL in param:
        finfo = f"L{L}_kappa{kappa:.1f}_f{f:.2f}_gL{gL:.2f}"
        filename = f"{folder}/obs_{finfo}.csv"
        data = np.genfromtxt(filename, delimiter=',', skip_header=1)
        R2, Rg2, Rxx, Rxz = data[0, 11], data[0, 12], data[0, 13], data[0, 17]
        R2_err, Rg2_err, Rxx_err, Rxz_err = data[1, 11], data[1, 12], data[1, 13], data[1, 17]

        R2s.append(R2)
        Rg2s.append(Rg2)
        R2_errs.append(R2_err)
        Rg2_errs.append(Rg2_err)

        Rxxs.append(Rxx)
        Rxzs.append(Rxz)
        Rxx_errs.append(Rxx_err)
        Rxz_errs.append(Rxz_err)

    return np.array(R2s), np.array(Rg2s), np.array(Rxxs), np.array(Rxzs)


def get_tts_data(folder, param):
    tts = []
    tts_err = []
    spBs = []
    for L, kappa, f, gL in param:
        finfo = f"L{L}_kappa{kappa:.1f}_f{f:.2f}_gL{gL:.2f}"
        filename = f"{folder}/obs_{finfo}.csv"
        data = np.genfromtxt(filename, delimiter=',', skip_header=1)
        tt = data[3, 19:]
        tt_err = data[4, 19:]
        spB = data[5, 19:]
        tts.append(tt)
        tts_err.append(tt_err)
        spBs.append(spB)
    return tts, tts_err, spBs


def ax_fit(x, a):
    return a*x


def fit_l_persistence(spB, tts):
    popt, pcov = curve_fit(ax_fit, spB, np.log(tts))
    return -1/popt[0]


def calc_persistence_length(tts, spB):
    lp = fit_l_persistence(spB, tts)
    lp_theta = -1/np.log(tts[1])
    return lp, lp_theta


def get_lp_data(folder, param, fitn=5):
    tts, tts_err, spBs = get_tts_data(folder, param)
    lps = []
    lp_thetas = []
    for i in range(len(param)):
        lp, lp_theta = calc_persistence_length(tts[i][:fitn], spBs[i][:fitn])
        lps.append(lp)
        lp_thetas.append(lp_theta)
    return np.array(lps), np.array(lp_thetas)


def plot_obs_kappa(tex_lw=240.71031, ppi=72):

    fig = plt.figure(figsize=(tex_lw / ppi * 1, tex_lw / ppi * 0.9))
    plt.rc("text", usetex=True)
    plt.rc("text.latex", preamble=r"\usepackage{physics}")

    ax11 = fig.add_subplot(221)
    ax21 = fig.add_subplot(222)

    ax12 = fig.add_subplot(223)
    ax22 = fig.add_subplot(224)

    # plot tts vs. kappa
    ms = 4
    labelpad = -0.75
    folder = "../data/20240805"
    kappas = [2.0, 4.0, 8.0, 16.0]
    param = [(100, kappa, 0.0, 0.0) for kappa in kappas]

    tts, tts_err, spBs = get_tts_data(folder, param)
    markers = ["o", "x", "s", "+", "d"]
    pltn = 10
    for i in range(len(kappas)):
        ax11.semilogy(spBs[i][:pltn], tts[i][:pltn], marker=markers[i], ms=ms, mfc="None", ls="None", lw=1, label=fr"${kappas[i]:.0f}$")
    ax11.set_ylabel(r"$\left<\cos{\theta}(s)\right>$", fontsize=9, labelpad=labelpad)
    ax11.set_xlabel(r"$s$", fontsize=9, labelpad=labelpad)
    ax11.tick_params(which="both", direction="in", top="on", right="on", labelbottom=True, labelleft=True, labelsize=7)
    ax11.legend(title=r"$\kappa$", ncol=2, columnspacing=0.5, handlelength=0.5, handletextpad=0.2, frameon=False, fontsize=9)
    ax11.set_ylim(2e-2, 1.2)
    ax11.xaxis.set_major_locator(plt.MultipleLocator(2))
    ax11.xaxis.set_minor_locator(plt.MultipleLocator(1))
    # ax.yaxis.set_major_locator(plt.MultipleLocator(20))
    # ax.yaxis.set_minor_locator(plt.MultipleLocator(10))

    kappas = np.arange(2.0, 20.01, 2.0)
    L = 100
    param = [(L, kappa, 0.0, 0.0) for kappa in kappas]
    lps, lp_thetas = get_lp_data(folder, param)
    ax21.plot(kappas, lps, marker="x", ms=ms, mfc="None", ls="none", lw=1, label=r"$l_p$")
    ax21.plot(kappas, lp_thetas, marker="+", ms=ms, mfc="None", ls="none", lw=1, label=r"$l_{p,\theta}$")
    ax21.plot(kappas, kappas, marker="none", ms=ms, mfc="None", ls="-", lw=1, color="gray", label="theory")
    ax21.legend(ncol=2, columnspacing=0.5, handlelength=0.5, handletextpad=0.2, frameon=False, fontsize=9)
    ax21.set_ylabel(r"$l_p,l_{p,\theta}$", fontsize=9, labelpad=labelpad)
    ax21.set_xlabel(r"$\kappa$", fontsize=9, labelpad=labelpad)
    ax21.tick_params(which="both", direction="in", top="on", right="on", labelbottom=True, labelleft=True, labelsize=7)
    ax21.xaxis.set_major_locator(plt.MultipleLocator(4))
    ax21.xaxis.set_minor_locator(plt.MultipleLocator(2))
    ax21.yaxis.set_major_locator(plt.MultipleLocator(4))
    ax21.yaxis.set_minor_locator(plt.MultipleLocator(2))

    # plot R2 and Rg vs. kappa
    kappas = np.arange(2.0, 20.01, 2.0)
    for L in [100, 200]:
        param = [(L, kappa, 0.0, 0.0) for kappa in kappas]
        R2s, Rgs, Rxxs, Rxzs = get_obs_data(folder, param)
        # ax12.errorbar(kappas, np.array(R2s)/L, yerr=np.array(R2_errs)/L, marker="s", ls="none", ms=ms, mfc="None")
        ax12.plot(kappas, np.array(R2s)/L, marker="s", ls="none", ms=ms, mfc="None", label=fr"${L}$")

        lps = np.array(lps)
        t = np.exp(-1/kappas)
        R2_pL_theo = (1+t)/(1-t) + 2*t/L * (np.power(t, L)-1)/(1-t)**2
        if (L == 200):
            ax12.plot(kappas, R2_pL_theo, marker="none", ms=ms, mfc="None", ls="-", lw=1, color="gray", label="theory")
        else:
            ax12.plot(kappas, R2_pL_theo, marker="none", ms=ms, mfc="None", ls="-", lw=1, color="gray")
        ax22.plot(kappas, Rgs/L, ls="none", marker="v", ms=ms, mfc="None", label=fr"${L}$")

    # ax12.plot(kappas, 2*lps*(1-lps/L*(1-np.exp(-L/lps))), marker="none", ms=ms, mfc="None", ls="-", lw=1, color="gray", label=r"$2l_p(1-\frac{l_p}{L}(1-e^{-L/l_p}))$")
    # above equation is for R^2
    ax12.set_ylabel(r"$R^2/L$", fontsize=9, labelpad=labelpad)
    ax12.set_xlabel(r"$\kappa$", fontsize=9, labelpad=labelpad)
    ax12.legend(title=r"$L$", ncol=2, columnspacing=0.5, handlelength=0.5, handletextpad=0.2, frameon=False, fontsize=9)
    ax12.tick_params(which="both", direction="in", top="on", right="on", labelbottom=True, labelleft=True, labelsize=7)
    ax12.xaxis.set_major_locator(plt.MultipleLocator(4))
    ax12.xaxis.set_minor_locator(plt.MultipleLocator(2))
    ax22.yaxis.set_major_locator(plt.MultipleLocator(10))
    ax22.yaxis.set_minor_locator(plt.MultipleLocator(5))

    # ax22.errorbar(kappas, Rgs, yerr=Rg_errs, ls="none", marker="v", ms=ms, mfc="None")

    ax22.set_ylabel(r"$R_g^2/L$", fontsize=9, labelpad=labelpad)
    ax22.set_xlabel(r"$\kappa$", fontsize=9, labelpad=labelpad)
    ax22.legend(title=r"$L$", ncol=1, columnspacing=0.5, handlelength=0.5, handletextpad=0.2, frameon=False, fontsize=9)
    ax22.tick_params(which="both", direction="in", top="on", right="on", labelbottom=True, labelleft=True, labelsize=7)
    ax22.xaxis.set_major_locator(plt.MultipleLocator(4))
    ax22.xaxis.set_minor_locator(plt.MultipleLocator(2))
    ax22.yaxis.set_major_locator(plt.MultipleLocator(1))
    ax22.yaxis.set_minor_locator(plt.MultipleLocator(0.5))

    # add a,b,c,d,e,f,g,h
    annotation = [r"$(a)$", r"$(b)$", r"$(c)$", r"$(d)$"]
    for ax in [ax11, ax21, ax12, ax22]:
        ax.text(0.8, 0.1, annotation.pop(0), fontsize=9, transform=ax.transAxes)

    plt.tight_layout(pad=0.05)
    plt.savefig("./figures/obs_kappa.pdf", format="pdf")
    plt.show()
    plt.close()


def plot_obs_f_g(tex_lw=240.71031, ppi=72):
    # TODO add config and make it double column

    fig = plt.figure(figsize=(tex_lw / ppi * 2, tex_lw / ppi * 0.9))
    plt.rc("text", usetex=True)
    plt.rc("text.latex", preamble=r"\usepackage{physics}")

    ax11 = fig.add_subplot(241, projection="3d")
    ax12 = fig.add_subplot(242)
    ax13 = fig.add_subplot(243, sharex=ax12)
    ax14 = fig.add_subplot(244, sharex=ax12)

    ax21 = fig.add_subplot(245, projection="3d")
    ax22 = fig.add_subplot(246)
    ax23 = fig.add_subplot(247, sharex=ax22)
    ax24 = fig.add_subplot(248, sharex=ax22)

    ax11_2d = fig.add_subplot(241)
    ax11_2d.set_axis_off()
    ax21_2d = fig.add_subplot(245)
    ax21_2d.set_axis_off()

    ms = 3.5
    labelpad = 0.1
    # plot lp vs f
    folder = "../data/20240807"
    L = 200

    # plot config for f

    for f in [0.04, 0.12, 0.20, 0.28]:
        ax_plot_config(ax11, "../data/20240730", [100, 10.0, f, 0.00], -10, fr"${f:.2f}$")
    ax11.view_init(elev=32., azim=-75)
    ax11.quiver(40, 5, 5, 20, 0, 0, color="black", arrow_length_ratio=0.4)
    ax11.text(60, 5, 5, r"$\vu{x}$", fontsize=9)
    ax11.legend(title=r"$f$", loc="lower center", ncol=2, columnspacing=0.5, handlelength=0.5, handletextpad=0.2, frameon=False, fontsize=9)
    ax11.set_axis_off()

    kappas = [10.0]
    lss = ['-', "--"]
    markers = ["s", "o", "v"]
    color = ["tomato", "royalblue"]
    for i in range(len(kappas)):
        kappa = kappas[i]
        ls = lss[i]
        marker = markers[i]
        fs = np.arange(0.00, 0.301, 0.02)
        param = [(L, kappa, f, 0.0) for f in fs]
        lps, lp_thetas = get_lp_data(folder, param)
        if i == 0:
            ax12.plot(fs, lps, color="tomato", ls="none", marker=marker, ms=ms, mfc="None", label=r"$l_p$")
            ax12.plot(fs, lp_thetas, color="royalblue", ls="none", marker=marker, ms=ms, mfc="None", label=r"$l_{p,\theta}$")
        else:
            ax12.plot(fs, lps, color="tomato", ls="none", marker=marker, ms=ms, mfc="None")
            ax12.plot(fs, lp_thetas, color="royalblue", ls="none", marker=marker, ms=ms, mfc="None")
        #ax12.text(0, kappa-1.7, fr"$\kappa={kappa:.0f}$", fontsize=9)

    ax12.legend(title=rf"$\kappa={kappa:.0f}$",ncol=2, columnspacing=0.5, handlelength=0.5, handletextpad=0.2, frameon=False, fontsize=9)
    ax12.set_ylabel(r"$l_p,l_{p,\theta}$", fontsize=9, labelpad=labelpad)
    ax12.set_xlabel(r"$f$", fontsize=9, labelpad=labelpad)
    ax12.tick_params(which="both", direction="in", top="on", right="on", labelbottom=True, labelleft=True, labelsize=7)

    ax12.xaxis.set_major_locator(plt.MultipleLocator(0.2))
    ax12.xaxis.set_minor_locator(plt.MultipleLocator(0.1))
    ax12.yaxis.set_major_locator(plt.MultipleLocator(1))
    ax12.yaxis.set_minor_locator(plt.MultipleLocator(0.5))
    #ax12.set_ylim(2, 13)

    # plot R2 and Rg vs. f
    kappas = [5.0, 10.0]
    for i in range(len(kappas)):
        kappa = kappas[i]
        ls = lss[i]
        marker = markers[i]
        ls = "None"
        fs = np.arange(0.00, 0.301, 0.02)
        param = [(L, kappa, f, 0.0) for f in fs]
        R2, Rg2, Rxx, Rxz = get_obs_data(folder, param)
        # ax21.errorbar(fs, R2/L, yerr=R2_err/L, ls=ls, marker=marker, ms=ms, mfc="None", label=fr"${kappa:.0f}$")
        ax13.plot(fs, R2/L, ls=ls, marker=marker, ms=ms, mfc="None", label=fr"${kappa:.0f}$")
        # ax31.errorbar(fs, Rg, yerr=Rg_err, ls=ls, marker=marker, ms=ms, mfc="None", label=fr"${kappa:.0f}$")
        if kappa == 10:
            #ax14.plot(fs, Rg2/L, ls=ls, marker=marker, ms=ms, mfc="None", label=r"$R_g^2$")
            ax14.plot(fs, Rxx/Rg2, ls=ls, marker=marker, ms=ms, mfc="None", label=r"$xx$")
            ax14.plot(fs, Rxz/Rg2, ls=ls, marker=marker, ms=ms, mfc="None", label=r"$xz$")

    ax13.legend(title=r"$\kappa$", ncol=1, columnspacing=0.5, handlelength=0.5, handletextpad=0.2, frameon=False, fontsize=9)
    ax13.set_ylabel(r"$R^2/L$", fontsize=9, labelpad=labelpad)
    ax13.set_xlabel(r"$f$", fontsize=9, labelpad=labelpad)
    ax13.tick_params(which="both", direction="in", top="on", right="on", labelbottom=True, labelleft=True, labelsize=7)
    # ax13.xaxis.set_major_locator(plt.MultipleLocator(0.2))
    # ax13.xaxis.set_minor_locator(plt.MultipleLocator(0.1))
    ax13.yaxis.set_major_locator(plt.MultipleLocator(20))
    ax13.yaxis.set_minor_locator(plt.MultipleLocator(10))

    #title=r"$\kappa=10$",
    ax14.legend(title=r"$\mu$", ncol=1, columnspacing=0.5, handlelength=0.5, handletextpad=0.2, frameon=False, fontsize=9)
    ax14.set_ylabel(r"$R_{\mu}/R_g^2$", fontsize=9, labelpad=labelpad)
    ax14.set_xlabel(r"$f$", fontsize=9, labelpad=labelpad)
    ax14.tick_params(which="both", direction="in", top="on", right="on", labelbottom=True, labelleft=True, labelsize=7)

    ax14.xaxis.set_major_locator(plt.MultipleLocator(0.1))
    ax14.xaxis.set_minor_locator(plt.MultipleLocator(0.05))
    ax14.yaxis.set_major_locator(plt.MultipleLocator(0.2))
    ax14.yaxis.set_minor_locator(plt.MultipleLocator(0.1))
    #ax14.set_ylim(-0.25,6)

    # plot config vs g
    for gL in [0.40, 0.80, 1.20]:
        ax_plot_config(ax21, "../data/20240730", [100, 10.0, 0.00, gL], -20, fr"${gL:.1f}$")
    # ax11.view_init(elev=32., azim=-75)
    ax21.quiver(40, 5, 5, 20, 0, 0, color="black", arrow_length_ratio=0.4)
    ax21.text(60, 5, 5, r"$\vu{x}$", fontsize=9)
    ax21.quiver(40, 5, 5, 0, 0, 20, color="black", arrow_length_ratio=0.4)
    ax21.text(40, 5, 25, r"$\vu{z}$", fontsize=9)

    ax21.legend(title=r"$gL$", loc="lower center", ncol=2, columnspacing=0.5, handlelength=0.5, handletextpad=0.2, frameon=False, fontsize=9)
    ax21.set_axis_off()

    # plot lp vs g
    L = 200
    kappa = 10
    fs = [0.00]
    gLs = np.arange(0.00, 2.001, 0.20)
    lss = ['-', "--"]
    markers = ["s", "o", "v"]
    color = ["tomato", "royalblue"]
    for i in range(len(fs)):
        f = fs[i]
        ls = lss[i]
        ls = "none"
        marker = markers[i]
        param = [(L, kappa, f, gL) for gL in gLs]
        lps, lp_thetas = get_lp_data(folder, param)
        if i == 0:
            ax22.plot(gLs, lps, color="tomato", ls=ls, marker=marker, ms=ms, mfc="None", label=r"$l_p$")
            ax22.plot(gLs, lp_thetas, color="royalblue", ls=ls, marker=marker, ms=ms, mfc="None", label=r"$l_{p,\theta}$")
        else:
            ax22.plot(gLs, lps, color="tomato", marker=marker, ls=ls, ms=ms, mfc="None")
            ax22.plot(gLs, lp_thetas, color="royalblue", marker=marker, ls=ls, ms=ms, mfc="None")
        # ax22.text(0, kappa-1.5, fr"$f={f}$", fontsize=9)

    ax22.legend(title=r"$f=0$", ncol=2, columnspacing=0.5, handlelength=0.5, handletextpad=0.2, frameon=False, fontsize=9)
    ax22.set_ylabel(r"$l_p,l_{p,\theta}$", fontsize=9, labelpad=labelpad)
    ax22.set_xlabel(r"$gL$", fontsize=9, labelpad=labelpad)
    ax22.tick_params(which="both", direction="in", top="on", right="on", labelbottom=True, labelleft=True, labelsize=7)

    # ax22.xaxis.set_major_locator(plt.MultipleLocator(0.2))
    # ax22.xaxis.set_minor_locator(plt.MultipleLocator(0.1))
    ax22.yaxis.set_major_locator(plt.MultipleLocator(2))
    ax22.yaxis.set_minor_locator(plt.MultipleLocator(1))
    # ax22.set_ylim(2, 13)

    # plot R2 and Rg vs. g
    fs = [0.00, 0.20]
    for i in range(len(fs)):
        f = fs[i]
        ls = lss[i]
        ls = "None"
        marker = markers[i]
        param = [(L, kappa, f, gL) for gL in gLs]
        R2, Rg2, Rxx, Rxz = get_obs_data(folder, param)
        # ax23.errorbar(gLs, R2/L, yerr=R2_err/L, ls=ls, marker=marker, ms=ms, mfc="None", label=fr"${f:.1f}$")
        ax23.plot(gLs, R2/L, ls=ls, marker=marker, ms=ms, mfc="None", label=fr"${f:.1f}$")
        # ax24.errorbar(gLs, Rg, yerr=Rg_err, ls=ls,  marker=marker, ms=ms, mfc="None", label=fr"${f:.1f}$")
        if f == 0:
            #ax24.plot(gLs, Rg2/L, ls=ls, marker=marker, ms=ms, mfc="None", label=r"$R_g^2$")
            ax24.plot(gLs, Rxx/Rg2, ls=ls, marker=marker, ms=ms, mfc="None", label=r"$xx$")
            ax24.plot(gLs, Rxz/Rg2, ls=ls, marker=marker, ms=ms, mfc="None", label=r"$xz$")

    ax23.legend(title=r"$f$", ncol=1, columnspacing=0.5, handlelength=0.5, handletextpad=0.2, frameon=False, fontsize=9)
    ax23.set_ylabel(r"$R^2/L$", fontsize=9, labelpad=labelpad)
    ax23.set_xlabel(r"$gL$", fontsize=9, labelpad=labelpad)
    ax23.tick_params(which="both", direction="in", top="on", right="on", labelbottom=True, labelleft=True, labelsize=7)
    # ax23.xaxis.set_major_locator(plt.MultipleLocator(0.2))
    # ax23.xaxis.set_minor_locator(plt.MultipleLocator(0.1))
    ax23.yaxis.set_major_locator(plt.MultipleLocator(20))
    ax23.yaxis.set_minor_locator(plt.MultipleLocator(10))
    #title=r"$f=0$"
    ax24.legend(title=r"$\mu$", ncol=1, columnspacing=0.5, handlelength=0.5, handletextpad=0.2, frameon=False, fontsize=9)
    ax24.set_ylabel(r"$R_{\mu}/R_g^2$", fontsize=9, labelpad=labelpad)
    ax24.set_xlabel(r"$gL$", fontsize=9, labelpad=labelpad)
    ax24.tick_params(which="both", direction="in", top="on", right="on", labelbottom=True, labelleft=True, labelsize=7)
    #ax24.set_ylim(-0.25,7)
    ax24.xaxis.set_major_locator(plt.MultipleLocator(0.5))
    ax24.xaxis.set_minor_locator(plt.MultipleLocator(0.25))
    ax24.yaxis.set_major_locator(plt.MultipleLocator(0.2))
    ax24.yaxis.set_minor_locator(plt.MultipleLocator(0.1))

    annotation = [r"$(a)$", r"$(b)$", r"$(c)$", r"$(d)$", r"$(e)$", r"$(f)$", r"$(g)$", r"$(h)$"]

    axs = [ax11_2d, ax12, ax13, ax14, ax21_2d, ax22, ax23, ax24]
    for i in range(len(axs)):
        ax = axs[i]
        ax.text(0.82, 0.12, annotation[i], fontsize=9, transform=ax.transAxes)

    # ax12.text(0.82, 0.07, annotation[i], fontsize=9, transform=ax12.transAxes)
    # ax22.text(0.82, 0.07, annotation[i], fontsize=9, transform=ax22.transAxes)

    plt.tight_layout(pad=0.05)
    plt.savefig("./figures/obs_f_g.pdf", format="pdf")
    plt.show()
    plt.close()
