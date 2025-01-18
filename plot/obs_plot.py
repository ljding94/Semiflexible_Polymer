import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches
from scipy.optimize import curve_fit
from matplotlib import rc
from config_plot import *

from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from Sq_plot import *

def get_obs_data(folder, param):
    Xs, Ys, Zs, R2s, Rg2s, R2_errs, Rg2_errs = [], [], [], [], [], [], []
    X_errs, Rxxs, Rxzs, Rxx_errs, Rxz_errs = [], [], [], [], []
    for L, kappa, f, gL in param:
        finfo = f"L{L}_kappa{kappa:.1f}_f{f:.2f}_gL{gL:.2f}"
        filename = f"{folder}/obs_{finfo}.csv"
        data = np.genfromtxt(filename, delimiter=',', skip_header=1)
        X, Y, Z, R2, Rg2, Rxx, Rxz = data[0, 7], data[0, 8], data[0, 9], data[0, 11], data[0, 12], data[0, 13], data[0, 17]
        X_err, R2_err, Rg2_err, Rxx_err, Rxz_err = data[1, 7], data[1, 11], data[1, 12], data[1, 13], data[1, 17]

        Xs.append(X)
        Ys.append(Y)
        Zs.append(Z)
        R2s.append(R2)
        Rg2s.append(Rg2)
        R2_errs.append(R2_err)
        Rg2_errs.append(Rg2_err)

        X_errs.append(X_err)
        Rxxs.append(Rxx)
        Rxzs.append(Rxz)
        Rxx_errs.append(Rxx_err)
        Rxz_errs.append(Rxz_err)

    return np.array(Xs), np.array(Ys), np.array(Zs), np.array(R2s), np.array(Rg2s), np.array(Rxxs), np.array(Rxzs)


def get_obs_data_g(folder, param):
    Xs, Ys, Zs, XsignZs, ZsignXs, R2s, Rg2s, R2_errs, Rg2_errs = [], [], [], [], [], [], [], [], []
    X_errs, Rxxs, Rxzs, Rxx_errs, Rxz_errs = [], [], [], [], []
    for L, kappa, f, gL in param:
        finfo = f"L{L}_kappa{kappa:.1f}_f{f:.2f}_gL{gL:.2f}"
        filename = f"{folder}/obs_{finfo}.csv"
        data = np.genfromtxt(filename, delimiter=',', skip_header=1)
        X, Y, Z, XsignZ, ZsignX, R2, Rg2, Rxx, Rxz = data[0, 7], data[0, 8], data[0, 9], data[0, 10], data[0, 11], data[0, 13], data[0, 14], data[0, 15], data[0, 19]
        X_err, Y_err, Z_err, XsignZ_err, ZsignX_err, R2_err, Rg2_err, Rxx_err, Rxz_err = data[1, 7], data[1, 8], data[1, 9], data[1, 10], data[1, 11], data[1, 13], data[1, 14], data[1, 15], data[1, 19]

        Xs.append(X)
        Ys.append(Y)
        Zs.append(Z)
        XsignZs.append(XsignZ)
        ZsignXs.append(ZsignX)

        R2s.append(R2)
        Rg2s.append(Rg2)
        R2_errs.append(R2_err)
        Rg2_errs.append(Rg2_err)

        X_errs.append(X_err)
        Rxxs.append(Rxx)
        Rxzs.append(Rxz)
        Rxx_errs.append(Rxx_err)
        Rxz_errs.append(Rxz_err)

    return np.array(Xs), np.array(Ys), np.array(Zs), np.array(XsignZs), np.array(ZsignXs), np.array(R2s), np.array(Rg2s), np.array(Rxxs), np.array(Rxzs)


def get_obs_data_err_g(folder, param):
    X_errs, Y_errs, Z_errs, XsignZ_errs, ZsignX_errs, R2_errs, Rg2_errs, Rxx_errs, Rxz_errs = [], [], [], [], [], [], [], [], []
    for L, kappa, f, gL in param:
        finfo = f"L{L}_kappa{kappa:.1f}_f{f:.2f}_gL{gL:.2f}"
        filename = f"{folder}/obs_{finfo}.csv"
        data = np.genfromtxt(filename, delimiter=',', skip_header=1)
        X_err, Y_err, Z_err, XsignZ_err, ZsignX_err, R2_err, Rg2_err, Rxx_err, Rxz_err = data[1, 7], data[1, 8], data[1, 9], data[1, 10], data[1, 11], data[1, 13], data[1, 14], data[1, 15], data[1, 19]

        X_errs.append(X_err)
        Y_errs.append(Y_err)
        Z_errs.append(Z_err)
        XsignZ_errs.append(XsignZ_err)
        ZsignX_errs.append(ZsignX_err)
        R2_errs.append(R2_err)
        Rg2_errs.append(Rg2_err)
        Rxx_errs.append(Rxx_err)
        Rxz_errs.append(Rxz_err)

    return np.array(X_errs), np.array(Y_errs), np.array(Z_errs), np.array(XsignZ_errs), np.array(ZsignX_errs), np.array(R2_errs), np.array(Rg2_errs), np.array(Rxx_errs), np.array(Rxz_errs)


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


def get_tts_data_g(folder, param):
    tts = []
    tts_err = []
    spBs = []
    for L, kappa, f, gL in param:
        finfo = f"L{L}_kappa{kappa:.1f}_f{f:.2f}_gL{gL:.2f}"
        filename = f"{folder}/obs_{finfo}.csv"
        data = np.genfromtxt(filename, delimiter=',', skip_header=1)
        tt = data[3, 21:]
        tt_err = data[4, 21:]
        spB = data[5, 21:]
        tts.append(tt)
        tts_err.append(tt_err)
        spBs.append(spB)
    return tts, tts_err, spBs


def ax_fit(x, a):
    return a*x


def fit_l_persistence(spB, tts):
    popt, pcov = curve_fit(ax_fit, spB, np.log(tts))
    perr = np.sqrt(np.diag(pcov))
    return -1/popt[0], popt[0]**2*perr[0]


def calc_persistence_length(tts, spB):
    lp, lp_err = fit_l_persistence(spB, tts)
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


def get_lp_data_g(folder, param, fitn=5):
    tts, tts_err, spBs = get_tts_data_g(folder, param)
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
    labelpad = 0.0
    folder = "../data/20240807"
    kappas = [2.0, 4.0, 8.0, 16.0]
    param = [(200, kappa, 0.0, 0.0) for kappa in kappas]

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
    L = 200
    param = [(L, kappa, 0.0, 0.0) for kappa in kappas]
    lps, lp_thetas = get_lp_data(folder, param)
    ax21.plot(kappas, lps, marker="x", ms=ms, mfc="None", ls="none", lw=1, label=r"$l_p$")
    ax21.plot(kappas, lp_thetas, marker="+", ms=ms, mfc="None", ls="none", lw=1, label=r"$l_{p,\theta}$")
    ax21.plot(kappas, kappas, marker="none", ms=ms, mfc="None", ls="-", lw=1, color="gray")  # label="theory")
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
    for L in [100, 200, 400]:
        param = [(L, kappa, 0.0, 0.0) for kappa in kappas]
        X, Y, Z, R2s, Rgs, Rxxs, Rxzs = get_obs_data(folder, param)
        # ax12.errorbar(kappas, np.array(R2s)/L, yerr=np.array(R2_errs)/L, marker="s", ls="none", ms=ms, mfc="None")
        ax12.plot(kappas, np.array(R2s)/L, marker="s", ls="none", ms=ms, mfc="None", label=fr"${L}$")

        lps = np.array(lps)
        t = np.exp(-1/kappas)
        R2_pL_theo = (1+t)/(1-t) + 2*t/L * (np.power(t, L)-1)/(1-t)**2
        if (L == 200):
            line_theory, = ax12.plot(kappas, R2_pL_theo, marker="none", ms=ms, mfc="None", ls="-", lw=1, color="gray")
        else:
            ax12.plot(kappas, R2_pL_theo, marker="none", ms=ms, mfc="None", ls="-", lw=1, color="gray")
        ax22.plot(kappas, Rgs/L, ls="none", marker="v", ms=ms, mfc="None", label=fr"${L}$")

    # ax12.plot(kappas, 2*lps*(1-lps/L*(1-np.exp(-L/lps))), marker="none", ms=ms, mfc="None", ls="-", lw=1, color="gray", label=r"$2l_p(1-\frac{l_p}{L}(1-e^{-L/l_p}))$")
    # above equation is for R^2
    ax12.set_ylabel(r"$R^2/L$", fontsize=9, labelpad=labelpad)
    ax12.set_xlabel(r"$\kappa$", fontsize=9, labelpad=labelpad)
    # theory_legend = ax12.legend(handles=[line_theory], labels=[r"theory"], loc="upper center", ncol=1, columnspacing=0.5, handlelength=0.5, handletextpad=0.2, frameon=False, fontsize=9)
    # ax12.add_artist(theory_legend)
    ax12.legend(title=r"$L$", ncol=1, columnspacing=0.5, handlelength=0.5, handletextpad=0.2, frameon=False, fontsize=9)

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

    plt.tight_layout(pad=0.2)
    #plt.subplots_adjust(wspace=0.2, hspace=0.25)
    plt.savefig("./figures/obs_kappa.pdf", format="pdf")
    plt.savefig("./figures/obs_kappa.png", format="png", dpi=300)
    plt.show()
    plt.close()

# stretching/tension


def plot_obs_f(tex_lw=240.71031, ppi=72):
    fig = plt.figure(figsize=(tex_lw / ppi * 1.0, tex_lw / ppi * 0.9))
    plt.rc("text", usetex=True)
    plt.rc("text.latex", preamble=r"\usepackage{physics}")

    ax11 = fig.add_subplot(221, projection="3d")
    ax12 = fig.add_subplot(222)
    ax21 = fig.add_subplot(223, sharex=ax12)
    ax22 = fig.add_subplot(224, sharex=ax12)

    ax11_2d = fig.add_subplot(221)  # just for plotting annotation (a)
    ax11_2d.set_axis_off()

    ms = 3.5
    labelpad = 0.1
    # plot lp vs f
    folder = "../data/20240821"
    L = 200

    # plot config for f
    for f in [0.04, 0.12, 0.20, 0.28]:
        ax_plot_config(ax11, "../data/20240730", [100, 10.0, f, 0.00], -10, fr"${f:.2f}$")
        # ax_plot_config(ax11, folder, [200, 10.0, f, 0.00], -10, fr"${f:.2f}$")
    # for f in [0.08, 0.14, 0.28]:
        # ax_plot_config(ax11, folder, [200, 10.0, f, 0.00], -10, fr"${f:.2f}$")
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
        lps, lp_thetas = get_lp_data_g(folder, param)
        if i == 0:
            ax12.plot(fs, lps, color="tomato", ls="none", marker=marker, ms=ms, mfc="None", label=r"$l_p$")
            ax12.plot(fs, lp_thetas, color="royalblue", ls="none", marker=marker, ms=ms, mfc="None", label=r"$l_{p,\theta}$")
        else:
            ax12.plot(fs, lps, color="tomato", ls="none", marker=marker, ms=ms, mfc="None")
            ax12.plot(fs, lp_thetas, color="royalblue", ls="none", marker=marker, ms=ms, mfc="None")
        # ax12.text(0, kappa-1.7, fr"$\kappa={kappa:.0f}$", fontsize=9)

    ax12.legend(title=rf"$\kappa={kappa:.0f}$", ncol=2, columnspacing=0.5, handlelength=0.5, handletextpad=0.2, frameon=False, fontsize=9)
    ax12.set_ylabel(r"$l_p,l_{p,\theta}$", fontsize=9, labelpad=labelpad)
    ax12.set_xlabel(r"$f$", fontsize=9, labelpad=labelpad)
    ax12.tick_params(which="both", direction="in", top="on", right="on", labelbottom=True, labelleft=True, labelsize=7)

    ax12.xaxis.set_major_locator(plt.MultipleLocator(0.2))
    ax12.xaxis.set_minor_locator(plt.MultipleLocator(0.1))
    ax12.yaxis.set_major_locator(plt.MultipleLocator(1))
    ax12.yaxis.set_minor_locator(plt.MultipleLocator(0.5))
    # ax12.set_ylim(2, 13)

    # plot R2 and Rg vs. f
    kappas = [5.0, 10.0]
    for i in range(len(kappas)):
        kappa = kappas[i]
        ls = lss[i]
        marker = markers[i]
        ls = "None"
        fs = np.arange(0.00, 0.301, 0.02)
        param = [(L, kappa, f, 0.0) for f in fs]
        X, Y, Z, XsignZ, ZsignX, R2, Rg2, Rxx, Rxz = get_obs_data_g(folder, param)
        X_errs, Y_errs, Z_errs, XsignZ_errs, ZsignX_errs, R2_errs, Rg2_errs, Rxx_errs, Rxz_errs = get_obs_data_err_g(folder, param)
        # ax21.errorbar(fs, R2/L, yerr=R2_err/L, ls=ls, marker=marker, ms=ms, mfc="None", label=fr"${kappa:.0f}$")
        # if (kappa == 5):
        # ax21.plot((X/L+0.25/(1-X/L)**2 - 0.25)/kappa, X/L, "-", color="gray") #, label="theory")
        ax21.plot((X/L+0.25/(1-X/L)**2 - 0.25)/kappa, X/L, "-", color="gray")

        ax21.plot(fs, X/L, ls=ls, marker=marker, ms=ms, mfc="None", label=fr"${kappa:.0f}$")
        # ax21.errorbar(fs, X/L, yerr=X_errs/L,ls=ls, marker=marker, ms=ms, mfc="None", label=fr"${kappa:.0f}$")

        if kappa == 10:
            ax22.plot(fs, Rxx/Rg2, ls=ls, marker=marker, ms=ms, mfc="None", label=r"$xx$")
            ax22.plot(fs, Rxz/Rg2, ls=ls, marker=marker, ms=ms, mfc="None", label=r"$xz$")
            # ax22.errorbar(fs, Rxx/Rg2, Rxx_errs/Rg2, ls=ls, marker=marker, ms=ms, mfc="None", label=r"$xx$")
            # ax22.errorbar(fs, Rxz/Rg2, Rxz_errs/Rg2, ls=ls, marker=marker, ms=ms, mfc="None", label=r"$xz$")

    ax21.legend(title=r"$\kappa$", loc="center right", ncol=1, columnspacing=0.5, handlelength=0.5, handletextpad=0.2, frameon=False, fontsize=9)
    ax21.set_ylabel(r"$X/L$", fontsize=9, labelpad=labelpad)
    ax21.set_xlabel(r"$f$", fontsize=9, labelpad=labelpad)
    ax21.tick_params(which="both", direction="in", top="on", right="on", labelbottom=True, labelleft=True, labelsize=7)
    # ax21.xaxis.set_major_locator(plt.MultipleLocator(0.2))
    # ax21.xaxis.set_minor_locator(plt.MultipleLocator(0.1))
    ax21.yaxis.set_major_locator(plt.MultipleLocator(0.2))
    ax21.yaxis.set_minor_locator(plt.MultipleLocator(0.1))

    # title=r"$\kappa=10$",
    ax22.legend(title=r"$\mu$", ncol=1, columnspacing=0.5, handlelength=0.5, handletextpad=0.2, frameon=False, fontsize=9)
    ax22.set_ylabel(r"$R_{\mu}/R_g^2$", fontsize=9, labelpad=labelpad)
    ax22.set_xlabel(r"$f$", fontsize=9, labelpad=labelpad)
    ax22.tick_params(which="both", direction="in", top="on", right="on", labelbottom=True, labelleft=True, labelsize=7)

    ax22.xaxis.set_major_locator(plt.MultipleLocator(0.1))
    ax22.xaxis.set_minor_locator(plt.MultipleLocator(0.05))
    ax22.yaxis.set_major_locator(plt.MultipleLocator(0.2))
    ax22.yaxis.set_minor_locator(plt.MultipleLocator(0.1))
    # ax14.set_ylim(-0.25,6)

    annotation = [r"$(a)$", r"$(b)$", r"$(c)$", r"$(d)$"]

    axs = [ax11_2d, ax12, ax21, ax22]
    for i in range(len(axs)):
        ax = axs[i]
        ax.text(0.82, 0.125, annotation[i], fontsize=9, transform=ax.transAxes)

    plt.tight_layout(pad=0.1)
    '''
    plt.subplots_adjust(left=0.08,
                    bottom=0.08,
                    right=0.99,
                    top=0.99,
                    wspace=0.3,
                    hspace=0.3)
    '''
    plt.subplots_adjust(wspace=0.3, hspace=0.25)
    plt.savefig("./figures/obs_f.pdf", format="pdf")
    plt.savefig("./figures/obs_f.png", format="png", dpi=300)
    # plt.subplot_tool()
    plt.show()
    plt.close()


def plot_obs_gamma(tex_lw=240.71031, ppi=72):

    fig = plt.figure(figsize=(tex_lw / ppi * 1, tex_lw / ppi * 0.9))
    plt.rc("text", usetex=True)
    plt.rc("text.latex", preamble=r"\usepackage{physics}")

    ax11 = fig.add_subplot(221, projection="3d")
    ax12 = fig.add_subplot(222)
    ax21 = fig.add_subplot(223, sharex=ax12)
    ax22 = fig.add_subplot(224, sharex=ax12)

    ax11_2d = fig.add_subplot(221)
    ax11_2d.set_axis_off()

    ms = 3.5
    labelpad = 0.1
    # plot lp vs f
    folder = "../data/20240820"
    L = 200

    # plot config vs g
    for gL in [0.40, 0.80, 1.20]:
        ax_plot_config(ax11, "../data/20240730", [100, 10.0, 0.00, gL], -30, fr"${gL:.1f}$")
    # ax11.view_init(elev=32., azim=-75)
    ax11.quiver(40, 5, 5, 20, 0, 0, color="black", arrow_length_ratio=0.4)
    ax11.text(60, 5, 5, r"$\vu{x}$", fontsize=9)
    ax11.quiver(40, 5, 5, 0, 0, 20, color="black", arrow_length_ratio=0.4)
    ax11.text(40, 5, 25, r"$\vu{z}$", fontsize=9)

    ax11.legend(title=r"$\gamma L$", loc="lower center", ncol=2, columnspacing=0.5, handlelength=0.5, handletextpad=0.2, frameon=False, fontsize=9)
    ax11.set_axis_off()

    # plot lp vs g
    L = 200
    kappa = 10
    fs = [0.00]
    gLs = np.arange(0.00, 1.501, 0.10)
    lss = ['-', "--"]
    markers = ["s", "o", "v", "d"]
    color = ["tomato", "royalblue"]
    for i in range(len(fs)):
        f = fs[i]
        ls = lss[i]
        ls = "none"
        marker = markers[i]
        param = [(L, kappa, f, gL) for gL in gLs]
        lps, lp_thetas = get_lp_data_g(folder, param)
        if i == 0:
            ax12.plot(gLs, lps, color="tomato", ls=ls, marker=marker, ms=ms, mfc="None", label=r"$l_p$")
            ax12.plot(gLs, lp_thetas, color="royalblue", ls=ls, marker=marker, ms=ms, mfc="None", label=r"$l_{p,\theta}$")
        else:
            ax12.plot(gLs, lps, color="tomato", marker=marker, ls=ls, ms=ms, mfc="None")
            ax12.plot(gLs, lp_thetas, color="royalblue", marker=marker, ls=ls, ms=ms, mfc="None")
        # ax12.text(0, kappa-1.5, fr"$f={f}$", fontsize=9)

    ax12.legend(title=rf"$\kappa={kappa:.0f}$", ncol=2, columnspacing=0.5, handlelength=0.5, handletextpad=0.2, frameon=False, fontsize=9)
    ax12.set_ylabel(r"$l_p,l_{p,\theta}$", fontsize=9, labelpad=labelpad)
    ax12.set_xlabel(r"$\gamma L$", fontsize=9, labelpad=labelpad)
    ax12.tick_params(which="both", direction="in", top="on", right="on", labelbottom=True, labelleft=True, labelsize=7)

    # ax12.xaxis.set_major_locator(plt.MultipleLocator(0.2))
    # ax12.xaxis.set_minor_locator(plt.MultipleLocator(0.1))
    # ax12.yaxis.set_major_locator(plt.MultipleLocator(2))
    # ax12.yaxis.set_minor_locator(plt.MultipleLocator(1))
    # ax12.set_ylim(2, 13)

    # plot X/L and Rg vs. g  , f is f_t, g is fs = =
    fs = [0.00]
    f = 0.0
    kappas = [5.0, 10.0]
    folder = "../data/20240820"
    for i in range(len(kappas)):
        # f = fs[i]
        kappa = kappas[i]
        # ls = lss[i]
        ls = "None"
        marker = markers[i]
        param = [(L, kappa, f, gL) for gL in gLs]
        # X, Y, Z, R2, Rg2, Rxx, Rxz = get_obs_data(folder, param)
        X, Y, Z, XsignZ, ZsignX, R2, Rg2, Rxx, Rxz = get_obs_data_g(folder, param)
        # ax23.errorbar(gLs, R2/L, yerr=R2_err/L, ls=ls, marker=marker, ms=ms, mfc="None", label=fr"${f:.1f}$")
        # ax21.plot(gLs, R2/L, ls=ls, marker=marker, ms=ms, mfc="None", label=fr"${f:.1f}$")
        ax21.plot(gLs, XsignZ/L, ls=ls, marker=marker, ms=ms+1, mfc="None", label=rf"${kappa:.0f}$")
        # ax21.plot(gLs, ZsignX/L, ls=ls, marker="2", ms=ms+2, mfc="None", label=r"$Z\theta(X)$")
        # ax21.plot(gLs, np.abs(Y)/L, ls=ls, marker=marker, ms=ms, mfc="None", label=fr"Y${f:.1f}$")
        # ax21.plot(gLs, np.abs(Z)/L, ls=ls, marker=marker, ms=ms, mfc="None", label=fr"Z${f:.1f}$")
        # ax22.errorbar(gLs, Rg, yerr=Rg_err, ls=ls,  marker=marker, ms=ms, mfc="None", label=fr"${f:.1f}$")
        if kappa == 10:
            # ax22.plot(gLs, Rg2/L, ls=ls, marker=marker, ms=ms, mfc="None", label=r"$R_g^2$")
            ax22.plot(gLs, Rxx/Rg2, ls=ls, marker=marker, ms=ms, mfc="None", label=r"$xx$")
            ax22.plot(gLs, Rxz/Rg2, ls=ls, marker=marker, ms=ms, mfc="None", label=r"$xz$")

    ax21.legend(title=r"$\kappa$", ncol=1, columnspacing=0.5, handlelength=0.5, handletextpad=0.2, frameon=False, fontsize=9)
    ax21.set_ylabel(r"$X\frac{Z}{|Z|}/L$", fontsize=9, labelpad=labelpad)
    ax21.set_xlabel(r"$\gamma L$", fontsize=9, labelpad=labelpad)
    ax21.tick_params(which="both", direction="in", top="on", right="on", labelbottom=True, labelleft=True, labelsize=7)
    ax21.xaxis.set_major_locator(plt.MultipleLocator(0.2))
    ax21.xaxis.set_minor_locator(plt.MultipleLocator(0.1))
    ax21.yaxis.set_major_locator(plt.MultipleLocator(0.2))
    ax21.yaxis.set_minor_locator(plt.MultipleLocator(0.1))
    # title=r"$f=0$"
    ax22.legend(title=r"$\mu$", ncol=1, columnspacing=0.5, handlelength=0.5, handletextpad=0.2, frameon=False, fontsize=9)
    ax22.set_ylabel(r"$R_{\mu}/R_g^2$", fontsize=9, labelpad=labelpad)
    ax22.set_xlabel(r"$\gamma L$", fontsize=9, labelpad=labelpad)
    ax22.tick_params(which="both", direction="in", top="on", right="on", labelbottom=True, labelleft=True, labelsize=7)
    # ax22.set_ylim(-0.25,7)
    ax22.xaxis.set_major_locator(plt.MultipleLocator(0.5))
    ax22.xaxis.set_minor_locator(plt.MultipleLocator(0.25))
    ax22.yaxis.set_major_locator(plt.MultipleLocator(0.2))
    ax22.yaxis.set_minor_locator(plt.MultipleLocator(0.1))

    annotation = [r"$(a)$", r"$(b)$", r"$(c)$", r"$(d)$"]

    axs = [ax11_2d, ax12, ax21, ax22]
    for i in range(len(axs)):
        ax = axs[i]
        ax.text(0.82, 0.125, annotation[i], fontsize=9, transform=ax.transAxes)

    plt.tight_layout(pad=0.1)
    plt.subplots_adjust(wspace=0.3, hspace=0.25)
    plt.savefig("./figures/obs_gamma.pdf", format="pdf")
    plt.savefig("./figures/obs_gamma.png", format="png", dpi=300)
    plt.show()
    plt.close()


def plot_obs_f_3panel(tex_lw=240.71031, ppi=72):
    fig = plt.figure(figsize=(tex_lw / ppi * 1, tex_lw / ppi * 0.75))
    plt.rc("text", usetex=True)
    plt.rc("text.latex", preamble=r"\usepackage{physics}")

    ax11 = plt.subplot2grid((5, 2), (0, 0), colspan=2, rowspan=2)
    ax21 = plt.subplot2grid((5, 2), (2, 0), rowspan=3)
    ax22 = plt.subplot2grid((5, 2), (2, 1), rowspan=3)

    ms = 3.5
    labelpad = 0.1
    # plot lp vs f
    folder = "../data/20240828"
    L = 200
    kappa = 10
    gL = 0.0
    fs = np.arange(0.00, 0.301, 0.02)
    fi2plot = [0, 5, 10, 15]
    xshifts = [0, 100, 150, 175]
    xshifts = np.cumsum(xshifts)
    xtextshifts = [0, 150, 300, 475]
    for i in range(len(fi2plot)):
        f = fs[fi2plot[i]]
        finfo = f"L{L}_kappa{kappa:.1f}_f{f:.2f}_gL{gL:.2f}"
        colors = plt.cm.rainbow(np.linspace(0, 1, 30))
        for j in range(30):
            filename = f"{folder}/{finfo}/config_{j}.csv"
            ax_plot_2dconfig_from_file(ax11, filename, colors[j], 120, 0, xshifts[i], 50)
        ax11.text(xtextshifts[i]-50, -50, fr"$f={f:.1f}$", fontsize=9)

    alen=40
    x0,y0=-100,-20
    ax11.arrow(x0, y0, alen, 0, color=(0,0,1), lw=1, head_length=alen*0.3, head_width=alen*0.3)  # , overhang=1.0) #, length_includes_head=True)
    ax11.text(x0+alen, y0+0.1*alen, r"$\vu{x}$", fontsize=9)
    ax11.arrow(x0, y0, 0, alen, color=(1,0,0), lw=1, head_length=alen*0.3, head_width=alen*0.3)  # , overhang=1.0) #, length_includes_head=True)
    ax11.text(x0+0.1*alen, y0+alen, r"$\vu{z}$", fontsize=9)

    ax11.set_axis_off()

    # plot R2 and Rg vs. f
    kappas = [5.0, 10.0]
    lss = ['-', "--"]
    markers = ["s", "o", "v", "d"]
    for i in range(len(kappas)):
        kappa = kappas[i]
        ls = lss[i]
        marker = markers[i]
        ls = "None"
        fs = np.arange(0.00, 0.301, 0.02)
        param = [(L, kappa, f, 0.0) for f in fs]
        X, Y, Z, XsignZ, ZsignX, R2, Rg2, Rxx, Rxz = get_obs_data_g(folder, param)
        X_errs, Y_errs, Z_errs, XsignZ_errs, ZsignX_errs, R2_errs, Rg2_errs, Rxx_errs, Rxz_errs = get_obs_data_err_g(folder, param)
        # ax21.errorbar(fs, R2/L, yerr=R2_err/L, ls=ls, marker=marker, ms=ms, mfc="None", label=fr"${kappa:.0f}$")
        # if (kappa == 5):
        # ax21.plot((X/L+0.25/(1-X/L)**2 - 0.25)/kappa, X/L, "-", color="gray") #, label="theory")
        ax21.plot((X/L+0.25/(1-X/L)**2 - 0.25)/kappa, X/L, "-", color="gray")

        ax21.plot(fs, X/L, ls=ls, marker=marker, ms=ms, mfc="None", label=fr"${kappa:.0f}$")
        # ax21.errorbar(fs, X/L, yerr=X_errs/L,ls=ls, marker=marker, ms=ms, mfc="None", label=fr"${kappa:.0f}$")

        if kappa == 10:
            ax22.plot(fs, Rxx/Rg2, ls=ls, marker=marker, ms=ms, mfc="None", label=r"$xx$")
            ax22.plot(fs, Rxz/Rg2, ls=ls, marker=marker, ms=ms, mfc="None", label=r"$xz$")
            # ax22.errorbar(fs, Rxx/Rg2, Rxx_errs/Rg2, ls=ls, marker=marker, ms=ms, mfc="None", label=r"$xx$")
            # ax22.errorbar(fs, Rxz/Rg2, Rxz_errs/Rg2, ls=ls, marker=marker, ms=ms, mfc="None", label=r"$xz$")

    ax21.legend(title=r"$\kappa$", loc="center right", ncol=1, columnspacing=0.5, handlelength=0.5, handletextpad=0.2, frameon=False, fontsize=9)
    ax21.set_ylabel(r"$X/L$", fontsize=9, labelpad=labelpad)
    ax21.set_xlabel(r"$f$", fontsize=9, labelpad=labelpad)
    ax21.tick_params(which="both", direction="in", top="on", right="on", labelbottom=True, labelleft=True, labelsize=7)
    ax21.xaxis.set_major_locator(plt.MultipleLocator(0.1))
    ax21.xaxis.set_minor_locator(plt.MultipleLocator(0.05))
    ax21.yaxis.set_major_locator(plt.MultipleLocator(0.2))
    ax21.yaxis.set_minor_locator(plt.MultipleLocator(0.1))

    # title=r"$\kappa=10$",
    ax22.legend(title=r"$\mu$", ncol=1, columnspacing=0.5, handlelength=0.5, handletextpad=0.2, frameon=False, fontsize=9)
    ax22.set_ylabel(r"$R_{\mu}/R_g^2$", fontsize=9, labelpad=labelpad)
    ax22.set_xlabel(r"$f$", fontsize=9, labelpad=labelpad)
    ax22.tick_params(which="both", direction="in", top="on", right="on", labelbottom=True, labelleft=True, labelsize=7)

    ax22.xaxis.set_major_locator(plt.MultipleLocator(0.1))
    ax22.xaxis.set_minor_locator(plt.MultipleLocator(0.05))
    ax22.yaxis.set_major_locator(plt.MultipleLocator(0.2))
    ax22.yaxis.set_minor_locator(plt.MultipleLocator(0.1))
    # ax14.set_ylim(-0.25,6)

    ax11.text(0.91, 0.0, r"$(a)$", fontsize=9, transform=ax11.transAxes)
    ax21.text(0.82, 0.125, r"$(b)$", fontsize=9, transform=ax21.transAxes)
    ax22.text(0.82, 0.125, r"$(c)$", fontsize=9, transform=ax22.transAxes)
    #ax11.set_rasterized(True)
    #ax21.set_rasterized(True)
    #ax22.set_rasterized(True)
    plt.tight_layout(pad=0.1)

    plt.subplots_adjust(wspace=0.3, hspace=0.25)
    plt.savefig("./figures/obs_f_3panel.pdf", format="pdf", dpi=300)
    plt.savefig("./figures/obs_f_3panel.png", format="png", dpi=300)
    plt.show()
    plt.close()


def plot_obs_gamma_3panel(tex_lw=240.71031, ppi=72):

    fig = plt.figure(figsize=(tex_lw / ppi * 1, tex_lw / ppi * 0.75))
    plt.rc("text", usetex=True)
    plt.rc("text.latex", preamble=r"\usepackage{physics}")

    ax11 = plt.subplot2grid((5, 2), (0, 0), colspan=2, rowspan=2)
    ax21 = plt.subplot2grid((5, 2), (2, 0), rowspan=3)
    ax22 = plt.subplot2grid((5, 2), (2, 1), rowspan=3)

    # ax11_2d = plt.subplot2grid((5, 2), (0, 0), colspan=2)
    # ax11_2d.set_axis_off()

    ms = 3.5
    labelpad = 0.1
    # plot lp vs f
    folder = "../data/20240820"
    L = 200
    gLs = np.arange(0.00, 1.501, 0.10)
    f = 0.0
    kappa = 10.0

    xlim = 120
    xshift = 220
    gLi2plot = [0, 3, 9, 15]
    for i in range(len(gLi2plot)):
        gL = gLs[gLi2plot[i]]
        finfo = f"L{L}_kappa{kappa:.1f}_f{f:.2f}_gL{gL:.2f}"
        colors = plt.cm.rainbow(np.linspace(0, 1, 30))
        for j in range(30):
            filename = f"{folder}/{finfo}/config_{j}.csv"
            ax_plot_2dconfig_from_file(ax11, filename, colors[j], xlim, 1, i*xshift, 50)
        ax11.text(i*xshift-100, -100, fr"$\gamma L={gL:.1f}$", fontsize=9)

    alen=40
    x0,y0=-100,-50
    ax11.arrow(x0, y0, alen, 0, color=(0,0,1), lw=1, head_length=alen*0.3, head_width=alen*0.3)  # , overhang=1.0) #, length_includes_head=True)
    ax11.text(x0+alen, y0+0.1*alen, r"$\vu{x}$", fontsize=9)
    ax11.arrow(x0, y0, 0, alen, color=(1,0,0), lw=1, head_length=alen*0.3, head_width=alen*0.3)  # , overhang=1.0) #, length_includes_head=True)
    ax11.text(x0+0.1*alen, y0+alen, r"$\vu{z}$", fontsize=9)

    ax11.set_axis_off()

    lss = ['-', "--"]
    markers = ["s", "o", "v", "d"]
    color = ["tomato", "royalblue"]
    # plot X/L and Rg vs. g  , f is f_t, g is fs = =
    f = 0.0
    kappas = [5.0, 10.0]
    folder = "../data/20240820"
    for i in range(len(kappas)):
        kappa = kappas[i]
        ls = "None"
        marker = markers[i]
        param = [(L, kappa, f, gL) for gL in gLs]
        X, Y, Z, XsignZ, ZsignX, R2, Rg2, Rxx, Rxz = get_obs_data_g(folder, param)
        ax21.plot(gLs, XsignZ/L, ls=ls, marker=marker, ms=ms+1, mfc="None", label=rf"${kappa:.0f}$")
        if kappa == 10:
            # ax22.plot(gLs, Rg2/L, ls=ls, marker=marker, ms=ms, mfc="None", label=r"$R_g^2$")
            ax22.plot(gLs, Rxx/Rg2, ls=ls, marker=marker, ms=ms, mfc="None", label=r"$xx$")
            ax22.plot(gLs, Rxz/Rg2, ls=ls, marker=marker, ms=ms, mfc="None", label=r"$xz$")

    ax21.legend(title=r"$\kappa$", ncol=1, columnspacing=0.5, handlelength=0.5, handletextpad=0.2, frameon=False, fontsize=9)
    ax21.set_ylabel(r"$X\frac{Z}{|Z|}/L$", fontsize=9, labelpad=labelpad)
    ax21.set_xlabel(r"$\gamma L$", fontsize=9, labelpad=labelpad)
    ax21.tick_params(which="both", direction="in", top="on", right="on", labelbottom=True, labelleft=True, labelsize=7)
    ax21.xaxis.set_major_locator(plt.MultipleLocator(0.5))
    ax21.xaxis.set_minor_locator(plt.MultipleLocator(0.25))
    ax21.yaxis.set_major_locator(plt.MultipleLocator(0.2))
    ax21.yaxis.set_minor_locator(plt.MultipleLocator(0.1))
    # title=r"$f=0$"
    ax22.legend(title=r"$\mu$", ncol=1, columnspacing=0.5, handlelength=0.5, handletextpad=0.2, frameon=False, fontsize=9)
    ax22.set_ylabel(r"$R_{\mu}/R_g^2$", fontsize=9, labelpad=labelpad)
    ax22.set_xlabel(r"$\gamma L$", fontsize=9, labelpad=labelpad)
    ax22.tick_params(which="both", direction="in", top="on", right="on", labelbottom=True, labelleft=True, labelsize=7)
    # ax22.set_ylim(-0.25,7)
    ax22.xaxis.set_major_locator(plt.MultipleLocator(0.5))
    ax22.xaxis.set_minor_locator(plt.MultipleLocator(0.25))
    ax22.yaxis.set_major_locator(plt.MultipleLocator(0.2))
    ax22.yaxis.set_minor_locator(plt.MultipleLocator(0.1))

    ax11.text(0.91, 0.125, r"$(a)$", fontsize=9, transform=ax11.transAxes)
    ax21.text(0.82, 0.125, r"$(b)$", fontsize=9, transform=ax21.transAxes)
    ax22.text(0.82, 0.125, r"$(c)$", fontsize=9, transform=ax22.transAxes)
    #ax11.set_rasterized(True)
    #ax21.set_rasterized(True)
    #ax22.set_rasterized(True)

    plt.tight_layout(pad=0.1)
    plt.subplots_adjust(wspace=0.3, hspace=0.25)
    plt.savefig("./figures/obs_gamma_3panel.pdf", format="pdf", dpi=300)
    plt.savefig("./figures/obs_gamma_3panel.png", format="png", dpi=300)
    plt.show()
    plt.close()


def plot_obs_f_3panel_sq(tex_lw=240.71031, ppi=72):
    fig = plt.figure(figsize=(tex_lw / ppi * 1, tex_lw / ppi * 1.2))
    plt.rc("text", usetex=True)
    plt.rc("text.latex", preamble=r"\usepackage{physics}")

    #axconf = plt.subplot2grid((5, 2), (0, 0), colspan=2, rowspan=2)
    #axX = plt.subplot2grid((5, 2), (2, 0), rowspan=3)
    #axRg = plt.subplot2grid((5, 2), (2, 1), rowspan=3)

    axconf = plt.subplot2grid((7, 4), (0, 0), colspan=4, rowspan=2)
    axsq0 = plt.subplot2grid((7, 4), (2, 0), rowspan=2)
    axsq1 = plt.subplot2grid((7, 4), (2, 1), rowspan=2)
    axsq2 = plt.subplot2grid((7, 4), (2, 2), rowspan=2)
    axsq3 = plt.subplot2grid((7, 4), (2, 3), rowspan=2)
    axX = plt.subplot2grid((7, 4), (4, 0), colspan=2, rowspan=3)
    axRg = plt.subplot2grid((7, 4), (4, 2), colspan=2, rowspan=3)

    ms = 3.5
    labelpad = 0.1
    # plot lp vs f
    folder = "../data/20240828"
    L = 200
    kappa = 10
    gL = 0.0
    fs = np.arange(0.00, 0.301, 0.02)
    fi2plot = [0, 5, 10, 15]
    xshifts = [0, 100, 150, 175]
    xshifts = np.cumsum(xshifts)
    xtextshifts = [0, 150, 300, 475]
    for i in range(len(fi2plot)):
        f = fs[fi2plot[i]]
        finfo = f"L{L}_kappa{kappa:.1f}_f{f:.2f}_gL{gL:.2f}"
        colors = plt.cm.rainbow(np.linspace(0, 1, 30))
        for j in range(30):
            filename = f"{folder}/{finfo}/config_{j}.csv"
            ax_plot_2dconfig_from_file(axconf, filename, colors[j], 120, 0, xshifts[i], 50)
        axconf.text(xtextshifts[i]-50, -50, fr"$f={f:.1f}$", fontsize=9)

    alen=40
    x0,y0=-100,-20
    axconf.arrow(x0, y0, alen, 0, color=(0,0,1), lw=1, head_length=alen*0.3, head_width=alen*0.3)  # , overhang=1.0) #, length_includes_head=True)
    axconf.text(x0+alen, y0+0.1*alen, r"$\vu{x}$", fontsize=9)
    axconf.arrow(x0, y0, 0, alen, color=(1,0,0), lw=1, head_length=alen*0.3, head_width=alen*0.3)  # , overhang=1.0) #, length_includes_head=True)
    axconf.text(x0+0.1*alen, y0+alen, r"$\vu{z}$", fontsize=9)

    axconf.set_axis_off()


    kappa = 10.0
    f = [0.00, 0.10, 0.20, 0.30]
    # plot Sq2D for various kappa
    axs = [axsq0, axsq1, axsq2, axsq3]
    for i in range(len(f)):
        filename = f"../data/20240920_precision/obs_L{L:.0f}_kappa{kappa:.0f}_f{f[i]:.2f}_gL{gL:.2f}.csv"
        Sq2D, qB = get_Sq2D_data(filename)
        qBx, qBz = np.meshgrid(qB, qB)
        print(np.min(Sq2D), np.max(Sq2D))
        axs[i].pcolormesh(qBx, qBz, np.log10(Sq2D), vmax=0, vmin=-3, cmap="rainbow", shading='gouraud')
        Cs = axs[i].contour(qBx, qBz, np.log10(Sq2D), vmax=0, vmin=-3, levels=np.linspace(-3, 0, 7), colors="gray", linewidths=0.5, linestyle=":")
        axs[i].clabel(Cs, Cs.levels, inline=True, fontsize=7, fmt="%1.1f", colors="black")
        axs[i].set_xlabel(r"$Q_x$", fontsize=9, labelpad=-0.0)
        axs[i].tick_params(which="both", direction="in", top="on", right="on", labelbottom=True, labelleft=False, labelsize=7)
        axs[i].xaxis.set_major_locator(plt.MultipleLocator(0.5))
        axs[i].xaxis.set_minor_locator(plt.MultipleLocator(0.25))
        axs[i].yaxis.set_major_locator(plt.MultipleLocator(0.5))
        axs[i].yaxis.set_minor_locator(plt.MultipleLocator(0.25))
        axs[i].set_title(r"$f = $"+f"{f[i]:.1f}", fontsize=9, pad=2.5)
        axs[i].set_aspect("equal")

    axs[0].set_ylabel(r"$Q_z$", fontsize=9, labelpad=-0.0)
    axs[0].tick_params(which="both", direction="in", top="on", right="on", labelbottom=True, labelleft=True, labelsize=7)



    folder = "../data/20240828"

    # plot R2 and Rg vs. f
    kappas = [5.0, 10.0]
    lss = ['-', "--"]
    markers = ["s", "o", "v", "d"]
    for i in range(len(kappas)):
        kappa = kappas[i]
        ls = lss[i]
        marker = markers[i]
        ls = "None"
        fs = np.arange(0.00, 0.301, 0.02)
        param = [(L, kappa, f, 0.0) for f in fs]
        X, Y, Z, XsignZ, ZsignX, R2, Rg2, Rxx, Rxz = get_obs_data_g(folder, param)
        X_errs, Y_errs, Z_errs, XsignZ_errs, ZsignX_errs, R2_errs, Rg2_errs, Rxx_errs, Rxz_errs = get_obs_data_err_g(folder, param)
        # axX.errorbar(fs, R2/L, yerr=R2_err/L, ls=ls, marker=marker, ms=ms, mfc="None", label=fr"${kappa:.0f}$")
        # if (kappa == 5):
        # axX.plot((X/L+0.25/(1-X/L)**2 - 0.25)/kappa, X/L, "-", color="gray") #, label="theory")
        axX.plot((X/L+0.25/(1-X/L)**2 - 0.25)/kappa, X/L, "-", color="gray")

        axX.plot(fs, X/L, ls=ls, marker=marker, ms=ms, mfc="None", label=fr"${kappa:.0f}$")
        # axX.errorbar(fs, X/L, yerr=X_errs/L,ls=ls, marker=marker, ms=ms, mfc="None", label=fr"${kappa:.0f}$")

        if kappa == 10:
            axRg.plot(fs, Rxx/Rg2, ls=ls, marker=marker, ms=ms, mfc="None", label=r"$xx$")
            axRg.plot(fs, Rxz/Rg2, ls=ls, marker=marker, ms=ms, mfc="None", label=r"$xz$")
            # axRg.errorbar(fs, Rxx/Rg2, Rxx_errs/Rg2, ls=ls, marker=marker, ms=ms, mfc="None", label=r"$xx$")
            # axRg.errorbar(fs, Rxz/Rg2, Rxz_errs/Rg2, ls=ls, marker=marker, ms=ms, mfc="None", label=r"$xz$")

    axX.legend(title=r"$\kappa$", loc="center right", ncol=1, columnspacing=0.5, handlelength=0.5, handletextpad=0.2, frameon=False, fontsize=9)
    axX.set_ylabel(r"$X/L$", fontsize=9, labelpad=labelpad)
    axX.set_xlabel(r"$f$", fontsize=9, labelpad=labelpad)
    axX.tick_params(which="both", direction="in", top="on", right="on", labelbottom=True, labelleft=True, labelsize=7)
    axX.xaxis.set_major_locator(plt.MultipleLocator(0.1))
    axX.xaxis.set_minor_locator(plt.MultipleLocator(0.05))
    axX.yaxis.set_major_locator(plt.MultipleLocator(0.2))
    axX.yaxis.set_minor_locator(plt.MultipleLocator(0.1))

    # title=r"$\kappa=10$",
    axRg.legend(title=r"$\mu$", ncol=1, columnspacing=0.5, handlelength=0.5, handletextpad=0.2, frameon=False, fontsize=9)
    axRg.set_ylabel(r"$R_{\mu}/R_g^2$", fontsize=9, labelpad=labelpad)
    axRg.set_xlabel(r"$f$", fontsize=9, labelpad=labelpad)
    axRg.tick_params(which="both", direction="in", top="on", right="on", labelbottom=True, labelleft=True, labelsize=7)

    axRg.xaxis.set_major_locator(plt.MultipleLocator(0.1))
    axRg.xaxis.set_minor_locator(plt.MultipleLocator(0.05))
    axRg.yaxis.set_major_locator(plt.MultipleLocator(0.2))
    axRg.yaxis.set_minor_locator(plt.MultipleLocator(0.1))
    # ax14.set_ylim(-0.25,6)



    axconf.text(0.91, 0.0, r"$(a)$", fontsize=9, transform=axconf.transAxes)
    axX.text(0.82, 0.125, r"$(b)$", fontsize=9, transform=axX.transAxes)
    axRg.text(0.82, 0.125, r"$(c)$", fontsize=9, transform=axRg.transAxes)
    #ax11.set_rasterized(True)
    #axX.set_rasterized(True)
    #axRg.set_rasterized(True)
    plt.tight_layout(pad=0.1)

    plt.subplots_adjust(wspace=0.3, hspace=0.25)
    plt.savefig("./figures/obs_f_3panel.pdf", format="pdf", dpi=300)
    plt.savefig("./figures/obs_f_3panel.png", format="png", dpi=300)
    plt.show()
    plt.close()


def plot_obs_gamma_3panel_sq(tex_lw=240.71031, ppi=72):

    fig = plt.figure(figsize=(tex_lw / ppi * 1, tex_lw / ppi * 0.75))
    plt.rc("text", usetex=True)
    plt.rc("text.latex", preamble=r"\usepackage{physics}")

    ax11 = plt.subplot2grid((5, 2), (0, 0), colspan=2, rowspan=2)
    ax21 = plt.subplot2grid((5, 2), (2, 0), rowspan=3)
    ax22 = plt.subplot2grid((5, 2), (2, 1), rowspan=3)

    # ax11_2d = plt.subplot2grid((5, 2), (0, 0), colspan=2)
    # ax11_2d.set_axis_off()

    ms = 3.5
    labelpad = 0.1
    # plot lp vs f
    folder = "../data/20240820"
    L = 200
    gLs = np.arange(0.00, 1.501, 0.10)
    f = 0.0
    kappa = 10.0

    xlim = 120
    xshift = 220
    gLi2plot = [0, 3, 9, 15]
    for i in range(len(gLi2plot)):
        gL = gLs[gLi2plot[i]]
        finfo = f"L{L}_kappa{kappa:.1f}_f{f:.2f}_gL{gL:.2f}"
        colors = plt.cm.rainbow(np.linspace(0, 1, 30))
        for j in range(30):
            filename = f"{folder}/{finfo}/config_{j}.csv"
            ax_plot_2dconfig_from_file(ax11, filename, colors[j], xlim, 1, i*xshift, 50)
        ax11.text(i*xshift-100, -100, fr"$\gamma L={gL:.1f}$", fontsize=9)

    alen=40
    x0,y0=-100,-50
    ax11.arrow(x0, y0, alen, 0, color=(0,0,1), lw=1, head_length=alen*0.3, head_width=alen*0.3)  # , overhang=1.0) #, length_includes_head=True)
    ax11.text(x0+alen, y0+0.1*alen, r"$\vu{x}$", fontsize=9)
    ax11.arrow(x0, y0, 0, alen, color=(1,0,0), lw=1, head_length=alen*0.3, head_width=alen*0.3)  # , overhang=1.0) #, length_includes_head=True)
    ax11.text(x0+0.1*alen, y0+alen, r"$\vu{z}$", fontsize=9)

    ax11.set_axis_off()

    lss = ['-', "--"]
    markers = ["s", "o", "v", "d"]
    color = ["tomato", "royalblue"]
    # plot X/L and Rg vs. g  , f is f_t, g is fs = =
    f = 0.0
    kappas = [5.0, 10.0]
    folder = "../data/20240820"
    for i in range(len(kappas)):
        kappa = kappas[i]
        ls = "None"
        marker = markers[i]
        param = [(L, kappa, f, gL) for gL in gLs]
        X, Y, Z, XsignZ, ZsignX, R2, Rg2, Rxx, Rxz = get_obs_data_g(folder, param)
        ax21.plot(gLs, XsignZ/L, ls=ls, marker=marker, ms=ms+1, mfc="None", label=rf"${kappa:.0f}$")
        if kappa == 10:
            # ax22.plot(gLs, Rg2/L, ls=ls, marker=marker, ms=ms, mfc="None", label=r"$R_g^2$")
            ax22.plot(gLs, Rxx/Rg2, ls=ls, marker=marker, ms=ms, mfc="None", label=r"$xx$")
            ax22.plot(gLs, Rxz/Rg2, ls=ls, marker=marker, ms=ms, mfc="None", label=r"$xz$")

    ax21.legend(title=r"$\kappa$", ncol=1, columnspacing=0.5, handlelength=0.5, handletextpad=0.2, frameon=False, fontsize=9)
    ax21.set_ylabel(r"$X\frac{Z}{|Z|}/L$", fontsize=9, labelpad=labelpad)
    ax21.set_xlabel(r"$\gamma L$", fontsize=9, labelpad=labelpad)
    ax21.tick_params(which="both", direction="in", top="on", right="on", labelbottom=True, labelleft=True, labelsize=7)
    ax21.xaxis.set_major_locator(plt.MultipleLocator(0.5))
    ax21.xaxis.set_minor_locator(plt.MultipleLocator(0.25))
    ax21.yaxis.set_major_locator(plt.MultipleLocator(0.2))
    ax21.yaxis.set_minor_locator(plt.MultipleLocator(0.1))
    # title=r"$f=0$"
    ax22.legend(title=r"$\mu$", ncol=1, columnspacing=0.5, handlelength=0.5, handletextpad=0.2, frameon=False, fontsize=9)
    ax22.set_ylabel(r"$R_{\mu}/R_g^2$", fontsize=9, labelpad=labelpad)
    ax22.set_xlabel(r"$\gamma L$", fontsize=9, labelpad=labelpad)
    ax22.tick_params(which="both", direction="in", top="on", right="on", labelbottom=True, labelleft=True, labelsize=7)
    # ax22.set_ylim(-0.25,7)
    ax22.xaxis.set_major_locator(plt.MultipleLocator(0.5))
    ax22.xaxis.set_minor_locator(plt.MultipleLocator(0.25))
    ax22.yaxis.set_major_locator(plt.MultipleLocator(0.2))
    ax22.yaxis.set_minor_locator(plt.MultipleLocator(0.1))

    ax11.text(0.91, 0.125, r"$(a)$", fontsize=9, transform=ax11.transAxes)
    ax21.text(0.82, 0.125, r"$(b)$", fontsize=9, transform=ax21.transAxes)
    ax22.text(0.82, 0.125, r"$(c)$", fontsize=9, transform=ax22.transAxes)
    #ax11.set_rasterized(True)
    #ax21.set_rasterized(True)
    #ax22.set_rasterized(True)

    plt.tight_layout(pad=0.1)
    plt.subplots_adjust(wspace=0.3, hspace=0.25)
    plt.savefig("./figures/obs_gamma_3panel.pdf", format="pdf", dpi=300)
    plt.savefig("./figures/obs_gamma_3panel.png", format="png", dpi=300)
    plt.show()
    plt.close()





def plot_TOC(tex_lw=240.71031, ppi=72):
    cm = 1/2.54
    # fig = plt.figure(figsize=(9*cm, 3.5*cm))
    fig = plt.figure(figsize=(8.25*cm, 4.45*cm))
    plt.rc("text", usetex=True)
    plt.rc("text.latex", preamble=r"\usepackage{physics}")

    ax = fig.add_subplot(111)
    ms = 3.5
    labelpad = 0.1
    # plot lp vs f
    folder = "../data/20240820"
    L = 200
    f = 0.0
    kappa = 10.0
    gLs = np.arange(0.00, 1.501, 0.10)
    folder = "../data/20240820"
    param = [(L, kappa, f, gL) for gL in gLs]
    X, Y, Z, XsignZ, ZsignX, R2, Rg2, Rxx, Rxz = get_obs_data_g(folder, param)
    ax.plot(gLs, XsignZ/L, "o-", lw=1, mfc="None", zorder=10)
    ax.set_ylabel(r"$X\frac{Z}{|Z|}/L$", fontsize=10, labelpad=labelpad)
    ax.set_xlabel(r"$\gamma L$ (Shear)", fontsize=10, labelpad=labelpad)
    ax.tick_params(which="both", direction="in", top="on", right="on", labelbottom=True, labelleft=True, labelsize=8)
    ax.xaxis.set_major_locator(plt.MultipleLocator(0.2))
    ax.xaxis.set_minor_locator(plt.MultipleLocator(0.1))
    ax.yaxis.set_major_locator(plt.MultipleLocator(0.2))
    ax.yaxis.set_minor_locator(plt.MultipleLocator(0.1))

    # plot config
    # g=0
    xlim = 130
    gL = 0.
    inset_size = 0.5
    inset_axs = []  # [ax1, ax2, ax3]
    inset_axs_pos = [[0.0, 0.37], [0.3, 0.35], [0.6, 0.35]]
    inset_axs_abs_pos = [[0.1, 0.2], [0.6, 0.22], [1.3, 0.3]]
    gLi2plot = [0, 4, 12]
    for i in range(len(gLi2plot)):
        gL = gLs[gLi2plot[i]]
        XsignZp = XsignZ[gLi2plot[i]]
        inset_axs.append(fig.add_axes([inset_axs_pos[i][0], inset_axs_pos[i][1], inset_size, inset_size]))
        ax.annotate(" ", xytext=inset_axs_abs_pos[i], xy=[gL, XsignZp/L], arrowprops=dict(arrowstyle="->", lw=1, color="gray"), annotation_clip=False)
        finfo = f"L{L}_kappa{kappa:.1f}_f{f:.2f}_gL{gL:.2f}"
        colors = plt.cm.rainbow(np.linspace(0, 1, 30))
        for j in range(30):
            filename = f"{folder}/{finfo}/config_{j}.csv"
            ax_plot_2dconfig_from_file(inset_axs[i], filename, colors[j], xlim)
        inset_axs[i].set_axis_off()
        inset_axs[i].set_xlim(-xlim, xlim)
        inset_axs[i].set_ylim(-xlim, xlim)
        # inset_ax1 = ax.inset_axes([-0.2, 0.25, inset_size, inset_size], projection='3d')
        # inset_ax1 = ax.inset_axes([0.1, 0.25, inset_size, inset_size])
    # inset_ax1 = fig.add_axes([.0, .35, inset_size, inset_size])
    # inset_ax1.text(0.1, 0.9, rf"$f_sL={gL:.1f}$", fontsize=9, transform=inset_ax1.transAxes)
    # ax.annotate(" ", xytext=(0.0, 0.35), xy=[gLs[0], XsignZ[0]/L], arrowprops=dict(arrowstyle="->", lw=1, color="gray"), annotation_clip=False)

    # inset_ax1.view_init(elev=32., azim=-60)

    inset_axs[2].arrow(xlim*0.3, -xlim*0.5, xlim*0.5, 0, color=(0,0,1), lw=1, head_length=xlim*0.1, head_width=xlim*0.1)  # , overhang=1.0) #, length_includes_head=True)
    inset_axs[2].text(xlim*0.8, -xlim*0.4, r"$\vu{x}$", fontsize=9)
    inset_axs[2].arrow(xlim*0.3, -xlim*0.5, 0, xlim*0.5, color=(1,0,0), lw=1, head_length=xlim*0.1, head_width=xlim*0.1)  # , overhang=1.0) #, length_includes_head=True)
    inset_axs[2].text(xlim*0.4, 0, r"$\vu{z}$", fontsize=9)

    # inset_ax3.view_init(elev=32., azim=-60)
    #inset_axs[0].set_rasterized(True)
    #inset_axs[1].set_rasterized(True)
    #inset_axs[2].set_rasterized(True)

    plt.tight_layout(pad=0.1)
    plt.savefig("./figures/obs_TOC.pdf", format="pdf", dpi=300)
    plt.savefig("./figures/obs_TOC.png", format="png", dpi=300)
    plt.show()
    plt.close()


def plot_obs_x_f_g(tex_lw=240.71031, ppi=72):
    fig = plt.figure(figsize=(tex_lw / ppi * 1, tex_lw / ppi * 0.8))
    plt.rc("text", usetex=True)
    plt.rc("text.latex", preamble=r"\usepackage{physics}")
    ax11 = fig.add_subplot(111)
    # ax12 = fig.add_subplot(122)

    ms = 3.5
    labelpad = 0.1
    # plot lp vs f
    folder = "../data/20240807"
    L = 200
    lss = ['-', "--"]
    markers = ["s", "o", "v"]
    color = ["tomato", "royalblue"]
    kappas = [5.0, 10.0]
    for i in range(len(kappas)):
        kappa = kappas[i]
        ls = lss[i]
        marker = markers[i]
        ls = "None"
        fs = np.arange(0.00, 0.301, 0.02)
        param = [(L, kappa, f, 0.0) for f in fs]
        X, Y, Z, R2, Rg2, Rxx, Rxz = get_obs_data(folder, param)
        ax11.plot(fs, X/L, ls=ls, marker=marker, ms=ms, mfc="None", label=fr"${kappa:.0f}$")
        if (kappa == 10):
            ax11.plot((X/L+0.25/(1-X/L)**2 - 0.25)/kappa, X/L, "-", color="gray", label="theory")
        ax11.plot((X/L+0.25/(1-X/L)**2 - 0.25)/kappa, X/L, "-", color="gray")
        # ax11.plot(fs[2:], 1-0.5/np.sqrt(fs[2:]*kappa), "-", color="gray")

    ax11.legend(title=r"$\kappa$", ncol=1, columnspacing=0.5, handlelength=0.5, handletextpad=0.2, frameon=False, fontsize=9)
    ax11.set_ylabel(r"$X/L$", fontsize=9, labelpad=labelpad)
    ax11.set_xlabel(r"$f$", fontsize=9, labelpad=labelpad)
    ax11.tick_params(which="both", direction="in", top="on", right="on", labelbottom=True, labelleft=True, labelsize=7)
    # ax13.xaxis.set_major_locator(plt.MultipleLocator(0.2))
    # ax13.xaxis.set_minor_locator(plt.MultipleLocator(0.1))
    # ax13.yaxis.set_major_locator(plt.MultipleLocator(20))
    # ax13.yaxis.set_minor_locator(plt.MultipleLocator(10))
    plt.tight_layout(pad=0.05)
    plt.savefig("./figures/obs_X_f.pdf", format="pdf")
    plt.savefig("./figures/obs_X_f.png", format="png")
    plt.show()
    plt.close()


def plot_obs_f_g(tex_lw=240.71031, ppi=72):

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
        # ax12.text(0, kappa-1.7, fr"$\kappa={kappa:.0f}$", fontsize=9)

    ax12.legend(title=rf"$\kappa={kappa:.0f}$", ncol=2, columnspacing=0.5, handlelength=0.5, handletextpad=0.2, frameon=False, fontsize=9)
    ax12.set_ylabel(r"$l_p,l_{p,\theta}$", fontsize=9, labelpad=labelpad)
    ax12.set_xlabel(r"$f$", fontsize=9, labelpad=labelpad)
    ax12.tick_params(which="both", direction="in", top="on", right="on", labelbottom=True, labelleft=True, labelsize=7)

    ax12.xaxis.set_major_locator(plt.MultipleLocator(0.2))
    ax12.xaxis.set_minor_locator(plt.MultipleLocator(0.1))
    ax12.yaxis.set_major_locator(plt.MultipleLocator(1))
    ax12.yaxis.set_minor_locator(plt.MultipleLocator(0.5))
    # ax12.set_ylim(2, 13)

    # plot R2 and Rg vs. f
    kappas = [5.0, 10.0]
    for i in range(len(kappas)):
        kappa = kappas[i]
        ls = lss[i]
        marker = markers[i]
        ls = "None"
        fs = np.arange(0.00, 0.301, 0.02)
        param = [(L, kappa, f, 0.0) for f in fs]
        X, Y, Z, R2, Rg2, Rxx, Rxz = get_obs_data(folder, param)
        # ax21.errorbar(fs, R2/L, yerr=R2_err/L, ls=ls, marker=marker, ms=ms, mfc="None", label=fr"${kappa:.0f}$")
        ax13.plot(fs, R2/L, ls=ls, marker=marker, ms=ms, mfc="None", label=fr"${kappa:.0f}$")
        # ax31.errorbar(fs, Rg, yerr=Rg_err, ls=ls, marker=marker, ms=ms, mfc="None", label=fr"${kappa:.0f}$")
        if kappa == 10:
            # ax14.plot(fs, Rg2/L, ls=ls, marker=marker, ms=ms, mfc="None", label=r"$R_g^2$")
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

    # title=r"$\kappa=10$",
    ax14.legend(title=r"$\mu$", ncol=1, columnspacing=0.5, handlelength=0.5, handletextpad=0.2, frameon=False, fontsize=9)
    ax14.set_ylabel(r"$R_{\mu}/R_g^2$", fontsize=9, labelpad=labelpad)
    ax14.set_xlabel(r"$f$", fontsize=9, labelpad=labelpad)
    ax14.tick_params(which="both", direction="in", top="on", right="on", labelbottom=True, labelleft=True, labelsize=7)

    ax14.xaxis.set_major_locator(plt.MultipleLocator(0.1))
    ax14.xaxis.set_minor_locator(plt.MultipleLocator(0.05))
    ax14.yaxis.set_major_locator(plt.MultipleLocator(0.2))
    ax14.yaxis.set_minor_locator(plt.MultipleLocator(0.1))
    # ax14.set_ylim(-0.25,6)

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
        X, Y, Z, R2, Rg2, Rxx, Rxz = get_obs_data(folder, param)
        # ax23.errorbar(gLs, R2/L, yerr=R2_err/L, ls=ls, marker=marker, ms=ms, mfc="None", label=fr"${f:.1f}$")
        ax23.plot(gLs, R2/L, ls=ls, marker=marker, ms=ms, mfc="None", label=fr"${f:.1f}$")
        # ax24.errorbar(gLs, Rg, yerr=Rg_err, ls=ls,  marker=marker, ms=ms, mfc="None", label=fr"${f:.1f}$")
        if f == 0:
            # ax24.plot(gLs, Rg2/L, ls=ls, marker=marker, ms=ms, mfc="None", label=r"$R_g^2$")
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
    # title=r"$f=0$"
    ax24.legend(title=r"$\mu$", ncol=1, columnspacing=0.5, handlelength=0.5, handletextpad=0.2, frameon=False, fontsize=9)
    ax24.set_ylabel(r"$R_{\mu}/R_g^2$", fontsize=9, labelpad=labelpad)
    ax24.set_xlabel(r"$gL$", fontsize=9, labelpad=labelpad)
    ax24.tick_params(which="both", direction="in", top="on", right="on", labelbottom=True, labelleft=True, labelsize=7)
    # ax24.set_ylim(-0.25,7)
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
    plt.savefig("./figures/obs_f_g.png", format="png")
    plt.show()
    plt.close()
