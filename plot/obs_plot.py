import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches
from scipy.optimize import curve_fit
from matplotlib import rc


def get_obs_data(folder, param):
    Rs, Rgs, R_errs, Rg_errs = [], [], [], []
    for L, kappa, f, g in param:
        finfo = f"L{L}_kappa{kappa:.1f}_f{f:.2f}_g{g:.2f}"
        filename = f"{folder}/obs_{finfo}.csv"
        data = np.genfromtxt(filename, delimiter=',', skip_header=1)
        R, Rg = data[0, 10], data[0, 11]
        R_err, Rg_err = data[1, 10], data[1, 11]
        Rs.append(R)
        Rgs.append(Rg)
        R_errs.append(R_err)
        Rg_errs.append(Rg_err)
    return np.array(Rs), np.array(Rgs), np.array(R_errs), np.array(Rg_errs)


def get_tts_data(folder, param):
    tts = []
    tts_err = []
    spBs = []
    for L, kappa, f, g in param:
        finfo = f"L{L}_kappa{kappa:.1f}_f{f:.2f}_g{g:.2f}"
        filename = f"{folder}/obs_{finfo}.csv"
        data = np.genfromtxt(filename, delimiter=',', skip_header=1)
        tt = data[3, 18:]
        tt_err = data[4, 18:]
        spB = data[5, 18:]
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
    # plt.rc("text", usetex=True)
    # plt.rc("text.latex", preamble=r"\usepackage{physics}")

    ax11 = fig.add_subplot(221)
    ax21 = fig.add_subplot(223)

    ax12 = fig.add_subplot(222)
    ax22 = fig.add_subplot(224)

    # plot tts vs. kappa
    ms = 4
    labelpad = -0.75
    folder = "../data/20240730"
    kappas = [2.0, 4.0, 8.0, 16.0]
    param = [(100, kappa, 0.0, 0.0) for kappa in kappas]

    tts, tts_err, spBs = get_tts_data(folder, param)
    markers = ["o", "x", "s", "+", "d"]
    pltn = 10
    for i in range(len(kappas)):
        ax11.semilogy(spBs[i][:pltn], tts[i][:pltn], marker=markers[i], ms=ms, mfc="None", ls="None", label=fr"${kappas[i]}$")
    ax11.set_ylabel(r"$\left<\cos{\theta}(s)\right>$", fontsize=9, labelpad=labelpad)
    ax11.set_xlabel(r"$s/l_b$", fontsize=9, labelpad=labelpad)
    ax11.tick_params(which="both", direction="in", top="on", right="on", labelbottom=True, labelleft=True, labelsize=7)
    ax11.legend(title=r"$\kappa$", ncol=2, columnspacing=0.5, handlelength=0.5, handletextpad=0.1, frameon=False, fontsize=9)
    ax11.set_ylim(2e-2, 1)
    ax11.xaxis.set_major_locator(plt.MultipleLocator(2))
    ax11.xaxis.set_minor_locator(plt.MultipleLocator(1))
    # ax.yaxis.set_major_locator(plt.MultipleLocator(20))
    # ax.yaxis.set_minor_locator(plt.MultipleLocator(10))

    kappas = np.arange(0.0, 20.01, 1.0)
    L = 100
    param = [(L, kappa, 0.0, 0.0) for kappa in kappas]
    lps, lp_thetas = get_lp_data(folder, param)
    ax21.plot(kappas, lps, marker="x", ms=ms, mfc="None", ls="none", label=r"$l_p$")
    ax21.plot(kappas, lp_thetas, marker="+", ms=ms, mfc="None", ls="none", label=r"$l_{p,\theta}$")
    ax21.plot(kappas, kappas, marker="none", ms=ms, mfc="None", ls="-", lw=1, color="gray")
    ax21.legend(ncol=1, columnspacing=0.5, handlelength=0.5, handletextpad=0.1, frameon=False, fontsize=9)
    ax21.set_ylabel(r"$l_p/l_b$", fontsize=9, labelpad=labelpad)
    ax21.set_xlabel(r"$\kappa/(k_B T)$", fontsize=9, labelpad=labelpad)
    ax21.tick_params(which="both", direction="in", top="on", right="on", labelbottom=True, labelleft=True, labelsize=7)
    ax21.xaxis.set_major_locator(plt.MultipleLocator(4))
    ax21.xaxis.set_minor_locator(plt.MultipleLocator(2))
    ax21.yaxis.set_major_locator(plt.MultipleLocator(4))
    ax21.yaxis.set_minor_locator(plt.MultipleLocator(2))

    # plot R and Rg vs. kappa
    kappas = np.arange(0.0, 20.01, 1.0)
    L = 100
    param = [(L, kappa, 0.0, 0.0) for kappa in kappas]
    Rs, Rgs, R_errs, Rg_errs = get_obs_data(folder, param)
    ax12.errorbar(kappas, np.array(Rs), yerr=np.array(R_errs), ls="-", ms=ms, mfc="None")
    lps = np.array(lps)
    # ax12.plot(kappas, 2*lps*(1-lps/L*(1-np.exp(-L/lps))), marker="none", ms=ms, mfc="None", ls="-", lw=1, color="gray", label=r"$2l_p(1-\frac{l_p}{L}(1-e^{-L/l_p}))$")
    # above equation is for R^2
    ax12.set_ylabel(r"$R/l_b$", fontsize=9, labelpad=labelpad)
    ax12.set_xlabel(r"$\kappa/(k_B T)$", fontsize=9, labelpad=labelpad)
    ax12.tick_params(which="both", direction="in", top="on", right="on", labelbottom=True, labelleft=True, labelsize=7)
    ax12.xaxis.set_major_locator(plt.MultipleLocator(4))
    ax12.xaxis.set_minor_locator(plt.MultipleLocator(2))
    ax22.yaxis.set_major_locator(plt.MultipleLocator(10))
    ax22.yaxis.set_minor_locator(plt.MultipleLocator(5))

    ax22.errorbar(kappas, Rgs, yerr=Rg_errs, ls="-", ms=ms, mfc="None")
    ax22.set_ylabel(r"$R_g/l_b$", fontsize=9, labelpad=labelpad)
    ax22.set_xlabel(r"$\kappa/(k_B T)$", fontsize=9, labelpad=labelpad)
    ax22.tick_params(which="both", direction="in", top="on", right="on", labelbottom=True, labelleft=True, labelsize=7)
    ax22.xaxis.set_major_locator(plt.MultipleLocator(4))
    ax22.xaxis.set_minor_locator(plt.MultipleLocator(2))
    ax22.yaxis.set_major_locator(plt.MultipleLocator(4))
    ax22.yaxis.set_minor_locator(plt.MultipleLocator(2))

    plt.tight_layout(pad=0.05)
    plt.savefig("./figures/obs_kappa.pdf", format="pdf")
    plt.show()
    plt.close()


def plot_obs_f_g(tex_lw=240.71031, ppi=72):
    fig = plt.figure(figsize=(tex_lw / ppi * 1, tex_lw / ppi * 1.2))
    # plt.rc("text", usetex=True)
    # plt.rc("text.latex", preamble=r"\usepackage{physics}")

    ax11 = fig.add_subplot(321)
    ax21 = fig.add_subplot(323, sharex=ax11)
    ax31 = fig.add_subplot(325, sharex=ax11)

    ax12 = fig.add_subplot(322)
    ax22 = fig.add_subplot(324, sharex=ax12)
    ax32 = fig.add_subplot(326, sharex=ax12)

    ms = 4
    labelpad = -0.75
    # plot lp vs f
    folder = "../data/20240730"
    L = 100

    kappas = [5.0, 10.0]
    lss = ['-', "--"]
    color = ["tomato", "royalblue"]
    for i in range(len(kappas)):
        kappa = kappas[i]
        ls = lss[i]
        fs = np.arange(0.00, 0.401, 0.02)
        param = [(L, kappa, f, 0.0) for f in fs]
        lps, lp_thetas = get_lp_data(folder, param)
        if i == 0:
            ax11.plot(fs, lps, color="tomato", ls=ls, label=fr"$l_p$")
            ax11.plot(fs, lp_thetas, color="royalblue", ls=ls, label=r"$l_{p,\theta}$")
        else:
            ax11.plot(fs, lps, color="tomato", ls=ls)
            ax11.plot(fs, lp_thetas, color="royalblue", ls=ls)
        ax11.text(0, kappa-1.5, fr"$\kappa/(k_B T)={kappa}$", fontsize=9)

    ax11.legend(ncol=2, columnspacing=0.5, handlelength=0.5, handletextpad=0.1, frameon=False, fontsize=9)
    ax11.set_ylabel(r"$l_p/l_b$", fontsize=9, labelpad=labelpad)
    # ax11.set_xlabel(r"$fl_b/(k_B T)$", fontsize=9, labelpad=labelpad)
    ax11.tick_params(which="both", direction="in", top="on", right="on", labelbottom=False, labelleft=True, labelsize=7)

    # ax11.xaxis.set_major_locator(plt.MultipleLocator(0.2))
    # ax11.xaxis.set_minor_locator(plt.MultipleLocator(0.1))
    ax11.yaxis.set_major_locator(plt.MultipleLocator(4))
    ax11.yaxis.set_minor_locator(plt.MultipleLocator(2))
    ax11.set_ylim(2, 13)

    # plot R and Rg vs. f
    for i in range(len(kappas)):
        kappa = kappas[i]
        ls = lss[i]
        fs = np.arange(0.00, 0.401, 0.02)
        param = [(L, kappa, f, 0.0) for f in fs]
        R, Rg, R_err, Rg_err = get_obs_data(folder, param)
        ax21.errorbar(fs, R, yerr=R_err, ls=ls, label=fr"${kappa}$")
        ax31.errorbar(fs, Rg, yerr=Rg_err, ls=ls, label=fr"${kappa}$")

    ax21.legend(title=r"$\kappa/(k_B T)$", ncol=1, columnspacing=0.5, handlelength=0.5, handletextpad=0.1, frameon=False, fontsize=9)
    ax21.set_ylabel(r"$R/l_b$", fontsize=9, labelpad=labelpad)
    # ax21.set_xlabel(r"$fl_b/(k_B T)$", fontsize=9, labelpad=labelpad)
    ax21.tick_params(which="both", direction="in", top="on", right="on", labelbottom=False, labelleft=True, labelsize=7)
    # ax21.xaxis.set_major_locator(plt.MultipleLocator(0.2))
    # ax21.xaxis.set_minor_locator(plt.MultipleLocator(0.1))
    ax21.yaxis.set_major_locator(plt.MultipleLocator(10))
    ax21.yaxis.set_minor_locator(plt.MultipleLocator(5))

    ax31.legend(title=r"$\kappa/(k_B T)$", ncol=1, columnspacing=0.5, handlelength=0.5, handletextpad=0.1, frameon=False, fontsize=9)
    ax31.set_ylabel(r"$R_g/l_b$", fontsize=9, labelpad=labelpad)
    ax31.set_xlabel(r"$fl_b/(k_B T)$", fontsize=9, labelpad=labelpad)
    ax31.tick_params(which="both", direction="in", top="on", right="on", labelbottom=True, labelleft=True, labelsize=7)
    ax31.xaxis.set_major_locator(plt.MultipleLocator(0.1))
    ax31.xaxis.set_minor_locator(plt.MultipleLocator(0.05))
    ax31.yaxis.set_major_locator(plt.MultipleLocator(4))
    ax31.yaxis.set_minor_locator(plt.MultipleLocator(2))

    # plot lp vs g
    L=100
    kappa = 10
    fs = [0.00]
    gs = np.arange(0.00, 0.201, 0.01)
    lss = ['-', "--"]
    color = ["tomato", "royalblue"]
    for i in range(len(fs)):
        f = fs[i]
        ls = lss[i]
        param = [(L, kappa, f, g) for g in gs]
        lps, lp_thetas = get_lp_data(folder, param)
        if i == 0:
            ax12.plot(gs, lps, color="tomato", ls=ls, label=fr"$l_p$")
            ax12.plot(gs, lp_thetas, color="royalblue", ls=ls, label=r"$l_{p,\theta}$")
        else:
            ax12.plot(gs, lps, color="tomato", ls=ls)
            ax12.plot(gs, lp_thetas, color="royalblue", ls=ls)
        #ax12.text(0, kappa-1.5, fr"$f={f}$", fontsize=9)

    ax12.legend(ncol=2, columnspacing=0.5, handlelength=0.5, handletextpad=0.1, frameon=False, fontsize=9)
    ax12.set_ylabel(r"$l_p/l_b$", fontsize=9, labelpad=labelpad)
    # ax12.set_xlabel(r"$fl_b/(k_B T)$", fontsize=9, labelpad=labelpad)
    ax12.tick_params(which="both", direction="in", top="on", right="on", labelbottom=False, labelleft=True, labelsize=7)

    # ax12.xaxis.set_major_locator(plt.MultipleLocator(0.2))
    # ax12.xaxis.set_minor_locator(plt.MultipleLocator(0.1))
    ax12.yaxis.set_major_locator(plt.MultipleLocator(10))
    ax12.yaxis.set_minor_locator(plt.MultipleLocator(5))
    #ax12.set_ylim(2, 13)

    # plot R and Rg vs. g
    fs = [0.00, 0.20]
    for i in range(len(fs)):
        f = fs[i]
        ls = lss[i]
        param = [(L, kappa, f, g) for g in gs]
        R, Rg, R_err, Rg_err = get_obs_data(folder, param)
        ax22.errorbar(gs, R, yerr=R_err, ls=ls, label=fr"${f}$")
        ax32.errorbar(gs, Rg, yerr=Rg_err, ls=ls, label=fr"${f}$")

    ax22.legend(title=r"$f l_b/(k_B T)$", ncol=1, columnspacing=0.5, handlelength=0.5, handletextpad=0.1, frameon=False, fontsize=9)
    ax22.set_ylabel(r"$R/l_b$", fontsize=9, labelpad=labelpad)
    # ax22.set_xlabel(r"$fl_b/(k_B T)$", fontsize=9, labelpad=labelpad)
    ax22.tick_params(which="both", direction="in", top="on", right="on", labelbottom=False, labelleft=True, labelsize=7)
    # ax22.xaxis.set_major_locator(plt.MultipleLocator(0.2))
    # ax22.xaxis.set_minor_locator(plt.MultipleLocator(0.1))
    ax22.yaxis.set_major_locator(plt.MultipleLocator(10))
    ax22.yaxis.set_minor_locator(plt.MultipleLocator(5))

    ax32.legend(title=r"$f l_b/(k_B T)$", ncol=1, columnspacing=0.5, handlelength=0.5, handletextpad=0.1, frameon=False, fontsize=9)
    ax32.set_ylabel(r"$R_g/l_b$", fontsize=9, labelpad=labelpad)
    ax32.set_xlabel(r"$gl^2_b/(k_B T)$", fontsize=9, labelpad=labelpad)
    ax32.tick_params(which="both", direction="in", top="on", right="on", labelbottom=True, labelleft=True, labelsize=7)
    ax32.xaxis.set_major_locator(plt.MultipleLocator(0.05))
    ax32.xaxis.set_minor_locator(plt.MultipleLocator(0.025))
    ax32.yaxis.set_major_locator(plt.MultipleLocator(4))
    ax32.yaxis.set_minor_locator(plt.MultipleLocator(2))


    plt.tight_layout(pad=0.05)
    plt.savefig("./figures/obs_f_g.pdf", format="pdf")
    plt.show()
    plt.close()
