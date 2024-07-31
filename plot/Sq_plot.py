import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches
from scipy.optimize import curve_fit


def calc_Sq_discrete_infinite_thin_rod(q, L):
    # numereical calculation
    Sq = [1.0/L for i in range(len(q))]
    for k in range(len(q)):
        Sqk = 0
        qk = q[k]
        for i in range(L-1):
            for j in range(i+1, L):
                Sqk += 2.0*np.sin(qk*(i-j))/(qk*(i-j))/(L*L)
        Sq[k] += Sqk
    return np.array(Sq)


def get_Sq_data(folder, param):
    Sqs = []
    qBs = []
    for L, kappa, f, g in param:
        finfo = f"L{L}_kappa{kappa:.1f}_f{f:.2f}_g{g:.2f}"
        filename = f"{folder}/obs_{finfo}.csv"
        data = np.genfromtxt(filename, delimiter=',', skip_header=1)
        Sq = data[0, 18:]
        qB = data[2, 18:]
        Sqs.append(Sq)
        qBs.append(qB)
    return np.array(Sqs), np.array(qBs)


def plot_Sq(tex_lw=240.71031, ppi=72):
    fig = plt.figure(figsize=(tex_lw / ppi * 1, tex_lw / ppi * 0.85))
    # plt.rc("text", usetex=True)
    # plt.rc("text.latex", preamble=r"\usepackage{physics}")
    ax00 = plt.subplot2grid((2, 2), (0, 0))
    ax10 = plt.subplot2grid((2, 2), (0, 1),sharex=ax00)
    ax11 = plt.subplot2grid((2, 2), (1, 0),sharex=ax00)
    ax12 = plt.subplot2grid((2, 2), (1, 1),sharex=ax00)

    ms = 4
    labelpad = -0.75
    # plot Sq for various kappa
    folder = "../data/20240730"
    L = 100
    kappas = [2.0, 4.0, 8.0, 16.0]
    param = [(L, kappa, 0.0, 0.0) for kappa in kappas]
    Sqs, qBs = get_Sq_data(folder, param)
    Sq_rod = calc_Sq_discrete_infinite_thin_rod(qBs[0], L)
    for i in range(len(kappas)):
        ax00.loglog(qBs[i], Sqs[i], "-", label=f"{kappas[i]}")
    ax00.loglog(qBs[0], Sq_rod, "k--", label="rod")
    #ax00.set_xlabel(r"$ql_b$", fontsize=9, labelpad=labelpad)
    ax00.set_ylabel(r"$S(ql_b)$", fontsize=9, labelpad=labelpad)
    ax00.legend(title=r"$\kappa/(k_B T)$", ncol=2, columnspacing=0.5, handlelength=0.5, handletextpad=0.1, frameon=False, fontsize=9)
    ax00.tick_params(which="both", direction="in", top="on", right="on", labelbottom=False, labelleft=True, labelsize=7)
    #ax00.xaxis.set_major_locator(plt.MultipleLocator(0.1))
    #ax00.xaxis.set_minor_locator(plt.MultipleLocator(0.05))
    #ax00.yaxis.set_major_locator(plt.MultipleLocator(4))
    #ax00.yaxis.set_minor_locator(plt.MultipleLocator(2))


    # plot Del Sq for various kappa
    for i in range(len(kappas)):
        ax10.semilogx(qBs[i], Sqs[i]/Sq_rod, "-", label=f"{kappas[i]}")
    #ax10.set_xlabel(r"$ql_b$", fontsize=9, labelpad=labelpad)
    ax10.set_ylabel(r"$\Delta S(ql_b)$", fontsize=9, labelpad=labelpad)
    ax10.legend(title=r"$\kappa/(k_B T)$", ncol=1, columnspacing=0.5, handlelength=0.5, handletextpad=0.1, frameon=False, fontsize=9)
    ax10.tick_params(which="both", direction="in", top="on", right="on", labelbottom=False, labelleft=True, labelsize=7)

    # plot Del Sq for various f
    fs = [0.00, 0.10, 0.20, 0.30]
    param = [(L, 10.0, f, 0.0) for f in fs]
    Sqs, qBs = get_Sq_data(folder, param)
    for i in range(len(fs)):
        ax11.semilogx(qBs[i], Sqs[i]/Sq_rod, "-", label=f"{fs[i]:.2f}")
    ax11.set_xlabel(r"$ql_b$", fontsize=9, labelpad=labelpad)
    ax11.set_ylabel(r"$\Delta S(ql_b)$", fontsize=9, labelpad=labelpad)
    ax11.legend(title=r"$fl_b/(k_B T)$", ncol=1, columnspacing=0.5, handlelength=0.5, handletextpad=0.1, frameon=False, fontsize=9)
    ax11.tick_params(which="both", direction="in", top="on", right="on", labelbottom=True, labelleft=True, labelsize=7)

    # plot Del Sq for various g
    gs = np.arange(0.00, 0.201, 0.01)
    gs = [0.00, 0.02, 0.04, 0.06]
    param = [(L, 10.0, 0.0, g) for g in gs]
    Sqs, qBs = get_Sq_data(folder, param)
    for i in range(len(gs)):
        ax12.semilogx(qBs[i], Sqs[i]/Sq_rod, "-", label=f"{gs[i]:.2f}")
    ax12.set_xlabel(r"$ql_b$", fontsize=9, labelpad=labelpad)
    ax12.set_ylabel(r"$\Delta S(ql_b)$", fontsize=9, labelpad=labelpad)
    ax12.legend(title=r"$gl_b^2(k_B T)$", ncol=1, columnspacing=0.5, handlelength=0.5, handletextpad=0.1, frameon=False, fontsize=9)
    ax12.tick_params(which="both", direction="in", top="on", right="on", labelbottom=True, labelleft=True, labelsize=7)

    plt.tight_layout(pad=0.05)

    plt.savefig("./figures/Sq.pdf", format="pdf")
    plt.show()
    plt.close()
