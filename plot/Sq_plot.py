import matplotlib.pyplot as plt
import numpy as np


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
    for L, kappa, f, gL in param:
        finfo = f"L{L}_kappa{kappa:.1f}_f{f:.2f}_gL{gL:.2f}"
        filename = f"{folder}/obs_{finfo}.csv"
        data = np.genfromtxt(filename, delimiter=',', skip_header=1)
        Sq = data[0, 21:]
        qB = data[2, 21:]
        Sqs.append(Sq)
        qBs.append(qB)
    return np.array(Sqs), np.array(qBs)


def plot_Sq(tex_lw=240.71031, ppi=72):
    fig = plt.figure(figsize=(tex_lw / ppi * 1, tex_lw / ppi * 0.8))
    plt.rc("text", usetex=True)
    plt.rc("text.latex", preamble=r"\usepackage{physics}")
    ax00 = plt.subplot2grid((2, 2), (0, 0))
    ax10 = plt.subplot2grid((2, 2), (0, 1), sharex=ax00)
    ax01 = plt.subplot2grid((2, 2), (1, 0), sharex=ax00)
    ax11 = plt.subplot2grid((2, 2), (1, 1), sharex=ax00)

    ms = 4
    labelpad = -0.0
    # plot Sq for various kappa
    folder = "../data/20240807"
    L = 200
    kappas = [2.0, 4.0, 8.0, 16.0]
    # kappas = [4.0, 8.0, 14.0, 20.0]
    param = [(L, kappa, 0.0, 0.0) for kappa in kappas]
    Sqs, qBs = get_Sq_data(folder, param)
    Sq_rod = calc_Sq_discrete_infinite_thin_rod(qBs[0], L)
    for i in range(len(kappas)):
        ax00.loglog(qBs[i], Sqs[i], "-", lw=1, label=f"{kappas[i]:.0f}")
    ax00.loglog(qBs[0], Sq_rod, "k--", lw=1, label="rod")
    # ax00.set_xlabel(r"$q$", fontsize=9, labelpad=labelpad)
    ax00.set_ylabel(r"$S(Q)$", fontsize=9, labelpad=labelpad)
    ax00.legend(title=r"$\kappa$", ncol=2, columnspacing=0.5, handlelength=0.5, handletextpad=0.2, frameon=False, fontsize=9)
    ax00.tick_params(which="both", direction="in", top="on", right="on", labelbottom=False, labelleft=True, labelsize=7)
    # ax00.xaxis.set_major_locator(plt.MultipleLocator(0.1))
    # ax00.xaxis.set_minor_locator(plt.MultipleLocator(0.05))
    # ax00.yaxis.set_major_locator(plt.MultipleLocator(4))
    # ax00.yaxis.set_minor_locator(plt.MultipleLocator(2))

    # plot Del Sq for various kappa
    for i in range(len(kappas)):
        ax10.semilogx(qBs[i], np.log(Sqs[i]/Sq_rod), "-", lw=1, label=f"{kappas[i]:.0f}")
    # ax10.set_xlabel(r"$q$", fontsize=9, labelpad=labelpad)
    ax10.set_ylabel(r"$\Delta S(Q)$", fontsize=9, labelpad=labelpad)
    ax10.legend(title=r"$\kappa$", ncol=1, columnspacing=0.5, handlelength=0.5, handletextpad=0.2, frameon=False, fontsize=9)
    ax10.tick_params(which="both", direction="in", top="on", right="on", labelbottom=False, labelleft=True, labelsize=7)
    ax10.yaxis.set_major_locator(plt.MultipleLocator(0.4))
    ax10.yaxis.set_minor_locator(plt.MultipleLocator(0.2))

    folder = "../data/20240821"
    # plot Del Sq for various f
    fs = [0.00, 0.10, 0.20, 0.30]
    param = [(L, 10.0, f, 0.0) for f in fs]
    Sqs, qBs = get_Sq_data(folder, param)
    Sq_rod = calc_Sq_discrete_infinite_thin_rod(qBs[0], L)
    for i in range(len(fs)):
        ax01.semilogx(qBs[i], np.log(Sqs[i]/Sq_rod), "-", lw=1, label=f"{fs[i]:.1f}")
    ax01.set_xlabel(r"$Q$", fontsize=9, labelpad=labelpad)
    ax01.set_ylabel(r"$\Delta S(Q)$", fontsize=9, labelpad=labelpad)
    ax01.legend(title=r"$f$", ncol=1, columnspacing=0.5, handlelength=0.5, handletextpad=0.2, frameon=False, fontsize=9)
    ax01.tick_params(which="both", direction="in", top="on", right="on", labelbottom=True, labelleft=True, labelsize=7)
    ax01.yaxis.set_major_locator(plt.MultipleLocator(0.4))
    ax01.yaxis.set_minor_locator(plt.MultipleLocator(0.2))

    # plot Del Sq for various g
    # gs = np.arange(0.00, 0.201, 0.01)
    folder = "../data/20240820"
    gLs = [0.00, 0.50, 1.0, 1.50]
    param = [(L, 10.0, 0.0, gL) for gL in gLs]
    Sqs, qBs = get_Sq_data(folder, param)
    for i in range(len(gLs)):
        ax11.semilogx(qBs[i], np.log(Sqs[i]/Sq_rod), "-", lw=1, label=f"{gLs[i]:.1f}")
    ax11.set_xlabel(r"$Q$", fontsize=9, labelpad=labelpad)
    ax11.set_ylabel(r"$\Delta S(Q)$", fontsize=9, labelpad=labelpad)
    ax11.legend(title=r"$\gamma L$", ncol=1, columnspacing=0.5, handlelength=0.5, handletextpad=0.2, frameon=False, fontsize=9)
    ax11.tick_params(which="both", direction="in", top="on", right="on", labelbottom=True, labelleft=True, labelsize=7)
    ax11.yaxis.set_major_locator(plt.MultipleLocator(0.4))
    ax11.yaxis.set_minor_locator(plt.MultipleLocator(0.2))

    annotation = [r"$(a)$", r"$(b)$", r"$(c)$", r"$(d)$"]
    for ax in [ax00, ax10, ax01, ax11]:
        ax.text(0.8, 0.85, annotation.pop(0), fontsize=9, transform=ax.transAxes)

    plt.tight_layout(pad=0.05)

    plt.savefig("./figures/Sq.pdf", format="pdf")
    plt.show()
    plt.close()


def get_Sq2D_data(filename, CS=False):
    if CS:
        data = np.genfromtxt(filename, delimiter=',', skip_header=0)
        qB = np.linspace(-0.25*np.pi, 0.25*np.pi, 51)
        Sq2D = data[:,:]
    else:
        data = np.genfromtxt(filename, delimiter=',', skip_header=1)
        qB = data[2, 21:]
        Sq2D = data[6:, 21:]
    return np.array(Sq2D), np.array(qB)


def plot_Sq2D(tex_lw=240.71031, ppi=72):
    fig = plt.figure(figsize=(tex_lw / ppi * 1, tex_lw / ppi * 0.45))
    plt.rc("text", usetex=True)
    plt.rc("text.latex", preamble=r"\usepackage{physics}")
    ax00 = plt.subplot2grid((1, 2), (0, 0))
    ax01 = plt.subplot2grid((1, 2), (0, 1))

    # on-lattice model
    Sq2D, qB = get_Sq2D_data("../data/scratch_local/CH_spect/L_200_kappa_10_SC.csv", True)
    print("Sq2D.shape",Sq2D.shape)
    print("qB.shape",qB.shape)
    qBx, qBy = np.meshgrid(qB, qB)
    print(np.min(Sq2D), np.max(Sq2D))
    ax00.pcolormesh(qBx, qBy, np.log10(Sq2D), vmax=0, vmin=-3, cmap="rainbow", shading='gouraud')
    ax00.contour(qBx, qBy, np.log10(Sq2D), vmax=0, vmin=-3, levels=np.linspace(-3, 0, 6), colors="white", linewidths=1, linestyle=":")
    #ax00.pcolormesh(qqx*N_backbone,qqy*N_backbone,np.log(S_q_2D_list[1,:,:,0]).T, vmax=0,vmin=-6, shading='gouraud')
    ax00.set_xlabel(r"$Q_xB$", fontsize=9, labelpad=-0.0)
    ax00.set_ylabel(r"$Q_yB$", fontsize=9, labelpad=-1)
    ax00.tick_params(which="both", direction="in", top="on", right="on", labelbottom=True, labelleft=True, labelsize=7)
    ax00.xaxis.set_major_locator(plt.MultipleLocator(0.5))
    ax00.xaxis.set_minor_locator(plt.MultipleLocator(0.25))
    ax00.yaxis.set_major_locator(plt.MultipleLocator(0.5))
    ax00.yaxis.set_minor_locator(plt.MultipleLocator(0.25))

    # off-lattice model
    #Sq2D, qB = get_Sq2D_data("../data/scratch_local/20240830/obs_L50_kappa2.0_f0.00_gL0.00.csv")
    Sq2D, qB = get_Sq2D_data("../data/scratch_local/20240830/obs_L50_kappa5.0_f0.00_gL0.00.csv")
    #Sq2D, qB = get_Sq2D_data("../data/scratch_local/20240829/obs_L200_kappa10.0_f0.00_gL0.00.csv")
    #Sq2D, qB = get_Sq2D_data("../data/20240829/obs_L200_kappa5.0_f0.00_gL0.00.csv")
    qBx, qBy = np.meshgrid(qB, qB)
    print("qBx.shape", qBx.shape)
    #Sq2D = np.sqrt(Sq2D.T)
    Sq2D = Sq2D.T
    print("Sq2D.shape",Sq2D.shape)
    #ax01.contourf(qBx, qBy, Sq2D, levels=20, cmap="jet", norm="log")
    print(np.min(Sq2D), np.max(Sq2D))
    ax01.pcolormesh(qBx, qBy, np.log10(Sq2D), vmax=0, vmin=-3, cmap="rainbow", shading='gouraud')
    #ax01.pcolormesh(qBx, qBy, Sq2D, cmap="rainbow", shading='gouraud')
    ax01.contour(qBx, qBy, np.log10(Sq2D), vmax=0, vmin=-3, levels=np.linspace(-3, 0, 6), colors="white", linewidths=1, linestyle=":")
    ax01.set_xlabel(r"$Q_xB$", fontsize=9, labelpad=-0.0)
    ax01.set_ylabel(r"$Q_yB$", fontsize=9, labelpad=-1)
    ax01.tick_params(which="both", direction="in", top="on", right="on", labelbottom=True, labelleft=True, labelsize=7)
    ax01.xaxis.set_major_locator(plt.MultipleLocator(0.5))
    ax01.xaxis.set_minor_locator(plt.MultipleLocator(0.25))
    ax01.yaxis.set_major_locator(plt.MultipleLocator(0.5))
    ax01.yaxis.set_minor_locator(plt.MultipleLocator(0.25))

    ax00.text(0.8, 0.15, r"$(a)$", fontsize=9, transform=ax00.transAxes, color='white')
    ax01.text(0.8, 0.15, r"$(b)$", fontsize=9, transform=ax01.transAxes, color='white')

    plt.tight_layout(pad=0.1)
    plt.subplots_adjust(wspace=0.35)
    plt.savefig("./figures/Sq2D.pdf", format="pdf")
    plt.savefig("./figures/Sq2D.png", format="png", dpi=300)
    plt.show()
    plt.close()
