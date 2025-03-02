from scipy import stats
import numpy as np
import matplotlib.pyplot as plt
import os


def get_all_feature_Sq2D_data(folder, parameters):
    all_kappa, all_f, all_gL = [], [], []  # force related
    all_R2, all_Rg2 = [], []  # R2, Rg2 related
    all_Sxx, all_Syy, all_Szz, all_Sxy, all_Sxz, all_Syz = [], [], [], [], [], []  # S related

    all_Sq2D_flatten = []
    qB = []
    for L, run_num in parameters:
        filename = f"{folder}/obs_L{L}_random_run{run_num}.csv"
        if not os.path.exists(filename):
            print(f"File not found: {filename}")
            continue
        data = np.genfromtxt(filename, delimiter=",", skip_header=1)
        kappa, f, g = data[0, 2:5]
        R2, Rg2 = data[0, 13:15]
        Sxx, Syy, Szz, Sxy, Sxz, Syz = data[0, 15:21]
        all_kappa.append(kappa)
        all_f.append(f)
        all_gL.append(g * L)
        all_R2.append(R2)
        all_Rg2.append(Rg2)
        all_Sxx.append(Sxx / Rg2)
        all_Syy.append(Syy / Rg2)
        all_Szz.append(Szz / Rg2)
        all_Sxy.append(Sxy / Rg2)
        all_Sxz.append(Sxz / Rg2)
        all_Syz.append(Syz / Rg2)

        qB = data[2, 21:]
        Sq2D = data[6:, 21:]
        all_Sq2D_flatten.append(Sq2D.flatten())

    all_feature = np.array([all_kappa, all_f, all_gL, all_R2, all_Rg2, all_Sxx, all_Syy, all_Szz, all_Sxy, all_Sxz, all_Syz])
    all_feature_name = ["kappa", "f", "gL", "R2", "Rg2", "Sxx", "Syy", "Szz", "Sxy", "Sxz", "Syz"]
    all_Sq2D_flatten = np.array(all_Sq2D_flatten)
    qB = np.array(qB)
    return all_feature.T, all_feature_name, all_Sq2D_flatten, qB


def read_SVD_data(folder):
    data = np.loadtxt(f"{folder}/data_svd_projection.txt", skiprows=1, delimiter=",", unpack=True)
    print("data.shape", data.shape)
    print("data[:,-3:].shape", data[-3:, :].shape)
    sqv0, sqv1, sqv2 = data[-3:, :]

    all_feature = data[:-3, :]
    return all_feature, sqv0, sqv1, sqv2


def calc_sum_log_likelihood(x, k):
    mu1 = np.mean(x[:k])
    mu2 = np.mean(x[k:])
    sigma1 = np.std(x[:k])
    sigma2 = np.std(x[k:])
    q = len(x)
    sigma = np.sqrt(((k - 1) * sigma1**2 + (q - k - 1) * sigma2**2) / (q - 2))
    log_norm_x1 = stats.norm.logpdf(x[:k], mu1, sigma)
    log_norm_x2 = stats.norm.logpdf(x[k:], mu2, sigma)
    return np.sum(log_norm_x1) + np.sum(log_norm_x2)


def calc_maximum_profile_likelihood(svdS):
    k = range(1, 20)
    pll = np.zeros(len(k))
    for i in range(len(k)):
        pll[i] = calc_sum_log_likelihood(svdS, k[i])
    return k, pll


def plot_SVD_data(tex_lw=240.71031, ppi=72):
    folder = "../data/20240924_random"
    #folder = "../data/20240918_random" # 18 used the wrong convention
    rand_max = 1024
    L = 200
    parameters = [[L, rand_num] for rand_num in range(rand_max)]
    all_features, all_feature_names, all_Sq2D_flatten, qB = get_all_feature_Sq2D_data(folder, parameters)
    svd = np.linalg.svd(all_Sq2D_flatten)

    plt.figure(figsize=(tex_lw / ppi * 1, tex_lw / ppi * 0.8))
    plt.rc("text", usetex=True)
    plt.rc("text.latex", preamble=r"\usepackage{physics}")
    ax00 = plt.subplot2grid((2, 6), (0, 0), colspan=6)
    # ax01 = plt.subplot2grid((2, 6), (0, 3), colspan=3)
    ax10 = plt.subplot2grid((2, 6), (1, 0), colspan=2)
    ax11 = plt.subplot2grid((2, 6), (1, 2), colspan=2, sharex=ax10, sharey=ax10)
    ax12 = plt.subplot2grid((2, 6), (1, 4), colspan=2, sharex=ax10, sharey=ax10)

    # Subplot for svd.S
    # ax00.plot(range(len(svd.S[:51])), svd.S[:51], "x", markersize=5, markerfacecolor='none', label=r"$\Sigma$")
    ax00.semilogx(range(1, len(svd.S) + 1), svd.S, "x", markersize=5, markerfacecolor="none", label=r"$\Sigma$")
    ax00.plot(range(1, len(svd.S) + 1)[:3], svd.S[:3], "ro", markersize=5, markerfacecolor="none")
    ax00.set_xlabel("SVR", fontsize=9, labelpad=0)
    ax00.set_ylabel(r"$\Sigma$", fontsize=9, labelpad=0)
    ax00.tick_params(which="both", direction="in", top="on", right="on", labelbottom=True, labelleft=True, labelsize=7)
    # ax00.xaxis.set_major_locator(plt.MultipleLocator(20))
    # ax00.xaxis.set_minor_locator(plt.MultipleLocator(10))

    # subplot for profile likelihood
    """
    k, pll = calc_maximum_profile_likelihood(svd.S)
    ax01.plot(k, pll, "x", markersize=5, markerfacecolor='none')
    ax01.plot(k[:3], pll[:3], "x", markersize=5, markerfacecolor='none')
    ax01.set_xlabel("SVR", fontsize=9, labelpad=0)
    ax01.set_ylabel("Profile Likelihood", fontsize=9, labelpad=0)
    ax01.tick_params(which="both", direction="in", top="on", right="on", labelbottom=True, labelleft=True, labelsize=7)
    """

    # Subplot for svd.U
    qBx, qBy = np.meshgrid(qB, qB)
    vmin, vmax = -0.2, 0.2
    print("np.minimum(svd.Vh[0]), np.maximum(svd.Vh[0])", svd.Vh[0].min(), svd.Vh[0].max())
    cf = ax10.contourf(qBx, qBy, svd.Vh[0].reshape(np.shape(qBx)), vmin=vmin, vmax=vmax, levels=np.linspace(vmin, vmax, 10), cmap="rainbow")
    ax10.text(0.6, 0.8, r"$V0$", fontsize=9, transform=ax10.transAxes)
    ax10.set_xlabel(r"$Q_x$", fontsize=9, labelpad=0)
    ax10.set_ylabel(r"$Q_z$", fontsize=9, labelpad=0)
    ax10.tick_params(which="both", direction="in", top="on", right="on", labelbottom=True, labelleft=True, labelsize=7)

    print("np.minimum(svd.Vh[1]), np.maximum(svd.Vh[1])", svd.Vh[1].min(), svd.Vh[1].max())
    ax11.contourf(qBx, qBy, svd.Vh[1].reshape(np.shape(qBx)), vmin=vmin, vmax=vmax, levels=np.linspace(vmin, vmax, 10), cmap="rainbow")
    ax11.set_xlabel(r"$Q_x$", fontsize=9, labelpad=0)
    ax11.text(0.6, 0.8, r"$V1$", fontsize=9, transform=ax11.transAxes)
    ax11.tick_params(which="both", direction="in", top="on", right="on", labelbottom=True, labelleft=False, labelsize=7)

    print("np.minimum(svd.Vh[2]), np.maximum(svd.Vh[2])", svd.Vh[2].min(), svd.Vh[2].max())
    cf12 = ax12.contourf(qBx, qBy, svd.Vh[2].reshape(np.shape(qBx)), vmin=vmin, vmax=vmax, levels=np.linspace(vmin, vmax, 10), cmap="rainbow")
    ax12.set_xlabel(r"$Q_x$", fontsize=9, labelpad=0)
    ax12.text(0.6, 0.8, r"$V2$", fontsize=9, transform=ax12.transAxes)
    ax12.tick_params(which="both", direction="in", top="on", right="on", labelbottom=True, labelleft=False, labelsize=7)

    cbar = plt.colorbar(cf12, ax=ax12, fraction=0.046, pad=0.02)
    cbar.set_ticks(np.linspace(vmin, vmax, 5))
    cbar.ax.tick_params(labelsize=7, direction="in", rotation=45)

    ax00.text(0.9, 0.15, r"$(a)$", transform=ax00.transAxes, fontsize=9)
    ax10.text(0.7, 0.15, r"$(b)$", transform=ax10.transAxes, fontsize=9)
    ax11.text(0.7, 0.15, r"$(c)$", transform=ax11.transAxes, fontsize=9)
    ax12.text(0.7, 0.15, r"$(d)$", transform=ax12.transAxes, fontsize=9)

    plt.tight_layout(pad=0.1)
    plt.savefig("figures/SVD.pdf", dpi=300)
    plt.savefig("figures/SVD.png", dpi=300)  # png figure for slides making

    plt.show()


def plot_SVD_feature_data(tex_lw=240.71031, ppi=72):
    folder = "../data/20240924_random"
    #folder = "../data/20240918_random"
    all_feature, sqv0, sqv1, sqv2 = read_SVD_data(folder)
    all_feature_name = ["kappa", "f", "gL", "R2", "Rg2", "Sxx", "Syy", "Szz", "Sxy", "Sxz", "Syz"]
    feature_to_tex = {"kappa": r"$\kappa$", "f": r"$f$", "gL": r"$\gamma L$", "R2": r"$R^2/L^2$", "Rg2": r"$R_g^2/L$", "Sxz": r"$R_{xz}/R_g^2$"}

    fig = plt.figure(figsize=(tex_lw / ppi * 1, tex_lw / ppi * 1.2))
    plt.rc("text", usetex=True)
    plt.rc("text.latex", preamble=r"\usepackage{physics}")

    # kappa, f, gL
    ax11 = fig.add_subplot(321, projection="3d")
    ax12 = fig.add_subplot(322, projection="3d")
    ax13 = fig.add_subplot(323, projection="3d")
    # R2, Rg2, Sxz
    ax21 = fig.add_subplot(324, projection="3d")
    ax22 = fig.add_subplot(325, projection="3d")
    ax23 = fig.add_subplot(326, projection="3d")

    axs = [ax11, ax12, ax13, ax21, ax22, ax23]

    i = -1
    for feature_name, feature_tex in feature_to_tex.items():
        i += 1
        # cbar_major_locator = [0.4, 0.03, 5, 5, 40, 0.2, 0.1, 40]
        if feature_name not in all_feature_name:
            continue
        print("plotting: ", feature_name)

        feature_index = all_feature_name.index(feature_name)
        mu = all_feature[feature_index, :]
        L = 200
        if feature_name == "R2":
            mu /= L * L
        if feature_name == "Rg2":
            mu /= L
        ax = axs[i]
        scatter = ax.scatter(sqv0, sqv1, sqv2, s=0.5, c=mu, cmap="rainbow", rasterized=True)
        ax.view_init(elev=25.0, azim=-100)

        ax.set_xlabel(r"$FV0$", fontsize=9, labelpad=-12, rotation=0)
        ax.set_ylabel(r"$FV1$", fontsize=9, labelpad=-5, rotation=0)
        ax.set_zlabel(r"$FV2$", fontsize=9, labelpad=-7, rotation=0)
        # ax.tick_params(labelsize=7, pad=0)
        ax.tick_params("x", labelsize=7, pad=-7)
        ax.tick_params("y", labelsize=7, pad=-2)
        ax.tick_params("z", labelsize=7, pad=-2)

        # ax.set_title(features_tex[i])
        cbar = fig.colorbar(scatter, ax=ax, fraction=0.025, pad=-0.02)  # , location="top", orientation='horizontal')
        cbar.ax.tick_params(labelsize=7, pad=0)

        # cbar.ax.yaxis.set_major_locator(plt.MultipleLocator(cbar_major_locator[i]))
        # cbar.ax.yaxis.set_minor_locator(plt.MultipleLocator(cbar_major_locator[i]*0.5))

        # cbar.ax.yaxis.set_major_locator(matplotlib.ticker.AutoLocator())
        # cbar.ax.set_axes_locator(plt.MultipleLocator(cbar_major_locator[i]))
        # cbar.ax.xaxis.set_minor_locator(plt.MultipleLocator(cbar_major_locator[i]*0.5))
        cbar.ax.set_title(feature_tex, fontsize=9)

        # cbar = .colorbar(axs.collections[0])
        # cbar.set_label(mu, fontsize=9)

        # ax.xaxis.set_major_locator(plt.MultipleLocator(0.5))
        # ax.xaxis.set_minor_locator(plt.MultipleLocator(0.25))
        # ax.yaxis.set_major_locator(plt.MultipleLocator(0.5))
        # ax.yaxis.set_minor_locator(plt.MultipleLocator(0.25))
        # ax.zaxis.set_major_locator(plt.MultipleLocator(0.4))
        # ax.zaxis.set_minor_locator(plt.MultipleLocator(0.2))

        ax.grid(True, which="minor")

    annotation = [r"$(a)$", r"$(b)$", r"$(c)$", r"$(d)$", r"$(e)$", r"$(f)$"]
    for i in range(len(annotation)):
        axs[i].text2D(0.7, 0.8, annotation[i], transform=axs[i].transAxes, fontsize=9)

    # plt.tight_layout( h_pad=0.2, w_pad=1.7) #, h_pad=-3, w_pad=2)
    plt.tight_layout(pad=0.5, w_pad=1.2)
    plt.savefig("figures/SVD_feature.pdf", format="pdf", dpi=300)
    plt.savefig("figures/SVD_feature.png", format="png", dpi=300)  # png figure for slides making
    plt.show()
    plt.close()


def get_LML_date(folder, feature):
    data = np.loadtxt(f"{folder}/data_{feature}_LML.txt", skiprows=1, delimiter=",", unpack=True)
    gp_theta0, gp_theta1, theta0, theta1, LML = data[0], data[1], data[2], data[3], data[4:]
    return gp_theta0, gp_theta1, theta0, theta1, LML


def plot_LML_contour(tex_lw=240.71031, ppi=72):
    #folder = "../data/20240916_random"
    folder = "../data/20240918_random" # the data used fo submission arxiv
    #folder = "../data/20240924_random"
    feature_to_tex = {"kappa": r"$\kappa$", "f": r"$f$", "gL": r"$\gamma L$", "R2": r"$R^2/L^2$", "Rg2": r"$R_g^2/L$", "Sxz": r"$R_{xz}/R_g^2$"}
    fig = plt.figure(figsize=(tex_lw / ppi * 1, tex_lw / ppi * 1.15))
    plt.rc("text", usetex=True)
    plt.rc("text.latex", preamble=r"\usepackage{physics}")
    # kappa, f, gL
    ax11 = fig.add_subplot(321)
    ax12 = fig.add_subplot(322)
    ax13 = fig.add_subplot(323)
    # R2, Rg2, Sxz
    ax21 = fig.add_subplot(324)
    ax22 = fig.add_subplot(325)
    ax23 = fig.add_subplot(326)

    axs = [ax11, ax12, ax13, ax21, ax22, ax23]
    grid_size = 2
    ticks = [
        (np.logspace(-1, 0, grid_size), np.logspace(-3, -2, grid_size)),
        (np.logspace(-1, 0, grid_size), np.logspace(-3, -2, grid_size)),
        (np.logspace(-1, 0, grid_size), np.logspace(-3, -2, grid_size)),
        (np.logspace(0, 0.5, grid_size), np.logspace(-5, -4, grid_size)),
        (np.logspace(0.2, 0.4, grid_size), np.logspace(-7, -5, grid_size)),
        (np.logspace(0.3, 0.6, grid_size), np.logspace(-7, -5, grid_size)),
    ]
    ticklabels = [
        ((r"$10^{-1}$", r"$10^{0}$"), (r"$10^{-3}$", r"$10^{-2}$")),
        ((r"$10^{-1}$", r"$10^{0}$"), (r"$10^{-3}$", r"$10^{-2}$")),
        ((r"$10^{-1}$", r"$10^{0}$"), (r"$10^{-3}$", r"$10^{-2}$")),
        ((r"$10^{0}$", r"$10^{0.5}$"), (r"$10^{-5}$", r"$10^{-4}$")),
        ((r"$10^{0.2}$", r"$10^{0.4}$"), (r"$10^{-7}$", r"$10^{-5}$")),
        ((r"$10^{0.3}$", r"$10^{0.6}$"), (r"$10^{-7}$", r"$10^{-5}$")),
    ]

    # "kappa": (np.logspace(-1, 0, grid_size), np.logspace(-3, -2, grid_size)),
    # "f": (np.logspace(-1, 0, grid_size), np.logspace(-3, -2, grid_size)),
    # "gL": (np.logspace(-1, 0, grid_size), np.logspace(-3, -2, grid_size)),
    # "R2": (np.logspace(0, 0.5, grid_size), np.logspace(-5, -4, grid_size)),
    # "Rg2": (np.logspace(0.2, 0.5, grid_size), np.logspace(-7, -6, grid_size)),
    # "Sxz": (np.logspace(0.3, 0.6, grid_size), np.logspace(-7, -6, grid_size)),
    i = -1
    for feature_name, feature_tex in feature_to_tex.items():
        i += 1
        gp_theta0, gp_theta1, theta0, theta1, LML = get_LML_date(folder, feature_name)
        Theta0, Theta1 = np.meshgrid(theta0, theta1)
        axs[i].contour(Theta0, Theta1, LML, levels=50)  # , cmap="summer")
        axs[i].plot([gp_theta0[0]], [gp_theta1[0]], "x", color="black", markersize=5, markeredgewidth=1)  # , label=r"l=%.2e, $\sigma$=%.2e" % (gp_theta0[0], gp_theta1[0]))
        axs[i].set_xscale("log")
        axs[i].set_yscale("log")
        axs[i].tick_params(which="both", direction="in", top="on", right="on", labelbottom=True, labelleft=True, labelsize=7)
        # axs[i].set_xlabel(r'theta0: l')
        # axs[i].set_ylabel(r'theta1: $\sigma$')
        axs[i].set_xticks(ticks[i][0])
        axs[i].set_yticks(ticks[i][1])
        axs[i].set_xticklabels(ticklabels[i][0])
        axs[i].set_yticklabels(ticklabels[i][1])
        axs[i].minorticks_off()

        """
        yticks = axs[i].yaxis.get_minor_ticks()
        for tick in yticks:
            if np.log10(tick.get_loc())%1 != 0:
                tick.label1.set_visible(False)
        xticks = axs[i].xaxis.get_minor_ticks()
        for tick in xticks:
            if np.log10(tick.get_loc())%1 != 0:
                tick.label1.set_visible(False)
        """
        axs[i].legend(title=feature_tex, fontsize=9, frameon=False)

    axall = fig.add_subplot(111, frameon=False)
    axall.tick_params(labelcolor="none", which="both", top=False, bottom=False, left=False, right=False)
    axall.set_xlabel(r"$l$", fontsize=9, labelpad=-3)
    axall.set_ylabel(r"$\sigma$", fontsize=9, labelpad=-3)
    # ax13.set_ylabel("ML Inversion", fontsize=9, labelpad=0)

    annotation = [r"$(a)$", r"$(b)$", r"$(c)$", r"$(d)$", r"$(e)$", r"$(f)$"]
    for i in range(len(annotation)):
        axs[i].text(0.8, 0.15, annotation[i], transform=axs[i].transAxes, fontsize=9)

    plt.tight_layout(pad=0)
    plt.subplots_adjust(left=0.12, bottom=0.08)

    plt.savefig("figures/LML_contour.pdf", dpi=300)
    plt.savefig("figures/LML_contour.png", dpi=300)
    plt.show()
    plt.close()


def get_GRP_data(folder, feature):
    data = np.loadtxt(f"{folder}/data_{feature}_prediction.txt", skiprows=1, delimiter=",", unpack=True)
    return data


def plot_GPR_prediction(tex_lw=240.71031, ppi=72):
    folder = "../data/20240918_random"
    #folder = "../data/20240924_random"

    all_feature, sqv0, sqv1, sqv2 = read_SVD_data(folder)
    all_feature_name = ["kappa", "f", "gL", "R2", "Rg2", "Sxx", "Syy", "Szz", "Sxy", "Sxz", "Syz"]
    feature_to_tex = {"kappa": r"$\kappa$", "f": r"$f$", "gL": r"$\gamma L$", "R2": r"$R^2/L^2$", "Rg2": r"$R_g^2/L$", "Sxz": r"$R_{xz}/R_g^2$"}

    fig = plt.figure(figsize=(tex_lw / ppi * 1, tex_lw / ppi * 1.15))
    plt.rc("text", usetex=True)
    plt.rc("text.latex", preamble=r"\usepackage{physics}")

    # kappa, f, gL
    ax11 = fig.add_subplot(321)
    ax12 = fig.add_subplot(322)
    ax13 = fig.add_subplot(323)
    # R2, Rg2, Sxz
    ax21 = fig.add_subplot(324)
    ax22 = fig.add_subplot(325)
    ax23 = fig.add_subplot(326)

    axs = [ax11, ax12, ax13, ax21, ax22, ax23]
    major_locator = [5, 0.2, 0.5, 0.2, 5, 0.2]
    minor_locator = [2.5, 0.1, 0.25, 0.1, 2.5, 0.1]
    i = -1
    for feature_name, feature_tex in feature_to_tex.items():
        i += 1
        mu, mu_pred, mu_err = get_GRP_data(folder, feature_name)
        L = 200
        if feature_name == "R2":
            mu /= L * L
            mu_pred /= L * L
        if feature_name == "Rg2":
            mu /= L
            mu_pred /= L
        relative_err = np.abs(mu - mu_pred) / mu
        norm = plt.Normalize(0, 1)
        #axs[i].scatter(mu, mu_pred, s=0.5, marker=".", c=relative_err, cmap="rainbow", norm=norm)
        axs[i].scatter(mu, mu_pred, s=0.5, marker=".", c="royalblue")
        axs[i].plot(mu, mu, color="gray", linestyle="--", lw=0.25, alpha=0.5)
        axs[i].xaxis.set_major_locator(plt.MultipleLocator(major_locator[i]))
        axs[i].xaxis.set_minor_locator(plt.MultipleLocator(minor_locator[i]))
        axs[i].yaxis.set_major_locator(plt.MultipleLocator(major_locator[i]))
        axs[i].yaxis.set_minor_locator(plt.MultipleLocator(minor_locator[i]))
        axs[i].grid(True, which="both", linestyle="--", linewidth=0.5)
        axs[i].tick_params(which="both", direction="in", top="on", right="on", labelbottom=True, labelleft=True, labelsize=7)
        # axs[i].legend(title = feature_tex, fontsize=9, loc="upper left")
        axs[i].text(0.2, 0.7, feature_tex, transform=axs[i].transAxes, fontsize=9)
        r_squ = np.corrcoef(mu, mu_pred)[0, 1]**2
        axs[i].text(0.4, 0.1, r"$r^2$=%.5f" % r_squ, transform=axs[i].transAxes, fontsize=9)

    axall = fig.add_subplot(111, frameon=False)
    axall.tick_params(labelcolor="none", which="both", top=False, bottom=False, left=False, right=False)
    axall.set_xlabel("MC References", fontsize=9, labelpad=-3)
    axall.set_ylabel("ML Inversion", fontsize=9, labelpad=-3)
    # ax13.set_ylabel("ML Inversion", fontsize=9, labelpad=0)

    annotation = [r"$(a)$", r"$(b)$", r"$(c)$", r"$(d)$", r"$(e)$", r"$(f)$"]
    for i in range(len(annotation)):
        axs[i].text(0.7, 0.3, annotation[i], transform=axs[i].transAxes, fontsize=9)

    plt.tight_layout(pad=0)
    plt.subplots_adjust(left=0.12, bottom=0.08)
    plt.savefig("figures/GPR_prediction.pdf", format="pdf", dpi=300)
    plt.savefig("figures/GPR_prediction.png", format="png", dpi=300)
    plt.show()
    plt.close()
