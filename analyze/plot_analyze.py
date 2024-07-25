import matplotlib.pyplot as plt
import numpy as np
from matplotlib.cm import ScalarMappable
from matplotlib import colormaps as cm
from matplotlib.colors import Normalize


def plot_polymer_config(filename, finfo, show=False):
    data = np.genfromtxt(filename, delimiter=',', skip_header=1)
    x, y, z = data[:, 0], data[:, 1], data[:, 2]
    tx, ty, tz = data[:, 3], data[:, 4], data[:, 5]

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    #ax.plot(x, y, z, "-o", color="gray", mfc="black", mec="black")
    ax.plot(x, y, z, "-", color="gray")
    ax.plot([x[0]], [y[0]], [z[0]], "o", color="red", label="start")
    ax.plot([x[-1]], [y[-1]], [z[-1]], "o", color="blue", label="end")

    #ymax = max(np.max(y), -np.min(y))
    #zmax = max(np.max(z), -np.min(z))
    ax.set_xlim(np.min(x)-1, np.max(x)+1)
    ax.set_ylim(np.min(y)-1, np.max(y)+1)
    ax.set_zlim(np.min(z)-1, np.max(z)+1)

    #ax.set_ylim(-ymax-1, ymax+1)
    #ax.set_zlim(-zmax-1, zmax+1)

    ax.legend()
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title(finfo)
    '''
    cmap = cm.get_cmap("jet_r")
    norm = Normalize(vmin=0, vmax=0.5*np.pi)
    deg = np.arccos(np.abs(n[:, 2]))
    for i in range(len(r)):
        # ax.plot(np.append(beadsrx[1:], beadsrx[1]), np.append(beadsry[1:], beadsry[1]), np.append(beadsrz[1:], beadsrz[1]), "-", color="gray")
        pass
    cbar = plt.colorbar(ScalarMappable(norm=Normalize(vmin=0, vmax=0.5*np.pi), cmap=cm.get_cmap("jet_r")), ax=ax, ticks=[0, np.pi/6, np.pi/3, np.pi/2])
    cbar.ax.set_yticklabels([r"$0$", r"$\pi/6$", r"$\pi/3$", r"$\pi/2$"])
    cbar.ax.tick_params(direction="in")
    cbar.ax.set_title(r"$\arccos{|\hat{\mathrm{n}}\cdot\hat{\mathrm{z}}|}$")

    cbar.set_label('Z')
    '''
    '''
    ax.xaxis.set_major_locator(plt.MultipleLocator(1))
    ax.xaxis.set_minor_locator(plt.MultipleLocator(0.5))
    ax.yaxis.set_major_locator(plt.MultipleLocator(1))
    ax.yaxis.set_minor_locator(plt.MultipleLocator(0.5))
    ax.zaxis.set_major_locator(plt.MultipleLocator(1))
    ax.zaxis.set_minor_locator(plt.MultipleLocator(0.5))
    '''
    ax.set_aspect("equal")
    plt.savefig(filename.replace(".csv", ".png"))
    if show:
        plt.show()
    plt.close()


def plot_MC_step(filename, finfo):
    data = np.genfromtxt(filename, delimiter=',', skip_header=4)
    E, Tb, X, Y, Z, R = data[:, 0], data[:, 1], data[:, 2], data[:, 3], data[:, 4], data[:, 5]

    fig, axs = plt.subplots(7, 1, figsize=(3, 15))

    axs[1].plot(range(len(Tb)), Tb, "*")
    axs[1].set_xlabel("MC sweep (1000L step per sweep)")
    axs[1].set_ylabel("Tb")

    axs[2].plot(range(len(X)), X, "*")
    axs[2].set_xlabel("MC sweep (1000L step per sweep)")
    axs[2].set_ylabel("X")

    axs[3].plot(range(len(Y)), Y, "*")
    axs[3].set_xlabel("MC sweep (1000L step per sweep)")
    axs[3].set_ylabel("Y")

    axs[4].plot(range(len(Z)), Z, "*")
    axs[4].set_xlabel("MC sweep (1000L step per sweep)")
    axs[4].set_ylabel("Z")

    axs[5].plot(range(len(R)), R, "*")
    axs[5].set_xlabel("MC sweep (1000L step per sweep)")
    axs[5].set_ylabel("R")

    axs[6].plot(range(len(E)), E, "*")
    axs[6].set_xlabel("MC sweep (1000L step per sweep)")
    axs[6].set_ylabel("Energy")
    plt.tight_layout()
    plt.savefig(filename.replace(".csv", ".png"))
    plt.close()


def plot_obs(folder, finfos, parameters, xparam):
    fig, axs = plt.subplots(1, 3, figsize=(12, 3))

    xpar = []
    all_X, all_X_err = [], []
    all_Y, all_Y_err = [], []
    all_Z, all_Z_err = [], []
    all_R, all_R_err = [], []

    all_Sq = []
    all_qB = []
    all_tts = []
    all_spB = []

    for i in range(len(finfos)):
        finfo = finfos[i]
        L, kappa, f, g = parameters[i]
        filename = f"{folder}/obs_{finfo}.csv"
        data = np.genfromtxt(filename, delimiter=',', skip_header=1)
        E, Tb, X, Y, Z, R = data[0, 5], data[0, 6], data[0, 7], data[0, 8], data[0, 9], data[0, 10]
        E_err, Tb_err, X_err, Y_err, Z_err, R_err = data[1, 5], data[1, 6], data[1, 7], data[1, 8], data[1, 9], data[1, 10]
        Sq = data[0, 11:]
        qB = data[2, 11:]
        tts = data[3, 11:]
        spB = data[5, 11:]

        if xparam == "L":
            xpar.append(L)
        elif xparam == "kappa":
            xpar.append(kappa)
        elif xparam == "f":
            xpar.append(f)
        elif xparam == "g":
            xpar.append(g)
        else:
            raise ValueError(f"Unknown xparam={xparam}")

        all_X.append(X)
        all_X_err.append(X_err)
        all_Y.append(Y)
        all_Y_err.append(Y_err)
        all_Z.append(Z)
        all_Z_err.append(Z_err)
        all_R.append(R)
        all_R_err.append(R_err)

        all_Sq.append(Sq)
        all_qB.append(qB)
        all_tts.append(tts)
        all_spB.append(spB)

    axs[0].errorbar(xpar, all_X, yerr=all_X_err, marker="o", label="X")
    axs[0].errorbar(xpar, all_Y, yerr=all_Y_err, marker="x", label="Y")
    axs[0].errorbar(xpar, all_Z, yerr=all_Z_err, marker="s", label="Z")
    axs[0].errorbar(xpar, all_R, yerr=all_R_err, marker="+", label="R")

    axs[0].set_xlabel(xparam)
    axs[0].set_ylabel("end-end distance")
    axs[0].legend()

    for i in range(len(all_Sq)):
        axs[1].loglog(all_qB[i], all_Sq[i], "-", label=f"{xparam}={xpar[i]}")
    axs[1].set_xlabel("qB")
    axs[1].set_ylabel("S(qB)")
    axs[1].legend()

    for i in range(len(all_tts)):
        print("all_tts[i][1]", all_tts[i][1])
        print("-1/np.log(all_tts[i][1])", -1/np.log(all_tts[i][1]))
        lp = -1/np.log(all_tts[i][1])

        axs[2].semilogy(all_spB[i][:20], all_tts[i][:20], "-", label=f"{xparam}={xpar[i]}")
        #axs[2].plot(all_spB[i], all_tts[i], "-", label=f"{xparam}={xpar[i]}")
    axs[2].set_xlabel(r"$s/B$")
    axs[2].set_ylabel(r"$<\cos{\theta}(s)>$")
    axs[2].legend()
    plt.tight_layout()
    plt.savefig(f"{folder}/obs_{xparam}.png")
    plt.close()