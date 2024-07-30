import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches


def ax_plot_config(ax, start, filename, color="gray"):
    data = np.genfromtxt(filename, delimiter=',', skip_header=1)
    x, y, z = data[start:, 0], data[start:, 1], data[start:, 2]
    # tx, ty, tz = data[:, 3], data[:, 4], data[:, 5]
    # ax.plot(x, y, z, "-o", color="gray", mfc="black", mec="black")
    ax.plot(x, y, z, "-o", color=color, mfc="gray", ms=3)

    ax.set_xlim(np.min(x)-1, np.max(x)+1)
    ax.set_ylim(np.min(y)-1, np.max(y)+1)
    ax.set_zlim(np.min(z)-1, np.max(z)+1)

    # ax.set_ylim(-ymax-1, ymax+1)
    # ax.set_zlim(-zmax-1, zmax+1)

    ax.legend()
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')


def ax_plot_config_xz(ax, filename, color="black"):
    data = np.genfromtxt(filename, delimiter=',', skip_header=1)
    x, y, z = data[:, 0], data[:, 1], data[:, 2]
    # tx, ty, tz = data[:, 3], data[:, 4], data[:, 5]
    # ax.plot(x, y, z, "-o", color="gray", mfc="black", mec="black")
    ax.plot(x, z, "-o", color=color, mfc="black", ms=3)

    ax.set_xlim(np.min(x)-1, np.max(x)+1)
    # ax.set_ylim(np.min(y)-1, np.max(y)+1)
    ax.set_ylim(np.min(z)-1, np.max(z)+1)

    # ax.set_ylim(-ymax-1, ymax+1)
    # ax.set_zlim(-zmax-1, zmax+1)


def ax_plot_con_rot_update(ax):  # , start, filename_pre, filename_post):
    dx = np.array([0.8, 0.7, 0.9, 0.8, 0.65, 0.75, 0.85])
    dzmax = np.sqrt(1.0 - dx**2)
    dz = [0.8, 0.9, -0.65, 0.7, -0.5, -0.6, 0.7]
    dz = np.array(dz)*dzmax
    x = [0.0]
    z = [0.0]
    for i in range(len(dx)):
        x.append(x[-1]+dx[i])
        z.append(z[-1]+dz[i])
    i, j = 1, 6
    z1 = []
    for k in range(len(x)):
        if k > i and k < j:
            z1.append(z[i]-0.6*(z[k]-z[i]))
        else:
            z1.append(z[k])
    # z1 = [1, 1, 0.6, 0.7, 0.6, 0.7, 1.0, 1.
    # z1 = np.array(z1)*np.array(z)
    ms = 15
    ax.plot(x[i:j+1], z1[i:j+1], "o:", color="gray", mfc="white", mec="gray", ms=ms)
    ax.plot(x, z, "o-", color="black", mfc="white", mec="black", ms=ms)

    ax.plot([x[i], x[j]], [z[i], z[j]], "o", color="red", mfc="none", mec="red", ms=ms)
    xcorr = -0.05
    zcorr = -0.1
    ax.annotate(r"$i$", xy=(x[i]+xcorr, z[i]+zcorr), fontsize=9)
    ax.annotate(r"$j$", xy=(x[j]+xcorr, z[j]+zcorr), fontsize=9)

    # a3 = mpatches.FancyArrowPatch((0.5*(x[i]+x[i+1]), 0.5*(z[i]+z[i+1])), (0.5*(x[i]+x[i+1]), 0.5*(z1[i]+z1[i+1])), arrowstyle="Fancy", connectionstyle="arc3,rad=-.5")
    # ax.add_patch(a3)

    # ax.set_xlim(-0.5,6)
    # ax.set_ylim(-1.5,2)

    # ax.view_init(elev=30., azim=-140)
    # ax.set_aspect('equal')
    # ax.set_axis_off()


def ax_plot_tan_rot_update(ax, yshift=-4):
    dx = np.array([0.8, 0.7, 0.9, 0.8, 0.65, 0.75, 0.85])
    dzmax = np.sqrt(1.0 - dx**2)
    dz = [0.8, 0.9, -0.65, 0.7, -0.5, -0.6, 0.7]
    dz = np.array(dz)*dzmax
    x = [0.0]
    z = [0.0]
    for i in range(len(dx)):
        x.append(x[-1]+dx[i])
        z.append(z[-1]+dz[i])
    cos_theta = np.cos(-np.pi/6)
    sin_theta = np.sin(-np.pi/6)
    print("cos_theta, sin_theta", cos_theta, sin_theta)
    k = 3
    x1 = x[:k+1]
    z1 = z[:k+1]
    for i in range(len(x))[k+1:]:
        x1.append((x[i]-x[k])*cos_theta - (z[i]-z[k])*sin_theta + x[k])
        z1.append((x[i]-x[k])*sin_theta + (z[i]-z[k])*cos_theta + z[k])
    z = np.array(z)+yshift
    z1 = np.array(z1)+yshift

    ms = 15
    ax.plot(x1[k:], z1[k:], "o:", color="gray", mfc="white", mec="gray", ms=ms)
    ax.plot(x, z, "o-", color="black", mfc="white", mec="black", ms=ms)

    ax.plot([x[k]], [z[k]], "o", color="red", mfc="none", mec="red", ms=ms)

    xcorr = -0.05
    zcorr = -0.1
    ax.annotate(r"$k$", xy=(x[k]+xcorr, z[k]+zcorr), fontsize=9)
    ax.annotate(r"$l$", xy=(x[k+1]+xcorr, z[k+1]+zcorr), fontsize=9)

    # ax.set_xlim(-0.5,6)
    # ax.set_ylim(-1.5,2)
    # ax.axis("equal")
    # ax.set_axis_off()


def plot_config_update_demo(tex_lw=240.71031, ppi=72):

    fig = plt.figure(figsize=(tex_lw / ppi * 1, tex_lw / ppi * 0.9))
    # plt.rc("text", usetex=True)
    # plt.rc("text.latex", preamble=r"\usepackage{physics}")
    # for outofplane twist (2D CANAL)
    ax = fig.add_subplot(111)  # , projection='3d')
    # ax2 = fig.add_subplot(212)  # , projection='3d')

    # con_rot_file_pre = "config_data/config_con_rot_pre.csv"
    # con_rot_file_post = "config_data/config_con_rot_post.csv"
    ax_plot_con_rot_update(ax)  # , 35, con_rot_file_pre, con_rot_file_post)

    ax_plot_tan_rot_update(ax, -2.5)
    ax.set_xlim(-0.5, 6)
    ax.axis("equal")
    ax.set_axis_off()

    # axs = [ax1, ax2]
    bbox = ax.get_position()
    # annotation
    x = bbox.x0 + bbox.width*(0.0)
    y = bbox.y0 + bbox.height*(0.9)
    fig.text(x, y, r"$(a)$", ha='center', va='center', fontsize=9)
    x = bbox.x0 + bbox.width*(0.0)
    y = bbox.y0 + bbox.height*(0.4)
    fig.text(x, y, r"$(b)$", ha='center', va='center', fontsize=9)

    plt.tight_layout()
    # plt.tight_layout(pad=0.1, h_pad=-1.4, w_pad=-1)
    # plt.tight_layout(h_pad=-2, w_pad=-6)
    # plt.show()
    plt.savefig("figures/config_update_demo.pdf", format="pdf")
    plt.close()


