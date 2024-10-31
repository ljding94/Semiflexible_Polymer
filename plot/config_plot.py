import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches
import random


def ax_plot_config(ax, folder, param, zlim_shift, label):
    L, kappa, f, gL = param
    finfo = f"L{L}_kappa{kappa:.1f}_f{f:.2f}_gL{gL:.2f}"
    filename = f"{folder}/config_{finfo}.csv"
    data = np.genfromtxt(filename, delimiter=',', skip_header=1)
    x, y, z = data[:, 0], data[:, 1], data[:, 2]
    # tx, ty, tz = data[:, 3], data[:, 4], data[:, 5]
    ax.plot(x, y, z, "-", lw=2, label=label)
    # y_min = np.min(y)
    # y_max = np.max(y)

    # for i in range(len(x)-1):
    # alpha_xz = 0.8*(y[i]-y_min+0.01)/(y_max-y_min+0.01)+0.1
    # ax.plot([x[i], x[i+1]], [z[i], z[i+1]], "-", alpha=alpha_xz, color=color, lw=1, label=label)

    # ax.set_xlim(np.min(x)-1, np.max(x)+1)
    ax.set_ylim(np.min(y)-1, np.max(y)+1)
    ax.set_zlim(np.min(z)-1+zlim_shift, np.max(z)+1)
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$y$")
    ax.set_zlabel(r"$z$")
    ax.set_aspect("equal")


def ax_plot_2dconfig_from_file(ax, filename, color, xlim, randflip=1, flip=0, xshift=0, yshift=0):
    data = np.genfromtxt(filename, delimiter=',', skip_header=1)
    x, y, z = data[:, 0], data[:, 1], data[:, 2]
    '''
    if (-xlim > np.min(x) or xlim < np.max(x)):
        print("np.min(x), np.max(x)", np.min(x), np.max(x))
    if (-xlim > np.min(y) or xlim < np.max(y)):
        print("np.min(y), np.max(y)", np.min(y), np.max(y))
    if (-xlim > np.min(z) or xlim < np.max(z)):
        print("np.min(z), np.max(z)", np.min(z), np.max(z))
    '''
    if (randflip and random.random() < 0.5):
        x, z = -x, -z
    if(flip and (z[-1]<0 or x[-1]<0)):
        x, z = -x, -z

    ymin = np.min(y)
    ymax = np.max(y)
    # ax.plot([x[0], x[1]], [z[0], z[1]], "-", alpha=0.8*(y[0]-ymin)/(ymax-ymin)+0.1, lw=1)
    # color = plt.gca().lines[-1].get_color()
    # color = (np.random.random(), np.random.random(), np.random.random())
    X, Y, Z = x[-1], y[-1], z[-1]
    R = np.sqrt(X**2 + Y**2 + Z**2)
    #color = (0.4*X/R+0.6, (0.4*Y/R+0.6)*1.0, 0.4*Z/R+0.6)
    color =  (np.abs(Z)/R, np.abs(Y)/R, np.abs(X)/R)

    #ax.plot(x+xshift, z+yshift+0.15, "-", color="white", lw=1, solid_capstyle='round') #, rasterized=True)
    #ax.plot(x+xshift, z+yshift-0.25, "-", color="gray", lw=1, solid_capstyle='round') #, rasterized=True)
    ax.plot(x+xshift, z+yshift, "-", color=color, lw=0.75, solid_capstyle='round') #, rasterized=True)
    #for i in range(len(x)-1):
    #    ax.plot([x[i]+xshift, x[i+1]+xshift], [z[i] + yshift, z[i+1]+yshift], "-", alpha=0.6*(y[i]-ymin)/(ymax-ymin)+0.2, color=color, lw=1, solid_capstyle='round', rasterized=True)
    # ax.set_xlim(-xlim, xlim)
    # ax.set_ylim(-xlim, xlim)
    # ax.set_zlim(-xlim, xlim)
    ax.set_aspect("equal")
    # ax.set_axis_off()


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
    # ax.annotate(r"$k'$", xy=(x[k+1]+xcorr, z[k+1]+zcorr), fontsize=9)

    # ax.set_xlim(-0.5,6)
    # ax.set_ylim(-1.5,2)
    # ax.axis("equal")
    # ax.set_axis_off()


def plot_config_update_demo(tex_lw=240.71031, ppi=72):

    fig = plt.figure(figsize=(tex_lw / ppi * 1, tex_lw / ppi * 0.75))
    plt.rc("text", usetex=True)
    plt.rc("text.latex", preamble=r"\usepackage{physics}")
    # for outofplane twist (2D CANAL)
    ax = fig.add_subplot(111)  # , projection='3d')
    # ax2 = fig.add_subplot(212)  # , projection='3d')

    # con_rot_file_pre = "config_data/config_con_rot_pre.csv"
    # con_rot_file_post = "config_data/config_con_rot_post.csv"
    ax_plot_con_rot_update(ax)  # , 35, con_rot_file_pre, con_rot_file_post)

    ax_plot_tan_rot_update(ax, -1.8)
    ax.set_xlim(-0.5, 6)
    ax.axis("equal")
    ax.set_axis_off()

    # axs = [ax1, ax2]
    bbox = ax.get_position()
    # annotation
    x = bbox.x0 + bbox.width*(0.0)
    y = bbox.y0 + bbox.height*(0.95)
    fig.text(x, y, r"$(a)$", ha='center', va='center', fontsize=9)
    x = bbox.x0 + bbox.width*(0.0)
    y = bbox.y0 + bbox.height*(0.4)
    fig.text(x, y, r"$(b)$", ha='center', va='center', fontsize=9)

    plt.tight_layout(pad=0.2)
    # plt.tight_layout(pad=0.1, h_pad=-1.4, w_pad=-1)
    # plt.tight_layout(h_pad=-2, w_pad=-6)

    plt.savefig("figures/config_update_demo.pdf", format="pdf")
    plt.show()
    plt.close()
