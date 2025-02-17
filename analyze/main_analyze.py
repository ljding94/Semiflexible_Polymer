#!/opt/homebrew/bin/python3
from plot_analyze import *
import numpy as np


def main():
    print("hello!")
    #plot_polymer_config("../data/20240828/config_L200_kappa10.0_f0.00_gL0.00.csv", "", True)
    #plot_R_distribution("../data/20240828/obs_MC_L200_kappa10.0_f0.00_gL0.00.csv")
    #plot_polymer_config("../data/scratch_local/20240829/config_L10_kappa10.0_f0.0_gL0.0.csv", "", True)
    #plot_R_distribution("../data/scratch_local/20240830/20240830/obs_MC_L200_kappa10.0_f0.00_gL0.00.csv")
    #plot_multi_config("../data/scratch_local/20240829/L2_kappa1.0_f0.0_gL0.0")
    #return 0
    #folder = "../data/scratch_local/20240725"
    # plot_polymer_config(folder+"/config_L50_kappa0.0_f0.0_g0.0.csv", "", True)
    # plot_MC_step(folder+"/obs_MC_L50_kappa0.0_f0.0_g0.0.csv", "", True)
    #plot_polymer_config(folder+"/config_L50_kappa10.0_f0.0_g0.0.csv", "", True)
    # plot_MC_step(folder+"/obs_MC_L50_kappa10.0_f0.0_g0.0.csv", "", True)
    # return 0

    folder = "../data/20240821"
    L = 200
    kappa = 10.0
    gL = 0.0
    parameters = []
    finfos = []
    fs = np.arange(0.00, 0.301, 0.02)
    for f in fs:
        finfo = f"L{L}_kappa{kappa:.1f}_f{f:.2f}_gL{gL:.2f}"
        parameters.append((L, kappa, f, gL))
        finfos.append(finfo)
        #plot_polymer_config(folder+f"/config_{finfo}.csv", finfo, True)
        plot_polymer_config(folder+f"/config_{finfo}.csv", finfo)
        #plot_MC_step(folder+f"/obs_MC_{finfo}.csv", finfo)
    plot_obs(folder, finfos, parameters, "f")
    #return 0

    folder = "../data/20240820"
    L = 200
    kappa = 10.0
    f = 0.0
    parameters = []
    finfos = []
    gLs = np.arange(0.00, 1.501, 0.10)
    for gL in gLs:
        finfo = f"L{L}_kappa{kappa:.1f}_f{f:.2f}_gL{gL:.2f}"
        parameters.append((L, kappa, f, gL))
        finfos.append(finfo)
        # plot_polymer_config(folder+f"/config_{finfo}.csv", finfo, True)
        plot_polymer_config(folder+f"/config_{finfo}.csv", finfo)
        #plot_MC_step(folder+f"/obs_MC_{finfo}.csv", finfo)
    #plot_obs(folder, finfos, parameters, "g")


if __name__ == "__main__":
    main()
