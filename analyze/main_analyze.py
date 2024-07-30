#!/opt/homebrew/bin/python3
from plot_analyze import *
import numpy as np


def main():
    print("hello!")
    folder = "../data/scratch_local/20240725"
    # plot_polymer_config(folder+"/config_L50_kappa0.0_f0.0_g0.0.csv", "", True)
    # plot_MC_step(folder+"/obs_MC_L50_kappa0.0_f0.0_g0.0.csv", "", True)
    # plot_polymer_config(folder+"/config_L50_kappa10.0_f0.0_g0.0.csv", "", True)
    # plot_MC_step(folder+"/obs_MC_L50_kappa10.0_f0.0_g0.0.csv", "", True)
    # return 0
    L = 50
    kappa = 10.0
    f = 0.0
    g = 0.0
    parameters = []
    finfos = []
    for kappa in [1.0, 2.0, 4.0, 8.0, 10.0]:
        finfo = f"L{L}_kappa{kappa:.1f}_f{f:.1f}_g{g:.1f}"
        parameters.append((L, kappa, f, g))
        finfos.append(finfo)
        # plot_polymer_config(folder+f"/config_{finfo}.csv", finfo, True)
        # plot_polymer_config(folder+f"/config_{finfo}.csv", finfo)
        # plot_MC_step(folder+f"/obs_MC_{finfo}.csv", finfo)

    plot_obs(folder, finfos, parameters, "kappa")
    # return 0

    kappa = 10.0
    g = 0.0
    parameters = []
    finfos = []
    for f in [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]:
        finfo = f"L{L}_kappa{kappa:.1f}_f{f:.1f}_g{g:.1f}"
        parameters.append((L, kappa, f, g))
        finfos.append(finfo)
        # plot_polymer_config(folder+f"/config_{finfo}.csv", finfo, True)
        # plot_polymer_config(folder+f"/config_{finfo}.csv", finfo)
        plot_MC_step(folder+f"/obs_MC_{finfo}.csv", finfo)
    plot_obs(folder, finfos, parameters, "f")
    # return 0

    kappa = 10.0
    f = 0.0
    parameters = []
    finfos = []
    for g in [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]:
        finfo = f"L{L}_kappa{kappa:.1f}_f{f:.1f}_g{g:.1f}"
        parameters.append((L, kappa, f, g))
        finfos.append(finfo)
        # plot_polymer_config(folder+f"/config_{finfo}.csv", finfo, True)
        plot_polymer_config(folder+f"/config_{finfo}.csv", finfo)
        plot_MC_step(folder+f"/obs_MC_{finfo}.csv", finfo)
    plot_obs(folder, finfos, parameters, "g")


if __name__ == "__main__":
    main()
