#!/opt/homebrew/bin/python3
from config_plot import *
from obs_plot import *
from Sq_plot import *
from ML_plot import *
def main():

    # for model paper
    #plot_config_update_demo()
    #plot_obs_kappa()
    #plot_obs_f()
    #plot_obs_gamma()

    #plot_obs_f_3panel()
    #plot_obs_gamma_3panel()
    #plot_obs_f_3panel_sq()

    #plot_Sq()
    #plot_Sq2D()
    #plot_Sq2D_f_conf()
    #plot_Sq2D_g()
    #plot_TOC()

    # for ML paper
    #plot_Sq2D_kappa()

    # need some fix since I lost some code that failed to save
    plot_config_fg()
    #

    #plot_Sq2D_fg()

    #plot_SVD_data()
    #plot_SVD_feature_data()
    #plot_LML_contour()
    #plot_GPR_prediction()


# todo: 1. increas hspace between subplot
# 2. remove legend for theory line
# 3. use Q instead of q
# 4. add kappa=10 in Sq plot
# 5. wor on TOC, use exact size as instructed



if __name__ == '__main__':
    main()