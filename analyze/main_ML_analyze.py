#!/opt/homebrew/bin/python3
from plot_analyze import *
import numpy as np
from ML_analyze import *
import sys
import random
import time


def main():

    print("analyzing data using ML model")
    folder = "../data/20240913_random"
    rand_max = 96*6
    L = 200
    parameters = [[L, rand_num] for rand_num in range(rand_max)]

    print("parameters", parameters)
    print("total number of parameters", len(parameters))

    calc_svd(folder, parameters)
    # plot_pddf_acf(folder, parameters, max_z=15, n_bin=100)

    random.shuffle(parameters)
    parameters_train = parameters[:int(0.7*len(parameters))]
    parameters_test = parameters[int(0.7*len(parameters)):]

    # all_feature_mean, all_feature_std, all_gp_per_feature = GaussianProcess_optimization(folder, parameters_train, all_feature_names)
    # GaussianProcess_prediction(folder, parameters_test, all_feature_names, all_feature_mean, all_feature_std, all_gp_per_feature)

    calc_Sq_fitted_Rg2(folder, parameters_test, all_feature_names)


if __name__ == "__main__":
    start_time = time.time()
    main()
    end_time = time.time()
    execution_time = end_time - start_time
    print(f"Execution time: {execution_time} seconds")
