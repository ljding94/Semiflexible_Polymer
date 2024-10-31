#!/opt/homebrew/bin/python3
from plot_analyze import *
import numpy as np
from ML_analyze import *
import sys
import random
import time
from KAN_analyze import *

def main():

    print("analyzing data using ML model")
    folder = "../data/20240924_random"
    rand_num = 5500
    rand_max = 1600
    L = 200
    parameters = []
    for i in range(rand_num):
        filename = f"{folder}/obs_L{L}_random_run{i}.csv"
        if os.path.exists(filename):
            parameters.append([L, i])
        if len(parameters) >= rand_max:
            break
    #print("parameters", parameters)
    print("total number of parameters", len(parameters))

    get_data_for_KAN(folder, parameters)


if __name__ == "__main__":
    start_time = time.time()
    main()
    end_time = time.time()
    execution_time = end_time - start_time
    print(f"Execution time: {execution_time} seconds")
