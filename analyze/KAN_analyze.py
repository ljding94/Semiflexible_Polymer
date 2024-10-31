import numpy as np
import matplotlib.pyplot as plt
import pickle
from kan import *
from ML_analyze import get_all_feature_Sq2D_data
import torch
from kan.utils import create_dataset_from_data


def get_data_for_KAN(folder, parameters):
    all_feature, all_feature_name, all_Sq2D_flatten, qB = get_all_feature_Sq2D_data(folder, parameters)
    qx, qy = np.meshgrid(qB, qB)
    qx = qx.flatten()
    qy = qy.flatten()
    qr = np.sqrt(qx**2 + qy**2)
    qtheta = np.arctan2(qy, qx)
    kappa, f, gamma = all_feature[:, 0], all_feature[:, 1], all_feature[:, 2]
    x_KAN = []
    y_KAN = []
    print()
    for i in range(len(kappa)):
        for j in range(len(qr)):
            #x_KAN.append([kappa[i], f[i], gamma[i], qr[j], qtheta[j]])
            x_KAN.append([kappa[i], f[i], gamma[i], qx[j], qy[j]])
            y_KAN.append(all_Sq2D_flatten[i, j])
    print("qr.shape, qtheta.shape, kappa.shape, f.shape, gamma.shape", qr.shape, qtheta.shape, kappa.shape, f.shape, gamma.shape)
    print("all_Sq2D_flatten.shape", all_Sq2D_flatten.shape)
    x_KAN = np.array(x_KAN)
    y_KAN = np.array(y_KAN)

    print("x_KAN.shape, y_KAN.shape", x_KAN.shape, y_KAN.shape)
    return x_KAN, y_KAN


def KAN_optimize(folder, parameters):
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(device)
    # get data for optimization/training
    x_kan, y_kan = get_data_for_KAN(folder, parameters)
    dataset = create_dataset_from_data(x_kan, y_kan, device=device)
    print("dataset['train_input'].shape, dataset['train_label'].shape", dataset['train_input'].shape, dataset['train_label'].shape)

    model = KAN(width=[5, 5, 1], grid=3, k=3, seed=42, device=device)

    # plot KAN at initialization
    model(dataset['train_input'])
    model.plot()

    # train the model
    model.fit(dataset, opt="LBFGS", steps=50, lamb=0.001)
