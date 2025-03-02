# some function to visualize the NN data
import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # Required for 3D plotting
import torch


def plot_loss(folder):
    # Load the loss data
    data = np.loadtxt(os.path.join(folder, "vae_losses.txt"), skiprows=1, delimiter=",")  # Skip the header row
    vae_train_loss = data[:, 0]
    vae_test_loss = data[:, 1]

    data = np.loadtxt(os.path.join(folder, "converterP2L_losses.txt"), skiprows=1, delimiter=",")
    converterP2L_train_loss = data[:, 0]
    converterP2L_test_loss = data[:, 1]

    data = np.loadtxt(os.path.join(folder, "converterL2F_losses.txt"), skiprows=1, delimiter=",")
    converterL2F_train_loss = data[:, 0]
    converterL2F_test_loss = data[:, 1]

    # Create a new figure
    plt.figure(figsize=(15, 6))
    ax1 = plt.subplot(131)
    ax2 = plt.subplot(132)
    ax3 = plt.subplot(133)
    # Plot the VAE losses
    ax1.loglog(vae_train_loss, label="VAE Train Loss", color="blue")
    ax1.loglog(vae_test_loss, label="VAE Test Loss", color="orange")
    ax1.set_xlabel("Epochs")
    ax1.set_ylabel("Loss")
    ax1.set_title("VAE Training and Test Loss")
    ax1.legend()
    ax1.grid(True)
    # Plot the Converter losses
    ax2.loglog(converterP2L_train_loss, label="Converter Train Loss", color="blue")
    ax2.loglog(converterP2L_test_loss, label="Converter Test Loss", color="orange")
    ax2.set_xlabel("Epochs")
    ax2.set_ylabel("Loss")
    ax2.set_title("Converter Training and Test Loss")
    ax2.legend()
    ax2.grid(True)

    # Plot the Converter losses
    ax3.loglog(converterL2F_train_loss, label="Converter Train Loss", color="blue")
    ax3.loglog(converterL2F_test_loss, label="Converter Test Loss", color="orange")
    ax3.set_xlabel("Epochs")
    ax3.set_ylabel("Loss")
    ax3.set_title("Converter Training and Test Loss")
    ax3.legend()
    ax3.grid(True)

    # Adjust layout to prevent overlap
    plt.tight_layout()
    # Save the plot
    plt.savefig(os.path.join(folder, "loss_plot.png"), dpi=300)
    # Show the plot
    plt.show()


def plot_laten_variable(folder, vae, converterP2L, train_loader, features_mean, features_std, device="cpu"):
    # Aggregate latent representations from both the VAE encoder and the Converter.
    all_mu_encoder = []  # Latent μ from VAE encoder.
    all_mu_converter = []  # Latent μ from Converter.
    all_logvar_encoder = []  # Latent logvar from VAE encoder.
    all_logvar_converter = []  # Latent logvar from Converter.
    all_params = []  # Polymer parameters (for coloring).

    vae.eval()
    converterP2L.eval()
    with torch.no_grad():
        for params, features, images in train_loader:
            images = images.to(device)
            # Get latent μ from the encoder.
            mu, logvar = vae.encoder(images)  # shape: [batch, latent_dim]
            all_mu_encoder.append(mu.detach().cpu().numpy())
            all_logvar_encoder.append(logvar.detach().cpu().numpy())

            # Get converter output, which now returns (z, μ, logvar).
            features_device = features.to(device)
            _, mu_conv, logvar_conv = converterP2L(features_device)  # shape: [batch, latent_dim]
            all_mu_converter.append(mu_conv.detach().cpu().numpy())
            all_logvar_converter.append(logvar_conv.detach().cpu().numpy())

            all_params.append(params.detach().cpu().numpy() * np.array(features_std) + np.array(features_mean))  # Denormalize features.

    # Concatenate arrays across batches.
    all_mu_encoder = np.concatenate(all_mu_encoder, axis=0)  # Shape: [N, 3]
    all_mu_converter = np.concatenate(all_mu_converter, axis=0)  # Shape: [N, 3]
    all_logvar_encoder = np.concatenate(all_logvar_encoder, axis=0)  # Shape: [N, 3]
    all_logvar_converter = np.concatenate(all_logvar_converter, axis=0)  # Shape: [N, 3]
    all_params = np.concatenate(all_params, axis=0)  # Shape: [N, 3]

    # Define feature names for coloring.
    feature_names = ["kappa", "f", "gL"]

    # Create a 2x3 grid of subplots.
    fig, axs = plt.subplots(4, 3, figsize=(18, 24), subplot_kw={"projection": "3d"})

    # Top row: latent μ from the VAE encoder.
    for i in range(3):
        sc = axs[0, i].scatter(all_mu_encoder[:, 0], all_mu_encoder[:, 1], all_mu_encoder[:, 2], c=all_params[:, i], cmap="rainbow", alpha=0.7)
        axs[0, i].set_xlabel("μ[0]")
        axs[0, i].set_ylabel("μ[1]")
        axs[0, i].set_zlabel("μ[2]")
        axs[0, i].set_title(f"Encoder Latent (colored by {feature_names[i]})")
        fig.colorbar(sc, ax=axs[0, i], shrink=0.5, label=feature_names[i])

    # second row: latent μ from the Converter.
    for i in range(3):
        sc = axs[1, i].scatter(all_mu_converter[:, 0], all_mu_converter[:, 1], all_mu_converter[:, 2], c=all_params[:, i], cmap="rainbow", alpha=0.7)
        axs[1, i].set_xlabel("μ_conv[0]")
        axs[1, i].set_ylabel("μ_conv[1]")
        axs[1, i].set_zlabel("μ_conv[2]")
        axs[1, i].set_title(f"Converter Latent (colored by {feature_names[i]})")
        fig.colorbar(sc, ax=axs[1, i], shrink=0.5, label=feature_names[i])

    # third row: latent logvar from the VAE encoder.
    for i in range(3):
        sc = axs[2, i].scatter(all_logvar_encoder[:, 0], all_logvar_encoder[:, 1], all_logvar_encoder[:, 2], c=all_params[:, i], cmap="rainbow", alpha=0.7)
        axs[2, i].set_xlabel("logvar[0]")
        axs[2, i].set_ylabel("logvar[1]")
        axs[2, i].set_zlabel("logvar[2]")
        axs[2, i].set_title(f"Encoder Logvar (colored by {feature_names[i]})")
        fig.colorbar(sc, ax=axs[2, i], shrink=0.5, label=feature_names[i])

    # fourth row: latent logvar from the Converter.
    for i in range(3):
        sc = axs[3, i].scatter(all_logvar_converter[:, 0], all_logvar_converter[:, 1], all_logvar_converter[:, 2], c=all_params[:, i], cmap="rainbow", alpha=0.7)
        axs[3, i].set_xlabel("logvar_conv[0]")
        axs[3, i].set_ylabel("logvar_conv[1]")
        axs[3, i].set_zlabel("logvar_conv[2]")
        axs[3, i].set_title(f"Converter Logvar (colored by {feature_names[i]})")
        fig.colorbar(sc, ax=axs[3, i], shrink=0.5, label=feature_names[i])

    # Save the latent representations and polymer params to a text file.
    combined_data = np.hstack((all_mu_encoder, all_mu_converter, all_logvar_encoder, all_logvar_converter, all_params))
    header = "encoder_mu0,encoder_mu1,encoder_mu2,converter_mu0,converter_mu1,converter_mu2,encoder_logvar0,encoder_logvar1,encoder_logvar2,converter_logvar0,converter_logvar1,converter_logvar2, param_kappa,param_f,param_gL"
    np.savetxt(os.path.join(folder, "latent_mu_logvar_data.txt"), combined_data, header=header, delimiter=",", comments="")
    print(f"Latent variables saved to {folder}/latent_mu_logvar_data.txt")

    plt.tight_layout()
    plt.savefig(os.path.join(folder, "latent_variable_plot.png"), dpi=300)
    # Show the plot
    plt.show()


def plot_LS_fitting(folder):
    # Load the fitted and ground truth parameters.
    data = np.loadtxt(os.path.join(folder, "ls_fitted_data.txt"), delimiter=",", skiprows=1)
    fitted_params = data[:, :3]  # Fitted parameters (first 3 columns).
    ground_truth_params = data[:, 3:]  # Ground truth parameters (last 3 columns).

    # Extract individual parameter arrays for easier plotting.
    fitted_kappa = fitted_params[:, 0]
    fitted_f = fitted_params[:, 1]
    fitted_gL = fitted_params[:, 2]

    gt_kappa = ground_truth_params[:, 0]
    gt_f = ground_truth_params[:, 1]
    gt_gL = ground_truth_params[:, 2]

    # Create a new figure
    fig, axs = plt.subplots(1, 3, figsize=(15, 5))

    # Plot kappa
    axs[0].scatter(gt_kappa, fitted_kappa, alpha=0.7)
    axs[0].plot([min(gt_kappa), max(gt_kappa)], [min(gt_kappa), max(gt_kappa)], 'r--')  # Diagonal line
    axs[0].set_xlabel("Ground Truth Kappa")
    axs[0].set_ylabel("Fitted Kappa")
    axs[0].set_title("Kappa Fitting")
    axs[0].set_xlim([min(gt_kappa), max(gt_kappa)])
    axs[0].set_ylim([min(gt_kappa), max(gt_kappa)])

    # Plot f
    axs[1].scatter(gt_f, fitted_f, alpha=0.7)
    axs[1].plot([min(gt_f), max(gt_f)], [min(gt_f), max(gt_f)], 'r--')  # Diagonal line
    axs[1].set_xlabel("Ground Truth f")
    axs[1].set_ylabel("Fitted f")
    axs[1].set_title("f Fitting")
    axs[1].set_xlim([min(gt_f), max(gt_f)])
    axs[1].set_ylim([min(gt_f), max(gt_f)])

    # Plot gL
    axs[2].scatter(gt_gL, fitted_gL, alpha=0.7)
    axs[2].plot([min(gt_gL), max(gt_gL)], [min(gt_gL), max(gt_gL)], 'r--')  # Diagonal line
    axs[2].set_xlabel("Ground Truth gL")
    axs[2].set_ylabel("Fitted gL")
    axs[2].set_title("gL Fitting")
    axs[2].set_xlim([min(gt_gL), max(gt_gL)])
    axs[2].set_ylim([min(gt_gL), max(gt_gL)])

    # Adjust layout to prevent overlap
    plt.tight_layout()
    # Save the plot
    plt.savefig(os.path.join(folder, "LS_fitting_plot.png"), dpi=300)
    # Show the plot
    plt.show()
    plt.close()


def plot_infer_feature(folder):
    # Move the model to the correct device and set it to evaluation mode.
    data = np.loadtxt(os.path.join(folder, "inferred_features.txt"), delimiter=",", skiprows=1)
    inferred_kappa, inferred_f, inferred_gL, gt_kappa, gt_f, gt_gL = data[:, 0], data[:, 1], data[:, 2], data[:, 3], data[:, 4], data[:, 5]

    # Create a new figure
    fig, axs = plt.subplots(1, 3, figsize=(15, 5))
    # Plot kappa
    axs[0].scatter(gt_kappa, inferred_kappa, alpha=0.7)
    axs[0].plot([min(gt_kappa), max(gt_kappa)], [min(gt_kappa), max(gt_kappa)], 'r--')  # Diagonal line
    axs[0].set_xlabel("Ground Truth Kappa")
    axs[0].set_ylabel("Inferred Kappa")
    axs[0].set_title("Kappa Inference")
    axs[0].set_xlim([min(gt_kappa), max(gt_kappa)])
    axs[0].set_ylim([min(gt_kappa), max(gt_kappa)])
    # Plot f
    axs[1].scatter(gt_f, inferred_f, alpha=0.7)
    axs[1].plot([min(gt_f), max(gt_f)], [min(gt_f), max(gt_f)], 'r--')  # Diagonal line
    axs[1].set_xlabel("Ground Truth f")
    axs[1].set_ylabel("Inferred f")
    axs[1].set_title("f Inference")
    axs[1].set_xlim([min(gt_f), max(gt_f)])
    axs[1].set_ylim([min(gt_f), max(gt_f)])
    # Plot gL
    axs[2].scatter(gt_gL, inferred_gL, alpha=0.7)
    axs[2].plot([min(gt_gL), max(gt_gL)], [min(gt_gL), max(gt_gL)], 'r--')  # Diagonal line
    axs[2].set_xlabel("Ground Truth gL")
    axs[2].set_ylabel("Inferred gL")
    axs[2].set_title("gL Inference")
    axs[2].set_xlim([min(gt_gL), max(gt_gL)])
    axs[2].set_ylim([min(gt_gL), max(gt_gL)])
    # Adjust layout to prevent overlap
    plt.tight_layout()
    # Save the plot
    plt.savefig(os.path.join(folder, "LS_inference_plot.png"), dpi=300)
    # Show the plot
    plt.show()
    plt.close()

