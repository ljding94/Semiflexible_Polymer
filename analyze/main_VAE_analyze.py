# analyze 2D Iq using VAE
from VAE_NN_analyze import *
from VAE_LS_analyze import *
from VAE_plot_analyze import *
from torch.utils.data import Subset


def main():

    folder = "../data/20250217"

    if 1:
        print("loading the data")
        L = 200
        rand_num = 5500
        rand_max = 4500
        parameters = generate_parameter_list(folder, L, rand_num, rand_max)
        print("Number of parameters:", len(parameters))

        batch_size = 32
        full_loader, dataset = create_dataloader(folder, parameters, batch_size=batch_size, shuffle=True)

        # save the normalization parameters for the data
        # Save normalization parameters separately as image_mean and image_std are likely 2D arrays
        normalization_data = {
            "param_mean": dataset.param_mean,
            "param_std": dataset.param_std,
            "features_mean": dataset.features_mean,
            "features_std": dataset.features_std,
            "images_mean": dataset.images_mean,
            "images_std": dataset.images_std,
            "qB": dataset.qB,
        }

        print("qB:", dataset.qB)

        np.savez(os.path.join(folder, "normalization_params.npz"), **normalization_data)
        print(f"Normalization parameters saved to {folder}/normalization_params.npz")

        # Split the dataset into training and test sets (e.g., 80% train, 20% test).
        dataset_size = len(dataset)
        test_size = int(0.1 * dataset_size)
        train_size = dataset_size - test_size
        indices = list(range(dataset_size))
        train_indices = indices[:train_size]
        test_indices = indices[train_size:]
        train_dataset = Subset(dataset, train_indices)
        test_dataset = Subset(dataset, test_indices)
        fit_indices = indices[train_size + test_size // 2 :]
        fit_dataset = Subset(dataset, fit_indices)

        train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
        test_loader = DataLoader(test_dataset, batch_size=batch_size, shuffle=False)
        fit_loader = DataLoader(fit_dataset, batch_size=batch_size, shuffle=False)  # used for least suqare fitting

        sample_params, sample_features, sample_images = next(iter(train_loader))
        nq = sample_images.shape[-1]  # Assuming images are [B, 1, H, W] and square.
        print("Detected image size:", nq, "x", nq)
        latent_dim = 3  # Use 3 to directly correspond to the 3 polymer parameters.
        device = "cuda" if torch.cuda.is_available() else "cpu"

        print("successfully loaded the data")
        print("train_loader length", len(train_loader))
        print("test_loader length", len(test_loader))

    if 0:
        print("traiing the NN")
        # train the model and save the model
        # train_and_save_model(folder, train_loader, test_loader, latent_dim, nq, vae_epoch=100, converter_epoch=200, device=device)
        train_and_save_model(folder, train_loader, test_loader, latent_dim, nq, vae_epoch=100, converter_epoch=200, end_to_end_epoch=100, device=device)

    if 0:
        # illustrate the training
        # plot the loss
        plot_loss(folder)

        # plot the latent variable
        features_mean = dataset.features_mean  # shape: (3,)
        features_std = dataset.features_std  # shape: (3,)
        vae, converterP2L, converterL2F = load_models(latent_dim, nq, folder)
        plot_laten_variable(folder, vae, converterP2L, train_loader, features_mean, features_std, device=device)

    if 0:
        print("testing the NN: comparing the laten variable")
        # Load the trained models.
        vae, converterP2L, converterL2F = load_models(latent_dim, nq, folder)

        s_params = torch.tensor([[0.8, 1.2, 0.5]], dtype=torch.float32)  # raw polymer parameters
        generated_image = generate_scattering(converterP2L, vae.decoder, s_params, device=device)
        print("Generated image shape:", generated_image.shape)  # Expected: [1, 1, nq, nq]

    if 0:
        # test fitting parameters
        # Load the trained models.
        vae, converterP2L, converterL2F = load_models(latent_dim, nq, folder)
        decoder = vae.decoder

        # Get one sample from the test set.
        fit_iter = iter(fit_loader)
        params_batch, features_batch, images_batch = next(fit_iter)  # params, features shape: [B, 3], images shape: [B, 1, H, W]
        # Pick the first sample.

        target_image = images_batch[1:2].to(device)  # normalized target image.
        gt_normalized = params_batch[1:2]  # normalized ground truth parameters.

        # Fit parameters for the target image.
        fitted_raw = fit_parameters(
            target_image,
            converterP2L,
            decoder,
            target_loss=2e-4,
            max_steps=10000,
            lr=1e-2,
            device=device,
        )

        # Denormalize the ground truth polymer parameters.
        gt_raw = gt_normalized.cpu().numpy()

        print("Fitted raw polymer parameters:", fitted_raw)
        print("Ground truth raw polymer parameters:", gt_raw)

    if 0:
        # fit all params
        # Load the trained models.
        vae, converterP2L, converterL2F = load_models(latent_dim, nq, folder)
        decoder = vae.decoder

        fitted_params, ground_truth_params = fit_all_parameters(
            fit_loader,
            converterP2L,
            decoder,  # use the pretrained (and frozen) decoder
            target_loss=5e-4,
            max_steps=5000,
            lr=1e-2,
            device=device,
        )

        # Denormalize the fitted parameters.
        fitted_params = fitted_params * np.array(dataset.features_std) + np.array(dataset.features_mean)
        # Denormalize the ground truth polymer parameters.
        ground_truth_params = ground_truth_params * np.array(dataset.features_std) + np.array(dataset.features_mean)

        combined_data = np.hstack((fitted_params, ground_truth_params))
        header = "fitted_kappa, fitted_f, fitted_gL, gt_kappa, gt_f, gt_gL"
        # Save the fitted and ground truth parameters to a text file.
        np.savetxt(os.path.join(folder, "ls_fitted_data.txt"), combined_data, header=header, delimiter=",", comments="")
        print(f"Latent variables saved to {folder}/ls_fitted_data.txt")

    if 0:
        # Plot the fitted parameters vs. ground truth.
        plot_LS_fitting(folder)

    if 0:
        # test infer feature
        # Load the trained models.
        vae, converterP2L, converterL2F = load_models(latent_dim, nq, folder)
        encoder = vae.encoder
        # Get one sample from the test set.
        fit_iter = iter(fit_loader)
        params_batch, features_batch, images_batch = next(fit_iter)
        # Pick the first sample.
        input_image = images_batch[1:2].to(device)  # normalized input image.
        gt_normalized = features_batch[1:2]  # normalized ground truth parameters.
        # Infer the polymer features from the target image.
        inferred_features = infer_feature(encoder, converterL2F, input_image, device=device)
        # Denormalize the ground truth polymer parameters.
        gt_raw = gt_normalized.cpu().numpy()
        print("Inferred polymer features:", inferred_features)
        print("Ground truth polymer features:", gt_raw)
        print("---------------------------------------------------")

    if 0:
        # inferred all features
        # Load the trained models.
        vae, converterP2L, converterL2F = load_models(latent_dim, nq, folder)
        encoder = vae.encoder
        inferred_features, ground_truth_features = infer_all_features(
            fit_loader,
            encoder,
            converterL2F,
            device=device,
        )

        # Denormalize the ground truth polymer parameters.
        ground_truth_features = ground_truth_features * np.array(dataset.features_std) + np.array(dataset.features_mean)
        # Denormalize the inferred features.
        inferred_features = inferred_features * np.array(dataset.features_std) + np.array(dataset.features_mean)

        # Save the inferred features and ground truth features to a text file.
        combined_data = np.hstack((inferred_features, ground_truth_features))
        header = "inferred_kappa, inferred_f, inferred_gL, gt_kappa, gt_f, gt_gL"
        np.savetxt(os.path.join(folder, "inferred_features.txt"), combined_data, header=header, delimiter=",", comments="")
        print(f"Inferred features saved to {folder}/inferred_features.txt")

    if 0:
        plot_infer_feature(folder)


if __name__ == "__main__":
    main()
