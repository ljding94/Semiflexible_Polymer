import torch
import torch.nn.functional as F
import numpy as np
from tqdm import tqdm


def generate_scattering(converterP2L, decoder, params, device="cpu", enable_grad=False):
    """
    Generate scattering image given raw polymer parameters.

    The function first normalizes the raw polymer parameters using
    features_mean and features_std, passes them through the converter
    and decoder, and finally (optionally) denormalizes the output image
    using images_mean and images_std.

    Args:
        polymer_params (torch.Tensor): Raw polymer parameters [B,3].
        features_mean, features_std: Normalization stats for polymer parameters.
        images_mean, images_std: Normalization stats for images.
        enable_grad (bool): If True, gradients are computed (for fitting); otherwise, no_grad is used.
        denormalize_output (bool): If True, output is returned in raw scale; if False, output remains normalized.

    Returns:
        output (torch.Tensor): Generated scattering image.
    """
    converterP2L.to(device)
    decoder.to(device)
    converterP2L.eval()
    decoder.eval()

    def _forward():
        z, mu_conv, logvar_conv = converterP2L(params)
        norm_output = decoder(z)
        return norm_output

    if enable_grad:
        return _forward()
    else:
        with torch.no_grad():
            return _forward()


def fit_parameters(target_image, converterP2L, decoder, target_loss=1e-3, max_steps=1000, lr=1e-2, device="cpu"):
    """
    Fit raw polymer parameters (kappa, f, gL) using the Adam optimizer so that the generated
    scattering image best matches the target_image. The optimization stops when the loss is below
    a desired threshold (target_loss) or when max_steps is reached.

    **Note:** The target_image is assumed to be normalized.
    The generation function is called with denormalize_output=False so that the generated image
    is in normalized space, matching the target image.

    Boundaries:
      - kappa: [-2, 2]
      - f: [-2, 2]
      - gL: [-2, 2]
      # since are all normalized

    An improved initial guess is used (midpoint of allowed ranges).

    Args:
        target_image (torch.Tensor): Target scattering image [1, 1, H, W] (normalized).
        converter (nn.Module): Trained Converter network.
        decoder (nn.Module): Trained (and frozen) Decoder network.
        features_mean, features_std: Normalization stats for polymer parameters.
        images_mean, images_std: Normalization stats for images.
        target_loss (float): Desired loss threshold for early stopping.
        max_steps (int): Maximum number of iterations.
        lr (float): Learning rate.
        device (str): "cpu" or "cuda".

    Returns:
        fitted_raw (np.ndarray): Fitted polymer parameters in raw scale (shape [1, 3]).
    """
    # Move models to device, set to eval mode, and freeze their parameters.
    converterP2L.to(device)
    decoder.to(device)
    converterP2L.eval()
    decoder.eval()
    for param in converterP2L.parameters():
        param.requires_grad = False
    for param in decoder.parameters():
        param.requires_grad = False

    # Initialize raw polymer parameters at the midpoint of the allowed ranges:
    # kappa: [2,20] -> 11, f: [0,0.5] -> 0.25, gL: [0,2.0] -> 1.0.
    params = torch.tensor([[0.0, 0.0, 0.0]], device=device, requires_grad=True)

    optimizer = torch.optim.Adam([params], lr=lr)

    step = 0
    loss_val = float("inf")
    while step < max_steps and loss_val > target_loss:
        optimizer.zero_grad()
        # Generate scattering image in normalized space.
        gen = generate_scattering(converterP2L, decoder, params, device=device, enable_grad=True)
        loss = F.mse_loss(gen, target_image)
        loss.backward()
        optimizer.step()

        # Enforce boundaries for each parameter.
        with torch.no_grad():
            params[0, 0].clamp_(min=-2, max=2)  # kappa
            params[0, 1].clamp_(min=-2, max=2)  # f
            params[0, 2].clamp_(min=-2, max=2)  # gL

        loss_val = loss.item()
        step += 1

        if step % 100 == 0 or loss_val <= target_loss:
            print(f"Iteration {step}, loss: {loss_val:.6f}")

    if loss_val <= target_loss:
        print(f"Target loss achieved at step {step}.")
    else:
        print(f"Maximum steps reached ({max_steps}). Final loss: {loss_val:.6f}")

    fitted_raw = params.detach().cpu().numpy()
    return fitted_raw


# Fit parameters for all samples in the test_loader.
def fit_all_parameters(fit_loader, converterP2L, decoder, target_loss=1e-3, max_steps=100, lr=1e-2, device="cpu"):
    """
    For each sample in the fit_loader, optimize raw polymer parameters so that
    the generated scattering image best matches the target image.
    After fitting, denormalize the fitted parameters using features_mean and features_std.

    Returns:
        fitted_all (np.ndarray): Array of fitted raw polymer parameters, shape [N, 3].
        gt_all (np.ndarray): Array of ground-truth raw polymer parameters, shape [N, 3].
    """
    fitted_list = []
    gt_list = []

    # Move models to the correct device and set them to evaluation mode.
    converterP2L.to(device)
    decoder.to(device)
    converterP2L.eval()
    decoder.eval()
    # Freeze the parameters of the models.
    for param in converterP2L.parameters():
        param.requires_grad = False
    for param in decoder.parameters():
        param.requires_grad = False

    # Iterate over batches in the test loader.
    for batch_idx, (params, features, images) in enumerate(fit_loader):
        print(f"Fitting parameters for batch {batch_idx+1}/{len(fit_loader)}")
        # Iterate over each sample in the batch.
        for i in range(params.size(0)):
            print(f"Fitting parameters for sample {i+1}/{params.size(0)} in batch {batch_idx+1}/{len(fit_loader)}")
            # Extract one image at a time and move it to the device.
            target_image = images[i : i + 1].to(device)  # shape: [1, 1, H, W] in normalized space.
            # Fit the parameters for the current image.
            fitted_raw = fit_parameters(target_image, converterP2L, decoder, target_loss, max_steps, device=device)
            fitted_list.append(fitted_raw)

            # Denormalize ground truth: raw = normalized * std + mean.
            gt_raw = params[i : i + 1].cpu().numpy()
            gt_list.append(gt_raw)
            print("Fitted raw polymer parameters:", fitted_raw)
            print("Ground truth raw polymer parameters:", gt_raw)
            print("---------------------------------------------------")

    # Concatenate all the individual fitted and ground-truth parameters.
    fitted_all = np.concatenate(fitted_list, axis=0)
    gt_all = np.concatenate(gt_list, axis=0)
    return fitted_all, gt_all


def infer_feature(encoder, converterL2F, images, device="cpu"):
    """
    Given raw images, this function normalizes them using the provided image stats,
    passes them through the encoder to obtain a latent representation, and then uses
    the ConverterL2F module to infer the full (N,6) features. The output features are then
    optionally denormalized using the provided feature stats.
    """
    encoder.to(device)
    converterL2F.to(device)
    encoder.eval()
    converterL2F.eval()
    with torch.no_grad():
        mu, logvar = encoder(images)
        # For deterministic inference, we use mu as the latent code.
        z = mu
        norm_features = converterL2F(z)

    return norm_features


def infer_all_features(fit_loader, encoder, converterL2F,  device="cpu"):
    """
    Given a DataLoader containing images, this function iterates over the batches,
    infers the full features for each image using the encoder and converterL2F,
    and returns all the inferred features as a single array.
    """
    encoder.to(device)
    converterL2F.to(device)
    encoder.eval()
    converterL2F.eval()

    # freeze the models
    for param in encoder.parameters():
        param.requires_grad = False
    for param in converterL2F.parameters():
        param.requires_grad = False

    inferred_list = []
    gt_list = []

    for batch_idx, (params, features, images) in enumerate(fit_loader):
        print(f"Fitting parameters for batch {batch_idx+1}/{len(fit_loader)}")
        # Iterate over each sample in the batch.
        for i in range(features.size(0)):
            print(f"Fitting parameters for sample {i+1}/{features.size(0)} in batch {batch_idx+1}/{len(fit_loader)}")
            input_image = images[i : i + 1].to(device)  # shape: [1, 1, H, W] in normalized space.
            # infer the features for the current image.
            norm_features = infer_feature(encoder, converterL2F, input_image, device=device)
            inferred_list.append(norm_features.cpu().numpy())
            gt_list.append(features[i : i + 1].cpu().numpy())

            print("Inferred feature", norm_features.cpu().numpy())
            print("Ground truth feature", features[i : i + 1].cpu().numpy())
            print("---------------------------------------------------")
    # Concatenate all the individual inferred features.
    # Concatenate all the individual ground-truth features.
    gt_list = np.concatenate(gt_list, axis=0)
    inferred_list = np.concatenate(inferred_list, axis=0)
    return inferred_list, gt_list
