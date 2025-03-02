# all NN function for VAE

import os
import numpy as np
import torch
from torch.utils.data import Dataset, DataLoader
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim


# --------------------------
# Data Loading Functions
# --------------------------
def generate_parameter_list(folder, L, rand_num, rand_max):
    parameters = []
    for i in range(rand_num):
        filename = f"{folder}/obs_L{L}_random_run{i}.csv"
        if os.path.exists(filename):
            parameters.append([L, i])
        if len(parameters) >= rand_max:
            break
    return parameters


def load_data(folder, parameters):
    all_kappa, all_f, all_gL = [], [], []
    all_R2, all_Rg2, all_Rxz = [], [], []
    all_Sq2D = []
    qB = None
    for L, run_num in parameters:
        filename = f"{folder}/obs_L{L}_random_run{run_num}.csv"
        if not os.path.exists(filename):
            continue
        data = np.genfromtxt(filename, delimiter=",", skip_header=1)
        # Extract polymer parameters
        kappa, f, g = data[0, 2:5]
        R2, Rg2, Rxz = data[0, 14], data[0, 15], data[0, 20]

        all_kappa.append(kappa)
        all_f.append(f)
        all_gL.append(g * L)

        all_R2.append(R2 / L**2)
        all_Rg2.append(Rg2 / L**2)
        all_Rxz.append(Rxz / Rg2)

        qB = data[3, 21:]
        Sq2D = data[7:, 21:]
        all_Sq2D.append(Sq2D)
    # Create the (N,3) params array and the (N,6) features (for later use)
    params = np.array([all_kappa, all_f, all_gL]).T  # shape: (N, 3)
    features = np.array([all_kappa, all_f, all_gL, all_R2, all_Rg2, all_Rxz]).T  # shape: (N, 6)
    features = params.copy()  # (N,3) params
    qB = np.array(qB)
    # Apply log-transform on images.
    images = -np.log(all_Sq2D)
    return params, features, images, qB


class ScatteringDataset(Dataset):
    def __init__(self, folder, parameters, transform=None, normalize_data=True, normalize_images=True):
        # Unpack raw data from load_data.
        raw_params, raw_features, raw_images, self.qB = load_data(folder, parameters)

        if normalize_data:
            # Normalize the (N,3) polymer parameters.
            self.param_mean = raw_params.mean(axis=0)
            self.param_std = raw_params.std(axis=0)
            param_norm = (raw_params - self.param_mean) / self.param_std

            # Normalize the (N,6) features.
            self.features_mean = raw_features.mean(axis=0)
            self.features_std = raw_features.std(axis=0)
            features_norm = (raw_features - self.features_mean) / self.features_std
        else:
            param_norm = raw_params
            features_norm = raw_features
            self.param_mean, self.param_std = None, None
            self.features_mean, self.features_std = None, None

        # Normalize images globally.
        if normalize_images:
            stacked = np.stack(raw_images, axis=0)  # shape: (N, H, W)
            self.images_mean = stacked.mean()
            self.images_std = stacked.std()
            images_norm = [(img - self.images_mean) / self.images_std for img in raw_images]
        else:
            images_norm = raw_images
            self.images_mean, self.images_std = None, None

        # Save the normalized arrays.
        self.params = param_norm  # (N, 3)
        self.features = features_norm  # (N, 6)
        self.images = images_norm  # list of (H, W) arrays
        self.transform = transform

    def __len__(self):
        return len(self.params)

    def __getitem__(self, idx):
        params = self.params[idx]  # (3,)
        feature = self.features[idx]  # (6,)
        image = self.images[idx]  # (H, W)
        if image.ndim == 2:
            image = image[np.newaxis, :, :]  # (1, H, W)
        params = torch.tensor(params, dtype=torch.float32)
        feature = torch.tensor(feature, dtype=torch.float32)
        image = torch.tensor(image, dtype=torch.float32)
        if self.transform:
            image = self.transform(image)
        return params, feature, image


def create_dataloader(folder, parameters, batch_size=32, shuffle=True, transform=None):
    dataset = ScatteringDataset(folder, parameters, transform=transform)
    return DataLoader(dataset, batch_size=batch_size, shuffle=shuffle), dataset


# --------------------------
# Model Definitions
# --------------------------
class Encoder(nn.Module):
    def __init__(self, latent_dim, nq):
        super(Encoder, self).__init__()
        self.conv = nn.Sequential(
            nn.Conv2d(1, 32, kernel_size=3, stride=2, padding=1), nn.ReLU(), nn.Conv2d(32, 64, kernel_size=3, stride=2, padding=1), nn.ReLU()  # [B,1,nq,nq] -> [B,32,nq/2,nq/2] -> [B,64,nq/4,nq/4]
        )
        dummy = torch.zeros(1, 1, nq, nq)
        conv_out = self.conv(dummy)
        self.flatten_dim = conv_out.numel()
        self.fc_mu = nn.Linear(self.flatten_dim, latent_dim)
        self.fc_logvar = nn.Linear(self.flatten_dim, latent_dim)

    def forward(self, x):
        batch_size = x.size(0)
        x = self.conv(x)
        x = x.view(batch_size, -1)
        return self.fc_mu(x), self.fc_logvar(x)


class Decoder(nn.Module):
    def __init__(self, latent_dim, nq):
        super(Decoder, self).__init__()
        self.nq = nq
        dummy = torch.zeros(1, 1, nq, nq)
        conv_out = nn.Sequential(nn.Conv2d(1, 32, kernel_size=3, stride=2, padding=1), nn.ReLU(), nn.Conv2d(32, 64, kernel_size=3, stride=2, padding=1), nn.ReLU())(dummy)
        self.channels = conv_out.size(1)
        self.height = conv_out.size(2)
        self.width = conv_out.size(3)
        self.flatten_dim = self.channels * self.height * self.width
        self.fc = nn.Linear(latent_dim, self.flatten_dim)
        self.deconv = nn.Sequential(
            nn.ReLU(),
            nn.ConvTranspose2d(self.channels, 32, kernel_size=4, stride=2, padding=1),
            nn.ReLU(),
            nn.ConvTranspose2d(32, 1, kernel_size=4, stride=2, padding=1),
        )

    def forward(self, z):
        batch_size = z.size(0)
        x = self.fc(z)
        x = x.view(batch_size, self.channels, self.height, self.width)
        x = self.deconv(x)
        if x.shape[-1] != self.nq:
            x = F.interpolate(x, size=(self.nq, self.nq), mode="bilinear", align_corners=False)
        return x


class VAE(nn.Module):
    def __init__(self, latent_dim, nq):
        super(VAE, self).__init__()
        self.encoder = Encoder(latent_dim, nq)
        self.decoder = Decoder(latent_dim, nq)

    def reparameterize(self, mu, logvar):
        std = torch.exp(0.5 * logvar)
        eps = torch.randn_like(std)
        return mu + eps * std

    def forward(self, x):
        mu, logvar = self.encoder(x)
        z = self.reparameterize(mu, logvar)
        return self.decoder(z), mu, logvar


class ConverterP2L(nn.Module):
    def __init__(self, latent_dim):
        super(ConverterP2L, self).__init__()  # params to latent
        # Map params (size 3) to a hidden representation.
        self.shared = nn.Sequential(
            nn.Linear(3, 9),
            nn.BatchNorm1d(9),
            nn.Tanh(),
        )
        # Branch for μ and logvar.
        self.fc_mu = nn.Linear(9, latent_dim)
        self.fc_logvar = nn.Linear(9, latent_dim)

    def reparameterize(self, mu, logvar):
        std = torch.exp(0.5 * logvar)
        eps = torch.randn_like(std)
        return mu + eps * std

    def forward(self, params):
        h = self.shared(params)
        mu = self.fc_mu(h)
        logvar = self.fc_logvar(h)
        z = self.reparameterize(mu, logvar)
        return z, mu, logvar


class ConverterL2F(nn.Module):
    # Convert latent variables to the full set of (N,6) features.
    def __init__(self, latent_dim):
        super(ConverterL2F, self).__init__()  # latent to features
        self.fc = nn.Sequential(
            nn.Linear(latent_dim, 9),
            nn.BatchNorm1d(9),
            nn.Tanh(),
            nn.Linear(9, 3),  # output: (N,6) features
        )

    def forward(self, z):
        return self.fc(z)


# -------------------------
# save and load data
# -------------------------
def save_trained_models(vae, converterP2L, converterL2F, folder):
    """
    Save the trained VAE and Converter model state dictionaries.
    """
    torch.save(vae.state_dict(), f"{folder}/vae_model_state.pth")
    torch.save(converterP2L.state_dict(), f"{folder}/converterP2L_model_state.pth")
    torch.save(converterL2F.state_dict(), f"{folder}/converterL2F_model_state.pth")
    print(f"VAE state_dict saved to {folder}/vae_model_state.pth")
    print(f"ConverterP2L state_dict saved to {folder}/converterP2L_model_state.pth")
    print(f"ConverterL2F state_dict saved to {folder}/converterL2F_model_state.pth")


def load_models(latent_dim, nq, folder, device="cpu"):
    """
    Instantiate the VAE and Converter models, load their state dictionaries,
    and set them to evaluation mode.

    Args:
        latent_dim (int): The dimension of the latent space.
        nq (int): The input image size (assumed square, nq x nq).
        device (str): The device to load the models on ("cpu" or "cuda").
        vae_filename (str): Filename for the VAE state dictionary.
        converter_filename (str): Filename for the Converter state dictionary.

    Returns:
        vae, converter: The loaded model instances.
    """
    # Instantiate the models with the same hyperparameters used during training.
    vae = VAE(latent_dim, nq)
    converterP2L = ConverterP2L(latent_dim)
    converterL2F = ConverterL2F(latent_dim)

    # Load state dictionaries.
    vae.load_state_dict(torch.load(f"{folder}/vae_model_state.pth", map_location=device))
    converterP2L.load_state_dict(torch.load(f"{folder}/converterP2L_model_state.pth", map_location=device))
    converterL2F.load_state_dict(torch.load(f"{folder}/converterL2F_model_state.pth", map_location=device))

    # Move models to the device and set to evaluation mode.
    vae.to(device)
    converterP2L.to(device)
    converterL2F.to(device)
    vae.eval()
    converterP2L.eval()
    converterL2F.eval()

    print(f"Loaded VAE state_dict from {folder}/vae_model_state.pth")
    print(f"Loaded ConverterP2L state_dict from {folder}/converterP2L_model_state.pth")
    print(f"Loaded ConverterL2F state_dict from {folder}/converterL2F_model_state.pth")
    print("Models are set to evaluation mode.")
    return vae, converterP2L, converterL2F


def save_losses(folder, vae_train_loss, vae_test_loss, converterP2L_train_loss, converterP2L_test_loss, converterL2F_train_loss, converterL2F_test_loss):
    data = np.column_stack((vae_train_loss, vae_test_loss))
    column_names = ["VAE Train Loss", "VAE Test Loss"]
    np.savetxt(f"{folder}/vae_losses.txt", data, delimiter=",", header=",".join(column_names), comments="")

    data = np.column_stack((converterP2L_train_loss, converterP2L_test_loss))
    column_names = ["ConverterP2L Train Loss", "ConverterP2L Test Loss"]
    np.savetxt(f"{folder}/converterP2L_losses.txt", data, delimiter=",", header=",".join(column_names), comments="")

    data = np.column_stack((converterL2F_train_loss, converterL2F_test_loss))
    column_names = ["ConverterL2F Train Loss", "ConverterL2F Test Loss"]
    np.savetxt(f"{folder}/converterL2F_losses.txt", data, delimiter=",", header=",".join(column_names), comments="")

    print(f"VAE losses saved to {folder}/vae_losses.txt")
    print(f"ConverterP2L losses saved to {folder}/converterP2L_losses.txt")
    print(f"ConverterL2F losses saved to {folder}/converterL2F_losses.txt")


# --------------------------
# Training Functions
# --------------------------
def train_vae(vae, train_loader, test_loader, num_epochs=50, lr=1e-3, device="cpu"):
    optimizer = optim.Adam(vae.parameters(), lr=lr)
    vae.to(device)
    vae.train()
    all_train_loss, all_test_loss = [], []
    for epoch in range(num_epochs):
        train_loss = 0
        test_loss = 0
        for _, _, images in train_loader:  # VAE training uses only images.
            images = images.to(device)
            optimizer.zero_grad()
            recon, mu, logvar = vae(images)
            recon_loss = F.mse_loss(recon, images, reduction="mean")
            loss = recon_loss  # (Optional KL loss term can be added)
            loss.backward()
            optimizer.step()
            train_loss += loss.item()
        for _, _, images in test_loader:
            images = images.to(device)
            recon, mu, logvar = vae(images)
            recon_loss = F.mse_loss(recon, images, reduction="mean")
            loss = recon_loss
            test_loss += loss.item()
        train_loss /= len(train_loader.dataset)
        test_loss /= len(test_loader.dataset)
        all_train_loss.append(train_loss)
        all_test_loss.append(test_loss)

        print(f"VAE Epoch {epoch+1}, train Loss: {train_loss:.6f}, test Loss: {test_loss:.6f}")
    return all_train_loss, all_test_loss


def train_converterP2L(converterP2L, decoder, train_loader, test_loader, num_epochs=30, lr=1e-3, device="cpu"):
    """
    Train the ConverterP2L module to map polymer parameters to latent codes.
    The decoder is frozen and used to compute the image reconstruction loss.
    """
    # Freeze the decoder.
    decoder.to(device)
    for p in decoder.parameters():
        p.requires_grad = False
    decoder.eval()

    optimizer = optim.Adam(converterP2L.parameters(), lr=lr)
    converterP2L.to(device)
    converterP2L.train()

    all_train_loss, all_test_loss = [], []

    for epoch in range(num_epochs):
        train_loss = 0
        test_loss = 0
        for params, features, images in test_loader:
            params = params.to(device)
            images = images.to(device)
            optimizer.zero_grad()
            # Map polymer parameters (N,3) to latent space.
            z, mu, logvar = converterP2L(params)
            # Generate images using the frozen decoder.
            image_output = decoder(z)
            loss = F.mse_loss(image_output, images)
            loss.backward()
            optimizer.step()
            train_loss += loss.item()
        for params, features, images in test_loader:
            params = params.to(device)
            images = images.to(device)
            # Map polymer parameters (N,3) to latent space.
            z, mu, logvar = converterP2L(params)
            # Generate images using the frozen decoder.
            image_output = decoder(z)
            loss = F.mse_loss(image_output, images)
            test_loss += loss.item()
        train_loss /= len(train_loader.dataset)
        test_loss /= len(test_loader.dataset)
        all_train_loss.append(train_loss)
        all_test_loss.append(test_loss)
        print(f"[ConverterP2L] Epoch {epoch+1}, Train Loss: {train_loss:.6f}, Test Loss: {test_loss:.6f}")

    # unfreeze the decoder after training
    for p in decoder.parameters():
        p.requires_grad = True
    decoder.train()

    return all_train_loss, all_test_loss


def train_converterL2F(encoder, converterL2F, train_loader, test_loader, num_epochs=30, lr=1e-3, device="cpu"):
    """
    Train the ConverterL2F module to map latent codes to full features.
    The ConverterP2L module is frozen so that it provides a stable latent representation.
    """

    # freeze the vae
    encoder.to(device)
    for p in encoder.parameters():
        p.requires_grad = False
    encoder.eval()

    optimizer = optim.Adam(converterL2F.parameters(), lr=lr)
    converterL2F.to(device)
    converterL2F.train()

    all_train_loss, all_test_loss = [], []

    for epoch in range(num_epochs):
        train_loss = 0
        test_loss = 0
        for params, features, images in test_loader:
            params = params.to(device)
            features = features.to(device)
            optimizer.zero_grad()

            # Obtain latent representation from the frozen encoder
            images = images.to(device)
            mu, logvar = encoder(images)
            z = mu + torch.randn_like(mu) * torch.exp(0.5 * logvar)  # reparameterization
            # Map latent code to features.
            features_output = converterL2F(z)
            loss = F.mse_loss(features_output, features)
            loss.backward()
            optimizer.step()
            train_loss += loss.item()

        for params, features, images in test_loader:
            params = params.to(device)
            features = features.to(device)
            # Obtain latent representation from the frozen encoder
            images = images.to(device)
            mu, logvar = encoder(images)
            z = mu + torch.randn_like(mu) * torch.exp(0.5 * logvar)
            # Map latent code to features.
            features_output = converterL2F(z)
            loss = F.mse_loss(features_output, features)
            test_loss += loss.item()
        train_loss /= len(train_loader.dataset)
        test_loss /= len(test_loader.dataset)
        all_train_loss.append(train_loss)
        all_test_loss.append(test_loss)
        print(f"[ConverterL2F] Epoch {epoch+1}, Train Loss: {train_loss:.6f}, Test Loss: {test_loss:.6f}")

    # unfreeze converterP2L after training if further joint fine-tuning is desired.
    for p in encoder.parameters():
        p.requires_grad = True
    encoder.train()

    return all_train_loss, all_test_loss


def train_end_to_end(encoder, decoder, converterP2L, converterL2F, dataloader, num_epochs=50, lr=1e-4, device="cpu"):
    """
    End-to-end fine-tuning of the entire network:
    encoder, decoder, converterP2L, and converterL2F are trained together.
    The combined loss includes:
      - VAE reconstruction loss (using the encoder & decoder)
      - Converter image loss (using converterP2L & decoder)
      - Converter feature loss (using converterP2L & converterL2F)
    """
    # Combine parameters from all modules.
    all_params = list(encoder.parameters()) + list(decoder.parameters()) + list(converterP2L.parameters()) + list(converterL2F.parameters())
    optimizer = optim.Adam(all_params, lr=lr)

    # Move all modules to device and set to train mode.
    encoder.to(device)
    decoder.to(device)
    converterP2L.to(device)
    converterL2F.to(device)
    encoder.train()
    decoder.train()
    converterP2L.train()
    converterL2F.train()

    for epoch in range(num_epochs):
        total_loss = 0
        for params, features, images in dataloader:
            params = params.to(device)
            features = features.to(device)
            images = images.to(device)
            optimizer.zero_grad()

            # VAE branch: Encode images and reconstruct.
            mu_enc, logvar_enc = encoder(images)
            z_enc = mu_enc + torch.randn_like(mu_enc) * torch.exp(0.5 * logvar_enc)  # reparameterization
            recon_images = decoder(z_enc)
            vae_loss = F.mse_loss(recon_images, images)

            # Converter branch: Map polymer params to latent space.
            z_conv, mu_conv, logvar_conv = converterP2L(params)
            conv_image = decoder(z_conv)
            conv_features = converterL2F(z_conv)
            conv_image_loss = F.mse_loss(conv_image, images)
            conv_feature_loss = F.mse_loss(conv_features, features)

            # Combined loss.
            loss = vae_loss + conv_image_loss + conv_feature_loss
            loss.backward()
            optimizer.step()
            total_loss += loss.item()
        print(f"End-to-End Epoch {epoch+1}, Loss: {total_loss/len(dataloader):.6f}")


def train_and_save_model(folder, train_loader, test_loader, latent_dim, nq, vae_epoch=50, converter_epoch=200, end_to_end_epoch=100, device="cpu"):
    # Automatically determine image size (nq) from the dataset.
    sample_params, sample_features, sample_images = next(iter(train_loader))
    nq = sample_images.shape[-1]  # Assuming images are [B, 1, H, W] and square.
    print("Detected image size:", nq, "x", nq)

    vae = VAE(latent_dim, nq)
    converterP2L = ConverterP2L(latent_dim)
    converterL2F = ConverterL2F(latent_dim)

    print("Training VAE...")
    vae_train_loss, vae_test_loss = train_vae(vae, train_loader, test_loader, num_epochs=vae_epoch, lr=1e-3, device=device)

    print("Training ConverterP2L...")
    converterP2L_train_loss, converterP2L_test_loss = train_converterP2L(converterP2L, vae.decoder, train_loader, test_loader, num_epochs=converter_epoch, lr=1e-3, device=device)

    print("Training ConverterL2F...")
    converterL2F_train_loss, converterL2F_test_loss = train_converterL2F(vae.encoder, converterL2F, train_loader, test_loader, num_epochs=converter_epoch, lr=1e-3, device=device)

    # train end to end
    print("Training end-to-end...")
    train_end_to_end(vae.encoder, vae.decoder, converterP2L, converterL2F, train_loader, num_epochs=end_to_end_epoch, lr=1e-4, device=device)

    # Save the trained models
    save_trained_models(vae, converterP2L, converterL2F, folder)
    print("Models saved successfully.")

    # save the training and test loss
    save_losses(folder, vae_train_loss, vae_test_loss, converterP2L_train_loss, converterP2L_test_loss, converterL2F_train_loss, converterL2F_test_loss)
    print("Losses saved successfully.")
