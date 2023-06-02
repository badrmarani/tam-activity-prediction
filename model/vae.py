import torch
from torch import nn, distributions

class Encoder(nn.Module):
    def __init__(self, data_dim, embedding_dim):
        super(Encoder, self).__init__()
        seq = []
        self.data_dim = data_dim
        self.seq = torch.nn.Sequential(
            torch.nn.Conv2d(3, 64, kernel_size=5, stride=2, padding=2),
            torch.nn.BatchNorm2d(64),
            torch.nn.ReLU(),
            torch.nn.Conv2d(64, 128, kernel_size=5, stride=2, padding=2),
            torch.nn.BatchNorm2d(128),
            torch.nn.ReLU(),
            torch.nn.Conv2d(128, 256, kernel_size=5, stride=2, padding=2),
            torch.nn.BatchNorm2d(256),
            torch.nn.ReLU(),
            torch.nn.Conv2d(256, 512, kernel_size=5, stride=2, padding=2),
            torch.nn.BatchNorm2d(512),
            torch.nn.ReLU(),
            torch.nn.Conv2d(512, 400, kernel_size=5, stride=2, padding=2),
            torch.nn.BatchNorm2d(400),
            torch.nn.ReLU(),
        )

        for i in range(5):
            mesh = self.output_features(data_dim, 5, 2, 2)
        self.fc1 = nn.Linear(400*mesh**2, embedding_dim)
        self.fc2 = nn.Linear(400*mesh**2, embedding_dim)
        self.seq = nn.Sequential(*seq)


    def output_features(self, size, kernel_size, stride, padding):
        return (((size - kernel_size) + (2*padding)) // stride) + 1

    def forward(self, x):
        x = x.view(x.size(0), 3, self.data_dim, self.data_dim)
        feat = self.seq(x)
        mu = self.fc1(feat)
        logvar = self.fc2(feat)
        return distributions.Normal(mu, logvar.mul(0.5).exp())
        

class Decoder(nn.Module):
    def __init__(self, data_dim):
        super(Decoder, self).__init__()

        padding = (0,0,0)
        if self.data_dim == 32:
            padding = (1,1,1)

        seq = torch.nn.Sequential(
            torch.nn.Conv3d(1, self.data_dim, kernel_size=4, stride=2, padding=(1, 1, 1)),
            torch.nn.BatchNorm3d(self.data_dim),
            torch.nn.LeakyReLU(self.args.leak_value),
            torch.nn.Conv3d(self.data_dim, self.data_dim*2, kernel_size=4, stride=2, padding=(1, 1, 1)),
            torch.nn.BatchNorm3d(self.data_dim*2),
            torch.nn.LeakyReLU(self.args.leak_value),
            torch.nn.Conv3d(self.data_dim*2, self.data_dim*4, kernel_size=4, stride=2, padding=(1, 1, 1)),
            torch.nn.BatchNorm3d(self.data_dim*4),
            torch.nn.LeakyReLU(self.args.leak_value),
            torch.nn.Conv3d(self.data_dim*4, self.data_dim*8, kernel_size=4, stride=2, padding=(1, 1, 1)),
            torch.nn.BatchNorm3d(self.data_dim*8),
            torch.nn.LeakyReLU(self.args.leak_value),
            torch.nn.Conv3d(self.data_dim*8, 1, kernel_size=4, stride=2, padding=padding),
            torch.nn.Sigmoid()
        )

        self.seq = nn.Sequential(*seq)
        self.scale = nn.Parameter(torch.ones(data_dim))

    def forward(self, x):
        x = x.view(-1, 1, self.data_dim, self.data_dim, self.data_dim)
        mu = self.seq(x)
        return distributions.Normal(mu, self.scale)

class VAE(nn.Module):
    def __init__(self, encoder, decoder):
        super(VAE, self).__init__()
        self.encoder = encoder
        self.decoder = decoder
    
    def forward(self, x):
        qzx = self.encoder(x)
        z = qzx.rsample()
        pxz = self.decoder(z)
        return qzx, pxz

def rsample(mu, logvar):
    std = logvar.mul(0.5).exp()
    return mu+std*torch.randn_like(std)
