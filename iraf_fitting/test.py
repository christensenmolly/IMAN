import numpy as np
from astropy.modeling.models import Gaussian2D
from photutils.datasets import make_noise_image
g = Gaussian2D(100., 75, 75, 20, 12, theta=40.*np.pi/180.)
ny = nx = 150
y, x = np.mgrid[0:ny, 0:nx]
noise = make_noise_image((ny, nx), type='gaussian', mean=0.,
                         stddev=2., random_state=12345)
data = g(x, y) + noise
#print(data.shape)
#exit()
from photutils.isophote import EllipseGeometry
geometry = EllipseGeometry(x0=75, y0=75, sma=20, eps=0.5,
                           pa=20.*np.pi/180.)

from photutils.isophote import Ellipse
ellipse = Ellipse(data, geometry)

isolist = ellipse.fit_image()

print(isolist)

from photutils.isophote import build_ellipse_model
model_image = build_ellipse_model(data.shape, isolist)
residual = data - model_image