import numpy as np
import matplotlib.pyplot as plt
import scipy.ndimage as ndi
from astropy.visualization import simple_norm
from astropy.modeling import models
import photutils
import time
import statmorph
from statmorph.utils.image_diagnostics import make_figure



ny, nx = 2400, 2400
y, x = np.mgrid[0:ny, 0:nx]
sersic_model = models.Sersic2D(amplitude=1, r_eff=200, n=1.5, x_0=0.5*nx, y_0=0.4*ny,
                               ellip=0.5, theta=0.5)
image = sersic_model(x, y)

size = 20  # on each side from the center
sigma_psf = 2.0
y, x = np.mgrid[-size:size+1, -size:size+1]
psf = np.exp(-(x**2 + y**2)/(2.0*sigma_psf**2))
psf /= np.sum(psf)

image = ndi.convolve(image, psf)

np.random.seed(1)
snp = 100.0
image += (1.0 / snp) * np.random.standard_normal(size=(ny, nx))

gain = 1000.0


threshold = photutils.detect_threshold(image, 1.5)
npixels = 5  # minimum number of connected pixels
segm = photutils.detect_sources(image, threshold, npixels)

label = np.argmax(segm.areas) + 1
segmap = segm.data == label

segmap_float = ndi.uniform_filter(np.float64(segmap), size=10)
segmap = segmap_float > 0.5

source_morphs = statmorph.source_morphology(image, segmap, gain=gain, psf=psf)

morph = source_morphs[0]
print(morph.sersic_rhalf)

fig = make_figure(morph)