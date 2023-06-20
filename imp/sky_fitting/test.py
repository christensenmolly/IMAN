#!/bin/python3.7

"""

Use spline interpolation to on grid of x,y,z value where z is either xdiff or ydiff for use as imagemagick 2D displacement maps

"""

import numpy as np
from scipy.interpolate import SmoothBivariateSpline
from skimage import io


# python lists of x,y dst control points and zx=xsrc-xdiff, zy=ysrc-ydiff offsets to be interpolated over full image of size 129x129
x = [8.5, 20.5, 33.5, 48.5, 64.5, 80.5, 95.5, 109.5, 121.5, 5.5, 17.5, 31.5, 46.5, 64.5, 81.5, 97.5, 111.5, 123.5, 2.5, 14.5, 29.5, 45.5, 64.5, 83.5, 99.5, 113.5, 125.5, 1.5, 12.5, 26.5, 43.5, 64.5, 85.5, 103.5, 116.5, 127.5, 0.5, 11.5, 24.5, 41.5, 64.5, 87.5, 103.5, 117.5, 128.5, 1.5, 12.5, 25.5, 42.5, 64.5, 86.5, 103.5, 116.5, 127.5, 2.5, 14.5, 28.5, 45.5, 64.5, 83.5, 100.5, 114.5, 125.5, 5.5, 17.5, 30.5, 46.5, 64.5, 81.5, 97.5, 111.5, 123.5, 8.5, 19.5, 33.5, 48.5, 64.5, 80.5, 95.5, 109.5, 121.5]
y = [7.5, 5.5, 3.5, 1.5, 1.5, 1.5, 3.5, 5.5, 7.5, 20.5, 16.5, 14.5, 12.5, 11.5, 12.5, 15.5, 16.5, 19.5, 33.5, 31.5, 28.5, 26.5, 24.5, 26.5, 28.5, 31.5, 33.5, 48.5, 47.5, 45.5, 42.5, 40.5, 42.5, 45.5, 46.5, 48.5, 64.5, 64.5, 64.5, 64.5, 64.5, 64.5, 64.5, 64.5, 64.5, 80.5, 81.5, 83.5, 86.5, 87.5, 86.5, 83.5, 81.5, 80.5, 95.5, 97.5, 100.5, 103.5, 104.5, 102.5, 100.5, 97.5, 95.5, 109.5, 111.5, 114.5, 116.5, 117.5, 116.5, 114.5, 111.5, 109.5, 121.5, 123.5, 125.5, 127.5, 127.5, 127.5, 125.5, 123.5, 120.5]
zx = [119.5, 123.5, 126.5, 127.5, 127.5, 127.5, 128.5, 130.5, 134.5, 122.5, 126.5, 128.5, 129.5, 127.5, 126.5, 126.5, 128.5, 132.5, 125.5, 129.5, 130.5, 130.5, 127.5, 124.5, 124.5, 126.5, 130.5, 126.5, 131.5, 133.5, 132.5, 127.5, 122.5, 120.5, 123.5, 128.5, 127.5, 132.5, 135.5, 134.5, 127.5, 120.5, 120.5, 122.5, 127.5, 126.5, 131.5, 134.5, 133.5, 127.5, 121.5, 120.5, 123.5, 128.5, 125.5, 129.5, 131.5, 130.5, 127.5, 124.5, 123.5, 125.5, 130.5, 122.5, 126.5, 129.5, 129.5, 127.5, 126.5, 126.5, 128.5, 132.5, 119.5, 124.5, 126.5, 127.5, 127.5, 127.5, 128.5, 130.5, 134.5]
zy = [120.5, 122.5, 124.5, 126.5, 126.5, 126.5, 124.5, 122.5, 120.5, 123.5, 127.5, 129.5, 131.5, 132.5, 131.5, 128.5, 127.5, 124.5, 126.5, 128.5, 131.5, 133.5, 135.5, 133.5, 131.5, 128.5, 126.5, 127.5, 128.5, 130.5, 133.5, 135.5, 133.5, 130.5, 129.5, 127.5, 127.5, 127.5, 127.5, 127.5, 127.5, 127.5, 127.5, 127.5, 127.5, 127.5, 126.5, 124.5, 121.5, 120.5, 121.5, 124.5, 126.5, 127.5, 128.5, 126.5, 123.5, 120.5, 119.5, 121.5, 123.5, 126.5, 128.5, 130.5, 128.5, 125.5, 123.5, 122.5, 123.5, 125.5, 128.5, 130.5, 134.5, 132.5, 130.5, 128.5, 128.5, 128.5, 130.5, 132.5, 135.5]

# convert python lists to numpy arrays
ax = np.asarray(x)
ay = np.asarray(y)
azx = np.asarray(zx)
azy = np.asarray(zy)

print(azx)
exit()

# define bbox of interpolated data
# bbox=[minx, maxx, miny, maxy]
bbox=[0, 129, 0, 129]

# convert bbox to numpy array
abbox = np.asarray(bbox)

# do interpolations
xd = SmoothBivariateSpline(ax, ay, azx, w=None, bbox=abbox, kx=3, ky=3)
yd = SmoothBivariateSpline(ax, ay, azy, w=None, bbox=abbox, kx=3, ky=3)

# define integer grid onto which to interpolate
grid_x=np.linspace(0, 129, num=129)
grid_y=np.linspace(0, 129, num=129)


# evaluate at grid points
xdisplace = xd.__call__(grid_x, grid_y, grid=True)
ydisplace = yd.__call__(grid_x, grid_y, grid=True)

# save output using skimage
io.imsave("xdimgs.png", xdisplace.astype('uint8'))
io.imsave("ydimgs.png", ydisplace.astype('uint8'))

# view output using skimage
io.imshow(xdisplace.astype('uint8')) 
io.show()
io.imshow(ydisplace.astype('uint8')) 
io.show()