import matplotlib.pyplot as plt
from astropy.visualization import make_lupton_rgb
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
import argparse

import plot_smoothed_image

from PIL import Image
import os

def main(g_name, r_name, i_name, output_file, text=None, pix2sec=1., kpc_per_arc=1., hor_pos=0.03, vert_pos=0.8, L_bar=60.):
    g = fits.open(g_name)[0].data
    r = fits.open(r_name)[0].data
    i = fits.open(i_name)[0].data

    rgb_default = make_lupton_rgb(i, r, g, filename='tmp.png')
    
    img = Image.open('tmp.png')


    # Plot final image
    fig =  plt.figure(0)
    ax = plt.subplot()
    fsize = 13
    nx = img.size[0]
    ny = img.size[1]

    ax.imshow(img, origin="upper", aspect='equal')
    ax.set_frame_on(False)
    ax.set_xticks([]); ax.set_yticks([])
    plt.axis('off')    
    ax.text(hor_pos, vert_pos, text, fontsize=fsize, color='black',transform=ax.transAxes, horizontalalignment='left',verticalalignment='baseline',backgroundcolor='whitesmoke')
    
    plot_smoothed_image.add_scale_bar(ax, nx, ny, pix2sec, L_bar=L_bar, kpc_per_arc=kpc_per_arc)
    plt.draw() 

    plt.savefig(output_file, bbox_inches='tight', pad_inches=0.01, dpi = 300)
    plt.clf()
    plt.close()
    os.remove('tmp.png')
    


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--g', required=True)
    parser.add_argument('--r', required=True)
    parser.add_argument('--i', required=True)
    parser.add_argument('--o', default='gri.jpg', required=True)
    args = parser.parse_args()

    g_name = args.g
    r_name = args.r
    i_name = args.i
    output_file = args.o

    main(g_name, r_name, i_name, output_file)
