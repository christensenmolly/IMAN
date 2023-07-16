# I put any useful functions in this module
# They make the code readable and more understandable. (I hope...)
#===================================================================

# Loading the required packages

# Astropy
from astroquery.vizier import Vizier 
import astropy.coordinates as coord
from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as u

from matplotlib import pyplot as plt
from numpy import sin, cos, arctan, degrees, radians, sqrt
import numpy as np
from pathlib import Path

from photutils.aperture import aperture_photometry
from photutils.aperture import CircularAperture


#-------------------------------------------------------------------

# My functions with some comments

def get_file(fname):
    # Gets data & headef from FITS
    hdu    = fits.open(fname)
    data   = hdu[0].data
    header = hdu[0].header
    return data, header

#-------------------------------------------------------------------
def angular_distance(ra, dec, RA, DEC):
    # INPUT AND OUTPUT DATA IN RADIANS!
    vec1 = [cos(ra)*cos(dec), sin(ra)*cos(dec), sin(dec)]
    vec2 = [cos(RA)*cos(DEC), sin(RA)*cos(DEC), sin(DEC)]
    return np.arccos(np.dot(vec1, vec2))

#-------------------------------------------------------------------
def get_image_center(filelist):
    # returns the eq. coordinates (RA & DEC) of the center of the image
    # A list of images is submitted to the input. 
    # This can be useful in the case of several close frames
    # OUTPUT DATA IN DEGREES!
    
    max_radius = 0
    ras    = []
    decs   = []

    for fname in filelist:
        #print(fname)
        imdata, header = get_file(fname)
        try:
            w = header['IMAGEW']
            h = header['IMAGEH']
        except:
            h, w = imdata.shape
        W = WCS(header)

        RAc, DECc = W.wcs_pix2world(np.array([[w/2, h / 2]], np.float_), 1)[0]
        RA1, DEC1 = W.wcs_pix2world(np.array([[0, 0]], np.float_), 1)[0]
        RA2, DEC2 = W.wcs_pix2world(np.array([[w, h]], np.float_), 1)[0]
        radius = degrees(angular_distance(radians(RA1), radians(DEC1), 
                                            radians(RA2), radians(DEC2)))/2
        if radius > max_radius:
            max_radius = radius
        decs.append(DECc)
        ras.append(RAc)

    #print('All Racs = ', ras)
    #print('All Decs = ', decs)
    RAc_mean  = np.mean(ras)
    DECc_mean = np.mean(decs)
    
    print('DECc = ', DECc_mean)
    print('RAc  = ', RAc_mean)
    print('radius = %s arcmin' %(max_radius*60))
    return RAc_mean, DECc_mean, radius

#-----------------------------------------------------------------
def get_stars_coords(df, W, bdr, w, h, catalog):
    # Extracts data about stars from the pandas DataFrame
    # bdr is the width of the border to cut off the stars close to the edge

    pxpos = W.wcs_world2pix(np.stack((df['ra'], df['dec']), axis=-1),1)

    p_indx = np.where((pxpos[:,0]>bdr) & (pxpos[:,1]>bdr) & (pxpos[:,0]<w-bdr) & (pxpos[:,1]<h-bdr))

    xy = pxpos[p_indx]
    mags   = np.asarray(df['phot_g_mean_mag'])[p_indx]
    
    ra, dec = np.asarray(df.iloc[p_indx]['ra']),np.asarray(df.iloc[p_indx]['dec'])
    
    return xy, mags, ra, dec
#-------------------------------------------------------------------

def load_stars_from_Vizier(RAc, DECc, rc, low_mag, up_mag, filt, 
                          out_file=None, catalog=['NOMAD']):
    # Load stars from NOMAD catalog with astroquery,Vizir
    # low_mag < Vmag < up_mag
    #out_file = 'test.csv'
    Vizier.ROW_LIMIT = -1
    RAc = RAc * u.degree
    DECc = DECc * u.degree
    radius = rc * u.degree

    if (catalog[0] == 'NOMAD1') or (catalog[0] == 'NOMAD'):
        filters={'%smag' %filt: '%s : %s' %(low_mag, up_mag)}  
    elif catalog[0] == 'SDSS': 
        filters={'%smag' %filt: '%s : %s' %(low_mag, up_mag), 'cl': '6'}  
        #filters={'class': '=6'}  
    elif catalog[0] == 'PS1':
        filters={'%smag' %filt: '(>%s) & (<%s)' %(low_mag, up_mag)}  



    result = Vizier.query_region(coord.SkyCoord(ra=RAc, dec=DECc, frame='icrs'),radius=radius,catalog=catalog,  column_filters=filters )
    df = result[0].to_pandas()
    if out_file is None:
        pass
    else:
        out_file = out_file
        df.to_csv(out_file)
    return df
    
#--------------------------------------------------------------------
def get_stars_for_phot(df, W, bdr, w, h, catalog, filters):
    # Extracts data about stars from the pandas DataFrame
    # bdr is the width of the border to cut off the stars close to the edge


    if (catalog == 'NOMAD1') or (catalog == 'NOMAD'):
        pxpos = W.wcs_world2pix(np.stack((df['RAJ2000'], df['DEJ2000']), axis=-1),1)
        p_indx = np.where((pxpos[:,0]>bdr) & (pxpos[:,1]>bdr) & (pxpos[:,0]<w-bdr) & (pxpos[:,1]<h-bdr))
    
    elif catalog == 'SDSS': 
        pxpos = W.wcs_world2pix(np.stack((df['RA_ICRS'], df['DE_ICRS']), axis=-1),1)
        p_indx = np.where((pxpos[:,0]>bdr) & (pxpos[:,1]>bdr) & (pxpos[:,0]<w-bdr) & (pxpos[:,1]<h-bdr))
        #Bmag   = np.asarray(df['umag'])[p_indx]
        #Vmag   = np.asarray(df['gmag'])[p_indx]
        #Rmag   = np.asarray(df['rmag'])[p_indx]
    elif catalog == 'PS1':
        pxpos = W.wcs_world2pix(np.stack((df['RAJ2000'], df['DEJ2000']), axis=-1),1)
        p_indx = np.where((pxpos[:,0]>bdr) & (pxpos[:,1]>bdr) & (pxpos[:,0]<w-bdr) & (pxpos[:,1]<h-bdr))
        #Bmag   = np.asarray(df['gmag'])[p_indx]
        #Vmag   = np.asarray(df['rmag'])[p_indx]
        #Rmag   = np.asarray(df['imag'])[p_indx]
         
     
    Bmag   = np.asarray(df['%smag' %filters[0]])[p_indx]
    Vmag   = np.asarray(df['%smag' %filters[1]])[p_indx]
    Rmag   = np.asarray(df['%smag' %filters[2]])[p_indx]
     #pxpos = W.wcs_world2pix(np.stack((df['ra'], df['dec']), axis=-1),1)


    xy = pxpos[p_indx]
    #Bmag   = np.asarray(df['phot_bp_mean_mag'])[p_indx]#np.asarray(df['Bmag'])[p_indx]
    #Vmag   = np.asarray(df['phot_g_mean_mag'])[p_indx]#np.asarray(df['Vmag'])[p_indx]
    #Rmag   = np.asarray(df['phot_rp_mean_mag'])[p_indx]#np.asarray(df['Rmag'])[p_indx]
    
    
   
    # Deletes stars witn no data 
    inds = np.where(~ (np.isnan(Bmag) + np.isnan(Vmag) + np.isnan(Rmag)))
    xy = xy[inds]
    Bmag   = Bmag[inds]
    Vmag   = Vmag[inds]
    Rmag   = Rmag[inds]
        
    return xy, Bmag, Vmag, Rmag

#------------------------------------------------------------------------------
def plot_sky(ima, cf1=0.1, cf2=1.6, xy=None, asteroid=None):
    # Plot object image with cuts cf1, cf2:
    # vmin = mean - cf1*std
    # vmax = mean + cf2*std
    # xy - pix coord of the objectthat you want to additionally display on the frame

    cf1 = cf1
    cf2 = cf2

    fig, ax = plt.subplots(figsize=(3.0,3.0), dpi=300)
    ax.xaxis.set_tick_params(labelsize=5)
    ax.yaxis.set_tick_params(labelsize=5)
    plt.xlabel('X, pix', fontsize = 4)
    plt.ylabel('Y, pix', fontsize=4)

    plt.imshow(ima, vmin=np.median(ima)- cf1*np.std(ima) ,
           vmax=np.median(ima)  + cf2*np.std(ima), origin='lower', cmap='gray')
    if xy is None:
        pass
    else:
        for p in xy:
            plt.scatter(p[0], p[1], s=7, facecolors='none', edgecolors='g', linewidth=0.4)
        #plt.scatter(y_a, x_a, s=7, facecolors='none', edgecolors='r', linewidth=0.4)    
    plt.tight_layout()
    plt.show()

#------------------------------------------------------------------------------
    
def del_nearest(xs, ys, ap):
    goods = []
    i = 0
    ind = []
    for x, y in zip(xs, ys):
        dist = np.sort(np.hypot(xs - x, ys - y))[1]
        if dist > ap:
            goods.append([x, y])
            ind.append(i)
        i += 1
    return goods, ind    
#------------------------------------------------------------------------------

from astroquery.gaia import Gaia
import ssl
import pandas as pd
def get_stars_from_GAIA(RAc, DECc, limmag1, limmag2, radius, filename=None):
    # Similar to the function above, but extracts stars from GAIA

    limmag1,limmag2 = limmag1, limmag2
    ssl._create_default_https_context = ssl._create_unverified_context
    fov = radius # это размер области в градусах
    if filename is None:
        filename = 'Gaia_cur.csv'
    # список звезд сохраняется в файл BP_PSC_Gaia.csv
    job = Gaia.launch_job_async("SELECT * FROM gaiadr3.gaia_source \
        WHERE CONTAINS(POINT('ICRS',gaiadr3.gaia_source.ra,gaiadr3.gaia_source.dec),CIRCLE('ICRS',%f,%f,%f))=1\
        AND  phot_g_mean_mag>%f AND  phot_g_mean_mag<%f;"\
                %(RAc,DECc,fov,limmag1,limmag2), dump_to_file=True, output_format='csv',
                            output_file=filename)
    df = pd.read_csv(filename)
    return df
#---------------------------------------------------------------------------------

def echo(args):
    dic = args.__dict__

    for key in dic:
        print('%s: %s' %(key, dic[key]))
    return

#---------------------------------------------------------------------------------
def save(cols, header, out_dir):

    out_path = Path(out_dir, 'photometry_result.dat')
    with open(out_path, 'w') as out:
        print('#' + ' '.join(header), file=out)
        n_cols = len(cols)
        n_rows = len(cols[0])

        for i in range(n_rows):
            for j in range(n_cols):
                print('%8.3f' %(np.round(cols[j][i], 3)), end='    ', file=out)
            print('',file=out)
    return


