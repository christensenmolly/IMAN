import move_phot_sys    as mps
import star_photometry  as sph
import useful_functions as ufoo
import settings         

from pathlib import Path
from astropy.wcs import WCS
import subprocess
import numpy as np
import argparse
import os

def main(fnames, bdr=0, low_mag=14, up_mag=17, fwhm=9.0, k0=2, k1=2.5, k2=3, 
            path_to_mask=None, catalog='NOMAD', manual=True, exptime=None, filters='BVR', out_dir='.', app=None, out_file='final_image_corr.fits', show=False):

    sets = settings.Settings()

    RAc, DECc, rc = ufoo.get_image_center(fnames)
    if len(fnames) == 1:
        f = filters[0]
        filters = 3*filters
    else:
        f = filters[1]

    if catalog in sets.known_catalogs: 
        print('Getting data from the star catalogues. This may take some time...')
        df = ufoo.load_stars_from_Vizier(RAc, DECc, rc, low_mag=low_mag, up_mag=up_mag, filt=f,
                          out_file=Path(out_dir, 'input_catalog.csv'), catalog=[catalog])
    else:        
        print('I do not know such a catalog... Exit.')
        os._exit(1)

    mags   = []
    fwhms  = []
    coords = []
    first  = True
    
    for fname in fnames:
    
        data, header = ufoo.get_file(fname)
        gal = fname.stem #fname.split('-')[0].split('/')[-1]
        
        if path_to_mask is None:
            pass
        else:
            mask, _ = ufoo.get_file(path_to_mask)
            mask = mask >= 1
            data[mask] = float('nan')
        
        if first:
            h, w = data.shape
            W = WCS(header)
            xy, Bmag, Vmag, Rmag = sph.get_stars_for_phot(df, W, bdr, w, h, catalog, filters)
            print('Objects found: %i' %len(Bmag)) 
                
            if app is None:
                sph.create_ds9_apertures(xy, fwhm, k0, k1, k2)
                app = 'aper.reg'
            app = Path(app)
            
            if manual:
                #sph.grow_curve(data, xs, ys, gal)

                print('Open ds9')
                sph.call_ds9(fname, app)
                apertures, annulus_apertures, inds = sph.load_ds9_apertures(app)
                
                # !!!!
                xy   = xy[inds]
                Bmag = Bmag[inds]
                Vmag = Vmag[inds]
                Rmag = Rmag[inds]

            else:
                apertures, annulus_apertures, inds = sph.load_ds9_apertures(app)
                xy   = xy[inds]
                Bmag = Bmag[inds]
                Vmag = Vmag[inds]
                Rmag = Rmag[inds]
                
                #apertures         = None
                #annulus_apertures = None
                #inds              = np.arange(len(Bmag))
            first = False

        # mag, fwhm_, xy = sph.aperture(data, xy, fwhm=fwhm, k0=k0, k1=k1, k2=k2, apertures=apertures, annulus_apertures=annulus_apertures)
        mag, fwhm_, xy_ = sph.aperture(data, apertures=apertures, annulus_apertures=annulus_apertures)
        mags.append(mag)
        fwhms. append(fwhm_)
        coords.append(xy_)

    if len(fnames) == 1:

        Best  = mags[0]
        Bfwhm = fwhms[0]
        Bxy   = coords[0]

        m_fwhm = np.median(Bfwhm[~np.isnan(Bfwhm)])
        mad_fwhm  = sph.mad(Bfwhm[~np.isnan(Bfwhm)])
        print('FWHM:', m_fwhm, '+-', mad_fwhm)
        Binds = np.where(abs(Bfwhm - m_fwhm) <= 3*mad_fwhm)
        Bfwhm = Bfwhm[Binds] 


        Best = Best[Binds]
        Bmag = Bmag[Binds]
        Bxy   = Bxy[Binds]
        inds = inds[Binds]

        print('Objects extracted: %i' %len(Bmag))
        if len(Bmag) < 3:
            print('There are too few points to establish a dependency. Exit....')
            os._exit(1)
       
        apertures = apertures[Binds]
        annulus_apertures = annulus_apertures[Binds]
        

        #sph.plot_apers(data, apertures, annulus_apertures, out_dir)
        sph.update_ds9_apertures(region=app, inds=inds, outdir=out_dir)
        
        dm = sph.grow_curve(data, Bxy, gal, fwhm, k0, k2)
        #Best = Best + dm
        print('Delta m:', dm)

        xs = [x[0] for x in Bxy]
        ys = [x[1] for x in Bxy]
        if show:
            ufoo.save([inds, xs, ys,  Best, Bmag, Bfwhm], ['index', 'x_center', 'y_center', '%s_instrumental' %filters[0], '%s_catalog' %filters[0], 'FWHM'], out_dir) 
    
        mps.get_equals_solo(cat_B=Bmag, est_B=Best, fnameB=fnames[0], filt=filters[0], out_dir=out_dir, inds0=inds, out_file=out_file)
        return

    else:
        Best, Vest, Rest    = mags
        Bfwhm, Vfwhm, Rfwhm = fwhms
        Bxy, Vxy, Rxy       = coords
       
        
        Bm_fwhm = np.median(Bfwhm[~np.isnan(Bfwhm)])
        Bmad_fwhm  = sph.mad(Bfwhm[~np.isnan(Bfwhm)])

        print('FWHM B:', Bm_fwhm, '+-', Bmad_fwhm)

        Vm_fwhm = np.median(Vfwhm[~np.isnan(Vfwhm)])
        Vmad_fwhm  = sph.mad(Vfwhm[~np.isnan(Vfwhm)])
        
        print('FWHM V:', Vm_fwhm, '+-', Vmad_fwhm)

        Rm_fwhm = np.median(Rfwhm[~np.isnan(Rfwhm)])
        Rmad_fwhm  = sph.mad(Rfwhm[~np.isnan(Rfwhm)])

        print('FWHM R:', Rm_fwhm, '+-', Rmad_fwhm)

        Binds = np.where((abs(Bfwhm - Bm_fwhm) <= 3*Bmad_fwhm)*(abs(Vfwhm - Vm_fwhm) <= 3*Vmad_fwhm)*(abs(Rfwhm - Rm_fwhm) <= 3*Rmad_fwhm))
        Bfwhm = Bfwhm[Binds]
        Vfwhm = Vfwhm[Binds]
        Rfwhm = Rfwhm[Binds]

        

        #xy   = xy[Binds]
        #xs = [x[0] for x in xy]
        #ys = [x[1] for x in xy]

        Best = Best[Binds]
        Bmag = Bmag[Binds]
        Bxy  = Bxy[Binds]
        data, header = ufoo.get_file(fnames[0])
        dm = sph.grow_curve(data, Bxy, gal, fwhm, k0, k2)
        Best += dm
        print('Delta Bm:', dm)


        Vest = Vest[Binds]
        Vmag = Vmag[Binds]
        Vxy  = Vxy[Binds]
        data, header = ufoo.get_file(fnames[1])
        dm = sph.grow_curve(data, Vxy, gal, fwhm, k0, k2)
        Vest += dm
        print('Delta Vm:', dm)

        Rest = Rest[Binds]
        Rmag = Rmag[Binds]
        Rxy  = Rxy[Binds]
        data, header = ufoo.get_file(fnames[2])
        dm = sph.grow_curve(data, Rxy, gal, fwhm, k0, k2)
        Rest += dm
        print('Delta Rm:', dm)

        
        print('Objects extracted: %i' %len(Bmag))
        if len(Bmag) < 3:
            print('There are too few points to establish a dependency. Exit....')
            os._exit(1)
        apertures = apertures[Binds]
        annulus_apertures = annulus_apertures[Binds]

        #sph.plot_apers(data, apertures, annulus_apertures, out_dir)

        inds = inds[Binds]
        sph.update_ds9_apertures(region=app, inds=inds, outdir=out_dir)

        if exptime == None:
            pass
        else:
            Best += 2.5*np.log10(exptime)
            Vest += 2.5*np.log10(exptime)
            Rest += 2.5*np.log10(exptime)
        xs = [x[0] for x in np.mean([Bxy, Vxy, Rxy], axis=0)]
        ys = [x[1] for x in np.mean([Bxy, Vxy, Rxy], axis=0)]
        ufoo.save([inds, xs, ys,  Best, Bmag, Bfwhm, Vest, Vmag, Vfwhm,  Rest, Rmag, Rfwhm], ['index', 'x_center', 'y_center', '%s_instrumental' %filters[0], '%s_catalog' %filters[0], '%s_FWHM' %filters[0], '%s_instrumental' %filters[1], '%s_catalog' %filters[1], '%s_FWHM' %filters[1], '%s_instrumental' %filters[2], '%s_catalog' %filters[2], '%s_FWHM' %filters[2]], out_dir) 
    
        mps.get_equals(cat_B=Bmag, cat_V=Vmag, cat_R=Rmag, est_B=Best, est_V=Vest, est_R=Rest, fnameB=fnames[0], fnameV=fnames[1], fnameR=fnames[2], filters=filters, out_dir=out_dir, inds0=inds, out_file=out_file)
    
        return
    

def run(args):
   
    fnames = []
    file1 = Path(args.first_file)
    fnames.append(file1)
    
    file2 = args.second_file
    if file2 == 'N':
        pass
    else:
        file2 = Path(args.second_file)
        fnames.append(file2)
    file3 = args.third_file
    if file3 == 'N':
        pass
    else:
        file3 = Path(args.third_file)
        fnames.append(file3)
    #print(fnames) 
    #fnames  = [file1, file2, file3]
    catalog = args.catalog
    filters = args.filters
    low_mag = args.low_mag
    up_mag  = args.up_mag 
    bdr     = args.bdr
    fwhm    = args.fwhm
    k0      = args.k0
    k1      = args.k1
    k2      = args.k2
    manual  = args.manual
    p2m     = args.path_to_mask
    app     = args.aperture_file
    if p2m is None:
        pass
    else:
        p2m = Path(p2m)
    exptime = None
    out_dir = args.path_to_out_dir
    if out_dir is None:
        result_dir = Path('.')
    else:
        if os.path.exists(out_dir):
            out_dir = Path(args.path_to_out_dir)
        else:
            key = True
            while key:
                ans = input("No path %s. Do you want to create path? (y/n) \n" %out_dir)

                if (ans == "y"):
                    subprocess.run("mkdir -p %s" %out_dir, shell=True)  
                    key = not(key)
                    out_dir = Path(args.path_to_out_dir)
                elif (ans == "n"):
                    print("Exit...")
                    key = not(key)
                    return
    echo = True    
    print('===INPUT_PARAMETERS===')
    if echo:
        ufoo.echo(args)
    print('===START===')

    main(fnames=fnames, filters=filters, bdr=bdr, low_mag=low_mag, up_mag=up_mag, fwhm=fwhm, k0=k0, k1=k1, k2=k2, 
            path_to_mask=p2m, catalog=catalog, manual=manual, exptime=exptime, out_dir=out_dir, app=app)



#-------------------------------------------------------------------------------------------

#if __name__ == '__main__':
#    os.chdir('../')
#    gals = ['M100']#['UGC9560','NGC7743','NGC6181','NGC5660','NGC5430','NGC5371','NGC3583','NGC895' ,'M100']
#    for gal in gals:
#        main(['WMO/bkg_est_result/%s-B.fits' %gal, 'WMO/bkg_est_result/%s-V.fits' %gal,'WMO/bkg_est_result/%s-R.fits' %gal])

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    
    parser.add_argument("first_file", type=str, help='First fits file')
    parser.add_argument("second_file", type=str, help='Second fits file')
    parser.add_argument("third_file", type=str,  help='Third fits file') 
    parser.add_argument("-c", "--catalog", type=str, default="NOMAD", help='Catalog of standard stars. Currently available: "SDSS", "NOMAD", "PS1" ')
    parser.add_argument("-f", "--filters", type=str, default="BVR",  help='Enter three filters for each image, respectively. Enter together, for example, "BVR" for the "NOMAD" catalog')
    parser.add_argument("-b", "--bdr", type=int, default=100, help="Indent from the edge of the image [pix]")
    parser.add_argument("-lm", "--low_mag", type=float, default=14.0, help="Filtering stars by stellar magnitude. low_mag < V (g) < up_mag")
    parser.add_argument("-um", "--up_mag", type=float, default=17.0, help="Filtering stars by stellar magnitude. low_mag < V (g) < up_mag")
    parser.add_argument( "--fwhm", type=float, default=9.0, help="FWHM of stars [pix]")
    parser.add_argument("--k0", type=float, default=2.0, help="The radius of the aperture is k0*fwhm")
    parser.add_argument("--k1", type=float, default=2.5, help="The inner radius of the annulus is k1*fwhm")
    parser.add_argument("--k2", type=float, default=3.0, help="The outer radius of the annulus is k2*fwhm")
    parser.add_argument("--manual", action='store_true',  help='Enabling manual mode')
    parser.add_argument("-p2m", "--path_to_mask", default=None,  help='Path to mask file')
    parser.add_argument("-out", "--path_to_out_dir", type=str, default='.',  help='Path to output directory')
    parser.add_argument('-app', '--aperture_file', type=str, default=None, help='ds9.reg file with apertures')


    
    args = parser.parse_args()
    run(args)
