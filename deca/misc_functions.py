#!/usr/bin/env python
# Module with misc functions
# *****************************************************************
# **         DECA -- DECompoistion Analysis of galaxies          **
# **         Astronomical Observatory, Ghent University         **
# *****************************************************************

# Import standard packages
import os
import shutil
import sys
import math
import numpy as np
from scipy.interpolate import interp1d

# -----------------------------------------------------------------
# CLASS TO HIGHLLIGHT THE OUTPUT TEXT
class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'

    def disable(self):
        self.HEADER = ''
        self.OKBLUE = ''
        self.OKGREEN = ''
        self.WARNING = ''
        self.FAIL = ''
        self.ENDC = ''
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# FUNCTIONS FOR PRINTING MESSAGES
def OUT_FAIL(message):
    print(bcolors.FAIL + message + bcolors.ENDC)
    
def OUT_WARNING(message):
    print(bcolors.WARNING + message + bcolors.ENDC)
    
def OUT_HEADER(message):
    print(bcolors.HEADER + message + bcolors.ENDC)
    
def OUT_INFO(message):
    print(bcolors.OKGREEN + message + bcolors.ENDC)

def OUT_TITLE(message):
    print(bcolors.OKBLUE + message + bcolors.ENDC)
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# FUNCTIONS TO RETURN BANDS DECA CAN WORK WITH
def FILTERS_LIST(out=None):
  arr = ('FUV, 153', 'NUV, 227', 'u, 354', 'U, 365','B, 445','g, 475','V, 551','r, 622','R, 658','i, 763','I, 806','Z, 900','z, 905','J, 1220','H, 1630','Ks, 2159','K, 2190','W1, 3400','L, 3450','I1, 3600', 'I2, 4500', 'W2, 4600')
  if out==None:
    return arr
  elif out=='band':
    a = []
    for k in range(len(arr)):
      a.append(arr[k].split(',')[0])
    return a
  elif out=='wave':
    a = []
    for k in range(len(arr)):
      a.append(float(arr[k].split(', ')[-1]))
    return a    
# -----------------------------------------------------------------

        
# -----------------------------------------------------------------
# FUNCTION TO RETURN WAVELENGTH AND MSUN BY A GIVEN FILTER NAME
def filters(Filter):
        # from http://mips.as.arizona.edu/~cnaw/sun.html
        # B&M 1998
        # [Vega,AB]
        filter_names = ['FUV', 'NUV', 'u', 'U', 'B', 'g', 'V', 'r', 'R', 'i', 'I', 'Z', 'z', 'J', 'H', 'Ks', 'K', 'W1', 'L', 'I1', 'I2', 'W2']
        wavelength = [153., 227., 354., 365., 445., 475., 551., 622., 658., 763., 806., 900., 905., 1220., 1630., 2159., 2190., 3400., 3450., 3600., 4500., 4600]
        Msun_Vega = [13.97, 8.45, 5.46, 5.61, 5.48, 5.22, 4.78, 4.50, 4.42, 4.16, 4.08, 4.00, 4.01, 3.64, 3.32, 3.27, 3.28, 3.24, 3.25, 3.24, 3.27, 3.28]
        Msun_AB = [16.42, 10.31, 6.45, 6.32, 5.36, 5.14, 4.79, 4.65, 4.65, 4.54, 4.55, 4.52, 4.52, 4.57, 4.71, 5.18, 5.19, 5.0, 5.0, 4.9, 5.2, 5.3]

        if Filter in filter_names:
            ind = filter_names.index(Filter)
            if deca_setup.mag_system=='Vega':
                return wavelength[ind], Msun_Vega[ind]
            elif deca_setup.mag_system=='AB':
                return wavelength[ind], Msun_AB[ind]
        else:
            try:
                # If Filter is wavelength in nm
                w = float(Filter)
                if deca_setup.mag_system=='Vega':
                    f2 = interp1d(np.array(wavelength),np.array(Msun_Vega))
                    return w, f2(w)
                elif deca_setup.mag_system=='AB':
                    f2 = interp1d(np.array(wavelength),np.array(Msun_AB))
                    return w, f2(w)
            except:
                try:
                    # Could not find Msun
                    return float(Filter), None
                except:
                    OUT_FAIL( '%s is not in the list of wavelengths! Please specify correct filter!' % (Filter) )
                    exit()
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# FUNCTION TO READ THE HEADER OF AN INPUT FITS-FILE
def read_header(galaxy_image, mask_image=None, psf_image=None, model=None, log_txt = ''):
  if '.' not in model:   
    hdulist = pyfits.open(galaxy_image)
    prihdr = hdulist[0].header
    data = hdulist[0].data
    nx, ny = data.shape[1], data.shape[0]

    if os.path.exists(mask_image):
        hdulist_mask = pyfits.open(mask_image)
        data_mask = hdulist_mask[0].data

        # Create astropy mask
        mask_astropy = np.zeros_like(np.array(data_mask,dtype=float),dtype=bool)

        for k in range(ny):
            for i in range(nx):
                if data_mask[k,i]!=0.:
                    mask_astropy[k,i] = True

    else:
        mask_astropy = None

    scale,note = GET_RESOLUTION(galaxy_image)        
    scale = float(scale)

    if 'SCALE' in prihdr:
        scale,log_txt,status_tmp = READ_KEYWORD(prihdr, 'SCALE', default_value=scale, log_txt=log_txt)

    fwhm,log_txt,status_tmp = READ_KEYWORD(prihdr, 'FWHM', default_value=float('nan'), log_txt=log_txt)

    if np.isnan(fwhm):
      if psf_image!=None:
        # Fit it
        hdulist_psf = pyfits.open(psf_image)
        prihdr_psf = hdulist_psf[0].header
        data_psf = hdulist_psf[0].data  
        nx_psf, ny_psf = data_psf.shape[1], data_psf.shape[0]
        xc_psf = nx_psf/2.
        yc_psf = ny_psf/2.
        fwhm = psf_fit.main(psf_image, xc_psf, yc_psf, nx_psf, ny_psf, 'gauss')
        fwhm = float(fwhm)
        log_txt = log_txt + 'FWHM is %.3f pix\n' % (fwhm)
      else:
        OUT_FAIL('Nor FWHM in header nor PSF image was found! Set to 3.0.')
        log_txt = log_txt + 'Nor FWHM in header nor PSF image was found! Set to 3.0.\n'
        fwhm = float(fwhm)


    if psf_image!=None:
        hdulist_psf = pyfits.open(psf_image)
        prihdr_psf = hdulist_psf[0].header
        if "SAMPLING" in prihdr_psf:
            sampling = int(prihdr_psf["SAMPLING"])
        else:
            sampling = 1

        if "CONVOLUTION_BOX" in prihdr_psf:
            convolution_box = int(prihdr_psf["CONVOLUTION_BOX"])
        else:
            convolution_box = None           
    else:
        sampling = 1
        convolution_box = int(deca_setup.conv_size)
        
    gain,log_txt,status_tmp = READ_KEYWORD(prihdr, 'GAIN', default_value=7., log_txt=log_txt)
    read_out_noise,log_txt,status_tmp = READ_KEYWORD(prihdr, 'RDNOISE', default_value=1., log_txt=log_txt)
    ncombine,log_txt,status_tmp = READ_KEYWORD(prihdr, 'NCOMBINE', default_value=1, log_txt=log_txt)
    exptime,log_txt,status_tmp = READ_KEYWORD(prihdr, 'EXPTIME', default_value=1, log_txt=log_txt)
    m0,log_txt,status_tmp = READ_KEYWORD(prihdr, 'M0', default_value=None, log_txt=log_txt)
    Sky_level,log_txt,status_tmp = READ_KEYWORD(prihdr, 'SKY_LEVEL', default_value=0., log_txt=log_txt)
    SkySubtr,log_txt,status_tmp = READ_KEYWORD(prihdr, 'SKY_SUBTR', default_value=1, log_txt=log_txt)
    xc,log_txt,status_tmp = READ_KEYWORD(prihdr,'XC', default_value=nx/2., log_txt=log_txt)
    yc,log_txt,status_tmp = READ_KEYWORD(prihdr,'YC', default_value=ny/2., log_txt=log_txt)
    std,log_txt,status_tmp = READ_KEYWORD(prihdr,'SKY_STD',default_value=float('nan'), log_txt=log_txt)
    if np.isnan(std):
          I_mean, I_median, std = sigma_clipped_stats(data, mask=mask_astropy, sigma=3.0, iters=5)

    read_out_noise = std/gain

    hdulist.close()
    observation_info_non_full = [nx,ny,m0,scale,gain,ncombine,exptime,read_out_noise,fwhm,Sky_level,SkySubtr,sampling,convolution_box]

    return observation_info_non_full,xc,yc,log_txt
  else:
    hdulist = pyfits.open(galaxy_image)
    prihdr = hdulist[0].header
    data = hdulist[0].data
    ny,nx = np.shape(data)
    m0 = 20.
    scale = 1.
    gain = 4. 
    ncombine = 1
    exptime = 1.
    read_out_noise = 10.
    fwhm = 3.
    Sky_level = 0.
    SkySubtr = 1.
    sampling = 1
    convolution_box = 51
    xc = nx/2.
    yc = ny/2.


    ff = open(model, 'r')
    for line in ff:
        if 'J)' in line and '# Magnitude photometric zeropoint' in line:
            MagZP = float(line.split()[1])
            if 'EXPTIME' in prihdr:
                exptime = float(prihdr['EXPTIME'])
                m0 = MagZP + 2.5*log10(exptime)
            else:
                m0 = MagZP
        if 'K)' in line and '# Plate scale (dx dy)   [arcsec per pixel]' in line:
            scale = float(line.split()[1])
    ff.close()
    return [nx,ny,m0,scale,gain,ncombine,exptime,read_out_noise,fwhm,Sky_level,SkySubtr,sampling,convolution_box],xc,yc,log_txt            
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# # FUNCTION TO CONVERT UNITS             
def units_converter(value,par_units,input_geom_units,input_SB_units,input_lum_units,output_geom_units,
                    output_SB_units,output_lum_units,m0,arcsec_per_pixel,Distance,kpc_per_arcsec,Filter,Aext=0.,Kcorr=0.):
  Msun = filters(Filter)
  
  if par_units=='geom':
        if input_geom_units=='pix' and input_geom_units!=output_geom_units:
          if output_geom_units=='arcsec':
            value = value * arcsec_per_pixel
          if output_geom_units=='kpc':
            #print(value, arcsec_per_pixel, kpc_per_arcsec)
            value = value * arcsec_per_pixel * kpc_per_arcsec
            
        if input_geom_units=='arcsec' and input_geom_units!=output_geom_units:
          if output_geom_units=='pix':
            value = value / arcsec_per_pixel
          if output_geom_units=='kpc':
            value = value * kpc_per_arcsec    
            
        if input_geom_units=='kpc' and input_geom_units!=output_geom_units:
          if output_geom_units=='pix':
            value = value / kpc_per_arcsec / arcsec_per_pixel
          if output_geom_units=='arcsec':
            value = value / kpc_per_arcsec
            
  elif par_units=='SB':
        if input_SB_units=='mag/arcsec2' and input_SB_units!=output_SB_units:
          if output_SB_units=='ADU/pix2':
            value = 10**(0.4*(m0 - value + 5.*log10(arcsec_per_pixel) + Aext + Kcorr))
          if output_SB_units=='Lsun/pc2':
            value = 10**(0.4*(Msun - value + 21.572 + Aext + Kcorr))
            
        if input_SB_units=='ADU/pix2' and input_SB_units!=output_SB_units:
          if output_SB_units=='mag/arcsec2':
            value = m0 -2.5*log10(value) + 5.*log10(arcsec_per_pixel) - Aext - Kcorr
          if output_SB_units=='Lsun/pc2':
            value = 10**(0.4*(Msun+21.572-m0 -2.5*log10(value) + 5.*log10(arcsec_per_pixel) - Aext - Kcorr))
        
        if input_SB_units=='Lsun/pc2' and input_SB_units!=output_SB_units:
          if output_SB_units=='mag/arcsec2':
            value = Msun + 21.572 - 2.5*log10(value) - Aext - Kcorr
          if output_SB_units=='ADU/pix2':
            value = 10**(0.4*(m0 - Msun + 21.572 - 2.5*log10(value) - Aext - Kcorr + 5.*log10(arcsec_per_pixel)))

    
  elif par_units=='lum':
        if input_lum_units=='mag' and input_lum_units!=output_lum_units:
          if output_lum_units=='ADU':
            value = 10**(0.4*(m0 - value + Aext + Kcorr))
          if output_lum_units=='Lsun':
            value = value - 5.*log10(Distance) - 25.
            value = 10**(0.4*(Msun - value))
            
        if input_lum_units=='ADU' and input_lum_units!=output_lum_units:
          if output_lum_units=='mag':
            value = m0 -2.5*log10(value) - Aext - Kcorr
          if output_lum_units=='Lsun':
            value = m0 -2.5*log10(value) - Aext - Kcorr
            value = value - 5.*log10(Distance) - 25.
            value = 10**(0.4*(Msun - value))
            
        if input_lum_units=='Lsun' and input_lum_units!=output_lum_units:
          if output_lum_units=='ADU':
            value = Msun - 2.5*log10(value)
            value = value + 5.*log10(Distance) + 25.
            value = 10**(0.4*(m0 - value + Aext + Kcorr))
          if output_lum_units=='mag':    
            value = Msun - 2.5*log10(value)
            value = value + 5.*log10(Distance) + 25.

  else:
        value = value
  return value
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# FUNCTION TO CHECK IF FILE EXISTS AND GET PATH TO IT
def CHECK_FILE_EXISTANCE(filename, log_text):
  import httplib2
  from urllib.parse import urlparse
  if filename!='none':
    try:
        '''
        request = requests.get(filename)
        if request.status_code == 200:
            return filename
        else:
            OUT_FAIL('%s does not exist! Exiting ...' % (filename))
            exit(1)
        '''
        p = urlparse(filename)
        conn = httplib2.HTTPConnection(p.netloc)
        conn.request('HEAD', p.path)
        resp = conn.getresponse()
        if resp.status < 400:
            return filename,0,log_text
        else:
            OUT_FAIL('%s does not exist!' % (filename))
            log_text = log_text + '%s does not exist!\n' % (filename)
            return None,1,log_text            
    except:
        if os.path.exists(filename):
            filename = os.path.abspath(filename)
            return filename,0,log_text
        else:
            OUT_FAIL('%s does not exist!' % (filename))
            log_text = log_text + '%s does not exist!\n' % (filename)
            return None,1,log_text  
  else:
      return 'none',0,log_text
# -----------------------------------------------------------------
  
  
# -----------------------------------------------------------------
# FUNCTION TO FIND TYPE OF THE GIVEN COMPRESSED FILE
def FILE_TYPE(filename):
    magic_dict = {
        "\x1f\x8b\x08": "gz",
        "\x42\x5a\x68": "bz2",
        "\x50\x4b\x03\x04": "zip"
        }

    max_len = max(len(x) for x in magic_dict)
    with open(filename) as f:
        file_start = f.read(max_len)
    for magic, filetype in magic_dict.items():
        if file_start.startswith(magic):
            return filetype
    return "no match"
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# FUNCTION TO DECOMPRESS FILE
def UNCOMPRESS_COPY_FILE(filename, outFilePath):
    import urllib2
    import StringIO
    import gzip

    # Download the file and save it
    downfile = os.path.dirname(os.path.abspath(outFilePath)) + '/' + filename.split('/')[-1]  

    #response = urllib2.urlopen(filename)
    
    opener = urllib2.build_opener(urllib2.ProxyHandler({}))
    #opener.open('http://anything.com/').read()


    with open(downfile,'wb') as output:
        #output.write(response.read())
        output.write(opener.open(filename).read())

    OUT_INFO('File %s has been successfully downloaded!'% (filename))
    
    filetype = FILE_TYPE(downfile)
    if filetype=='zip':
        import zipfile
        #zip_ref = zipfile.ZipFile(downfile, 'r')        
        directory_to_extract_to = os.path.dirname(os.path.abspath(outFilePath))

        zip_ref.extractall(directory_to_extract_to)
        zip_ref.close()
        shutil.move(directory_to_extract_to+'/'+outFilePath.split('/')[-1], outFilePath)

    elif filetype=='gz':
        try:
            f=gzip.open(downfile,'rb')
            compressedFile = StringIO.StringIO(f.read())
            decompressedFile = gzip.GzipFile(fileobj=compressedFile)

            with open(outFilePath, 'w') as outfile:
                outfile.write(decompressedFile.read())        
        except:
            os.rename(downfile, outFilePath)
        f.close()
    elif filetype=='bz2':
        import bz2
        zipfile = bz2.BZ2File(downfile) # open the file
        data = zipfile.read() # get the decompressed data
        open(outFilePath, 'wb').write(data) # write an uncompressed file
    else:
        os.rename(downfile, outFilePath)
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# FUNCTION TO COPY IMAGES FROM THE GIVEN DIRECTORY TO A NEW ONE
def COPY_IMAGES(IMAGES, directory, band=None):
    new_images = []
    status = 0
    for k in range(len(IMAGES)):
      if band==None:
        new_image = IMAGES[k].split('/')[-1]
        if '.gz' in new_image or '.bz2' in new_image or '.zip' in new_image:
            new_image = os.path.splitext(new_image)[0]
            
      else:
        new_image = IMAGES[k].split('/')[-1].split('.fits')[0]+'_'+band+'.fits'

      if IMAGES[k]!='none':# and os.path.exists(IMAGES[k]):   # Copying only if the input and output files are not identical! 
        if not os.path.exists(directory+new_image):     
          try:
            UNCOMPRESS_COPY_FILE(IMAGES[k], directory+new_image)
          except:
              try:
                shutil.copy(IMAGES[k], directory+new_image)
              except:
                status = 2
        
        
        #else:
        #    if not filecmp.cmp(IMAGES[k], directory+new_image):
        #       shutil.copy(IMAGES[k], directory+new_image) 
      new_images.append(new_image)

    return new_images, status
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# FUNCTION TO GET PIXEL SCALE OF THE GIVEN IMAGE
def GET_RESOLUTION(input_image):
    from astropy.io import fits
    hdulist = fits.open(input_image)#, ignore_missing_end=True)
    header = hdulist[0].header

    pixelscale = 1.
    note = ' '

    if 'PIXSCALE' in header:

        pixelscale = header['PIXSCALE']

    elif 'SECPIX' in header:

        pixelscale = header['SECPIX']

    elif 'PFOV' in header:

        pixelscale = header['PFOV']

    elif 'CD1_1' in header and 'CD1_2' in header:

        pixelscale = math.sqrt(header['CD1_1']**2 + header['CD1_2']**2 ) * 3600.0

    elif 'CD1_1' in header:

        pixelscale = abs(header['CD1_1']) * 3600.0

    elif 'CDELT1' in header:

        pixelscale = abs(header['CDELT1']) * 3600.0

    else:

        OUT_WARNING("Could not determine the pixel scale from the image header. Set to 1 pix/arcsec.")
        note = '*'

    # Return the pixel scale (in arcseconds)
    hdulist.close()
    return str(pixelscale),note
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# FUNCTION FOR READING KEYWORDS OF THE HEADER
def READ_KEYWORD(header, keyw, default_value=None, log_txt=''):
    if keyw in header:
      return header[keyw],log_txt,0
    else:
      if default_value!=None:
        log_txt = log_txt + 'This KEYWORD [%s] is not found in the input fits file with the galaxy image! Set to %f.\n' % (keyw, default_value)
        return default_value,log_txt,0
      else:
        OUT_FAIL('This KEYWORD [%s] is not found in the input fits file with the galaxy image!' % (keyw))
        log_txt = log_txt + 'This KEYWORD [%s] is not found in the input fits file with the galaxy image!' % (keyw)
        return None,log_txt,3
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# FUNCTION FOR EXITING
def EXIT_DECA():
    import deca_setup
    if deca_setup.stop_deca==True:
        exit()
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# OLD FUNCTION TO READ IN REGION FILE WITH OBJECTS
def READ_INFO_REGION_OBJECTS(line):
    name = 'none'
    model = 'none'
    redshift = float(nan)
    Aext = float(nan)
    Kcorr = float(nan)
    if 'text={' in line:
        pars = line.split('text={')[-1].split('}')[0].split(',')
        try:
            name = pars[0]
            model = pars[1]
            redshift = float(pars[2])
            Aext = float(pars[3])
            Kcorr = float(pars[4])
        except:
            z=1
    return [name,model,redshift,Aext,Kcorr]
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# FUNCTION TO READ IN REGION FILE WITH OBJECTS
def READ_REGION_OBJECTS(file_with_galaxies, input_image):
  # Read in the region file:
  f = open(file_with_galaxies, "r")
  
  # Determine DS9 system: 
  system = 'image'
  for line in f:
    if 'image' in line:
      system = 'image'
    if 'fk5' in line:
      system = 'fk5'
      hdulist = pyfits.open(input_image)
      w = wcs.WCS(hdulist[0].header)
  f.close()

  f = open(file_with_galaxies, "r")  
  xc = []; yc = []; ellA = []; ellB = []; ellPA = []
  Name = []; Model = []; Redshift = []; AExt = []; KCorr = []
  xc_stars = []; yc_stars = []
  for line in f:
    if system=='image':
        if "ellipse(" in line:
            params = line.split(",")
            xcc = float(params[0].split('(')[1])
            ycc = float(params[1])
            sma = float(params[2])
            smb = float(params[3])
            PA = float(params[4].split(')')[0])
            if sma < smb:
                sma, smb = smb, sma
                PA += 90
                
            [[xl,yl],[xr,yr]] = ellipse_borders([xcc,ycc],sma,smb,PA)
            if check_belong(xl,yl,xr,yr,nx,ny)==True:
                xc.append(xcc)
                yc.append(ycc)
                ellA.append(sma)
                ellB.append(smb)
                ellPA.append(PA)
                
            [name,model,redshift,Aext,Kcorr] = READ_INFO_REGION_OBJECTS(line)
            Name.append(name)
            Model.append(model)
            Redshift.append(redshift)
            AExt.append(Aext)
            KCorr.append(Kcorr)
            
        if "circle(" in line:
            params = line.split(",")
            xcc = float(params[0].split('(')[1])
            ycc = float(params[1])
            r = float(params[2].split(')')[0])
            [[xl,yl],[xr,yr]] = ellipse_borders([xcc,ycc],r,r,0.)
            if check_belong(xl,yl,xr,yr,nx,ny)==True:
                xc.append(xcc)
                yc.append(ycc)
                ellA.append(r)
                ellB.append(r)
                ellPA.append(0.)
            [name,model,redshift,Aext,Kcorr] = READ_INFO_REGION_OBJECTS(line)
            Name.append(name)
            Model.append(model)
            Redshift.append(redshift)
            AExt.append(Aext)
            KCorr.append(Kcorr)
            
        if "point(" in line:
            xc_stars.append(float(line.split('(')[1].split(',')[0]))
            yc_stars.append(float(line.split(',')[1].split(')')[0]))
            
    if system=='fk5':
        if "ellipse(" in line:
            params = line.split(",")
            xcc_fk5 = float(params[0].split('(')[1])    # RA, in deg
            ycc_fk5 = float(params[1])                  # DEC, in deg
            sma_fk5 = float(params[2])                  # Sma, in arcsec
            smb_fk5 = float(params[3])                  # Smb, in arcsec
            PA_fk5 = float(params[4].split(')')[0])     # HyperLeda position angle of the major axis of the isophote 25 mag/arcsec**2 in the B-band for galaxies. It is counted from the North (pa=0) toward East between 0 and 180.
            if sma_fk5 < smb_fk5:
                sma_fk5, smb_fk5 = smb_fk5, sma_fk5
                PA_fk5 += 90
                
            # Convert fk5 to pixels
            [[xcc,ycc]] = w.wcs_world2pix([[xcc_fk5,ycc_fk5]],1)
            skycoord = SkyCoord(ra=xcc_fk5*u.degree, dec=ycc_fk5*u.degree, frame=FK5)
            PA_x_N = angle_at_skycoord(skycoord, w)
            PA = PA_x_N - PA_fk5
            sma = sma/scale
            smb = smb/scale
            
            [[xl,yl],[xr,yr]] = ellipse_borders([xcc,ycc],sma,smb,PA)
            if check_belong(xl,yl,xr,yr,nx,ny)==True:
                xc.append(xcc)
                yc.append(ycc)
                ellA.append(sma)
                ellB.append(smb)
                ellPA.append(PA)
            [name,model,redshift,Aext,Kcorr] = READ_INFO_REGION_OBJECTS(line)
            Name.append(name)
            Model.append(model)
            Redshift.append(redshift)
            AExt.append(Aext)
            KCorr.append(Kcorr)
            
        if "circle(" in line:
            params = line.split(",")
            xcc_fk5 = float(params[0].split('(')[1])
            ycc_fk5 = float(params[1])
            r_fk5 = float(params[2].split(')')[0])

            # Convert fk5 to pixels
            [[xcc,ycc]] = w.wcs_world2pix([[xcc_fk5,ycc_fk5]],1)
            r = r_fk5/scale            
            
            [[xl,yl],[xr,yr]] = ellipse_borders([xcc,ycc],r,r,0.)
            if check_belong(xl,yl,xr,yr,nx,ny)==True:
                xc.append(xcc)
                yc.append(ycc)
                ellA.append(r)
                ellB.append(r)
                ellPA.append(0.)
            [name,model,redshift,Aext,Kcorr] = READ_INFO_REGION_OBJECTS(line)
            Name.append(name)
            Model.append(model)
            Redshift.append(redshift)
            AExt.append(Aext)
            KCorr.append(Kcorr)
            
        if "point(" in line:
            xc_star_fk5 = float(line.split('(')[1].split(',')[0])
            yc_star_fk5 = float(line.split(',')[1].split(')')[0])
            [[xc_star,yc_star]] = w.wcs_world2pix([[xc_star_fk5,yc_star_fk5]],1)
            
            xc_stars.append(xc_star)
            yc_stars.append(yc_star)

  f.close()
  return [xc,yc,ellA,ellB,ellPA,Name,Model,Redshift,AExt,KCorr],[xc_star,yc_star]
# -----------------------------------------------------------------


# -----------------------------------------------------------------
def remove_file(file):
    if os.path.exists(file):
        os.remove(file)
