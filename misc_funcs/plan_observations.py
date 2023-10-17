# Observing Run Preparation - GROWTH Astronomy School 2020
# https://www.eso.org/sci/observing/tools.html
# http://www.astro.wisc.edu/~moravec/data_reduction.html
# https://arxiv.org/pdf/1812.07963.pdf
# https://www.astro.rug.nl/~ndouglas/teaching/JAFFE/reduce.htm

# http://astronomy.nmsu.edu:8000/apo-wiki
# https://mwcraig.github.io/ccd-as-book/00-00-Preface.html   !!! PACKAGE
# https://image-analysis.readthedocs.io/en/latest/index.html
# http://astro.corlan.net/gcx/html/node7.html
# https://faculty.virginia.edu/ASTR5110/lectures/detectors/detectors_red.html !!! Good description of the data reduction
# EXAMPLE: python3 plan_observations.py Filaments_for_observations.dat --d=' '
# EXAMPLE: python3 ~/MyGit/IMAN/misc_funcs/plan_observations.py stripe82_prg_cands.dat 2023-07-15 01:00:00 2023-07-15 05:40:00 --d=' ' --s desi

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import EarthLocation
from astropy.coordinates import AltAz
from astropy.coordinates import ICRS, Galactic, FK4, FK5
from astropy.time import Time
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
from astropy.time import TimeDelta
from astroplan import Observer
from pytz import all_timezones_set, common_timezones_set
from astroplan import FixedTarget
from astroplan.plots import dark_style_sheet, plot_airmass
from astroplan.plots import plot_sky
from astroplan.plots import plot_finder_image
from astroplan import FixedTarget
from astropy.coordinates import Angle
import argparse

import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning) 
import pytz
import time
import os
import sys
import shutil

LOCAL_DIR = "/misc_funcs"
IMAN_DIR = os.path.dirname(__file__).rpartition(LOCAL_DIR)[0]


sys.path.append(os.path.join(IMAN_DIR, 'detect_objects'))
sys.path.append(os.path.join(IMAN_DIR, 'misc_funcs'))

import download_sdss_jpg
import download_panstarrs_jpg
import download_legacy_jpg_new
import read_data
import get_galaxy_center_ned

def get_time_zone(t, lat, lon):
    from timezonefinder import TimezoneFinder
    t_find = TimezoneFinder()
    timezone_info = t_find.certain_timezone_at(lng=lon, lat=lat)

    #day = datetime(2021, 4, 24)
    day = datetime.strptime(t, "%Y-%m-%d %H:%M:%S")

    #offset = day.utcoffset()#.total_seconds()/60/60

    tz = pytz.timezone(timezone_info)
    offset = tz.utcoffset(day).total_seconds()/3600.
    #print('Time zone:', offset)
    return int(offset)
    '''
    if timezone_info is None:
        print("Could not determine the time zone")
    else:
        # Display the current time in that time zone
        timezone = pytz.timezone(timezone_info)
        print(timezone_info)
    '''




def check_winter_summer_time(t):
    time1 = time.strptime(t, "%Y-%m-%d %H:%M:%S")
    ad_hour = time.localtime(time.mktime(time1)).tm_isdst
    return ad_hour

#check_winter_summer_time('2022-08-24 20:10:00')
#exit()


def get_ha(time, location, targets):
    lst = time.sidereal_time('apparent', longitude=location.lon)
    ha = lst - targets.ra
    
    # return negative values for rising targets
    return (ha + 12 * u.hourangle) % (24 * u.hourangle) - 12 * u.hourangle

def make_visibility_plot(time, location, targets, names):
    ha = get_ha(time, location, targets)
    for k in range(len(names)):
        plt.plot(ha[k].hour, targets[k].dec.deg, 'o', label=names[k])
    
    lst = time.sidereal_time('apparent', longitude=location.lon)
    #print('Local siderial time:',lst)

    az = np.linspace(0.1, 359.8, 200) * u.deg
    for alt in [45, 30, 20]:
        # convert (Az, Alt) to (HA, Dec)
        line = SkyCoord(az, alt * u.deg, frame=AltAz(obstime=time, location=location))
        radec = line.transform_to(ICRS)
        ha = get_ha(time, location, radec)
        X = 1./np.sin(np.radians(alt+244./(165.+47.*alt**(1.1))))
        #X = 1/np.sin(np.radians(alt))

        plt.plot(np.array(ha.hour), np.array(radec.dec.deg), '--', label='Alt=%i; X=%.2f' % (alt, X))
    plt.legend()
    plt.grid(True)
    plt.xlabel('Hour Angle (%s-R.A.)' % (str(lst)))
    plt.ylabel('Declination [deg]')
    plt.title("Visibility from APO at %s" % (str(time)))
    plt.show() 

def read_targets(input_file, delimiter=',', header_line=0, units_line=None, skip_lines = [], comment='#'):
    f = open(input_file, "r")
    cul_num = len(f.readlines()[header_line].split(delimiter))

    if cul_num == 1:
        # There is only one column - the name of the object
        names = np.loadtxt(input_file, usecols=[0], skiprows = 1, dtype=str)

        ra = []; dec = []
        for name in names:
            coords = get_galaxy_center_ned.get_coords_by_name(name)
            ra.append(coords[0])
            dec.append(coords[1])
        targets = SkyCoord(ra, dec, frame=FK5, unit="deg")

        return targets,names

    # Main columns should be: name, ra, dec
    try:
        target_list,units = read_data.main(input_file, delimiter=delimiter, header_line=header_line, units_line=units_line, skip_lines = skip_lines, comment=comment)


        if ':' in str(target_list['ra'][0]):
            targets = SkyCoord(target_list['ra'], target_list['dec'], unit=(u.hourangle, u.deg))
        elif 'h' in str(target_list['ra'][0]):
            targets = SkyCoord(target_list['ra'], target_list['dec'], frame=FK5)
        else:
            targets = SkyCoord(target_list['ra'], target_list['dec'], frame=FK5, unit="deg")

    except:
        dtype = []
        dtype.append( ('name', 'U18') )
        dtype.append( ('ra', 'U18') )
        dtype.append( ('dec', 'U18') )
        dtype.append( ('epoch', 'U18') )
        target_list = np.genfromtxt(input_file, dtype=dtype)



        targets = SkyCoord(target_list['ra'], target_list['dec'], frame=FK5, unit='hour, deg')
    return targets,target_list['name']






def main(table, targets, names, local_time_start='2023-05-17 00:05:00', local_time_end='2023-05-17 01:04:00', lon=-105.819644, lat=32.780359, elevation=2788., observatory='APO', bad_if_airmass_is_larger=2.0, pressure_bar=1.0, temperature_C=0.0, verbosity=False, plot_visibility=False, survey=None):

    Time_zone = get_time_zone(local_time_start, lat, lon) # This includes saving daytime hour!!!

    #ad_hour = check_winter_summer_time(local_time_start)
    #Time_zone = -7 + ad_hour #Time_zone = -7 #??? for APO ! 6 - summer time, 7 - winter time
    #Time_zone = -6

    # Set location and site conditions
    location = EarthLocation.from_geodetic(lon*u.deg, lat*u.deg, elevation*u.m)
    OBSERVATORY = Observer(location=location, name=observatory, timezone="Etc/GMT%s" % (str(Time_zone)), elevation=elevation*u.m, pressure=pressure_bar*u.bar, temperature=temperature_C*u.deg_C)


    # Convert local time to UTC
    dt = TimeDelta(3600.0*Time_zone, format='sec')
    t_start = Time(local_time_start, format='iso', scale='utc', location=location) - dt
    t_end = Time(local_time_end, format='iso', scale='utc', location=location) - dt   


    sunset_tonight = OBSERVATORY.sun_set_time(t_start, which='nearest')
    sunrise_tonight = OBSERVATORY.sun_rise_time(t_end, which='nearest')
    moon_rise = OBSERVATORY.moon_rise_time(t_start)
    moon_set = OBSERVATORY.moon_set_time(t_end)    
    
    if verbosity:
        print('\tUTC sunset time:', sunset_tonight.iso)
        print('\tUTC sunrise time:', sunrise_tonight.iso)

        print('\tLocal sunset time:', Time(sunset_tonight.iso, format='iso', scale='utc', location=location) + dt)
        print('\tLocal sunrise time:', Time(sunrise_tonight.iso, format='iso', scale='utc', location=location) + dt)

        print('\tUTC moon rise time:', moon_rise.iso)
        print('\tUTC moon set time:', moon_set.iso)

        print('\tLocal moon rise time:', Time(moon_rise.iso, format='iso', scale='utc', location=location) + dt)
        print('\tLocal moon set time:', Time(moon_set.iso, format='iso', scale='utc', location=location) + dt)

        print('\n\tLimits on our observing window:')
        print('\tApproximate range for dark time (after evening twilight and before morning twilight) UTC: %s --- %s' % (Time(sunset_tonight.iso, format='iso', scale='utc', location=location) + TimeDelta(3600.0*1.5, format='sec'), Time(sunrise_tonight.iso, format='iso', scale='utc', location=location) - TimeDelta(3600.0*1.5, format='sec')))

        print('\tApproximate range for dark time (after evening twilight and before morning twilight) Local time: %s --- %s' % (Time(sunset_tonight.iso, format='iso', scale='utc', location=location) + TimeDelta(3600.0*1.5, format='sec') + dt, Time(sunrise_tonight.iso, format='iso', scale='utc', location=location) - TimeDelta(3600.0*1.5, format='sec') + dt))


    print('\n======================================================================================================')
    print('Night %s - %s: OBJECTS WITH AIRMASS<%.2f' % (local_time_start, local_time_end, bad_if_airmass_is_larger))
    print('SUN set-rise: %s - %s' % (Time(sunset_tonight.iso, format='iso', scale='utc', location=location) + dt, Time(sunrise_tonight.iso, format='iso', scale='utc', location=location) + dt))
    print('MOON rise-set: %s - %s' % (Time(moon_rise.iso, format='iso', scale='utc', location=location) + dt, Time(moon_set.iso, format='iso', scale='utc', location=location) + dt))
    print('------------------------------------------------------------------------------------------------------')
    print('%s | %s | %s| %s|%s| %s |%s' % ('Ind'.ljust(3),'Object'.ljust(20), 'OBS_start'.ljust(24), 'OBS_end'.ljust(24), 'H_start'.ljust(7), 'H_end'.ljust(5), 'MOON_IL'.ljust(7) ))
    print('------------------------------------------------------------------------------------------------------')

    
    if plot_visibility:
        make_visibility_plot(t_start, location, targets, names)

    LST_start = t_start.sidereal_time('apparent', longitude=location.lon)
    LST_end = t_end.sidereal_time('apparent', longitude=location.lon)
    deltat = LST_end - LST_start
    times = LST_start + deltat * np.linspace(0.,1.,100)

    #print('\n================================================')
    #print("\n\tObjects which have airmass<%.1f over the observational period:" % (bad_if_airmass_is_larger))
    observe_time = np.arange(t_start,t_end,TimeDelta(60., format='sec')) # step = 1 minute

    sel_targets = []
    for i in range(len(targets)):
        target = targets[i]

        #delta_t = t_end - t_start
        #observe_time = t_start + delta_t*np.linspace(0, 1, 75)
        #print(t_start+TimeDelta(60., format='sec'))
        

        
        coordinates = SkyCoord(target.ra.deg*u.degree, target.dec.deg*u.degree, frame='fk5')
        TARGET = FixedTarget(coord=coordinates)  
        zenith_angle = Angle(OBSERVATORY.altaz(observe_time, TARGET).zen, unit=u.deg)
        zenith_angle = np.array(zenith_angle.hour*15., float)
        alt = 90.-zenith_angle

        airmasses = 1./np.sin(np.radians(alt+244./(165.+47.*alt**(1.1)))) # Pickering (2002): http://www.dioi.org/vols/wc0.pdf 
        #airmasses = (np.cos(np.radians(zenith_angle)) + 0.025*np.exp(-11.*np.cos(np.radians(zenith_angle))))**(-1) 
        inds = []
        for k in range(len(airmasses)):
            if airmasses[k]<=bad_if_airmass_is_larger:
                inds.append(k)
        

        if inds!=[]:
            inds = np.array(inds)
            ttime = np.array(observe_time)[inds]
            aairmasses = airmasses[inds]
            Moon_illumination = np.mean(OBSERVATORY.moon_illumination(ttime))
            #for tttime in ttime:
            #    print(OBSERVATORY.moon_illumination(tttime))
            '''
            if OBSERVATORY.moon_illumination(np.min(ttime)) > bad_if_moon_illumination_is_larger:
                print('\t',names[i],Time(np.min(ttime), format='iso', scale='utc', location=location) + dt, Time(np.max(ttime), format='iso', scale='utc', location=location) + dt, 'WARNING: MOON ILLUMINATION is LARGER THAN %.1f' % (bad_if_moon_illumination_is_larger))
            else:
                print('\t',names[i],Time(np.min(ttime), format='iso', scale='utc', location=location) + dt, Time(np.max(ttime), format='iso', scale='utc', location=location) + dt)
            '''
                
        
            # Alterantive estimate (they're roughly consistent)
            '''
            Alt = []
            for time in times:
                    Alt.append(float(np.degrees(np.arcsin(np.sin(location.lat.radian) *
                                    np.sin(target.dec) +
                                    np.cos(location.lat.radian) *
                                    np.cos(target.dec) *
                                    np.cos(time.radian - target.ra.radian))).value))
            print(names[i],'Alt(deg) = %.2f (start obs), %.2f (end obs), %.2f (max ALT), %.2f (min ALT)' % (Alt[0], Alt[-1], np.max(Alt), np.min(Alt)))
            '''
            #print('\t',names[i],'Alt(deg) = %.2f (start obs), %.2f (end obs), %.2f (max ALT), %.2f (min ALT)' % (alt[0], alt[-1], np.max(alt), np.min(alt)))
            
            #print('%s | %.2f | %.2f | %.2f | %.2f | %.2f | %s | %s ' % (names[i], alt[0], alt[-1], np.max(alt), np.min(alt), start_moon_illumination, Time(moon_rise.iso, format='iso', scale='utc', location=location) + dt, Time(moon_set.iso, format='iso', scale='utc', location=location) + dt))
            
            if (Time(moon_rise.iso, format='iso', scale='utc', location=location) + dt > Time(np.min(ttime), format='iso', scale='utc', location=location) + dt and Time(moon_rise.iso, format='iso', scale='utc', location=location) + dt < Time(np.max(ttime), format='iso', scale='utc', location=location) + dt) or (Time(moon_set.iso, format='iso', scale='utc', location=location) + dt > Time(np.min(ttime), format='iso', scale='utc', location=location) + dt and Time(moon_set.iso, format='iso', scale='utc', location=location) + dt < Time(np.max(ttime), format='iso', scale='utc', location=location) + dt) or (Time(moon_rise.iso, format='iso', scale='utc', location=location) + dt < Time(np.min(ttime), format='iso', scale='utc', location=location) + dt and Time(moon_set.iso, format='iso', scale='utc', location=location) + dt > Time(np.max(ttime), format='iso', scale='utc', location=location) + dt):
                # Moon is present during this time of observations
                moon_illumination = Moon_illumination
            else:
                moon_illumination = 0.0
            
            print('%s | %s | %s | %s | %.2f | %.2f | %.2f ' % (str(i).ljust(3), str(names[i]).ljust(20), Time(np.min(ttime), format='iso', scale='utc', location=location) + dt, Time(np.max(ttime), format='iso', scale='utc', location=location) + dt, alt[inds[0]], alt[inds[-1]], moon_illumination))
            sel_targets.append(i)
            
            
    print('======================================================================================================')

    if survey is not None:
        downl_targets = input('\n Select indices for which you want to download %s jpgs?' % (survey))

        width = 1.0
        width = float(input('\n Enter the angular width (in arcmin) of your output pictures (default: %.1f arcmin)?' % (width)) or width)

        if downl_targets=='':
            downl_targets = sel_targets
        else:
            downl_targets = list(map(int, downl_targets.split(',')))


        #print(downl_targets)
        #exit()
        #print(float(target.ra.deg))
        #exit()

        #os.mkdir('pics_%s' % (table.split('.')[0]))
        #

        if os.path.exists('pics_%s' % (table.split('.')[0])):
            shutil.rmtree('pics_%s' % (table.split('.')[0]))
        os.makedirs('pics_%s' % (table.split('.')[0]))
        os.chdir('pics_%s' % (table.split('.')[0]))
        #exit()
        for ii in downl_targets:
            #print(ii)
            #exit()
            if survey.lower()=='sdss':
                download_sdss_jpg.download_sdss(float(targets[ii].ra.deg), float(targets[ii].dec.deg), width=width, output_file='%i_%s.jpg' % (ii,names[ii]), add_bar=True, text=None)
            elif survey.lower()=='panstarrs':
                download_panstarrs_jpg.download_panstarrs(float(targets[ii].ra.deg), float(targets[ii].dec.deg), width=width, output_file='%i_%s.jpg' % (ii,names[ii]))
            elif survey.lower()=='desi':
                download_legacy_jpg_new.main(ii, names[ii], float(targets[ii].ra.deg), float(targets[ii].dec.deg), width/2.0, 0.0, None, pixscale=0.262, resolution=600, brightness_factor=4.0, contrast_factor=15.,sharpness_factor=0.01, invert=True, composite=True, output_dir='./',L_bar=30., output_file='%i_%s.png' % (ii,names[ii]))
        os.chdir('..')
 


# local_time_start='2023-05-17 00:05:00', local_time_end='2023-05-17 01:04:00'
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Plan Observations") 
    parser.add_argument("table", help="Input table with targets", type=str)
    
    parser.add_argument("start_date", help="Start date in the format yyyy-mm-dd, e.g. 2023-05-17")
    parser.add_argument("start_time", help="Start time in the format hh-mm-ss, e.g. 00:05:00")
    parser.add_argument("end_date", help="End date in the format hh-mm-ss, e.g. 01:04:00")
    parser.add_argument("end_time", help="End time")
    parser.add_argument("--d", help="Optional: delimiter", type=str, default=' ')
    parser.add_argument("--h", help="Optional: header line with the names of the columns", type=int, default=0) 
    parser.add_argument("--v", help="Optional: plot visibility", default=False, action="store_true")
    parser.add_argument("--s", help="Optional: survey from which you would like to get a pciture of this object [SDSS, PanSTARRS, DESI]", type=str, default=None)

    args = parser.parse_args()
    
    table = args.table
    start_date = args.start_date
    start_time = args.start_time
    end_date = args.end_date
    end_time = args.end_time

    delimiter = args.d
    header_line = args.h
    visibility = args.v
    survey = args.s
    
    targets, names = read_targets(table, delimiter=delimiter, header_line=header_line, units_line=None, skip_lines = [], comment='#')
    main(table, targets, names, local_time_start='%s %s' % (start_date,start_time), local_time_end='%s %s' % (end_date,end_time), plot_visibility=visibility, survey=survey)

'''
#targets,names = read_targets('apo1.lst')
targets,names = read_targets('apo1.lst', delimiter=',', header_line=0, units_line=None, skip_lines = [], comment='#')
print(targets,names)
exit()
main(targets,names)
'''



# HA = LST - RA (Angle west of the meredian (usually expressed in hours from -12 to +12)
# Targets on the meridian have LST = RA and thus HA=0
# Rising targets have HA<0
# setting targets have HA>0

# Airmass increases rapidly with zenith angles>60
# Best to observe targets at airmass < 2-3
# Astronomical twilight: 18 deg below the horizon

# https://www.esrl.noaa.gov/gmd/grad/solcalc/sunrise.html
# http://users.softlab.ntua.gr/~ipanag/fromnetmode/scripts/suntime.htm
# https://stardate.org/nightsky/moon
# https://www.mooncalc.org/

