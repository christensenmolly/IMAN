from astroquery.vizier import Vizier 
import astropy.units as u 
import astropy.coordinates as coord

def gaia_query(ra_deg, dec_deg, rad_deg, maxmag=20, 
               maxsources=10000): 
    """
    Query Gaia DR1 @ VizieR using astroquery.vizier
    :param ra_deg: RA in degrees
    :param dec_deg: Declination in degrees
    :param rad_deg: field radius in degrees
    :param maxmag: upper limit G magnitude (optional)
    :param maxsources: maximum number of sources
    :return: astropy.table object
    """
    vquery = Vizier(columns=['Source', 'RA_ICRS', 'DE_ICRS', 
                             'phot_g_mean_mag'], 
                    column_filters={"phot_g_mean_mag": 
                                    ("<%f" % maxmag)}, 
                    row_limit = maxsources) 
 
    field = coord.SkyCoord(ra=ra_deg, dec=dec_deg, 
                           unit=(u.deg, u.deg), 
                           frame='icrs')
    return vquery.query_region(field, 
                               width=("%fd" % rad_deg), 
                               catalog="I/345/gaia2")[0] 

# Example query
print(gaia_query(12.345, 67.89, 0.1))