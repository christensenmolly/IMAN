import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits as pyfits
from math import *
import sys
import os
import subprocess
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
from shapely.geometry import LineString
import scipy.ndimage as ndimage


LOCAL_DIR = "/fit_features/tidal_features"
IMAN_DIR = os.path.dirname(__file__).rpartition(LOCAL_DIR)[0]

sys.path.append(os.path.join(IMAN_DIR, 'misc_funcs'))
sys.path.append(os.path.join(IMAN_DIR, 'imp/masking'))
sys.path.append(os.path.join(IMAN_DIR, 'cosmology'))
import image_interpolation_astropy
import convert_reg_to_mask
import cosmo_calc_NED_func

def get_length_of_polygon(poly):
    # get minimum bounding box around polygon
    box = poly.minimum_rotated_rectangle

    # get coordinates of polygon vertices
    x, y = box.exterior.coords.xy
    
    # get length of bounding box edges
    edge_length = (Point(x[0], y[0]).distance(Point(x[1], y[1])), Point(x[1], y[1]).distance(Point(x[2], y[2])))

    # get length of polygon as the longest edge of the bounding box
    length = max(edge_length)

    # get width of polygon as the shortest edge of the bounding box
    width = min(edge_length)
    x,y = box.exterior.xy
    #plt.plot(x,y)
    #plt.show()
    return length


def sort_to_form_line(unsorted_list):
    """
    Given a list of neighboring points which forms a line, but in random order, 
    sort them to the correct order.
    IMPORTANT: Each point must be a neighbor (8-point sense) 
    to a least one other point!
    """
    sorted_list = [unsorted_list.pop(0)]

    while len(unsorted_list) > 0:
        i = 0
        while i < len(unsorted_list):
            if are_neighbours(sorted_list[0], unsorted_list[i]):
                #neighbours at front of list
                sorted_list.insert(0, unsorted_list.pop(i))
            elif are_neighbours(sorted_list[-1], unsorted_list[i]):
                #neighbours at rear of list
                sorted_list.append(unsorted_list.pop(i))
            else:
                i = i+1

    return sorted_list

def are_neighbours(pt1, pt2):
    """
    Check if pt1 and pt2 are neighbours, in the 8-point sense
    pt1 and pt2 has integer coordinates
    """
    return (np.abs(pt1[0]-pt2[0]) < 2) and (np.abs(pt1[1]-pt2[1]) < 2)

def order_points(x,y):
    from sklearn.neighbors import NearestNeighbors
    import networkx as nx

    points = np.c_[x, y]

    clf = NearestNeighbors(2).fit(points)
    G = clf.kneighbors_graph()

    T = nx.from_scipy_sparse_matrix(G)

    order = list(nx.dfs_preorder_nodes(T, 0))

    xx = x[order]
    yy = y[order]
    '''
    paths = [list(nx.dfs_preorder_nodes(T, i)) for i in range(len(points))]

    mindist = np.inf
    minidx = 0

    for i in range(len(points)):
        p = paths[i]           # order of nodes
        ordered = points[p]    # ordered nodes
        # find cost of that order by the sum of euclidean distances between points (i) and (i+1)
        cost = (((ordered[:-1] - ordered[1:])**2).sum(1)).sum()
        if cost < mindist:
            mindist = cost
            minidx = i

    opt_order = paths[minidx]

    xx = x[opt_order]
    yy = y[opt_order]
    '''
    return xx,yy


def distance(P1, P2):
    """
    This function computes the distance between 2 points defined by
    P1 = (x1,y1) and P2 = (x2,y2) 
    """

    return ((P1[0] - P2[0])**2 + (P1[1] - P2[1])**2) ** 0.5


def optimized_path(coords, start=None):
    """
    This function finds the nearest point to a point
    coords should be a list in this format coords = [ [x1, y1], [x2, y2] , ...] 

    """
    if start is None:
        start = coords[0]
    pass_by = coords
    path = [start]
    pass_by.remove(start)
    while pass_by:
        nearest = min(pass_by, key=lambda x: distance(path[-1], x))
        path.append(nearest)
        pass_by.remove(nearest)
    return path

def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth


def func(x, a, b, c):
    return a * np.exp(-b * x) + c

def read_centerline(region_file):
    f = open(region_file, 'r')
    lines = f.readlines()
    coords = []
    x = []
    y = []
    for line in lines:
        if 'circle(' in line:
            params = line.split("(")[-1].split(",")[0]
            xc = float(line.split("(")[-1].split(",")[0])
            yc = float(line.split("(")[-1].split(",")[1])
            coords.append([int(xc),int(yc)])
            x.append(int(xc))
            y.append(int(yc))
    x = np.array(x)
    y = np.array(y)
    plt.plot(x,y)        
    # Sort data according to x.
    #temp_data = sorted(zip(x, y))

    # Unpack sorted data.
    #x_sorted, y_sorted = zip(*temp_data)
    
    #xx,yy = order_points(x,y)
    
    
    
    #plt.plot(xx,yy)
    #plt.show()
    #my_coords = sort_to_form_line(coords)
    
    f = open('corrected_centerline.reg', 'w')
    f.write('image\n')
    for k in range(len(xx)):
        
        f.write('circle(%.1f,%.1f,1) # color=red\n' % (float(xx[k]),float(yy[k])))
    f.close()
    exit()
    

    #my_coords = optimized_path(coords,coords[0])
    
    #ZZ = np.polyfit(x, y, 11)
    #ff = np.poly1d(ZZ)
    #plt.plot(x,ff(x),'o',color='green')
    
    #from scipy.interpolate import splprep, splev
    #x = np.array([23, 24, 24, 25, 25])
    #y = np.array([13, 12, 13, 12, 13])
    #pts = np.vstack((xx,yy))
    #tck, u = splprep([x,y], s=0)
    #new_points = splev(u, tck)
    #plt.plot(new_points[0], new_points[1], 'r-')

    
    coords = []
    for k in range(len(xx)):
        coords.append((xx[k],yy[k]))

    Line = LineString(coords)
    x,y = Line.xy
    plt.plot(x,y)
    plt.show()
    return coords,Line

def read_line(region_file):
    f = open(region_file, 'r')
    lines = f.readlines()
    coords = []
    for line in lines:
        if 'line(' in line:
            x1 = float(line.split("(")[-1].split(",")[0])
            y1 = float(line.split("(")[-1].split(",")[1])
            x2 = float(line.split("(")[-1].split(",")[2])
            y2 = float(line.split("(")[-1].split(",")[3].split(")")[0])
    return (y2-y1)/(x2-x1)

def find_polygon_center(region_file):
    f = open(region_file, 'r')
    lines = f.readlines()
    new_lines = []
    point_numbers = []
    for line in lines:
        if 'polygon(' in line:
            point_numbers.append(len(line))
            new_lines.append(line)
    f.close()

    coords = new_lines[0].replace('(',',').replace(')',',').split(',')[1:-1]
    pol = []
    for kk in range(0,len(coords)-1,2):
        pol.append((int(float(coords[kk])),int(float(coords[kk+1]))))
    polygon = Polygon(pol)
    x,y = polygon.exterior.xy
    plt.plot(x,y)

    [(x,y)] = polygon.centroid.coords
    length = get_length_of_polygon(polygon)


    #plt.show()
    #exit()

    return x,y, length, polygon.area    



def fit_polygon(data, polygon_segm, average_sb, std_sb):
    ny,nx = np.shape(data)
    x = []
    y = []
    for k in range(ny):
        for i in range(nx):
            if polygon_segm[k,i]!=0 and data[k,i]>average_sb-0.5*std_sb and data[k,i]<average_sb+0.5*std_sb:
                x.append(i)
                y.append(k)
    x = np.array(x)
    y = np.array(y)
    
    p = np.polyfit(x,y,2) 
    f = np.poly1d(p)
    f_deriv = f.deriv(1)
    
    xx = np.arange(min(x),max(x),0.1)
    yy = f(xx)
    
    #plt.plot(x, y, 'bo', label="Data")
    #plt.plot(xx,f(xx), 'b-', color='red', label="Polyfit")
    #plt.show()
    
    coords = []
    for k in range(len(xx)):
        coords.append((xx[k],yy[k]))

    Line = LineString(coords)
    
    return Line.length,coords,f_deriv
    
    
def f(x, A,B,C):
    return A*x**2 + B*x + C
    
def fit_polygon1(data, polygon_segm):
    from scipy.optimize import curve_fit, fmin
    

    ny,nx = np.shape(data)
    x = []
    y = []
    noise_sigma = []
    for k in range(ny):
        for i in range(nx):
            if polygon_segm[k,i]!=0:
                x.append(i)
                y.append(k)
                noise_sigma.append(data[k,i])
    x = np.array(x)
    y = np.array(y)
    noise_sigma = np.array(noise_sigma)
    
    popt, pcov = curve_fit(f, x, y)#, sigma=1./noise_sigma, absolute_sigma=True)

    xx = np.arange(min(x),max(x),0.1)
    
    plt.plot(x, y, 'bo', label="Data")
    plt.plot(xx,f(xx, popt[0],popt[1],popt[2]), 'b-', color='red', label="Polyfit")
    plt.show()    


def find_flux(data, segm_data):
    ny,nx = np.shape(data)
    flux = np.nansum(data[segm_data!=0])
    return flux

def average_sb(data, segm_data):    
    ny,nx = np.shape(data)
    av_intensity = np.nanmean(data[segm_data!=0])
    std_intensity = np.nanstd(data[segm_data!=0])
    return av_intensity,std_intensity


def find_nearest_point_on_centerline_to_center(xc, yc, coords):
    D2 = []
    for k in range(len(coords)):
        D2.append((xc-coords[k][0])**2 + (yc-coords[k][1])**2)
    return coords[D2.index(min(D2))][0],coords[D2.index(min(D2))][1],D2.index(min(D2))


def main(input_image, galaxy_region_file, lsb_structure_region_file, lsb_structure_centerline_file, polar_line, main_line, mask_image, m0, pix2sec, z, sigma_smooth=3.):
    DA_Mpc,kpc_DA,DL_Mpc = cosmo_calc_NED_func.main(z, H0=70., WM=0.3, WV=0.7)
    Scale = kpc_DA
    print("Scale (kpc/arcsec): %.3f" % (Scale))
    
    # Convert galaxy polygon to segm
    convert_reg_to_mask.mask(input_image, galaxy_region_file, output_image=None, output_mask='tmp_galaxy.fits', mask_value=1, mask_DN=None, verbosity=False)
    
    # Convert LSB structure polygon to segm
    convert_reg_to_mask.mask(input_image, lsb_structure_region_file, output_image=None, output_mask='tmp_lsb_structure.fits', mask_value=2, mask_DN=None, verbosity=False)

    # Do interpolation of the image within the mask
    image_interpolation_astropy.astropy_smoothing(input_image, mask_image, output_image='interpolated.fits', sigma_smooth=5., sampling_factor=None, sigma_back=None)
    
    hdulist = pyfits.open('interpolated.fits')
    img = hdulist[0].data
    ny,nx = np.shape(img)
    
    # Smooth the image: this is done mainly for smoothing noisy LSB structure
    smoothed_img = ndimage.gaussian_filter(img, sigma=sigma_smooth, order=0) # Now we work with this smoothed data only!
    outHDU = pyfits.PrimaryHDU(smoothed_img)
    outHDU.writeto('smoothed_img.fits', overwrite=True)  

    
    
    
    

    hdulist1 = pyfits.open('tmp_galaxy.fits')
    galaxy_segm = hdulist1[0].data

    hdulist2 = pyfits.open('tmp_lsb_structure.fits')
    lsb_structure_segm = hdulist2[0].data

    # 1. Find average SB inside the LSB structure:
    av_intensity_lsb_structure,std_intensity_lsb_structure = average_sb(smoothed_img, lsb_structure_segm)
    print('Average surface brightness of the stream (mag/arcsec2): %.2f+/-%.2f' % (m0-2.5*log10(av_intensity_lsb_structure/(pix2sec**2)),2.5*std_intensity_lsb_structure/(av_intensity_lsb_structure*2.303)))
    
    
    # 2. Find the length of the SLB structure:
    #read_centerline(lsb_structure_centerline_file)
    #exit()
    
    
    LSB_length,coords_centerline,f_deriv = fit_polygon(smoothed_img, lsb_structure_segm,av_intensity_lsb_structure,std_intensity_lsb_structure)
    print('Length of the stream (pix): %.1f' % (LSB_length))
    print('Length of the stream (arcsec): %.1f' % (LSB_length*pix2sec))
    print('Length of the stream (kpc): %.2f' % (LSB_length*pix2sec*Scale))
    #exit()
    



    galaxy_no_lsb_segm = np.copy(galaxy_segm)

    for k in range(ny):
        for i in range(nx):
            if galaxy_segm[k,i]==1:
                if lsb_structure_segm[k,i]==1:
                    galaxy_no_lsb_segm[k,i]=0.

    outHDU = pyfits.PrimaryHDU(galaxy_no_lsb_segm)
    outHDU.writeto('tmp_galaxy_minus_lsb_structure.fits', overwrite=True)  

    





    # Determine parameters of the galaxy 
    galaxy_total_flux = find_flux(smoothed_img, galaxy_segm)

    lsb_structure_flux = find_flux(smoothed_img, lsb_structure_segm)
    f = lsb_structure_flux/galaxy_total_flux

    print('Total magnitude of the galaxy: %.2f' % (m0-2.5*log10(galaxy_total_flux)))
    print('Total magnitude of the stream: %.2f' % (m0-2.5*log10(lsb_structure_flux)))
    print('Fraction of the stream: %.3f' % (f))



    
    xc_lsb_structure,yc_lsb_structure, length, area = find_polygon_center(lsb_structure_region_file)
    print('Center of polygon: %.1f, %.1f' % (xc_lsb_structure,yc_lsb_structure))
    print('Length of the stream-box (pix): %.1f' % (length))
    print('Length of the stream-box (arcsec): %.1f' % (length*pix2sec))
    print('Length of the stream-box (kpc): %.2f' % (length*pix2sec*Scale))
    print('Area of the polygon (pix2): %.2f' % (area))
    print('Area of the polygon (arcsec2): %.2f' % (area*pix2sec*pix2sec))
    print('Area of the polygon (kpc2): %.2f' % (area*pix2sec*pix2sec*Scale*Scale))
    #coords_centerline, shapely_centerline = read_centerline(lsb_structure_centerline_file)

    #LSB_length = shapely_centerline.length
    
    #print('Length of the stream (pix): %.1f' % (LSB_length))
    #exit()

    # Find nearest point on centerline to the center of the polygon
    xc_line,yc_line,ind = find_nearest_point_on_centerline_to_center(xc_lsb_structure, yc_lsb_structure, coords_centerline)
    print('Nearest point on centerline to center: %.1f, %.1f' % (xc_line,yc_line))
    
    several_points_x = []
    several_points_y = []
    #for k in range(ind-2,ind+3):
    #    several_points_x.append(coords_centerline[k][0])
    #    several_points_y.append(coords_centerline[k][1])
    
    #slope = np.polyfit(several_points_x, several_points_y, 1)[0]
    
    slope = f_deriv(xc_line)
    
    slope_polar = read_line(polar_line)
    slope_main = read_line(main_line)

    tan_psi_with_polar = (slope-slope_main)/(1.+slope_polar*slope)
    psi_with_polar = np.degrees(atan(tan_psi_with_polar))
    tan_psi_with_main  = (slope-slope_main)/(1.+slope_main*slope)
    psi_with_main = np.degrees(atan(tan_psi_with_main))
    tan_pa_polar_main = (slope_main-slope_polar)/(1.+slope_polar*slope_main)
    pa_polar_main = np.degrees(atan(tan_pa_polar_main))

    print('Angle between polar ring and main body (deg): +/-%.1f' % (pa_polar_main))
    print('Angle between polar ring and stream (deg): +/-%.1f' % (psi_with_polar))
    print('Angle between main body and stream (deg): +/-%.1f' % (psi_with_main))
    
    print('Done!')
'''
input_image = '73_stacked.fits'
galaxy_region_file = 'total_galaxy.reg'
lsb_structure_region_file = 'stream.reg'
lsb_structure_centerline_file = 'centerline_stream.reg'
polar_line = 'polar_line.reg'
main_line = 'main_line.reg'
mask_image = 'stars_on_stream.fits'
m0 = 28.17
pix2sec = 0.396
z = 0.05909

main(input_image, galaxy_region_file, lsb_structure_region_file, lsb_structure_centerline_file, polar_line, main_line, mask_image, m0, pix2sec, z)
'''