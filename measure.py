#!/usr/bin/env python
"""
NOTE
----
Forked out from https://github.com/AXXE251/AGN-JET-RL that was originally created on 20230712 by Dr. Wei Zhao and Dr. Xiaofeng Li
"""
from astropy.io import fits
import numpy as np
from numpy import unravel_index
from scipy import ndimage
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, SymLogNorm
from matplotlib.patches import Ellipse
import matplotlib
matplotlib.rcParams.update({'font.size': 14})

def word2pix(h,xy):
	x,y = xy
	x = h['crpix1'] -1 + x/(h['cdelt1']*3.6E6)
	y = h['crpix2'] -1 + y/(h['cdelt2']*3.6E6)
	return x,y
def pix2word(h,xy):
	x,y = xy
	x = h['cdelt1']*3.6E6*(x-h['crpix1']+1)
	y = h['cdelt2']*3.6E6*(y-h['crpix2']+1)
	return x,y
def W2w(h,w=()):
	if len(w)==4:
		x0, x1, y0, y1 = w
	else:
		x0, y0 = 0, 0
		x1, y1 = h['naxis1']-1, h['naxis2']-1
	x0, y0 = pix2word(h,(x0,y0))
	x1, y1 = pix2word(h,(x1,y1))
	w = x0, x1, y0, y1
	return w

def w2W(h,w=()):
	if len(w)==4:
		x0, x1, y0, y1 = w
		x0, y0 = word2pix(h,(x0,y0))
		x1, y1 = word2pix(h,(x1,y1))
	else:
		x0, y0 = 0, 0
		x1, y1 = h['naxis1']-1, h['naxis2']-1
	w = x0, x1, y0, y1
	return w


def max_d(h,img,a,o,d,a1,a2):
    t1 = a1*(np.pi/180.0)
    t2 = a2*(np.pi/180.0)
    t=np.linspace(t1,t2,2000)
    ang=t*(180.0/np.pi)
    x = o[0] - d* np.cos(a*np.pi/180.0+t)/np.cos(t)                          #slicing with straight lines
    y = o[1] + d* np.sin(a*np.pi/180.0+t)/np.cos(t)                          #slicing with straight lines
    #x = o[0] - d* np.cos(a*np.pi/180.0+t)                                   #slicing with circular arcs
    #y = o[1] + d* np.sin(a*np.pi/180.0+t)                                   #slicing with circular arcs
    X, Y = word2pix(h,(x,y))
    flux = ndimage.map_coordinates(np.transpose(img),np.vstack((X,Y))) #determine the flux density along the slices
    flux_max=flux[np.argmax(flux)]                                           #determine the max flux density along the slices
    ang_max=ang[np.argmax(flux)]                                             #determine the position angle of the max 
    xmax=x[np.argmax(flux)]                                                  #determine the position of the max  in R.A.
    ymax=y[np.argmax(flux)]                                                  #determine the position of the max  in Dec
    doc=open("RL_data.txt",'a')                
    print>>doc,xmax,ymax,flux_max,ang_max
    doc.close() 

def AGN_jet_position_angle(fitsimage):
    """
    Input parameters
    ----------------
    fitsimage : str
        VLBI fits image, normally the AGN model well made for calibration.
    """
    hdu = fits.open(fitsimage)
    h = hdu[0]
    h.header
    return h
    sys.exit()

    a=-90        # assumed jet orientation at the starting point
    o=(0,0)      # starting point for slicing
    jl=1.6       # slicing range in radial direction in mas 
    a1=120       # slicing range in azimuthal dirction in degrees upper limit
    a2=-120      # slicing range in azimuthal dirction in degrees lower limit
    sl=0.01      # step length in radial direction  in mas 
    n=int(jl/0.01) # number of steps in radial direction


    for i in range(n):
        max_d(h,hdu[0].data[0,0,:,:],a,o,sl*i,a1,a2)
    f=np.loadtxt('RL_data.txt')


    fig = plt.figure()
    fig.set_size_inches((5,6))
    rms = 0.0005
    levs = 3* rms * np.array([1,2,4,8,16,32,64,128,256,512,1024,2048,4096,8192])
    w = (2.5,-2.5,-4,2)
    win = W2w(h)
    plt.plot(f[0:n-1,0],f[0:n-1,1],'r')
    plt.contour(hdu[0].data[0,0,:,:],levs,extent=win,linewidths=0.5,colors='k')
    plt.tick_params('both',direction='in',right=True,top=True,which='both')
    plt.minorticks_on()
    plt.xlim(w[0],w[1])
    plt.ylim(w[2],w[3])
    plt.xlabel('Relative RA (mas)')
    plt.ylabel('Relative Dec (mas)')
    fig.tight_layout()
    plt.savefig('hello_jet.eps')

def read_fits_image(fitsimage):
    """
    Input parameter
    ---------------
    fitsimage : str
        normally AGN model well made for calibration
    """
    hdus = fits.open(fitsimage) ## HDU --> Header Data Unit
    data = hdus[0].data
    header = hdus[0].header
    hdus.close()
    ndim = data.ndim
    while ndim > 2: ## normally, ndim==4 for VLBI fits image
        data = data[0]
        ndim -= 1
    if ndim != 2:
        print('not a 2D image; aborting')
        sys.exit()
    return data, header

def raster_fits_image(fitsimage):
    data, header = read_fits_image(fitsimage)
    vmin = 3*np.mean(data)
    plt.imshow(data, cmap='cividis', norm=LogNorm(vmin=vmin))
    plt.colorbar()
    plt.show()
    plt.clf()

def ____pixels_on_a_circle_around_the_reference_pixel(header, diffpix):
    """
    Input parameters
    ----------------
    header : 
        hdu[0].header
    diffpix : int (in pixel)
        offset from the reference pixel

    Output paramters
    ----------------
    pixels_RA : numpy array of float (in pixel)
        RAs of the positions on the circle
    pixels_Dec : numpy array of float (in pixel)
        Decs of the positions on the circle
    """
    refpix_RA  = header['crpix1'] ## reference pixel in RA direction
    refpix_Dec = header['crpix2'] ## reference pixel in Dec direction
    PAs = position_angles = np.linspace(0, 2*np.pi, 1000) ## east of north
    pixels_RA  = pixels_RA_on_circle = refpix_RA + diffpix * np.sin(PAs) ## numpy array of float
    pixels_Dec = pixels_Dec_on_circle = refpix_Dec + diffpix * np.cos(PAs) ## numpy array of float
    return pixels_RA, pixels_Dec, PAs
def resample_the_in_the_polar_coordinate(header, data=None, refpix_RA=None, refpix_Dec=None):
    """
    Input parameters
    ----------------
    header : 
        hdu[0].header

    Output paramters
    ----------------
    pixels_RA : 2D numpy array of float (in pixel)
        RAs of the positions on the circle
    pixels_Dec : 2D numpy array of float (in pixel)
        Decs of the positions on the circle
    Rs : 1D numpy array of float (in pixel)
        R values of the polar coordinate.
    PAs : 1D numpy array of float (in rad)
        position angle values of the polar coordinate.
    """
    if (refpix_RA==None) and (refpix_Dec==None):
        use_brightest_spot_as_refpix = True
        if data.all() == None:
            print('data has to be provided when use_brightest_spot_as_refpix; aborting...')
            sys.exit()
        refpix_RA, refpix_Dec = find_the_pixel_with_the_highest_flux_density(data)
    #refpix_RA  = header['CRPIX1'] ## reference pixel in RA direction
    #refpix_Dec = header['CRPIX2'] ## reference pixel in Dec direction
    #beamPA = header['BPA'] * np.pi / 180. ## in rad
    #beamRA = np.abs(np.sin(beamPA)) * header['BMAJ'] ## in deg
    #beamRA /= np.abs(header['CDELT1']) ## in pixel 
    #beamDec = np.abs(np.cos(beamPA)) * header['BMAJ'] ## in deg
    #beamDec /= np.abs(header['CDELT2']) ## in pixel 
    beam_major = header['BMAJ'] / min(np.abs([header['CDELT1'], header['CDELT2']])) ## in pixel
    
    #DelR = np.mean([beamRA, beamDec]) / 2. ## increment in R (of the polar coordinate); in pixel
    #DelR = max(beamRA, beamDec) / 2. ## increment in R (of the polar coordinate); in pixel
    DelR = beam_major / 2. ## increment in R (of the polar coordinate); in pixel
    DelR *= 2. ## angular broadening factor (highly uncertain)
    length_RA = header['NAXIS1'] ## in pixel
    length_Dec = header['NAXIS2'] ## in pixel
    R_max = min(length_RA-refpix_RA, refpix_RA, length_Dec-refpix_Dec, refpix_Dec) ## in pixel
    Rs = np.arange(0+DelR, R_max+DelR, DelR) ## in pixel
    PAs = position_angles = np.linspace(0, 2*np.pi, 500) ## east of north (anti-clockwise in the image)

    pixels_RA  = refpix_RA - np.array([Rs]).T * np.sin(PAs) ## 2D numpy array of float
    pixels_Dec = refpix_Dec + np.array([Rs]).T * np.cos(PAs) ## 2D numpy array of float
    return pixels_RA, pixels_Dec, Rs, PAs
def find_the_pixel_with_the_highest_flux_density(data):
    refpix_Dec, refpix_RA = unravel_index(data.argmax(), data.shape)
    return refpix_RA, refpix_Dec
def get_flux_densities_in_polar_coordinate(data, header):
    """
    Output parameters
    -----------------
    fluxes : 2D array of float
        in Jy/beam. fluxes as a map of R and PA.
    """
    pixels_RA, pixels_Dec, Rs, PAs = resample_the_in_the_polar_coordinate(header, data)
    fluxes = ndimage.map_coordinates(data, np.stack((pixels_Dec, pixels_RA))) ## determine the flux density along the slices
    return fluxes, Rs, PAs

def ____fluxes_on_a_circle_around_the_reference_pixel(data, header, diffpix):
    pixels_RA, pixels_Dec, PAs = ____pixels_on_a_circle_around_the_reference_pixel(header, diffpix)
    fluxes = get_flux_densities_at_a_subset_of_pixels(data, pixels_RA, pixels_Dec)
    return PAs, fluxes

def obtain_jet_ridgeline(fitsimage, how_many_sigma=7):
    data, header = read_fits_image(fitsimage)
    flux_threshold = derive_flux_threshold(data, how_many_sigma) ## in Jy/beam
    fluxes, Rs, PAs = get_flux_densities_in_polar_coordinate(data, header)
    
    fluxes_max = np.array([])
    Rs_max = np.array([])
    PAs_max = np.array([])
    for i in range(len(Rs)):
        R = Rs[i]
        fluxes_at_R = fluxes[i,:]
        max_flux = max(fluxes_at_R)
        if max_flux > flux_threshold:
            #print(fluxes_at_R)
            index = fluxes_at_R == max_flux
            PA_max = PAs[index]
            PAs_max = np.append(PAs_max, PA_max)
            fluxes_max = np.append(fluxes_max, max_flux)
            Rs_max = np.append(Rs_max, R)
    return PAs_max, Rs_max, fluxes_max


def derive_flux_threshold(data, auto_mask=True, how_many_sigma=7):
    data1 = data.flatten()
    if auto_mask: ## remove detected components in an iterative manner
        iterations = 10
        count = 1
        while count < iterations:
            std_flux = np.std(data1)
            mean_flux = np.mean(data1)
            mask_threshold = mean_flux + 7 * std_flux ## in Jy/beam
            #print('mask_threshold = %f Jy/beam' % mask_threshold)
            index = data1 < mask_threshold
            data1 = data1[index]
            count += 1
    
    std_flux = np.std(data1)
    mean_flux = np.mean(data1)
    flux_threshold = mean_flux + how_many_sigma * std_flux ## in Jy/beam
    return flux_threshold ## in Jy/beam

def plot_jet_ridgeline(fitsimage, how_many_sigma=7):
    data, header = read_fits_image(fitsimage)
    refpix_RA, refpix_Dec = find_the_pixel_with_the_highest_flux_density(data)
    #refpix_RA  = header['CRPIX1'] ## reference pixel in RA direction
    #refpix_Dec = header['CRPIX2'] ## reference pixel in Dec direction
    PAs_max, Rs_max, fluxes_max = obtain_jet_ridgeline(fitsimage, how_many_sigma)
    pixRAs_max  = refpix_RA - Rs_max * np.sin(PAs_max)
    pixDecs_max = refpix_Dec + Rs_max * np.cos(PAs_max)
    pixRAs_max = np.concatenate(([refpix_RA], pixRAs_max))
    pixDecs_max = np.concatenate(([refpix_Dec], pixDecs_max))
    
    fig = plt.figure()
    plt.imshow(data, cmap='rainbow', norm=SymLogNorm(1e-2, base=10), origin='lower')
    #plt.imshow(data, cmap='cividis', norm=LogNorm(), origin='lower')
    
    plt.plot(pixRAs_max[:-1], pixDecs_max[:-1], color='white', linewidth=1)
    plt.arrow(pixRAs_max[-2], pixDecs_max[-2], pixRAs_max[-1]-pixRAs_max[-2], pixDecs_max[-1]-pixDecs_max[-2], color='white', linewidth=0.3, width=1.5, head_width=6, head_length=12, fill=True, length_includes_head=True)
    
    cbar = plt.colorbar()
    cbar.set_label(r'flux density ($\mathrm{Jy~{beam}^{-1}}$)', rotation=-90, labelpad=15)
    #plt.plot(pixRAs_max, pixDecs_max, color='white')
    fig.tight_layout()
    fitsimagename = fitsimage.split('/')[-1]
    outputpdf = fitsimagename.replace('clean.fits', 'jet_ridgeline.pdf')
    plt.savefig(outputpdf)
    #plt.show()
    plt.clf()
