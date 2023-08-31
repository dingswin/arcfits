#!/usr/bin/env python
"""
NOTE
----
Forked out from https://github.com/AXXE251/AGN-JET-RL that was originally created on 20230712 by Dr. Wei Zhao and Dr. Xiaofeng Li
"""
import sys, os
from astropy.io import fits
from astropy.table import Table
import numpy as np
from numpy import unravel_index
from scipy import ndimage
from scipy import special as sp
from arcfits import others as _others
from arcfits import plot as _plot

def ____AGN_jet_position_angle(fitsimage):
    """
    All function started with '____' are deprecated.

    Input parameters
    ----------------
    fitsimage : str
        VLBI fits image, normally the AGN model well made for calibration.
    """
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

def beam_info(header):
    """
    Note
    ----
    Here it is assumed that header['CDELT1'] == header['CDELT2'].
    """
    beam_maj = header['BMAJ'] / min(np.abs([header['CDELT1'], header['CDELT2']])) ## in pixel
    beam_min = header['BMIN'] / min(np.abs([header['CDELT1'], header['CDELT2']])) ## in pixel
    beam_PA = header['BPA'] ## in deg
    return beam_maj, beam_min, beam_PA
    

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
def ____resample_the_in_the_polar_coordinate(header, data=None, refpix_RA=None, refpix_Dec=None):
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
def resample_the_in_a_coordinate_resembling_the_synthesized_beam(header, data=None, refpix_RA=None, refpix_Dec=None):
    """
    Input parameters
    ----------------
    header : 
        hdu[0].header

    Output paramters
    ----------------
    pixels_RA : 2D numpy array of float (in pixel)
        RAs of the positions on the ellipse
    pixels_Dec : 2D numpy array of float (in pixel)
        Decs of the positions on the ellipse
    Rs : 1D numpy array of float (in pixel)
        long axis R values of the beam-like elliptic coordinate.
    PAs : 1D numpy array of float (in rad)
        position angle values of the elliptic coordinate.
    """
    if (refpix_RA==None) and (refpix_Dec==None):
        use_brightest_spot_as_refpix = True
        if data.all() == None:
            print('data has to be provided when use_brightest_spot_as_refpix; aborting...')
            sys.exit()
        refpix_RA, refpix_Dec = find_the_pixel_with_the_highest_flux_density(data)
    
    beam_maj, beam_min, beam_PA = beam_info(header) ## in pixel and deg
    ratio_min_maj = float(beam_min) / float(beam_maj) 
    R_min = beam_maj / 2. ## smallest R to estimate jet direction (within R0 is the beam); in pixel
    length_RA = header['NAXIS1'] ## in pixel
    length_Dec = header['NAXIS2'] ## in pixel
    R_max = min(length_RA-refpix_RA, refpix_RA, length_Dec-refpix_Dec, refpix_Dec) ## in pixel
    DelR = 10
    Rs = np.arange(R_min, R_max+DelR, DelR) ## in pixel
    PAs = position_angles = np.linspace(0, 2*np.pi, 500) ## east of north (anti-clockwise in the image)

    pixels_RA, pixels_Dec = project_beam_like_elliptic_coordinate_to_Cartesian(Rs, PAs, ratio_min_maj, beam_PA, refpix_RA, refpix_Dec)
    return pixels_RA, pixels_Dec, Rs, PAs
def project_beam_like_elliptic_coordinate_to_Cartesian(Rs, PAs, ratio_min_maj, beam_PA, refpix_RA, refpix_Dec):
    k = ratio_min_maj
    PA0 = beam_PA * np.pi / 180. ## in rad
    y2s = np.array([Rs]).T * np.cos(PAs - PA0) ## 2D numpy array of float
    x2s = -k * np.array([Rs]).T * np.sin(PAs - PA0) ## 2D numpy array of float
    
    x1s = x2s * np.cos(PA0) - y2s * np.sin(PA0) ## 2D numpy array of float
    y1s = y2s * np.cos(PA0) + x2s * np.sin(PA0) ## 2D numpy array of float
    pixels_RA = x1s + refpix_RA ## 2D numpy array of float
    pixels_Dec = y1s + refpix_Dec ## 2D numpy array of float
    return pixels_RA, pixels_Dec
def project_beam_like_elliptic_coordinate_to_Cartesian_1D(Rs, PAs, ratio_min_maj, beam_PA, refpix_RA, refpix_Dec):
    if len(Rs) != len(PAs):
        print('The length of Rs needs to match the length of PAs! Aborting...')
        sys.exit()
    k = ratio_min_maj
    PA0 = beam_PA * np.pi / 180. ## in rad
    y2s = np.array(Rs) * np.cos(PAs - PA0) ## 1D numpy array of float
    x2s = -k * np.array(Rs) * np.sin(PAs - PA0) ## 1D numpy array of float
    
    x1s = x2s * np.cos(PA0) - y2s * np.sin(PA0) ## 1D numpy array of float
    y1s = y2s * np.cos(PA0) + x2s * np.sin(PA0) ## 1D numpy array of float
    pixels_RA = x1s + refpix_RA ## 1D numpy array of float
    pixels_Dec = y1s + refpix_Dec ## 1D numpy array of float
    return pixels_RA, pixels_Dec

    
def find_the_pixel_with_the_highest_flux_density(data):
    refpix_Dec, refpix_RA = unravel_index(data.argmax(), data.shape)
    return refpix_RA, refpix_Dec
def ____get_flux_densities_in_polar_coordinate(data, header):
    """
    Output parameters
    -----------------
    fluxes : 2D array of float
        in Jy/beam. fluxes as a map of R and PA.
    """
    pixels_RA, pixels_Dec, Rs, PAs = ____resample_the_in_the_polar_coordinate(header, data)
    fluxes = ndimage.map_coordinates(data, np.stack((pixels_Dec, pixels_RA))) ## determine the flux density along the slices
    return fluxes, Rs, PAs
def get_flux_densities_in_beam_like_elliptic_coordinate(data, header):
    """
    Output parameters
    -----------------
    fluxes : 2D array of float
        in Jy/beam. fluxes as a map of R and PA.
    """
    pixels_RA, pixels_Dec, Rs, PAs = resample_the_in_a_coordinate_resembling_the_synthesized_beam(header, data)
    fluxes = ndimage.map_coordinates(data, np.stack((pixels_Dec, pixels_RA))) ## determine the flux density along the slices
    return fluxes, Rs, PAs

def ____fluxes_on_a_circle_around_the_reference_pixel(data, header, diffpix):
    pixels_RA, pixels_Dec, PAs = ____pixels_on_a_circle_around_the_reference_pixel(header, diffpix)
    fluxes = get_flux_densities_at_a_subset_of_pixels(data, pixels_RA, pixels_Dec)
    return PAs, fluxes

def obtain_jet_ridgeline(fitsimage, how_many_rms=7, how_many_sigma=4, write_out_to_file=True):
    """
    Input parameters
    ----------------
    how_many_rms : int
        The significance of detection in the fits image.
    how_many_sigma : int
        This parameter is used to estimate the position angle range, in which the jet ridgeline resides at how_many_sigma significance.
    """
    data, header = read_fits_image(fitsimage)
    flux_threshold = derive_flux_threshold(data, True, how_many_rms) ## in Jy/beam
    flux_threshold_3rms = derive_flux_threshold(data, True, 3)
    #print(flux_threshold_3rms)
    flux_density_deficit, junk1, junk2 = calculate_maximum_flux_density_deficit(data, how_many_sigma)
    fluxes, Rs, PAs = get_flux_densities_in_beam_like_elliptic_coordinate(data, header)
    
    no_extended_radio_feature = True
    fluxes_max = np.array([])
    Rs_max = np.array([])
    PAs_max = np.array([])
    PAs_max_lower = np.array([])
    PAs_max_upper = np.array([])
    for i in range(len(Rs)):
        R = Rs[i]
        fluxes_at_R = fluxes[i,:] ## actually on an ellipse with the half major axis of R
        max_flux = max(fluxes_at_R)
        #flux_75th_percentile = np.percentile(fluxes_at_R, 75) ## the largest quarter of the fluxes_at_R
        #flux_25th_percentile = np.percentile(fluxes_at_R, 25) ## the largest quarter of the fluxes_at_R
        #print(flux_75th_percentile, flux_25th_percentile)
        median_flux = np.median(fluxes_at_R)
        if max_flux > flux_threshold:
            if median_flux > flux_threshold_3rms:
                print('Still inside the jet-base blob. Skip to the next step...')
                continue
            else:
                #print(fluxes_at_R)
                index = fluxes_at_R == max_flux
                PA_max = PAs[index]

                lower_flux = max_flux - flux_density_deficit ## in Jy/beam
                PAs_target = flux_to_position_angles(fluxes_at_R, PAs, lower_flux)
                if len(PAs_target) != 2:
                    print('The number of solutions is %d instead of 2! Possibly the solutions do not come from extended jet structures (instead belong to the unresolved region). Skip to the next step...' % len(PAs_target))
                    continue
                else:
                    no_extended_radio_feature = False
                    PAs_max = np.append(PAs_max, PA_max)
                    fluxes_max = np.append(fluxes_max, max_flux)
                    Rs_max = np.append(Rs_max, R)
                    _plot.flux_to_PA_relation(fluxes_at_R, PAs, PAs_max, PAs_target, flux_threshold, R)
                    PAs_max_lower = np.append(PAs_max_lower, min(PAs_target))
                    PAs_max_upper = np.append(PAs_max_upper, max(PAs_target))

    if no_extended_radio_feature:
        print('No radio component other than the jet core is detected in the image. In other words, the radio jet is too compact for determining the jet direction. Aborting...')
        sys.exit()
    if write_out_to_file:
        table = Table([PAs_max, PAs_max_lower, PAs_max_upper, Rs_max, fluxes_max], names = ['PA_max', 'PA_max_lower', 'PA_max_upper', 'R_max', 'flux_max'])
        fitsimagename = fitsimage.split('/')[-1]
        output = fitsimagename.replace('clean.fits', 'position_angles_at_jet_ridgeline.dat')
        table.write(output, format='ascii', overwrite=True)
    return table

def flux_to_position_angles(flux_chain, PAs, flux_target):
    fluxes = np.array(flux_chain)
    flux_diffs = fluxes - flux_target
    multiply = flux_diffs[1:] * flux_diffs[:-1]
    index = multiply <= 0
    PAs_target = PAs[:-1][index]
    return PAs_target

def derive_flux_threshold(data, auto_mask=True, how_many_rms=7):
    data1 = data.flatten()
    if auto_mask: ## remove detected components in an iterative manner
        data1 = remove_detected_components_iteratively(data1, 7, 10)
    
    std_flux = np.std(data1)
    mean_flux = np.mean(data1)
    flux_threshold = mean_flux + how_many_rms * std_flux ## in Jy/beam
    return flux_threshold ## in Jy/beam
def remove_detected_components_iteratively(chain, how_many_rms=7, iterations=10):
    count = 1
    while count < iterations:
        std_flux = np.std(chain)
        mean_flux = np.mean(chain)
        mask_threshold = mean_flux + how_many_rms * std_flux ## in Jy/beam
        index = chain < mask_threshold
        chain = chain[index]
        count += 1
    return chain
def count_negative_flux_densities(data, how_many_rms=3):
    data1 = data.flatten()
    data1 = remove_detected_components_iteratively(data1, 7, 10)
    number_noises = len(data1)
    
    std_flux = np.std(data1)
    mean_flux = np.mean(data1)
    threshold = mean_flux - how_many_rms * std_flux ## in Jy/beam
    index = data1 < threshold
    number_negative = len(data1[index])
    prob_negative = number_negative / number_noises
    return number_negative, number_noises, prob_negative, std_flux
def calculate_maximum_flux_density_deficit(data, how_many_sigma=4):
    prob_target = 1 - sp.erf(how_many_sigma / 2**0.5)
    junk1, number_noises, junk2, rms = count_negative_flux_densities(data, 1)
    min_prob = 1. / number_noises
    if prob_target < min_prob:
        print('The significance target is not reachable. Use a lower how_many_sigma. Aborting...')
        sys.exit()
    how_many_rms_low = 0
    how_many_rms_high = 10
    how_many_rms = 0
    prob = 1
    iterations = 0
    iteration_limit = 100
    while abs(prob - prob_target) > 0.1 * prob_target:
        how_many_rms = (how_many_rms_low + how_many_rms_high) / 2.
        junk1, junk2, prob, junk3 = count_negative_flux_densities(data, how_many_rms)
        if prob > prob_target:
            how_many_rms_low = how_many_rms
        else:
            how_many_rms_high = how_many_rms
        iterations += 1
        if iterations > iteration_limit:
            print('iteration limit reached. Use a lower how_many_sigma. Aborting...')
            sys.exit()

    how_many_rms_down = how_many_rms
    max_flux_density_deficit = how_many_rms_down * rms
    return max_flux_density_deficit, how_many_rms_down, iterations
