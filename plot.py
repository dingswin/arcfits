#!/usr/bin/env python
import sys, os
from astropy.table import Table
import numpy as np
from arcfits import others as _others
from arcfits import measure as _measure
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from matplotlib.colors import LogNorm, SymLogNorm
#matplotlib.rcParams.update({'font.size': 14})

def flux_to_PA_relation(fluxes_at_R, PAs, PAs_max, PAs_max_limit, flux_threshold, R):
    """
    make diagnostic plots for later check
    """
    print('making a diagnostic plot for later check...')
    fig = plt.figure()
    #PAs *= 180 / np.pi ## in deg
    #PAs_max_limit *= 180. / np.pi ## in deg
    plt.plot(PAs*180/np.pi, fluxes_at_R)
    for PA in PAs_max_limit*180./np.pi:
        plt.axvline(x=PA, ls='--') ## mark the uncertainty range of PA_max
    plt.axhline(y=flux_threshold, ls = '-.') ## mark adopted flux threshold
    plt.xlabel(r'jet position angle (deg)')
    plt.ylabel(r'flux density ($\mathrm{Jy~{beam}^{-1}}$)')
    fig.tight_layout()

    #srcname = fitsimage.split('/')[-1].split('.')[0]
    outputdir = 'diagnostic_plots'
    if not os.path.exists(outputdir):
        os.mkdir(outputdir)
    output = '%s/flux_densities_at_%dpixel_from_the_jet_base.pdf' % (outputdir, R) 
    plt.savefig(output)
    plt.clf()
    plt.close()

def ____prepare_locations_in_the_image(t_RL, refpix_RA, refpix_Dec):
    PAs_max = t_RL['PA_max']
    Rs_max = t_RL['R_max']
    PAs_max_lower = t_RL['PA_max_lower']
    PAs_max_upper = t_RL['PA_max_upper']
    pixRAs_max  = refpix_RA - Rs_max * np.sin(PAs_max)
    pixRAs_max_lower  = refpix_RA - Rs_max * np.sin(PAs_max_lower)
    pixRAs_max_upper  = refpix_RA - Rs_max * np.sin(PAs_max_upper)
    pixDecs_max = refpix_Dec + Rs_max * np.cos(PAs_max)
    pixDecs_max_lower = refpix_Dec + Rs_max * np.cos(PAs_max_lower)
    pixDecs_max_upper = refpix_Dec + Rs_max * np.cos(PAs_max_upper)
    pixRAs_max = np.concatenate(([refpix_RA], pixRAs_max))
    pixRAs_max_lower = np.concatenate(([refpix_RA], pixRAs_max_lower))
    pixRAs_max_upper = np.concatenate(([refpix_RA], pixRAs_max_upper))
    pixDecs_max = np.concatenate(([refpix_Dec], pixDecs_max))
    pixDecs_max_lower = np.concatenate(([refpix_Dec], pixDecs_max_lower))
    pixDecs_max_upper = np.concatenate(([refpix_Dec], pixDecs_max_upper))
    return pixRAs_max, pixDecs_max, pixRAs_max_lower, pixDecs_max_lower, pixRAs_max_upper, pixDecs_max_upper
def prepare_locations_in_the_image(t_RL, header, refpix_RA, refpix_Dec):
    beam_maj, beam_min, beam_PA = _measure.beam_info(header) ## in pixel and deg
    ratio_min_maj = float(beam_min) / float(beam_maj) 
    PAs_max = t_RL['PA_max']
    Rs_max = t_RL['R_max']
    PAs_max_lower = t_RL['PA_max_lower']
    PAs_max_upper = t_RL['PA_max_upper']
    #pixRAs_max  = refpix_RA - Rs_max * np.sin(PAs_max)
    #pixRAs_max_lower  = refpix_RA - Rs_max * np.sin(PAs_max_lower)
    #pixRAs_max_upper  = refpix_RA - Rs_max * np.sin(PAs_max_upper)
    #pixDecs_max = refpix_Dec + Rs_max * np.cos(PAs_max)
    #pixDecs_max_lower = refpix_Dec + Rs_max * np.cos(PAs_max_lower)
    #pixDecs_max_upper = refpix_Dec + Rs_max * np.cos(PAs_max_upper)
    pixRAs_max, pixDecs_max = _measure.project_beam_like_elliptic_coordinate_to_Cartesian_1D(Rs_max, PAs_max, ratio_min_maj, beam_PA, refpix_RA, refpix_Dec)
    pixRAs_max_lower, pixDecs_max_lower = _measure.project_beam_like_elliptic_coordinate_to_Cartesian_1D(Rs_max, PAs_max_lower, ratio_min_maj, beam_PA, refpix_RA, refpix_Dec)
    pixRAs_max_upper, pixDecs_max_upper = _measure.project_beam_like_elliptic_coordinate_to_Cartesian_1D(Rs_max, PAs_max_upper, ratio_min_maj, beam_PA, refpix_RA, refpix_Dec)

    pixRAs_max = np.concatenate(([refpix_RA], pixRAs_max))
    pixRAs_max_lower = np.concatenate(([refpix_RA], pixRAs_max_lower))
    pixRAs_max_upper = np.concatenate(([refpix_RA], pixRAs_max_upper))
    pixDecs_max = np.concatenate(([refpix_Dec], pixDecs_max))
    pixDecs_max_lower = np.concatenate(([refpix_Dec], pixDecs_max_lower))
    pixDecs_max_upper = np.concatenate(([refpix_Dec], pixDecs_max_upper))
    return pixRAs_max, pixDecs_max, pixRAs_max_lower, pixDecs_max_lower, pixRAs_max_upper, pixDecs_max_upper
    
def jet_ridgeline(fitsimage, how_many_rms=7, how_many_sigma=4):
    """
    Input parameters
    ----------------
    how_many_rms : int
        The significance of detection in the fits image.
    how_many_sigma : int
        This parameter is used to estimate the position angle range, in which the jet ridgeline resides at how_many_sigma significance.
    """
    data, header = _measure.read_fits_image(fitsimage)
    refpix_RA, refpix_Dec = _measure.find_the_pixel_with_the_highest_flux_density(data)
    #refpix_RA  = header['CRPIX1'] ## reference pixel in RA direction
    #refpix_Dec = header['CRPIX2'] ## reference pixel in Dec direction
    t_RL = t_ridgeline = _measure.obtain_jet_ridgeline(fitsimage, how_many_rms, how_many_sigma)
    pixRAs_max, pixDecs_max, pixRAs_max_lower, pixDecs_max_lower, pixRAs_max_upper, pixDecs_max_upper = prepare_locations_in_the_image(t_RL, header, refpix_RA, refpix_Dec)
    
    fig = plt.figure()
    plt.imshow(data, cmap='rainbow', norm=SymLogNorm(1e-2, base=10), origin='lower')
    #plt.imshow(data, cmap='cividis', norm=LogNorm(), origin='lower')
    
    plt.plot(pixRAs_max[:-1], pixDecs_max[:-1], color='white', linewidth=1)
    plt.arrow(pixRAs_max[-2], pixDecs_max[-2], pixRAs_max[-1]-pixRAs_max[-2], pixDecs_max[-1]-pixDecs_max[-2], color='white', linewidth=0.3, width=1.5, head_width=6, head_length=12, fill=True, length_includes_head=True)
    plt.plot(pixRAs_max_lower, pixDecs_max_lower, '--', color='white', linewidth=0.5)
    plt.plot(pixRAs_max_upper, pixDecs_max_upper, '--', color='white', linewidth=0.5)
    

    cbar = plt.colorbar()
    cbar.set_label(r'flux density ($\mathrm{Jy~{beam}^{-1}}$)', rotation=-90, labelpad=15)
    #plt.plot(pixRAs_max, pixDecs_max, color='white')
    
    ## >> now plot synthesized beam 
    beam_maj, beam_min, beam_PA = _measure.beam_info(header)
    ax = plt.gca()
    ellipse = Ellipse(xy=(header['NAXIS1']/10, header['NAXIS2']/10), width=beam_min, height=beam_maj, angle=beam_PA, fill=True, color='grey')
    ax.add_patch(ellipse)
    ## <<

    plt.xlabel(r'$-\Delta\alpha\,(\mathrm{pixel}; 1\,\mathrm{pixel}=%.2f\,\mathrm{mas}$)' % abs(float(header['CDELT1']*3.6e6)))
    plt.ylabel(r'$\Delta\delta\,(\mathrm{pixel}; 1\,\mathrm{pixel}=%.2f\,\mathrm{mas}$)' % float(header['CDELT2']*3.6e6))
    RA0 = _others.deg2dms(header['CRVAL1'] / 15., 3)
    Dec0 = _others.deg2dms(header['CRVAL2'], 2)
    pixRA0 = header['CRPIX1']
    pixDec0 = header['CRPIX2']
    plt.title(r'%s, %s at (%d, %d)' % (RA0, Dec0, pixRA0, pixDec0), fontsize=15, pad=10)

    fig.tight_layout()
    fitsimagename = fitsimage.split('/')[-1]
    outputpdf = fitsimagename.replace('clean.fits', 'jet_ridgeline.pdf')
    plt.savefig(outputpdf)
    #plt.show()
    plt.clf()
    plt.close()
