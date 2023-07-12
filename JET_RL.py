# -*- coding: utf-8 -*-
"""
Created on 20230712

@author: 
line 10-58 by Xiaofeng Li
the rest by Wei Zhao 
"""
from astropy.io import fits
import numpy as np
import scipy.ndimage
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import matplotlib
matplotlib.rcParams.update({'font.size': 14})

class Beam:
	def __init__(self,x0=0.0,y0=0.0,bmaj=0.0,bmin=0.0,bpa=0.0):
		self.maj = bmaj
		self.min = bmin
		self.pa = bpa
		self.x0 = x0
		self.y0 = y0
	def show(self):
		print(self.maj,self.min,self.pa)
		print(self.x0,self.y0)
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
	flux = scipy.ndimage.map_coordinates(np.transpose(img),np.vstack((X,Y))) #determine the flux density along the slices
        flux_max=flux[np.argmax(flux)]                                           #determine the max flux density along the slices
        ang_max=ang[np.argmax(flux)]                                             #determine the position angle of the max 
        xmax=x[np.argmax(flux)]                                                  #determine the position of the max  in R.A.
        ymax=y[np.argmax(flux)]                                                  #determine the position of the max  in Dec
        doc=open("RL_data.txt",'a')                
        print>>doc,xmax,ymax,flux_max,ang_max
        doc.close() 



fname = 'hello_jet.fits'
hdu = fits.open(fname)
h = hdu[0].header

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
