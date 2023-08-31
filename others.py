#!/usr/bin/env python
import numpy as np
import os, sys, math


def diffposition(RAs, refRA, Decs, refDec):
    diffDecs = dms2deg(Decs) - dms2deg(refDec)
    diffRAs = dms2deg(RAs) - dms2deg(refRA)
    diffDecs_mas = diffDecs*3600*1000
    diffRAs_ms = diffRAs*3600*1000
    refDec_rad = dms2deg(refDec)*math.pi/180
    diffRAs_mas = diffRAs_ms*15*math.cos(refDec_rad)
    return diffRAs_mas, diffDecs_mas

def deg2dms(array, decimal_place=7):
    if type(array)==float or type(array)==np.float64:
        dmss = deg_float2dms(array, decimal_place)
    else:
        dmss = np.array([])
        for number in array:
            dms = deg_float2dms(number, decimal_place)
            dmss = np.append(dmss, dms)
        if len(dmss) == 1:
            dmss = dmss[0]
    return dmss
def deg_float2dms(degree, decimal_place): #-- degree to dd:mm:ss.sssssss
    sign = np.sign(degree)
    degree = float(abs(degree)) 
    d = np.floor(degree)
    res = (degree - d) * 60 ## in arcmin
    m = np.floor(res)
    s = (res - m) * 60 ## in arcsec
    if 60 - s < 1e-7:
        s = 0
        m += 1
    dms="%02d:%02d:%0*.*f" % (d , m, int(decimal_place+3), int(decimal_place), s)
    if sign == -1:
        dms = '-' + dms
    return dms 

def dms_str2deg(string): #-- convert dd:mm:ss.ssss to dd.dddd
    if ':' in string:
        a = string.split(':')
    else:
        a = string.split(' ')
    sign1 = 1
    if a[0].strip().startswith('-'):
        sign1 = -1
    b = []
    for item in a:
        b.append(float(item))
    if b[0] < 0:
        sign = -1
    elif b[0] == 0 and sign1 == -1:
        sign = -1
    else:
        sign = 1
    b[0] = abs(b[0])
    i=len(b)
    degree = b[i-1]
    while i != 1:
        degree /= 60
        i -= 1
        degree += b[i-1]
    degree *= sign
    #degree = ((b[2]/60+b[1])/60+b[0])*sign
    return degree
def dms2deg(array):
    if type(array)==str:
        degrees = dms_str2deg(array)
    else:
        degrees = np.array([])
        for string in array:
            degree = dms_str2deg(string)
            degrees = np.append(degrees, degree)
    return degrees
def sample2median(array1):
    array = sorted(array1)
    length = len(array)
    if length % 2 == 0:
        median = 0.5*(array[length//2-1] + array[length//2])
    else:
        median = array[(length-1)//2]
    return median
def sample2median_range(array1, confidencelevel):
    array = sorted(array1)
    CL = confidencelevel
    if CL<1:
        SV = int(CL*len(array)) #SV -> significant volume
    elif CL>=1 and CL<10:
        CL = math.erf(CL/2**0.5)
        SV = int(CL*len(array))
    elif CL>=10:
        SV = CL
    index_start = int((len(array)-SV)/2-1)
    index_end = index_start + SV
    return array[index_start], array[index_end]
