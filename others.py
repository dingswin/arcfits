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
def periodic_sample2estimate(alist, period=360, confidencelevel=1):
    """
    periodic sample cannot use median as the estimate.
    instead, it should adopt the most compact symmetric confidence interval
    as the estimate and the associated uncertainty.

    Output parameters
    -----------------
    median : flaot
    upper-side error of median : float
    lower-side error of median : float
    """
    alist = np.sort(np.array(alist))
    alist = alist % period
    CL = confidencelevel
    error = float('inf')
    median = None
    value = None
    for threshold in np.arange(0,period,period/360.):
        reordered_list = move_elements_larger_than_a_threshold_to_the_head_of_a_list(alist, threshold)
        value1, error1, median1 = sample2estimate_and_median(reordered_list, CL)
        if error1 < error:
            error = error1
            value = value1
            median = median1
    return median % period, value+error-median, median-(value-error)
def move_elements_larger_than_a_threshold_to_the_head_of_a_list(a_sorted_list, threshold, period=360):
    """
    a_sorted_list should be a numpy array
    """
    aSL = a_sorted_list
    length = len(aSL)
    index_threshold = (np.abs(aSL - threshold)).argmin()
    new_head_indice = np.arange(index_threshold, length)
    list_new_head = aSL[new_head_indice] - period
    list_new = np.concatenate((list_new_head, aSL[np.arange(index_threshold)]))
    return list_new
def sample2estimate_and_median(array1,confidencelevel):
    """
    find the narrowest confidence interval and report the median of the this interval
    """
    array = sorted(array1)
    CL = confidencelevel
    if CL<1:
        SV = int(CL*len(array)) #SV -> significant volume
    elif CL>=1 and CL<10:
        CL = math.erf(CL/2**0.5)
        SV = int(CL*len(array))
    elif CL>=10:
        SV = CL    
    delta = float('inf')
    for i in range(len(array)-SV-1):
        diff = array[i+SV] - array[i]
        if diff < delta:
            j=i
            delta = diff
    confidence_min = array[j]
    confidence_max = array[j+SV]
    value = 0.5*(confidence_min+confidence_max)
    error = 0.5*(confidence_max-confidence_min)
    return value,error,array[int(j+SV/2.)]

class lsqfit: #y=ax+b
    def __init__(s):
        #Xs = np.array([1,2,4,7,10])
        #errXs = np.array([0.1,0.3,0.3,0.5,0.2])
        #Ys = np.array([1,3,6,2,4])
        #errYs = errXs
        #print s.linearfit2D(Xs, errXs, Ys, errYs)
        #[a, b, r_chi_sq, run] = s.linearfit2D(Xs, errXs, Ys, errYs)
        #s.plot_linearfit2D(a, b, Xs, Ys, errXs, errYs)
        pass

    def linearfit2D(s, Xs, errXs, Ys, errYs):
        DoF = len(Xs) - 2 #degree of freedom
        X1s = Xs
        chi_sq_old = float('inf') 
        chi_sq = 1e10
        run = 0
        while abs(chi_sq/chi_sq_old)<1-1e-5:
            run += 1
            chi_sq_old = chi_sq
            Delta = sum(1./errYs**2)*sum(X1s**2/errYs**2) 
            Delta -= (sum(X1s/errYs**2))**2
            b = sum(X1s**2/errYs**2)*sum(Ys/errYs**2) - sum(X1s/errYs**2)*sum(X1s*Ys/errYs**2)
            b /= Delta
            a = sum(1/errYs**2)*sum(X1s*Ys/errYs**2) - sum(X1s/errYs**2)*sum(Ys/errYs**2)
            a /= Delta
            chi_sq = sum((X1s-Xs)**2/errXs**2 + (Ys-a*X1s-b)**2/errYs**2)
            X1s = Xs/errXs**2 + a*(Ys-b)/errYs**2
            X1s /= 1/errXs**2 + a**2/errYs**2
        r_chi_sq = chi_sq/DoF
        return a, b, r_chi_sq, run
    def linearfit_no_errors(s, Xs, Ys): ## y=a*x+b
        n = len(Xs)
        a = (n * np.sum(Xs * Ys) - np.sum(Xs) * np.sum(Ys)) / (n * np.sum(Xs**2) - np.sum(Xs)**2)
        b = (np.sum(Ys) - a * np.sum(Xs)) / n
        return a, b
    def linearfit_with_asymmetric_X_errs(s, Xs, X_upper_errs, X_lower_errs, Ys, N=1e4, HowManySigma=1):
        n = len(Xs)
        N = int(N)
        X_sims = np.zeros((N,n))
        for i in range(n):
            X_sims[:,i] = split_normal_random_sampler(Xs[i], X_upper_errs[i], X_lower_errs[i], N)
        As, Bs = np.array([]), np.array([])
        for j in range(N):
            A, B = s.linearfit_no_errors(X_sims[j,:], Ys)
            As = np.append(As, A)
            Bs = np.append(Bs, B)
        A_median = np.median(As)
        A_lower, A_upper = sample2median_range(As, HowManySigma)
        B_median = np.median(Bs)
        B_lower, B_upper = sample2median_range(Bs, HowManySigma)
        print('A = %f + %f - %f (at %dsigma confidence).' % (A_median, A_upper-A_median, A_median-A_lower, HowManySigma))
        print('B = %f + %f - %f (at %dsigma confidence).' % (B_median, B_upper-B_median, B_median-B_lower, HowManySigma))
        return A_median, A_lower, A_upper, B_median, B_lower, B_upper
    
    def linearfit_with_asymmetric_errs(s, Xs, Ys, N=1e4, HowManySigma=1, **kwargs):
        n = len(Xs)
        N = int(N)
        having_X_errs, having_Y_errs = False, False
        try:
            X_upper_errs = kwargs['X_upper_errs']
            X_lower_errs = kwargs['X_lower_errs']
            having_X_errs = True
            X_sims = np.zeros((N,n))
            for i in range(n):
                X_sims[:,i] = split_normal_random_sampler(Xs[i], X_upper_errs[i], X_lower_errs[i], N)
        except KeyError:
            pass
       
        try:
            Y_upper_errs = kwargs['Y_upper_errs']
            Y_lower_errs = kwargs['Y_lower_errs']
            having_Y_errs = True
            Y_sims = np.zeros((N,n))
            for i in range(n):
                Y_sims[:,i] = split_normal_random_sampler(Ys[i], Y_upper_errs[i], Y_lower_errs[i], N)
        except KeyError:
            pass
        
        if having_X_errs==False and (having_Y_errs==False):
            print('carrying out linear fitting without any error...')
            A, B = s.linearfit_no_errors(Xs, Ys)
            return A, B
        else:
            As, Bs = np.array([]), np.array([])
            if having_X_errs and (not having_Y_errs):
                for j in range(N):
                    A, B = s.linearfit_no_errors(X_sims[j,:], Ys)
                    As = np.append(As, A)
                    Bs = np.append(Bs, B)
            elif having_Y_errs and (not having_X_errs):
                for j in range(N):
                    A, B = s.linearfit_no_errors(Xs, Y_sims[j,:])
                    As = np.append(As, A)
                    Bs = np.append(Bs, B)
            elif having_X_errs and having_Y_errs:
                for j in range(N):
                    A, B = s.linearfit_no_errors(X_sims[j,:], Y_sims[j,:])
                    As = np.append(As, A)
                    Bs = np.append(Bs, B)
                
            A_median = np.median(As)
            A_lower, A_upper = sample2median_range(As, HowManySigma)
            B_median = np.median(Bs)
            B_lower, B_upper = sample2median_range(Bs, HowManySigma)
            print('A = %f + %f - %f (at %dsigma confidence).' % (A_median, A_upper-A_median, A_median-A_lower, HowManySigma))
            print('B = %f + %f - %f (at %dsigma confidence).' % (B_median, B_upper-B_median, B_median-B_lower, HowManySigma))
            return A_median, A_lower, A_upper, B_median, B_lower, B_upper


def split_normal_random_sampler(mu, upper_sigma, lower_sigma, size):
    """
    input parameters
    ----------------
    mu : float 
        mode of the distributin
    upper_sigma : float
        upper 1-sigma error
    lower_sigma : float
        lower 1-sigma error, positive only!
    size : int
        size of the simulations

    outputs
    -------
    numbers : list of float
    """
    size = int(size)
    mu, lower_sigma, upper_sigma = np.float64([mu, lower_sigma, upper_sigma])
    numbers = np.abs(np.random.normal(0, 1, size))
    prob = np.array([upper_sigma, lower_sigma]) / (upper_sigma + lower_sigma) ## the PDF has to be smooth around the mode, which is ensured by prob.
    numbers *= np.random.choice([upper_sigma, -lower_sigma], size, p=prob)
    numbers += mu
    return numbers
