# -*- coding: utf-8 -*-
"""
Created on Fri Dec 18 19:21:19 2015

@author: thorsten
"""

import numpy as np

def DFT_time2freq( t, val, freq, signal_type='pulse'):
    assert len(t)==len(val)
    assert len(freq)>0
    f_val = np.zeros(len(freq))*1j
    for n_f in range(len(freq)):
        f_val[n_f] = np.sum( val*np.exp( -1j*2*np.pi*freq[n_f] * t ) )

    if signal_type == 'pulse':
        f_val *= t[1]-t[0]
    elif signal_type == 'periodic':
        f_val /= len(t)
    else:
        raise Exception('Unknown signal type: "{}"'.format(signal_type))

    return 2*f_val  # single-sided spectrum

def Check_Array_Equal(a,b, tol, relative=False):
    a = np.array(a)
    b = np.array(b)
    if a.shape!=b.shape:
        return False
    if tol==0:
        return (a==b).all()
    if relative:
        d = np.abs((a-b)/a)
    else:
        d = np.abs((a-b))
    return np.max(d)<tol

if __name__=="__main__":
    import pylab as plt

    t = np.linspace(0,2,201)

    s = np.sin(2*np.pi*2*t)
    plt.plot(t,s)

    f = np.linspace(0,3,101)
    sf = DFT_time2freq(t, s, f, 'periodic')

    plt.figure()
    plt.plot(f, np.abs(sf))

    plt.show()


