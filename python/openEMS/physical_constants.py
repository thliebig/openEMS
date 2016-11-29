# -*- coding: utf-8 -*-
"""
Created on Fri Dec 18 20:58:01 2015

@author: thorsten
"""

import numpy as np

C0 = 299792458 # m/s
MUE0 = 4e-7*np.pi # N/A^2
EPS0 = 1/(MUE0*C0**2) # F/m

# free space wave impedance
Z0 = np.sqrt(MUE0/EPS0) # Ohm
