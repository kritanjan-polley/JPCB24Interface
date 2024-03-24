import numpy as np
import scipy.sparse as sps
import math
from pandas import read_csv as pd
import os

eps = 1e-9
left_cutoff = 10.0
right_cutoff = 40.0

cwd = os.getcwd()
cwd = cwd + '/'
data_dir = cwd
temp_dir = cwd
shooting_data_dir = cwd + 'data_dir/'
# print(cwd)
# print(shooting_data_dir)

## read data file
def read_data_file(pos, left_cutoff = left_cutoff, right_cutoff = right_cutoff):
    data_file = pd(shooting_data_dir+'data.'+str(pos)+'.shooting',header=0, delim_whitespace=True).to_numpy()

    zvals = abs(data_file[:,1])
    dt = data_file[2,0] - data_file[1,0]
    tvals = np.arange(0.0, (len(zvals))*dt + eps, dt)[:len(zvals)]
    combined = np.vstack((tvals,zvals)).T

    int_array = combined[combined[:,1] < right_cutoff]
    reducted_array = int_array[int_array[:,1] > left_cutoff]

    return reducted_array

## find nearest value in an array
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

## normal distribution
def gauss(x, mu, sig):
    return 1.0/math.sqrt(2.0*math.pi*sig*sig)*np.exp( -0.5*((x-mu)/sig)**2 )

## one sided gaussian
def one_sided_gaussian(xArr, mu, sig, side = 'L'):
    if side == 'L':
        heaviside = np.heaviside(xArr-mu, 0.0)
        flip_heaviside = np.where( heaviside == 0.0, 1.0, 0.0)
        return gauss(xArr, mu, sig)*flip_heaviside
    elif side == 'R':
        return gauss(xArr, mu, sig)*np.heaviside(xArr-mu, 0.0)

## chop function, similar to Mathematica
def chop(expr, delta=10**-12):
    return np.ma.masked_inside(expr, -delta, delta).filled(0)

## moving average
def moving_avg(x, n):
    cumsum = np.cumsum(np.insert(x, 0, 0))
    return (cumsum[n:] - cumsum[:-n]) / float(n)