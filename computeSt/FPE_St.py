import numpy as np
import math
from pandas import read_csv as pd

from scipy.stats import norm
from scipy import linalg as LA
from scipy import interpolate
from scipy.signal import savgol_filter as sf

import torch
from torch import linalg as tLA
import os, sys

from time import time
from functions import *
from joblib import Parallel, delayed

print('Starting simulation')
start_time = time()
sys.stdout.flush()

filename = 'obelix_eq_st.txt'

eps = 1e-9
tag = 0
assert tag in [0, 1]

## the diffusion constant : https://www.degruyter.com/document/doi/10.1515/zna-1987-0915/pdf
dLvals = [1.85e-4, 2.9e-4]##diffusion constant, Angstrom^2/fs
dgvals = [1.4, 2.0]
mass = [3*16e-3, 17e-3]
factor = 4.184*1e-7

D = dLvals[tag]
dgas = dgvals[tag]
m = mass[tag]

kBT = 0.5962 ## kB*T @300K,  Kcal/mol
mg = D/kBT ## 1/(m*gamma)
gamma_liquid = (kBT/(m*D))*factor ### in 1/fs unit
gamma_gas = (kBT/(m*dgas))*factor ### in 1/fs unit

rt = 8.314*300
sig = 1e-5*math.sqrt(rt/m) ## in angstrom/fs unit
mu0 = sig

# %%
# mymol = ['O3', 'OH']
# mol = mymol[tag]

data = pd('freefile.txt', comment = '#', header = None, \
         delim_whitespace = True).to_numpy()
xdata = data[:,0]
xdata_shifted = xdata + abs(min(xdata))
ydata = data[:,1]

z = interpolate.interp1d( xdata_shifted, ydata)
fit = z

xpoints = int(1e5)
left_extend = 5
right_extend = 5
total_x_max = left_extend + 0.5*(max(xdata) - min(xdata)) + right_extend

xdata_in_use = np.linspace(0.0, total_x_max, xpoints)

L = max(abs(xdata_in_use)) - min(abs(xdata_in_use)) ##length of the box, in Angstrom
spi = math.sqrt(math.pi)
piL = math.pi/L
dx = abs(xdata_in_use[1]-xdata_in_use[0])
L2x = 2/L*dx

def one_sided_potential(x_val):
    if x_val <= left_extend:
        pot = fit(max(xdata_shifted)*0.5)
    elif x_val >= left_extend + 0.5*(max(xdata) - min(xdata)):
        pot = fit(max(xdata_shifted))
    else:
        pot = fit(x_val - left_extend + 0.5*max(xdata_shifted))
    return pot

fit_potential_in_use = np.array(list(map(one_sided_potential, xdata_in_use)))
v1p_in_use = np.gradient( fit_potential_in_use, xdata_in_use )

v1p_in_use = sf(v1p_in_use, int(xpoints/5)+1, 2)

nmax = 400
alpha = 4
alphamax = 2*alpha + 1
ntotal = nmax * alphamax

f1 = pd('optim_sm.txt', header = None, delim_whitespace=True).to_numpy()
min_list = np.argsort(f1[:,0])[:50]
width_x, shift_x, power = f1[min_list[0]][1:]
high_x = gamma_liquid

print('parameters for friction: width_x = %s, shift_x = %s, power = %s'%(width_x, shift_x, power))
sys.stdout.flush()
def gamma_one_sided(y):
    x = y - left_extend
    result = np.where(x > 0, (np.tanh(width_x * (-x + shift_x)) / 2 + 0.5)**power * (high_x - gamma_gas) + gamma_gas, (np.tanh(width_x * (x + shift_x)) / 2 + 0.5)**power * (high_x - gamma_gas) + gamma_gas)
    return result

gamma_arr_in_use = gamma_one_sided(xdata_in_use)
muarr = np.linspace(-2*mu0, 2*mu0+eps, alphamax)

## integrate( gauss(i)*gauss(j) )
def goverlap(i, j, sig = sig):
    return 1/2/sig/spi*math.exp(-0.25*((muarr[i]-muarr[j])/sig)**2)

def gi1(i, j, sig = sig):
    return 0.25/spi/sig*math.exp(-0.25*((muarr[i]-muarr[j])/sig)**2)*(muarr[i]+muarr[j])

def gi2(i, j, sig = sig):
    return 0.125/spi/(sig**3)*math.exp(-0.25*((muarr[i]-muarr[j])/sig)**2)*(-muarr[i]**2+muarr[j]**2+2*sig**2)

def gi3(i, j, sig = sig):
    return 0.25/spi/(sig**3)*math.exp(-0.25*((muarr[i]-muarr[j])/sig)**2)*(-muarr[i]+muarr[j])

def gi4(i, j, sig = sig):
    return 0.125/spi/(sig**5)*math.exp(-0.25*((muarr[i]-muarr[j])/sig)**2)*((muarr[i]-muarr[j])**2-2*sig**2)

## PIB wave function integrations
def integrateijp(i, j):
    return np.sum(np.sin(i*piL*xdata_in_use)*np.cos(j*piL*xdata_in_use))*L2x*j*piL

def integrateijV(i, j):
    ival = np.sum(np.sin(i*piL*xdata_in_use)*v1p_in_use*np.sin(j*piL*xdata_in_use))*L2x
    return ival


# testArr = np.linspace(26.0,62.0,10000)
check_pos = np.argmin(fit_potential_in_use)
mymu = xdata_in_use[check_pos] ## L/2
# mymu = left_extend + 35.0
sigma_initial = 0.3

print('initial position is', mymu)
sys.stdout.flush()
def initial_coeff_2(nmax_1, L = L, mu = mymu, sigma = sigma_initial):
    nmax_array = np.arange(1, nmax_1 + 1 + eps, 1)
    initial_coeff = np.zeros(nmax_1)
    gauss_x = np.asarray( gauss(xdata_in_use, mymu, sigma) )
    for i in range(len(initial_coeff)):
        initial_coeff[i] = np.sum(gauss_x*np.sin(nmax_array[i]*piL*np.asarray(xdata_in_use)))*dx*math.sqrt(2/L)
    return initial_coeff


initial_coeff_position = initial_coeff_2(nmax)
initial_coeff_velocity = np.zeros(alphamax)
initial_coeff_velocity[alpha] += 1
initial = np.kron( initial_coeff_position, initial_coeff_velocity)
# print('The non-zero initial matrix elements are',initial[np.where(initial.numpy()!=0.)]) ## print 
def test_initial_pos(x):
    checkin = nmax
    test1 = initial_coeff_2(checkin)
    temp1 = math.sqrt(2/L)*np.sin(math.pi/L*np.outer(np.arange(1,checkin+1,1),x))
    return np.sum(test1[:,None]*temp1, axis = 0)

vdata = np.linspace( -20*sig, 20*sig, 1000)
def test_initial_vel(v):
    checkin = alphamax
    test1 = initial_coeff_velocity
    temp1 = 0
    for i in range(len(muarr)):
        temp1 += gauss(vdata, muarr[i], sig)*test1[i]
    return temp1

# #### $S$-matrix and its inverse

# %%
overlap = np.zeros((alphamax, alphamax))
for i in range(alphamax):
    for j in range(alphamax):
        overlap[i, j] += goverlap(i, j)

smat = np.kron(np.eye(nmax), overlap)
sinv = LA.solve(smat, np.eye(len(smat)))

# %%
scale_potential = 1.0
def from_coeff_matrix(high_x = high_x, width_x = width_x, shift_x = shift_x):
    ## set up the A matrix : S \dot{C}(t)=A\cdot C(t)

    gamma_arr = gamma_one_sided(xdata_in_use)

    def integrate_gamma(i, j):
        return np.sum(np.sin(i*piL*xdata_in_use)*gamma_arr_in_use*np.sin(j*piL*xdata_in_use))*L2x

    amat = np.zeros((ntotal, ntotal))
    ## term 1, -v*D[p[x,v,t],x]
    xcoeff = np.zeros((nmax, nmax))
    vcoeff = np.zeros((alphamax, alphamax))
    for i in range(nmax):
        for j in range(nmax):
            xcoeff[i, j] += integrateijp(i+1, j+1)
    for i in range(alphamax):
        for j in range(alphamax):
            vcoeff[i, j] += gi1(i, j)
    amat -= np.kron(xcoeff, vcoeff)

    ## term 2, gamma(z)*D[v*p[x,v,t],v]
    xcoeff = np.zeros((nmax, nmax))
    vcoeff = np.zeros((alphamax, alphamax))
    for i in range(nmax):
        for j in range(nmax):
            xcoeff[i, j] += integrate_gamma(i+1, j+1)
    for i in range(alphamax):
        for j in range(alphamax):
            vcoeff[i, j] += gi2(i, j)
    amat += np.kron(xcoeff, vcoeff)

    # term 4, k_B*T*\gamma (z)/m *D[p[x,v,t],{v,2}]
    vcoeff = np.zeros(( alphamax, alphamax ))
    xcoeff = xcoeff*(kBT/m)*factor
    for i in range(alphamax):
        for j in range(alphamax):
            vcoeff[i, j] += gi4(i, j)
    amat += np.kron(xcoeff, vcoeff)

    # term 3, V'(x)/m*D[p[x,v,t],v]
    xcoeff = np.zeros((nmax, nmax))
    vcoeff = np.zeros((alphamax, alphamax))
    for i in range(nmax):
        for j in range(nmax):
            xcoeff[i, j] += integrateijV(i+1, j+1)*factor/m
    for i in range(alphamax):
        for j in range(alphamax):
            vcoeff[i, j] += gi3(i, j) 
    amat += np.kron(scale_potential*xcoeff, vcoeff)

    coeff_matrix = sinv@amat

    return coeff_matrix

# %%


# %%
coeff_matrix = chop(from_coeff_matrix(high_x = high_x, \
                                      width_x = width_x, shift_x = shift_x))

# %%
e_val, e_vec = LA.eig(coeff_matrix)
e_vec_inv = chop(LA.solve(e_vec, np.eye(e_vec.shape[0])))
def mymatexp(t):
    # return chop( e_vec @ np.diag(np.exp(e_val*t)) @ e_vec_inv )
    return torch.matrix_exp( torch.tensor(coeff_matrix*t) ).numpy()
    # return chop(LA.expm( coeff_matrix*t ))

# %%


# %%
# tmaxs = [1.5e7, 1e9]
tmaxs = [5e4,1e9]
tmax = tmaxs[tag]
numSpacing = 6000

vdata = np.linspace( -20*sig, 20*sig, len(xdata_in_use) )
checkx = xdata_in_use
dx = checkx[1] - checkx[0]

sin_integral = np.zeros(nmax)
for i in range(nmax):
    sin_integral[i] += np.sum(np.sin(piL*(i+1)*xdata_in_use))*math.sqrt(2/L)*dx

# tlist = np.append(0.0, np.geomspace(0.5, tmax, numSpacing))
tlist = np.linspace(0.0 + abs(np.random.normal()/100) , tmax + abs(np.random.normal()/100), numSpacing )

def cmn_part_vel():
    test_a = []
    for i in range(alphamax):
        test_a.append(gauss(vdata, muarr[i], sig))
    test_a = np.asarray(test_a)
    return np.kron( sin_integral, test_a.T )

cmn_int_1 = chop(np.kron(np.sin(piL*np.outer(np.arange(1,nmax+1,1), xdata_in_use)).T, np.ones(alphamax) ))
cmn_int_2 = chop(cmn_part_vel())

def checkPx(t) :
    ct = mymatexp(t).dot(initial) ## np.matmul(mymatexp(t), initial) ## 
    return np.real(cmn_int_1.dot(ct))*math.sqrt(2/L)

def checkPv(t):
    ct = mymatexp(t).dot(initial) ## mymatexp(t).dot(initial)
    return np.real(cmn_int_2.dot(ct))

def survival(t):
    return np.sum(checkPx(t))*dx

# %%
stationary_pt = find_nearest(e_val.real, 10.0)

print(e_val[stationary_pt])
sys.stdout.flush()
# %%

# %%
survival_prob_list = np.zeros( len(tlist) )
if os.path.isfile(filename) == False:
    f = open(filename, 'w')
    f.close()


def compute_st(ii):
   temp = survival(tlist[ii])
#    survival_prob_list[ii] += temp
   print('Done %s of %s calculations after %s seconds at t = %s fs where st = %s '%(ii+1, numSpacing, time()-start_time, tlist[ii], temp))
   sys.stdout.flush()
   f = open(filename, 'a')
   print('%20.8f %12.8f'%(tlist[ii], temp), file = f)
   # sys.stdout.flush()
   f.close()

cpu_count = 64
Parallel( n_jobs = cpu_count, verbose = len(tlist) )( delayed(compute_st)(ii) for ii in range(len(tlist)) )

# np.savetxt('st.txt', np.vstack( (tlist, survival_prob_list) ).T)

# f2.close()

print('It took %s seconds'%(time()-start_time))
