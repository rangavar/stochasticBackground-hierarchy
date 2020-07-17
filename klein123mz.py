import matplotlib
matplotlib.use('Agg')
import numpy as np
import os
import h5py
import operator
import pickle
import random
import matplotlib.pyplot as plot
from astropy.cosmology import Planck15 as cosmo
import numpy
import os
import PhenomA as pa
import LISA as li
import WaveformTools as wt

font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 17}

#reading the mass and redshifts from simulations and meanwhile calucalting the snr in filtering for snr>7.0
f = open('data/Klein16_PopIII.dat','r')
lines = f.readlines()
f.close()

#calculate the eventRates per redshfit bin per year
def eventRate(N,z1,z2):
    return (1./(z2-z1))*4.*numpy.pi*N*(3.064e-7)*cosmo.comoving_distance((z1+z2)/2.).value**2.

Z = []
M1 = []
M2 = []
detections = {}
for line in lines:
    z = float(line.split(' ')[1])
    Z.append(z)
    m1 = float(line.split(' ')[2])
    M1.append(m1)
    m2 = float(line.split(' ')[3])
    M2.append(m2)  
    
    if '%s-%s-%s'%(m1,m2,z) in detections.keys():
        detections['%s-%s-%s'%(m1,m2,z)] = detections['%s-%s-%s'%(m1,m2,z)] + float(line.split(' ')[-1].split('\n')[0])
    else:
        detections['%s-%s-%s'%(m1,m2,z)] = float(line.split(' ')[-1].split('\n')[0])

#binning z
zBinsSize = 1.0
zBins = numpy.arange(0.0,20.+zBinsSize,zBinsSize)

#didnt sort cause already sorted
M1 = numpy.unique(M1)
M2 = numpy.unique(M2)

#binning the M1 and M2
M1[-1] = M1[-1]+1.
M2[-1] = M2[-1]+1.

Mtotal = [M1[0]+M2[0],M1[-1]+M2[-1]]
Mtotallog = np.log10(Mtotal)

M1log = numpy.log10(M1)
M2log = numpy.log10(M2)

M1BinSize = (M1log[-1]-M1log[0])/20.
M1logBins = []
for i in range(21):
    M1logBins.append(M1log[0]+M1BinSize*i)

M2BinSize = (M2log[-1]-M2log[0])/20.
M2logBins = []
for i in range(21):
    M2logBins.append(M2log[0]+M2BinSize*i)

MTBinSize = (Mtotallog[-1]-Mtotallog[0])/10.
MTlogBins = []
MTlogs = []
for i in range(11):
    MTlogBins.append(Mtotallog[0]+MTBinSize*i)

MZ = [Mtotal[0]*(1.+zBins[0]),Mtotal[1]*(1.+zBins[-1])]
MZlog = np.log10(MZ)
MZBinSize = (MZlog[-1]-MZlog[0])/10.
MZlogBins = []
MZlogs = []
for i in range(11):
    MZlogBins.append(MZlog[0]+MZBinSize*i)

#binning finished, now start
###popIII
f = open('data/Klein16_PopIII.dat','r')
lines = f.readlines()
f.close()

detections = {}
for line in lines:
    z = float(line.split(' ')[1])
    bnZ = int((z-zBins[0])/zBinsSize)
    m1 = float(line.split(' ')[2])
    m2 = float(line.split(' ')[3])   
    
    #bn1 = int((numpy.log10(m1)-M1logBins[0])/M1BinSize)
    #bn2 = int((numpy.log10(m2)-M2logBins[0])/M2BinSize)
    bnZ = int((z-zBins[0])/zBinsSize)

    mt = m1 + m2
    bnt = int((np.log10(mt)-MTlogBins[0])/MTBinSize)
    
    mz = mt * (1.+z)
    bnmz = int((np.log10(mz)-MZlogBins[0])/MZBinSize)

    if '%s-%s'%(bnmz,bnZ) in detections.keys():
        detections['%s-%s'%(bnmz,bnZ)] = detections['%s-%s'%(bnmz,bnZ)] + float(line.split(' ')[-1].split('\n')[0])
    else:
        detections['%s-%s'%(bnmz,bnZ)] = float(line.split(' ')[-1].split('\n')[0])

#event rate
eventRates = {}
for mz in range(len(MZlogBins)):
    summ = 0.0
    for z in range(len(zBins)):
        if '%s-%s'%(mz,z) in detections.keys():
            summ = summ + eventRate(detections['%s-%s'%(mz,z)],zBins[z],zBins[z+1])
    eventRates[MZlogBins[mz]] = summ

#plotting
x = []
y = []
for key in eventRates.keys():
    x.append(float(key))
x = numpy.sort(x)
for key in x:
    y.append(numpy.log10(eventRates[key]))
plot.plot(x,y, label='popIII')

###Q3delay
f = open('data/Klein16_Q3delays.dat','r')
lines = f.readlines()
f.close()

detections = {}
for line in lines:
    z = float(line.split(' ')[1])
    bnZ = int((z-zBins[0])/zBinsSize)
    m1 = float(line.split(' ')[2])
    m2 = float(line.split(' ')[3])   
    
    #bn1 = int((numpy.log10(m1)-M1logBins[0])/M1BinSize)
    #bn2 = int((numpy.log10(m2)-M2logBins[0])/M2BinSize)
    bnZ = int((z-zBins[0])/zBinsSize)

    mt = m1 + m2
    bnt = int((np.log10(mt)-MTlogBins[0])/MTBinSize)
    
    mz = mt * (1.+z)
    bnmz = int((np.log10(mz)-MZlogBins[0])/MZBinSize)

    if '%s-%s'%(bnmz,bnZ) in detections.keys():
        detections['%s-%s'%(bnmz,bnZ)] = detections['%s-%s'%(bnmz,bnZ)] + float(line.split(' ')[-1].split('\n')[0])
    else:
        detections['%s-%s'%(bnmz,bnZ)] = float(line.split(' ')[-1].split('\n')[0])

#event rate
eventRates = {}
for mz in range(len(MZlogBins)):
    summ = 0.0
    for z in range(len(zBins)):
        if '%s-%s'%(mz,z) in detections.keys():
            summ = summ + eventRate(detections['%s-%s'%(mz,z)],zBins[z],zBins[z+1])
    eventRates[MZlogBins[mz]] = summ

#plotting
x = []
y = []
for key in eventRates.keys():
    x.append(float(key))
x = numpy.sort(x)
for key in x:
    y.append(numpy.log10(eventRates[key]))
plot.plot(x,y, label='Q3delay')

###Q3noDelay
f = open('data/Klein16_Q3nodelays.dat','r')
lines = f.readlines()
f.close()

detections = {}
for line in lines:
    z = float(line.split(' ')[1])
    bnZ = int((z-zBins[0])/zBinsSize)
    m1 = float(line.split(' ')[2])
    m2 = float(line.split(' ')[3]) 
    
    #bn1 = int((numpy.log10(m1)-M1logBins[0])/M1BinSize)
    #bn2 = int((numpy.log10(m2)-M2logBins[0])/M2BinSize)
    bnZ = int((z-zBins[0])/zBinsSize)

    mt = m1 + m2
    bnt = int((np.log10(mt)-MTlogBins[0])/MTBinSize)
    
    mz = mt * (1.+z)
    bnmz = int((np.log10(mz)-MZlogBins[0])/MZBinSize)

    if '%s-%s'%(bnmz,bnZ) in detections.keys():
        detections['%s-%s'%(bnmz,bnZ)] = detections['%s-%s'%(bnmz,bnZ)] + float(line.split(' ')[-1].split('\n')[0])
    else:
        detections['%s-%s'%(bnmz,bnZ)] = float(line.split(' ')[-1].split('\n')[0])

#event rate
eventRates = {}
for mz in range(len(MZlogBins)):
    summ = 0.0
    for z in range(len(zBins)):
        if '%s-%s'%(mz,z) in detections.keys():
            summ = summ + eventRate(detections['%s-%s'%(mz,z)],zBins[z],zBins[z+1])
    eventRates[MZlogBins[mz]] = summ
    
#plotting
x = []
y = []
for key in eventRates.keys():
    x.append(float(key))
x = numpy.sort(x)
for key in x:
    y.append(numpy.log10(eventRates[key]))
plot.plot(x,y, label='Q3noDelay')

plot.xlabel('$log(M_z) \;\;[M_{\odot}]$')
plot.ylabel('$log(d^2N/dtdM_z) \;\;[year]^{-1}$')
plot.legend()
plot.savefig('plots/eventRate_mz_klein.eps',format='eps',dpi=1000)
plot.close()
