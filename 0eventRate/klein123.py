import matplotlib
matplotlib.use('Agg')
import numpy
from astropy.cosmology import Planck15 as cosmo
import matplotlib.pyplot as plot

font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 14}

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
zBinsSize = 2.0
zBins = numpy.arange(0.0,20.+zBinsSize,zBinsSize)

#binning finished, now start
###popIII
f = open('data/Klein16_PopIII.dat','r')
lines = f.readlines()
f.close()

detections = {}
for line in lines:
    z = float(line.split(' ')[1])
    bnZ = int((z-zBins[0])/zBinsSize)
    if '%s'%(bnZ) in detections.keys():
        detections['%s'%(bnZ)] = detections['%s'%(bnZ)] + float(line.split(' ')[-1].split('\n')[0])
    else:
        detections['%s'%(bnZ)] = float(line.split(' ')[-1].split('\n')[0])

#event rate
eventRates = {}
for j in detections.keys():
    i = int(j)
    eventRates[zBins[i]] = eventRate(detections['%s'%i],zBins[i],zBins[i+1])

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
    if '%s'%(bnZ) in detections.keys():
        detections['%s'%(bnZ)] = detections['%s'%(bnZ)] + float(line.split(' ')[-1].split('\n')[0])
    else:
        detections['%s'%(bnZ)] = float(line.split(' ')[-1].split('\n')[0])

#event rate
eventRates = {}
for j in detections.keys():
    i = int(j)
    eventRates[zBins[i]] = eventRate(detections['%s'%i],zBins[i],zBins[i+1])

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
    if '%s'%(bnZ) in detections.keys():
        detections['%s'%(bnZ)] = detections['%s'%(bnZ)] + float(line.split(' ')[-1].split('\n')[0])
    else:
        detections['%s'%(bnZ)] = float(line.split(' ')[-1].split('\n')[0])

#event rate
eventRates = {}
for j in detections.keys():
    i = int(j)
    eventRates[zBins[i]] = eventRate(detections['%s'%i],zBins[i],zBins[i+1])

#plotting
x = []
y = []
for key in eventRates.keys():
    x.append(float(key))
x = numpy.sort(x)
for key in x:
    y.append(numpy.log10(eventRates[key]))
plot.plot(x,y, label='Q3noDelay')

plot.xlabel('$z$')
plot.ylabel('$log(d^2N/dtdz) \;\;[year]^{-1}$')
plot.legend()
plot.savefig('plots/eventRate_klein.eps',format='eps',dpi=1000)
plot.close()
