import matplotlib
matplotlib.use('Agg')
import numpy
from astropy.cosmology import Planck15 as cosmo
import matplotlib.pyplot as plot
fig_size = plot.rcParams["figure.figsize"]
fig_size[0] = 12
fig_size[1] = 4
plot.rcParams["figure.figsize"] = fig_size

##
f = open('data/phyprop_bhmergers_cdm_nouv_j30_snr.dat','r')
lines = f.readlines()
f.close()

data = {}
data7 = {}
data1  = {}
data2 = {}
data3 = {}
number = {}
for line in lines:
    z = float(line.split('\t')[3])
    ty = float(line.split('\t')[4])
    snr = float(line.split('\t')[-1].split('\n')[0])

    if z in data.keys():
        data[z] = data[z] + (10.**(float(line.split('\t')[2])))
    else:
        data[z] = (10.**(float(line.split('\t')[2])))
    
    if snr > 7.0:
        if z in data7.keys():
            data7[z] = data7[z] + (10.**(float(line.split('\t')[2])))
        else:
            data7[z] = (10.**(float(line.split('\t')[2])))

        if ty == 1.:
            if z in data1.keys():
                data1[z] = data1[z] + (10.**(float(line.split('\t')[2])))
            else:
                data1[z] = (10.**(float(line.split('\t')[2])))

        if ty == 2.:
            if z in data2.keys():
                data2[z] = data2[z] + (10.**(float(line.split('\t')[2])))
            else:
                data2[z] = (10.**(float(line.split('\t')[2])))
    
        if ty == 3.:
            if z in data3.keys():
                data3[z] = data3[z] + (10.**(float(line.split('\t')[2])))
            else:
                data3[z] = (10.**(float(line.split('\t')[2])))
    
    if z in number.keys():
        number[z] = number[z] + 1.
    else:    
        number[z] = 1.

def eventRate(N,z1,z2):
    return (1./(z2-z1))*0.01*4.*numpy.pi*N*(3.064e-7)*cosmo.comoving_distance(z1).value**2.

zBins = []
for key in data.keys():
    zBins.append(key)
zBins.append(18.43)
zBins = numpy.sort(zBins)

eventRates = {}
for i in range(len(data.keys())):
    eventRates[zBins[i]] = eventRate(data[zBins[i]],zBins[i],zBins[i+1])

eventRates7 = {}
for i in range(len(data.keys())):
    if zBins[i] in data7.keys():
        eventRates7[zBins[i]] = eventRate(data7[zBins[i]],zBins[i],zBins[i+1])
    else:
        eventRates7[zBins[i]] = 0.00000000000001

eventRates1 = {}
for i in range(len(data.keys())):
    if zBins[i] in data1.keys():
        eventRates1[zBins[i]] = eventRate(data1[zBins[i]],zBins[i],zBins[i+1])
    else:
        eventRates1[zBins[i]] = 0.00000000000001

eventRates2 = {}
for i in range(len(data.keys())):
    if zBins[i] in data2.keys():
        eventRates2[zBins[i]] = eventRate(data2[zBins[i]],zBins[i],zBins[i+1])
    else:
        eventRates2[zBins[i]] = 0.00000000000001

eventRates3 = {}
for i in range(len(data.keys())):
    if zBins[i] in data3.keys():
        eventRates3[zBins[i]] = eventRate(data3[zBins[i]],zBins[i],zBins[i+1])
    else:
        eventRates3[zBins[i]] = 0.00000000000001

x = []
y = []
for key in eventRates.keys():
    x.append(float(key))
x = numpy.sort(x)
for key in x:
    y.append(numpy.log10(eventRates[key]))

x7 = []
y7 = []
for key in eventRates7.keys():
    x7.append(float(key))
x7 = numpy.sort(x7)
for key in x7:
    y7.append(numpy.log10(eventRates7[key]))

x1 = []
y1 = []
for key in eventRates1.keys():
    x1.append(float(key))
x1 = numpy.sort(x1)
for key in x:
    y1.append(numpy.log10(eventRates1[key]))

x2 = []
y2 = []
for key in eventRates2.keys():
    x2.append(float(key))
x2 = numpy.sort(x2)
for key in x2:
    y2.append(numpy.log10(eventRates2[key]))

x3 = []
y3 = []
for key in eventRates3.keys():
    x3.append(float(key))
x3 = numpy.sort(x3)
for key in x3:
    y3.append(numpy.log10(eventRates3[key]))

plot.subplot(122)
plot.plot(x,y,label = 'all events')
plot.plot(x7,y7, label = 'snr > 7')
plot.plot(x1,y1, label = 'type1 snr > 7')
plot.plot(x2,y2, label = 'type2 snr > 7' )
plot.plot(x3,y3, label = 'type3 snr > 7')
plot.ylim(-3.,2.)
plot.xlim(3,20)
plot.title('nouv j30')
plot.xlabel('$z$')
plot.ylabel('$log(d^2N/dtdz) [year]^{-1}$')
plot.legend()

f = open('data/Klein16_PopIII.dat','r')
lines = f.readlines()
f.close()

data = {}
number = {}
for line in lines:
    z = float(line.split(' ')[1])

    for i in range(len(zBins)-1):
        if zBins[i]<z<zBins[i+1]:
            if zBins[i] in data.keys():
                data[zBins[i]] = data[zBins[i]] + float(line.split(' ')[-1].split('\n')[0])
            else:
                data[zBins[i]] = float(line.split(' ')[-1].split('\n')[0])


def eventRate(N,z1,z2):
    return (1./(z2-z1))*4.*numpy.pi*N*(3.064e-7)*cosmo.comoving_distance(z1).value**2.

eventRates = {}
for i in range(len(data.keys())):
    eventRates[zBins[i]] = eventRate(data[zBins[i]],zBins[i],zBins[i+1])

x = []
y = []
for key in eventRates.keys():
    x.append(float(key))
x = numpy.sort(x)
for key in x:
    y.append(numpy.log10(eventRates[key]))

plot.subplot(121)
plot.plot(x,y)
plot.ylim(-3.,2.)
plot.xlim(3,20)
plot.xlabel('$z$')
plot.ylabel('$log(d^2N/dtdz) [year]^{-1}$')
plot.title('Klein16 PopIII')

#import os
#os.mkdir('plots')

plot.savefig('plots/eventRate_nouv30_vs_klein.eps',format='eps',dpi=1000)
plot.close()
