import matplotlib
matplotlib.use('Agg')
import numpy as np
import os
import h5py
import operator
import pickle
import random
import matplotlib.pyplot as plot

#reading the pre-made templates (1,2,3)
templateStrain = {}
#1(popIII)
f = open('../0template/klein16_popIII.txt','r')
lines = f.readlines()
f.close()
for line in lines:
    freq = float(line.split(' ')[0])
    templateStrain['1-%1.8f'%freq] = float(line.split(' ')[1])
#2(Q3delay)
f = open('../0template/klein16_Q3delays.txt','r')
lines = f.readlines()
f.close()
for line in lines:
    freq = float(line.split(' ')[0])
    templateStrain['2-%1.8f'%freq] = float(line.split(' ')[1])
#3(Q3noDelay)
f = open('../0template/klein16_Q3nodelays.txt','r')
lines = f.readlines()
f.close()
for line in lines:
    freq = float(line.split(' ')[0])
    templateStrain['3-%1.8f'%freq] = float(line.split(' ')[1])

#reading a sample of strain to work with for now 
#misiing the noise model for now!
sampleFreq = []
sampleStrain = {}
f = open('../0template/klein16_Q3nodelays.txt','r')
lines = f.readlines()
f.close()
for line in lines:
    sampleFreq.append(float(line.split(' ')[0]))
    freq = float(line.split(' ')[0])
    sampleStrain['%1.8f'%freq] = (float(line.split(' ')[1]))

#somewhat of a noise model
noise = 10.**-18.
for f in sampleFreq:
    sampleStrain['%1.8f'%f] = sampleStrain['%1.8f'%f] + (np.random.rand(1)[0]-0.5)*noise

#defining the template numbers and frequency
templateNumber = [1,2,3]
frequency = sampleFreq

#defining chi squared
def xaiSquared(sObserved,sTemplate,scatter):
    y = ((sObserved-sTemplate)**2.)/(scatter**2.)
    return y

#define probability from xai2
def prob(xai):
    p = np.exp( -xai/2. )
    return p

#main ML algorithm
probes = {}
sumProbes = 0.
for tn in templateNumber:
    sumXai = 0.
    for f in frequency: 
        print(tn,f)

        xai2 = xaiSquared(sampleStrain['%1.8f'%f],templateStrain['%s-%1.8f'%(tn,f)],noise)
        sumXai = sumXai + xai2

    probability = prob(sumXai)
    probes['%s'%tn] = probability
    sumProbes = sumProbes + probability

nProbes = {}
plotProbes = []
for tn in templateNumber:
    nProbes['%s'%tn] = probes['%s'%tn]/sumProbes
    plotProbes.append(nProbes['%s'%tn])

lines = []
for tn in templateNumber:
    lines.append('%s-%s'%(tn,nProbes['%s'%tn]))

f = open('likelihood.txt','w')
f.writelines(lines)
f.close()

plot.plot(templateNumber,plotProbes)
#plot.axvline(x='3')
plot.xlabel('simulation')
plot.ylabel('probability')
plot.savefig('probPlots.png')
plot.close()
