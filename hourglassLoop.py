import matplotlib
matplotlib.use('Agg')
import numpy as np
import os
import h5py
import operator
import pickle
import random
import matplotlib.pyplot as plt
from astropy.cosmology import Planck15 as cosmo
import numpy
import PhenomA as pa
import LISA as li
import WaveformTools as wt
import os
import pickle as pickleRick

font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 15}

matplotlib.rc('font', **font)

#defining the error range and the loop size on each error as well
allErrors = np.arange(0.0,0.305,0.005)
loopSize = 500

#plotting the LISA sensitivity line for future use
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
lisa = li.LISA()
f  = np.logspace(np.log10(1.0e-5), np.log10(100.), 1000)
Sn = lisa.Sn(f)
li.PlotSensitivityCurve(f, Sn)
plt.plot(f,Sn)

#reading the pre-made templates (1,2,3)
templateStrain = {}
#1(popIII)
f = open('../0template/klein16_popIII_LISA.txt','r')
lines = f.readlines()
f.close()
fred = []
strain = []
for line in lines:
    freq = float(line.split(' ')[0])
    templateStrain['1-%1.8f'%freq] = float(line.split(' ')[1])
    fred.append(float(line.split(' ')[0]))
    strain.append(float(line.split(' ')[1]))
plt.plot(fred,strain,label='popIII(1)')
#2(Q3delay)
f = open('../0template/klein16_Q3delays_LISA.txt','r')
lines = f.readlines()
f.close()
fred = []
strain = []
for line in lines:
    freq = float(line.split(' ')[0])
    templateStrain['2-%1.8f'%freq] = float(line.split(' ')[1])
    fred.append(float(line.split(' ')[0]))
    strain.append(float(line.split(' ')[1]))
plt.plot(fred,strain,label='delay(2)')
#3(Q3noDelay)
f = open('../0template/klein16_Q3nodelays_LISA.txt','r')
lines = f.readlines()
f.close()
fred = []
strain = []
for line in lines:
    freq = float(line.split(' ')[0])
    templateStrain['3-%1.8f'%freq] = float(line.split(' ')[1])
    fred.append(float(line.split(' ')[0]))
    strain.append(float(line.split(' ')[1]))
plt.plot(fred,strain,label='noDelay(3)')

#defining chi squared (simple gaussian)
def xaiSquared(sObserved,sTemplate,scatter):
    y = ((sObserved-sTemplate)**2.)/(scatter**2.)
    return y

#define probability from xai2
def prob(xai):
    p = np.exp( -xai/2. )
    return p

#defining the template numbers and frequency
templateNumber = ['1','2','3']

#this is new: defining the function that compares the likelihood ratios between the true type and the second most probable given the data
def pRatioFunction(array,true):
    tempArray = np.sort(array)
    a = array[int(true)-1]
    b = tempArray[-1]
    if a == b:
        b = tempArray[-2]
    y = a/b
    return y

trueModelsAddress = ['klein16_popIII_LISA.txt','klein16_Q3delays_LISA.txt','klein16_Q3nodelays_LISA.txt']
pRatio = {}
for error in allErrors:
    #assigning which model is the true model number (the one we use to make the sims)
    for trueModelNumber in ['1','2','3']:

        for i in range(loopSize):
            #reading one of the strains to simulate/sample the observations
            sampleStrain = {}
            sampleFreq = []
            fred = []
            strain = []
            f = open('../0template/%s'%trueModelsAddress[int(trueModelNumber)-1],'r')
            lines = f.readlines()
            f.close()
            for line in lines:
                sampleStrain['%1.8f'%freq] = (float(line.split(' ')[1]))
                sampleFreq.append(float(line.split(' ')[0]))
                fred.append(float(line.split(' ')[0]))
                strain.append(float(line.split(' ')[1]))

            frequency = sampleFreq

            #adding the noise model
            count = 0
            for f in fred:
                workStrain = np.log10(strain[count])
                workStrain = np.random.normal(workStrain,-workStrain*error)
                #workStrain = workStrain + workStrain*error
                strain[count] = 10.**(workStrain)
                sampleStrain['%1.8f'%f] = strain[count]
                count = count + 1

            #we just want to use the data from LISA: e-4 t0 2e-2
            frequency = []
            for f in sampleFreq:
                if 0.0001 < f < 0.02:
                    frequency.append(f)

            #strain maximum likelihood algorithm
            probes = {}
            sumProbes = 0.
            for tn in templateNumber:
                sumXai = 0.
                for f in frequency:
                    #print(tn,f)

                    xai2 = xaiSquared(np.log10(sampleStrain['%1.8f'%f]),np.log10(templateStrain['%s-%1.8f'%(tn,f)]),np.log10(sampleStrain['%1.8f'%f])*error)
                    sumXai = sumXai + xai2

                probability = prob(sumXai)
                probes['%s'%tn] = probability
                sumProbes = sumProbes + probability

            #normalizing the probaility
            nProbes = {}
            plotProbes = []
            for tn in templateNumber:
                nProbes['%s'%tn] = probes['%s'%tn]/sumProbes
                plotProbes.append(nProbes['%s'%tn])

            #saving the normalized probability
            #lines = []
            #lines.append('probability distribution from strain only\n')
            #for tn in templateNumber:
            #    lines.append('%s-%s\n'%(tn,nProbes['%s'%tn]))

            #f = open('outputs/LStrain_loop_%s_%2.1f.txt'%(trueModelNumber,error*100),'w')
            #f.writelines(lines)
            #f.close()

            if '%s-%1.3f'%(trueModelNumber,error) not in pRatio.keys():
                pRatio['%s-%1.3f'%(trueModelNumber,error)] = pRatioFunction(plotProbes,trueModelNumber)
            else:
                pRatio['%s-%1.3f'%(trueModelNumber,error)] = pRatio['%s-%1.3f'%(trueModelNumber,error)] + pRatioFunction(plotProbes,trueModelNumber)
    
            print('%s:error-itteration = %1.3f - %s' %(trueModelNumber,error,i))

#now that the loop is over, normalization
for key in pRatio.keys():
    pRatio[key] = pRatio[key]/float(loopSize)

#saving the dictionary for tomorrow morning
pickle_out = open("outputs/loop.pickle","wb")
pickleRick.dump(pRatio, pickle_out)
pickle_out.close()

#load the dictionary
pickle_in = open("outputs/loop.pickle","rb")
pRatio = pickleRick.load(pickle_in)

#plotting
plotRatio1 = []
plotRatio2 = []
plotRatio3 = []
for error in allErrors:
    plotRatio1.append(pRatio['%s-%1.3f'%('1',error)])
    plotRatio2.append(pRatio['%s-%1.3f'%('2',error)])
    plotRatio3.append(pRatio['%s-%1.3f'%('3',error)])

plotRatio1[0] = plotRatio1[1]
for i in range(len(plotRatio1)):
    if numpy.isnan(plotRatio1[i]):
        print('yes') 
        plotRatio1[i] = 0.000000001

plotRatio2[0] = plotRatio2[1]
for i in range(len(plotRatio2)):
    if numpy.isnan(plotRatio2[i]):
        print('yes') 
        plotRatio2[i] = 0.000000001

plotRatio3[0] = plotRatio3[1]
for i in range(len(plotRatio3)):
    if numpy.isnan(plotRatio3[i]):
        print('yes') 
        plotRatio3[i] = 0.0000000001

plotRatio1 = np.log10(plotRatio1)
plotRatio2 = np.log10(plotRatio2)
plotRatio3 = np.log10(plotRatio3)

plt.plot(allErrors,plotRatio1,label = 'popIII')
plt.plot(allErrors,plotRatio2,label = 'Q3delay')
plt.plot(allErrors,plotRatio3,label = 'Q3noDelay')
plt.axhline(np.log10(5.0),c='m')
plt.legend()
plt.ylim(0,10)
plt.xlim(0,0.3)
plt.xlabel('noise ($\epsilon$)')
plt.ylabel('$log(R_p)$')
plt.savefig('plots/pRatio.eps',format='eps',dpi=1000)
plt.close()
