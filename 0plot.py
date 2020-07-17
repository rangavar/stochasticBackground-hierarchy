import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import PhenomA as pa
import LISA as li
import WaveformTools as wt
import os

#os.mkdir('plots')

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
lisa = li.LISA() 
f  = np.logspace(np.log10(1.0e-5), np.log10(100.), 1000)
Sn = lisa.Sn(f)
li.PlotSensitivityCurve(f, Sn)
plt.plot(f,Sn)

f = open('../1energyDensity/sn_ept1_j30/phyprop_bhmergers_cdm_nouv_j30.txt','r')
lines = f.readlines()
f.close()
freq = []
strain = []
for line in lines:
    freq.append(float(line.split(' ')[0]))
    strain.append(float(line.split(' ')[1]))
plt.plot(freq,strain,label='nouv j30')

f = open('../1energyDensity/sn_ept1_j30_tdf/phyprop_bhmergers_cdm_nouv_tdf_j30.txt','r')
lines = f.readlines()
f.close()
freq = []
strain = []
for line in lines:
    freq.append(float(line.split(' ')[0]))
    strain.append(float(line.split(' ')[1]))
plt.plot(freq,strain,label='nouv tdf j30')

f = open('../1energyDensity/snuv_ept1_j300/phyprop_bhmergers_cdm_uv_j300.txt','r')
lines = f.readlines()
f.close()
freq = []
strain = []
for line in lines:
    freq.append(float(line.split(' ')[0]))
    strain.append(float(line.split(' ')[1]))
plt.plot(freq,strain,label='uv j300')

f = open('../1energyDensity/snuv_ept1_j300_tdf/phyprop_bhmergers_cdm_uv_tdf_j300.txt','r')
lines = f.readlines()
f.close()
freq = []
strain = []
for line in lines:
    freq.append(float(line.split(' ')[0]))
    strain.append(float(line.split(' ')[1]))
plt.plot(freq,strain,label='uv tdf j300')

f = open('../1energyDensity/klein16_popIII/klein16_popIII.txt','r')
lines = f.readlines()
f.close()
freq = []
strain = []
for line in lines:
    freq.append(float(line.split(' ')[0]))
    strain.append(float(line.split(' ')[1]))
plt.plot(freq,strain,label='klein popIII')

f = open('../1energyDensity/klein16_popIII/klein16_Q3nodelays.txt','r')
lines = f.readlines()
f.close()
freq = []
strain = []
for line in lines:
    freq.append(float(line.split(' ')[0]))
    strain.append(float(line.split(' ')[1]))
plt.plot(freq,strain,label='klein Q3nodelays')

f = open('../1energyDensity/klein16_popIII/klein16_Q3delays.txt','r')
lines = f.readlines()
f.close()
freq = []
strain = []
for line in lines:
    freq.append(float(line.split(' ')[0]))
    strain.append(float(line.split(' ')[1]))
plt.plot(freq,strain,label='klein Q3delays')

f= open('decigo.csv','r')
lines = f.readlines()
f.close()

freq = []
strain = []
for line in lines[1:]:
    freq.append(10.**float(line.split(',')[0]))
    strain.append(10.**float(line.split(',')[1]))
plt.plot(freq,strain,label='decigo')

plt.xlim(10**(-5),100)
plt.ylim(10**(-26),10**(-15))
plt.legend()
plt.xlabel('log(f) - log(Hz)')
plt.ylabel('strain')
#plt.legend()
plt.savefig('plots/strain.eps',format='eps',dpi=1000)
