import matplotlib
matplotlib.use('Agg')
import numpy
from astropy.cosmology import Planck15 as cosmo
import matplotlib.pyplot as plot
from joblib import Parallel, delayed

#constants
c = 3.064e-7 #lightspeed in Mpc/yr  WA: c in Mpc/yr
G = 4.49E-33 #gravitational constant in Mpc, Mo, yr  WA: Gravitational constant in Megaparsec^3/solar mass/year^2

#functions to calculate min and max frequency of inspiral binary (1psc to ISCO)
def fmin(m1,m2,z):
    y = 100./(2.*numpy.pi*3.154*(1.+z))*numpy.sqrt(G*(m1+m2))
    return y
def fmax(m1,m2,z):
    y = 1./(6.*numpy.pi*numpy.sqrt(6.)*(1.+z)*3.154e+7)*c**3./(G*(m1+m2))
    return y

#########
###enrico
f = open('Klein16_PopIII.dat','r')
lines = f.readlines()
f.close()

Z = []
M1 = []
M2 = []
N = {}
#maxFreq = {}
#minFreq = {}
for line in lines:
    #redshifts
    z = float(line.split(' ')[1])
    Z.append(z)

    #Mass of each blackhole
    m1 = float(line.split(' ')[2])
    m2 = float(line.split(' ')[3])
    M1.append(m1)
    M2.append(m2)

    #calculate min and mass frequency for that specific total mass
    #maxFreq['%s-%s-%s'%(m1,m2,z)] = fmax(m1,m2,z)
    #minFreq['%s-%s-%s'%(m1,m2,z)] = fmin(m1,m2,z)

    #reading the comoving densities
    if '%s-%s-%s' %(m1,m2,z) in N.keys():
        N['%s-%s-%s'%(m1,m2,z)] = N['%s-%s-%s'%(m1,m2,z)] + float(line.split(' ')[-1].split('\n')[0])
    else:
        N['%s-%s-%s'%(m1,m2,z)] = float(line.split(' ')[-1].split('\n')[0])

#redshift bins with extrapolation for the last one
zBins = Z
zBins.append(20.00)
zBins = numpy.sort(zBins)
zBins = numpy.unique(zBins)

#unique m1,m2
M1 = numpy.unique(M1)
M2 = numpy.unique(M2)

#binning the M1 and M2
M1[-1] = M1[-1]+1.
M2[-1] = M2[-1]+1.

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

f = open('Klein16_PopIII.dat','r')
lines = f.readlines()
f.close()

NBinned = {}
#maxFreq = {}
#minFreq = {}
for line in lines:
    #redshifts
    z = float(line.split(' ')[1])

    #Mass of each blackhole
    m1 = float(line.split(' ')[2])
    m2 = float(line.split(' ')[3])

    #finding the bin numbers for m1 and m2
    bn1 = int((numpy.log10(m1)-M1logBins[0])/M1BinSize)
    m1 = 10.**((M1logBins[bn1]+M1logBins[bn1+1])/2.)

    bn2 = int((numpy.log10(m2)-M2logBins[0])/M2BinSize)
    m2 = 10.**((M2logBins[bn2]+M2logBins[bn2+1])/2.)

    #reading the comoving densities
    if '%s-%s-%2.2f'%(bn1,bn2,z) in NBinned.keys():
        NBinned['%s-%s-%2.2f'%(bn1,bn2,z)] = NBinned['%s-%s-%2.2f'%(bn1,bn2,z)] + float(line.split(' ')[-1].split('\n')[0])
    else:
        NBinned['%s-%s-%2.2f'%(bn1,bn2,z)] = float(line.split(' ')[-1].split('\n')[0])

#defining the inner intgral
def integralM(m1,m2,f):
    #c = counter
    bn1 = int((numpy.log10(m1)-M1logBins[0])/M1BinSize)
    bn2 = int((numpy.log10(m2)-M2logBins[0])/M2BinSize)
    summ = 0.
    for i in range(len(zBins)-1):
        if '%s-%s-%2.2f' %(bn1,bn2,zBins[i]) in NBinned.keys() and fmin(m1,m2,zBins[i]) < f < fmax(m1,m2,zBins[i]):
            summ = summ + NBinned['%s-%s-%2.2f'%(bn1,bn2,zBins[i])]/( (1.+(zBins[i]+zBins[i+1])/2.)**(1./3.) )
        #else: 
            #print(zBins[i])
    return summ

#defining the chirp mass
def mChirp(m1,m2):
    return (m1*m2)**(3./5.)/(m1+m2)**(1./5.)

#defining the h(f) for each frequency
def h(f):
    sum = 0.
    #i = 0
    for m1log in M1logBins:
        m1 = 10.**(m1log+M1BinSize/2.)
        m1start = 10.**(m1log)
        m1end = 10.**(m1log+M1BinSize)
        deltam1 = m1end - m1start
        for m2log in M2logBins:
            m2 = 10.**(m2log+M2BinSize/2.)
            m2start = 10.**(m2log)
            m2end = 10.**(m2log+M2BinSize)
            deltam2 = m2end - m2start
            sum = sum + integralM(m1,m2,f)*mChirp(m1,m2)**(5./3.)
            #print('%s-%s'%(i,sum))
        #i = i+1
    y = ( 4.*G**(5./3.)/(3.*numpy.pi**(1./3.)*c**2.) )**(0.5)*numpy.sqrt(sum)*f**(-2./3.) * (3.17e-8)**(2./3.)
    return f,y

#hc = {}
#freq = numpy.arange(0.001,0.1,0.005)
#for f in freq:
#    hc['%s'%f] = h(f)

freq = numpy.arange(-5.,4.2,0.2)
freq = 10.**(freq)

lines = []
for f in freq:
    x,y = h(f)
    lines.append('%s %s\n'%(x,y))

f = open('klein16_popIII.txt','w')
f.writelines(lines)
f.close()

################
################
################
###Q3-delay
f = open('Klein16_Q3delays.dat','r')
lines = f.readlines()
f.close()

Z = []
M1 = []
M2 = []
N = {}
#maxFreq = {}
#minFreq = {}
for line in lines:
    #redshifts
    z = float(line.split(' ')[1])
    Z.append(z)

    #Mass of each blackhole
    m1 = float(line.split(' ')[2])
    m2 = float(line.split(' ')[3])
    M1.append(m1)
    M2.append(m2)

    #calculate min and mass frequency for that specific total mass
    #maxFreq['%s-%s-%s'%(m1,m2,z)] = fmax(m1,m2,z)
    #minFreq['%s-%s-%s'%(m1,m2,z)] = fmin(m1,m2,z)

    #reading the comoving densities
    if '%s-%s-%s' %(m1,m2,z) in N.keys():
        N['%s-%s-%s'%(m1,m2,z)] = N['%s-%s-%s'%(m1,m2,z)] + float(line.split(' ')[-1].split('\n')[0])
    else:
        N['%s-%s-%s'%(m1,m2,z)] = float(line.split(' ')[-1].split('\n')[0])

#redshift bins with extrapolation for the last one
zBins = Z
zBins.append(20.00)
zBins = numpy.sort(zBins)
zBins = numpy.unique(zBins)

#unique m1,m2
M1 = numpy.unique(M1)
M2 = numpy.unique(M2)

#binning the M1 and M2
M1[-1] = M1[-1]+1.
M2[-1] = M2[-1]+1.

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

f = open('Klein16_Q3delays.dat','r')
lines = f.readlines()
f.close()

NBinned = {}
#maxFreq = {}
#minFreq = {}
for line in lines:
    #redshifts
    z = float(line.split(' ')[1])

    #Mass of each blackhole
    m1 = float(line.split(' ')[2])
    m2 = float(line.split(' ')[3])

    #finding the bin numbers for m1 and m2
    bn1 = int((numpy.log10(m1)-M1logBins[0])/M1BinSize)
    m1 = 10.**((M1logBins[bn1]+M1logBins[bn1+1])/2.)

    bn2 = int((numpy.log10(m2)-M2logBins[0])/M2BinSize)
    m2 = 10.**((M2logBins[bn2]+M2logBins[bn2+1])/2.)

    #reading the comoving densities
    if '%s-%s-%2.2f'%(bn1,bn2,z) in NBinned.keys():
        NBinned['%s-%s-%2.2f'%(bn1,bn2,z)] = NBinned['%s-%s-%2.2f'%(bn1,bn2,z)] + float(line.split(' ')[-1].split('\n')[0])
    else:
        NBinned['%s-%s-%2.2f'%(bn1,bn2,z)] = float(line.split(' ')[-1].split('\n')[0])

#hc = {}
#freq = numpy.arange(0.001,0.1,0.005)
#for f in freq:
#    hc['%s'%f] = h(f)

freq = numpy.arange(-5.,4.2,0.2)
freq = 10.**(freq)

lines = []
for f in freq:
    x,y = h(f)
    lines.append('%s %s\n'%(x,y))

f = open('klein16_Q3delays.txt','w')
f.writelines(lines)
f.close()

##################
##################
###Q3-noDelay
f = open('Klein16_Q3nodelays.dat','r')
lines = f.readlines()
f.close()

Z = []
M1 = []
M2 = []
N = {}
#maxFreq = {}
#minFreq = {}
for line in lines:
    #redshifts
    z = float(line.split(' ')[1])
    Z.append(z)

    #Mass of each blackhole
    m1 = float(line.split(' ')[2])
    m2 = float(line.split(' ')[3])
    M1.append(m1)
    M2.append(m2)

    #calculate min and mass frequency for that specific total mass
    #maxFreq['%s-%s-%s'%(m1,m2,z)] = fmax(m1,m2,z)
    #minFreq['%s-%s-%s'%(m1,m2,z)] = fmin(m1,m2,z)

    #reading the comoving densities
    if '%s-%s-%s' %(m1,m2,z) in N.keys():
        N['%s-%s-%s'%(m1,m2,z)] = N['%s-%s-%s'%(m1,m2,z)] + float(line.split(' ')[-1].split('\n')[0])
    else:
        N['%s-%s-%s'%(m1,m2,z)] = float(line.split(' ')[-1].split('\n')[0])

#redshift bins with extrapolation for the last one
zBins = Z
zBins.append(20.00)
zBins = numpy.sort(zBins)
zBins = numpy.unique(zBins)

#unique m1,m2
M1 = numpy.unique(M1)
M2 = numpy.unique(M2)

#binning the M1 and M2
M1[-1] = M1[-1]+1.
M2[-1] = M2[-1]+1.

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

f = open('Klein16_Q3nodelays.dat','r')
lines = f.readlines()
f.close()

NBinned = {}
#maxFreq = {}
#minFreq = {}
for line in lines:
    #redshifts
    z = float(line.split(' ')[1])

    #Mass of each blackhole
    m1 = float(line.split(' ')[2])
    m2 = float(line.split(' ')[3])

    #finding the bin numbers for m1 and m2
    bn1 = int((numpy.log10(m1)-M1logBins[0])/M1BinSize)
    m1 = 10.**((M1logBins[bn1]+M1logBins[bn1+1])/2.)

    bn2 = int((numpy.log10(m2)-M2logBins[0])/M2BinSize)
    m2 = 10.**((M2logBins[bn2]+M2logBins[bn2+1])/2.)

    #reading the comoving densities
    if '%s-%s-%2.2f'%(bn1,bn2,z) in NBinned.keys():
        NBinned['%s-%s-%2.2f'%(bn1,bn2,z)] = NBinned['%s-%s-%2.2f'%(bn1,bn2,z)] + float(line.split(' ')[-1].split('\n')[0])
    else:
        NBinned['%s-%s-%2.2f'%(bn1,bn2,z)] = float(line.split(' ')[-1].split('\n')[0])

#hc = {}
#freq = numpy.arange(0.001,0.1,0.005)
#for f in freq:
#    hc['%s'%f] = h(f)

freq = numpy.arange(-5.,4.2,0.2)
freq = 10.**(freq)

lines = []
for f in freq:
    x,y = h(f)
    lines.append('%s %s\n'%(x,y))

f = open('klein16_Q3nodelays.txt','w')
f.writelines(lines)
f.close()


