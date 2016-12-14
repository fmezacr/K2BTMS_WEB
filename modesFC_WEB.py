"""
Created on Mon Oct  3 11:39:09 2016
@author: fmeza
"""

# ojo: solo fits no .DS ademas mejorar lo de solo una freq.

# Libraries

from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import ascii
import os
import bisect
import argparse


from ivs.timeseries import pergrams
from ivs.timeseries import freqanalyse as fa

import matplotlib
matplotlib.rcParams.update({'font.size': 16})

if os.path.exists("outputWEB"):
	os.system("rmdir outputWEB")

cat0 = 0
cat1 = 0
cat2 = 0
cat3 = 0
cat4 = 0
total = 0 

parser = argparse.ArgumentParser()
parser.add_argument("filename")
args = parser.parse_args()

filename = args.filename

print "Processing candidate star #{0}...".format(total+1)
star = filename[0:9]	
#star = args.star
os.system("mkdir outputWEB")	
os.system("mkdir outputWEB/{0}".format(star))	

# Reading Files

hdulist = fits.open("input/{0}".format(filename));
data = hdulist[1].data
plotx = 'TIME'
ploty = 'FLUX'

cols = hdulist[1].columns

xdata = data.field(plotx)
ydatanan = data.field(ploty)

ydata = np.nan_to_num(ydatanan)	

# ASCII Creation, from .fits to .dat for Light Curve (lc)

ascii.write([xdata, ydata], 'outputWEB/{0}/{0}_lc.dat'.format(star),format='no_header')

# Display Light Curve

plt.figure(figsize=(18,6))
plt.xlabel('Time (d)')
plt.ylabel('Relative Flux (ppm)')
plt.title('LIGHT CURVE KASOC - EPIC {0}'.format(star))
plt.plot(xdata,ydata,'k')
#plt.xlim(xmin=0, xmax=80)
plt.savefig('outputWEB/{0}/{0}_lc.png'.format(star))

# Calculate Periodogram FOURIER

times,signal = data['TIME'],data['FLUX']
#times,signal = data['TIME'],data['FLUX']-np.average(data['FLUX'])
newsignal = np.nan_to_num(signal)

#freqs,ampls = pergrams.scargle(times,signal,f0=0, fn=1, df=0.1, norm='amplitude')
freqs,ampls = pergrams.scargle(times,newsignal)

# ASCII Creation, from .fits to .dat for Lomb-Scargle periodogram (lsp)

ascii.write([freqs, ampls], 'outputWEB/{0}/{0}_lsp.dat'.format(star),format='no_header')

# Display Periodogram FOURIER 

#plt.subplot(212)
plt.figure(figsize=(18,6))
plt.xlabel('Frequency [c d$^{-1}$]')
plt.ylabel('Amplitude [magnitude]')
plt.title('LOMB SCARGLE PERGRAM - EPIC {0}'.format(star))
#plt.axis([-10, 100, -0.000002, 0.000002])
plt.plot(freqs,ampls,'k')
plt.savefig('outputWEB/{0}/{0}_lsp.png'.format(star))

# Iterative Pre-whitening, SNR > 4

filename = 'outputWEB/{0}/{0}_lc.dat'.format(star)
maxiter = 1000
#f0=0.
#fn=20.

data = np.loadtxt(filename)
times,signal = data[:,0],data[:,1]-np.average(data[:,1])
#params,model = fa.iterative_prewhitening(times,signal,f0=f0,fn=fn,maxiter=maxiter,threads='max',full_outputWEB=True,stopcrit=(fa.stopcrit_scargle_snr,4.,1.))
#params = fa.iterative_prewhitening(freqs,ampls,f0=2,fn=7,maxiter=6,stopcrit=(fa.stopcrit_scargle_snr,4.,1,))
params,model = fa.iterative_prewhitening(times,signal,maxiter=maxiter,threads='max',full_output=True,stopcrit=(fa.stopcrit_scargle_snr,4.,1.),prewhiteningorder_snr=True,prewhiteningorder_snr_window=1.)
#print pl.mlab.rec2txt(params,precision=6)


paramsSave = np.vstack((params['const'],params['ampl'],params['e_ampl'],params['freq'],params['e_freq'],params['phase'],params['e_phase'],params['stopcrit']))
paramsSave = paramsSave.T
ascii.write(paramsSave, 'outputWEB/{0}/{0}_lc_pars.dat'.format(star),format='no_header')

# Aca va un IF con una salida para el caso de solo una freq. If freq = 0 then make frequecies cero y sale del loop
if len(paramsSave) > 1 :
# Frequency Combination Avoiding	
	snrlimit = 4
	frequlimit_min = 0.5
	frequlimit_max = 25
	parentfrequlimit_min = 0.5
	parentfrequlimit_max = 25
	
	combo_maxlenght = 2
	#combo_multiplicator_limits = [-1,1]
	combo_multiplicator_limits = [0,8]
	
	data = np.loadtxt('outputWEB/{0}/{0}_lc.dat'.format(star))
	times = data[:,0]
	signal = data[:,1]
	
	timerange = times[-1]-times[0]
	rayleigh = 1./timerange
	
	params = np.loadtxt('outputWEB/{0}/{0}_lc_pars.dat'.format(star))
	params = params[:-1,(0,1,3,5,7)] #because the last one is not matching the SNR crit anymore
	cons_ = params[:,0]
	ampl_ = params[:,1]
	freq_ = params[:,2]
	phas_ = params[:,3]
	snrs_ = params[:,4]
	cons = cons_[(snrs_>=snrlimit) & (freq_>=frequlimit_min) & (freq_<=frequlimit_max)]
	amplitudes = ampl_[(snrs_>=snrlimit) & (freq_>=frequlimit_min) & (freq_<=frequlimit_max)]
	frequencies = freq_[(snrs_>=snrlimit) & (freq_>=frequlimit_min) & (freq_<=frequlimit_max)]
	phas = phas_[(snrs_>=snrlimit) & (freq_>=frequlimit_min) & (freq_<=frequlimit_max)]
	snrs = snrs_[(snrs_>=snrlimit) & (freq_>=frequlimit_min) & (freq_<=frequlimit_max)]
	#print 'Total number of frequencies:',len(snrs)
	
	c_2 = 0 #combination frequency counters
	c_3 = 0
	c_4 = 0
	c_5 = 0
	c_5p= 0

else: 
	frequencies = [] 	

# Creation of the mainfreqs.txt file.

#	txt = open("outputWEB/{0}/{0}_mainfreqs.txt".format(star),"w")
#	print >> txt, plt.mlab.rec2txt(frequencies,precision=6)
#	txt.close()

# Display main frequancies, SNR > 4

plt.figure()
plt.figure(figsize=(18,6))
plt.vlines(frequencies,0,amplitudes,color='k',linestyle='-')
plt.xlabel('Frequency [c d$^{-1}$]')
plt.ylabel('Amplitude [magnitude]')
plt.title('ALL MODES (SNR > 4) - EPIC {0}'.format(star))
#plt.annotate('$\ell=2$',(6.46,0.015),xycoords='data',va='bottom',ha='center',fontsize='large')
#plt.annotate('$\ell=0$',(3,0.294),xycoords='data',va='bottom',ha='center',fontsize='large')
#plt.annotate('$\ell=1$',(6.99,0.015),xycoords='data',va='bottom',ha='center',fontsize='large')
plt.xlim(0,25)
#plt.ylim(0,1500)
plt.savefig('outputWEB/{0}/{0}_modes.png'.format(star))

plt.figure()
plt.vlines(frequencies,0,amplitudes,color='k',linestyle='-')
plt.xlabel('Frequency [c d$^{-1}$]')
plt.ylabel('Amplitude [magnitude]')
plt.title('HO G-MODES - EPIC {0}'.format(star))
#plt.annotate('$\ell=2$',(6.46,0.015),xycoords='data',va='bottom',ha='center',fontsize='large')
#plt.annotate('$\ell=0$',(3,0.294),xycoords='data',va='bottom',ha='center',fontsize='large')
#plt.annotate('$\ell=1$',(6.99,0.015),xycoords='data',va='bottom',ha='center',fontsize='large')
plt.xlim(0,2.99)
#plt.ylim(0,1500)
plt.savefig('outputWEB/{0}/{0}_hogmodes.png'.format(star))

plt.figure()
plt.vlines(frequencies,0,amplitudes,color='k',linestyle='-')
plt.xlabel('Frequency [c d$^{-1}$]')
plt.ylabel('Amplitude [magnitude]')
plt.title('LO G/P-MODES - EPIC {0}'.format(star))
#plt.annotate('$\ell=2$',(6.46,0.015),xycoords='data',va='bottom',ha='center',fontsize='large')
#plt.annotate('$\ell=0$',(3,0.294),xycoords='data',va='bottom',ha='center',fontsize='large')
#plt.annotate('$\ell=1$',(6.99,0.015),xycoords='data',va='bottom',ha='center',fontsize='large')
plt.xlim(3,5.99)
#plt.ylim(0,1500)
plt.savefig('outputWEB/{0}/{0}_logpmodes.png'.format(star))

plt.figure()
plt.vlines(frequencies,0,amplitudes,color='k',linestyle='-')
plt.xlabel('Frequency [c d$^{-1}$]')
plt.ylabel('Amplitude [magnitude]')
plt.title('HO P-MODES - EPIC {0}'.format(star))
#plt.annotate('$\ell=2$',(6.46,0.015),xycoords='data',va='bottom',ha='center',fontsize='large')
#plt.annotate('$\ell=0$',(3,0.294),xycoords='data',va='bottom',ha='center',fontsize='large')
#plt.annotate('$\ell=1$',(6.99,0.015),xycoords='data',va='bottom',ha='center',fontsize='large')
plt.xlim(6,25)
#plt.ylim(0,25)
plt.savefig('outputWEB/{0}/{0}_hopmodes.png'.format(star))

plt.close("all")

# Modes counting

modes = frequencies 

source = sorted(modes) # problema cuando modes = 0 

poles = (0,3,6,25)

outputWEB = [source[bisect.bisect_left(source,poles[i]):bisect.bisect_right(source,poles[i+1])] for i in range(len(poles)-1)]

hog = len(outputWEB[0])
logp = len(outputWEB[1])
hop = len(outputWEB[2])

B1 = "0" 
B2 = "0"
B3 = "0"

for f in modes:
	if 0 <= f <3 : 
		B1 = "1"
	if 3 <= f < 6:
		B2 = "1"
	if 6 <= f < 25 :
		B3 = "1"

comb = B1 + B2 + B3

if comb == "000" :
	cat = "no-modes"
	cat4 += 1
elif comb == "100" :
	cat = "SPB"
	cat1 += 1
elif comb == "011" :
	cat = "beta Cep"
	cat2 += 1
elif comb == "111" :
	cat = "hybrid"
	cat3 += 1
else: 
	cat = "OTHERS" 	
	cat0 += 1
	
os.system("touch outputWEB/{0}/{1}".format(star,cat))
total += 1

# Journal PDF Creation
 	
from reportlab.pdfgen import canvas
c = canvas.Canvas('outputWEB/{0}/{0}.pdf'.format(star))
c.drawString(250,790, "EPIC {0}".format(star))
c.drawImage('outputWEB/{0}/{0}_lc.png'.format(star),0,620,600,150)
c.drawImage('outputWEB/{0}/{0}_lsp.png'.format(star),0,460,600,150)
c.drawImage('outputWEB/{0}/{0}_modes.png'.format(star),0,300,600,150)
c.drawImage('outputWEB/{0}/{0}_hogmodes.png'.format(star),60,140,170,150)
c.drawImage('outputWEB/{0}/{0}_logpmodes.png'.format(star),230,140,170,150)
c.drawImage('outputWEB/{0}/{0}_hopmodes.png'.format(star),400,140,170,150)
c.drawString(90,110, "EPIC {0} has ({1}) HO g-modes, ({2}) LO g- and p-modes and ({3}) HO p-modes.".format(star,hog,logp,hop) )
c.drawString(90,90, "EPIC {0} ===> {1} candidate.".format(star,cat) )
c.showPage()
c.save()





		
