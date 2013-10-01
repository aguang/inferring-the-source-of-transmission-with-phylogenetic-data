"""
OUTPUT FOR PHI=10PC, ST 1
number of simulations: 527
pdonor_isdonor Wed Aug 28 17:32:09 2013
number of pairs: 13368 Wed Aug 28 17:32:10 2013
number of transmissions: 825 Wed Aug 28 17:32:10 2013
expected number of transmissions: 800 Wed Aug 28 17:32:10 2013
(array(-0.9681825836275303), 0.33297070847435867)
local regression data Wed Aug 28 17:32:11 2013
slope, intercept, r_value, p_value, slope_std_error
1.01011477629 0.00119479710972 0.448790564001 0.0 0.0173975224172
[ 1.01488066]


OUTPUT FOR PHI=10PC, ST 0
number of simulations: 631
pdonor_isdonor Thu Aug 29 10:11:53 2013
number of pairs: 45826 Thu Aug 29 10:12:06 2013
number of transmissions: 2343 Thu Aug 29 10:12:06 2013
expected number of transmissions: 2310 Thu Aug 29 10:12:06 2013
(array(-0.7857331374982653), 0.43202790782957123)
local regression data Thu Aug 29 10:12:07 2013
slope, intercept, r_value, p_value, slope_std_error
1.0270028229 -0.000648055165958 0.470882256054 0.0 0.00898830769319
Optimization terminated successfully.
         Current function value: 1730.270244
         Iterations: 1
         Function evaluations: 9
         Gradient evaluations: 3
[ 1.02443841]


"""
from pylab import *
import cPickle, csv, pdb, time, os, sys

print 'start', time.ctime()


phiflag, stflag = sys.argv[1:3]
ostem = 'phi%s_st%s_' % ( phiflag, stflag)
infns = [ fn for fn in os.listdir('output') if ostem in fn]
ftdata = list()
for fn in infns:
	d = cPickle.load(open('output/%s' % fn))
	ftdata.extend(d)
	print 'loaded cpickle', time.ctime()
#

radius = 25#50

print 'number of simulations:', len(ftdata)

pdonor_isdonor = list()
for d in ftdata:
	#~ nwk, sampleTimes, sampleStates, tip2donortip,  Cs, Ps, P, tipnames, o_prevalence = d
	nwk, sampleTimes, sampleStates, tip2donortip,  W,  o_prevalence = d
	
	if len(sys.argv)>3:
		#make a plot
		t = linspace(0., sampleTimes.values()[0] , len(o_prevalence))
		plot(t, o_prevalence[:,0], 'k-')
	#
	
	P = W
	tipnames = W.keys()
	tip2donortip = tip2donortip.items()
	#~ print 'loaded data; %i tipnames' % len(tipnames), time.ctime()
	for u in tipnames:
		for v in tipnames:
			if u==v:
				continue
			#
			p = P[u][v] # prob u transm v
			i = 0
			if (v,u) in tip2donortip:
				i = 1
			#
			pdonor_isdonor.append((p,i))
#

if len(sys.argv) > 3:
	xlabel('Time', fontsize = 16)
	ylabel('Number infected', fontsize = 16)
	savefig('figures/%s_prevalence.png' % ostem, dpi = 300)
	savefig('figures/%s_prevalence.eps' % ostem, dpi = 300)
	savefig('figures/%s_prevalence.tiff' % ostem, dpi = 300)
	#~ show()

print 'pdonor_isdonor', time.ctime()

pdis = array(pdonor_isdonor)
pi_order = argsort(pdis[:,0])
pdis = pdis[pi_order]




#~ downsample; get rid of a lot of zeros
pdis = pdis[pdis[:,0]>1e-4,:]
#~ pdis = pdis[pdis[:,0]>5*1e-2,:]
#~ pdis = pdis[-int(.1*len(pdis)):] # top 10 pc

print 'number of pairs: %i' % len(pdis), time.ctime()
print 'number of transmissions: %i' % sum(pdis[:,1]), time.ctime()
print 'expected number of transmissions: %i' % sum(pdis[:,0]), time.ctime()
#~ assert False

from scipy.stats import ttest_1samp
print ttest_1samp( pdis[:,0] - pdis[:,1], 0.)#distribution of residuals should have mean equal 0


#~ local average
lrp = pdis[radius:-radius,0] # local regression p
lrisd = list() # local regression isdonor
for i in range(radius, len(pdis)-radius):
	lrisd.append( mean(pdis[(i-radius):(i+radius),1]) )
#

print 'local regression data', time.ctime()

# linear regression
from scipy.stats import linregress
slope, intercept, r_value, p_value, slope_std_error = linregress(pdis[:,0], pdis[:,1])
print 'slope, intercept, r_value, p_value, slope_std_error'
print slope, intercept, r_value, p_value, slope_std_error

from scipy.optimize import fmin_bfgs
rss = lambda a: sum( (pdis[:,0]*a - pdis[:,1])**2.)
fmin = fmin_bfgs( rss, 1.)
print fmin
intercept = 0.;  slope = fmin[0]


#~ make figures
rc('legend',**{'fontsize':8})

f1 = figure(1, figsize=(3.5,3.0), dpi=600)
f1.subplots_adjust(bottom=.13, top=.98)
plot(pdis[:,0], pdis[:,1], 'k|', markersize=6.)
plot(lrp, lrisd, 'r-', label='Moving avg(radius=25)')
plot( [0,1], [0,1], 'k-', label='y=x')
plot( pdis[:,0], intercept + slope * pdis[:,0], 'g--', label='Linear regression' )
ylim(-.05, 1.05)
xlim(-.05, 1.05)
#~ legend(loc='lower left', bbox_to_anchor=(0., 1.))
#~ legend(loc='center left')
legend(loc='lower left', bbox_to_anchor=(0., .65))
xlabel('Estimated transmission probability')
savefig('figures/%s_lr1.png' % ostem, dpi=600)
savefig('figures/%s_lr1.tiff' % ostem, dpi=600)
savefig('figures/%s_lr1.svg' % ostem, dpi=600)
savefig('figures/%s_lr1.eps' % ostem, dpi=600)

#~ figure(2) 
#~ # histogram
#~ x = pdis[pdis[:,1]==1.,0]
#~ hist(x, bins = 50)


#~ figure(3)
#lowess (doesn't work well for binary data)
#~ from rpy2.robjects import r as R 
#~ lowess = R['lowess'](pdis[:,0].tolist(), pdis[:,1].tolist(), f=.05)
#~ plot( array( lowess[0] ), array( lowess[1] ), 'r-')
#~ import statsmodels.api as sm
#~ pdis2 = pdis[pdis[:,0]>1e-8,:]
#~ lowess = sm.nonparametric.lowess(pdis2[:,0], pdis2[:,1], frac=.1)
#~ print 'lowess', time.ctime()
#~ plot( lowess[:,0], lowess[:,1], 'r-')
#~ plot( [0,1], [0,1], 'k-')

f4 = figure(4, figsize=(3.5,3.0), dpi=600)
f4.subplots_adjust(bottom=.13, top=.98)
#different sort of local average
x = linspace(0., 1., 20)
xradius = .05
y = [ mean( pdis[ abs(pdis[:,0]-xx)< xradius,1] ) for xx in x]
plot(x, y, 'r-', label=r'Moving avg(radius=5%)')
plot( [0,1], [0,1], 'k-', label='x=y')
plot(pdis[:,0], pdis[:,1], 'k|', markersize=6.)
plot( pdis[:,0], intercept + slope * pdis[:,0], 'g--', label='Linear regression' )
ylim(-.05, 1.05)
xlim(-.05, 1.05)
legend(loc='lower left', bbox_to_anchor=(0., .65))
xlabel('Estimated transmission probability')
savefig('figures/%s_lr2.png' % ostem, dpi=600)
savefig('figures/%s_lr2.tiff' % ostem, dpi=600)
savefig('figures/%s_lr2.svg' % ostem, dpi=600)
savefig('figures/%s_lr2.eps' % ostem, dpi=600)




#~ fancy version of 5
from mpl_toolkits.axes_grid1 import make_axes_locatable
f5 = figure(5, figsize=(3.5,5.5), dpi=600)
f5.subplots_adjust(bottom=.09, top=.98, right=.96, left=.16)
ax = subplot(111)
#different sort of local average
x = linspace(0., 1., 20)
xradius = .05
y = [ mean( pdis[ abs(pdis[:,0]-xx)< xradius,1] ) for xx in x]
plot(x, y, 'r-', label=r'Moving avg(radius=5%)')
plot( [0,1], [0,1], 'k-', label='x=y')
plot(pdis[:,0], pdis[:,1], 'k|', markersize=6.)
plot( pdis[:,0], intercept + slope * pdis[:,0], 'g--', label='Linear regression' )
ylim(-.05, 1.05)
xlim(-.05, 1.05)
legend(loc='lower left', bbox_to_anchor=(0., .65))
#~ xlabel('Estimated transmission probability')

divider = make_axes_locatable(ax)
axHistU = divider.append_axes("top", 1.2, pad=0.1, sharex=ax)
axHistB = divider.append_axes("bottom", 1.2, pad=0.1, sharex=ax)
setp(axHistU.get_xticklabels() + ax.get_xticklabels(), visible=False) #axHistB.get_yticklabels()
axHistB.set_xlabel('Estimated transmission probability')
pU = pdis[pdis[:,1]==1,0]
pB = pdis[pdis[:,1]==0,0]
axHistU.hist(pU, bins=20)
axHistB.hist(pB, bins=20)

savefig('figures/%s_lr2f.png' % ostem, dpi=600)
savefig('figures/%s_lr2f.tiff' % ostem, dpi=600)
savefig('figures/%s_lr2f.svg' % ostem, dpi=600)
savefig('figures/%s_lr2f.eps' % ostem, dpi=600)





#~ show()
print time.ctime(), 'complete'


"""
for v in tipnames:
	print v, sum( [P[u][v] for u in tipnames if u!=v] )

 [sum( [P[u][v] for u in tipnames if u!=v] ) for v in tipnames ]
for v in tipnames:
	print v, sum( [P[v][u] for u in tipnames if u!=v] )
"""
