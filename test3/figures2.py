"""
SIMNUM = 0
n2load = 50
radius = 25
slope, intercept, r_value, p_value, slope_std_error
1.03667617149 0.00572063240026 0.422754180801 7.74363664361e-108 0.0446967869781

"""

calcMeanW = False #True
SIMNUM = 0
n2load = 50
radius = 25

from pylab import *
import cPickle, csv, pdb, time, colgem2



true_nwk, sampleTimes, sampleStates, Fs, Gs, Ys, transmissions = cPickle.load(open( 'sim%i/data.cPickle' % SIMNUM, 'r') )
t_incidence = loadtxt('fgy/coalescent_taxis_mle1.csv', delimiter=',')


Ws = list()
for i in range(n2load):
	Ws.append( cPickle.load(open('sim%i/W_nwk%i.cPickle' % (SIMNUM,i), 'r')) ) #
	print time.ctime(), 'loaded W', i
Ps = Ws

keyTranslator = kt = dict.fromkeys(Ws[-1].keys())
kt_ = dict()
for k in kt.keys():
	kt[k] = '\"%s\"' % k.split('_')[0]
	kt_['\"%s\"' % k.split('_')[0]] = k

#~ average P
if calcMeanW:
	W = P = dict.fromkeys(kt.values())
	for u in P.keys():
		P[u] = dict()
	#
	for u in P.keys():
		for v in P.keys():
			#~ try:
			P[u][v] = mean([PP[kt_[u]][kt_[v]] for PP in Ps] )
			#~ except:
				#~ continue
	cPickle.dump(W, open('sim%i/Wmean.cPickle' % SIMNUM, 'w') )
else:
	W = P = cPickle.load(open('sim%i/Wmean.cPickle' % SIMNUM, 'r') )
#
print time.ctime(), 'meanW'

pdonor_isdonor = list()
tipnames = sampleTimes.keys()
tip2donortip = dict()
for u,v in transmissions:
	tip2donortip['\"%i\"' % v] = '\"%i\"' % u
for u in tipnames:
	for v in tipnames:
		if u==v:
			continue
		#
		i = 0
		#~ if (u,v) in transmissions:
		if tip2donortip[v] == u:
			i = 1
		#
		try:
			pdonor_isdonor.append((P[u][v],i))
		except:
			pdb.set_trace()
#
print 'pdonor_isdonor', time.ctime()

pdis = array(pdonor_isdonor)
pi_order = argsort(pdis[:,0])
pdis = pdis[pi_order]
#~ downsample; get rid of a lot of zeros
pdis = pdis[pdis[:,0]>1e-8,:]

print 'pdis; len %i' % len(pdis), time.ctime()


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


#~ make figures
rc('legend',**{'fontsize':8})





#~ fancy version of 5
from mpl_toolkits.axes_grid1 import make_axes_locatable
f5 = figure(5, figsize=(3.5,5.5), dpi=600)
f5.subplots_adjust(bottom=.09, top=.98, right=.96, left=.16)
ax = subplot(111)
#different sort of local average
x = linspace(0., 1., 20)
xradius = .05
#~ y = [ mean( pdis[ abs(pdis[:,0]-xx)< xradius,1] ) for xx in x]
#~ plot(x, y, 'r-', label=r'Moving avg(radius=5%)')
plot(lrp, lrisd, 'r-', label='Moving avg(radius=25)')
plot( [0,1], [0,1], 'k-', label='x=y')
plot(pdis[:,0], pdis[:,1], 'k|', markersize=6.)
plot( pdis[:,0], intercept + slope * pdis[:,0], 'g--', label='Linear regression' )
ylim(-.05, 1.05)
xlim(-.01, max(lrp) +.05)
legend(loc='lower left', bbox_to_anchor=(0., .45))
#~ xlabel('Estimated transmission probability')

divider = make_axes_locatable(ax)
axHistU = divider.append_axes("top", 1.2, pad=0.1, sharex=ax)
axHistB = divider.append_axes("bottom", 1.2, pad=0.1, sharex=ax)
setp(axHistU.get_xticklabels() + ax.get_xticklabels(), visible=False) #axHistB.get_yticklabels()
axHistB.set_xlabel('Estimated infector probability')
pU = pdis[pdis[:,1]==1,0]
pB = pdis[pdis[:,1]==0,0]
axHistU.hist(pU, bins=20)
axHistB.hist(pB, bins=20)
axHistB.set_xticks(arange(0, max(lrp) +.05, .1))
xlim(-.01, max(lrp) +.05)

savefig('figures/sim2_hiv_lr2f.png', dpi=600)
savefig('figures/sim2_hiv_lr2f.tiff', dpi=600)
savefig('figures/sim2_hiv_lr2f.svg', dpi=600)
savefig('figures/sim2_hiv_lr2f.eps', dpi=600)
