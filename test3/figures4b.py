"""208 cherries
auc est trees 0.847117565764
auc true tree 0.848208815399
"""

SIMNUM = 0

from pylab import *
import cPickle, csv, pdb, time, colgem2, re, copy
from ete2 import *

print 'start', time.ctime()
true_nwk, sampleTimes, sampleStates, Fs, Gs, Ys, transmissions = cPickle.load(open( 'sim%i/data.cPickle' % SIMNUM, 'r') )
print 'loaded simdata', time.ctime()
t_incidence = loadtxt('fgy/coalescent_taxis_mle1.csv', delimiter=',')

#~ find the set of cherries in the true tree
tree = Tree(true_nwk) 
cherries =  [u.get_leaf_names() for u in tree.get_descendants() if len(u.get_leaves())==2] 
getname = lambda taxon: '\"%s\"' % re.search('([0-9]+)_', taxon).groups()[0]
cherries = [ (getname(u), getname(v))  for u,v in cherries  ]
########################################################################
# ROC fig
spid = 'hiv_sim_beast'

# true W
tW = tP = cPickle.load(open('sim%i/trueW.cPickle' % SIMNUM, 'r'))

#~ average P
W = P = cPickle.load(open('sim%i/Wmean.cPickle' % SIMNUM, 'r') )
print 'averaged W', time.ctime()


pdonor_isdonor = list()
pdonor_isdonor_tp = list()
tipnames = sampleTimes.keys()
tip2donortip = dict()
#~ pdb.set_trace()
for u,v in transmissions:
	tip2donortip['\"%i\"' % v] = '\"%i\"' % u

for u,v in cherries:
	i = 0
	#~ if (u,v) in transmissions:
	if tip2donortip[v] == u:
		i = 1
	#
	pdonor_isdonor.append( (P[u][v], i) )
	pdonor_isdonor_tp.append( (tP[u][v], i) )
#
print 'pdonor_isdonor', time.ctime()


#####################

pdis = array(pdonor_isdonor)
pi_order = argsort(pdis[:,0])
pdis = pdis[pi_order]
#~ downsample; get rid of a lot of zeros
#~ pdis = pdis[pdis[:,0]>1e-8,:]
#~ minpdis = min(pdis[pdis[:,0]>0,0])
#~ pdis = pdis[pdis[:,0]>=minpdis,:]

pdistp = array(pdonor_isdonor_tp)
pi_order = argsort(pdistp[:,0])
pdistp = pdistp[pi_order]
#~ downsample; get rid of a lot of zeros
#~ pdistp = pdistp[pdistp[:,0]>1e-8,:]
#~ pdistp = pdistp[pdistp[:,0]>minpdis,:]

print 'pdis; len %i' % len(pdis), time.ctime()


##
fn_rates = list()
tp_rates = list()
fp_rates = list()
all_pos = float ( sum(pdis[:,1]) )
all_neg = float( len(pdis) - all_pos ) 
#~ qrange = linspace(1e-8, 1., 100)
qrange = linspace(0, 1., 100)
for q in qrange:
	#~ posmax = max( (pdis[:,0]<=q) * arange(len(pdis)) ) 
	posmax = sum(pdis[:,0]<q) 
	tp = sum(pdis[posmax:,1])
	#~ tp_rate = tp / float(len(pdis) - posmax)
	tp_rate = tp / all_pos
	
	#~ fp = sum(pdis[posmax:,1]==0)
	#~ fp_rate = float(len(pdis) - posmax)
	fn_rate = sum(pdis[:posmax,1]) / all_pos
	
	fp_rate = sum(pdis[posmax:,1]==0) / all_neg
	
	#~ pdb.set_trace()
	
	fn_rates.append(fn_rate)
	tp_rates.append(tp_rate)
	fp_rates.append( fp_rate)
#

fn_ratestp = list()
tp_ratestp = list()
fp_ratestp = list()
all_postp = float ( sum(pdistp[:,1]) )
all_negtp = float( len(pdistp) - all_postp ) 
for q in qrange:
	#~ posmax = max( (pdistp[:,0]<=q) * arange(len(pdistp)) ) 
	posmax = sum(pdis[:,0]<q) 
	tp = sum(pdistp[posmax:,1])
	#~ tp_rate = tp / float(len(pdis) - posmax)
	tp_rate = tp / all_postp
	
	#~ fp = sum(pdis[posmax:,1]==0)
	#~ fp_rate = float(len(pdis) - posmax)
	fn_rate = sum(pdistp[:posmax,1]) / all_postp
	
	fp_rate = sum(pdistp[posmax:,1]==0) / all_negtp #?
	
	fn_ratestp.append(fn_rate)
	tp_ratestp.append(tp_rate)
	fp_ratestp.append( fp_rate)
#
##

#~ figure()
#~ plot(qrange, c_[fn_rates, tp_rates], 'r-')
#~ plot([0,1], [0,1], 'k-')


f = figure(figsize=(3.5,3.0), dpi=600) 
ax = subplot(111)
f.subplots_adjust(bottom=.13, top=.96, right=.96, left=.18)
plot( fp_ratestp, tp_ratestp, 'k-', label='True tree')
plot( fp_rates, tp_rates, 'k--', label='Est. trees')
plot([0,1], [0,1], 'r-')
xlabel('False positive rate')
ylabel('True positive rate')
legend(loc='best')

savefig('figures/%s_roc.png' % spid, dpi=600)
savefig('figures/%s_roc.tiff' % spid, dpi=600)
savefig('figures/%s_roc.svg' % spid, dpi=600)
savefig('figures/%s_roc.eps' % spid, dpi=600)


#~ show()

#~ calculate auc
from scipy.interpolate import interp1d
i = interp1d(fp_rates[::-1], tp_rates[::-1], kind='linear', fill_value=0., bounds_error=False)
auc = mean( i(linspace(0., 1., 1000)[1:] ) ) 
print 'auc est trees', auc

i = interp1d(fp_ratestp[::-1], tp_ratestp[::-1], kind='linear', fill_value=0., bounds_error=False)
auc = mean( i(linspace(0., 1., 1000)[1:] ) ) 
print 'auc true tree', auc
