"""
[[ 1.          0.83206661]
 [ 0.83206661  1.        ]]
"""

should_calculate_P_truenwk = False#True
SIMNUM = 0
n2load = 50

from pylab import *
import cPickle, csv, pdb, time, colgem2, re, copy

print 'start', time.ctime()
#~ format of data: cPickle.dump( (nwk, oSampleTimes, oSampleStates, Fs, Gs, Ys, t_incidence, transmissions),  open( 'transmProbTest/simulationsHIV.cPickle', 'w') )
true_nwk, sampleTimes, sampleStates, Fs, Gs, Ys, transmissions = cPickle.load(open( 'sim%i/data.cPickle' % SIMNUM, 'r') )
print 'loaded simdata', time.ctime()
t_incidence = loadtxt('fgy/coalescent_taxis_mle1.csv', delimiter=',')


# get W for true tree
if should_calculate_P_truenwk:
	replace_string = lambda inmatch: '%s' % inmatch.group().split('_')[0]
	true_nwk = re.sub('([0-9]+_[0-9]+)',replace_string, true_nwk)
	
	fgy = colgem2.FGY(t_incidence, Fs, Gs, Ys)
	ggs= colgem2.GeneGenealogyStates(true_nwk, sampleTimes, sampleStates, fgy, repair_negative_branch_lengths=True)
	W = P = colgem2.transmission_probabilities(ggs) 
	cPickle.dump(W, open('sim%i/trueW.cPickle' % SIMNUM, 'w'))
else:
	W = P = cPickle.load(open('sim%i/trueW.cPickle' % SIMNUM, 'r'))
#

print time.ctime(), 'true W'



# load W from beast posterior samples
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





PvsP = list()
for u in P.keys():
	for v in P.keys():
		if u==v or P[u][v] < 1e-8:
			continue
		#~ pdb.set_trace()
		PvsP += [(P[u][v], PP[kt_[u]][kt_[v]]) for PP in Ps]
PvsP = array(PvsP)

#~ make the figure
f = figure(figsize=(3.5, 3.1), dpi = 600)
f.subplots_adjust(bottom = .18, left = .18, top = .95, right=.95)
plot(PvsP[:,0], PvsP[:,1], 'ko', markersize=.5)
plot(PvsP[:,0], PvsP[:,0], 'r-')
xlabel('Infector probability for true tree', fontsize=10)
ylabel('Infector probability for estimated trees', fontsize=10)

savefig('figures/hiv_PvsP.png', dpi=600)
savefig('figures/hiv_PvsP.tiff', dpi=600)
savefig('figures/hiv_PvsP.svg', dpi=600)
savefig('figures/hiv_PvsP.eps', dpi=600)


#~ corrcoef
print corrcoef(PvsP[:,0],  PvsP[:,1])


