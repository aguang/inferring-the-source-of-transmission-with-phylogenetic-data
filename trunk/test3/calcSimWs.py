"""

python calcSimWs.py 0 0
python calcSimWs.py 0 1 2 3 4 5 6 7 8 9 10 11 12 &
python calcSimWs.py 0 13 14 15 16 17 18 19 20 21 22 23 24 &
python calcSimWs.py 0 25 26 27 28 29 30 31 32 33 34 35 36 &
python calcSimWs.py 0 36 37 38 39 40 41 42 43 44 45 46 47 48 49 &
"""

from pylab import *
import pdb, cPickle, time, colgem2, sys

whichsim = sys.argv[1]
whichnwks = [int(arg) for arg in sys.argv[2:]]

nwks=[nwk for nwk in open('sim%s.nwk' % whichsim, 'r').read().split() if len(nwk) > 0]
nwks=[nwk for i,nwk in enumerate(nwks) if i in whichnwks]

simdata = cPickle.load(open('sim%s/data.cPickle' % whichsim, 'r'))
true_nwk, sampleTimes_, sampleStates_, Fs, Gs, Ys, transmissions = simdata
sampleTimes = dict([( taxon.strip('\"')+'_%i' % int(st), st ) for taxon, st in sampleTimes_.items()])
sampleStates = dict([(taxon.strip('\"')+'_%i' % int(sampleTimes_[taxon]), ss) for taxon, ss in sampleStates_.items()])
t_incidence = loadtxt('fgy/coalescent_taxis_mle1.csv', delimiter=',')


fgy = colgem2.FGY(t_incidence, Fs, Gs, Ys)
for i,nwk in zip(whichnwks,nwks):
	print 'start P-calc', time.ctime()
	try:
		ggs= colgem2.GeneGenealogyStates(nwk, sampleTimes, sampleStates, fgy, repair_negative_branch_lengths=True)
		W = colgem2.transmission_probabilities(ggs) 
		cPickle.dump(W, open('sim%s/W_nwk%i.cPickle' % (whichsim, i), 'w'))
	except colgem2.NegativeBranchLengthException:
		
		continue
	print 'completed a P-calc', time.ctime()
#


