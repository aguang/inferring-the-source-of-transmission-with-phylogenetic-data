""" aug 21 
 add flags for 2x2 analysis (sample time x sample fraction)
"""
#~ from modelGenerator import *
from colgem2 import *
from simulationsHeterochronous3x2 import *
import os, sys
from time import ctime
from scipy.integrate import odeint

phiflag, stflag = sys.argv[1:]

STIMES = [16.65, 150.] 

phi = float( phiflag)/100. 
tfin = STIMES[ int(stflag ) ]
stem = 'output/%s_phi%s_st%s_' % (repr(os.getpid()), phiflag, stflag)

#~ parameters
m = 2 # constant
target_n = 150*75 #1000
#~ ntrees = 2

stLimit = .999 #gives essentially homochronous sample
finiteSizeCorrection = True


#~ generate the model
class DeterministicSIRS: 
	N = 1e+04#5e+3
	gamma = 1.
	beta = 1.5
	mu = .1
	t = linspace(0., tfin, 100)
	def _F(self, x):
		s, i, r = x
		return array( [[self.beta * i * s / self.N, 0.],[0., 0.]] )
	def _G(self, x):
		s, i, r = x
		return array( [[0., self.gamma*i],[0., 0.]] )
	#
	def xdot(self, x, t):
		s, i, r = x
		F = sum( self._F(x) )
		G = sum( self._G(x) )
		ds = -F + self.mu * r
		di = F - G
		dr = G - self.mu * r
		return ds, di, dr
	def __init__(self,**kwargs):
		self.__dict__.update(**kwargs)
		x0 = (self.N-1., 1., 0.) # one infection
		self.x = odeint(self.xdot, x0, self.t)
		self.s, self.i, self.r = self.x.T
		self.Y = self.x[:,1:3]
		h = self.t[1]
		self.F = [ self._F(x) * h for t,x in zip( self.t[1:], self.x )]
		self.G = [ self._G(x) * h for t,x in zip( self.t[1:], self.x )]
#

print ctime(), 'begin'
dsirs = DeterministicSIRS()

if len(sys.argv)>3:
	# make a plot
	assert False

#~ plot(dsirs.t, dsirs.Y); show()

#~ generate the trees, transm events, transm probabilities
#~ pdb.set_trace()
fgy = FGY(dsirs.t, dsirs.F, dsirs.G, dsirs.Y)
ftdata = list() #ft sims
#~ for i in range(ntrees):
nsampled = 0
while nsampled < target_n:
	deaths = zeros((len(dsirs.t), 2)) #not applicable
	nwk, sampleTimes, sampleStates, o_births, o_migrations, o_prevalence,tip2donortip  = ft_simulate(dsirs.F, dsirs.G, dsirs.Y, deaths, dsirs.t, phi=phi, stLimit = stLimit)
	
	if sum(sampleStates.values(), axis = 0)[0] != len(sampleStates):
		#check only samples I state
		pdb.set_trace()
	
	print ctime(), 'simulation complete'
	
	nsampled += len(sampleTimes)
	
	o_births = [ array( o_births[i:i+m]) for i in range(m, len(o_births), m) ]
	o_migrations = [ array( o_migrations[i:i+m]) for i in range(m, len(o_migrations), m) ]
	o_prevalence = array(o_prevalence, dtype=float)
	
	#~ pdb.set_trace() #check dimensions
	
	#~ transm events
	fgy2  = FGY(dsirs.t, o_births, o_migrations, o_prevalence)
	ggs = GeneGenealogyStates(nwk, sampleTimes, sampleStates, fgy2, finiteSizeCorrection=finiteSizeCorrection, repair_negative_branch_lengths = True)
	#~ Cs, Ps, P, rterminals = transmission_probabilities(ggs)
	#~ tipnames = [rt().name for rt in rterminals]
	#~ ftdata.append((nwk, sampleTimes, sampleStates, tip2donortip,  Cs, Ps, P, tipnames, o_prevalence))
	
	print ctime(), 'calculated GGS'
	#~ pdb.set_trace()
	
	W= transmission_probabilities(ggs, tol=1e-12)
	ftdata.append((nwk, sampleTimes, sampleStates, tip2donortip,  W,  o_prevalence))
	
#




#~ output
import cPickle
cPickle.dump( ftdata, open(stem + 'ft.cPickle', 'w'))
print ctime(), 'done'

#~ print [P[recip][donor] for recip, donor in tip2donortip.items() if P.has_key(recip) and P.has_key(donor)]
