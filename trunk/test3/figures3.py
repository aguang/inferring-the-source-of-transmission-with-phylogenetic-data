from pylab import *
import cPickle
import pdb, time

#~ def figure3():
#~ make a nice figure of prevalence over time

true_nwk, sampleTimes, sampleStates, Fs, Gs, Ys, transmissions = cPickle.load(open( 'sim0/data.cPickle', 'r') )
t_incidence = 1975. + loadtxt('fgy/coalescent_taxis_mle1.csv', delimiter=',') / 365.
Y = array(Ys)
acute = Y[:,0] + Y[:,5]
aids = Y[:,4] + Y[:,-1]
chronic = sum(Y[:,1:4], axis=1) + sum( Y[:,1:4], axis = 1)
diagnosed = sum(Y[:,5:], axis = 1)
undiagnosed = sum(Y[:,:5], axis = 1)
xticks([])

f = figure(figsize = (3.5, 5.8), dpi = 600)
f.subplots_adjust(left = .14, bottom = .1, top = .97, right = .97)
topaxes = f.add_subplot(211)
#~ subplot(211)
plot(t_incidence, acute, label='EHI')
plot(t_incidence, chronic, label='Chronic')
plot(t_incidence, aids, label='AIDS')
legend(loc='best')
topaxes.xaxis.set_ticklabels([''])
#~ topaxes.xaxis.set_visible(False)
subplot(212)
plot(t_incidence, undiagnosed, label='Undiag.')
plot(t_incidence, diagnosed, label='Diag.')
legend(loc='lower left')
xlabel('Year')
xticks([1975, 1985, 1995, 2005])

savefig('figures/simulationsHIV-prevalence.png', dpi = 600)
savefig('figures/simulationsHIV-prevalence.tiff', dpi = 600)
savefig('figures/simulationsHIV-prevalence.svg', dpi = 600)
savefig('figures/simulationsHIV-prevalence.eps', dpi = 600)
