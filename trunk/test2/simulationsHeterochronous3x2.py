"""  12 12 12 OO version
transm indicator 0..1 
dict keyed by tip_%i 

#TODO neversamplestates 1
"""

from pylab import *
from likelihoodC import *
import networkx as nx
from scipy.stats.distributions import rv_discrete
import pdb, csv, time
from copy import deepcopy



def ft_simulate(births, migrations, prevalence, deaths, taxis, phi = .20, stateToSample=0, stLimit = .95, forgiveNoneSampleable = False):
	n =  int(phi * prevalence[-1][stateToSample]) # new in 3x2: sample only from one state
	print 'n', n
	m = len(prevalence[0])
	#~ random sampling
	sampleTimes = linspace(max(taxis)*stLimit, max(taxis), n)
	#~ sampleStates = [eye(m)[randint(0,m)] for i in range(n) ]
	sampleStates = [eye(m)[stateToSample] for i in range(n) ] #only sample one state
	
	sampleTimes_pop = copy(sampleTimes).tolist() # temporary object to keep track of which have been reached
	sampleStates_pop = copy(sampleStates).tolist()
	
	
	state_labels = [ 'I%i' % i for i in range(m/2) ] + ['D%i' % i for i in range(m/2) ]
	
	
	c = 0 # counter to keep track of names of new dudes
	cc =0 # counter for sample
	prevalence0 = [int(round(p)) for p in prevalence[0]]
	sid2state = dict()
	sid2branchState = dict() # for ancestor functions
	sids_set = set() # list of infected
	sampled = set() # those who are sampled
	state_sids = dict([(i, set()) for i in range(m)])
	transmission_graph = nx.DiGraph()
	transmission_genealogy = nx.DiGraph() # includes nodes that correspond to transmission events
	sid2node = dict() # will map individual infected to a node in the genealogy
	genealogy_node_time= dict()
	for i,p in enumerate(prevalence0):
		for pp in range(p):
			sid2state[c] = i
			sid2branchState[c] = [(0., i)]
			sids_set.add(c)
			state_sids[i].add(c)
			sid2node[c] = repr(c)
			genealogy_node_time[repr(c)] = 0.
			transmission_genealogy.add_node(sid2node[c])
			c+=1
			
	#
	
	
	
	
	tips = list()
	tip_states = list()
	tip_times = list()
	tip2state = dict()
	sample=list()
	sample_states = list()
	tip2sid = dict() # donor is sid
	sid2tip = dict()
	
	
	#~ t = 0
	#~ interval =0
	interval = 1
	t = taxis[interval-1] 
	t_next_interval = taxis[interval]
	t_delta = taxis[interval]- taxis[interval-1]
	
	#
	state_timeSeries = o_prevalence = [[len(state_sids[i]) for i in range(m)]]
	o_migrations = zeros((m,m)).tolist(); o_migration = zeros((m,m))
	o_births = zeros((m,m)).tolist(); o_birth = zeros((m,m))
	o_deaths = list(); 
	
	rounded_births = 0.#TODO
	
	done = False
	while not done:
		try:
			birth = round_(births[interval-1]).astype(int)
		except:
			pdb.set_trace()
		migration = round_(migrations[interval-1]).astype(int)
		death = round_(deaths[interval]).astype(int)
		
		summigrations = sum(migrations[interval-1])
		sumbirths = sum(births[interval-1])
		sumdeaths = sum(deaths[interval])
		
		db =  sum(births[interval-1]) - sum(birth) 
		#~ pdb.set_trace()
		if rand() < db: #(db - int(db)):
			db = ceil(db)
		else:
			db = floor(db);
		if db >= 1:
			whichbirth = rv_discrete(name='whichbirthrvs', values = [arange(m**2), births[interval - 1]/sumbirths]).rvs(size=int(db))
			for wb in whichbirth:
				fromstate = int(wb) /  m
				tostate = mod(int(wb), m)
				birth[fromstate,tostate] += 1.
				#~ print 'print fromstate, tostate',  fromstate, tostate
			#
		#
		
		db = ceil( sum(migrations[interval - 1]) - sum(migration) )
		if rand() < db: #(db - int(db)):
			db = ceil(db)
		else:
			db = floor(db);
		if db >= 1:
			whichmig = rv_discrete(name='whichrvs', values = [arange(m**2), migrations[interval - 1]/summigrations]).rvs(size=int(db))
			for wb in whichmig:
				fromstate = int(wb) /  m
				tostate = mod(int(wb), m)
				migration[fromstate,tostate] += 1.
			#
		#
		
		db = ceil( sum(deaths[interval]) - sum(death) )
		if rand() < db:
			db = ceil(db)
		if db >= 1:
			whichdeath = rv_discrete(name='whichd', values = [arange(m), deaths[interval] / sumdeaths]).rvs(size=int(db))
			for wb in whichdeath:
				death[wb] += 1.
		#
		
		rounded_births += sum(birth) #TODO 
		
		events = list()
		for k in range(m):
			for l in range(m):
				events += [('b', k, l)]*birth[k,l]
				events += [('m', k, l)]*migration[k,l]
			#
			events += [('d', k, k)]*death[k]
		#
		shuffle(events)
		tevents = sort(  uniform(t, t_next_interval, size=len(events)) )
		for tt, event in zip(tevents, events):
			if len(sampleTimes_pop)==0:
				done=True
				continue
			#
			while sampleTimes_pop[0] < tt: 
				ttt = sampleTimes_pop.pop(0)
				state = argmax( sampleStates_pop.pop(0) )
				
				# simple random sample
				#~ sid = list(sids_set)[randint(0,len(sids_set))] # sampled guy
				
				
				#~ sample specific state
				sampleable = state_sids[state].difference(sampled) 
				#~ pdb.set_trace()
				if len(sampleable) > 0:
					sid = list(sampleable)[randint(0,len(sampleable) ) ]  
				elif not forgiveNoneSampleable:
					print 'could not sample state', state
					continue
				elif forgiveNoneSampleable:
					sampleable = sids_set.difference(sampled)
					sid = list(sampleable)[randint(0,len(sampleable))] # sampled guy
					print 'could not sample state', state
				#
				sampled.add(sid)
				
				#~ make leaf nodes
				new_genealogy_tip = 'tip_' + repr(cc) #'(%s,%s)' % (sid2node[sid], 'tip_' + repr(cc)) #
				new_genealogy_node = '(%s,%s)' % (sid2node[sid], 'sample_' + repr(cc))
				tips.append(new_genealogy_tip)
				tip_states.append( eye(m)[sid2state[sid]].tolist() )
				tip_times.append( ttt )
				cc+=1
				
				sid2tip[sid] = deepcopy( new_genealogy_tip )
				tip2sid[new_genealogy_tip] = sid # sid2tip #new_genealogy_node#deepcopy( transmission_graph.predecessors(sid)[0] )
				
				transmission_genealogy.add_edge(sid2node[sid], new_genealogy_node)
				transmission_genealogy.add_edge(new_genealogy_node, new_genealogy_tip)
				genealogy_node_time[new_genealogy_tip] = ttt
				genealogy_node_time[new_genealogy_node] = ttt
				transmission_genealogy[new_genealogy_node][new_genealogy_tip]['interval'] = (ttt, ttt)
				transmission_genealogy[sid2node[sid]][new_genealogy_node]['interval'] = (genealogy_node_time[sid2node[sid]], ttt)
				transmission_genealogy[new_genealogy_node][new_genealogy_tip]['state'] = sid2branchState[sid] #[(ttt, sid2state[sid])]
				transmission_genealogy[sid2node[sid]][new_genealogy_node]['branch_length'] = ttt - genealogy_node_time[sid2node[sid]]
				
				sid2node[sid] = new_genealogy_node #new_genealogy_tip #does this break it? 
				
				tip2state[new_genealogy_tip] = tip_states[-1]
				
				#~ make sample, record states
				sample.append(new_genealogy_tip)
				sample_states.append(tip_states[-1])
				#~ continue
				
				if len(sampleTimes_pop)==0:
					done = True
					break
			#			
			if event[0]=='b':
				#~ birthhappens
				fromstate = event[1]
				tostate = event[2]
				cwhile = 0; 
				donorfound = False
				try:
					donor = list(state_sids[fromstate])[randint(0,len(state_sids[fromstate]))]
				except:
					#~ pdb.set_trace()
					while len(state_sids[fromstate]) ==0:
						cwhile+=1
						if cwhile > 100:
							fromstate = argmax([len(state_sids[k]) for k in range(m)])
							print 'birth fail', event, tt, [len(state_sids[k]) for k in range(m)]
							#~ pdb.set_trace()
							donor=sid2node.keys()[0]
							break
						rbirths = birth.ravel()
						sumbirths2 = sum(birth)
						whichbirth = rv_discrete(name='whichbirthrvs', values = [arange(m**2), rbirths/sumbirths2]).rvs()
						fromstate = int(whichbirth) /  m
						tostate = mod(int(whichbirth), m)
					#
				recipient = c
				c+=1
				state_sids[tostate].add(recipient)
				sid2state[recipient] = tostate
				sids_set.add(recipient)
				transmission_graph.add_edge(donor, recipient)
				
				o_birth[fromstate,tostate] += 1.
				
				new_genealogy_node = '(' + sid2node[donor] + ',' + repr(recipient) + ')'
				transmission_genealogy.add_edge(sid2node[donor], new_genealogy_node)
				
				genealogy_node_time[new_genealogy_node] = tt
				transmission_genealogy[sid2node[donor]][new_genealogy_node]['branch_length'] = tt - genealogy_node_time[sid2node[donor]]
				transmission_genealogy[sid2node[donor]][new_genealogy_node]['interval'] = (genealogy_node_time[sid2node[donor]], tt)
				transmission_genealogy[sid2node[donor]][new_genealogy_node]['state'] = sid2branchState[donor] #[ (tt, tostate) ]
				
				sid2node[donor]= new_genealogy_node
				sid2node[recipient]= new_genealogy_node
				sid2branchState[recipient] = [(t, sid2state[recipient])]
			elif event[0]=='m':
				#~ migration happens
				fromstate = event[1]
				tostate = event[2]
				if len(state_sids[fromstate]) > 1:#0: #3x2, protect state 0 from being empty
					donor = list(state_sids[fromstate])[randint(0,len(state_sids[fromstate]))]
					state_sids[fromstate].remove(donor)
					state_sids[tostate].add(donor)
					sid2state[donor] = tostate
					sid2branchState[donor].append( (t, tostate) )
					o_migration[fromstate,tostate]+=1.
			elif event[0]=='d':
				#~ death happens
				whichdeath = event[1]
				if len(state_sids[whichdeath]) > 0:
					deadguy = list(state_sids[whichdeath])[randint(0,len(state_sids[whichdeath]))]
					state_sids[whichdeath].remove(deadguy)
					sids_set.remove(deadguy)
			#
		#
		
		interval += 1
		t = taxis[interval-1]
		
		# for validation 
		state_timeSeries.append([len(state_sids[i]) for i in range(m)])
		o_migrations.extend(o_migration)
		o_births.extend(o_birth) 
		o_birth = zeros((m,m))
		o_migration = zeros((m,m))
		#~ print state_timeSeries[-1]
		#~ print t, time.ctime()
		
		#~ for ancestor functions
		for sid in sids_set:
			"""v = sid2node[sid]
			try:
				u = transmission_genealogy.predecessors(v)[0]
			except:
				#~ no parent
				continue
			try:
				transmission_genealogy[u][v]['state'].append( (t, sid2state[sid]) )
			except:
				transmission_genealogy[u][v]['state'] =  [(t, sid2state[sid]) ]
			"""
			sid2branchState[sid].append( (t, sid2state[sid]) )
		
		try:
			t_next_interval = taxis[interval ]
		except IndexError:
			done = True
			t = taxis[-1]
	#
	#
	print time.ctime(), 'make genealogy'
	#~ make sample genealogy
	include = set(copy(sample))
	rtg = transmission_genealogy.reverse()
	for ego in sample:
		include = include.union( nx.dfs_tree(rtg, source=ego).nodes() )
	#
	transmission_genealogy_sample = phylo = transmission_genealogy.subgraph(include)
	
	#~ TODO
	rnodes = [node for node in phylo.nodes() if phylo.out_degree(node) >=3 ]
	
	
	roots = list()
	for node in copy(transmission_genealogy_sample.nodes()):
		if phylo.out_degree(node)==1 and phylo.in_degree(node)==1: # do root nodes separately
			# merge branches
			parent = phylo.predecessors(node)[0]
			child = phylo.successors(node)[0]
			phylo.add_edge(parent, child)
			#~ phylo[parent][child]['branch length'] = phylo[node][child]['branch length'] + phylo[parent][node]['branch length']
			phylo.remove_node(node)
		elif phylo.in_degree(node)==0:
			roots.append(node)
	# 
	roots = set(roots)
	print roots
	done = False
	while not done:
		done = True
		for root in roots:
			if phylo.out_degree(root)<2:
				roots.remove(root)
				roots.add(phylo.successors(root)[0])
				phylo.remove_node(root)
				done = False
			#
	#
	print time.ctime(), 'made roots, make branch lengths'
	#~ make branch lengths
	for parent,child in transmission_genealogy_sample.edges():
		transmission_genealogy_sample[parent][child]['branch length'] = genealogy_node_time[child] - genealogy_node_time[parent]
	
	
	
	
	#~ make the transmission indicator #TODO
	#~ tips -> sids -> donors (sids) -> tips (if avail)
	tipsids = [tip2sid[tip] for tip in tips]
	tipsids_donors = [transmission_graph.predecessors(sid)[0] for sid in tipsids]
	tip2donortip = dict([(sid2tip[tipsid], sid2tip[donortipsid])  for tipsid, donortipsid in zip(tipsids, tipsids_donors) if sid2tip.has_key(donortipsid)  ])
	
	#~ donortips = [sid2tip[donor] for donor in [tip2donor[tn] for tn in tips] if sid2tip.has_key(donor)]
	#~ tipdonors = [tip2donor[tn] for tn in tips] 
	#~ tip2donortip = dict([ (tip, sid2tip[tip2donor[tip]]) for tip in tips if tip2donor[tip] in tipdonors and sid2tip.has_key(tip2donor[tip])] )
	
	strays = [tip for tip in tips if tip not in phylo.nodes()] #TODO 
	
	
	#~ make the newick
	print time.ctime(), 'make newick'
	node2newick = dict.fromkeys(phylo.nodes(), None)
	done = False
	while not done:
		done = True
		for node in phylo.nodes():
			if node2newick[node]==None:
				done = False
				children = phylo.successors(node)
				if len(children)==1 or len(children)>2:
					pdb.set_trace() # shouldn't happen
				elif len(children)==2  and node2newick[children[0]]!=None and node2newick[children[1]]!=None:
					if node in roots: # is root
						node2newick[node] = '(' + node2newick[children[0]] + ',' + node2newick[children[1]] + ');'
					else:
						parent = phylo.predecessors(node)[0]
						bl = phylo[parent][node]['branch length']
						node2newick[node] = '(' + node2newick[children[0]] + ',' + node2newick[children[1]] + '):' + repr(bl)
				elif len(children)==0: # tip
					parent = phylo.predecessors(node)[0]
					bl = phylo[parent][node]['branch length']
					node2newick[node] = node + ':' + repr(bl)
	#
	nwk = ''
	for root in roots:
		nwk+= node2newick[root] + '\n'
		nwk = nwk[:-1]
	#
	
	#~ pdb.set_trace()
	return nwk, dict(zip(tips, tip_times)), dict(zip(tips, tip_states)), o_births, o_migrations, o_prevalence, tip2donortip
#
