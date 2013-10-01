#!/bin/bash
#~ http://bioweb2.pasteur.fr/docs/seq-gen/
#~ meanRate 4.446774e-06, ucld.mean 4.4245e-6
#~ april 3 2012, I think this parameter is actually half of what beast reports... -t9.468

date

for i in {0..19}
do
	STEM=sim$i
	seq-gen -mHKY -t4.734 \
	 -f0.416,0.196,0.189,0.199 \
	 -i0.472 \
	 -a0.714 \
	 -g4 \
	 -s4.446774e-06 \
	 -l1200 \
	 -n1 \
	 -on \
	 -q \
	 $STEM/tree.nwk > $STEM/msa.nexus
	date
done

#~ use -op for interleaved format (phylip)
