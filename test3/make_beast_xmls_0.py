""" april 27 '13
take fasta, gen xml from stem

replace these at top of file: 
<?xml version="1.0" standalone="yes"?>

<beast>
	<taxa id="taxa">
		<taxon id="21910_8835">
			<date value="8835.0" direction="forwards" units="days"/>
		</taxon>
	</taxa>
	
	
		<alignment id="alignment" dataType="nucleotide">
		<sequence>
			<taxon idref="17362_10066"/>
			ACTG
		</sequence>
	</alignment>
	
	
	edit this: 
	<log id="fileLog" logEvery="50000" fileName="stem.log" overwrite="false">
	
	edit this:
	<logTree id="treeFileLog" logEvery="50000" nexusFormat="true" fileName="stem.trees" sortTranslationTable="true">
"""
import pdb, os, re
from Bio import SeqIO
STEMXML = 'stem.xml'

stemstring = open( STEMXML, 'r').read()
for i in range(20):
	DIR = 'sim%i' % i
	
	top = '''<?xml version="1.0" standalone="yes"?>
<beast>
	<taxa id="taxa">'''
	
	date_from_sid = lambda sid: sid.split('_')[-1]
	for sr in SeqIO.parse( '%s/msa.nexus' % DIR , 'nexus'):
		top += '''
			<taxon id="%s">
				<date value="%s" direction="forwards" units="days"/>
			</taxon>''' % (sr.id, date_from_sid(sr.id))
	#
	
	top += '''
	</taxa>

	<alignment id="alignment" dataType="nucleotide">
	'''
	
	nseq = 0
	for sr in SeqIO.parse( '%s/msa.nexus' % DIR , 'nexus'):
		nseq += 1
		top += '''
			<sequence>
				<taxon idref="%s"/>
				%s
			</sequence>
			''' % (sr.id, sr.seq.tostring() )
	top += '''
	</alignment>
	'''
	
	
	bottom = stemstring
	bottom = re.sub( 'stem\.log', DIR + '.log', bottom )
	bottom = re.sub( 'stem\.trees', DIR + '.trees', bottom)
	#~ if re.match('skyride', bottom):
	if 'skyride' in bottom:
		bottom = re.sub('102', repr(nseq-1), bottom) #102 is the skyride dimension in stem.xml
		bottom = re.sub('204', 2*repr(nseq-1), bottom) #102 is the skyride dimension in stem.xml
	
	open( '%s/%s.xml' %(DIR,DIR), 'w').write(top + bottom)
#
