#!/usr/bin/python
# December 2013
# Nate Tippens

import sys, os, re
from ProteinClasses import *

DNABPs = set() # DNA Binding Proteins
DNABDs = set() # DNA Binding Domains
PBDs = set() # Protein Binding Domains


with open("./datafiles/Expanded_DNABDs.txt", 'r') as DNABDs_file:
	for textline in DNABDs_file:
		DNABDs.add( textline.strip() )


with open("../datafiles/hSIN_ProtBindDomains.txt", 'r') as PBDs_file:
#with open("../datafiles/hSIN_cc_ProtBindDomains.txt", 'r') as PBDs_file:
	for tabline in PBDs_file:
		(uniprot, pfam) = tabline.strip().split("\t")
		PBDs.add( (uniprot, pfam) )


with open("./datafiles/DNABPs_Entrez.txt", "r") as DNABPs_file:
#with open("../datafiles/DNAcc_TFs.txt", 'r') as DNABPs_file:
	for tabline in DNABPs_file:
		#(etz, name) = tabline.strip().split("\t")
		etz = tabline.strip()
		DNABPs.add( int(etz) )




Proteins_Load( DNABPs )
#print "Loaded %s Proteins" % ( len(Protein.byEntrez) )



TotalMuts =		0

# NOTE: You must subscribe to HGMD and re-format its mutation datafile to run this code.
# Specifically, you must map all mutations onto Uniprot IDs and AA positions.
# Michael Meyer has published a webtool (Bisque) to perform these mappings.
MutPosRE = re.compile(r"c.\d+")
with open("../datafiles/HGMD2013.4_indels_Uniprots.txt", 'r') as hgmd_f:
	for tabline in hgmd_f:
		(uniprot, mut, mut_type, disease, inframe) = tabline.strip().split("\t")
		
		if mut_type != "DM" or mut == "" or inframe == "FALSE":
			continue
		
		try:
			prot = Protein.byUniprot[uniprot]
		except KeyError:
			continue
		
		try:
			mut_pos = int(mut)
		except ValueError:
			mut_pos = int( MutPosRE.search(mut).group()[2:] )
			mut = str(mut_pos)
		if mut_pos <= prot.Length:# and mut_pos not in prot.Mutations:
			prot.Mutations.append( mut_pos )


# compute non-overlapping domain coverage
ResInDNABDs =	0
ResInPBDs =		0
ResInBoth =		0
ResInDoms =		0
TotalRes =		0

# compute mutation type
MutsInDNABDs =	0
MutsInPBDs =	0
MutsInBoth =	0
MutsInDoms =	0

out_str = "Entrez\tSymbol\tDoubleMuts\tDNABMuts\tPBMuts\tDomMuts\tAllMuts\n"
for prot in Protein.byUniprot.values():
	if len(prot.Mutations) == 0:
		continue
	
	DNABres = [0] * prot.Length
	PBres = [0] * prot.Length
	Domres = [0] * prot.Length
	for dom in prot.Domains:
		Domres[dom.Start : dom.Stop+1] = [1]*(dom.Stop-dom.Start+1)
		if dom.Pfam in DNABDs:
			DNABres[dom.Start : dom.Stop+1] = [1]*(dom.Stop-dom.Start+1)
		if (prot.Uniprot, dom.Pfam) in PBDs:
			PBres[dom.Start : dom.Stop+1] = [1]*(dom.Stop-dom.Start+1)
	
	if(sum(DNABres) == 0):
		continue
	
	for i in range(prot.Length):
		if DNABres[i] + PBres[i] == 2:
			ResInBoth += 1
	ResInDNABDs += sum(DNABres)
	ResInPBDs += sum(PBres)
	ResInDoms += sum(Domres)
	TotalRes += prot.Length
	if sum(Domres) > prot.Length:
		print "%s\t%s\t%s" % (prot.Uniprot, sum(Domres), prot.Length)
	
	for mut in prot.Mutations:
		TotalMuts += 1
		if DNABres[mut-1] == 1:
			MutsInDNABDs += 1
		if PBres[mut-1] == 1:
			MutsInPBDs += 1
		if (DNABres[mut-1] + PBres[mut-1]) == 2:
			MutsInBoth += 1
		if Domres[mut-1] == 1:
			MutsInDoms += 1


print "AMINO ACIDS"
print "DNABD:\t%s" % ResInDNABDs
print "PBD:\t%s" % ResInPBDs
print "Double:\t%s" % ResInBoth
print "Domain:\t%s" % ResInDoms
print "Total:\t%s" % TotalRes
print "\n"

print "MUTATIONS"
print "DNABD:\t%s" % MutsInDNABDs
print "PBD:\t%s" % MutsInPBDs
print "Double:\t%s" % MutsInBoth
print "Domain:\t%s" % MutsInDoms
print "Total:\t%s" % TotalMuts


