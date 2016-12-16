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
	#print "Using full hSIN"
	for tabline in PBDs_file:
		(uniprot, pfam) = tabline.strip().split("\t")
		PBDs.add( (uniprot, pfam) )


with open("./datafiles/DNABPs_Entrez.txt", "r") as DNABPs_file:
#with open("../datafiles/DNAcc_TFs.txt", "r") as DNABPs_file:
	for tabline in DNABPs_file:
		#(etz, name) = tabline.strip().split("\t")
		etz = tabline.strip()
		DNABPs.add( int(etz) )




Proteins_Load( DNABPs )
#print "Loaded %s Proteins" % ( len(Protein.byEntrez) )



TotalMuts =		0

MutPosRE = re.compile(r"\d+")
with open("../datafiles/ESP6500_v2.txt", 'r') as esp_f:
	for tabline in esp_f:
		columns = tabline.strip().split("\t")
		
		try:
			entrez = int(columns[3])
			mut = int( MutPosRE.search(columns[8]).group() )
			freq = float(columns[5])
		except:
			continue
		muttype = columns[9]
		
		# Common: > 1%		Rare: < 0.1%
		if muttype != "missense" or freq <= 0.01:
			continue
		
		try:
			prot = Protein.byEntrez[entrez]
		except KeyError:
			continue
		
		if mut > prot.Length:
			continue
		
		prot.Mutations.append(mut)
		


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
		if DNABres[i] == 1 and PBres[i] == 1:
			ResInBoth += 1
	ResInDNABDs += sum(DNABres)
	ResInPBDs += sum(PBres)
	ResInDoms += sum(Domres)
	TotalRes += prot.Length
	
	DNABMuts = 0
	PBMuts = 0
	BothMuts = 0
	DomMuts = 0
	
	for mut in prot.Mutations:
		TotalMuts += 1
		if DNABres[mut-1] == 1:
			DNABMuts += 1
		if PBres[mut-1] == 1:
			PBMuts += 1
		if (DNABres[mut-1] + PBres[mut-1]) == 2:
			BothMuts += 1
		if Domres[mut-1] == 1:
			DomMuts += 1
	MutsInDNABDs += DNABMuts
	MutsInPBDs += PBMuts
	MutsInBoth += BothMuts
	MutsInDoms += DomMuts
	if prot.Entrez in hORF:
		out_str += "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (prot.Entrez, prot.Name, BothMuts, DNABMuts, PBMuts, DomMuts, len(prot.Mutations))
	
with open("./Vaq_ESP_Screen.txt", "w") as out_file:
	out_file.write( out_str )

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
