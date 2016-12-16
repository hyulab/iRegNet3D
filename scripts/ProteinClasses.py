#!/usr/bin/python
# February 2014
# Nate Tippens

import re;

hORF = set()
Vidal = set()


class Protein:
	byUniprot = dict()
	byEntrez = dict()

	def __init__(self):
		self.Uniprot = ""
		self.Entrez = -1
		self.Name = ""
		self.Length = None
		self.Domains = list()
		self.Mutations = list()
		self.GOs = set()
		self.PPIs = set()
		self.SCAPairs = list()



class Domain:
	def __init__(self, pfam, start, stop):
		self.Pfam = pfam
		self.Start = int(start)
		self.Stop = int(stop)
	
	def contains(self, aapos):
		return (aapos >= self.Start and aapos <= self.Stop)
	
	def __len__(self):
		return self.Stop - self.Start



class Pfam:
	get = dict()

	def __init__(self, pfam):
		if pfam in Pfam.get:
			raise ValueError("%s already exists!" % pfam)
		Pfam.get[pfam] = self
		self.Description = ""
		self.GOs = set()



class GOterm:
	get = dict()

	def __init__(self, GO):
		if GO in GOterm.get:
			raise ValueError("%s already exists!" % GO)
		self.GO = GO
		GOterm.get[GO] = self
		self.Parents = list()
		self.Children = list()








def Proteins_Load(idlist = None):
	subset = set()
	try:
		for pid in idlist:
			subset.add( int(pid) )
	except:
		subset = idlist
		
	with open("/home/nate/Rmain/resources/entrez/entrez_human2longest_uniprot.txt", 'r') as convers_f:
		for tabline in convers_f:
			(etz, uniprot) = tabline.strip().split("\t")
			etz = int(etz)
			if subset == None or etz in subset or uniprot in subset:
				prot = Protein()
				prot.Uniprot = uniprot
				prot.Entrez = etz
				Protein.byUniprot[uniprot] = prot
				Protein.byEntrez[etz] = prot


	with open("/home/nate/Rmain/resources/uniprot/9606_SPROT_domains.txt", 'r') as pfam_f:
		for tabline in pfam_f:
			(uniprot, pfam, start, stop) = tabline.strip().split("\t")
			if uniprot not in Protein.byUniprot:
				continue
			
			prot = Protein.byUniprot[uniprot]
			try:
				dom = Domain(pfam, start, stop)
			except ValueError:
				continue
			prot.Domains.append(dom)
			if pfam not in Pfam.get:
				Pfam(pfam)


	with open("/home/nate/Rmain/resources/uniprot/9606_SPROT_domain_info.txt", 'r') as dominfo_file:
		for tabline in dominfo_file:
			(pfam, id, name, type) = tabline.strip().split("\t")
			if pfam in Pfam.get:
				Pfam.get[pfam].Description = name


	with open("/home/nate/Rmain/resources/uniprot/uniprot_human2length.txt", 'r') as length_f:
		for tabline in length_f:
			(uniprot, length) = tabline.strip().split("\t")
			if uniprot in Protein.byUniprot:
				Protein.byUniprot[uniprot].Length = int(length)

	with open("/home/nate/Rmain/resources/entrez/Homo_sapiens.gene_info", 'r') as ginfo_f:
		for tabline in ginfo_f:
			cols = tabline.strip().split("\t")
			etz = int(cols[1])
			if etz in Protein.byEntrez:
				Protein.byEntrez[etz].Name = cols[2]



def hORF_Load():
	with open("/home/nate/Rmain/resources/yulab/hORFemone_V8.1_Entry_DNASU.txt", 'r') as horf_f:
		for tabline in horf_f:
			cols = tabline.strip().split("\t")
			try:
				etz = int(cols[5])
				hORF.add(etz)
			except:
				continue



def Y2HVidal_Load(unique = False):
	with open("/home/nate/Rmain/resources/yulab/HumanY2HVidal.txt", 'r') as yth_f:
		for tabline in yth_f:
			cols = tabline.strip().split("\t")
			protA = int(cols[0])
			protB = int(cols[1])
			if protA in Protein.byEntrez:
				Protein.byEntrez[protA].PPIs.add(protB)
			if protB in Protein.byEntrez:
				Protein.byEntrez[protB].PPIs.add(protA)
			if unique:
				Vidal.add(protA)
				Vidal.add(protB)




def Gene_GO_Load(ontology):
	ontname = None
	if ontology in ("molecular_function", "molecular", "mf", "m"):
		ontname = "Function"
	elif ontology in ("biological_process", "process", "bp", "b", "p"):
		ontname = "Process"
	elif ontology in ("cellular_component", "cell", "component", "cc", "c"):
		ontname = "Component"
	else:
		ontname = "Function"

	with open("/home/nate/Rmain/resources/entrez/Human_entrez2go.txt", 'r') as UniGO_file:
		for tabline in UniGO_file:
			cols = tabline.strip().split(",")
			try:
				etz = int(cols[1])
				GOid = str(int(cols[2][3:]))
			except:
				continue
			
			evidence = cols[3]
			category = cols[7]
			
			if category != ontname:
				continue
				
			# see http://www.geneontology.org/GO.evidence.shtml
			if evidence not in ("EXP", "IDA", "IPI", "IC", "IBA"):
				continue

			if etz in Protein.byEntrez:
				Protein.byEntrez[etz].GOs.add( GOid )



def Pfam_GO_Load():
	with open("/home/nate/Rmain/resources/GO/pfam2go.txt", 'r') as PFGO_file:
		for tabline in PFGO_file:
			(pfam, name, description, GOid) = tabline.strip().split("\t")
			GOu = GOterm[GOid]
			GOu.Description = description
			GOu.GO = GOid
			Pfam.get[pfam].GOs.add(GOu)



def GO_Heirarchy_Load(ontology):
	ontname = None
	if ontology in ("molecular_function", "molecular", "mf", "m", "f"):
		ontname = "molecular_function"
	elif ontology in ("biological_process", "process", "bp", "b", "p"):
		ontname = "biological_process"
	elif ontology in ("cellular_component", "cell", "component", "cc", "c"):
		ontname = "cellular_component"
	else:
		ontname = "molecular_function"
	
	with open("/home/nate/Rmain/resources/GO/GO_%s_1.2.txt" % ontname, 'r') as GO_f:
		for tabline in GO_f:
			(cstr, pstr) = tabline.strip().split("\t")
			child = str(int(cstr[4:]))
			parent = str(int(pstr[4:]))
			if child not in GOterm.get:
				GOterm(child)
			if parent not in GOterm.get:
				GOterm(parent)
			GOterm.get[child].Parents.append( GOterm.get[parent] )
			GOterm.get[parent].Children.append( GOterm.get[child] )

