def readCounts(contig1,contig2,linksfile):
	norm=len([x for x in linksfile if (contig1 in x) and (contig2 in x) and (x[4]=='1' and x[7]=='0')])
	rev=len([x for x in linksfile if (contig1 in x) and (contig2 in x) and (x[4]=='0' and x[7]=='1')])
	posinvc2=len([x for x in linksfile if (contig1 in x) and (contig2 in x) and (x[4]=='1' and x[7]=='1')])
	posinvc1=len([x for x in linksfile if (contig1 in x) and (contig2 in x) and (x[4]=='0' and x[7]=='0')])
	total=norm+rev+posinvc2+posinvc1
	return [norm,rev,posinvc2,posinvc1,total]

def checklinks(contig1,linksfile):
	links=self.getlinks(contig1,linksfile)
	notcontig1=set([x[0] for x in linksfile if x[0]!=contig1])
	nottig1=set([x[1] for x in linksfile if x[1]!=contig1])
	return list(nottig1|notcontig1)

def completecheck(self,cleanedlinks):
	connections={}
	for contig in self.contigNames:
		connections={contig: checklinks(contig,cleanedlinks)}
	return connections

def cleanlinks(linksfile,cutoff=4):
''' Trims of the incredibly unrealisticly gapped reads
 before use downstream - should save computation and prevent
erroneous computation'''
	bamNames=self.bamNames
	meaninsert={}
	stdinsert={}
	bamnameindex=8
	for name in bamNames:
		#Should give unique insert mean and insert std
		meaninstert[name]=[x for x in coveragefile if name in x][5]
		stdinstert[name]=[x for x in coveragefile if name in x][6]
		
		####Need to mod this list compre to get b
	cleanedlinks=[x for x in linksfile if lowerdist(x)<\
	meaninsert[x[bamnameindex].strip('.bam')]+cutoff*stdinsert[x[banameindex].strip('.bam')]]
	return cleanedlinks

def bamfreecleanlinks(linksfile,cutoff=4):
	''' Trims of the incredibly unrealisticly gapped reads
	before use downstream - should save computation and prevent
	erroneous computation'''
	try:
		linksfile.readline()
		linksfile[0]
	except TypeError:
		notlist=True
	if notlist:
		cleanedlinks=[x.strip().split('\t') for x in linksfile if lowerdist(x)<300+cutoff*30]
	else:
		cleanedlinks=[x for x in linksfile[1:] if lowerdist(x)<300+cutoff*30]
	return cleanedlinks


def lowerdist(onelink):
###severe lower bound on gap between reads (ignores gap between contigs)
###Probably need to implement gap estimation 
###Without gap estimate should be sever understimate
#Need to check bamm if pos is end of read or middle - if middle need to get length
	orientation1=int(onelink[4])
	orientation2=int(onelink[7])
	if orientation1:
		dist1=int(onelink[3])
	else:
		dist1=int(onelink[2])-int(onelink[3])
	if orientation2:
		dist2=int(onelink[6])
	else:
		dist2=int(onelink[5])-int(onelink[6])
	return (dist1+dist2)
		
	
def makegraphs(checklist,cleanedlinks,threshhold):
''' Intent is to get the linksfile which has been scrubbed of reads
with far to large insert, the dictionary of contigs and all tigs with at least one read.
This will be used to construct a graph, each contigs will be a node, it will join other contigs
if it passed the threshold number
'''

	

def makevertices(graph):
	vertice=[]
	for key in graph:
		for point in graph[key]:
			vertice.append((key,point))
	return vertice

def lonetigs(graph):
	singletigs=[]
	for node in graph:
		if graph[node]==[]:
			singletigs.append(node)
		else:
			pass
	return singletigs


def collapse(graph,cutoff):
''' Using the sspace-like decision process to simply the graph.
This result in the remove of certain paths where contigs had more than 1
pairing
'''



def algorithm(stuff):
# calls all the graph functions to go from start to finish in constructing the scaffolds


def gapest(stuff):
#implementation of sahlini et al 2013 algorithm.
import os
import sys
from scipy.special import erf
from scipy.stats import norm
from scipy.constants import pi
from math import exp



def retroply(cleanedgraph):
#Apply gapest to all linked members of the graph

def extractpath(graph):
	



	return path



