#!/usr/bin/env python
###############################################################################
#                                                                             #
#    dataparser.py                                                            #
#                                                                             #
#    Class for storing and printing scaffolds and associated summary info     #
#                                                                             #
#    Copyright (C) Alexander Baker                                            #
#                                                                             #
###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

__author__ = "Alexander Baker"
__copyright__ = "Copyright 2015"
__credits__ = ["Alexander Baker"]
__license__ = "GPLv3"
__maintainer__ = "Alexander Baker"
__email__ = "Alexander.baker@uqconnect.edu.au"

###############################################################################
###############################################################################
###############################################################################
import dataloader
import scaffold
###############################################################################
'''
Class which does all the precalculations needed for giving input to
both the gapestimator and the scaffolding algorithm


'''
class dataparser(object):
    def __init__(self,
                 links,
                 coverages,
                 inserts,
                 bamnames,
                 contignames):
        ###DataParser - takes data from DataLoader instance
        ###Process - scrub completely unrealistic links (Huge insert for library)
        ###Get dictionary of contig and possibly links
        ###Make an initial graph - connections for any no. reads >threshold
        ###Take graph - Isolated contigs - write them out as Scaffold
        ###If tigs not isolated - look at paths
        ###Use decision process to keep certain paths (sspace like acceptance)
        ###Some sorta visualisation tool
        ###Hoping for a sparsely connected graph
        self.links=links
        self.coverages=coverages
        self.inserts=inserts
        self.bamnames=bamnames
        self.contignames=contignames
        self.scaffoldset={}

    def readCounts(self,contig1,contig2,linksfile=None):
        if linksfile==None:
            linksfile=self.cleanedlinks
        else:
            pass
        norm=len([x for x in linksfile if (contig1 in x) and (contig2 in x) and (x[4]=='1' and x[7]=='0')])
        rev=len([x for x in linksfile if (contig1 in x) and (contig2 in x) and (x[4]=='0' and x[7]=='1')])
        posinvc2=len([x for x in linksfile if (contig1 in x) and (contig2 in x) and (x[4]=='1' and x[7]=='1')])
        posinvc1=len([x for x in linksfile if (contig1 in x) and (contig2 in x) and (x[4]=='0' and x[7]=='0')])
        total=norm+rev+posinvc2+posinvc1
        return [norm,rev,posinvc2,posinvc1,total]
        
    def getlinks(self,contig1,contig2=False,linksfile=None):
        if linksfile==None:
            linksfile=self.links
        try:
            if contig2==False:
                Links=[x for x in linksfile if contig1 in x]
            else:
                Links=[row for row in linksfile if (contig1 in row) and (contig2 in row)]
            if len(Links)>0:
                return Links
            else:
                return None
        except TypeError:
            print "Contig1 must be string, contig2 can be ither False or a string"
            
    def checklinks(self,contig1,linksfile=None):
        if linksfile==None:
            linksfile=self.cleanedlinks
        links=self.getlinks(contig1,False,linksfile)
        notcontig1=set([x[0] for x in links if x[0]!=contig1])
        nottig1=set([x[1] for x in links if x[1]!=contig1])
        return list(nottig1|notcontig1)

    def completecheck(self,cleanedlinks=None):
        if cleanedlinks==None:
            cleanedlinks=self.cleanedlinks
        connections={}
        for contig in self.contignames:
            connections[contig]= self.checklinks(contig,cleanedlinks)
        return connections

    def makegraph(self,connections=None,cleanedlinks=None,threshold=5):
        ''' Intent is to get the linksfile which has been scrubbed of reads
        with far to large insert, the dictionary of contigs and all tigs with at least one read.
        This will be used to construct a graph, each contigs will be a node, it will join other contigs
        if it passed the threshold number
        '''
        if cleanedlinks==None:
            cleanedlinks=self.cleanedlinks
        if connections==None:
            connections=self.completecheck()
        OrientedGraph={}
        Graph={}
        for contig1 in connections:
            Graph[contig1]=[]
            OrientedGraph[contig1]={}
            for contig2 in connections[contig1]:
                Counts=self.readCounts(contig1,contig2)
                index=[i for i,x in enumerate(Counts) if (x==max(Counts)) and (x>=threshold)]
                OrientedGraph[contig1][contig2]=[self.linkorientation(i) for i in index if self.linkorientation(i)!=None]
                Graph[contig1]=Graph[contig1]+[contig2]
        return [OrientedGraph,Graph]

    def makescaffolds(self,Graph,OrientedGraph):
        Scaffolds={}
        for contig1 in Graph:
            for contig2 in Graph[contig1]:
                Scaffolds["Scaffold"+str(len(Scaffolds))]={}
        return None
    def arrangenplace(self,linkdirs):
        ''' Return tuple indicating how a pair of linked
        contigs should be arranged based on their read orientations'''
        flipplace=False
        flipsequence=False
        or1,or2=linkdirs
        if or1==1 and or2==0:
            flipplace=True
            flipsequence=False
        elif or1==0 and or2==1:
            flipplace=False
            flipsequence=False
        elif or1==1 and or2==1:
            flipplace=None
            flipsequence=True
        elif or1==0 and or2==0:
            flipplace=None
            flipsequence=True
        else:
            pass
        return (flipplace,flisequence)
                
    def linkorientation(self,index):
        try:
            if index==0:
                return (1,0)
            elif index==1:
                return (0,1)
            elif index==2:
                return (1,1)
            elif index==3:
                return (0,0)
            else:
                pass
        except:
            print "The index needs to be integer"
            
    def cleanlinks(self,linksfile=None,insertfile=None,cutoff=3):
        ''' Trims of the incredibly unrealisticly gapped reads
         before use downstream - should save computation and prevent
        erroneous computation'''
        if linksfile==None:
            linksfile=self.links
        else:
            pass
        if insertfile==None:
            insertfile=self.inserts
        bamNames=self.bamnames
        meaninsert={}
        stdinsert={}
        bamnameindex=8
        for name in bamNames:
            #Should give unique insert mean and insert std
            meaninsert[name]=[x[1] for x in insertfile if name+".bam" in x]
            stdinsert[name]=[x[2] for x in insertfile if name+".bam" in x]
            ####Need to mod this list compre to get b
        print meaninsert
        print stdinsert
        cleanedlinks=[x for x in linksfile if self.lowerdist(x)<\
        meaninsert[x[bamnameindex].strip('.bam')]+cutoff*stdinsert[x[bamnameindex].strip('.bam')]]
        self.cleanedlinks=cleanedlinks
        self.mean=meaninsert
        self.std=stdinsert
        return
        
    def parse(self):
        '''links all functions to start from the DataParser input and to end
        with the processed scaffolds(with gap estimate and proper orientation).
        It will also call all of the scaffolds to print to the scaffold file.'''
        self.cleanlinks() #Cleans the obviously erroneous mappings
        
        
        return
    def lowerdist(self,onelink):
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

    def makevertices(self,graph):
        vertice=[]
        for key in graph:
            for point in graph[key]:
                vertice.append((key,point))
        return vertice

    def lonetigs(self,graph):
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
        return



    def algorithm(stuff):
    # calls all the graph functions to go from start to finish in constructing the scaffolds
    #Will return the dictionary of scaffolds suitable for scaffold class - or
    #Do final step in here as well - make dictionary of Scaffolds objects
        return

    def gapest(self,contig1=False,contig2=False):
    #implementation of sahlini et al 2013 algorithm.
        import os
        import sys
        from scipy.special import erf
        from scipy.stats import norm
        from scipy.constants import pi
        from math import exp
        
        if contig2==False:
            return 0
        else:
            #Stuff happens
            pass
        return 200

    def extractpath(self,graph):
        '''Scaffolds are equivalent to paths through the
        connected contig graph'''
        return path
