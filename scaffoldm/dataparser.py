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
class DataParser(object):
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
        '''Gets counts for each orientation pair of links'''
        #Might be able to force proper orientation here.
        if linksfile==None:
            linksfile=self.cleanedlinks
        else:
            pass
        #Changes from needing to maintain relative arrangement
        tig1ind=linksfile[0].find(contig1)
        tig2ind=linksfile[0].find(contig2)
        #Now maintains appropiate order, if refering to contigs in seperate arrangements
        #This is since index determines if x[4] refers to read of tig 1 or tig 2.
        if tig1ind<tig2ind:
            norm=len([x for x in linksfile if (contig1 in x) and (contig2 in x) and (x[4]=='1' and x[7]=='0')])
            rev=len([x for x in linksfile if (contig1 in x) and (contig2 in x) and (x[4]=='0' and x[7]=='1')])
        elif tig2ind<tig1ind:
            norm=len([x for x in linksfile if (contig1 in x) and (contig2 in x) and (x[4]=='0' and x[7]=='1')])
            rev=len([x for x in linksfile if (contig1 in x) and (contig2 in x) and (x[4]=='1' and x[7]=='0')])
        else:
            print "How is this possible?"
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
                if index==[]:
                    pass
                elif len(index)>0:
                    for i in index:
                        if self.linkorientation(i)!=None:
                            OrientedGraph[contig1][contig2]=[self.linkorientation(i)]
                    Graph[contig1]=Graph[contig1]+[contig2]
        print Graph
        print OrientedGraph
        return [OrientedGraph,Graph]

    def makescaffolds(self,Graph,OrientedGraph):
        '''Takes Graph of all connections and one with orientations,
        makes a dictionary of Scaffold objects for scaffoldset'''
        Scaffolds={}
        isolated=self.lonetigs(Graph)
        for tig in isolated:
            #Add all contigs with no links as individual scaffolds
            #Also remove from both Graph and OrientedGraph
            Scaffolds["Scaffold"+str(len(Scaffolds))]={tig:[0,0,0,0,0]}
            Graph.pop(tig,None)
            OrientedGraph.pop(tig,None)
        print Scaffolds
        #Now need to retrieve each scaffold with the order,gapsize,orientation
        paths=makepaths(Graph,OrientedGraph)
        return None
        
    def makepaths(self,Graph,OrientedGraph):
    '''Made with assumption that isolated tigs have been removed'''
        edges=self.makeedges(Graph)
        badedges={}
        for edge in edges:
            #Check orientation, expunge all illegal arrangements
            swap,flip=arrangenplace(OrientedGraph[edge[0]][edge[1]])
            if swap and not flip:
                badedges[(edge[0],edge[1])]=[0]
            else:
                pass
        trueedges=[x for x in edges if x not in badedges.keys()]
        connections={}
        for edge in trueedges:
            path[]
            connections[edge]=[(1,x) if edge[1]==x[0] else (0,x) if edge[0]==x[1] for x in trueedges]
        for key in connections:
            paths=[]
            for edges in connections[key]:
                if edge[0]==0:
                    list(key
                elif edge[1]==1:
                else:
                    pass
    #Is there some way to memoize this
    
    def findpath(Graph,verystart=None,start=None,end=None,path=[]):
        if start==None: #initialise specific path
            pass
        elif start not in Graph:
            return []
        else:
            path=path+[start]
        if start==end and start!=None:
            return [path]
        paths=[]
        if start==None:
            for key in Graph:
                for edge in Graph[key]:
                    if edge not in path:
                        print verystart,edge,"HELLLLLLLO"
                        extrapaths=findpath(Graph,key,key,edge,path)
                        for pathz in extrapaths:
                            paths.append(pathz)
        else:
            for edge in Graph[start]:
                if edge not in path:
                    print verystart,start,edge
                    extrapaths=findpath(Graph,verystart,end,edge,path)
                    for pathz in extrapaths:
                        paths.append(pathz)
        return paths

    def rearrange(self,lists,item1,item2):
        '''swap two items in a list'''
        a, b = i.index(item1), i.index(item2)
        lists[b], lists[a] = lists[a], lists[b]
        return lists
        
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
        '''Index of counts gives an idea of read orientation'''
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
        self.makegraph() #Makes initial connections with no.reads>threshold
        #self.makescaffolds() #Makes set of scaffolds
        #self.printscaffolds()
    def printscaffold(self,filename='testScaffold'):
        try:
            for key in self.scaffoldset:
                scaffoldset[key].printscaffold('testScaffold')
        except AttributeError:
            print 'These scaffolds do not appear to be the scaffold class'
        
        
        return
    def lowerdist(self,onelink):
        '''severe lower bound on gap between
         reads (ignores gap between contigs)'''
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

    def makeedges(self,graph):
        '''edges as (node1,node2) tuples.'''
        edges=[]
        for key in graph:
            for point in graph[key]:
                if (point,key) not in edges and (key,point):
                    edges.append((key,point))
        return edges
        
    def makepaths(self,graph):
        '''Scaffolds are equivalent to paths along the vertices of the graph.
        Here, scaffoldM uses the vertices to construct paths and returns them.
        The legality of paths is not checked'''
        vertices=makeedges(graph)
        
    def lonetigs(self,graph):
        '''Identifies all contigs which are not linked
        to any other contig.'''
        singletigs=[]
        for node in graph:
            if graph[node]==[]:
                singletigs.append(node)
            else:
                pass
        return singletigs


    def collapse(self,graph,cutoff):
        ''' Using the sspace-like decision process to simply the graph.
        This result in the remove of certain paths where contigs had more than 1
        pairing
        '''
        return None



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
