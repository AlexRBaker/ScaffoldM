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
import scaffold
import numpy as np
import pandas as pd
import datetime
###############################################################################
#Notes for possible conversion to numpy
class DataParser(object):
    '''
    Class which does all the precalculations needed for giving input to
    both the gapestimator and the scaffolding algorithm
    '''
    def __init__(self,
                 links,
                 coverages,
                 inserts,
                 bamnames,
                 contignames,
                 contigloc,
                 tiglen,
                 readsize=100,
                 Verbose=True):
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
        self.bamnames=[bam.split(".bam")[0] for bam in bamnames]
        self.contignames=contignames
        self.scaffoldset={}
        self.contigloc=contigloc
        self.tiglen=tiglen
        self.readlen=readsize
        self.gaps={}
    def tigin(self,contig1,linksfile):
        '''Uses numpy vectorised functions to check if a contig
        is present in a linksfile'''
        return np.logical_or(linksfile[:,0]==contig1,linksfile[:,1]==contig1)
        
    def readCounts(self,contig1,contig2,linksfile=None,clean=False):
        '''Gets the number of links between two contigs
        in each of 4 orientations, (0,1),(1,0),(1,1),(0,0)
        Here 0 is a read in the same direction as the contig'''
        linksfile=self.getlinks(contig1,contig2,clean=clean)
        if isinstance(linksfile,type(None)):
            return [0,0,0,0]
        #Changes from needing to maintain relative arrangement
        #tig1ind=linksfile[0].index(contig1)
        #tig2ind=linksfile[0].index(contig2)
        tig1ind=np.where(linksfile[0]==contig1)[0][0]
        tig2ind=np.where(linksfile[0]==contig2)[0][0]
        #Now maintains appropiate order, if refering to contigs in seperate arrangements
        #This is since index determines if x[4] refers to read of tig 1 or tig 2.
        norm=len(linksfile[np.logical_and(linksfile[:,4]=='1',linksfile[:,7]=='0'),0])
        if tig1ind<tig2ind: #Checks how the contig1 versus contig2 are place in readsfile
            norm=len(linksfile[np.logical_and(linksfile[:,4]=='1',linksfile[:,7]=='0'),0])
            #norm=len([x for x in linksfile if (contig1 in x) and (contig2 in x) and (x[4]=='1' and x[7]=='0')])
            rev=len(linksfile[np.logical_and(linksfile[:,4]=='0',linksfile[:,7]=='1'),0])
            #rev=len([x for x in linksfile if (contig1 in x) and (contig2 in x) and (x[4]=='0' and x[7]=='1')])
        elif tig2ind<tig1ind:
            norm=len(linksfile[np.logical_and(linksfile[:,4]=='0',linksfile[:,7]=='1'),0])
            #norm=len([x for x in linksfile if (contig1 in x) and (contig2 in x) and (x[4]=='0' and x[7]=='1')])
            rev=len(linksfile[np.logical_and(linksfile[:,4]=='1',linksfile[:,7]=='0'),0])
            #rev=len([x for x in linksfile if (contig1 in x) and (contig2 in x) and (x[4]=='1' and x[7]=='0')])
        else:
            print tig1ind
            print tig2ind
            print contig1
            print contig2
            norm=0
            rev=0
            print "How is this possible? contig1!=contig2"
        posinvc2=len(linksfile[np.logical_and(linksfile[:,4]=='1',linksfile[:,7]=='1'),0])
        posinvc1=len(linksfile[np.logical_and(linksfile[:,4]=='0',linksfile[:,7]=='0'),0])
        #posinvc2=len([x for x in linksfile if (contig1 in x) and (contig2 in x) and (x[4]=='1' and x[7]=='1')])
        #posinvc1=len([x for x in linksfile if (contig1 in x) and (contig2 in x) and (x[4]=='0' and x[7]=='0')])
        return np.array([norm,rev,posinvc2,posinvc1],dtype=int)
    def getlinks(self,contig1,contig2=False,linksfile=None,clean=False,bam=False):
        '''Retrieves every link entry for one contig for every link 
        between two contigs'''
        if isinstance(linksfile,type(None)):
            if clean:
                linksfile=self.cleanedlinks
            else:
                linksfile=self.links
        try:
            if contig2==False:
                if bam==False:
                    Links=linksfile[np.logical_or(linksfile[:,0]==contig1,linksfile[:,1]==contig1),:]
                else:
                    Links=linksfile[np.logical_and(np.logical_and(np.logical_or(linksfile[:,0]==contig1,linksfile[:,1]==contig1)\
                    ),linksfile[:,8]==bam),:]
            else:
                if bam==False:
                    Links=linksfile[np.logical_and(np.logical_or(linksfile[:,0]==contig1,linksfile[:,1]==contig1),\
                    np.logical_or(linksfile[:,0]==contig2,linksfile[:,1]==contig2)),:]
                else:
                    #Links=[row for row in linksfile if (contig1 in row) and (contig2 in row) and (bam in row)]
                    Links=linksfile[np.logical_and(np.logical_and(np.logical_or(linksfile[:,0]==contig1,linksfile[:,1]==contig1),\
                    np.logical_or(linksfile[:,0]==contig2,linksfile[:,1]==contig2)),linksfile[:,8]==bam),:]
                    #check if both contigs and desired bamfile
            if len(Links)>0:
                return Links
            else:
                return None
        except TypeError:
            print "Contig1 must be string, contig2 can be ither False or a string"
            
    #~ def checklinks(self,contig1,linksfile=None,clean=True):
        #~ '''Checks all contigs linked to the specified contig
        #~ and return a list of these contigs.
        #~ '''
        #~ links=self.getlinks(contig1,False,linksfile,clean)
        #~ if isinstance(links,type(None)):
            #~ return None
        #~ return list(np.unique1d(pd.unique(links[links[:,1]!=contig1,1]),pd.unique(links[links[:,0]!=contig1,0])))

    def completecheck(self,cleanedlinks=None,clean=True):
        '''Uses numpy goodness. Sorts by each column of contigs. THen, find indices of first
        instance of the unique and splits the links into subarrays by each column. Then it is
        simply a matter of looping over the subarrays to create the dictionary. About 400 times faster than old method'''
        ###Appears to be drastically faster
        connections={}
        if clean and isinstance(cleanedlinks,type(None)):
            links=self.cleanedlinks
        elif not clean and isinstance(cleanedlinks,type(None)):
            links=self.links
        else:
            links=cleanedlinks
        tigs1=links[links[:,0].argsort(),0:2] #indices to order first column
        tigs2=links[links[:,1].argsort(),0:2] #Indices to order second column
        splittigs1=np.vsplit(tigs1,np.unique(tigs1[:,0],return_index=True)[1]) #array split by first col
        splittigs2=np.vsplit(tigs2,np.unique(tigs2[:,1],return_index=True)[1]) #Array split by second col
        for array in splittigs1: #Should be ~38000 in ratdata
            if np.size(array)>0:
                key=pd.unique(array[:,0]) #Should be one contig
                if len(key)==1:
                    key=key[0] #Should now be the string
                else:
                    print "Dunno what happened here"
                if key not in connections:
                    connections[key]=list(pd.unique(array[:,1])) #Get all unique linked tigs
                elif key in connections:
                    connections[key]=list(np.union1d(connections[key],pd.unique(array[:,1])))
        for array in splittigs2:
            if np.size(array)>0:
                key=pd.unique(array[:,1])
                if len(key)==1:
                    key=key[0] #Should now be the string
                else:
                    print "Dunno what happened here",key
                if key not in connections:
                    connections[key]=list(pd.unique(array[:,0])) #Get all unique linked tigs
                elif key in connections:
                    connections[key]=list(np.union1d(connections[key],pd.unique(array[:,0])))
        return connections
        
    def graphtosif(self,graph,graphname,removed=True,dirty=False,final=False):
        #print graphname, "This is the supposed graph being parsed"
        #print "\n",graph
        header="Node1\tRelationship\tNode2\tRemoved\n"
        done=set([])
        with self.tryopen("{0}{1}".format(graphname,"_links"),"Contig1\tRelationship\tContig2\tSupportingLinks\tRemoved\n",".txt",True) as network2:
            for contig1,connected in graph.iteritems():
                for contig2 in connected:
                    if (contig1,contig2) not in done and (contig2,contig1) not in done:
                        if dirty:
                            Counts=self.readCounts(contig1,contig2,clean=False) #The 4 orientations
                            for i,links in enumerate(Counts):
                                if links!=0:
                                    network2.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(contig1,\
                                    self.linkorientation(i),contig2,links,removed[contig1][contig2]))
                        elif not dirty and not final:
                            Counts=self.readCounts(contig1,contig2,clean=True)
                            index=[i for i,x in enumerate(Counts) if (x==max(Counts))]
                            if len(index)>0:                               
                                network2.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(contig1,\
                                self.linkorientation(index[0]),contig2,max(Counts),removed[contig1][contig2]))
                        if final:
                            Counts=self.readCounts(contig1,contig2,clean=True)
                            network2.write("{0}\t{1}\t{2}\t{3}\n".format(contig1,\
                            '(0,1)',contig2,max(Counts)))
                        done|=set([(contig1,contig2)]) #Add current pair
                        done|=set([(contig1,contig2)[::-1]]) #Reverse of current pair
        if not dirty:
            with self.tryopen("{0}{1}".format("all","_nodes"),"Contig1\tContiglength\t\n",".txt",True) as nodes:
                #print self.tiglen
                for tig,length in self.tiglen.iteritems():
                    nodes.write("{0}\t{1}\t\n".format(tig,length))
        return
    def tryopen(self,filename,header,filetype,expunge=False):
        '''Looks for a file, if its not there, it makes the file, if
        it is there then it returns the opened file. Remember to close the file if you call this
        function or use: with tryopen(stuff) as morestuff.'''
        import os
        try:
            if not os.path.isfile(filename+filetype):
                with open(filename+filetype,'w') as newfile:
                    newfile.write(header)
            elif os.path.isfile(filename+filetype):
                if expunge:
                    temp=open(filename+filetype,'w+')
                    temp.write(header)
                    temp.close()
                return open(filename+filetype,'a+')
            return open(filename+filetype,'a+')
        except:
            print "Either could not create or open the file"
            raise
        
    def stagestosif(self):
        removed=self.notingraph(self.totalgraph,self.graph,"Initial","Threshold")
        self.graphtosif(self.totalgraph,"Initial",removed,dirty=True)
        removed=self.notingraph(self.graph,self.finalgraph,"Threshold","Cov_Links")
        self.graphtosif(self.graph,"Threshold",removed)
        self.graphtosif(self.finalgraph,"Cov_Links",final=True)
        return
        
    def notingraph(self,graph1,graph2,name1,name2):
        removed={}
        for key,item in graph1.iteritems():
            removed[key]={}
            for key2 in item:
                if key in graph2:
                    if key2 in graph2[key]:
                        removed[key][key2]="False"
                    else:
                        removed[key][key2]="True"
                else:
                    removed[key][key2]="True"
        return removed
                        
                            
        
    def makegraph(self,connections=None,cleanedlinks=None,threshold=5,total=False,dirty=False):
        ''' Intent is to get the linksfile which has been scrubbed of reads
        with far to large insert, the dictionary of contigs and all tigs with at least one read.
        This will be used to construct a graph, each contigs will be a node, it will join other contigs
        if it passed the threshold number
        '''
        if isinstance(connections,type(None)):
            print "Creating connection"
            connections=self.completecheck(clean=not dirty)
            print "Finished making connections"
            print datetime.datetime.now().time()
        if not dirty and isinstance(cleanedlinks,type(None)):
            links=self.cleanedlinks
        elif dirty and isinstance(cleanedlinks,type(None)):
            links=self.links
        else:
            links=cleanedlinks
        Graph={}
        tigs1=links[links[:,0].argsort(),0:2] #indices to order first column
        tigs2=links[links[:,1].argsort(),0:2] #Indices to order second column
        double=np.vstack((tigs1,np.fliplr(tigs2))) #flipped array and went left to right
        double=double[double[:,0].argsort(),0:2] #Sorts by the first column
        splitdouble=np.vsplit(double,np.unique(double[:,0],return_index=True)[1]) #array split by first col
        self.splitdouble=self.splitdouble
        for array in double:
            if np.size(array)>max(0,threshold*2-1): 
                #Condition ensures don't waste time on arrays
                #too small to have threshold supporting links
                if np.size(array)==2:
                    key=array[0] #Get first entry
                    Graph[key]=[]
                else:
                    key=pd.unique(array[:,0]) #Should be one contig
                if len(key)==1:
                    key=key[0] #Should now be the string
                    Graph[key]=[]
                elif isinstance(key,str):
                    pass
                else:
                    print "Dunno what happened here"
                if np.size(array)==2:
                    counts=np.unique(array[1],return_counts=True)
                else:
                    counts=np.unique(array[:,1],return_counts=True)
                candidates=counts[0][counts[1]>=threshold]
                if len(candidates)>0:
                    Graph[key]=np.concatenate((Graph[key],candidates),axis=0)
        #Currently has duplicate eg contig1 - contig2 and contig2 to contig 1 as entries
        #Shouldn't matter - handled in scaffold creation
        if not dirty:
            self.graph=Graph
        return Graph
    
    def makescaffolds(self,Graph):
        '''Takes Graph of all connections and one with orientations,
        makes a dictionary of Scaffold objects for scaffoldset'''
        Scaffolds={}
        isolated=self.lonetigs(Graph)
        self.isolated=isolated
        for tig in isolated:
            #Add all contigs with no links as individual scaffolds
            #Also remove from both Graph and OrientedGraph
            Scaffolds["Scaffold"+str(len(Scaffolds))]={tig:[0,0,0,0,0]}
            del Graph[tig]
        #print Scaffolds
        #Now need to retrieve each scaffold with the order,gapsize,orientation
        paths,pathz=self.makepathss(Graph)
        print "Paths have been completed",datetime.datetime.now().time()
        self.finalgraph=self.pathstograph(pathz) #Turns final paths back into a graph
        #Now unpacks the paths to make scaffolds
        for path in paths:
            scafname="Scaffold"+str(len(Scaffolds))
            Scaffolds[scafname]={}
            startnode=None
            prevnode=None
            curnode=None
            i=0
            for node in path:
                if startnode==None:
                    startnode=node
                    curnode=node            
                else:
                    prevnode=curnode
                    curnode=node
                    Scaffolds[scafname][prevnode]=\
                        [0,0,self.arrange(prevnode,curnode),self.gapest(prevnode,curnode,final=True),i]
                    i+=1
            finalnode=curnode
            Scaffolds[scafname][finalnode]=[0,0,self.arrange(prevnode,curnode),self.gapest(finalnode),i]
        for scaf in Scaffolds: #Loop over the constructed scaffolds
            scaf1={} #Make a dict for storing  just one scaffold
            scaf1[scaf]=Scaffolds[scaf] #Dict with one entry -structure expected by scaffold.scaf
            self.scaffoldset[scaf]=scaffold.Scaffold(scaf1,scaf,linesize=70,contigloc=self.contigloc) #Makes one scaffold structure
            #and stores it in the dictionary for later calling
        return

    def arrange(self,prevnode,curnode,front=True):
        '''Takes in the previous node and currrent node to check whether
        a contig needs to be flipped. The flipping flag will trigger a down stream process
        to write a flag file with possible inversions. It is checked elsewhere if the flip will ensure
        the 0,1 read orientation. For, for the righthand side this would be needing a 0,0 to flip and LHS
        a 1,1. Note that arrange is decided for each pair and is such that if it needs to flip then the prevnode will be flipped'''
        links=self.readCounts(prevnode,curnode,clean=True)
        linkind=links.index(max(links)) #Assumes the maximum linked direction is used 
        #This assumption is throughout ScaffoldM
        orientation=sum(self.linkorientation(linkind)) #Check states
        #eg (1,0) and (0,1) don't need an inversion, (0,0) and (1,1) do.
        if orientation==2:
            return 1 #INdicates a need to invert a sequence
        elif orientation==0: #Issue - need to flip curnode not prevnode
            ##Might nto be picked up when considering next pair
            return 1 #INdicates a need to invert a sequence
        elif orientation==1:
            return 0 #Do nothing
        return None
        
    def reversecompliment(self,seq):
        compliment={'A':'T','G':'C','T':'A','C':'G','a':'t','g':'c','t':'a','c':'g'}
        try:
            newseq="".join([compliment[char] for char in reversed(seq)])
            return newseq
        except:
            print "This sequence has illegal characters"
            raise
    def makepathss(self,Graph):
        '''Scaffolds are equivalent to paths along the vertices of the graph.
        Here, scaffoldM uses the vertices to construct paths and returns them.
        The legality of paths is not checked
        Made with assumption that isolated tigs have been removed'''
        edges=self.makeedges(Graph)
        badedges={}
        connections={}
        paths=self.findpath(edges) #Makes the paths
        pathz=paths #Storing for later use
        try:
            if type(paths[0])==list:
                paths=[self.tuplecollapse(x) for x in paths if x!=None]
        except:
            print paths, "This is the paths"
            pass
        #print paths
        return paths,pathz
    def mergedicts(self,x,y):
        '''Shouldn't overwrite values in x with those in y as
        they shouldn't share any contigs and therefore do not share keys'''
        z=x.copy()
        z.update(y)
        return z
    def pathstograph(self,listofpaths):
        ''' Loops over list of paths, makes a graph for each path and joins
        with the final graph'''
        finalgraph={}
        for x in listofpaths:
            finalgraph.update(self.edgestograph(x))
        #print finalgraph
        return finalgraph
        
    def edgestograph(self,edges):
        '''Creates a graph of each edge which access all linked edges'''
        Graph={}
        for edge in edges:
            if edge[0] not in Graph and edge[1] not in Graph:
                Graph[edge[0]]=[edge[1]]
                Graph[edge[1]]=[edge[0]]
            elif edge[1] not in Graph:
                Graph[edge[1]]=[edge[0]]
            elif edge[0] not in Graph:
                Graph[edge[0]]=[edge[1]]
            else:
                pass
            for edge2 in edges:
                if edge2==edge:
                    pass
                elif edge2[1] in edge:
                    if edge2[1]==edge[0]:
                        if edge2[0] not in Graph[edge[0]]:
                            Graph[edge[0]]+=[edge2[0]]
                    elif edge2[1]==edge[1]:
                        if edge2[0] not in Graph[edge[1]]:
                            Graph[edge[1]]+=[edge2[0]]
                elif edge2[0] in edge:
                    if edge2[0]==edge[0]:
                        if edge2[1] not in Graph[edge[0]]:
                            Graph[edge[0]]+=[edge2[1]]
                    elif edge2[0]==edge[1]:
                        if edge2[1] not in Graph[edge[1]]:
                            Graph[edge[1]]+=[edge2[1]]
                else:
                    pass     
        return Graph
        
    def findpath(self,Graph):
        '''Keeps extending path until no more possible edges'''
        paths=[] #The current path being considered
        doneedges=set([]) #All previously considered edges
        donetigs=set([])  #All previously considered contigs
        for edge in Graph:
            done=False
            if edge not in doneedges and edge[0] not in donetigs and edge[1] not in donetigs:
                #Check if the edge has already been considered in this round and if the
                #component contigs have been completely considered
                path=[edge]
                doneedges|=set([edge]) #Consider current tuple
                doneedges|=set([edge[::-1]]) #Consider the reversed tuple
                if len(path)==1: #Force a orientation consideration on starting pair
                    counts=self.readCounts(edge[0],edge[1],clean=True)
                    index=np.arange(len(Counts))
                    index=index[np.logical_and(Counts>=threshold,Counts==max(Counts))]
                    orientation=self.linkorientation(index)
                    start=True
                    if orientation==(1,0):
                        path[0]=path[0][::-1] #Flips tuples
                    elif orientation==(1,1):
                        pass
                    elif orientation==(0,0):
                        pass
                while not done: #Extend the path
                    frontedge=[x for x in Graph if path[-1][1]==x[0] and (x not in doneedges) and (x[::-1] not in doneedges) ]#and x[::-1]!=path[-1]
                    backedge=[x for x in Graph if path[0][0]==x[1] and (x not in doneedges) and (x[::-1] not in doneedges) ]#and x[::-1]!=path[0]
                    newfront=self.joinedges(path,frontedge,front=True)
                    newback=self.joinedges(path,backedge,front=False)
                    if isinstance(newback,tuple) and isinstance(newfront,tuple):
                        path=[newback]+path+[newfront]
                        doneedges=doneedges | set([newback])| set([newfront])
                        doneedges=doneedges | set([newback[::-1]]) | set([newfront[::-1]])
                    elif isinstance(newback,tuple):
                        path=[newback]+path
                        doneedges=doneedges | set([newback])
                        doneedges=doneedges | set([newback[::-1]])
                    elif isinstance(newfront,tuple):
                        path=path+[newfront]
                        doneedges=doneedges | set([newfront])
                        doneedges=doneedges | set([newfront[::-1]])
                    elif isinstance(newfront,type(None)) and isinstance(newback,type(None)):
                        done=True
                #After finishing a path i.e no possible extensions
                for tup in path:
                    for tig in tup:
                        donetigs|=set([tig])
                #Added checks for sensible path
                if all(isinstance(edge, tuple) for edge in path):
                    if all(not isinstance(ele,bool) for edge in path for ele in edge):
                        paths.append(path)
                    #all elements are tuples
                    #print path
                    
                else: #Only append if appears like a valid path
                    print path, "Possible bugs to fix here"
                    #Ignore paths where there is a non-tuple object
                    pass
        return paths

    def joinedges(self,path,posedges=[],threshold=0.8,front=False):
        '''Includes sspace-like decision for joining edges.
        Takes a path containing ordered edges and decides which, if any,
        of the possible edges should be joined tot the path.
        Now also has a coverage based test. This does not form links but can aid
        in accepting or rejecting links which have passed earlier tests. In the cases
        of ambiguity no links are made, this should reduce false positives.'''
        newedge=None
        #best can =1 since earlier threshold for edges was 5
        maxcounts=[]
        best=0
        nextcheck=True #Whether or not next linking test is needed
        test1=True #Whether sspace-like decision 1 was passed
        test2=True #whether sspace-like decision 2 was passed
        #print posedges
        if posedges==[]: ##Consider if there were no possible edges
            #print path,"This is the Path"
            return None
        #Need to add check against those reads in the path
        covtest={} #Dictionary of coverage test results
        for edge in posedges:
            covtest[edge]=[]
            if front:
                #Note path[-1][1]==edge[0] by how the possible edges were decided
                counts=self.readCounts(path[-1][0],edge[0],clean=True)
                ori=self.linkorientation(counts.index(max(counts)))
                covtest[edge]+=[self.metabatprobtest(path[-1][0],edge[0])]
            elif not front:
                #path[0][0]==edge[1] by how possible edges are evaluated
                counts=self.readCounts(edge[0],path[0][0],clean=True) #Reorder ensures orientation for later
                ori=self.linkorientation(counts.index(max(counts)))
                covtest[edge]+=[self.metabatprobtest(path[0][0],edge[0])]
            else:
                pass
            #Legal joins for each position
            if (front and ori==(0,1)) or (front==False and ori==(0,1) or (front and ori==(0,0)) or (not front and ori==(1,1))):
                maxcounts.append((max(counts),edge))
        #Whether any of the coverages are suitable for use in decision making
        covaccept={key: [True  if x=='Accept'  else False for x in covtest[key]] for key in covtest}
        covreject={key: [True if x=='Reject' else False for x in covtest[key]] for key in covtest}
        #maxcounts.append((max(self.readCounts(path,edge),
        if maxcounts!=[]:
            best=max(maxcounts)
        if type(best)!=tuple:
            #print best, "This is the path for typeerror"
            #print path,"This is the Path for typeerror" ###PAY ATTENTION TO THIS ERROR - IT COULD BE QUITE IMPORTANT
            return None
        #print maxcounts, "Here I am bugfixing"
        for count,tup in maxcounts:
            if float(count)/best[0]>=threshold and (count,tup)!=best and not covaccept[tup][0] and not covaccept[best[1]][0]:
                #print best, "Compared to",(count,tup), "Inconclusive for both"
                nextcheck=False #Don't do the next check - failed first test
                test1=False #Failed first threshold
            elif float(count)/best[0]>=threshold and (count,tup)!=best and covaccept[tup][0] and covaccept[best[1]][0]:
                #print best, "Compared to",(count,tup), "Support for both"
                nextcheck=False #Don't do the next check - failed first test
                test1=False #Failed first threshold
            elif float(count)/best[0]>=threshold and (count,tup)!=best and not covaccept[tup][0] and covaccept[best[1]][0]:
                #If the current best option passes the covtest but possible edge does not
                print "Support for best but not tup"
                pass 
            elif float(count)/best[0]>=threshold and (count,tup)!=best and covreject[tup][0] and not covreject[best[1]][0]:
                #If coverage indicates the current pair shouldn't be considered, don't consider breaking the threshold a failure
                #print "Support against tup but not best"
                pass
            else:
                #print "It does not appear that coverage brought me here"
                pass
        if nextcheck: #Sequence space analysis - SSPACE-like decision 2
            spaceperlink=[]
            for edge in posedges:
                if front:
                    spaceperlink+=[self.seqspace(path[-1][0],edge[0])]
                elif front==False:
                    spaceperlink+=[self.seqspace(path[0][0],edge[0])]
                else:
                    pass
            biggest=max(spaceperlink)
            bigind=spaceperlink.index(biggest)
            ratios=[x/float(biggest) for i,x in enumerate(spaceperlink) if i!=bigind] #Look at all ratios but biggest/biggest
            
            if ratios!=[]:
                if max(ratios)>=threshold: #Second decision point
                    test2=False
        #Now for coverage based method
        #merge in above loop for 2nd sspace later

        #print "These are the accept values",accept
        #print "These are the reject values", reject
        if test1 and test2: #Pass both SSPACE tests
            #print "The coverage accept values", covaccept
            #print "The coverage reject values", covreject
            #print "SSPACE-like accept"
            return best[1]
        else:
            #print "SSPACE-like rejection"
            #print posedges
            #print path
            #print "test1",test1,"test2",test2
            return None
        ###Now need to check if edge satisfies coverage - red false +ve
        #If one is rejection Could recalc for the next biggest no. links
        #Could increase true +ve and dec false -ve - might be a mistake to include
        #Could also inc false +ve
        #print path,"This is the Path for too many competing options"
    def seqspace(self,tig1,tig2):
        ''' The second SSPACE-like step. Here the 'sequence space' is formed in the manner described
        by SSPACE and comparisons are made. This works by deeming the searched sequence space as being
        the same size as the insert for the library. Then, it works out how much of this sequence space is
        occupied by each contig (insertsize-gap between contigs) and the amount of sequence space
        per link. The sequence space per link is then compared and if the ratio is > threshold not links
        are made.'''
        ###Need to improve gapestimator for this to be effective.
        seqspaces=self.mean ##Need to check inserts for each set of links and weight mean by each one
        cover=[]
        N_rpbam=[]
        tempgap=self.gapest(tig1,tig2)
        for key in seqspaces:
            temp=self.getlinks(tig1,tig2,bam=key+".bam")
            if not isinstance(temp,type(None)):
                N_rpbam+=[len(self.getlinks(tig1,tig2,bam=key+".bam"))] #Number of links per bam library
                cover+=[min(self.tiglen[tig2],seqspaces[key][0]-tempgap)]
            #Using arithmetic mean here.
        mean=sum([x/float(N_rpbam[i]) for i,x in enumerate(cover)])/len(seqspaces) #Average seqspace per link over bams
        return mean

    
    def tuplecollapse(self,tuplist):
        '''Takes a list of tuples and then joins them into one longer tuple
        if they share edges'''
        x=None
        y=None
        orientation=0
        fuser=None
        if len (tuplist)==1:
            return tuplist[0]
        for tup in tuplist:
            y=x
            x=tup
            if not isinstance(x,type('h')) or not isinstance(y,type('h')):
                #Both have to be legal strings, i.e, contig names
                pass
            else:
                if isinstance(fuser,type(None)):
                    fuser=self.dubtuple(y,x)
                if all(not isinstance(ele,bool) for ele in y) and all(not isinstance(ele2,bool) for ele2 in x):
                    fuser=self.dubtuple(fuser,x)
                else:
                    pass
        return fuser
                
    def dubtuple(self,tup1,tup2):
        ''' Takes in a pair of tuples and then decides how to join them
        on a shared edge. Eg, join (B,A) and (D,B) => (D,A,B)'''
        eletup1=set(tup1)
        eletup2=set(tup2)
        if tup1[-1]==tup2[0]:
            new=tup1+tuple(eletup2-eletup1)
        elif tup1[-1]==tup2[-1]:
            new=tup1+tuple(eletup2-eletup1)
        elif tup1[0]==tup2[0]:
            new=tuple(eletup2-eletup1)+tup1
        elif tup1[0]==tup2[-1]:
            new=tuple(eletup2-eletup1)+tup1
        else:
            new=False
        return new

    def rearrange(self,lists,item1,item2):
        '''swap two items in a list'''
        a, b = i.index(item1), i.index(item2)
        lists[b], lists[a] = lists[a], lists[b]
        return lists
        
    def flipdec(self,linkdirs):
        ''' Return tuple indicating how a pair of linked
        contigs should be arranged based on their read orientations'''
        flipplace=False
        or1,or2=linkdirs
        or1=int(or1)
        or2=int(or2)
        if (or1+or2)!=1:
            flipplace=True
        return flipplace
                
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
        ''' Trims off reads with unrealistic large reads for the library
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
            #print name
            #Should give unique insert mean and insert std
            meaninsert[name]=np.atleast_2d(insertfile)[np.atleast_2d(insertfile)[:,0]==name+".bam",1].astype(float)
            stdinsert[name]=np.atleast_2d(insertfile)[np.atleast_2d(insertfile)[:,0]==name+".bam",2].astype(float)
            #meaninsert[name]=[float(x[1]) for x in insertfile if name+".bam" in x]
            #stdinsert[name]=[float(x[2]) for x in insertfile if name+".bam" in x]
        self.mean=meaninsert
        self.std=stdinsert
        cleanedlinks=linksfile[self.veclowerdist(linksfile)<self.veccutoff(linksfile),:]
        self.cleanedlinks=cleanedlinks
        print np.size(cleanedlinks), "The new Size"
        print np.size(linksfile), "the old size"
        return
        
    def parse(self):
        import datetime
        import sys
        '''links all functions to start from the DataParser input and to end
        with the processed scaffolds(with gap estimate and proper orientation).
        It will also call all of the scaffolds to print to the scaffold file.'''
        self.cleanlinks() #Cleans the obviously erroneous mappings
        print "Links have been cleaned"
        print datetime.datetime.now().time()
        Graph=self.makegraph() #Makes initial connections with no.reads>threshold
        print "This initial graph has been made"
        print datetime.datetime.now().time()
        sys.stdout.flush()
        self.makescaffolds(Graph) #Makes set of scaffolds
        print "The scaffolds are now done"
        print datetime.datetime.now().time()
        self.printscaffolds()
        print "The scaffolds have been printed"
        print datetime.datetime.now().time()
        self.gapprint()
        print "The gaps have been printed"
        print datetime.datetime.now().time()
        self.totalgraph=self.makegraph(connections=self.completecheck(self.links),\
        cleanedlinks=self.links,threshold=0,dirty=True)
        self.stagestosif()
        print "The sif/network type files for cytoscape have been completed"
        print datetime.datetime.now().time()

    def printscaffolds(self,filename='testScaffold'):
        try:
            for key in self.scaffoldset:
                self.scaffoldset[key].printscaffold(filename+".fasta")
        except AttributeError:
            print 'These scaffolds do not appear to be the scaffold class'
        return

    def lowerdist(self,onelink):
        '''severe lower bound on gap between
         reads (ignores gap between contigs)'''
        orientation1=int(onelink[4])
        orientation2=int(onelink[7])
        if orientation1==1:
            dist1=int(onelink[3])
        else:
            dist1=int(onelink[2])-int(onelink[3])
        if orientation2==1:
            dist2=int(onelink[6])
        else:
            dist2=int(onelink[5])-int(onelink[6])
        return (dist1+dist2)
    def veclowerdist(self,alllinks):
        '''severe lower bound on gap between
         reads (ignores gap between contigs)'''
        size=len(alllinks[:,1])
        dist1=np.zeros(size)
        dist2=np.zeros(size)
        #Get the indices for each orientation pair
        ori11=np.arange(size)[alllinks[:,4].astype(int)==1]
        ori10=np.arange(size)[alllinks[:,4].astype(int)==0]
        ori21=np.arange(size)[alllinks[:,7].astype(int)==1]
        ori20=np.arange(size)[alllinks[:,7].astype(int)==1]
        dist1[ori11]=alllinks[ori11,3].astype(int)
        dist1[ori10]=alllinks[ori10,2].astype(int)-alllinks[ori10,3].astype(int)
        dist2[ori21]=alllinks[ori21,6].astype(int)
        dist2[ori20]=alllinks[ori20,5].astype(int)-alllinks[ori20,6].astype(int)
        return (dist1+dist2)
        
    def veccutoff(self,alllinks,cutoff=3):
        size=len(alllinks[:,1])
        thresh=np.zeros(size)
        indices=np.arange(size)
        bamindices={bam:indices[alllinks[:,8]==bam+".bam"] for bam in self.bamnames}
        for name in self.bamnames:
            thresh[bamindices[name]]=self.mean[name][0]\
            +cutoff*self.std[name][0]
        return thresh
        
    def makeedges(self,graph):
        '''edges as (node1,node2) tuples.'''
        edges=[]
        for key in graph:
            for point in graph[key]:
                if (key,point) not in edges:
                    edges.append((key,point))
        return edges
        
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

    def gapest(self,contig1=False,contig2=False,cleanedlinks=None,bamnames=None,default=False,final=False):
        '''Takes in the cleaned links and then uses them to estimate the gap between two contigs
        in the manner presented in Sahlini et al. 2012. There are issues where if the gap
        is too large compared to the insert size then the gapestimator is not effective. This
        is covered by testing for conditions where the estimator would fail and putting in a default value.'''
        if bamnames==None:
            bamnames=self.bamnames
        if default:
            if final:
                self.gaps[(contig1,contig2)]=200
            return 200
        if not contig2:
            return 0
    #implementation of sahlini et al 2013 algorithm.
        import os
        import sys
        from scipy.special import erf
        from scipy.stats import norm
        from scipy.constants import pi
        from math import exp
        from scaffold import Scaffold
        ##Given sufficiently large contigs then
        ##Can use binary search on the g(d) function to find 
        ##THe ML estimate
        #Data
        #observations=fsomething(cleanedlinks)
        nullscaf=Scaffold({},'Null',self.contigloc)
        meaninsert=self.mean
        stdinsert=self.std
        observations=self.veclowerdist(self.getlinks(contig1,contig2,clean=True))
        #observations=[self.lowerdist(x) for x in cleanedlinks if (contig1 in x) and (contig2 in x)]
        #print observations
        r=self.readlen #read length
        m=3 #Arbitrary tolerance level
        sigma=[float(stdinsert[x][0]) for x in stdinsert] #sd of insert library
        mu=[float(meaninsert[x][0]) for x in meaninsert]#Mean of insert library
        l=10 #smallest correct gap expected to see
        if contig1!=False and contig2!=False:
            c1=int(self.tiglen[contig1])#length of contig 1
            c2=int(self.tiglen[contig2])#Length of contig 2
            cmin=min(c1,c2) #Minimum length
            cmax=max(c1,c2) #Maximum length
            if len(observations)==0:
                print "There has been an error"
                return 1
            distance=self.MLsearch(cmin,cmax,r,mu[0],sigma[0],observations)
            #print "Tig 1:",contig1,"Tig2",contig2
            #print "The ML distance is {0}".format(distance)
            if final:
                self.gaps[(contig1,contig2)]=distance
        return distance

    ##Estimator functions from below are almost exact replicates of those seen in GapCalculator.py as
    ##Released by Sahlin et al 2014 on https://github.com/GapEst
    def Denominator(self,d,c_min,c_max,readLen,mean,stdDev):
        from scipy.special import erf
        from scipy.stats import norm
        from scipy.constants import pi
        from math import exp
        scale=(2**0.5*float(stdDev))
        #term 1,2 and 3 denodes what part of the function we are integrating term1 for first (ascending), etc...
        term2=(c_min-readLen+1)/2.0*(erf((c_max+d+readLen-mean)/(scale))- erf((c_min+d+readLen-mean)/(scale))   )
        #term2=((pi/2)**0.5)*(c_min-readLen+1)*(erf((c_max+d+readLen-mean)/(scale))- erf((c_min+d+readLen-mean)/(scale))   )
        first=-(d+2*readLen-mean-1)/2.0*( erf((c_min+d+readLen-mean)/(scale)) - erf((d+2*readLen-1-mean)/(scale))  )
        #first=-((pi/2)**0.5)*(d+2*readLen-mean-1)*( erf((c_min+d+readLen-mean)/(scale)) - erf((d+2*readLen-1-mean)/(scale))  )
        second=(stdDev/((2*pi)**0.5))*( exp(-( (d+2*readLen-1-mean)**2)/(scale**2)) - exp(-( (c_min+d+readLen-mean)**2)/(scale**2)))
        #second=stdDev*( exp(-( (d+2*readLen-1-mean)**2)/(scale**2)) - exp(-( (c_min+d+readLen-mean)**2)/(scale**2))) 
        term1=first+second

        first=(c_min+c_max+d-mean+1)/2.0*( erf((c_min+c_max+d-mean+1)/(scale)) - erf((c_max+readLen+d-mean)/(scale))  )
        #first=((pi/2)**0.5)*(c_min+c_max+d-mean+1)*( erf((c_min+c_max+d-mean)/(scale)) - erf((c_max+readLen+d-mean)/(scale))  )
        #print 'First: ',first
        second=(stdDev/((2*pi)**0.5))*( exp(-( (c_min+c_max+d-mean+1)**2)/(scale**2)) - exp(-( (c_max+readLen+d-mean)**2)/(scale**2)))
        #second=stdDev*( exp(-( (c_min+c_max+d-mean)**2)/(scale**2)) - exp(-( (c_max+readLen+d-mean)**2)/(scale**2)))
        #print 'Second: ',second
        term3=first+second
        denom=term1+term2+term3
        #print denom, "Denominator"
        #print term1,term2,term3
        return denom
    def gprimed(self,cmin,cmax,r,d,mu,sigma):
        from scipy.special import erf
        from scipy.stats import norm
        from scipy.constants import pi
        from math import exp
        scale=(2**0.5*float(sigma))
        #Straight from GapCalculator.py in Sahlin et al's code - rename
        #changes
        #erf((mu-d-cmax-r-1)/(2**0.5*float(sigma)))) to erf((mu-d-cmax-r)/(2**0.5*float(sigma))))
        #erf((cmin+d+r+1-mu)/(2**0.5*float(sigma))) to erf((cmin+d+r-mu)/(2**0.5*float(sigma)))
        #num1=( erf((cmin+cmax+d+1-mu)/(scale)) +erf((mu-d-cmax-r-1)/(scale)))*(pi/2)**0.5
        #num2=-(erf((cmin+d+r-mu+1)/(scale))+erf((mu-d-2*r+1)/(scale)))*(pi/2)**0.5
        num1=1/2.0*(erf((cmin+cmax+d+1-mu)/float(scale))+erf((d+2*r-1-mu)/(scale)))
        num2=-1/2.0*(erf((cmax+d+r-mu)/(scale))+erf((cmin+d+r-mu)/(scale)))
        num=num1+num2 #Complete g prime
        return num
        
        
    def fd(self,cmin,cmax,r,d,mu,sigma):
        numer=self.gprimed(cmin,cmax,r,d,mu,sigma)
        denom=self.Denominator(d,cmin,cmax,r,mu,sigma)
        #~ denom=self.gofd(d,cmin,cmax,r,mu,sigma)
        fofd=d+(numer/float(denom))*sigma**2
        #print "This is the quotient", fofd-d, "Will it be super small?"
        #print "This is the current gapsize being considered", d
        return fofd
        
    def MLsearch(self,cmin,cmax,r,mu,sigma,observations,m=3,bayesian=True):
        noLinks=len(observations)
        #~ print observations
        #~ print noLinks
        #~ print mu
        l=10
        #Heuristics cut off from Sahlin et al.
        #mean of 10 largest - mean of 10 smallest <6*stdDev
        #At least 10 links mapping
        obsval=(noLinks*mu-int(sum(observations)))/float(noLinks)-r #offset for 1 read from RF pair
        #obsval=(int(sum(observations)))/float(noLinks)
        #print "This is the sum of the observed distance", obsval
        #d_lower=max(-l,mu-cmin-cmax-m*sigma+2*r)
        #d_lower=-1000
        #d_upper=mu+m*sigma+2*r
        if not bayesian:
            d_upper=max(int(mu+4*sigma-2*r), 2*r) #from sahlini directly
            d_lower=min(-int(mu+4*sigma-2*r), -2*r)
            #Can do a binary search since f_d should be monotically increasing
            #As mentioned by Sahlin et al.
            while d_upper-d_lower>1:
                d_ML=(d_upper+d_lower)/2.0
                #print d_ML
                f_d=self.fd(cmin,cmax,r,d_ML,mu,sigma)
                if f_d>obsval:
                    d_upper=d_ML
                else:
                    d_lower=d_ML
            Gapsize=int(round((d_upper+d_lower)/2.0,0))
            return Gapsize
        elif bayesian:
            d_upper=max(int(mu+4*sigma-2*r), 2*r) #from sahlini directly
            d_lower=min(-int(mu+4*sigma-2*r), -2*r)
            #Can do a binary search since f_d should be monotically increasing
            #As mentioned by Sahlin et al.
            while d_upper-d_lower>1:
                d_MAP=(d_upper+d_lower)/2.0
                #print d_MAP
                f_d=self.bayesianeqn(cmin,cmax,r,d_MAP,mu,sigma)
                if f_d>obsval:
                    d_upper=d_MAP
                else:
                    d_lower=d_MAP
            Gapsize=int(round((d_upper+d_lower)/2.0,0))
            return Gapsize
        
    def bayesianeqn(self,cmin,cmax,r,d,mu,sigma):
        from scipy.stats import norm
        #Again from Sahlini
        ML_value=self.fd(cmin,cmax,r,d,mu,sigma)
        prior_value=sigma**2*(1/(1-norm.cdf(d,mu,sigma))*norm.pdf(d,mu,sigma))
        return ML_value+prior_value

    def metabatprobtest(self,contig1,contig2,lower=0.01,upper=0.95):
        pairs=self.gettigcov(contig1,contig2)
        prob=1
        j=0
        for i,pair in enumerate(pairs):
            repack=zip(*pair)
            #print repack
            if 0 not in repack[1]:
                prob*=self.distancecalc(repack[0],repack[1])
                j+=1
            if 0 in repack[1]: #Conservative method
                #For dealing with probabilities
                prob*=1
                j+=1 
        prob=prob**(1/float(j)) #Geometric mean of probabilities
        #Make a decision based on the probabilties
        if prob<=lower:
            return "Accept"
        elif prob>=upper:
            return "Reject"
        else:
            return "Uninformative"
        
    def gettigcov(self,contig1,contig2):
        covs={}
        covs[contig1]=[]
        covs[contig2]=[]
        for (key,val) in self.coverages.iteritems():
            covs[contig1]+=[val[contig1]]
            covs[contig2]+=[val[contig2]]
        pairs=zip(covs[contig1],covs[contig2])
        return pairs
            
    def gapprint(self):
        #print "CAN YOU HEAR ME UP THERE"
        import os
        if not os.path.isfile("Gapdata.txt"):
            with open("Gapdata.txt",'w') as Gaps:
                Gaps.write("Contig1,Contig2,Orientation,gapestimate\n")
        with open("Gapdata.txt",'a+') as Gaps:
            for tigtuple in self.gaps:
                Gaps.write("{0},{1},Default,{2}\n".format(tigtuple[0],tigtuple[1],self.gaps[tigtuple]))
        return
        
    def flagprint(self):
        
        return

    def scafsum(self):
        self.gapprint()
        return

    def output(self):
        self.gapprint()
        self.scafsum()
        return
        
    def distancecalc(self,means,var,tol=10**(-9)):
        ''' sums the difference in probability density between theorectical
        coverage distributions untill changes are less than some predefined tolerance''' 
        import random 
        from scipy.stats import norm
        import numpy as np
        from math import sqrt
        start=0
        end=0
        dist=0
        i=0
        change=1
        while abs((change))>tol or i<max(means)+100: # or i<max(means)
            #if i==max(means)-1:
                #print "Leaving the loop soon"
            change=abs(norm.pdf(i,loc=means[0],scale=sqrt(var[0]))-norm.pdf(i,loc=means[1],scale=sqrt(var[1])))
            dist+=change
            i+=1
        #print i
        return dist/float(2)
        
    def writesometigs(self):
        with open('ExtractionList.fna','w') as extract:
            for tig in self.contignames:
                extract.write(">{0}\n".format(tig))
        return
