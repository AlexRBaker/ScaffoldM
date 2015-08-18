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
###############################################################################

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
                 contigloc):
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
        self.contigloc=contigloc
        self.readlen=100

    def readCounts(self,contig1,contig2,linksfile=None):
        '''Gets the number of links between two contigs
        in each of 4 orientations, (0,1),(1,0),(1,1),(0,0)
        Here 0 is a read in the same direction as the contig'''
        linksfile=self.getlinks(contig1,contig2)
        #Changes from needing to maintain relative arrangement
        tig1ind=linksfile[0].index(contig1)
        tig2ind=linksfile[0].index(contig2)
        #Now maintains appropiate order, if refering to contigs in seperate arrangements
        #This is since index determines if x[4] refers to read of tig 1 or tig 2.
        if tig1ind<tig2ind:
            norm=len([x for x in linksfile if (contig1 in x) and (contig2 in x) and (x[4]=='1' and x[7]=='0')])
            rev=len([x for x in linksfile if (contig1 in x) and (contig2 in x) and (x[4]=='0' and x[7]=='1')])
        elif tig2ind<tig1ind:
            norm=len([x for x in linksfile if (contig1 in x) and (contig2 in x) and (x[4]=='0' and x[7]=='1')])
            rev=len([x for x in linksfile if (contig1 in x) and (contig2 in x) and (x[4]=='1' and x[7]=='0')])
        else:
            norm=0
            rev=0
            print "How is this possible?"
        posinvc2=len([x for x in linksfile if (contig1 in x) and (contig2 in x) and (x[4]=='1' and x[7]=='1')])
        posinvc1=len([x for x in linksfile if (contig1 in x) and (contig2 in x) and (x[4]=='0' and x[7]=='0')])
        total=norm+rev+posinvc2+posinvc1
        return [norm,rev,posinvc2,posinvc1,total]
        
    def getlinks(self,contig1,contig2=False,linksfile=None):
        '''Retrieves every link entry for one contig for every link 
        between two contigs'''
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
        '''Checks all contigs linked to the specified contig
        and return a list of these contigs.
        '''
        if linksfile==None:
            linksfile=self.cleanedlinks
        links=self.getlinks(contig1,False,linksfile)
        notcontig1=set([x[0] for x in links if x[0]!=contig1])
        nottig1=set([x[1] for x in links if x[1]!=contig1])
        return list(nottig1|notcontig1)

    def completecheck(self,cleanedlinks=None):
        '''Uses CheckLinks to Construct all a dictionary with a key for each contigs
        which contains all contigs to which they are linked.'''
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
        #print Graph
        #print self.makeedges(Graph)
        #print OrientedGraph
        return [OrientedGraph,Graph]

    def makescaffolds(self,Graph,OrientedGraph):
        '''Takes Graph of all connections and one with orientations,
        makes a dictionary of Scaffold objects for scaffoldset'''
        Scaffolds={}
        isolated=self.lonetigs(Graph)
        for tig in isolated:
            #Add all contigs with no links as individual scaffolds
            #Also remove from both Graph and OrientedGraph
            Scaffolds["Scaffold"+str(len(Scaffolds))]={tig:[0,0,0,self.gapest(tig),0]}
            Graph.pop(tig,None)
            OrientedGraph.pop(tig,None)
        print Scaffolds
        #Now need to retrieve each scaffold with the order,gapsize,orientation
        paths=self.makepathss(Graph,OrientedGraph)
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
                        [0,0,0,self.gapest(prevnode,curnode),i]
                    i+=1
            finalnode=curnode
            Scaffolds[scafname][finalnode]=[0,0,0,self.gapest(finalnode),i]
        for scaf in Scaffolds:
            scaf1={}
            scaf1[scaf]=Scaffolds[scaf]
            self.scaffoldset[scaf]=scaffold.Scaffold(scaf1,scaf)
        return

    
    def makepathss(self,Graph,OrientedGraph):
        '''Scaffolds are equivalent to paths along the vertices of the graph.
        Here, scaffoldM uses the vertices to construct paths and returns them.
        The legality of paths is not checked
        Made with assumption that isolated tigs have been removed'''
        edges=self.makeedges(Graph)
        #Nedges=self.edgestograph(edges)
        badedges={}
        connections={}
        paths=self.findpath(edges) #Makes the paths
        if type(paths[0])==list:
            paths=[self.tuplecollapse(x) for x in paths if x!=None]
        print paths
        return paths
        
    def edgestograph(self,edges):
        '''Creates a graph of each edge which access all linked edges'''
        Graph={}
        for edge in edges:
            Graph[edge]=[]
            for edge2 in edges:
                if edge2==edge:
                    pass
                elif edge2[1] in edge or edge2[0] in edge:
                    Graph[edge]=Graph[edge]+[edge2]
                else:
                    pass
        return Graph
        
    def findpath(self,Graph):
        '''Keeps extending path until no more possible edges'''
        paths=[] #The current path being considered
        doneedges=set([])
        for edge in Graph:
            path=[edge]
            done=False
            if edge not in doneedges:
                doneedges|=set([edge])
                doneedges|=set([edge[::-1]])
                if len(path)==1: #Force a orientation consideration on starting pair
                    counts=self.readCounts(edge[0],edge[1])
                    index=counts.index(max(counts))
                    orientation=self.linkorientation(index)
                    if orientation==(1,0):
                        path[0]=path[0][::-1] #Flips tuples
                while not done: #Extend the path
                    frontedge=[x for x in Graph if path[-1][1]==x[0] and (x not in doneedges) and (x[::-1] not in doneedges) ]#and x[::-1]!=path[-1]
                    backedge=[x for x in Graph if path[0][0]==x[1] and (x not in doneedges) and (x[::-1] not in doneedges) ]#and x[::-1]!=path[0]
                    newfront=self.joinedges(path,frontedge,front=True)
                    newback=self.joinedges(path,backedge,front=False)
                    if newback!=None and newfront!=None:
                        path=[newback]+path+[newfront]
                        doneedges=doneedges | set([newback])| set([newfront])
                        doneedges=doneedges | set([newback[::-1]]) | set([newfront[::-1]])
                    elif newback!=None:
                        path=[newback]+path
                        doneedges=doneedges | set([newback])
                        doneedges=doneedges | set([newback[::-1]])
                    elif newfront!=None:
                        path=path+[newfront]
                        doneedges=doneedges | set([newfront])
                        doneedges=doneedges | set([newfront[::-1]])
                    else:
                        done=True
                
                paths.append(path)
        #print paths
        return paths

    def joinedges(self,path,posedges=[],threshold=0.8,front=False):
        '''Includes sspace-like decision for joining edges.
        Takes a path containing ordered edges and decides which, if any,
        of the possible edges should be joined tot the path'''
        newedge=None
        #best can =1 since earlier threshold for edges was 5
        maxcounts=[]
        best=0
        print posedges
        if posedges==[]: ##Consider if there were no possible edges
            print path,"This is the Path"
            return None
        #Need to add check against those reads in the path
        for edge in posedges:
            if front:
                counts=self.readCounts(path[-1][0],edge[0])
                ori=self.linkorientation(counts.index(max(counts)))
            elif front==False:
                counts=self.readCounts(edge[0],path[0][0]) #Reorder ensures orientation for later
                ori=self.linkorientation(counts.index(max(counts)))
            else:
                pass
            if (front and ori==(0,1)) or (front==False and ori==(0,1)):
                maxcounts.append((max(counts),edge))
        #maxcounts.append((max(self.readCounts(path,edge),
        if maxcounts!=[]:
            best=max(maxcounts)
        if type(best)!=tuple:
            print path,"This is the Path for typeerror"
            return None
        for count,tup in maxcounts:
            if float(count)/best[0]>threshold and (count,tup)!=best:
                print path,"This is the Path for too many competing options"
                return None
        return best[1]

    def tuplecollapse(self,tuplist):
        x=None
        y=None
        orientation=0
        fuser=None
        if len (tuplist)==1:
            return tuplist[0]
        for tup in tuplist:
            y=x
            x=tup
            if x==None or y==None:
                pass
            elif fuser==None:
                fuser=self.dubtuple(y,x)
            else:
                fuser=self.dubtuple(fuser,x)
        return fuser
                
    def dubtuple(self,tup1,tup2):
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
        
    def arrangenplace(self,linkdirs):
        ''' Return tuple indicating how a pair of linked
        contigs should be arranged based on their read orientations'''
        flipplace=False
        flipsequence=False
        or1,or2=linkdirs
        or1=int(or1)
        or2=int(or2)
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
        return (meaninsert,stdinsert)
        
    def parse(self):
        '''links all functions to start from the DataParser input and to end
        with the processed scaffolds(with gap estimate and proper orientation).
        It will also call all of the scaffolds to print to the scaffold file.'''
        self.cleanlinks() #Cleans the obviously erroneous mappings
        OrientedGraph,Graph=self.makegraph() #Makes initial connections with no.reads>threshold
        
        self.makescaffolds(Graph,OrientedGraph) #Makes set of scaffolds
        self.printscaffolds()

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

    def gapest(self,contig1=False,contig2=False,cleanedlinks=None,bamnames=None,default=False):
        if cleanedlinks==None:
            cleanedlinks=self.cleanedlinks
        if bamnames==None:
            bamnames=self.bamnames
        if default:
            return 200
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
        meaninsert={}
        stdinsert={}
        for name in bamnames:
            #Should give unique insert mean and insert std
            meaninsert[name]=[x[1] for x in self.inserts if name+".bam" in x]
            stdinsert[name]=[x[2] for x in self.inserts if name+".bam" in x]
        observations=[self.lowerdist(x) for x in cleanedlinks if (contig1 in x) and (contig2 in x)]
        #print observations
        r=self.readlen #read length
        m=3 #Arbitrary tolerance level
        sigma=[float(stdinsert[x][0]) for x in stdinsert] #sd of insert library
        mu=[float(meaninsert[x][0]) for x in meaninsert]#Mean of insert library
        l=10 #smallest correct gap expected to see
        if contig1!=False and contig2!=False:
            c1=len(nullscaf.extractcontigs(contig1,self.contigloc,header=False).translate(None,'\n'))#length of contig 1
            c2=len(nullscaf.extractcontigs(contig2,self.contigloc,header=False).translate(None,'\n'))#Length of contig 2
            cmin=min(c1,c2) #Minimum length
            cmax=max(c1,c2) #Maximum length
            if len(observations)==0:
                print "There has been an error"
                return 1
            distance=self.MLsearch(cmin,cmax,r,mu[0],sigma[0],observations)
            print "Tig 1:",contig1,"Tig2",contig2
            print "The ML distance is {0}".format(distance)
        if contig2==False:
            distance=0
        return distance

    #~ def gofd(self,d,cmin,cmax,r,mu,sigma):
        #~ #Couldn't resolve difference between my and Sahlin's implementation - did they make further steps from paper?
        #~ #Some differences didn't match paper - check later
        #~ # The g (d) function from Sahlin et al 2013. This function
        #~ from scipy.special import erf
        #~ from scipy.stats import norm
        #~ from scipy.constants import pi
        #~ from math import exp
        #~ c_min=cmin
        #~ c_max=cmax
        #~ readLen=r
        #~ mean=mu
        #~ stdDev=sigma
        #~ # Terms are grouped by brackets in Sahlin et al. 2013
        #~ Term1=((c_min-readLen+1)/2.0)*(erf((c_max+d+readLen-mean)/((2**0.5)*stdDev))-erf((c_min+d+readLen-mean)/((2**0.5)*stdDev)))
        #~ Term2=((c_min+c_max+d+1-mean)/2.0)*(erf((c_min+c_max+d+1-mean)/((2**0.5)*stdDev))-erf((c_max+d+readLen-mean)/((2**0.5)*stdDev)))
        #~ Term3=((d+2*readLen-1-mean)/2.0)*(erf((d+2*readLen-1-mean)/((2**0.5)*stdDev))-erf((c_min+d+readLen-mean)/((2**0.5)*stdDev)))
        #~ Term4=(stdDev/((2*pi)**0.5))*(exp(-(cmax+cmin+d+1-mean)**2/(2*float(stdDev**2)))+exp(-(d+2*readLen-1-mean)**2/(2*float(stdDev**2)))-exp(-(c_max+d+readLen-mean)**2/(2*float(stdDev**2)))-exp(-(c_min+d+readLen-mean)**2/(2*float(stdDev**2))))
        #~ denom=Term1+Term2+Term3+Term4
        #~ print denom, "Denominator"
        #~ return denom
    ##Estimator functions from below are almost exact replicates of those seen in GapCalculator.py as
    ##Released by Sahlin et al 2014 on https://github.com/GapEst
    def Denominator(self,d,c_min,c_max,readLen,mean,stdDev):
        from scipy.special import erf
        from scipy.stats import norm
        from scipy.constants import pi
        from math import exp
        #term 1,2 and 3 denodes what part of the function we are integrating term1 for first (ascending), etc...
        term2=(c_min-readLen+1)/2.0*(erf((c_max+d+readLen-mean)/((2**0.5)*stdDev))- erf((c_min+d+readLen-mean)/((2**0.5)*stdDev))   )

        first=-((pi/2)**0.5)*(d+2*readLen-mean-1)*( erf((c_min+d+readLen-mean)/(2**0.5*float(stdDev))) - erf((d+2*readLen-1-mean)/(2**0.5*float(stdDev)))  )
        second=stdDev*( exp(-( (d+2*readLen-1-mean)**2)/(float(2*stdDev**2))) - exp(-( (c_min+d+readLen-mean)**2)/(float(2*stdDev**2)))) 
        term1=first+second

        first=((pi/2)**0.5)*(c_min+c_max+d-mean+1)*( erf((c_min+c_max+d-mean)/(2**0.5*float(stdDev))) - erf((c_max+readLen+d-mean)/(2**0.5*float(stdDev)))  )
        #print 'First: ',first
        second=stdDev*( exp(-( (c_min+c_max+d-mean)**2)/(float(2*stdDev**2))) - exp(-( (c_max+readLen+d-mean)**2)/(float(2*stdDev**2))))
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
        #Straight from GapCalculator.py in Sahlin et al's code - rename
        #num1=( erf((cmin+cmax+d+1-mu)/(2**0.5*float(sigma))) +erf((mu-d-cmax-r-1)/(2**0.5*float(sigma))))*(pi/2)**0.5
        #num2=-(erf((cmin+d+r+1-mu)/(2**0.5*float(sigma)))+erf((mu-d-2*r+1)/(2**0.5*float(sigma))))*(pi/2)**0.5
        num1=1/2.0*(erf((cmin+cmax+d+1-mu)/float(2**0.5*sigma))+erf((d+2*r-1-mu)/(2**0.5*float(sigma))))
        num2=-1/2.0*(erf((cmax+d+r-mu)/(2**0.5*sigma))+erf((cmin+d+r-mu)/(2**0.5*pi)))
        num=num1+num2 #Complete g prime
        return num
        
    def fd(self,cmin,cmax,r,d,mu,sigma):
        numer=self.gprimed(cmin,cmax,r,d,mu,sigma)
        denom=self.Denominator(d,cmin,cmax,r,mu,sigma)
        #~ denom=self.gofd(d,cmin,cmax,r,mu,sigma)
        fofd=d+(numer/float(denom))*sigma**2
        return fofd
        
    def MLsearch(self,cmin,cmax,r,mu,sigma,observations,m=3):
        noLinks=len(observations)
        #~ print observations
        #~ print noLinks
        #~ print mu
        l=10
        #Heuristics cut off from Sahlin et al.
        #mean of 10 largest - mean of 10 smallest <6*stdDev
        #At least 10 links mapping
        # 
        obsval=(noLinks*mu-sum(observations))/float(noLinks)
        d_lower=max(-l,mu-cmin-cmax-m*sigma+2*r)
        d_upper=mu+m*sigma+2*r
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
