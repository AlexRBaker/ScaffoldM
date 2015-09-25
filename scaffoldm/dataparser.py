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
        self.bamnames=bamnames
        self.contignames=contignames
        self.scaffoldset={}
        self.contigloc=contigloc
        self.tiglen=tiglen
        self.readlen=readsize
        self.gaps={}

    def readCounts(self,contig1,contig2,linksfile=None,clean=False):
        '''Gets the number of links between two contigs
        in each of 4 orientations, (0,1),(1,0),(1,1),(0,0)
        Here 0 is a read in the same direction as the contig'''
        linksfile=self.getlinks(contig1,contig2,clean=clean)
        #Changes from needing to maintain relative arrangement
        tig1ind=linksfile[0].index(contig1)
        tig2ind=linksfile[0].index(contig2)
        #Now maintains appropiate order, if refering to contigs in seperate arrangements
        #This is since index determines if x[4] refers to read of tig 1 or tig 2.
        if tig1ind<tig2ind: #Checks how the contig1 versus contig2 are place in readsfile
            norm=len([x for x in linksfile if (contig1 in x) and (contig2 in x) and (x[4]=='1' and x[7]=='0')])
            rev=len([x for x in linksfile if (contig1 in x) and (contig2 in x) and (x[4]=='0' and x[7]=='1')])
        elif tig2ind<tig1ind:
            norm=len([x for x in linksfile if (contig1 in x) and (contig2 in x) and (x[4]=='0' and x[7]=='1')])
            rev=len([x for x in linksfile if (contig1 in x) and (contig2 in x) and (x[4]=='1' and x[7]=='0')])
        else:
            norm=0
            rev=0
            print "How is this possible? contig1!=contig2"
        posinvc2=len([x for x in linksfile if (contig1 in x) and (contig2 in x) and (x[4]=='1' and x[7]=='1')])
        posinvc1=len([x for x in linksfile if (contig1 in x) and (contig2 in x) and (x[4]=='0' and x[7]=='0')])
        total=norm+rev+posinvc2+posinvc1
        return [norm,rev,posinvc2,posinvc1,total]
    def getlinks(self,contig1,contig2=False,linksfile=None,clean=False,bam=False):
        '''Retrieves every link entry for one contig for every link 
        between two contigs'''
        if linksfile==None:
            if clean:
                linksfile=self.cleanedlinks
            else:
                linksfile=self.links
        try:
            if contig2==False:
                if bam==False:
                    Links=[row for row in linksfile if contig1 in row]
                else:
                    Links=[row for row in linksfile if (contig1 in row) and (bam in row) ]
            else:
                if bam==False:
                    Links=[row for row in linksfile if (contig1 in row) and (contig2 in row)]
                else:
                    Links=[row for row in linksfile if (contig1 in row) and (contig2 in row) and (bam in row)]
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
        links=self.getlinks(contig1,False,linksfile,clean=True)
        if links==None:
            return None
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
        
    def graphtosif(self,graph,graphname,dirty=False,final=False):
        header="Node1   Relationship  Node2 links\n"
        with self.tryopen(graphname,header,".sif") as network:
            for contig1,connected in graph.iteritems():
                for contig2 in connected:
                    if dirty:
                        Counts=self.readCounts(contig1,contig2,clean=False)
                        for i,links in enumerate(Counts):
                            if links!=0:
                                network.write("{0}  {1} {2} {3}\n".format(contig1,\
                                self.linkorientation(i),contig2,links)
                    elif not dirty:
                        Counts=self.readCounts(contig1,contig2,clean=True)
                        network.write("{0}  {1} {2} {3}\n".format(contig1,\
                        self.orientedgraph[contig1][contig2][0],contig2,max(Counts)))
                    if final:
                        network.write("{0}  {1} {2} {3}\n".format(contig1,\
                        ,contig2,max(Counts)))
        return
    def tryopen(self,filename,header,filetype):
        '''Looks for a file, if its not there, it makes the file, if
        it is there then it returns the opened file. Remember to close the file if you call this
        function or use: with tryopen(stuff) as morestuff.'''
        import os
        try:
            if not os.path.isfile(filename+filetype):
                with open(filename+filetype,'w') as newfile:
                    newfile.write(header)
            elif os.path.isfile(filename+filetype):
                return open(filename+filetype,'a+')
            return open(filename+filetype,'a+')
        except:
            print "Either could not create or open the file"
            raise
        
    def stagestosif(self):
        self.graphtosif(self.totalgraph,"Initial",dirty=True)
        self.graphtosif(self.graph,"Threshold")
        #self.graphtosif(self.finalgraph,"Cov_Links",final=True)
        return

    def makegraph(self,connections=None,cleanedlinks=None,threshold=5,total=False):
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
            if connections[contig1]==None:
                pass
            else:
                for contig2 in connections[contig1]:
                    Counts=self.readCounts(contig1,contig2,clean=True)
                    #x==max(Counts) covers case of competing link orientations between two of the same
                    #contigs - default is to the the one with most in that orientation.
                    index=[i for i,x in enumerate(Counts) if (x==max(Counts)) and (x>=threshold)]
                    if index==[]:
                        pass
                    elif len(index)>0:
                        for i in index:
                            if self.linkorientation(i)!=None:
                                OrientedGraph[contig1][contig2]=[self.linkorientation(i)]
                        Graph[contig1]=Graph[contig1]+[contig2]
        self.graph=Graph
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
            Scaffolds["Scaffold"+str(len(Scaffolds))]={tig:[0,0,0,0,0]}
            del Graph[tig]
            del OrientedGraph[tig]
        #print Scaffolds
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
        #print paths
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
            done=False
            if edge not in doneedges:
                path=[edge]
                doneedges|=set([edge]) #Consider current tuple
                doneedges|=set([edge[::-1]]) #Consider the reversed tuple
                if len(path)==1: #Force a orientation consideration on starting pair
                    counts=self.readCounts(edge[0],edge[1],clean=True)
                    index=counts.index(max(counts))
                    orientation=self.linkorientation(index)
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
                        print edge
                        done=True
                paths.append(path)
        #print paths
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
        for edge in posedges:
            if front:
                #Note path[-1][1]==edge[0] by how the possible edges were decided
                counts=self.readCounts(path[-1][0],edge[0],clean=True)
                ori=self.linkorientation(counts.index(max(counts)))
            elif front==False:
                #path[0][0]==edge[1] by how possible edges are evaluated
                counts=self.readCounts(edge[0],path[0][0],clean=True) #Reorder ensures orientation for later
                ori=self.linkorientation(counts.index(max(counts)))
            else:
                pass
            if (front and ori==(0,1)) or (front==False and ori==(0,1)):
                maxcounts.append((max(counts),edge))
        #maxcounts.append((max(self.readCounts(path,edge),
        if maxcounts!=[]:
            best=max(maxcounts)
        if type(best)!=tuple:
            print best, "This is the path for typeerror"
            #print path,"This is the Path for typeerror" ###PAY ATTENTION TO THIS ERROR - IT COULD BE QUITE IMPORTANT
            return None
        print maxcounts, "Here I am bugfixing"
        for count,tup in maxcounts:
            if float(count)/best[0]>=threshold and (count,tup)!=best:
                print best, "Compared to",(count,tup)
                nextcheck=False #Don't do the next check - failed first test
                test1=False #Failed first threshold
                #Loops over pssible edges and use the second metric to evaluate acceptability
        #Won't include this in decision making for the moment
        if nextcheck:
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
            else:
                test2=False
        elif not nextcheck:
            test2=False
        #Now for coverage based method
        #merge in above loop for 2nd sspace later
        for edge in posedges:
            covtests=[]
            if front:
                covtests+=[self.metabatprobtest(path[-1][0],edge[0])]
            elif front==False:
                covtests+=[self.metabatprobtest(path[0][0],edge[0])]
        accept=[True  if x=='Accept'  else False for x in covtests]
        reject=[True if x=='Reject' else False for x in covtests]
        print "These are the accept values",accept
        print "These are the reject values", reject
        if test1 and test2: #Pass both SSPACE tests
            return best[1]
            #print "SSPACE-like accept"
        else:
            #print "SSPACE-like rejection"
            pass
        ###Now need to check if edge satisfies coverage - red false +ve
        #If one is rejection Could recalc for the next biggest no. links
        #Could increase true +ve and dec false -ve - might be a mistake to include
        #Could also inc false +ve
        #print path,"This is the Path for too many competing options"
        return best[1]
        
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
            N_rpbam+=[len(self.getlinks(tig1,tig2,bam=key+".bam"))] #Number of links per bam library
            cover+=[min(self.tiglen[tig2],seqspaces[key][0]-tempgap)]
            #Using arithmetic mean here.
        mean=sum([x/float(N_rpbam[i]) for i,x in enumerate(cover)])/3 #Average seqspace per link over bams
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
            if x==None or y==None:
                pass
            elif fuser==None:
                fuser=self.dubtuple(y,x)
            else:
                fuser=self.dubtuple(fuser,x)
        return fuser
                
    def dubtuple(self,tup1,tup2):
        ''' Takes in a pair of tuples and then decides how to join them
        on a shared edge. Eg, join (B,A) and (D,B) => (D,A,B)'''
        #print "This is the first tuple for the bug fix", tup1
        #print "This is the second tuple for the bug fix", tup2
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
        flipsequence=False
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
            #Should give unique insert mean and insert std
            meaninsert[name]=[float(x[1]) for x in insertfile if name+".bam" in x]
            stdinsert[name]=[float(x[2]) for x in insertfile if name+".bam" in x]
            ####Need to mod this list compre to get b
        cleanedlinks=[x for x in linksfile if self.lowerdist(x)<\
        (meaninsert[x[bamnameindex].strip('.bam')][0]+cutoff*stdinsert[x[bamnameindex].strip('.bam')][0])]
        for x in cleanedlinks:
            print self.lowerdist(x),
        #print meaninsert
        #print stdinsert
        for name in bamNames:
            print (meaninsert[name][0]+cutoff*stdinsert[name][0]), "This was the cutoff"
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
        self.writesometigs()
        self.makescaffolds(Graph,OrientedGraph) #Makes set of scaffolds
        self.printscaffolds()
        self.gapprint()
        self.totalgraph=makegraph(connections=self.completecheck(self.links),self.links,threshold=0)
        self.stagestosif()

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
        if cleanedlinks==None:
            cleanedlinks=self.cleanedlinks
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
            c1=int(self.tiglen[contig1])#length of contig 1
            c2=int(self.tiglen[contig2])#Length of contig 2
            cmin=min(c1,c2) #Minimum length
            cmax=max(c1,c2) #Maximum length
            if len(observations)==0:
                print "There has been an error"
                return 1
            distance=self.MLsearch(cmin,cmax,r,mu[0],sigma[0],observations)
            #print "Tig 1:",contig1,"Tig2",contig2
            print "The ML distance is {0}".format(distance)
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
        print "This is the sum of the observed distance", obsval
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
    
    def CoverageCheck(self,contig1,contig2):
        from scipy.special import erf
        from scipy.stats import norm
        from scipy.constants import pi
        from math import exp
        #Approximate the coverage as a normal distribution with mean mu
        #and stdDev sigma.
        #Then the probability that two given contigs are not from the same genome
        #based on covrage is the |norm.cdf(x)-
        #Cumulative Density function for the normal distribution norm.cdf(x)
        return


    def metabatprobtest(self,contig1,contig2,lower=0.01,upper=0.95):
        pairs=self.gettigcov(contig1,contig2)
        prob=1
        for i,pair in enumerate(pairs):
            repack=zip(*pair) 
            #print repack
            prob*=self.distancecalc(repack[0],repack[1])
        prob=prob**(1/float(i+1)) #Geometric mean of probabilities
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
        print "CAN YOU HEAR ME UP THERE"
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
        
    def networkvis(self):
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
            if i==max(means)-1:
                print "Leaving the loop soon"
            change=abs(norm.pdf(i,loc=means[0],scale=sqrt(var[0]))-norm.pdf(i,loc=means[1],scale=sqrt(var[1])))
            dist+=change
            i+=1
        #print i
        return dist/float(2)
        
    def writesometigs(self):
        with open('ExtractionList.fna','w') as extract:
            for tig in self.contignames[len(self.contignames)/2:]:
                extract.write(">{0}\n".format(tig))
        return
