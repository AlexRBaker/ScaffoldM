#!/usr/bin/env python
###############################################################################
#                                                                             #
#    compare.py                                                               #
#                                                                             #
#    Makes comparisons between ScaffoldM and SSPACE for improvement           #
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
import os
import argparse
#from dataloader import DataLoader
#from dataparser import DataParser
#from scaffold import Scaffold
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
###############################################################################
#First step - parse input to make library file for sspace_basic
#Library name readsfile1 readsfile2 insert size tolerated error, read orientation
#of form Lib1 file.1.1.fasta file.1.2.fasta 400 0.25 FR
# the above indicates that library 1 has first set of reads in file1.1.1 and second set of pair
#in file.1.2 and insert size of 400 and is willing to tolerate a 25% error (100), and reads
#map F(------>)R(<---------) onto the contigs.
def makelibrary(filename,libnames,pairedend1, pairedend2, insertsize,error,orientation,tab=False):
    import os
    import sys
    if not os.path.isfile(filename):
        scaffile=open(filename+".txt",'w')
        scaffile.close()
    with open(filename+".txt",'a+') as library:
        for i,libname in enumerate(libnames):
            if tab==False:
                library.write("{0} bwa {1} {2} {3} {4} {5}\n".\
                format(libname,pairedend1[i]+".fasta", pairedend2[i]+".fasta",\
                insertsize[i],error[i],orientation[i]))
            else:
                library.write("{0} TAB {1} {2} {3} {4}\n".\
                format(libname,TABfile+".tab",\
                insertsize[i],error[i],orientation[i]))
                

def splitter(interleavedreads):
    import sys
    try:
        with open(interleavedreads+".fna",'r') as reads:
            head=reads.readline()
            if not head.startswith('>'):
                raise TypeError("Not a FASTA file:")
            reads.seek(0)
            firstreads=[]
            secondreads=[]
            firstread=False
            secondread=False
            for line in reads:
                if line.startswith('>'):
                    firstread=False
                    secondread=False
                    if ".1" in line:
                        firstread=True
                    elif ".2" in line:
                        secondread=True
                    if firstread:
                        firstreads.append(line)
                    elif secondread:
                        secondreads.append(line)
                elif not line.startswith('>'):
                    if firstread:
                        firstreads.append(line)
                    elif secondread:
                        secondreads.append(line)
            read1=''.join(firstreads)
            read2=''.join(secondreads)
        with open(interleavedreads+"_1.fasta",'a+') as reads:
            reads.write(read1)
        with open(interleavedreads+"_2.fasta",'a+') as reads:
            reads.write(read2)
        return interleavedreads+"_1",interleavedreads+"_2"
    except:
        print "Error opening file:", interleavedreads, sys.exc_info()[0]
        raise
        
def chunker(string,chunksize,end):
    ''' Creates chunks from string and appends an end term
    Note this return generator - should be iterated over'''
    try:
        stringmod=string.translate(None,end)
        for i in xrange(0,len(stringmod),chunksize):
            if len(stringmod)>=i+chunksize:
                yield stringmod[i:i+chunksize]+end
            else:
                yield stringmod[i:i+chunksize]
    except TypeError:
        print "end must be concatenable to string, \
        intended use is type(str) for both"

def getfastalen(fastaname):
    try:
        with open(fastaname+".fasta",'r+') as fasta:
            fasta.readline()
            linelen=len(fasta.readline().strip('\n'))
            fasta.seek(0)
            for i, l in enumerate(fasta):
                pass
            return (i)*linelen #lose 1 due to header_assuming one  header for this type
    except IOError:
            with open(fastaname,'r+') as fasta:
                fasta.readline()
                linelen=len(fasta.readline().strip('\n'))
                fasta.seek(0)
                for i, l in enumerate(fasta):
                    pass
                return (i)*linelen #lose 1 due to header_assuming one  header for this type
        
        
def slicer(slices,filename):
    ''' slices are of form start, end, reps, orientation. This function will take those
    slices out of a specified fasta file and then print them as contigs'''
    import os
    import sys
    slices=[int(ele) for ele in slices]
    #print slices
    try:
        with open(filename,'r+') as Genome:
            linelen=len(Genome.readlines()[3].strip("\n"))
            start=min([slices[i] for i in range(len(slices)) if i%4==0])
            ends=max(slices) #assumes no reps or orientation greater than max contig position - reasonable
            seqslice=[]
            Genome.seek(0)
            header=Genome.readline().strip('>').rstrip('\n')
            refhead=header[0:min(24,len(header))]
            sequence=''.join([line.translate(None,"\n") for line in Genome.readlines() if not line.startswith('>')])
            for i in range(0,len(slices),4):
                if int(slices[i+3])==-1:
                    seqslice.append(reversecompliment(sequence[int(slices[i]):int(slices[i+1])]*int(slices[i+2])))
                elif int(slices[i+3])==1:
                    seqslice.append(sequence[int(slices[i]):int(slices[i+1])]*int(slices[i+2]))
    except IOError: #If not openable try for file in folder one level up
        with open('..'+os.sep+filename,'r+') as Genome:
            linelen=len(Genome.readlines()[3].strip("\n"))
            start=min([slices[i] for i in range(len(slices)) if i%4==0])
            ends=max(slices) #assumes no reps or orientation greater than max contig position - reasonable
            seqslice=[]
            Genome.seek(0)
            header=Genome.readline().strip('>').rstrip('\n')
            refhead=header[0:min(24,len(header))]
            sequence=''.join([line.translate(None,"\n") for line in Genome.readlines() if not line.startswith('>')])
            for i in range(0,len(slices),4):
                seqslice.append(sequence[int(slices[i]):int(slices[i+1])][::int(slices[i+2])]*int(slices[i+3]))
    parts=filename.split(os.sep)
    fileend=parts[-1].split(".fasta")[0]
    slicefilename=fileend+"slices.fna"
    completefilename="{0}S:{1}_E:{2}".format(fileend,start,ends)+"complete.fna"
    if not os.path.isfile(slicefilename):
        tigfile=open(slicefilename,'w')
        tigfile=tigfile.close()
    if not os.path.isfile(completefilename):
        tigfile=open(completefilename,'w')
        tigfile=tigfile.close()
    with open(slicefilename,'a+') as tigfile:
        for i,seq in enumerate(seqslice):
            tigname=refhead+"contig"+str(i+1)+"|"+header[int(min(24,len(header))):]
            tigfile.write(">{0}, S:{1}:E:{2}:R:{3}:OR:{4}\n".format(tigname,slices[i*4],slices[i*4+1],slices[i*4+3],slices[i*4+2]))
            for chunk in chunker(seq,linelen,"\n"):
                tigfile.write(chunk)
            tigfile.write('\n')
    with open(completefilename,'a+') as tigfile:
        tigname=refhead+"complete|"+header[min(24,len(header)):]
        tigfile.write(">{0}, S:{1}:E:{2}\n".format(tigname,start,ends))
        for chunk in chunker(sequence[start:ends],linelen,"\n"):
            tigfile.write(chunk)
        tigfile.write('\n') 
    return slicefilename,completefilename
    
def reversecompliment(seq):
    compliment={'A':'T','G':'C','T':'A','C':'G','a':'t','g':'c','t':'a','c':'g'}
    try:
        newseq="".join([compliment[char] for char in reversed(seq)])
        return newseq
    except:
        print "This sequence has illegal characters"
        raise

def randcuts(gap,seqlen,noslices=10,rep=False,ori=False,steps=100,gapvar=False,randgaps=True):
    import random
    starts=[]
    ends=[]
    orientation=[]
    reps=[]
    m=[1]*90+10*[-1]
    replist=[1,2,3,4]
    starts.append(random.randint(0,seqlen/4))
    upperlen=(seqlen-starts[0])//noslices
    for i in range(noslices):
        ends.append(random.randint(starts[i]+steps,starts[i]+upperlen+steps))
        if i!=(noslices-1):
            if randgaps:
                starts.append(ends[i]+int(random.gauss(gap,gap/5)))
            else:
                starts.append(ends[i]+gap)
        if ori:
            orientation.append(random.sample(m,1)[0])
        else:
            orientation.append(1)
        if rep:
            reps.append(random.sample(replist,1)[0])
        else:
            reps.append(1)
    weave=zip(starts,ends,reps,orientation)
    return [int(zipdat) for zipped in weave for zipdat in zipped]

def makereads(readnumber,readlength,\
               meaninsert,stdinsert,filename,readmaker="metasim"):
    if readmaker=="metasim":
        statoscom=("/home/baker/Packages/metasim/MetaSim cmd -r {0} -m \
            -g /home/baker/Packages/metasim/examples/errormodel-100bp.mconf \
            -2 /home/baker/Packages/metasim/examples/errormodel-100bp.mconf \
            --empirical-pe-probability {1} --clones-mean {2} --clones-param2 {3} \
            {4}").format(readnumber,readlength,meaninsert,stdinsert,\
           filename)
    elif readmaker=="gemsim":
        pass
    os.system(statoscom)
    
def makereadswrap(readnumber,readlength,\
               meaninsert,stdinsert,filename,readmaker="metasim"):
    makereads(readnumber,readlength,\
               meaninsert,stdinsert,filename,readmaker="metasim")
    if readmaker=='metasim':
        readname=filename.split(".fna")[0]+"-Empirical.fna"
        print "This is the read name", readname
        setreads(readname)
    
def setreads(filename,clean=True):
    import os
    tempname=filename.split(".fna")[0]+"redone.fna"
    print tempname, "This is the temporary name"
    os.rename(filename,tempname)
    #Make a new file for writing output
    with open(filename,'w') as make:
        pass
    with open(tempname) as reads:
        seq=[]
        Ind=False
        with open(filename,'a+') as oldfile:
            for line in reads:
                if line.startswith('>'):
                    Ind=False
                    if seq!=[]:
                        oldfile.write("".join(seq)+"\n")
                    if clean:
                        oldfile.write(line.split(" ")[0]+"\n")
                    else:
                        oldfile.write(line)
                    seq=[]
                else:
                    seq.append(line.translate(None,'\n'))
            oldfile.write("".join(seq)+"\n")
            seq=[]
    os.remove(tempname)
    return
                
def makesummary(cuts):
    ''' Takes the cuts used for genome slicing and makes
    a summary file detailing changes'''
    return
    
def scaffoldparse(scaffoldloc,scaffoldloc2,truegaps,trueorientations):
    '''Given two scaffold locations of predefined structure extract information from it and compare'''
    
    return
def Falsejoins(type1errors):
    with open("ScafMFalseJoins.fasta",'a+') as Joins:
        Joins.write("Mistake1,Mistake2\n")
        for tup in type1errors:
            ind=False
            for tig in tup:
                if ind==False:
                    Joins.write(str(tig)+",")
                elif ind:
                    Joins.write(str(tig)+"\n")
                ind=True
def trackdecisions(Truepos,Falsepos,falseneg,notigs=False,N_joins=False):
    import os
    if not os.path.isfile("../Results.txt"):
        with open("../Results.txt",'w') as data:
            data.write("{0},{1},{2}\n".format("TruePositive","FalsePositive","FalseNegative"))
    with open("../Results.txt",'a+') as data:
        data.write("{0},{1},{2}\n".format(len(Truepos),len(Falsepos),len(falseneg)))
    return
    
def Visualise(scaffoldnames,gaps,contigloc,covplot=False):
    import matplotlib.pyplot as plt
    import numpy as np
    import scipy as sp
    import pandas as pd
    scaffoldMgap,scaffoldMTjoins,scaffoldMFjoins,scaffoldMFNeg=validcheck(contigloc=contigloc)
    data={}
    Falsejoins(scaffoldMFjoins)
    trackdecisions(scaffoldMTjoins,scaffoldMFjoins,scaffoldMFNeg)
    sortMgap=sorted(scaffoldMgap,key=lambda x: min(x[0]))
    pairs=[(min(x[0])-1,x[1]) for x in sortMgap]
    minscafind,gap2=zip(*pairs)
    pairedgaps=[gap for i,gap in enumerate(gaps) if i in minscafind]
    for scaffold in scaffoldnames:
        data[scaffold]=[contiglen(scaffold)]
        data[scaffold]+=[[NXcalc(X*0.1,data[scaffold][0]) for X in range(1,11)]]
        data[scaffold]+=[[len(data[scaffold])-1]]
    #print data
    standardplot(gap2,pairedgaps,"The predicted gapsize(nt)","The actual gapsize(nt)","","ScaffoldMVTrueGap")
    multiplot([[X*0.1 for X in range(1,11)] for i in range(0,len(scaffoldnames))],[data[scaffold][-2] \
    for scaffold in scaffoldnames],"X","NX Value for the scaffold",\
    "The NX metric for various scaffolds", ["Contigs","ScaffoldM","SSPACE"],"N50Metric_Scaffolds")
    if covplot:
        for scaffold in scaffoldnames:
            plotcoverage(scaffold)
    return

def contigmap(evidencefile,coveragefile='covs.tsv'):
    Scaffolds={}
    with open(evidencefile,'r+') as SSPACE:
        for line in SSPACE:
            if line.startswith(">"):
                parts=line.split('|')
                scaffold=parts[0].strip('>')
                Scaffolds[scaffold]=[]
            else:
                parts=line.split('|')[0] #tig name - in form f_tign
                if parts!='\n':
                    if parts.startswith('r'):
                        parts='f'+parts[1:]
                    Scaffolds[scaffold]+=[parts]
    with open(coveragefile,'r+') as covs:
        covs.readline() #Move paste header
        orderedtigs=[]
        for line in covs:
            orderedtigs.append(line.split('\t')[0]) #Contig name
    N_tigs=len(orderedtigs) #Number of contigs
    #map f_tigi to orderedtigs[i] in dictionary
    Swapdict={"{0}{1}".format('f_tig',i):orderedtigs[i-1] for i in range(1,N_tigs+1)}
    Mapped={scaffold:[Swapdict[contig] for contig in contigs] for scaffold,contigs in Scaffolds.iteritems()}
    print Mapped
    return Mapped
    
def scaffoldtoedges(mapped):
    SS_data=np.zeros((1,3))
    for i,(scaffold,contigs) in enumerate(mapped.iteritems()):
        for j,contig in enumerate(contigs):
            if j<len(contigs)-1:
                SS_data=np.vstack((SS_data,np.array([contig,'(0,1)',contigs[j+1]])))
    SS_data=SS_data[1:,:] #Remove initial dummy row
    np.savetxt('SSPACE_Edges.txt',SS_data,fmt='%s',delimiter='\t',newline='\n', header='Edge1\trel\tEdge2\n')
    return SS_data
    
def graphtosif(self,graph,graphname,removed=False):
    #print graphname, "This is the supposed graph being parsed"
    #print "\n",graph
    done=set([])
    with tryopen("{0}{1}".format(graphname,"_links"),"Contig1\tRelationship\tContig2\n",".txt",True) as network2:
        for contig1,connected in graph.iteritems():
            for contig2 in connected:
                if (contig1,contig2) not in done and (contig2,contig1) not in done:
                    network2.write("{0}\t{1}\t{3}\n".format())
                    done|=set([(contig1,contig2)]) #Add current pair
                    done|=set([(contig1,contig2)[::-1]]) #Reverse of current pair
    with tryopen("{0}{1}".format(graphname,"_contigs"),"Contig\n",".txt") as contigs:
        for contig in graph:
            contigs.write("{0}".format(contig))
    return
    
def addcol(filename,column_s,header,d='\t'):
    '''Takes a text file containing tab separated columns and adds tab-separated columns
    to the end. This is primarily for updating the .txt files using in Cytoscape with additional
    info such as bin allocation etc. Loads whole file into memory.'''
    with tryopen(filename,'','') as oldfile:
        olf=oldfile.readlines()
        Newfile=[]
        for i,line in enumerate(olf):
            if i==0:
                processedline=[x.rstrip('\n') for x in line.split('\t')]
                processedline+=[head for head in header]
                processedline[-1]=processedline[-1]+"\n"
                Newfile+=processedline
            else:
                processedline=[x.rstrip('\n') for x in line.split('\t')]
                processedline+=[column[i] for column in column_s]
                processedline[-1]=processedline[-1]+"\n" #Add endline
                Newfile+=processedline
    with tryopen(filename,'','',True) as newfile:
        for line in Newfile:
            newfile.write(("{0}".format(d)).join(line))
    return
        
def sspaceconvert(evidencefile,coveragefile='covs.tsv'):
    Mapped=scaffoldtoedges(contigmap(evidencefile,coveragefile='covs.tsv'))
    return Mapped
    
def plotcoverage(name):
    return
#Simply comparison - one scaffolder and preprocessed dataset
def standardplot(x,y,xname,yname,title,saveloc,log=False):
    import matplotlib.pyplot as plt
    plt.gca().set_color_cycle(['blue', 'black'])
    plt.plot(x,y,'o')
    plt.plot(y,y,'-')
    plt.xlabel(xname)
    plt.ylabel(yname)
    plt.title(title)
    plt.savefig('./graphs/'+saveloc+'.png',bbox_inches='tight')
    plt.close()
    return
#More complicated comparisions - Likely between multiple scaffolders

def multiplot(x,y,xname,yname,title,legend,saveloc,log=False):
    '''For plotting lists of lists for x and y, along with an appropiate legend'''
    import matplotlib.pyplot as plt
    plt.gca().set_color_cycle(['red', 'blue', 'black'])
    print x, "This is the x variable"
    print y, "This is the y variable"
    for i,xdat in enumerate(x):
        plt.plot(x[i],y[i])
    plt.xlabel(xname)
    plt.ylabel(yname)
    plt.title(title)
    plt.xticks([0.1*X for X in range(0,11)])
    plt.axis([min(xv for xval in x for xv in xval),max(xv for xval in x for xv in xval),0,1.05*max(yv for yval in y for yv in yval)])
    leg=plt.legend(legend, loc='upper right',title='Scaffolder')
    leg.get_frame().set_alpha(0)
    plt.savefig('./graphs/'+saveloc+'.png',bbox_inches='tight')
    plt.close()

def multibar(x,y,xlab,ylab,title,saveloc,legend):
    #Sourced from :http://matplotlib.org/examples/api/barchart_demo.html
    #To be modified heavily later
    import numpy as np
    import matplotlib.pyplot as plt
    N = 5
    menMeans = (20, 35, 30, 35, 27)
    menStd =   (2, 3, 4, 1, 2)
    ind = np.arange(N)  # the x locations for the groups
    width = 0.35       # the width of the bars
    fig, ax = plt.subplots()
    rects1 = ax.bar(ind, menMeans, width, color='r', yerr=menStd)
    womenMeans = (25, 32, 34, 20, 25)
    womenStd =   (3, 5, 2, 3, 3)
    rects2 = ax.bar(ind+width, womenMeans, width, color='y', yerr=womenStd)
    # add some text for labels, title and axes ticks
    ax.set_ylabel('Scores')
    ax.set_title('Scores by group and gender')
    ax.set_xticks(ind+width)
    ax.set_xticklabels( ('G1', 'G2', 'G3', 'G4', 'G5') )
    ax.legend( (rects1[0], rects2[0]), ('Men', 'Women') )
    def autolabel(rects):
        # attach some text labels
        for rect in rects:
            height = rect.get_height()
            ax.text(rect.get_x()+rect.get_width()/2., 1.05*height, '%d'%int(height),
                    ha='center', va='bottom')
    autolabel(rects1)
    autolabel(rects2)
    plt.show()
    
def validcheck(gapdataloc="Gapdata.txt",contigloc='MG1655refslices.fna'):
    '''Uses Gapdata.txt to compare observed scaffolds to known scaffold'''  
    validgaps=[]
    truejoins=[]
    falsejoins=[]
    falsenegs=[]
    truepairs=[(i,i+1) for i in range(1,len(contiglen(contigloc)))]
    with open(gapdataloc) as gaps:
        gaps.readline() #Move past header
        for line in gaps:
            compsplit=line.split(",")
            tig1=compsplit[0].rstrip("|").split("|")[-1]
            tig2=compsplit[1].rstrip("|").split("|")[-1]
            #print compsplit[-1]
            gap=int(compsplit[-1])
            tig1ind=int(tig1.split('g')[-1])
            tig2ind=int(tig2.split('g')[-1])
            if abs(tig1ind-tig2ind)==1:
                validgaps.append(((tig1ind,tig2ind),gap))
                truejoins.append((tig1,tig2))
            else:
                falsejoins.append((tig1,tig2))
    print "THESE ARE THE FALSENEGATIVES"
    falsenegs=[x for x in truepairs if x not in zip(*validgaps)[0]] 
    print falsenegs, "The Reject|True"
    return [validgaps,truejoins,falsejoins,falsenegs]

def NXcalc(X,tiglengths):
    tot=sum(tiglengths)
    S_tig=sorted(tiglengths)
    N50=0
    runtot=0
    if X==0:
        return S_tig[-1]
    for i in range(len(S_tig)):
        N50=S_tig[-(i+1)]
        runtot+=N50
        if runtot>tot*X:
            return N50
            
def parsetsv(filename='links.tsv',delim=',',header=False):
    parsed=[]
    with open(filename) as tsv:
        if header==False:
            tsv.readline()
        for line in tsv:
            parsed.append(line.split(delim))
    return parsed
    
def getlinks(contig1,contig2,filename='links.tsv'):
    links=parsetsv(filename,delim='\t')
    flags=[]
    for link in links[1:]:
        tig1=0
        tig2=0
        for col in link:
            #Need to change for simple contig names since I stripped the full name earlier
            if contig1 in col:
                tig1=1
            if contig2 in col:
                tig2=1
        if tig1+tig2==2:
            flags.append(link)
    #print flags
    return flags

def linkdist(onelink):
    '''distance from relevant edge for each contig in the link'''
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
        #Returns tuple of distances and contig names
    return  ((dist1,onelink[0]),(dist2,onelink[1]))

def extracttigs(filename='ScafMFalseJoins.fasta',outfile='mislinkseq.fasta',contigloc='MG1655refslices.fna'):
    missjoins=parsetsv(filename,delim=',',header=False)
    missjoins=[[y.rstrip('\n') for y in x] for x in missjoins]
    print "THe faulty contig pair",missjoins
    print "You are right before the loop"
    for line in missjoins:
        #print line
        #print line[0],
        tigs=getlinks(line[0],line[1]) #Extract missjoins
        #print "This is the links",tigs
        orientation=(tigs[0][4],tigs[0][7]) #Assuming only one orientation present - risky
        #print "This is the orientation",orientation
        #Fix this assumption later
        dists=[linkdist(x) for x in tigs] #Get links distances
        #print "This is the distance",dists
        #Works with assumption that joins are always in same order - they are
        maximum=[max(x) for x in zip(*dists)] #Unzips tuple into list for each contig
        #print "This is the maximum",maximum
        #GEts maximum distance from edge
        W_faultylink(maximum[0][0],maximum[1][0],maximum[0][1],maximum[1][1],orientation,outfile,contigloc)

def W_faultylink(distance1,distance2,contig1,contig2,orientation,filename='mislinkseq.fasta',contigloc='MG1655refslices.fna'):
    import os
    import sys
    print "Did I make it this far"
    try:
        if not os.path.isfile("../{0}".format(filename)):
            with open("../{0}".format(filename),'w') as test:
                pass
        with open("../{0}".format(filename),'a+') as mislink:
            mislink.write(">{0}|Distance:{1}bp_from_edge\n".format(contig1,distance1))
            for chunk in chunker(cut(extractcontigs(contig1,contigloc,header=False).translate(None,'\n'),distance1,orientation,0),70,'\n'):
                mislink.write(chunk)
            mislink.write('\n')
            mislink.write(">{0}|Distance:{1}bp_from_edge\n".format(contig2,distance2))
            for chunk in chunker(cut(extractcontigs(contig2,contigloc,header=False).translate(None,'\n'),distance2,orientation,1),70,'\n'):
                mislink.write(chunk)
            mislink.write('\n')
    except:
        print "Errors opening file or running stuff"
        raise ValueError
    
def cut(seq,slicesize,orientation,tigpairno):
    or1=int(orientation[0])
    or2=int(orientation[1])
    if tigpairno==0:
        if or1==1:
            return seq[:slicesize]
        else:
            return seq[-slicesize:]
    elif tigpairno==1:
        if or2==1:
            return seq[:slicesize]
        else:
            return seq[-slicesize:]

def extractcontigs(contigname,contigloc,header=True):
    '''Just assigns contigs file via contigloc.
    Temporary just for use when making scaffold
    Will extract the text for that contig'''
    import sys
    try:
        with open(contigloc,'r') as Contigs:
            head=Contigs.readline()
            if not head.startswith('>'):
                raise TypeError("Not a FASTA file:")
            Contigs.seek(0)
            title=head[1:].rstrip() ##Strips whitespace and >
            record=0
            contigseq=[]
            for line in Contigs:
                if line.startswith('>') and line.find(contigname)>=0:
                    record=1
                    if header:
                        contigseq.append(line)
                elif line.startswith('>') and contigname not in line:
                    record=0
                elif record==1:
                    contigseq.append(line)
                else:
                    pass
            seq=''.join(contigseq)
            return seq
    except:
        print "Error opening file:", contigloc,sys.exc_info()[0]
        raise        
    
        
            
def contiglen(contigloc):
    '''Just goes through multi-fasta file, and works out sequence length for eahc entry'''
    import sys
    try:
        lengths=[]
        with open(contigloc,'r+') as Contigs:
            head=Contigs.readline()
            if not head.startswith('>'):
                raise TypeError("Not a FASTA file:")
            Contigs.seek(0)
            contigseq=[]
            for line in Contigs:
                if line.startswith('>'):
                    if contigseq!=[]:
                        lengths.append(len("".join(contigseq)))
                    contigseq=[]
                else:
                    contigseq.append(line.translate(None,"\n"))
            lengths.append(len("".join(contigseq)))
        return lengths
    except:
        print "Error opening file:", contigloc,sys.exc_info()[0]
        raise      
def writeout(data):
    return
def postprocess(sifloc,trueloc,final=False):
    with tryopen(trueloc,'','.txt') as correct:
        correct.readline() #Move past header
        Truescaf={}
        for line in correct.readlines():
            curline=line.split('\t')
            if curline[0] not in Truescaf:
                Truescaf[curline[0]]=[]
            if curline[1].rstrip('\n') not in Truescaf:
                Truescaf[curline[1].rstrip('\n')]=[]
            if curline[1].rstrip('\n') not in Truescaf[curline[0]]:
                Truescaf[curline[0]]+=[curline[1].rstrip('\n')]
            if curline[0] not in Truescaf[curline[1].rstrip('\n')]:
                Truescaf[curline[1].rstrip('\n')]+=[curline[0]]
    Newdata=[]
    with tryopen(sifloc,'','.txt') as olddata:
        if not final:
            header=olddata.readline().translate(None,'\n')+"\tTrueEdge\tDecision\n"
        else:
            header=olddata.readline().translate(None,'\n')+"\tTrueEdge\n"
        for line in olddata.readlines():
            curline=line.split('\t')
            tig1=curline[0]
            tig2=curline[2]
            if tig1 in Truescaf:
                if tig2 in Truescaf[tig1]:
                    TrueEdge="True"
                else:
                    TrueEdge="False"
            else:
                TrueEdge="False"
            if not final:
                remove=curline[4].rstrip('\n')
                #print remove
                if remove=="True" and TrueEdge=="True":
                    Decision="FalseNeg"
                elif remove=="True" and TrueEdge=="False":
                    Decision="TrueNeg"
                elif remove=="False" and TrueEdge=="True":
                    Decision="TruePos"
                elif remove=="False" and TrueEdge=="False":
                    Decision="FalsePos"
                Newdata+=[[x.rstrip('\n') for x in curline]+[TrueEdge]+[Decision]]
            else:
                Newdata+=[[x.rstrip('\n') for x in curline]+[TrueEdge]]
            
    with tryopen(sifloc,header,'.txt',True) as final:
        for line in Newdata:
            final.write("{0}\n".format("\t".join(line)))
    return

def totprocess(sifloc1,sifloc2,sifloc3,trueloc="TrueEdges"):
    postprocess(sifloc1,trueloc,False)
    postprocess(sifloc2,trueloc,False)
    postprocess(sifloc3,trueloc,True)
    return
    
def maketrueedges():
    '''Assumes that contigslices file is both ordered by position withiin each species
    and by species'''
    TrueEdges=[]
    with tryopen("covs",'','.tsv') as covs:
        i=1
        prevtig=False
        curtig=False
        donetigs=[]
        TrueEdge=[]
        covs.readline() #Move past header
        for line in covs:
            prevtig=curtig
            curtig=line.split('\t')[0]
            #print "This is the current Contig", curtig
            #print "This is the current line", line.split('\t')[0]
            tignumber='{0}{1}'.format('tig',i)
            if tignumber in curtig:
                i+=1
                if tignumber in donetigs:
                    donetigs=[]
                    i=2
                elif prevtig!=False:
                    TrueEdge+=[(prevtig,curtig)]
                    donetigs+=[tignumber]
            else:
                i=2
                donetigs=[]
                #print tignumber
                #print curtig
            

    with tryopen("TrueEdges",'contig1\tcontig2\n','.txt',True) as Edge:
        for edge in TrueEdge:
            #print edge
            Edge.write("{0}\t{1}\n".format(edge[0],edge[1]))      
    return
def binstotxt(checkdir,fileend='.fa'):
    import os
    files=os.listdir(checkdir)
    bins=[File for File in files if File.endswith('.fa')] #Should get those with .fa in name
    binpaths=[os.path.join(checkdir,binfile) for binfile in bins]
    contigbinmap={}
    for i,binloc in enumerate(binpaths):
        contigbinmap[bins[i].strip(fileend)]=getfastaheaders(binloc) #Maps the set of contigs in
        #bin file to that bin in a graph
    return contigbinmap
    
def writelis(filename,header,rows):
    with tryopen(filename,header,'') as newfile:
        for row in rows:
            newfile.write("\t".join(row)+"\n")
    return
        
def binwrapper(writedir,filename,checkdir,fileend='.fa'):
    import os
    bins=binstotxt(checkdir,fileend)
    rows=[(item,key) for key,items in bins.iteritems() for item in items]
    columns=zip(rows)
    writelis(os.path.join(writedir,filename),"Contig\tBin\n",rows)

def getfastaheaders(filename):
    try:
        Names=[]
        with tryopen(filename,'','') as fasta:
            for line in fasta:
                if line.startswith('>'):
                    #print line
                    Names+=[line.strip('>').rstrip('\n')]
                else:
                    pass
        return Names
    except:
        raise
        
def screenbins(binsfile,contamlevel, completeness):
    
    return
def tryopen(filename,header,filetype,expunge=False):
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
        
def makeboolean(string):
        Val_T=("true",'t','1','yes')
        Val_F=("false","f",'0','no')
        try:
            if isinstance(string,bool):
                return string
            if string.lower() in Val_T:
                return True
            elif string.lower() in Val_F:
                return False
        except:
            raise TypeError("This does not appear to even be an attempt at a boolean")
if __name__ == "__main__": ###Check if arguments coming in from command line
    import matplotlib.pyplot as plt
    import numpy as np
    import scipy as sp
    import os
    import sys
    import datetime
    parser = argparse.ArgumentParser(description='Takes a Genome in Fasta format and then parse it into predefined chunks\
     These reads and gaps are then parsed to both SSPace and ScaffoldM which then have there output extracted and compared.\
     These comparisons have formed the basis for improvements on ScaffoldM.')

    parser.add_argument('-N','--name', type=str, nargs='?', \
		help='The name of the file',default='/home/baker/Documents/Geneslab/TestFunctions/ReferenceFasta/MG1655ref.fasta')
    parser.add_argument('-P','--path', type=str, nargs='?', \
	default='/home/baker/Packages/SSPACE-STANDARD-3.0/SSPACE_Standard_v3.0.pl', \
		help='The name of absolute path')
    parser.add_argument('-L','--lists', type=int, nargs='*', help='Test',default=False)
    parser.add_argument('-l','--readlength',type=int,nargs='?',
    help='The length of reads in the simulated library',default=100)
    parser.add_argument('-c','--coverage',type=int,nargs='?',
    help='Coverage in simulated library',default=30)
    parser.add_argument('-b','--bams',type=str,nargs='*',
    help='Coverage in simulated library',default="NA")
    parser.add_argument('-m','--meaninsert',type=int,nargs='?',
    help='The mean insert size of the library, this is the expected gap between two paired reads',default=300)
    parser.add_argument('-s','--stdinsert',type=int,nargs='?',
    help='The standard deviation of the insert size between paired reads',default=30)
    parser.add_argument('-g','--gap',type=int,nargs='?',
    help='The gap between contigs',default=50) 
    parser.add_argument('-lim','--linklimit',type=int,nargs='?',
    help='The number of linking reads needed to link contigs',default=5)
    parser.add_argument('-r','--ratio',type=int,nargs='?',
    help='If one contig is linked to multiple others a comparison is made between the number of linking reads. \
    this ratio is the value at which both links will be rejected to avoid false positives',default=0.7)
    parser.add_argument('-ln','--libno',type=int,nargs='?',
    help='The number of libraries expected to be found in a bamm file',default=[1])
    parser.add_argument('-ns','--nslice',type=int,nargs='?',
    help='The number of slices to use for data simulation',default=100)
    parser.add_argument('-t','--trim',type=int,nargs='?',
    help='Whether or not',default=0.2)
    parser.add_argument('-e','--error',type=int,nargs='?',
    help='Whether or not',default=0.75)
    parser.add_argument('-w','--wrapperp',type=str,nargs='?',
    help='The path to ScaffoldM wrapper',default="~/Documents/Geneslab/ScaffoldM/scaffoldm/")
    parser.add_argument('-si','--sim',type=str,nargs='?',
    help='Whether or not to simulate',default=True)
    parser.add_argument('-C','--contiglocation',type=str,nargs='?',
    help='The path to contigs',default="")
    parser.add_argument('-O','--randominversions',type=bool,nargs='?',
    help='The path to contigs',default=False)
    parser.add_argument('-rep','--rep',type=bool,nargs='?',
    help='Boolean- whether to repeat sequences',default=False)
    parser.add_argument('-Tig','--tigname',type=str,nargs='?',
    help='Name of Fasta file for mapping',default="mergedslices.fasta")
    args = parser.parse_args()
    parser.add_argument('-Rsim','--simreads',type=str,nargs='?',
    help='Name of Fasta file for mapping',default=False)
    args = parser.parse_args()
    ###Stuff for Simulation
    sim=makeboolean(args.sim) #Whether or not to simulate
    path=args.path #path to sspace
    Name=args.name #Path to reference genome
    contigloc=args.contiglocation #path to contigs
    prename=os.sep.join(Name.split(os.sep)[:-1]) #Strips last layer of path
    if len(Name.split(os.sep)[-1].split(".fasta"))>1:
        postname=Name.split(os.sep)[-1].split(".fasta")[0] #Should be name of reference file - path and file type
    elif len(Name.split(os.sep)[-1].split(".fna"))>1:
        postname=Name.split(os.sep)[-1].split(".fna")[0]
    refpath="../{0}".format(Name.split(os.sep)[-2])
    newname="../{0}/{1}".format(Name.split(os.sep)[-2],Name.split(os.sep)[-1].split(".fasta")[0]) #For moving up and into reference folder
    coverage=args.coverage # A specified amount of coverage for simulation
    readlength=args.readlength #Simulated read length
    meaninsert=args.meaninsert #Mean insert size
    stdinsert=args.stdinsert #Std deviation in insert size
    lists=args.lists #THe list of cuts to be made to the reference genome
    linklimit=args.linklimit #The lower limit to accpet a pairing - eg ignore all pairs with less than k links
    ratio=args.ratio #Ratio for SSPACE algorithm
    gap=args.gap #Mean value for simulated gap size
    libno=args.libno #THe number of libraries BAMM should search for amongst the reads
    N_slices=args.nslice #THe number of slices to make
    error=args.error #The error for SSPACE to accept read inserts
    trim=args.trim #Can't remember what this does
    ori=args.randominversions #Whether or not to randomly take the reverse compliment
    rep=args.rep  #Whether or not to randomly repeat sequence in the simulation
    wrapperpath=args.wrapperp #The path to the python wrapper
    seqlen=getfastalen(Name) #Approximate length of the genome
    bams=args.bams #name of bams
    contigname=args.tigname #name of fasta file
    simreads=makeboolean(args.simreads)
    #STuff for all comparisons
    
    if sim:
        if trim==True:
            #Trim overall length randomly for more variability
            pass
        #Turn long list of cut locations into a string suitable for os.system
        if lists==False:
            slices=randcuts(gap,seqlen,N_slices,ori=ori)
            #print slices
            slices=[str(i) for i in slices]
            readcuts=" ".join(slices)
        else:
            slices=lists
            readcuts=" ".join(lists)
        start=min([int(j) for i,j in enumerate(slices) if i%4==0]) #Starting position of cuts
        end=max([int(j) for i,j in enumerate(slices) if i%4==1]) #Starting position of cuts
        os.mkdir(postname+"cuts_S:{0}_E:{1}".format(start,end))
        os.chdir(postname+"cuts_S:{0}_E:{1}".format(start,end))
        slicename,completename=slicer(slices,Name)
        #Make reads via metasim
        if simreads:
            readnumber=coverage*(end-start)/readlength
            makereadswrap(readnumber,readlength,meaninsert,stdinsert,\
            completename)
            #splits the reads into a format suitable for
            completename=completename.split(".fna")[0]
            file1,file2=splitter(completename+"-Empirical")
        else:
            #print os.getcwd()
            file1='{0}/{1}-Empirical_1'.format(refpath,postname)
            file2='{0}/{1}-Empirical_2'.format(refpath,postname)
            completename=newname
        contigloc=slicename
        #Make libraries.txts
        makelibrary("library",['Lib1'], [file1],[file2],[meaninsert],[error],orientation=['FR'])
    else:
        contigloc=contigname
        pass
    #To separate from SIM - needs a library file
    #perl SSPACE_Basic.pl -l libraries.txt -s contigs.fasta -x 0 -m 32 -o 20 -t 0 -k 5 -a 0.70 -n 15 -p 0 -v 0 -z 0 -g 0 -T 1 -b standard_out
    mapreads=False
    if mapreads:
        print "SSPACE", contigloc, "The contig file"
        print "perl {4} -l {0} -s {1} -x 0 \
         -k {2} -a {3} -b standard_out"\
         .format("library.txt",contigloc,linklimit,ratio,path)
        os.system("perl {4} -l {0} -s {1} -x 0 \
         -k {2} -a {3} -b standard_out"\
         .format("library.txt",contigloc,linklimit,ratio,path))
        print "Onto BamM"
        if sim:
            os.system("bamm make -d {0} -i {1} --quiet".format(slicename,completename+"-Empirical.fna"))
            libno=[str(ele) for ele in libno]
            librarynumbers=' '.join(libno)
            bamname="{0}{1}{2}".format(slicename.split(".fna")[0],".",postname)+"-Empirical"
            contigname=slicename
    print "The current time: ", datetime.datetime.now().time().isoformat()
    if type(bams)!=str: #Check if default is being used
        libno=[str(ele) for ele in libno]
        librarynumbers=' '.join(libno)
        print "The Bams", ' '.join(bams)
        os.system("python {0}wrapper.py -b {1} -f {2} -n {3}".format(wrapperpath,' '.join(bams),contigname,librarynumbers))
        real=True
        if not real:
            maketrueedges()
            totprocess("Initial_links","Threshold_links","Cov_Links_links")
        #print os.getcwd()
        if os.path.isdir('./standard_out'): #Only go if SSPACE worked
            SSPACEgraph=sspaceconvert("./standard_out/standard_out.final")
            graphtosif(SSPACEgraph)
        binned=True
        if binned:
            binwrapper('.','Node_BinClass.txt','bins_raw/')
            print "You made it pass checking the nodes"
        #Separate those with high enough completeness/quality scores
        #Visualise the remaining
        #Work out how to colour based on this in cytoscape
        #Bam, done, can compare into and out of binning occurrences
        print "This is the real end now"
    elif sim: #Should only occur on defaults
        os.system("python {0}wrapper.py -b {1} -f {2} -n {3}".format(wrapperpath,bamname,contigname,librarynumbers))
        maketrueedges()
        totprocess("Initial_links","Threshold_links","Cov_Links_links")
        print "This is the real end now"
    comparisons=False
    if comparisons:
        #Make some comparisons between SSPACE and ScaffoldM
        os.mkdir('graphs')
        Visualise([slicename,"testScaffold.fasta","./standard_out/standard_out.final.scaffolds.fasta"],\
        [int(slices[i+4])-int(slices[i+1]) for i in range(0,len(slices)-4,4)],contigloc)
        if sim:
            print "You made it to the link error extractions"
            extracttigs(contigloc=slicename)
        else:
            pass

        
        
