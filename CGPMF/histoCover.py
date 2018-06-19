# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 14:14:31 2016

@author: rachael
This file will open a histo file returned from Gromacs g_wham output and analyze
the overlap of the histograms given
Should also be given list of xvg pull files and then can return a list of 
frames needed to (a) continue and (b) start new to fully cover histogram
"""
import numpy as np,re,sys,argparse

def parseTFile(tprname):
    #this takes a tprfile and returns a list of integers corresponding
    #to the runfile indices
    f = open(tprname,'r')
    lines = f.readlines()
    f.close()
    finds = []
    for line in lines:
        m = re.search('\d+',line)
        ind = m.group(0)
        finds.append(int(ind))
    return finds

def parseHFile(fname,tprname):
    #given a file name, assume it is a histo.xvg file from g_wham
    #read in data and return in the form of a 2d numpy array
    #order of columns in hist.xvg: same as order in tpr-files.dat (I'm pretty
    # sure)
    #this function also takes the tpr-files.dat fname and then reorders the
    #columns to put neighboring files together
    f = open(fname,'r')
    lines = f.readlines()
    f.close()
    cinds = parseTFile(tprname)
    #ordinds = np.argsort(cinds)
    #lines 0-11 are comments
    #create array
    histos = np.zeros([len(lines)-12,len(lines[12].split())])
    for i in range(12,len(lines)):
        spline = lines[i].split()
        for s in range(len(spline)):
            
            histos[i-12,s] = float(spline[s])
    return (histos,cinds)
       

def sortedHistos(histos,cinds):
    #assuming we have an array as returned by the parsing function, we want to 
    #turn it into a numpy array with 3 columns: [0: where is the highest point
    #of the histogram, 1: where is the lowest val it covers, 2: where is the 
    #greatest val it covers
    #then each line is a histogram
    #and we want to sort these based on column 0
    hs = np.zeros([np.shape(histos)[1]-1,3])
    for i in range(1,np.shape(histos)[1]):
        h = histos[:,i]
        dinds = np.argwhere(h)
        low = min(dinds)[0]
        high = max(dinds)[0]
        maxpoint = np.where(h==max(h))[0][0]
        hs[i-1,0] = histos[maxpoint,0]
        hs[i-1,1] = histos[low,0]
        hs[i-1,2] = histos[high,0]
    sinds = hs[:,0].argsort()
    runs = np.array(cinds)
    runs = runs[sinds]
    hs = hs[sinds]
    return (hs,runs)

def checkOverlap(h1,h2):
    #just checks to see how much two histograms are overlapping, given their min
    #and max vals
    if (h1[0] <= h2[0]) and (h1[1] >= h2[1]):
        #complete overlap
        return 1.0
    elif h1[0] >= h2[0] and h1[1] <= h2[1]:
        #complete overlap the other way
        return 1.0
    elif h1[0] <= h2[0]:
        #return minimum percentage overlap
        olap = h1[1]-h2[0]
        polap1 = olap/(h1[1]-h1[0])
        polap2 = olap/(h2[1]-h2[0])
        return (min(polap1,polap2))
    else:
        olap = h2[1]-h1[0]
        polap1 = olap/(h1[1]-h1[0])
        polap2 = olap/(h2[1]-h2[0])
        return (min(polap1,polap2))

def pmHistos(hs,runs,ncut,xcut,nxcut):
    #now we've got those histos from sortedHistos and we want to check how our
    #overlap is
    #so what we do is run from L -> R and check overlap with R neighbor
    #if overlap < ncut (too small), we add a location midway between the
    #neighbors' highest points to newsampledists
    #else if the overlap > xcut (too big), we check if curr overlaps with R
    #neighbor's neighbor; if it does, we throw out R, but this is a more conservative cutoff than 
    #to add a neighbor (harder to throw away something already existing)
    print "ncut: {0}".format(ncut)
    print "xcut: {0}".format(xcut)
    print "nxcut: {0}".format(nxcut)
    newsampledists = []
    loseruns = []
    keepruns = list(runs)
    for i in range(np.shape(hs)[0]-1):
        olap = checkOverlap([hs[i,1],hs[i,2]],[hs[i+1,1],hs[i+1,2]])
        if olap < ncut:
            midpoint = (hs[i,0]+hs[i+1,0])/2.0
            newsampledists.append(midpoint)
        elif olap > xcut:
            if (i+2) < np.shape(hs)[0]:
                nolap = checkOverlap([hs[i,1],hs[i,2]],[hs[i+2,1],hs[i+2,2]])
                if nolap > nxcut:
                    loseruns.append(keepruns[i+1])
        else:
            continue

    for l in loseruns:
        keepruns.remove(l)
    return (newsampledists,keepruns,loseruns)
    
def findFrame(newd,dtabname,runs):
    #given a location that isn't well covered, find the frame closest to that
    #location from a summary distances type table which hasn't already been run
    d = open(dtabname,'r')
    dlines = d.readlines()
    d.close()
    dtable = np.zeros([len(dlines)-len(runs),2])
    i = 0
    for line in dlines:
        spline = line.split()
        r = int(spline[0])
        if r not in runs:
            dtable[i,0] = r
            dtable[i,1] = float(spline[1])
            i+=1
    dloc = np.argmin(abs(dtable[:,1]-newd))
    dframe = dtable[dloc,0]
    return int(dframe)
    
def newHistos(histos,keepruns,cinds,outname):
    #write out the new histograms for sanity check    
    o = open(outname,"w")
    #newhistos = np.zeros([np.shape(histos)[0],len(keepruns)])
    #k = 0
    for i in range(np.shape(histos)[0]):
        o.write('\t{0}'.format(histos[i,0]))
        for j in range(len(cinds)):
            if cinds[j] in keepruns:
                o.write('\t{0}'.format(histos[i,j+1]))
        o.write('\n')
    o.close()
    
def outputFiles(newdists,keepruns,loseruns,ndf,kf,lf,summaryd,framefolder,extendsteps,topolname,runprefix):
    #write out new tprfiles.dat, pullx-files.dat, pullf-files.dat
    #also write out file telling what new distances need to be calculated and what runs are being kept and lost
    ndfile = open(ndf,'w')
    for nd in newdists:
        ndfile.write('{0}\n'.format(nd))
    krfile = open(kf,'w')
    for kr in keepruns:
        krfile.write('{0}\n'.format(kr))
    
    lrfile = open(lf,'w')
    for lr in loseruns:
        lrfile.write('{0}\n'.format(lr))
    ndfile.close()
    krfile.close()
    lrfile.close()
    
    newruns = computeNewRuns(newdists,summaryd,keepruns,loseruns)
    tprf = open('tpr-files-cont.dat','w')
    pullxf = open('pullx-files-cont.dat','w')
    pullff = open('pullf-files-cont.dat','w')
    
    for k in keepruns:
        tprf.write('run{0}/umbrella{0}.tpr\n'.format(k))
        pullxf.write('run{0}/pullx-umbrella{0}.xvg\n'.format(k))
        pullff.write('run{0}/pullf-umbrella{0}.xvg\n'.format(k))
        
    for n in newruns:
        tprf.write('run{0}/umbrella{0}.tpr\n'.format(n))
        pullxf.write('run{0}/pullx-umbrella{0}.xvg\n'.format(n))
        pullff.write('run{0}/pullf-umbrella{0}.xvg\n'.format(n))
    
    tprf.close()
    pullxf.close()
    pullff.close()
    
    #finally set up sge cont files and sge full run files
    setupsge(keepruns,newruns,framefolder,extendsteps,topolname,runprefix)

def setupsge(keepruns,newruns,framefolder,extendsteps,topolname,runprefix):

    #for new runs, create initial frame file from frame-0 files by replacing all instances of
    # umbrella0,run0,conf0,npt0 with the 0 replaced by the number in question
    replaceterms = ['umbrella','run','conf','npt']
    frame0 = open(framefolder+'/frame-0_umbrella.sge','r')
    framedummy = frame0.readlines()
    frame0.close()
    for run in newruns:
        framename = '{0}/frame-{1}_umbrella.sge'.format(framefolder,run)
        framenew = open(framename,'w')
        for line in framedummy:
            for r in replaceterms:
                line = line.replace(r+'0',r+str(run))
            framenew.write(line)
        framenew.close()
    #next write cont files for all the keepruns and the newruns frames
    for run in newruns+keepruns:
        framename = '{0}/cont-frame-{1}_umbrella.sge'.format(framefolder,run)
        framecont = open(framename,'w')
        for line in framedummy:
            writeflag = 1
            if 'grompp -v -f npt_umbrella.mdp' in line or '-s npt0.tpr -o npt0.trr -c npt0.gro -g npt0.log -e npt0.edr -x npt0.xtc -cpo npt0.cpt' in line:
                writeflag = 0
            if writeflag:
                if 'grompp -v -f umbrella.mdp' in line:
                    line = '{2}tpbconv -s umbrella{0}.tpr -o umbrella{0}-cont.tpr -extend {1}\n'.format(run,extendsteps,runprefix)
                elif '-s umbrella0.tpr -o umbrella0.trr -c umbrella0.gro' in line:
                    line = line.replace('-cpo umbrella0.cpt',' ')
                    line = line.replace('-s umbrella0.tpr','-s umbrella0-cont.tpr')
                    line = line.strip() + ' -cpo umbrella0cont.cpt -cpi umbrella0.cpt -append\n'
                    for r in replaceterms:
                        line = line.replace(r+'0',r+str(run))
                else:
                    for r in replaceterms:
                        line = line.replace(r+'0',r+str(run))
                framecont.write(line)
        framecont.close()
    #finally, write a short bash script that will submit new files, continuations, and continuations of new files, with the last being held until the original new files are already run
    subsh = open(framefolder+'/submit.sh','w')
    subsh.write('#!/bin/bash\n')
    subsh.write('#This bash submission script was written out by histoCover.py.\n')
    subsh.write('#It should be set up to submit a mix of continuations and new runs for continued umbrella sampling runs.\n')
    for run in newruns:
        subsh.write('JID=$(qsub frame-{0}_umbrella.sge | cut -f3 -d" ")\n'.format(run))
        subsh.write('qsub -hold_jid $JID cont-frame-{0}_umbrella.sge\n'.format(run))
    for run in keepruns:
        subsh.write('qsub cont-frame-{0}_umbrella.sge\n'.format(run))
    subsh.close()
            
    
def computeNewRuns(newdists,summarydists,keepruns,loseruns):
    #takes a list of new places histograms need to be centered and the name of the summarydists file that lists what the distances for each of the runs is
    #returns the names of the runs that will place you closest to the desired distances that are not in the loseruns or keepruns list
    sd = open(summarydists,'r')
    lines = sd.readlines()
    sd.close()
    dmap = np.zeros([len(lines),2])
    newruns = []
    ind = 0
    for line in lines:
        spline = line.split()
        #only bother to check distances that aren't already kept or lost
        if not (keepruns.count(int(spline[0])) or loseruns.count(int(spline[0]))):
            dmap[ind,:] = map(float,spline)
            ind += 1
    dmap = dmap[0:ind,:]
    #now find the closest distance for each newdist
    print "Newly discovered distances with accompanying runs.  May be checked with {0}:".format(summarydists)
    for nd in newdists:
        
        run = int(dmap[abs(dmap[:,1]-nd).argmin(),0])
        print '{0} {1}\n'.format(nd,run)
        newruns.append(run)
    
    return newruns
    
def main():
    #take a set of command line arguments including input file names and write out list of new files
    print "Take care because this function is set up to use default values of all the possible parameters, but they may not match your usage needs.\n"
    parser = argparse.ArgumentParser()
    parser.add_argument("-f","--inputfile",help="name of input histogram xvg file")
    parser.add_argument("-o","--outputfile",help="name of output histogram xvg file")
    parser.add_argument("-n","--newdistsfile",help="file containing distances that are not covered")
    parser.add_argument("-k","--keeprunsfile",help="file containing the names of runs to be kept")
    parser.add_argument("-l","--loserunsfile",help="file containint the names of runs to be lost")
    parser.add_argument("-nc","--mincut",help="minimum cutoff")
    parser.add_argument("-xc","--maxcut",help="maximum cutoff")
    parser.add_argument("-nx","--nextcut",help="cutoff for next neighbor")
    parser.add_argument("-t","--tprlist",help="file containing list of tpr files")
    parser.add_argument("-s","--summarydists",help="file output in 4_pull by setupUmbrella.py")
    parser.add_argument("-ff","--framefolder",help="folder where .sge files are found")
    parser.add_argument("-ex","--extend",help="number of steps to extend continuing runs")
    parser.add_argument("-tn","--topolname",help="topology file name")
    parser.add_argument("-rp","--runprefix",help="prefix for grompp running")
    args = parser.parse_args()
    if args.inputfile:
        f = args.inputfile
    else:
        f = "histinit.xvg"
    if args.outputfile:
        o = args.outputfile
    else:
        o = "newhistinit.xvg"
    if args.newdistsfile:
        n = args.newdistsfile
    else:
        n = "newdists.dat"
    if args.keeprunsfile:
        kr = args.keeprunsfile
    else:
        kr = "keepruns.dat"
    if args.loserunsfile:
        lr = args.loserunsfile
    else:
        lr = "loseruns.dat"
    if args.mincut:
        nc = float(args.mincut)
    else:
        nc = 0.5
    if args.maxcut:
        xc = float(args.maxcut)
    else:
        xc = 0.75
    if args.nextcut:
        nx = float(args.nextcut)
    else:
        nx = 0.5
    if args.tprlist:
        t = args.tprlist
    else:
        t = "tpr-files.dat"
    if args.summarydists:
        summaryd = args.summarydists
    else:
        summaryd = "summary_distances.txt"
    if args.framefolder:
        ff = args.framefolder
    else:
        ff = '/home/rachael/cluster/home/mansbac2/coarsegraining/orient_restrain_mon_PMF/arom_line_restrain/end-end'
    if args.extend:
        ex = int(args.extend)
    else:
        ex = 5000
    if args.topolname:
        topolname = args.topolname
    else:
        topolname = "CG_DFAG"
    if args.runprefix:
        runprefix = args.runprefix
    else:
        runprefix = ""
    
    (histos,cinds) = parseHFile(f,t)
    (hs,runs) = sortedHistos(histos,cinds)
    (nd,k,l) = pmHistos(hs,runs,nc,xc,nx)
    print "\nKeeping {0} runs.".format(len(k))
    print "Losing {0} runs.".format(len(l))
    print "Adding {0} runs.\n".format(len(nd))
    newHistos(histos,k,cinds,o)
    outputFiles(nd,k,l,n,kr,lr,summaryd,ff,ex,topolname,runprefix)
    
if __name__ == "__main__":
    main()
    '''    
    mincutoff = 0.5
    maxcutoff = 0.75
    nxcutoff = 0.5
    (histos,cinds) = parseHFile("histinit.xvg","tpr-files.dat")
    (hs,runs) = sortedHistos(histos,cinds)
    (n,k,l) = pmHistos(hs,runs,mincutoff,maxcutoff,nxcutoff)
    newHistos(histos,k,cinds,'newhistinit.xvg')
     '''
