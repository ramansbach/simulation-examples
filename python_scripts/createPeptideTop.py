# -*- coding: utf-8 -*-
"""
Created on Fri Apr 27 08:18:14 2018

@author: rachael

Creates a Martini model of the type DXXX-OPV3-XXXD when run

Currently hard-coded to use DFAG-OPV3-GAFD as the base
"""
import createMartiniModel as cmm
import argparse
from martini22_ff import martini22

letterDict = dict()

letterDict['d'] = 'PAS';
letterDict['g'] = 'GLY';
letterDict['a'] = 'ALA';
letterDict['v'] = 'VAL';
letterDict['i'] = 'ILE';
letterDict['l'] = 'LEU';
letterDict['m'] = 'MET';
letterDict['f'] = 'PHE';
letterDict['w'] = 'TRP';
letterDict['y'] = 'TYR';
letterDict['e'] = 'GLU';

def main():
    parser = argparse.ArgumentParser(description='get side chain')
    parser.add_argument('sidechain',metavar='S')
    args = parser.parse_args()
    sidechain = args.sidechain
     
    chemistry = []
     
    for letter in sidechain:
        
        chemistry.append(letterDict[letter])
     
    Top1 = cmm.DXXXTopology('DFAG.itp','DFAG.gro')
    resIDs = [(1,15),(2,14),(3,13),(4,12)]
    for i in range(len(chemistry)-1,-1,-1):
        Top1.resSwap(chemistry[i],'C',resIDs[i][0])
        Top1.resSwap(chemistry[i],'C',resIDs[i][1])
    Top1.write(sidechain)

if __name__ == '__main__':
    main()


