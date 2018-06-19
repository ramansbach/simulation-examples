import argparse
parser = argparse.ArgumentParser(description='fix dimer gro file')
parser.add_argument('fname',metavar='F')
args = parser.parse_args()
fname = args.fname
ofname = 'dimer.gro'
f = open(fname,'r')
of = open(ofname,'w')
lines = f.readlines()
monsize = int(lines[1].split()[0])
of.write(lines[0])
of.write('\t{0}\n'.format(2*monsize))
for l in range(2,len(lines)):
    if (l!=monsize+2) and (l!=monsize+3) and (l!=monsize+4):
        of.write(lines[l])
of.close()
