#python script that just simply updates a template file based on the given three inputs
import re,sys
from numpy.random import permutation
A = sys.argv[1]
SCA = sys.argv[2]
SIGA = sys.argv[3]
SCB = sys.argv[4]
SIGB = sys.argv[5]
RUN = sys.argv[6]
bfrac = float(sys.argv[7])

def gen_array(bfrac):
    """
    Creates a string of 'EA' and 'EB' in
    a specified fraction but a random
    order

    ----------
    Parameters
    ----------
    bfrac: float
        fraction of strings that should
        be 'EB'

    -------
    Returns
    -------
    mxarray: string
        string with randomly permuted
        'EA' and 'EB
    """
    tlen = 20
    nb = int(bfrac*tlen)
    na = tlen - nb
    mxlist = []
    for i in range(nb):
        mxlist.append('EB')
    for i in range(na):
        mxlist.append('EA')
    mxlist = permutation(mxlist)
    mxarray = '["'+mxlist[0]
    for i in range(1,len(mxlist)):
        mxarray += '","'+mxlist[i]
    mxarray += '"]'
    return mxarray


templatefile = "production_template.py"
template = open(templatefile)
tlines = template.readlines()
template.close()
outfile = "production_template_temp.py"
outtemp=open(outfile,'w')
Afloat = float(A)/100.0
SCfloatA = float(SCA)
SCfloatB = float(SCB)
if SCA == '02':
    SCfloatA = 0.2
if SCB == '02':
    SCfloatB = 0.2
if SCA == '09':
    SCfloatA = 0.9
if SCB == '09':
    SCfloatB = 0.9
if SCB == '001':
    SCfloatB = 0.01
SIGfloatA = float(SIGA)/100.0
SIGfloatB = float(SIGB)/100.0

mxarray = gen_array(bfrac)

for line in tlines:
    line = line.replace('SSS','12345')
    line = line.replace('=AAA','='+str(Afloat))
    line = line.replace('[SCSCSCA','['+str(SCfloatA))
    line = line.replace('SCSCSCB]',str(SCfloatB)+']')
    line = line.replace('[ERADA','['+str(SIGfloatA))
    line = line.replace(',ERADB',','+str(SIGfloatB))
    line = line.replace('_AAA','_'+A)
    line = line.replace('-SCSCSCA','-'+SCA)
    line = line.replace('-SCSCSCB','-'+SCB)
    line = line.replace('-ERADA','-'+SIGA)
    line = line.replace('-ERADB','-'+SIGB)
    line = line.replace('RUN',RUN)
    line = line.replace('MIXARRAY',mxarray)
    outtemp.write(line)
outtemp.close()

