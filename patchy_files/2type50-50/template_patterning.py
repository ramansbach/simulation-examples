#python script that just simply updates a template file based on the given three inputs
import re,sys
A = sys.argv[1]
SCA = sys.argv[2]
SIGA = sys.argv[3]
SCB = sys.argv[4]
SIGB = sys.argv[5]
RUN = sys.argv[6]
templatefile = "production_template.py"
template = open(templatefile)
tlines = template.readlines()
template.close()
outfile = "production_template_temp.py"
outtemp=open(outfile,'w')
Afloat = float(A)/100.0
SCAfloat = float(SCA)
SCBfloat = float(SCB)
if SCA == '02':
    SCAfloat = 0.2
if SCB == '02':
    SCBfloat = 0.2
if SCA == '09':
    SCAfloat = 0.9
if SCB == '09':
    SCBfloat = 0.9
SIGAfloat = float(SIGA)/100.0
SIGBfloat = float(SIGB)/100.0
for line in tlines:
    line = line.replace('SSS','12345')
    line = line.replace('=AAA','='+str(Afloat))
    line = line.replace('[SCSCSCA','['+str(SCfloatA))
    line = line.replace('SCSCSCB]',str(SCfloatB)+']')
    line = line.replace('[ERADA','['+str(SIGfloatA))
    line = line.replace('ERADB]',str(SIGfloatB)+']')
    line = line.replace('_AAA','_'+A)
    line = line.replace('-SCSCSCA','-'+SCA)
    line = line.repalce('-SCSCSCB','-'+SCB)
    line = line.replace('-ERADA','-'+SIGA)
    line = line.replace('-ERADB','-'+SIGB)
    line = line.replace('RUN',RUN)
    outtemp.write(line)
outtemp.close()

