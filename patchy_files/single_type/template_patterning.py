#python script that just simply updates a template file based on the given three inputs
import re,sys
A = sys.argv[1]
SC = sys.argv[2]
SIG = sys.argv[3]
RUN = sys.argv[4]
templatefile = "production_template.py"
template = open(templatefile)
tlines = template.readlines()
template.close()
outfile = "production_template_temp.py"
outtemp=open(outfile,'w')
Afloat = float(A)/100.0
SCfloat = float(SC)
if SC == '02':
    SCfloat = 0.2
if SC == '09':
    SCfloat = 0.9
SIGfloat = float(SIG)/100.0
for line in tlines:
    line = line.replace('SSS','12345')
    line = line.replace('=AAA','='+str(Afloat))
    line = line.replace('=SCSCSC','='+str(SCfloat))
    line = line.replace('=ERAD','='+str(SIGfloat))
    line = line.replace('_AAA','_'+A)
    line = line.replace('-SCSCSC','-'+SC)
    line = line.replace('-ERAD','-'+SIG)
    line = line.replace('RUN',RUN)
    outtemp.write(line)
outtemp.close()

