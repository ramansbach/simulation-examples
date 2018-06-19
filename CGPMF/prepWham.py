#python script to prepare tprfiles and pullf and pullx file lists for use with g_wham
import glob
fnames = glob.glob("run*")
tprf = open("tpr-files.dat","w")
pullff = open("pullf-files.dat","w")
pullxf = open("pullx-files.dat","w")
for name in fnames:
	run = name.split("n")[1]
	tprf.write("run{0}/umbrella{0}.tpr\n".format(run))
	pullff.write("run{0}/pullf-umbrella{0}.xvg\n".format(run))
	pullxf.write("run{0}/pullx-umbrella{0}.xvg\n".format(run))

tprf.close()
pullff.close()
pullxf.close()
