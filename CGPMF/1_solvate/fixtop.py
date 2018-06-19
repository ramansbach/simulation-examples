"""
A simple script that writes out a CG topology file containing solvated water
"""
import argparse
def main():
    f=open("genbox.out","r")
    lines = f.readlines()
    f.close()
    npw=0
    for line in lines:
        if len(line) > 0:
            spline = line.split()
            if len(spline) > 0:
                if spline[0] == "Added":
                    npw=spline[1]
                    break
                
    parser = argparse.ArgumentParser(description='fix topology')
    parser.add_argument('nmols',metavar='M')
    parser.add_argument('name',metavar='N')
    args = parser.parse_args()

    nmols=args.nmols
    name=args.name
    t=open("CG_{0}.top".format(name),"w")
    t.write("#include \"martini.itp\"\n")
    t.write("#include \"{0}.itp\"\n\n".format(name))
    t.write("[ system ]\n\n")
    t.write("; name\n")
    t.write("Martini system from {0}.gro\n\n".format(name))
    t.write("[ molecules ]\n\n")
    t.write("; name\t\tnumber\n\n")
    t.write("Protein\t\t{0}\n".format(nmols))
    t.write("W\t\t{0}\n".format(npw))
    t.close()

if __name__ == "__main__":
    main()