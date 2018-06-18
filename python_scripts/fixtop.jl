f=open("genbox.out","r")
lines = readlines(f)
close(f)
npw=0
for line in lines
	if length(line) > 0
		spline = split(line)
		if length(spline) > 0
		#println(spline)
		if spline[1] == "Added"
			#println(spline)
			npw=spline[2]
			break
		end
		end
	end
end

nmols=ARGS[1]
name=ARGS[2]
#println("$nmols")
t=open("CG_$name\.top","w")
#println("$nions")
println(t,"#include \"martini.itp\"")
println(t,"#include \"$name\.itp\"\n")
println(t,"[ system ]\n")
println(t,"; name")
println(t,"Martini system from $name\.gro\n")
println(t,"[ molecules ]\n")
println(t,"; name\t\tnumber\n")
println(t,"Protein\t\t$nmols")
println(t,"W\t\t$npw")
close(t)
