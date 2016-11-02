#!/usr/bin/python

# 12/30/2015, Badri Adhikari
# Print R script to plot Chord diagram of a given contacts file
# Contact sources must be defined in an extra column in the RR file

import sys
import os
import math

if len(sys.argv) != 5:
	print(" Usage: " + sys.argv[0] + " contact-RR-file output-R-file output-R-png-file-name image-label")
	sys.exit()

contact  = sys.argv[1]
rfile    = sys.argv[2]
rimage   = sys.argv[3]
imagelbl = sys.argv[4]

if os.path.isfile(contact):
	print("Printing R chord diagram script..")
else:
	print("Contact file " + contact + " does not exist!")
	sys.exit()

aa2color = {
	"-" : "black",
	"A" : "cyan",
	"R" : "gold",
	"N" : "coral",
	"D" : "darkred",
	"B" : "dimgray",
	"C" : "dimgrey",
	"E" : "cornsilk",
	"Q" : "darkblue",
	"Z" : "darkcyan",
	"G" : "darkgray",
	"H" : "darkgrey",
	"I" : "deeppink",
	"L" : "chocolate",
	"K" : "darkgreen",
	"M" : "darkkhaki",
	"F" : "firebrick",
	"P" : "gainsboro",
	"S" : "goldenrod",
	"T" : "chartreuse",
	"W" : "darkorange",
	"Y" : "darkorchid",
	"V" : "darksalmon"
}

grpcolorlist = ["#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#A30059",
	"#7A4900", "#0000A6", "#4A3B53", "#D16100",
	"#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
	"#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
	"#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
	"#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",
	"#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
	"#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
	"#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
	"#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329",
	"#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C"
]

grpltylist = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24]

# compute the length of the sequence
L = 0
seq = ""
fCon = open(contact, "r")
for line in fCon:
	if line[0].isalpha() or line[0] == "-":
		L = len(line)
		seq = line
		break
	else:
		print("ERROR! First line must be a sequence!")
		sys.exit()
fCon.close()

# assign color to each residue number
rnum2color = {}
for i in range(1, L):
	rnum2color[i] = aa2color[seq[i-1]]

# get the list of all contact groups
groups = []
fCon = open(contact, "r")
for line in fCon:
	if line[0].isalpha() or line[0] == "-":
		continue
	cols = line.split()
	if len(cols) < 6:
		cols.append("all-contacts")
	if cols[5] in groups:
		continue;
	if not cols[5]:
			print("Column 5 must be defined for each row")
			sys.exit()
	groups.append(cols[5])
fCon.close()

# assign color to each group
group2color = {}
i = 0
for grp in groups:
	group2color[grp]= grpcolorlist[i]
	i = i + 1

# assign lty type to each group
group2lty = {}
i = 0
for grp in groups:
	group2lty[grp]= grpltylist[i]
	i = i + 1

# map residue number to angle1/angle2 and x/y coordinate in the plot
rnum2angle = {}
rnum2xy = {}
rnum2indexxy = {}
for i in range(1, L):
	a1 = round((i-1) * 2 * math.pi / L, 5)
	a2 = round(i * 2 * math.pi / L, 5)
	x  = round(0.5 + 0.5 * math.cos((a1+a2)/2), 5)
	y  = round(0.5 + 0.5 * math.sin((a1+a2)/2), 5)
	rnum2angle[str(i) + "a1"] = str(a1)
	rnum2angle[str(i) + "a2"] = str(a2)
	rnum2xy[str(i) + "x"] = x
	rnum2xy[str(i) + "y"] = y
	lx  = round(0.5 + 0.5 * math.cos((a1+a2)/2) + 0.02 * math.cos(i * 2 * math.pi / L), 5)
	ly  = round(0.5 + 0.5 * math.sin((a1+a2)/2) + 0.02 * math.sin(i * 2 * math.pi / L), 5)
	rnum2indexxy[str(i) + "x"] = lx
	rnum2indexxy[str(i) + "y"] = ly

# Initialize R script
fileR = open(rfile, 'w')
fileR.write("require(plotrix)\n")
fileR.write("png(file = \"" + rimage + "\", width=6.6, height=5, units=\"in\", res=500);\n")
fileR.write("par(mar=c(0, 0, 2, 10), xpd=TRUE);\n");
fileR.write("plot(0.5, 0.5, xlim = c(0, 1), ylim = c(0, 1), col = \"white\", main=\"" + imagelbl + "\", axes = FALSE);\n")
fileR.close()

# draw segments
fCon = open(contact, "r")
fileR = open(rfile, 'a')
for line in fCon:
	if line[0].isalpha() or line[0] == "-":
		continue
	cols = line.split()
	if len(cols) < 6:
		cols.append("all-contacts")
	i = cols[0]
	j = cols[1]
	fileR.write("segments(" + str(rnum2xy[str(i)+ "x"]) + ", " + str(rnum2xy[str(i)+ "y"]) + ", " + str(rnum2xy[str(j)+ "x"]) + ", " + str(rnum2xy[str(j)+ "y"]) + ", col = \"" + group2color[cols[5]] + "\", lty = " + str(group2lty[cols[5]]) + ", lwd = 1)\n")
fCon.close()
fileR.close()

# draw circle using arcs
fileR = open(rfile, 'a')
for i in range(1, L):
	# add the arc
	fileR.write("draw.arc(0.5, 0.5, radius = 0.5, angle1 = " + rnum2angle[str(i) + "a1"] + ", angle2 = " + rnum2angle[str(i) + "a2"] + ", col = \"" + rnum2color[i] + "\", lwd = 4);\n")
	# add residue number
	if L <= 100:
		fileR.write("text(" + str(rnum2indexxy[str(i)+ "x"]) + ", " + str(rnum2indexxy[str(i)+ "y"]) + ", cex = 0.5, " + str(i) + ");\n")
	elif L <= 200:
		if i % 2 == 1:
			fileR.write("text(" + str(rnum2indexxy[str(i)+ "x"]) + ", " + str(rnum2indexxy[str(i)+ "y"]) + ", cex = 0.5, " + str(i) + ");\n")
	else:
		if i % 3 == 1:
			fileR.write("text(" + str(rnum2indexxy[str(i)+ "x"]) + ", " + str(rnum2indexxy[str(i)+ "y"]) + ", cex = 0.5, " + str(i) + ");\n")

fileR.close()


# Add legend using the group color when there are at least 2 groups
if len(groups) > 1:
	fileR = open(rfile, 'a')
	fileR.write("legend(1.05, 1, c(")
	for i in range(0, len(groups)):
		if i == len(groups)-1:
			fileR.write("\"" + groups[i] + "\"")
		else:
			fileR.write("\"" + groups[i] + "\", ")
	fileR.write("), col = c(")
	for i in range(0, len(groups)):
		if i == len(groups)-1:
			fileR.write("\"" + group2color[groups[i]] + "\"")
		else:
			fileR.write("\"" + group2color[groups[i]] + "\", ")
	fileR.write("), lty = c(")
	for i in range(0, len(groups)):
		if i == len(groups)-1:
			fileR.write(str(group2lty[groups[i]]) + "")
		else:
			fileR.write(str(group2lty[groups[i]]) + ", ")
	fileR.write("))\n")
	fileR.close()
