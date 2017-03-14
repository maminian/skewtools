import sys
import glob

import scripts.skewtools as skewtools

# Read the input argument. We're expecting something with a wildcard, 
# but a single file is okay too.
files = glob.glob(sys.argv[1])
files.sort()

# Gather the data from the files.
times = skewtools.gatherDataFromFiles(files,'Time')
skews = skewtools.gatherDataFromFiles(files,'Avgd_Skewness')
aratios = skewtools.gatherDataFromFiles(files,'aratio')
peclets = skewtools.gatherDataFromFiles(files,'Peclet')

#
# Process the data from the files and spit it out.
#
print ""
print "Peclet       Aspect ratio  Time of minimum"
print "---------    ------------  ---------------"

for i in range(len(files)):

     # For each simulation, 
     # search for the minimum skewness within a specified time range.
     tlow = times[i][1]
     thigh = 0.1
     
     skmin,imin = skewtools.constrainedMinimum(times[i],skews[i],tlow,thigh)

     tmin = times[i][imin[0]]
     aratio = aratios[i][0]
     peclet = peclets[i][0]
     
     print "%.3e    %.3e     %.3e" % (peclet, aratio, tmin)
# end for

print ""

