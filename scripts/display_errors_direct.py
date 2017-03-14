#!/usr/bin/python

# Code similar to display_results.py.
# Instead of just plotting all the points,
# assumes all the inputs have the same evaluation points,
# and displays errors of second, third, ... files compared
# to the first.

from numpy import *
import sys



# Read the reference data to plot against.
allRef = {};

fo = open(sys.argv[1],'r');
ref_data = fo.readlines();

allRef["maxM"] = int(ref_data[0]);
allRef["maxN"] = int(ref_data[1]);
allRef["flow"] = ref_data[2].split()[0];

nAratios = int(ref_data[3]);
allRef["skews"] = zeros( (nAratios,2) );


i=0;
for line in ref_data[4:]:
     nums = line.split();
     allRef["skews"][i,0] = float(nums[0]);
     allRef["skews"][i,1] = float(nums[1]);
     i += 1;
     
fo.close();

# end for

# Read the rest of the data.
allData = {};

for file in sys.argv[2:]:
     
     fo = open(file,'r');
     data = fo.readlines();
          

     allData[file] = {};
     allData[file]["maxM"] = int(data[0]);
     allData[file]["maxN"] = int(data[1]);
     allData[file]["flow"] = data[2].split()[0];
     nAratios = int(data[3]);
     allData[file]["skews"] = zeros( (nAratios,2) );

     i=0;
     
     for line in data[4:]:
          nums = line.split();
          allData[file]["skews"][i,0] = float(nums[0]);
          allData[file]["skews"][i,1] = float(nums[1]);
          i += 1;

     # end if
     fo.close();
     
# end for

# moo

for data in allData:

     if (size(allData[data]["skews"]) == 2):
          y = log10(abs(allRef["skews"][:,1] - allData[data]["skews"][:,1]));
          
          plot(allRef["skews"][:,0], y, label=file, marker='.', markersize=10);
     else:
          y = log10(abs(allRef["skews"][:,1] - allData[data]["skews"][:,1]));
          
          plot(allRef["skews"][:,0], y, label=file, marker='.', markersize=5);


grid()
#legend(loc=4)
show(block=False)
savefig('errors.png')
