#!/usr/bin/python

from pylab import *

# Read the input arguments as a list of files to pull data from.
# Then pull the information from each of those files.

allData = {};

for file in sys.argv[1:]:
     
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
     
     
     
     if (size(allData[file]["skews"]) == 2):
          plot(allData[file]["skews"][:,0],allData[file]["skews"][:,1], label=file, marker='.', markersize=10);
     else:
          allData[file]["plot"] = plot(allData[file]["skews"][:,0],allData[file]["skews"][:,1], label=file, marker='.', markersize=5);
          
     # end if
     fo.close();
     
# end for

# For now, just produce the figure that will come as a result.
grid()
legend()
ax = gca()
ax.set_xscale('log')
show(block=False)
savefig('plots.png')
