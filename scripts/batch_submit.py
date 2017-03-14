# Script to submit batch jobs to Kure/Killdevil.
#
# Format as of 2015/01/09:
#
# bsub -n 1 -o outfile -e errfile ./duct_mc space_mc.txt time_mc.txt fname
#
# Where fname is whatever you want the output file name to be (what will be
# iterated over).

import subprocess

n = 200
fname_prefix = "dt"
folder = "dtest/"

exe_loc = "./duct_mc"
command_prefix = "bsub -n 1 -o out."

space_loc = "space_mc.txt"
time_loc = "time_mc.txt"
file_prefix = folder+"dt"


for i in range(1,n+1):
     stridx = str(i).zfill(4)
     
     exe_command = [exe_loc,space_loc,time_loc,file_prefix+stridx]
     bsub_prefix = ["bsub","-n","1","-o",folder+"out/"+"out."+stridx,"-e",folder+"err/"+"err."+stridx]
     
     bsub_command = bsub_prefix + exe_command
 
     subprocess.call(bsub_command)
# end for
