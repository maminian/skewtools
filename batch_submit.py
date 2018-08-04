# Script to execute a single execute locally.
#

import subprocess,os,multiprocessing,random
import numpy as np

# -------------------------
# Parameters specifying the number of runs,
# and location and names for output files.
# -------------------

execute = "katrina"   # if "local" then run directly, sequentially.
                # "bsub" submits jobs to Kure/Killdevil.
                # "longleaf" submits jobs to Longleaf.
                # "katrina" behaves like "local" but uses the
                # multiprocessing package to submit nprocs
                # jobs at a time.

n = 96            # number of trials to run.
fname_prefix = "tr_"
parent = "./triangle-rev1/"
sim_folder = ""

exe_loc = "./triangle_mc"

# ------------------
# More bsub options

nprocs = 48 # number of processes; currently only used by "katrina" execute.
memoryreq = 16384	# in MB
queue = "week"

# ---------------------------

folder = parent+sim_folder
file_prefix = folder+fname_prefix
out_suffix = '.h5'

# ------------------------

def process_parameter_file(time_loc,i):
# Looks at the parameter file, checks if the line
# specifying the RNG seed is filled in with an
# integer. If it is, leave it. If not, replace the line
# with a seed based on the operating system RNG
# (eg, /dev/urandom.) Python handles this part automatically.

    tfile = open(time_loc,'r')
    lines = tfile.readlines()

    # If thing is not an integer, replace the line
    # with a random seed.
    sidx = 31 # Line which should hold the seed.
    seed = lines[sidx].split()
    try:
       int(seed)
#       int(seed[0])
    except:

       #   int1 = int(os.urandom(10).decode(encode('hex'),16)
       #   int2 = i*int(os.urandom(10).encode('hex'),16)
       int1 = random.SystemRandom().randint(0,10**16)
       int2 = random.SystemRandom().randint(0,10**16)

       # Bitwise XOR. Why? I don't know. Just because.
       int3 = int1^int2

       lines[sidx] = str(int3)[-8:-1]+' \n'
    # end try
    tfile.close()

    tfile = open(time_loc,'w')
    tfile.writelines(lines)
    tfile.close()

# end def

# Creating simulation initialization files.
for i in range(n):
    # Make a four-digit index (0000 through 9999)
    stridx = str(i).zfill(4)

    param_loc = folder+"parameters_mc_"+fname_prefix+stridx+".txt"
    out_loc = file_prefix+stridx

    # Copy the parameter files to the appropriate place, and
    # generate seeds if necessary.
    copy_command1 = ["cp","parameters_mc.txt",param_loc]

    subprocess.call(copy_command1)

    process_parameter_file(param_loc,i)
# end if

# Check if the appropriate out/ and err/ folders
# exist; if not, make them.
stdoutfolder = folder+"out/"
stderrfolder = folder+"err/"
if not os.path.exists(stdoutfolder):
    os.makedirs(stdoutfolder)
#
if not os.path.exists(stderrfolder):
    os.makedirs(stderrfolder)
#

# Running simulations after the files are created.
def execute_job(j):
    stridx = str(j).zfill(4)

    param_loc = folder+"parameters_mc_"+fname_prefix+stridx+".txt"
    out_loc = file_prefix+stridx+out_suffix

    if (execute=="bsub"):
        bsub_prefix = ["bsub","-n","1","-o",folder+"out/"+"out."+stridx,"-e",folder+"err/"+"err."+stridx]

        if (queue == "week"):
            more = ["-q",queue]
        else:
            more = ["-M",str(int(memoryreq)),"-q",queue]
        # end if

        exe_command = [exe_loc,param_loc,out_loc]

        command = bsub_prefix + more + exe_command
    # elif (execute=="local"):
    elif (execute in ["local", "katrina"]):
        exe_command = [exe_loc,param_loc,out_loc]

        command = exe_command
    elif (execute=="longleaf"):
        longleaf_prefix = ["sbatch","-n","1","-o",folder+"out/"+"out."+stridx,"-e",folder+"err/"+"err."+stridx]

        if (queue == "day"):
            tlimit = ["-t","%i:00"%(24,)]
        elif (queue == "week"):
            tlimit = ["-t","%i-%i"%(10,23)]
        else:
            tlimit = ["-t","%i:00"%(1,)]
        # end if
        mem = ["--mem","%i"%memoryreq]

#        exe_command = ["--wrap=\' %s"%exe_loc,"%s"%param_loc, "%s \'"%out_loc]
        exe_command = ["--wrap=%s %s %s"%(exe_loc,param_loc,out_loc)]
#        exe_command = ["--wrap","=",exe_loc,param_loc,out_loc]

        command = longleaf_prefix + tlimit + mem + exe_command
#        print command
    else:
        print("Unrecognized platform; \"execute\" must be one of \"local\", \"bsub\", or \"longleaf\".")
    # end if


    subprocess.call(command)
#

if execute!='katrina':
    for i in range(n): execute_job(i)
else:
    p = multiprocessing.Pool(nprocs)
    results = p.map(execute_job, np.arange(n, dtype=int))
#
