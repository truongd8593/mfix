# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 11:08:05 2013

@author: rgarg
"""
#!pythonl
import sys 

import sys 
import inspect, os 
import os.path
import mfix_pyth_functions

import argparse
import subprocess
import shutil 
import comp_flags 
import time 

def hold():
    junk = raw_input("Enter to continue ") 
    return

vers="2013-1"
mfix_pyth_functions.mfix_print_version(vers)
    
run_dir = os.getcwd()

thisfile_fullpath = inspect.getfile(inspect.currentframe())

thisfilename = os.path.split(thisfile_fullpath)[1]
python_utils_base = os.path.split(thisfile_fullpath)[0]

if python_utils_base in "":
#if it is a null string, then that implies this script was called from 
#Python utility folder itself. So get the full path in this case. 
    python_utils_base = os.getcwd()

#do the following to get the python utils base directory name in 
#in full format. This is to prevent errors arising out of user 
#specifying relative paths 
    
os.chdir(python_utils_base)
python_utils_base = os.getcwd()
os.chdir(run_dir)

comp_options_file_loc = os.path.join(python_utils_base, "LAST_COMP_OPTIONS")

mfix = python_utils_base.replace("/MAKE_MFIX_PYTHON", "/model/")
mfix_base = python_utils_base.replace("/MAKE_MFIX_PYTHON", "/")

tools_dir = python_utils_base.replace("/MAKE_MFIX_PYTHON", "/tools/")

print " run_dir:{0:<s}\n python util directory:{1:<s} \n mfix directory:{2:<s}".format(run_dir, python_utils_base, mfix)
print " tools dir :{0:<s}".format(tools_dir)
if os.path.samefile(run_dir, mfix ):
    print "*** Execute this command from any directory other than the current directory!"
    print "*** It is usually executed from a run directory containing user defined files." 
    exit()


opsys    = subprocess.check_output(['uname', '-s']).rstrip()
proctype = subprocess.check_output(['uname', '-p']).rstrip()

print "Operating system detected as: ", opsys
print "Processor type   detected as: ", proctype 
print "Continuing with MFIX build ......."
time.sleep(2)# To allow users the chance to read this everything so far

optim_level_choices_tup = ('O1', 'O2', 'O3')
compiler_choices_tup = ('gfortran', 'ifort', 'pgi', 'pathscale', 'user1', 'user2')

parser_parent = argparse.ArgumentParser(add_help=True)

init_group = parser_parent.add_mutually_exclusive_group()

init_group.add_argument('-r', '--repeat', action="store_true", dest = 'comp_rep',
                    default = False, help = "Flag to read the compilation \
                    from last known flags")
init_group.add_argument('-c', '--clean', action="store_true", dest = 'comp_del', default = False, 
                    help =   "Flag to delete files created from last compilation")

parser_parent.add_argument('--comp', '--compiler', action = "store", dest = "compiler",
                           help = "Select the compiler out of the following options. There should be an appropriate \
                           flags file or the script will complain.", 
                           default = 'gfortran')
                           #choices = compiler_choices_tup)
#
#

optim_debug_group = parser_parent.add_mutually_exclusive_group()

optim_debug_group.add_argument('--optim', '--opt', action="store_true", dest = 'optim_mode',
                    default = False, help = "Flag for compiling in optim mode, default choice!")
                    
#                    choices=('O1', 'O2', 'O3'))

optim_debug_group.add_argument('--debug', '--dbg', action="store_true", dest = 'debug_mode',
                    default = False, help = "Flag for compiling in debug mode")

serial_smp_dmp_group = parser_parent.add_mutually_exclusive_group()

serial_smp_dmp_group.add_argument('--serial', action="store_true", dest = 'serial_mode',
                    default = False, help = "Flag for compiling in serial mode. \n \
                    Specify one between -serial, -smp, -dmp and -smp_dmp.")

serial_smp_dmp_group.add_argument('--smp', action="store_true", dest = 'smp_mode',
                    default = False, help = "Flag for compiling in shared memory acceleration mode. \
                    Specify one between -serial, -smp, -dmp and -smp_dmp.")
                    
serial_smp_dmp_group.add_argument('--dmp', action="store_true", dest = 'dmp_mode',
                    default = False, help = "Flag for compiling in distributed memory acceleration mode. \
                    Specify one between -serial, -smp, -dmp and -smp_dmp.")
                    
serial_smp_dmp_group.add_argument('--smp_dmp', action="store_true", dest = 'smp_dmp_mode',
                    default = False, help = "Flag for compiling in hybrid shared and distributed memory \
                    acceleration mode. Specify one between -serial, -smp, -dmp and -smp_dmp.")
                    
args = parser_parent.parse_args()

read_options = args.comp_rep
clean_mfix   = args.comp_del
compiler     = args.compiler
optim_mode  = args.optim_mode
debug_mode   = args.debug_mode
serial_mode  = args.serial_mode
smp_mode     = args.smp_mode
dmp_mode     = args.dmp_mode
smp_dmp_mode     = args.smp_dmp_mode 

if read_options | clean_mfix :
    if os.path.isfile(comp_options_file_loc):
        user_options_in = {}
        print "Reading the options from the user specified or last known flags"        
        user_options_in = mfix_pyth_functions.read_compilation_flags(comp_options_file_loc)
    else:
        print "The following file does not exist", comp_options_file_loc, "\n \t EXITTIG"   
        exit()
        
    compiler = mfix_pyth_functions.get_or_set_value_from_dictionary(user_options_in, 'compiler', 'gfortran')
    
    optim_mode0 = mfix_pyth_functions.get_or_set_value_from_dictionary(user_options_in, 'optim_mode', 'False')    
    debug_mode0 = mfix_pyth_functions.get_or_set_value_from_dictionary(user_options_in, 'debug_mode', 'False')
    
    serial_mode0 = mfix_pyth_functions.get_or_set_value_from_dictionary(user_options_in, 'serial_mode', 'True')
    smp_mode0 = mfix_pyth_functions.get_or_set_value_from_dictionary(user_options_in, 'smp_mode', 'False')    
    dmp_mode0 = mfix_pyth_functions.get_or_set_value_from_dictionary(user_options_in, 'dmp_mode', 'False')    
    smp_dmp_mode0 = mfix_pyth_functions.get_or_set_value_from_dictionary(user_options_in, 'smp_dmp_mode', 'False')
       
       
    optim_mode  = (optim_mode0.lower() in ("true", "t", "1", "yes", "y"))
    debug_mode  = (debug_mode0.lower() in ("true", "t", "1", "yes", "y"))
    serial_mode = (serial_mode0.lower() in ("true", "t", "1", "yes", "y"))
    smp_mode    = (smp_mode0.lower() in ("true", "t", "1", "yes", "y"))
    dmp_mode    = (dmp_mode0.lower() in ("true", "t", "1", "yes", "y"))
    smp_dmp_mode    = (smp_dmp_mode0.lower() in ("true", "t", "1", "yes", "y"))
    
    del(optim_mode0); del(debug_mode0); del(serial_mode0); del(smp_mode0); del(dmp_mode0); del(smp_dmp_mode0)
   
    
if (optim_mode & debug_mode):
    print " Both optimization level ({0:s}) and debug mode ({1:s}) set to True".format(str(optim_mode), str(debug_mode))
    print " Choose between optim and debug mode.\n Call script with --help to list options EXITTING!"
    exit()

if (not optim_mode) & (not debug_mode):
    print " Both optimization level ({0:s}) and debug mode ({1:s}) set to False".format(str(optim_mode), str(debug_mode))
    print " Choose between optim and debug mode.\n Call script with --help to list options EXITTING!"
    exit()



#since the parser does not formce presence of at least one mode of compilation option, 
#check for the same for both fresh and repeat flag cases. 
mfix_pyth_functions.check_mode_of_compilation(serial_mode, smp_mode, dmp_mode, smp_dmp_mode)

os.chdir(mfix)
#print "now in mfix directory ", os.getcwd()

DPO1 = mfix_pyth_functions.get_object_dirname(debug_mode, smp_mode, dmp_mode, smp_dmp_mode)
DPO1 = DPO1 + "_" + compiler.upper() + "/"
DPO_ABS  = os.path.join(mfix_base, DPO1)
DPO = "../" + DPO1

#print "DPO_ABS", DPO_ABS, DPO1
if os.path.isdir(DPO_ABS):
    pass
else:
    os.mkdir(DPO_ABS)
#print "DPO and DPOABS  = ", DPO, DPO_ABS

if clean_mfix: mfix_pyth_functions.delete_last_compiled_objects(DPO)

flags_file =  compiler + "_flags"
flags_file1 = "comp_flags/" + flags_file + ".py"
flags_file_loc = os.path.join(python_utils_base, flags_file1)

if os.path.isfile(flags_file_loc):
    print "Success!\nFound a flags file \n\"{0:s})\"\ncorresponding to the specified compiler ({1:s})".format(flags_file_loc, compiler)
    print "Moving on with setting up flags for MFIX build"
else:
    print "Error!\nCould not find a flags file \n\"{0:s})\"\ncorresponding to the specified compiler ({1:s})".format(flags_file_loc, compiler)
    print "See the documentation on adding a new compiler option. \nEXITTING"
    exit()
    

user_options = {}
user_options['compiler'] = compiler
user_options['optim_mode'] = optim_mode
user_options['debug_mode'] = debug_mode
user_options['serial_mode'] = serial_mode
user_options['smp_mode'] = smp_mode
user_options['dmp_mode'] = dmp_mode
user_options['smp_dmp_mode'] = smp_dmp_mode

#print on screen the user options 
mfix_pyth_functions.print_user_options(user_options)

with open(comp_options_file_loc, "w") as opt_file:
    #write the flags to the file for future builds
    for optname, value in user_options.iteritems() :
        opt_file.write(optname + "   " + str(value) + " " +  "\n")
    
machfile=mfix_pyth_functions.set_machinefile(opsys)
print "machine file overwritten with: ", machfile
time.sleep(1)

#print "mpi_include = ", mpi_include

name = "comp_flags." + flags_file
F = __import__(name, fromlist=["*"])
F.get_comp_flags()
mfix_pyth_functions.check_comp_flags(F)

if smp_mode | smp_dmp_mode :
    mp_flag = F.omp + " " 
else:
    mp_flag = " "


FORT_COMMAND = F.FORTRAN_CMD
mpi_include = " "

if (dmp_mode | smp_dmp_mode):
    FORT_COMMAND = "mpif90"
    mpi_include = subprocess.check_output(['mpif90', '-showme:incdirs'])     
    mpi_include = mpi_include.rstrip('\n')
else:
    mpi_include=mfix + "/dmp_modules/mpi_donothing"

cmd = "ln -sf " + mpi_include + "/mpif.h  ." 
#print "symbolic linking command = ", cmd     
#hold()  
subprocess.call(cmd, shell=True)
#Since the compiler looks for mpif.h in the model directory, 
#it is necessary to define a symbolic here. Another option could have 
#been to replace mpif.h by mpi_include/mpif.h but that would require 
#changing more than one source file. So going with this approach of 
#creating a symbolic link 

if debug_mode :
    pass    
    FORT_FLAGS  = mp_flag + F.compile_d1
    FORT_FLAGS2 = mp_flag + F.compile_d2
    FORT_FLAGS3 = mp_flag + F.compile_d3
    LINK_FLAGS  = mp_flag + F.link_d
    LIB_FLAGS   = F.libs_d
else:
    pass
    FORT_FLAGS  = mp_flag +  F.compile1
    FORT_FLAGS2 = mp_flag +  F.compile2
    FORT_FLAGS3 = mp_flag +  F.compile3
    LINK_FLAGS  = mp_flag +  F.link
    LIB_FLAGS   = F.libs

EXEC_FILE = "mfix.exe"

#Check if mfix_directory_path.inc has the correct path or not
mfix_pyth_functions.write_directory_path(mfix)

mfix_pyth_functions.backup_original_source(run_dir, mfix)
#the above will copy any new files to mfix directory and also 
#revert back any backups from previous build

#update the base mfix files, so that any dependencies are updated 
#since this is done everytime now, there is no need to force recompilation of run directory routines
#any new dependencies will be automatically captured
basemake_build_file_loc = os.path.join(tools_dir, "make_make_byPy.py")

if not os.path.isfile(basemake_build_file_loc):
    print "basemake build utility missing at: \"{0:<s}\" ".format(basemake_build_file_loc)
    print "Exitting"
    exit()
cmd ="cd " + tools_dir + " ; " + "python "+ basemake_build_file_loc + ";" + "cd " + os.getcwd()
print "renewing base mfix files by calling the following python utility: \"{0:<s}\"".format(basemake_build_file_loc)
subprocess.call(cmd, shell=True)


if F.makefile_upper_case:
    mfile = "mfix_u.make"
else:
    mfile = "mfix_l.make"

mfile =  mfile + F.makefile_append
mfile_loc = os.path.join(mfix, mfile)

print "Base mfix makefile used:", mfile

mfile_new_name = "mfix_temp_makefile"
mfile_new_loc = os.path.join(mfix, mfile_new_name)
mfile_new_fobj = open(mfile_new_loc, "w")

print "mfile_new_loc = ", mfile_new_loc

mfile_new_fobj.write("DPO=" +    DPO + "\n")
mfile_new_fobj.write("mpi_include=" + mpi_include.lstrip(" ").rstrip("\n") + "\n")
mfile_new_fobj.write("FORT_FLAGS=" + FORT_FLAGS.lstrip(" ") + "\n")
mfile_new_fobj.write("FORT_FLAGS2=" + FORT_FLAGS2.lstrip(" ")+ "\n")
mfile_new_fobj.write("FORT_FLAGS3=" + FORT_FLAGS3.lstrip(" ")+ "\n")
mfile_new_fobj.write("LIB_FLAGS=" + LIB_FLAGS.lstrip(" ")+ "\n")
mfile_new_fobj.write("OBJ_EXT=" + F.OBJ_EXT.lstrip(" ")+ "\n")
mfile_new_fobj.write("FORTRAN_CMD=" + FORT_COMMAND.lstrip(" ")+ "\n")
mfile_new_fobj.write("LINK_CMD=" + FORT_COMMAND.lstrip(" ")+ "\n")
mfile_new_fobj.write("EXEC_FILE=" + EXEC_FILE.lstrip(" ")+ "\n")

mfile_new_fobj.close()

cmd = "cat " + mfile_loc + " >> " + mfile_new_loc 
subprocess.call(cmd, shell=True)

mfile_final_name = "mfix_final_makefile"
mfile_final_loc = os.path.join(mfix, mfile_final_name)
mfile_final_fobj = open(mfile_final_loc, "w")

if serial_mode | smp_mode :
    mfix_pyth_functions.replace_dmp_with_donothing(mfile_new_loc, mfile_final_loc)
else:
    cmd  = "cp " + mfile_new_loc  + " " + mfile_final_loc
    subprocess.call(cmd, shell=True)




cmd = "make -f " + mfile_final_loc + "  "   + EXEC_FILE
print 'Compiling mfix. Final temporary make file is', mfile_final_loc
subprocess.call(cmd, shell=True)

cmd  = "cp " + EXEC_FILE + " "  + run_dir + "/"
subprocess.call(cmd, shell=True)


cmd  = "cp " + EXEC_FILE + " "  + DPO_ABS + "/"
subprocess.call(cmd, shell=True)
