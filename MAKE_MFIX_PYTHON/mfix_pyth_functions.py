# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 11:42:41 2013
make_mfix file related functions
@author: rgarg
"""
import sys 

import argparse
import re
import subprocess
import time 

def mfix_print_version(vers):
    print "=============================================================="
    print ""
    print "       MFIX: Multiphase Flow with Intephase eXchanges"
    print ""
    print "                      mfix.netl.doe.gov"
    print ""
    print "                       Version ", vers
    print ""
    print "=============================================================="
    return

def change_number(y):
    y = y*2
    return 
def set_machinefile(opsys):
    machfile="LINUX.F"
    if "OSF1" in opsys:
        machfile = "OSF1.f"
    elif "AIX" in opsys:
        machfile = "AIX.F"
    cmd = "chmod +w " +  "machine.f"
    #print cmd 
    subprocess.call(cmd, shell = True)
    cmd  = "/bin/cp -f " + machfile + " machine.f"
    subprocess.call(cmd , shell = True)
    return machfile
    
def read_compilation_flags(comp_options_file_loc):       
    """
    reads the old compilation file and returns a dictionary
    """
    old_options = {} # initialize it as dictionary
    with open(comp_options_file_loc, "r") as opt_file:
        #pass
        for line in opt_file :
            line = line.lstrip(" ") #remove any leading white space
            line= line.strip("\n")
            #print "line = ", line
            if line.startswith("#"):   
                pass # reserved for comments 
            else:           
                #print "line before = ", line
                line  = re.sub('\s+',' ',line)
                #the above will remove every multiple white spaces, tabs and newlines with 
                #one white space 
                #line = line.replace("\t", "  ") #no need for this 
                #print "line = ", line
                line_list = line.split(" ")
                line_list2 = []
                for items in line_list:
                    #re-make a list that does not contain any white space or hashes
                    if items in (""," ", "#"):
                        pass
                    else:
                        line_list2.append(items)  
                #print "line 2 = ", line_list2
                if len(line_list2) == 1 :
                    print "value for {0:s} not specified, so setting it to none.".format(line_list2[0])
                    print "Could lead to unforeseen errors. Recommended way is to specify something"
                    old_options[line_list2[0].lower()] = 'None'
                else:                        
                    old_options[line_list2[0].lower()] = line_list2[1]            
    return old_options


def print_user_options(user_options):
    print "\n******************************************************************\n"
    print "MFIX will be compiled for the following conditions"
    for opt, value in user_options.iteritems() :
        print "{0:<s} \t {1:>s}".format(opt.upper(), str(value))
    
    print "\n******************************************************************\n"
    return 

def get_or_set_value_from_dictionary(dictionary, dict_key, default_val):
    try:
        value    = dictionary[dict_key]
    except KeyError as e:
        print "*************************WARNING***************************************"
        print "Could not find a entry for {0:<s}".format(dict_key)
        print "Setting {0:<s} to default value of {1:>s}".format(dict_key, default_val)        
        print "*************************WARNING***************************************"
        value =  default_val
    return value 

def check_mode_of_compilation(serial_mode, smp_mode, dmp_mode, smp_dmp_mode):
    
    mode_count = 0 
        
    if serial_mode : mode_count = mode_count + 1
    if smp_mode :    mode_count = mode_count + 1
    if dmp_mode :    mode_count = mode_count + 1
    if smp_dmp_mode :    mode_count = mode_count + 1
    
    if (mode_count == 0):   
        print "Mode of compilation could not be determined unambiguously"
        print "serial_mode ({0:s}), smp_mode ({1:s}), dmp_mode ({2:s}), smp_dmp_mode ({3:s})".\
        format(str(serial_mode), str(smp_mode), str(dmp_mode), str(smp_dmp_mode))
        print "all specified as False", \
        "EXITTING"
        exit()
    elif (mode_count > 1):
        print "Mode of compilation could not be determined unambiguously"
        print "More then one out of the following:"
        print "serial_mode ({0:s}),smp_mode ({1:s}), dmp_mode ({2:s}), smp_dmp_mode ({3:s})".\
        format(str(serial_mode), str(smp_mode), str(dmp_mode), str(smp_dmp_mode))  
        print "specified as True. \nEXITTING"
        exit()
    
    return
        
def get_object_dirname(debug_mode, smp_mode, dmp_mode, smp_dmp_mode):
    DPO1 = "OBJECT_FILES"
    if debug_mode:
        DPO1 = DPO1 + '_DBG'
    else:
        DPO1 = DPO1 + '_OPT'
    
    if smp_mode:
        DPO1 = DPO1 + '_SMP'
    elif dmp_mode:
        DPO1 = DPO1 + '_DMP'
    elif smp_dmp_mode:
        DPO1 = DPO1 + '_HYB'
    
    return DPO1
    
    
def delete_last_compiled_objects(DPO):
    import os
    import os.path

    if os.path.exists(DPO):
        print "Directory {0:s} exists. Will be deleting the object and module files in this directory".format(DPO)
        
        for currdir, subdir, files in os.walk(DPO):
            for file in files :     
                #print file 
                if(os.path.splitext(file))[1].lower() in ('.o', '.mod', '.a'):
                    os.remove(os.path.join(currdir, file))
            
        print "Done!"
        print "Directory {0:s} is not deleted to prevent any inadvertent deletion of important files".format(DPO)
        
    else:
        print "Object Directory {0:s} from the last build does not exist.\nDoing nothing and exitting!".format(DPO)
        
    exit()
    return     
    
        
def replace_dmp_with_donothing(mfile_new_loc, mfile_final_loc):
    
    mfile_final_fobj = open(mfile_final_loc, "w")
    
    with open(mfile_new_loc, "r") as temp_make_file:        
        for oline in temp_make_file:
            line = oline.replace("dmp_modules", "dmp_modules/mpi_donothing")
            line2 = line.replace("des/desmpi_wrapper_mod", "des/mpi_donothing/desmpi_wrapper_mod")
            mfile_final_fobj.write(line2)       
            
    mfile_final_fobj.close()
    return 

def hold():
    junk = raw_input("Enter to continue ") 
    return
    
def backup_original_source(run_dir, mfix):
    
    force_recomp =  False
    check_for_changed_or_backup(run_dir, mfix, ".", force_recomp, 'model')
    check_for_changed_or_backup(run_dir, mfix, "cartesian_grid", force_recomp, 'cartesian_grid')
    check_for_changed_or_backup(run_dir, mfix, "chem", force_recomp, 'chem')
    check_for_changed_or_backup(run_dir, mfix, "cohesion", force_recomp, 'cohesion')
    check_for_changed_or_backup(run_dir, mfix, "des", force_recomp, 'des')
    check_for_changed_or_backup(run_dir, mfix, "des/mpi_donothing", force_recomp, "des/mpi_donothing")
    check_for_changed_or_backup(run_dir, mfix, "dmp_modules", force_recomp, "dmp_modules")
    check_for_changed_or_backup(run_dir, mfix, "dqmom", force_recomp, "dqmom")
    check_for_changed_or_backup(run_dir, mfix, "GhdTheory", force_recomp, "GhdTheory")
    check_for_changed_or_backup(run_dir, mfix, "qmomk", force_recomp, "qmomk")
    check_for_changed_or_backup(run_dir, mfix, "thermochemical", force_recomp, "thermochemical")

def check_for_changed_or_backup(run_dir, mfix, tdir, force_recomp, fold_name):
    import os 
    import os.path
    import filecmp
    import subprocess
    subdir = "/" + tdir 
    w_mfixdir = os.path.join(mfix, tdir)
    w_rundir  = os.path.join(run_dir, tdir)
        
    print "--------------------------------------------------------------------------------------"
    
    print "Backing up or restoring files in \"{0:<s}\" folder from the run directory".format(fold_name)
    
    if os.path.isdir(w_rundir):
        #print "Files in \"{0:<s}\" folder from the run directory used to build MFIX".format(fold_name)
        for file in os.listdir(w_rundir):
            filepath = os.path.join(w_rundir, file)
            if (os.path.isfile(filepath)) and  (os.path.splitext(file))[1] in ('.f', '.inc'):
                #print file 
                extension = os.path.splitext(file)[1]
                backup_file = os.path.splitext(file)[0] + ".0" + extension[1:]
        
                real_file_loc_in_src     = os.path.join(w_mfixdir, file)
                real_file_loc_in_run     = os.path.join(w_rundir , file)
                
                backup_file_loc_in_src = os.path.join(w_mfixdir, backup_file)
                #print "backup file = ", backup_file_loc_in_src, os.path.isfile(backup_file_loc_in_src), os.path.lexists(backup_file_loc_in_src)
                #print "source file = ", real_file_loc_in_src, os.path.isfile(real_file_loc_in_src), os.path.lexists(real_file_loc_in_src)
                #print "changed file = ", real_file_loc_in_run
                #os.isfile()
                if os.path.isfile(real_file_loc_in_src):
                    if not os.path.isfile(backup_file_loc_in_src) :
                        cmd = "cp " + real_file_loc_in_src + " " + backup_file_loc_in_src
                        #print "no backup yet"
                        #print "cmd ", cmd                 
                        subprocess.call(cmd, shell = True)
                    
    else:
        print "This directory \"{0:<s}\" does not exist in the run directory, so no new files will come from here".format(tdir)
        print "Will be reverting the earlier backed up files (if any)."
        
            #hold()
    for file in os.listdir(w_mfixdir):
        filepath = os.path.join(w_mfixdir, file)
        if (os.path.isfile(filepath)) and  (os.path.splitext(file))[1] in ('.0f', '.0inc'):
            #print file 
            extension = os.path.splitext(file)[1]
            real_file = os.path.splitext(file)[0] + "." + extension[2:]
            backedup_file_loc_in_src = os.path.join(w_mfixdir, file)
            real_file_loc_in_src     = os.path.join(w_mfixdir, real_file)
            real_file_loc_in_run     = os.path.join(w_rundir , real_file)
            #print "backed file = ", backedup_file_loc_in_src 
            #print "source file = ", real_file_loc_in_src 
        
            #print "changed file = ", real_file_loc_in_run
        
            #print real_file                 
            if os.path.isfile(real_file_loc_in_run) :
                
                
                if filecmp.cmp(real_file_loc_in_src, real_file_loc_in_run):
                    print "File in rundir, unchanged from last comp: \"{0:<s}\"".format(real_file)
                    if force_recomp:
                        cmd = "cp " + real_file_loc_in_run + " " + real_file_loc_in_src
                        subprocess.call(cmd, shell = True)                   
                        
                else:
                    print "File in rundir: changed from last comp  : \"{0:<s}\"".format(real_file)
                    cmd = "cp " + real_file_loc_in_run + " " + real_file_loc_in_src
                    subprocess.call(cmd, shell = True)
                    
                    #print filecmp.cmp(backedup_file_loc_in_src, real_file_loc_in_run), cmd
                    
                    #hold()  
                    
            else:
                if filecmp.cmp(real_file_loc_in_src, backedup_file_loc_in_src):
                    #if they are same, just remove the backed up file
                    
                    print "Restored file that didn't really change: \"{0:<s}\"".format(real_file)
                    os.remove(backedup_file_loc_in_src)
                else:
                    print "Restored file that did change: \"{0:<s}\"".format(real_file)
                    #Restore from the backed up file
                    cmd = "mv " + backedup_file_loc_in_src + " " + real_file_loc_in_src
                    subprocess.call(cmd, shell = True)
                    cmd = "touch " + real_file_loc_in_src
                    subprocess.call(cmd, shell = True)
                    
    #hold()
                    
    print "--------------------------------------------------------------------------------------"
    return
    
def write_directory_path(mfixpath):
    import os 
    import os.path
    import subprocess    
    filename = "mfix_directory_path.inc"
    file_loc = os.path.join(mfixpath, filename)
    cmd  = "tail -1 " + file_loc
    oline  = subprocess.check_output(cmd, shell = True)
    
    line_des = "     CHARACTER(len=132) :: MFIX_PATH = " + "\'" + mfixpath.rstrip('/') + "\'"
    
    please_write  = False
    if oline != line_des:
        print oline, line_des
        please_write = True
     
    if please_write :
        with open(file_loc, "w") as tfile:                
            tfile.write(line_des)
            
            print "--------------------------------------------------------------------------"
            print "filename", filename, "will be changed"
            print "Swapping {0:<s} \n with {1:>s}".format(oline, line_des)        
            print "This happens on first installation on a new machine or if you change the name of mfix root folder"
            print "--------------------------------------------------------------------------"
                
    return
            
def check_comp_flags(F):
    
    try:
        if F.libs        == ""   : raise_error("libs", F)
    except AttributeError as e:
        raise_flag_undefined_error("libs", F, e)
        
    try:
        if F.libs_d      == ""   : raise_error("libs_d", F)
    except AttributeError as e:
        raise_flag_undefined_error("libs_d", F, e)
    
    try:
        if F.link        == ""   : raise_error("link", F)
    except AttributeError as e:
        raise_flag_undefined_error("link", F, e)
    
    try:
        if F.link_d      == ""   : raise_error("link_d", F)        
    except AttributeError as e:
        raise_flag_undefined_error("link_d", F, e)
        
    try:
        if F.omp         == ""   : raise_error("omp", F)        
    except AttributeError as e:
        raise_flag_undefined_error("omp", F, e)
    
    try:
        if F.OBJ_EXT     == ""   : raise_error("OBJ_EXT", F)        
    except AttributeError as e:
        raise_flag_undefined_error("OBJ_EXT", F, e)
    
    try:
        if F.FORTRAN_EXT == ""   : raise_error("FORTRAN_EXT", F)
    except AttributeError as e:
        raise_flag_undefined_error("FORTRAN_EXT", F, e)
        
    try:
        if F.FORTRAN_CMD == ""   : raise_error("FORTRAN_CMD", F)
    except AttributeError as e:
        raise_flag_undefined_error("FORTRAN_CMD", F, e)
        
    try:
        if F.compile1    == ""   : raise_error("compile1", F)
    except AttributeError as e:
        raise_flag_undefined_error("compile1", F, e)
    
    try:
        if F.compile2    == ""   : raise_error("compile2", F)
    except AttributeError as e:
        raise_flag_undefined_error("compile2", F, e)
    
    try:
        if F.compile3    == ""   : raise_error("compile3", F)
    except AttributeError as e:
        raise_flag_undefined_error("compile3", F, e)
        
    try:
        if F.compile_d1  == ""   : raise_error("compile_d1", F)
    except AttributeError as e:
        raise_flag_undefined_error("compile_d1", F, e)
    
    try:
        if F.compile_d2  == ""   : raise_error("compile_d2", F)
    except AttributeError as e:
        raise_flag_undefined_error("compile_d2", F, e)
    
    try:        
        if F.compile_d3  == ""   : raise_error("compile_d3", F)
    except AttributeError as e:
        raise_flag_undefined_error("compile_d3", F, e)
    return 

def raise_flag_undefined_error(param_name, modname, error):
    print "************************ERROR BEGIN********************************************"
    print "The following error was raised for flag", param_name, "\n", error
    print "\nReferenced module for flags: \n", modname
    print "\nThis happens when a flag is not declared in the preamble of flags module"
    print "Edit the py file corresponding to the pyc file path to correct this error"
    print "EXITTING"
    
    print "************************ERROR END  ********************************************"

    exit()
def raise_error(param_name, modname):
    
    print "************************ERROR BEGIN********************************************"
    print "The following flag: ({0:<s}) was not initialized".format(param_name)
    print "This can happen either if the flag was not declared as global in the flags file"
    print "or if the user forgot to initialize it."    
    print "If you don't want it, initialize it as a whitespace"
    print "Referenced module for flags:", modname
    print "EXITTING"
    print "************************ERROR END  ********************************************"
    
    exit()
if __name__ == "__main__":
    import sys
    user_options = {}
    opsys = "Linux"
    machfile=set_machinefile(opsys)
    print "machine file overwritten with: ", machfile
    time.sleep(1)
    #fib(int(sys.argv[1]))
    #open('args.txt', 'w').write('-f\nbar')
    #parser = argparse.ArgumentParser(fromfile_prefix_chars='@')
    #parser.add_argument('-f')
    
    #print parser.parse_args(['@args.txt'])
    #parser.parse_args(['-f', 'foo', '@args.txt'])