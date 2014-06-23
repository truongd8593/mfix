# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 15:43:17 2013

@author: rgarg
Developed with v2.7
"""
import sys 
import os 
import os.path
import shutil

def call_get_use_inc_list(filename):

    """function that reads one line at a time and looks for use modules and 
    include files. Returns a dictionary with key values of uselist and inclist."""
    use_list = []
    inc_list = []
    with open(filename, "r") as modefile:
        for oline in modefile :
            line = oline.lstrip() #strips all the beginning white spaces
            if line.lower().startswith("use "):
                usename =  line.split()[1].strip(",").lower() #since fortran is case
                #insensitive. and filenames are conventionally named lower case
                try:
                    use_list.index(usename)
                    continue 
                except ValueError as e:
                    use_list.append(usename)
            
            if line.lower().startswith("include "):
                inc_name =  line.split()[1].strip("'").strip("\"")
                #to remove quotes around include statment 
                try:
                    inc_list.index(inc_name)
                    continue 
                except ValueError as e:
                    inc_list.append(inc_name)
    #return the lists as a collection in dict object
    return {'uselist':use_list, 'inclist':inc_list}
   
   
# Function definition is here
def check_file_repeat( listoftuples, filename ):
    "This checks if this file is already in the list or not"
    count = 0
    file_repeat = False
    for filefullpath, file in listoftuples:
        count  = count +  1 
        #print "fll..= ", filefullpath, file, filename, count

        if file == filename:
            file_repeat = True
            break 
    if file_repeat:
        return True
    else:
        return False
            
def readchar():
    junk = raw_input("Enter to continue ") 
    return

curr_dir = os.getcwd()
tools_dir = curr_dir
os.chdir("../model")
MODEL_PATH = os.getcwd()

files_list_loc = os.path.join(tools_dir, "mfix_files.list")
mod_files_list_loc = os.path.join(tools_dir, "mfix_mod_files.list")

files_list = open(files_list_loc, "w")
mod_files_list = open(mod_files_list_loc, "w")

mfix_l_make_loc = os.path.join(tools_dir, "mfix_l.make")
mfix_u_make_loc = os.path.join(tools_dir, "mfix_u.make")

mfix_l = open(mfix_l_make_loc, "w")
mfix_u = open(mfix_u_make_loc, "w")

preamble = ".$(FORTRAN_EXT).$(OBJ_EXT):\n\t$(FORTRAN_CMD) $(FORT_FLAGS) $<"
preamble = preamble + "\n\n" + "$(EXEC_FILE) : \\" + "\n"
mfix_l.write(preamble)
mfix_u.write(preamble)

all_mod_files = []
non_mod_files = []
all_fort_files = []

for currdir, subdir, files in os.walk(os.getcwd()):
    for file in files : 
        if(os.path.splitext(file))[1] == ".f":
            fullpath = os.path.join(currdir, file)
            fullpath = fullpath.replace(MODEL_PATH , ".")
            filetuple = (fullpath, file)
            
            if (check_file_repeat( all_fort_files, file)):
                #file already exists in the list. so move on to the next file
                continue
            else:
                pass
            all_fort_files.append(filetuple)        
            files_list.write( str(all_fort_files[-1]) + "\n")
            onlyfilename = os.path.splitext(file)[0]
            
            if onlyfilename.lower().endswith("_mod") :
                modname = onlyfilename.replace("_mod", ".mod")
                modtuple = (fullpath, file, modname)
                all_mod_files.append(modtuple)
                
                mod_files_list.write(str(all_mod_files[-1]) + "\n")
                new_lline = "    $(DPO)" + modname + " \\" + "\n"
                
                new_uline = (os.path.splitext(modname)[0]).upper() + (os.path.splitext(modname)[1])               
                new_uline = "    $(DPO)" + new_uline + " \\" + "\n" 
                
                mfix_l.write(new_lline)
                mfix_u.write(new_uline)
            else:
                non_mod_files.append(filetuple)        
    
    
    
for filefullpath, fortfile in non_mod_files:
    new_line = "    $(DPO)" + os.path.splitext(fortfile)[0] + ".$(OBJ_EXT) \\" + "\n"
    mfix_l.write(new_line)
    mfix_u.write(new_line)
    
    
new_line =  new_line + "    $(DPO)blas90.a $(DPO)odepack.a $(DPO)dgtsv90.a" + "\n"
new_line =  new_line + "\t" + "$(LINK_CMD) $(LINK_FLAGS) " + "\\" + "\n"

mfix_l.write(new_line)
mfix_u.write(new_line)

for filefullpath, fortfile in all_fort_files:
    new_line = "    $(DPO)" + os.path.splitext(fortfile)[0] + ".$(OBJ_EXT) \\" + "\n"
    mfix_l.write(new_line)
    mfix_u.write(new_line)
    
new_line  = ' -o $(EXEC_FILE) $(LIB_FLAGS)' + "\n\n"
new_line = new_line + '$(DPO)blas90.a : $(DPO)BLAS.o' + "\n"
new_line = new_line + "\t" + "ar cr $(DPO)blas90.a $(DPO)BLAS.o" + "\n"
new_line = new_line + "$(DPO)BLAS.o : BLAS.F" + "\n"
new_line = new_line + "\t" + "$(FORTRAN_CMD) $(FORT_FLAGS) BLAS.F -o $(DPO)BLAS.o" + "\n"
new_line = new_line + "$(DPO)dgtsv90.a : $(DPO)DGTSV.o" + "\n"
new_line = new_line + "\t" + "ar cr $(DPO)dgtsv90.a $(DPO)DGTSV.o" + "\n"
new_line = new_line + "$(DPO)DGTSV.o : DGTSV.F" + "\n"
new_line = new_line + "\t" + "$(FORTRAN_CMD) $(FORT_FLAGS) DGTSV.F -o $(DPO)DGTSV.o" + "\n"
new_line = new_line + "$(DPO)odepack.a : $(DPO)ODEPACK.o" + "\n"
new_line = new_line + "\t" + "ar cr $(DPO)odepack.a $(DPO)ODEPACK.o" + "\n"
new_line = new_line + "$(DPO)ODEPACK.o : ODEPACK.F" + "\n"
new_line = new_line + "\t" + "$(FORTRAN_CMD) $(FORT_FLAGS3) ODEPACK.F -o $(DPO)ODEPACK.o" + "\n"
    
mfix_l.write(new_line)
mfix_u.write(new_line)

new_lline = ""
new_uline = ""
#now the make files will be written out to only at the end. 
for filefullpath, filename, modname in all_mod_files:
    #use_list = []
    #inc_list = []
    
    useinc_dict = call_get_use_inc_list(filefullpath)        
    use_list = useinc_dict['uselist']
    inc_list = useinc_dict['inclist']
    
    len_dep = len(use_list) + len(inc_list)    
    BACK  = " \\"
    
    if(len_dep == 0) :      
        BACK = " "
        
    new_lline  = new_lline + "$(DPO)" + modname + " : " + filefullpath + BACK + "\n"    
    modname_u = (os.path.splitext(modname)[0]).upper() + (os.path.splitext(modname)[1])     
    new_uline  = new_uline + "$(DPO)" + modname_u + " : " + filefullpath + BACK + "\n"
    
    count = 0
    for usename in use_list:
        BACK  = " \\"
        count = count + 1
        if (count == len(use_list)) & (len(inc_list) == 0) :
            BACK = " "
            
        new_lline = new_lline + "\t" + "$(DPO)" + usename + ".mod"+ BACK +  "\n"    
        new_uline = new_uline + "\t" + "$(DPO)" + usename.upper() + ".mod"+ BACK +  "\n"    
        
    count = 0
    for inc_name in inc_list:
        BACK  = " \\"
        count = count + 1
        if (count == len(inc_list)):
            BACK = " "
            
        new_lline = new_lline + "\t" + inc_name + BACK +  "\n"  
        new_uline = new_uline + "\t" + inc_name + BACK +  "\n"          
        
    #now write the file linking
    new_lline = new_lline + "\t" + "$(FORTRAN_CMD) $(FORT_FLAGS) " + filefullpath + " -o "
    new_lline = new_lline + "$(DPO)" + os.path.splitext(filename)[0] + ".$(OBJ_EXT)"
    new_lline = new_lline + " $(MODDIRPREFIX)$(DPO) " + "\n"
    
    new_uline = new_uline + "\t" + "$(FORTRAN_CMD) $(FORT_FLAGS) " + filefullpath + " -o "
    new_uline = new_uline + "$(DPO)" + os.path.splitext(filename)[0] + ".$(OBJ_EXT)"
    new_uline = new_uline + " $(MODDIRPREFIX)$(DPO) " + "\n"
    
    
for filefullpath, filename in non_mod_files:
    #use_list = []
    #inc_list = []
    objname = os.path.splitext(filename)[0] + ".$(OBJ_EXT)"
    
    useinc_dict = call_get_use_inc_list(filefullpath)        
    use_list = useinc_dict['uselist']
    inc_list = useinc_dict['inclist']
                
    len_dep = len(use_list) + len(inc_list)
    BACK  = " \\"
    if(len_dep == 0) :
        BACK = " "
        
    new_lline  = new_lline + "$(DPO)" + objname + " : " + filefullpath + BACK + "\n"
    new_uline  = new_uline + "$(DPO)" + objname + " : " + filefullpath + BACK + "\n"
    
    count = 0
    for usename in use_list:
        BACK  = " \\"
        count = count + 1
        if (count == len(use_list)) & (len(inc_list) == 0) :
            BACK = " "
            
        new_lline = new_lline + "\t" + "$(DPO)" + usename + ".mod"+ BACK +  "\n"    
        new_uline = new_uline + "\t" + "$(DPO)" + usename.upper() + ".mod"+ BACK +  "\n" 
                           
    count = 0
    for inc_name in inc_list:
        BACK  = " \\"
        count = count + 1
        if (count == len(inc_list)):
            BACK = " "
        
        new_lline = new_lline + "\t" + inc_name + BACK +  "\n"   
        new_uline = new_uline + "\t" + inc_name + BACK +  "\n"   
        
        
    common_line = "\t" + "$(FORTRAN_CMD) $(FORT_FLAGS) " + filefullpath + " -o "
    common_line = common_line + "$(DPO)" + os.path.splitext(filename)[0] + ".$(OBJ_EXT)"
    common_line = common_line + " $(MODDIRPREFIX)$(DPO) " + "\n"
    
    new_lline = new_lline + common_line
    new_uline = new_uline + common_line
  
  
mfix_l.write(new_lline)
mfix_u.write(new_uline)
     
mfix_l.close()
mfix_u.close()

files_list.close() 
mod_files_list.close()

shutil.copy(mfix_l_make_loc, MODEL_PATH)
shutil.copy(mfix_u_make_loc, MODEL_PATH)

os.remove(mfix_l_make_loc)
os.remove(mfix_u_make_loc)
os.remove(files_list_loc)
os.remove(mod_files_list_loc)
