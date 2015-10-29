#!/bin/bsh -f

# Indicates if the preprocessor is in the reaction block @(RXNS
in_rxn_block="false"      # TFM
in_des_rxn_block="false"  # DPM

# Indicates if a reaction block was found in the deck file.
fnd_rxn="false"      # TFM
fnd_des_rxn="false"  # DPM

# Number of reaction constructs found
rx_cnt="0"      # TFM
des_rx_cnt="0"  # DPM

call_rxn_stop="false"
rxn_fmt_err="false"

rxn_dup_err="false"
all_aliases="_Blank"
alias_count=0
duplicate_aliases=""

format_err() {
  echo " "
  echo "  *********************************************************************"
  echo "  * Species aliases and reaction names must be FORTRAN variable name  *"
  echo "  * compliant: alphanumeric/underscore combinations with the first    *"
  echo "  * character being a letter. Please correct the variable name in the *"
  echo "  * mfix.dat file and usr_rates files as needed.                      *"
  echo "  *********************************************************************"
  echo " "
  echo " "
  return
}

duplicate_err() {
  echo " "
  echo "  *********************************************************************"
  echo "  * Species aliases and reaction names must be UNIQUE. Two or more    *"
  echo "  * entries were found to be identical. Please correct these entries  *"
  echo "  * in the mfix.dat and usr_rates files as needed.                    *"
  echo "  *                                                                   *"

  for dup in $duplicate_aliases; do
    lmsg="  * Duplicate entry: $dup"
    lc=$((21 + ${#dup}))
    while test $lc -lt 70; do
      lmsg=$lmsg" "
      lc=$((lc+1))
    done
    echo "$lmsg*"
  done

  echo "  *********************************************************************"
  echo " "

  return
}


check_val () {
# $1: Type of variable: Species Alias, Reaction Name
# $2: Value pulled from deck file.
# $3: Length constraint

# Clean up any trailing white spaces.
  v1=$(echo $2 | sed -e 's/^ $//g')
  v2=$(echo $2 | sed -e 's/[^[:alnum:]|_]//g')
# v3 is empty if the first character is a letter.
  v3=$(echo $v1 | cut -c 1 | grep '[^a-zA-Z]')
# Verify that there are no special characters in the name.
  if [ "$v1" != "$v2" ]; then
    echo " "
    echo "  * Format Error - Invalid $1: $2"
# Verify that the first character is a letter
    call_stop="true"
    rxn_fmt_err="true"
#  elif [ ! -z "$(echo ${v1:0:1} | grep '[^a-zA-Z]')" ]; then
  elif [ ! -z "$v3" ]; then
    echo " "
    echo "  * Format Error - Invalid $1: $2"
    call_stop="true"
    rxn_fmt_err="true"
  elif [ ${#v1} -gt $3 ]; then
    echo " "
    echo "  * Format Error - $1: $2 is too long! $3 character limit."
    call_stop="true"
    rxn_fmt_err="true"
  fi

# Verify that the alias is not a duplicate.
  for fndAlias in $all_aliases; do
    if [ "$fndAlias" = "$2" ]; then
# Duplicate entries flag an error and force make_mfix to exit
      call_stop="true"
      rxn_dup_err="true"
      duplicate_aliases=$duplicate_aliases" "$2
    fi
  done

# Unique entries are addeed to an array for later comparison.
  all_aliases=$all_aliases" "$2
  alias_count=$((alias_count+1))

  return
}


# If a species.inc file exits, delete it and start a new instance.
sfile="species.inc"
if [ -f "$sfile" ]; then
  rm "$sfile"
fi
touch $sfile

echo "! This file is automatically generated by make_mfix through processing" >> $sfile
echo "! species and reaction block input. Do not directly edit this file." >> $sfile

# Check to see if there is a deck file in the run directory.
if [ ! -f "mfix.dat" ]; then
  echo "  *********************************************************************"
  echo "  * Error: usr_rates.f located and mfix.dat is missing!               *"
  echo "  * ----------------------------------------------------------------- *"
  echo "  * User defined file for chemical reactions/phase changes was found  *"
  echo "  * in the compile directory however, no mfix.dat file was found.     *"
  echo "  * Pre-processing  of the mfix.dat is required for reacting flows.   *"
  echo "  *********************************************************************"
  echo ""
  echo ""
  exit
else

# Loop over each entry in the deck file.
  while read p ; do
# Clean up some control characters (Windows)
    p=$(echo $p | tr -d '\015')
# Remove white spaces and tabs.
    p=$(echo $p | sed -e 's/^[ \t]*//;s/[ \t]*$//')
# Remove comments with hash.
    p=$(echo $p | cut -d "#" -f1 | cut -d "!" -f1)

# Pull off SPECIES_ALIAS lines.
    p1g=$(echo $p | grep -i 'species_alias_g')
    p1g=$(echo $p1g | sed 's/.*\(species_alias_g.*\)/\1/gI')
    p1g=$(echo $p1g | sed 's/species_g.*//gI')

    p1s=$(echo $p | grep -i "species_alias_s.*") # solids phases
    p1s=$(echo $p1s | sed 's/.*\(species_alias_s.*\)/\1/gI')
    p1s=$(echo $p1s | sed 's/species_s.*//gI')  

    p1s_des=$(echo $p | grep -i "des_species_alias_s.*") # solids phases
    p1s_des=$(echo $p1s_des | sed 's/.*\(des_species_alias_s.*\)/\1/gI')
    p1s_des=$(echo $p1s_des | sed 's/des_species_s.*//gI')

# Process gas phase species aliases.
#-----------------------------------------------------------------------
    if [ ! -z "$p1g" ]; then
# Get the initial species index.
      sindex=$(echo $p1g | cut -d "(" -f2 | cut -d ")" -f1)
# Get the species alias.
      salias=$(echo $p1g | cut -d "=" -f2 | tr -d "'" | tr -d "\"")
# Replace any commas with spaces.
      salias=$(echo $salias | tr "," " ")
# Initialize the loop counter
      lc="0"
# Loop over compound entries
      for p1a in $salias; do
# Ensure that the speices aliases adhead to Fortran variable naming.
        check_val "Species Alias" $p1a 32
# Calculate the species index (for compound entries)
        nsindex=$((sindex+lc))
# Write the species alias link to the include file
        echo "      INTEGER, PARAMETER :: $p1a = $nsindex" >> $sfile
# Increment the loop counter
        lc=$((lc+1))
      done
# Process solids phase species aliases. (DEM/MP-PIC solids)
# Due to similarity in names, DPM alias checks have to come before
# any TFM solids checks.
#-----------------------------------------------------------------------
    elif [ ! -z "$p1s_des" ]; then
# Get the initial species index.
      sphase=$(echo $p1s_des | cut -d "(" -f2 | cut -d "," -f1)
# Get the initial species index.
      sindex=$(echo $p1s_des | cut -d "," -f2 | cut -d ")" -f1)
# Get the species alias.
      salias=$(echo $p1s_des | cut -d "=" -f2 | tr -d "'" | tr -d "\"")
# Replace any commas with spaces.
      salias=$(echo $salias | tr "," " ")
# Initialize the loop counter
      lc="0"
# Loop over compound entries
      for p1b in $salias; do
# Ensure that the speices aliases adhead to Fortran variable naming.
        check_val "Species Alias" $p1b 32
# Calculate the species index (for compound entries)
        nsphase=$((sphase+lc))
# Write the species alias link to the include file
        echo "      INTEGER, PARAMETER :: $p1b = $sindex" >> $sfile
# Increment the loop counter
        lc=$((lc+1))
      done
# Process solids phase species aliases. (TFM solids)
#-----------------------------------------------------------------------
    elif [ ! -z "$p1s" ]; then
# Get the initial species index.
      sphase=$(echo $p1s | cut -d "(" -f2 | cut -d "," -f1)
# Get the initial species index.
      sindex=$(echo $p1s | cut -d "," -f2 | cut -d ")" -f1)
# Get the species alias.
      salias=$(echo $p1s | cut -d "=" -f2 | tr -d "'" | tr -d "\"")
# Replace any commas with spaces.
      salias=$(echo $salias | tr "," " ")
# Initialize the loop counter
      lc="0"
# Loop over compound entries
      for p1b in $salias; do
# Ensure that the speices aliases adhead to Fortran variable naming.
        check_val "Species Alias" $p1b 32
# Calculate the species index (for compound entries)
        nsphase=$((sphase+lc))
# Write the species alias link to the include file
        echo "      INTEGER, PARAMETER :: $p1b = $sindex" >> $sfile
# Increment the loop counter
        lc=$((lc+1))
      done
    fi

# Pull off reaction blocks
#-----------------------------------------------------------------------
    p2END=$(echo $p | grep -i "@(END")
    p2TFM=$(echo $p | grep -i "@(RXNS")      # TFM
    p2DPM=$(echo $p | grep -i "@(DES_RXNS")  # DPM (DES/MP-PIC)

# Find the end of a reaction block.
    if [ ! -z "$p2END" ]; then
      if [ "$in_rxn_block" = "true" ]; then
        in_rxn_block="false"
      elif [ "$in_des_rxn_block" = "true" ]; then
        in_des_rxn_block="false"
      else
        echo "  Found the end of a reaction block without finding start"
        echo "  Fatal Error: Killing make_mfix"
        exit
      fi
    elif [ ! -z "$p2TFM" ]; then
      if [ "$in_des_rxn_block" = "true" ]; then
        echo "  Found the start of a TFM reaction block while"
        echo "  still inside a DPM reaction block!"
        echo "  Fatal Error: Killing make_mfix"
        exit
      else
        fnd_rxn="true"
        in_rxn_block="true"
      fi
    elif [ ! -z "$p2DPM" ]; then
      if [ "$in_rxn_block" = "true" ]; then
        echo "  Found the start of a DPM reaction block while"
        echo "  still inside a TFM reaction block!"
        echo "  Fatal Error: Killing make_mfix"
        exit
      else
        fnd_des_rxn="true"
        in_des_rxn_block="true"
      fi

# If already processing a TFM reaction block.
    elif [ "$in_rxn_block" = "true" ]; then
      if [ ! -z "$(echo $p | grep "{")" ]; then
# Take the name that precedes the start of the reaction construct {
        rx_name=$(echo $p | cut -d "{" -f1)
        check_val "Reaction Name" $rx_name 32
        rx_cnt=$((rx_cnt+1))
        echo "      INTEGER, PARAMETER :: $rx_name = $rx_cnt" >> $sfile
      fi
    elif [ "$in_des_rxn_block" = "true" ]; then
      if [ ! -z "$(echo $p | grep "{")" ]; then
# Take the name that precedes the start of the reaction construct {
        rx_name=$(echo $p | cut -d "{" -f1)
        check_val "Reaction Name" $rx_name 32
        des_rx_cnt=$((des_rx_cnt+1))
        echo "      INTEGER, PARAMETER :: $rx_name = $des_rx_cnt" >> $sfile
      fi
    fi
  done < "mfix.dat"
fi

# If there are no continuum phase reactions and there are no discrete
# phase reactions, delete sfile.
if test "$fnd_rxn" = "false" && test "$fnd_des_rxn" = "false" ; then
  rm "$sfile"
fi

# Report any format errors
if test "$rxn_fmt_err" = "true"; then
  format_err
fi
# Report any duplicate entries (species aliases/reaction names)
if test "$rxn_dup_err" = "true"; then
  duplicate_err
fi

if test "$call_stop" = "true"; then
  echo " "
  echo " "
  echo "  An input error for chemical reactions was detected."
  echo "  Please correct the error and execute make_mfix again."
  echo "  Exiting make_mfix."
  echo " "
  echo " "
  exit
fi

unset IFS

echo " "
echo "  Reaction data was successfully processed. "
echo " "
