#!/bin/bash
#
# This is a low level example of a script to run GENIE using
# captControl.  It uses the gevgen_capt.exe main program for GENIE
# which is lightly customized from gevgen.  

# The number of events to generate
EVENTS=5000

# The output file name.
GHEP_RUN=1
GHEP_PREFIX=captain-genie
GHEP_FILE="${GHEP_PREFIX}.${GHEP_RUN}.ghep.root"

# The flux to use.  GENIE assumes the energy is in units of GeV, so
# the peak energy is 5 GeV.
FLUX="text,14,../flux/fluka08_me000z200i_810km_0kmoa_flux.txt,0,2"
FLUX="${FLUX}:text,-14,../flux/fluka08_me000z200i_810km_0kmoa_flux.txt,0,5"

 gevgen_capt.exe -r ${GHEP_RUN} \
     -n ${EVENTS} \
     -f ${FLUX} 
  
gntpc -f rootracker -i ${GHEP_FILE} -o captain-genie.root

