#!/bin/bash
#
# This is a low level example of a script to run GENIE using
# captControl.  It uses the gevgen_capt.exe main program for GENIE
# which is lightly customized from gevgen.  

# The number of events to generate
EVENTS=5

# The output file name.
GHEP_RUN=1
GHEP_PREFIX=captain-genie
GHEP_FILE="captain-genie.${GHEP_RUN}.ghep.root"

echo $GHEP_PREFIX
echo $GHEP_FILE

# The flux to use.  GENIE assumes the energy is in units of GeV, so
# the peak energy is 5 GeV.
FLUX='x*x*exp(-(x/5.0)**2)'

# The neutrino flavor (PDG).  The default target is natural argon 
PDG=14

# Limit the amount of noise from GENIE
LOGLEVEL=${GENIE}/config/Messenger_laconic.xml

gevgen_capt.exe -r ${GHEP_RUN} \
    -o ${GHEP_PREFIX} \
    -n ${EVENTS} \
    --seed 0 \
    -e 0.1,15.0 -p ${PDG} -f ${FLUX} \
    --event-record-print-level 0 \
    --message-thresholds ${LOGLEVEL}

gntpc -f rootracker -i ${GHEP_FILE} -o captain-genie.root \
    --message-thresholds ${LOGLEVEL}

