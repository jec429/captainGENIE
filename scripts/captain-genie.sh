#!/bin/bash
#
# This is a low level example of a script to run GENIE using
# captControl.  It uses the gevgen_capt.exe main program for GENIE
# which is lightly customized from gevgen.  

source captain-control.bash

# Set the run number.
captain-experiment mc
captain-data-source gn
captain-run-type spl
captain-run-number 1

# The number of events to generate
EVENTS=10

# The output file name.
GHEP_PREFIX=$(basename $(captain-file "ghep") ".root")
GHEP_FILE="${GHEP_PREFIX}.$(captain-run-number).ghep.root"

echo $GHEP_PREFIX
echo $GHEP_FILE

# The flux to use.  GENIE assumes the energy is in units of GeV, so
# the peak energy is 5 GeV.
FLUX="x*x*exp(-(x/5.0)**2)"

# The neutrino flavor (PDG).  The default target is natural argon 
PDG=14

# Limit the amount of noise from GENIE
LOGLEVEL=${GENIE}/config/Messenger_laconic.xml

gevgen_capt.exe -r $(captain-run-number) \
    -o ${GHEP_PREFIX} \
    -n ${EVENTS} \
    -e 0.1,15.0 -p ${PDG} -f ${FLUX} \
    --message-thresholds ${LOGLEVEL}
gntpc -f rootracker -i ${GHEP_FILE} -o $(captain-file "gnmc") \
    --message-thresholds ${LOGLEVEL}

mv $GHEP_FILE $(captain-file "ghep")

# Write a GEANT4 macro file to process the output.
G4_MACRO=$(captain-file "g4in" "mac")
cat >> ${G4_MACRO} <<EOF
/dsim/control baseline 1.0
/dsim/update

/generator/kinematics/rooTracker/input $(captain-file "gnmc")
/generator/kinematics/set rooTracker

# Have exactly one interaction per event.
/generator/count/fixed/number 1
/generator/count/set fixed

# Choose the position based on the density (and only in the drift volume).
/generator/position/density/volume Drift
/generator/position/set density

/generator/add

/run/beamOn ${EVENTS}
EOF

#####################################################
# The is the meat of the script: Run the DETSIM, ELECSIM, calibration,
# and reconstruction.  The file names are generated based on the 
# commands above.
#####################################################
captain-process-detsim-macro $G4_MACRO
# captain-run-reconstruction
