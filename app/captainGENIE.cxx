//____________________________________________________________________________
/*!

\program gevgen_capt

\brief A simple GENIE v+A event generation driver (gevgen_capt) customized for
         CAPTAIN.  This needs to use the internal genie command line argument
         parser because it needs access to the built in arguments, but
         unfortunately, the parser is junk, so the command line controls are a
         kludge. 

         Syntax :
           gevgen_capt [-h] 
                   -n nev 
                  [-r run#] 
                  [-t target_pdg]
                  [-d directionCosines]
                  [-f flux_description] 
                  [-o file_prefix]
                  [--seed random_number_seed] 
                  [--cross-sections xml_file]
                  [--event-generator-list list_name]
                  [--message-thresholds xml_file]          
                  [--unphysical-event-mask mask]
                  [--event-record-print-level level]
                  [--mc-job-status-refresh-rate  rate]
                  [--cache-file root_file]

         Options :
           [] Denotes an optional argument.
           -h 
              Prints-out help on using gevgen_capt and exits.
           -n 
              Specifies the number of events to generate.
           -r 
              Specifies the MC run number.
           -o 
              Set the output file prefix (default is "captain-genie")
           -d
              Specify the neutrino direction cosines.  For example 
              '-d 0.1,0.2,0.3'.  The result will be normalized to a unit
              vector.  The default is (1,0,0)
           -t 
              Specifies the target PDG code (pdg format: 10LZZZAAAI) _or_ a
              target mix (pdg codes with corresponding weights) typed as a
              comma-separated list of pdg codes with the corresponding weight
              fractions in brackets, eg code1[fraction1],code2[fraction2],...
              For example, to use a target mix of 95% O16 and 5% H type: `-t
              1000080160[0.95],1000010010[0.05]'.
           -f 
              Specifies the neutrino flux spectrum.
              It can be any of:
              -- An energy in GeV.  This can't be combined with other flux
                 descriptions. eg ` -f mono,pdg,energy' 
              -- A function (x in GeV):
                 eg ` -f func,pdg,x*x+4*exp(-x),emin,emax' 
              -- A 1-D ROOT histogram (binned in GeV):
                 The general syntax is `-f hist,pdg,file.root,hist_name'
              -- A text file with columns of energy (in GeV) and flux:
                 The general syntax is `-f text,pdg,file.txt,Ecol,Fcol' with
                 ECol giving the column with the energy Fcol giving the column
                 with the flux.  The colunms are counted from 0.
              If multiple fluxes are specified, then they are separated by a
              ':', `-f hist,14,file.root,numu:hist,-14,file.root,numubar'
              would read the muon neutrino (pdg 14) flux from the histogram
              named numu in file.root, and the muon anti-neutrino (pdg -14)
              flux from the histogram named numubar in file.root.  This clunky
              interface is forced by the GENIE command line argument parser
              which can't handle multiple copies of an option.
           --seed
              Random number seed.  A value of -1 uses a default seed value. 
              A value of 0 uses a UUID based seed.
           --cross-sections
              Name (incl. full path) of an XML file with pre-computed
              cross-section values used for constructing splines.
           --event-generator-list            
              List of event generators to load in event generation drivers.
              [default: "Default"].
           --message-thresholds           
              Allows users to customize the message stream thresholds.
              The thresholds are specified using an XML file.
              See $GENIE/config/Messenger.xml for the XML schema.
           --unphysical-event-mask       
              Allows users to specify a 16-bit mask to allow certain types of
              unphysical events to be written in the output file.
              [default: all unphysical events are rejected]
           --event-record-print-level
              Allows users to set the level of information shown when the event
              record is printed in the screen. See GHepRecord::Print().
           --mc-job-status-refresh-rate   
              Allows users to customize the refresh rate of the status file.
           --cache-file                  
              Allows users to specify a cache file so that the cache can be
              re-used in subsequent MC jobs.

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory
\author  Modified by Clark McGrew <clark.mcgrew \at stonybrook.eduk>
         Stony Brook Univ.

\created October 05, 2014

\cpright Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
         GPL3: For the full text of the license visit 
         http://copyright.genie-mc.org or see $GENIE/LICENSE 
*/
//____________________________________________________________________________

#include <cstdlib>
#include <cassert>
#include <sstream>
#include <string>
#include <vector>
#include <map>

#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TVector3.h>
#include <TRandom.h>
#include <TH1D.h>
#include <TF1.h>
#include <Math/Interpolator.h>

#include "Conventions/XmlParserStatus.h"
#include "Conventions/GBuild.h"
#include "Conventions/GVersion.h"
#include "Conventions/Controls.h"
#include "EVGCore/EventRecord.h"
#include "EVGDrivers/GFluxI.h"
#include "EVGDrivers/GEVGDriver.h"
#include "EVGDrivers/GMCJDriver.h"
#include "EVGDrivers/GMCJMonitor.h"
#include "Interaction/Interaction.h"
#include "Ntuple/NtpWriter.h"
#include "Ntuple/NtpMCFormat.h"
#include "Numerical/RandomGen.h"
#include "Numerical/Spline.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "Utils/AppInit.h"
#include "Utils/RunOpt.h"
#include "Utils/XSecSplineList.h"
#include "Utils/StringUtils.h"
#include "Utils/PrintUtils.h"
#include "Utils/SystemUtils.h"
#include "Utils/CmdLnArgParser.h"

#include "FluxDrivers/GCylindTH1Flux.h"
#include "FluxDrivers/GMonoEnergeticFlux.h"
#include "Geo/PointGeomAnalyzer.h"

// User-specified options:
int             gOptNevents = 0;  // n-events to generate
std::string     gOptEvFilePrefix = "captain-genie"; // event file prefix
std::map<int,double> gOptTgtMix;  // target mix (each with its relative weight)
Long_t          gOptRunNu = 0.0;  // run number
std::string     gOptFlux;         // the flux to use.
long int        gOptRanSeed;      // random number seed
std::string     gOptInpXSecFile;  // cross-section splines
TVector3        gOptNuDir(1,0,0); // Neutrino direction in genie.

void Initialize() {
    // Initialization of random number generators, cross-section table, 
    // messenger thresholds, cache file
    genie::utils::app_init::MesgThresholds(
        genie::RunOpt::Instance()->MesgThresholdFiles());
    genie::utils::app_init::CacheFile(genie::RunOpt::Instance()->CacheFile());
    genie::utils::app_init::RandGen(gOptRanSeed);
    genie::utils::app_init::XSecTable(gOptInpXSecFile, false);
    
    // Set GHEP print level
    genie::GHepRecord::SetPrintLevel(
        genie::RunOpt::Instance()->EventRecordPrintLevel());
}

genie::GeomAnalyzerI * GeomDriver(void) {
    genie::GeomAnalyzerI * geom_driver 
        = new genie::geometry::PointGeomAnalyzer(gOptTgtMix);
    return geom_driver;
}

genie::GFluxI* MonoEnergeticFluxDriver(std::string description) {
    std::vector<std::string> splitFlux 
        = genie::utils::str::Split(description,",");

    // Check if this is a monoenergetic flux.
    if (splitFlux.size() > 0 || splitFlux[0] != "mono")  return NULL;

    if (splitFlux.size() != 3) {
        std::cout << "Invalid monoenergetic flux description: " << description
                  << std::endl;
        std::exit(1);
    }

    int pdg = std::atoi(splitFlux[1].c_str());
    double energy = std::atoi(splitFlux[2].c_str());

    std::cout << "Generate monoenergetic flux for " << pdg
              << " at " << energy << " GeV" << std::endl;

    genie::flux::GMonoEnergeticFlux* flux = 
        new genie::flux::GMonoEnergeticFlux(energy, pdg);

    TVector3 bspot(0,0,0);
    flux->SetNuDirection(gOptNuDir);
    flux->SetBeamSpot(bspot);

    return dynamic_cast<genie::GFluxI *>(flux);
}

// Read a histogram from a root file.
TH1D* HistogramFlux(int pdg, std::string file, std::string name) {
    std::cout << "Read histogram for " << pdg
              << " from " << file << " : " << name << std::endl;

    if (!gSystem->AccessPathName(file.c_str())) {
        std::cout << "Invalid flux file: " << file
                  << std::endl;
        std::exit(1);
    }
    
    // Open the file, get the flux histogram, and neutrino type.
    TFile flux_file(file.c_str(),"read");
    TH1* hist = dynamic_cast<TH1*>(flux_file.Get(name.c_str()));

    if (!hist) {
        std::cout << "Histogram not found: " 
                  << file << " " << name
                  << std::endl;
        std::exit(1);
    }

    std::ostringstream histName;
    histName << "fluxHistogram" << pdg;
    
    // copy all bins between emin,emax
    TH1D* spectrum  = new TH1D(histName.str().c_str(),"neutrino flux",
                               hist->GetNbinsX(),
                               hist->GetXaxis()->GetXmin(),
                               hist->GetXaxis()->GetXmax());

    for(int ibin = 1; ibin <= hist->GetNbinsX(); ibin++) {
        spectrum->SetBinContent(ibin, hist->GetBinContent(ibin));
    }
    
    flux_file.Close();        
    
    return spectrum;
}

// ** Generate the flux histogram from the input functional form
TH1D* FunctionalFlux(int pdg, std::string func, 
                     int bins, double emin, double emax) {
    std::cout << "Generate flux for " << pdg
              << " from " << func 
              << " w/ " << bins << " bins"
              << " b/w " << emin << " GeV and " << emax << " GeV"
              << std::endl;

    std::ostringstream histName;
    histName << "fluxHistogram" << pdg;
    TF1 input_func("input_func", func.c_str(), emin, emax);
    TH1D* spectrum  = new TH1D(histName.str().c_str(),func.c_str(), 
                               bins, emin, emax);
    for (int i=1; i<bins+1; ++i) {
        double x = spectrum->GetBinCenter(i);
        double f = input_func.Eval(x);
        spectrum->SetBinContent(i, f);
    }
    
    return spectrum;
}

TH1D* TextFlux(int pdg, std::string file, int energyColumn, int fluxColumn) {
    std::ifstream input(file.c_str());

    std::vector<double> energy;
    std::vector<double> flux;
    std::string line;
    while (std::getline(input,line)) {
        line = line.substr(0,line.find("#"));
        std::istringstream inputLine(line);
        double val;
        int column = 0;
        while (inputLine >> val) {
            if (column == energyColumn) {
                energy.push_back(val);
            }
            if (column == fluxColumn) {
                flux.push_back(val);
            }
            ++column;
        }
    }


    std::cout << "Generate flux for " << pdg
              << " from " << file 
              << " w/ " << energy.size() << " bins"
              << std::endl;

    if (flux.size() != energy.size()) {
        std::cout << "Mismatched number of energy and flux bins " << std::endl;
        std::exit(1);
    }

    std::ostringstream histName;
    histName << "fluxHistogram" << pdg;
    double delta = 1.0*(energy.back() - energy.front())/(energy.size()-1);
    TH1D* spectrum  = new TH1D(histName.str().c_str(),
                               file.c_str(),
                               energy.size(), 
                               energy.front()-0.5*delta, 
                               energy.back()+0.5*delta);

    int bins = spectrum->GetNbinsX();
    ROOT::Math::Interpolator func(energy,flux);
    for (int i=1; i<bins+1; ++i) {
        double x = spectrum->GetBinCenter(i);
        spectrum->SetBinContent(i+1, func.Eval(x));
    }

    return spectrum;
}

// create & configure one of the generic flux drivers
genie::GFluxI * FluxDriver(void) {

    std::vector<std::string> splitOption 
        = genie::utils::str::Split(gOptFlux,":");

    if (splitOption.size() == 1) {
        // Try the monoenergetic flux driver.  This is a special case, and
        // there can only be one.  Just return the created flux driver.
        genie::GFluxI* fluxDriver = MonoEnergeticFluxDriver(splitOption[0]);
        if (fluxDriver) return fluxDriver;
    }

    // All other flux types are described by a histogram.
    genie::flux::GCylindTH1Flux* fluxDriver = new genie::flux::GCylindTH1Flux;
    fluxDriver->SetNuDirection(gOptNuDir);
    TVector3 bspot(0,0,0);
    fluxDriver->SetBeamSpot(bspot);
    fluxDriver->SetTransverseRadius(-1);

    for (std::vector<std::string>::iterator s = splitOption.begin();
         s != splitOption.end(); ++s) {
        std::vector<std::string> splitFlux 
            = genie::utils::str::Split(*s,",");
        if (splitFlux.size() < 2) {
            std::cout << "Invalid flux specification: " << *s << std::endl;
            std::exit(1);
        }
        if (splitFlux[0] == "mono") {
            std::cout << "Cannot mix monoenergetic fluxes with other types."
                      << gOptFlux << std::endl;
            std::exit(1);
        }
        if (splitFlux[0] == "func") {
            if (splitFlux.size() != 6) {
                std::cout << "A functional flux needs 6 arguments."
                          << *s << std::endl;
                std::exit(1);
            }
            int pdg = std::atoi(splitFlux[1].c_str());
            int bins = std::atoi(splitFlux[3].c_str());
            double emin = std::atof(splitFlux[4].c_str());
            double emax = std::atof(splitFlux[5].c_str());
            TH1D* flux = FunctionalFlux(pdg,splitFlux[2], 
                                        bins, emin, emax);
            fluxDriver->AddEnergySpectrum(pdg, flux);
            continue;
        }
        if (splitFlux[0] == "hist") {
            if (splitFlux.size() != 4) {
                std::cout << "A histogram flux needs 4 arguments."
                          << *s << std::endl;
                std::exit(1);
            }
            int pdg = std::atoi(splitFlux[1].c_str());
            TH1D* flux = HistogramFlux(pdg, splitFlux[2], splitFlux[3]);
            fluxDriver->AddEnergySpectrum(pdg, flux);
            continue;
        }
        if (splitFlux[0] == "text") {
            if (splitFlux.size() != 5) {
                std::cout << "Functional flux needs 5 arguments."
                          << *s << std::endl;
                std::exit(1);
            }
            int pdg = std::atoi(splitFlux[1].c_str());
            int energyColumn = std::atoi(splitFlux[3].c_str());
            int fluxColumn = std::atoi(splitFlux[4].c_str());
            TH1D* flux = TextFlux(pdg,splitFlux[2], 
                                  energyColumn, fluxColumn);
            fluxDriver->AddEnergySpectrum(pdg, flux);
            continue;
        }
    }

    return fluxDriver;
}

void PrintSyntax(void) {
    std::cout << std::endl
              << std::endl
              << "Syntax:" 
              << "      gevgen_capt [-h]"
              << std::endl
              << "              [-r run#]"
              << std::endl 
              << "               -n nev"
              << std::endl 
              << "              [-t target_pdg] "
              << std::endl 
              << "              [-d dCosX,dCosY,dCosZ] (default: 1,0,0)"
              << std::endl 
              << "              [-o file_prefix]"
              << std::endl 
              << "              [-f flux_description]"
              << std::endl 
              << "              [-w]"
              << std::endl 
              << "              [--seed random_number_seed]"
              << std::endl 
              << "              [--cross-sections xml_file]"
              << std::endl
              << "              [--event-generator-list list_name]"
              << std::endl
              << "              [--message-thresholds xml_file]"
              << std::endl
              << "              [--unphysical-event-mask mask]"
              << std::endl 
              << "              [--event-record-print-level level]"
              << std::endl 
              << "              [--mc-job-status-refresh-rate  rate]"
              << std::endl
              << "              [--cache-file root_file]"
              << std::endl;
}

void GetCommandLineArgs(int argc, char ** argv) {
    // Common run options. Set defaults and read.
    genie::RunOpt::Instance()->EnableBareXSecPreCalc(true);
    genie::RunOpt::Instance()->ReadFromCommandLine(argc,argv);
    
    // Parse run options for this app
    genie::CmdLnArgParser parser(argc,argv);
    
    // help?
    bool help = parser.OptionExists('h');
    if (help) {
        PrintSyntax();
        exit(0);
    }
    
    // number of events
    if (parser.OptionExists('n')) {
        gOptNevents = parser.ArgAsInt('n');
    } 
    else {
        std::cout << "Must specify number of events to generate" << std::endl;
        PrintSyntax();
        exit(1);
    }

    // run number
    if (parser.OptionExists('r')) {
        gOptRunNu = parser.ArgAsLong('r');
    } 

    // event file prefix
    if (parser.OptionExists('o')) {
        gOptEvFilePrefix = parser.ArgAsString('o');
    } 

    // flux functional form
    if (parser.OptionExists('f')) {
        gOptFlux = parser.ArgAsString('f');
    }

    // target mix (their PDG codes with their corresponding weights)
    if (parser.OptionExists('t')) {
        std::string stgtmix = parser.ArgAsString('t');
        gOptTgtMix.clear();
        std::vector<std::string> tgtmix = genie::utils::str::Split(stgtmix,",");
        if (tgtmix.size()==1) {
            int    pdg = atoi(tgtmix[0].c_str());
            double wgt = 1.0;
            gOptTgtMix.insert(std::map<int, double>::value_type(pdg, wgt));
        } 
        else {
            for (std::vector<std::string>::const_iterator tgtmix_iter 
                     = tgtmix.begin(); 
                 tgtmix_iter != tgtmix.end(); ++tgtmix_iter) {
                std::string tgt_with_wgt = *tgtmix_iter;
                std::string::size_type open_bracket  = tgt_with_wgt.find("[");
                std::string::size_type close_bracket = tgt_with_wgt.find("]");
                std::string::size_type ibeg = 0;
                std::string::size_type iend = open_bracket;
                std::string::size_type jbeg = open_bracket+1;
                std::string::size_type jend = close_bracket-1;
                int    pdg = atoi(tgt_with_wgt.substr(ibeg,iend).c_str());
                double wgt = atof(tgt_with_wgt.substr(jbeg,jend).c_str());
                std::cout << "Adding to target mix: pdg = " << pdg 
                          << ", wgt = " << wgt
                          << std::endl;
                gOptTgtMix.insert(std::map<int, double>::value_type(pdg, wgt));
            }
        }
    }
    else {
        std::cout << "Simulating natural argon" << std::endl;
        gOptTgtMix[1000180400] = 0.9960;
        gOptTgtMix[1000180360] = 0.0034;
        gOptTgtMix[1000180380] = 0.0006;
    }

    // Set the neutrino direction.
    if (parser.OptionExists('d')) {
        std::vector<double> dcos = parser.ArgAsDoubleTokens('d',",");
        if (dcos.size() != 3) {
            std::cout << "Wrong number of direction cosines";
            exit(1);
        }
        gOptNuDir = TVector3(dcos[0], dcos[1], dcos[2]).Unit();
    }

    // random number seed
    if (parser.OptionExists("seed")) {
        gOptRanSeed = parser.ArgAsLong("seed");
        std::cout << "Setting random number seed to " 
                  << gOptRanSeed 
                  << std::endl;
    }
    else {
        gOptRanSeed = 0;
        std::cout << "Unspecified random number seed - Using time based seed."
                  << std::endl;
    }

    // input cross-section file
    if (parser.OptionExists("cross-sections")) {
        gOptInpXSecFile = parser.ArgAsString("cross-sections");
        std::cout << "Reading cross-section file"
                  << gOptInpXSecFile
                  << std::endl;
    } 
    else {
        const char* genie = gSystem->Getenv("GENIE");
        if (!genie) {
            std::cout << "GENIE source location not set.";
            std::exit(1);
        }
        gOptInpXSecFile = std::string(genie) 
            + "/data/cross-section/"
            + "gxspl-t2k-v" + __GENIE_RELEASE__ + ".xml";
        std::cout << "Reading cross-section file"
                  << gOptInpXSecFile
                  << std::endl;
    }

    //
    // print-out the command line options
    //
    std::cout << std::endl
              << genie::utils::print::PrintFramedMesg(
                  "gevgen_capt job configuration")
              << std::endl;

    std::cout << "MC Run Number: " << gOptRunNu << std::endl;
    if (gOptRanSeed != -1) {
        if (gOptRanSeed == 0) {
            gRandom->SetSeed(0);
            gOptRanSeed = gRandom->GetSeed();
        }
        std::cout << "Random number seed: " << gOptRanSeed 
                  << std::endl;
    } 
    else {
        std::cout << "Random number set to default value" 
                  << std::endl;
    }

    std::cout << "Number of events requested: " << gOptNevents 
              << std::endl;

    if (gOptInpXSecFile.size() > 0) {
        std::cout << "Using cross-section splines read from: " 
                  << gOptInpXSecFile
                  << std::endl;
    }
    else {
        std::cout << "No input cross-section spline file" << std::endl;
    }

    
    std::vector<std::string> fluxDescriptions 
        = genie::utils::str::Split(gOptFlux,":");
    for (std::vector<std::string>::iterator s = fluxDescriptions.begin();
         s != fluxDescriptions.end(); ++ s) {
        std::cout << "Flux: " << *s << std::endl;
    }

    std::cout << "Target code (PDG) & weight fraction:"
              << std::endl;

    std::map<int,double>::const_iterator iter;
    for(iter = gOptTgtMix.begin(); iter != gOptTgtMix.end(); ++iter) {
        int    tgtpdgc = iter->first;
        double wgt     = iter->second;
        std::cout << " >> " <<  tgtpdgc << " (weight fraction = " << wgt << ")"
                  << std::endl;
    }

}

void GenerateEvents(void) {
    // Get flux and geom drivers
    genie::GFluxI *        flux_driver = FluxDriver();
    genie::GeomAnalyzerI * geom_driver = GeomDriver();
    
    // Create the monte carlo job driver
    genie::GMCJDriver * mcj_driver = new genie::GMCJDriver;
    mcj_driver->SetEventGeneratorList(
        genie::RunOpt::Instance()->EventGeneratorList());
    mcj_driver->SetUnphysEventMask(
        *genie::RunOpt::Instance()->UnphysEventMask());
    mcj_driver->UseFluxDriver(flux_driver);
    mcj_driver->UseGeomAnalyzer(geom_driver);
    mcj_driver->Configure();
    mcj_driver->UseSplines();
    mcj_driver->ForceSingleProbScale();
    
    // Initialize an Ntuple Writer to save GHEP records into a TTree
    genie::NtpWriter ntpw(genie::kNFGHEP, gOptRunNu);
    ntpw.CustomizeFilenamePrefix(gOptEvFilePrefix);
    ntpw.Initialize();
    
    // Create an MC Job Monitor
    genie::GMCJMonitor mcjmonitor(gOptRunNu);
    mcjmonitor.SetRefreshRate(
        genie::RunOpt::Instance()->MCJobStatusRefreshRate());

    // Generate events / print the GHEP record / add it to the ntuple
    int ievent = 0;
    while ( ievent < gOptNevents) {
         // generate a single event for neutrinos coming from the specified
         // flux
        genie::EventRecord * event = mcj_driver->GenerateEvent();
         
        // add event at the output ntuple, refresh the mc job monitor and
        // clean-up
        ntpw.AddEventRecord(ievent, event);
        mcjmonitor.Update(ievent,event);
        delete event;
        ++ievent;
    }
    
    // Save the generated MC events
    ntpw.Save();
    
    delete flux_driver;
    delete geom_driver;
    delete mcj_driver;;
}

int main(int argc, char ** argv) {
    GetCommandLineArgs(argc,argv);

    Initialize();

    //
    // Generate neutrino events
    //
    GenerateEvents();

    return 0;
}


// Local Variables:
// mode:c++
// c-basic-offset: 4
// End:
