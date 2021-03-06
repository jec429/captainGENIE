//____________________________________________________________________________
/*!

\program gevgen_capt

\brief   A simple 'generic' GENIE v+A event generation driver (gevgen_capt).

         It handles:
         a) event generation for a fixed init state (v+A) at fixed energy, or
         b) event generation for simple fluxes (specified either via some
            functional form, tabular text file or a ROOT histogram) and for 
            simple 'geometries' (a target mix with its corresponding weights)

         For more complex event generation cases, for specific experiments, 
         using the outputs of detailed neutrino flux simulations and realistic 
         detector geometry descriptions see in $GENIE/src/support/.

         Syntax :
           gevgen_capt [-h] 
                   -n nev 
                   -e energy (or energy range) 
                   -p neutrino_pdg 
                  [-r run#] 
                  [-t target_pdg]
                  [-d directionCosines]
                  [-f flux_description] 
                  [-w] 
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
              Set the output file prefix (default is "gntp")
           -d
              Specify the neutrino direction cosines.  For example 
              '-d 0.1,0.2,0.3'.  The result will be normalized to a unit
              vector.
           -e 
              Specifies the neutrino energy.
              If what follows the -e option is a comma separated pair of values
              it will be interpreted as an energy range for the flux specified
              via the -f option (see below).
           -p 
              Specifies the neutrino PDG code.
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
              -- A function:
                 eg ` -f x*x+4*exp(-x)' 
              -- A vector file:
                 The vector file should contain 2 columns corresponding to 
                 energy,flux (see $GENIE/data/flux/ for few examples). 
              -- A 1-D ROOT histogram (TH1D):
                 The general syntax is `-f /full/path/file.root,object_name'
           -w 
              Forces generation of weighted events.
              This option is relevant only if a neutrino flux is specified.
              Note that 'weighted' refers to the selection of the primary
              flux neutrino + target that were forced to interact. A weighting
              scheme for the generated kinematics of individual processes can
              still be in effect if enabled..
              ** Only use that option if you understand what it means **
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

        ***  See the User Manual for more details and examples. ***

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory

\created October 05, 2004

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
#include <TH1.h>
#include <TF1.h>

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
#include "Messenger/Messenger.h"
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

#ifdef __GENIE_FLUX_DRIVERS_ENABLED__
#ifdef __GENIE_GEOM_DRIVERS_ENABLED__
#define __CAN_GENERATE_EVENTS_USING_A_FLUX_OR_TGTMIX__
#include "FluxDrivers/GCylindTH1Flux.h"
#include "FluxDrivers/GMonoEnergeticFlux.h"
#include "Geo/PointGeomAnalyzer.h"
#endif
#endif

using std::string;
using std::vector;
using std::map;
using std::ostringstream;

using namespace genie;
using namespace genie::controls;

void GetCommandLineArgs (int argc, char ** argv);
void Initialize         (void);
void PrintSyntax        (void);

#ifdef __CAN_GENERATE_EVENTS_USING_A_FLUX_OR_TGTMIX__
void            GenerateEventsUsingFluxOrTgtMix();
GeomAnalyzerI * GeomDriver              (void);
GFluxI *        FluxDriver              (void);
GFluxI *        MonoEnergeticFluxDriver (void);
GFluxI *        TH1FluxDriver           (void);
#endif

void GenerateEventsAtFixedInitState (void);

//Default options (override them using the command line arguments):
int           kDefOptNevents   = 0;         // n-events to generate
string        kDefOptEvFilePrefix = "gntp"; // default file prefix
NtpMCFormat_t kDefOptNtpFormat = kNFGHEP;   // ntuple format
Long_t        kDefOptRunNu     = 0;         // default run number

//User-specified options:
int             gOptNevents;      // n-events to generate
string          gOptEvFilePrefix; // event file prefix
double          gOptNuEnergy;     // neutrino E, or min neutrino E in spectrum
double          gOptNuEnergyRange;// energy range in input spectrum
int             gOptNuPdgCode;    // neutrino PDG code
map<int,double> gOptTgtMix;       // target mix (each with its relative weight)
Long_t          gOptRunNu;        // run number
string          gOptFlux;         // the flux to use.
bool            gOptWeighted;     
bool            gOptUsingFluxOrTgtMix = false;
long int        gOptRanSeed;      // random number seed
string          gOptInpXSecFile;  // cross-section splines
TVector3        gOptNuDir(1,0,0); // Neutrino direction in genie.

//____________________________________________________________________________
int main(int argc, char ** argv)
{
  GetCommandLineArgs(argc,argv);
  Initialize();

  //
  // Generate neutrino events
  //
  if(gOptUsingFluxOrTgtMix) {
    GenerateEventsUsingFluxOrTgtMix();
  } 
  else {
    GenerateEventsAtFixedInitState();
  }

  return 0;
}

//____________________________________________________________________________
void Initialize()
{
  // Initialization of random number generators, cross-section table, 
  // messenger thresholds, cache file
  utils::app_init::MesgThresholds(RunOpt::Instance()->MesgThresholdFiles());
  utils::app_init::CacheFile(RunOpt::Instance()->CacheFile());
  utils::app_init::RandGen(gOptRanSeed);
  utils::app_init::XSecTable(gOptInpXSecFile, false);

  // Set GHEP print level
  GHepRecord::SetPrintLevel(RunOpt::Instance()->EventRecordPrintLevel());
}

//____________________________________________________________________________
void GenerateEventsAtFixedInitState(void)
{
  int neutrino = gOptNuPdgCode;
  int target   = gOptTgtMix.begin()->first;
  double Ev    = gOptNuEnergy;
  TLorentzVector nu_p4(0.,0.,Ev,Ev); // px,py,pz,E (GeV)

  // Create init state
  InitialState init_state(target, neutrino);

  // Create/config event generation driver 
  GEVGDriver evg_driver;
  evg_driver.SetEventGeneratorList(RunOpt::Instance()->EventGeneratorList());
  evg_driver.SetUnphysEventMask(*RunOpt::Instance()->UnphysEventMask());
  evg_driver.Configure(init_state);

  // Initialize an Ntuple Writer
  NtpWriter ntpw(kDefOptNtpFormat, gOptRunNu);
  ntpw.CustomizeFilenamePrefix(gOptEvFilePrefix);
  ntpw.Initialize();

  // Create an MC Job Monitor
  GMCJMonitor mcjmonitor(gOptRunNu);
  mcjmonitor.SetRefreshRate(RunOpt::Instance()->MCJobStatusRefreshRate());

  LOG("gevgen_capt", pNOTICE) 
    << "\n ** Will generate " << gOptNevents << " events for \n" 
    << init_state << " at Ev = " << Ev << " GeV";

  // Generate events / print the GHEP record / add it to the ntuple
  int ievent = 0;
  while (ievent < gOptNevents) {
     LOG("gevgen_capt", pNOTICE) 
        << " *** Generating event............ " << ievent;

     // generate a single event
     EventRecord * event = evg_driver.GenerateEvent(nu_p4);

     if(!event) {
        LOG("gevgen_capt", pNOTICE) 
          << "Last attempt failed. Re-trying....";
        continue;
     }

     LOG("gevgen_capt", pNOTICE) 
        << "Generated Event GHEP Record: " << *event;

     // add event at the output ntuple, refresh the mc job monitor & clean up
     ntpw.AddEventRecord(ievent, event);
     mcjmonitor.Update(ievent,event);
     ievent++;
     delete event;
  }

  // Save the generated MC events
  ntpw.Save();
}
//____________________________________________________________________________

//............................................................................
void GenerateEventsUsingFluxOrTgtMix(void)
{
  // Get flux and geom drivers
  GFluxI *        flux_driver = FluxDriver();
  GeomAnalyzerI * geom_driver = GeomDriver();

  // Create the monte carlo job driver
  GMCJDriver * mcj_driver = new GMCJDriver;
  mcj_driver->SetEventGeneratorList(RunOpt::Instance()->EventGeneratorList());
  mcj_driver->SetUnphysEventMask(*RunOpt::Instance()->UnphysEventMask());
  mcj_driver->UseFluxDriver(flux_driver);
  mcj_driver->UseGeomAnalyzer(geom_driver);
  mcj_driver->Configure();
  mcj_driver->UseSplines();
  if(!gOptWeighted) 
        mcj_driver->ForceSingleProbScale();

  // Initialize an Ntuple Writer to save GHEP records into a TTree
  NtpWriter ntpw(kDefOptNtpFormat, gOptRunNu);
  ntpw.CustomizeFilenamePrefix(gOptEvFilePrefix);
  ntpw.Initialize();

  // Create an MC Job Monitor
  GMCJMonitor mcjmonitor(gOptRunNu);
  mcjmonitor.SetRefreshRate(RunOpt::Instance()->MCJobStatusRefreshRate());

  // Generate events / print the GHEP record / add it to the ntuple
  int ievent = 0;
  while ( ievent < gOptNevents) {

     LOG("gevgen_capt", pNOTICE) << " *** Generating event..... " << ievent;

     // generate a single event for neutrinos coming from the specified flux
     EventRecord * event = mcj_driver->GenerateEvent();

     LOG("gevgen_capt", pNOTICE) << "Generated Event GHEP Record: " << *event;

     // add event at the output ntuple, refresh the mc job monitor & clean-up
     ntpw.AddEventRecord(ievent, event);
     mcjmonitor.Update(ievent,event);
     ievent++;
     delete event;
  }

  // Save the generated MC events
  ntpw.Save();

  delete flux_driver;
  delete geom_driver;
  delete mcj_driver;;
}

//____________________________________________________________________________
GeomAnalyzerI * GeomDriver(void)
{
// create a trivial point geometry with the specified target or target mix

  GeomAnalyzerI * geom_driver = new geometry::PointGeomAnalyzer(gOptTgtMix);
  return geom_driver;
}

//____________________________________________________________________________
GFluxI * FluxDriver(void)
{
// create & configure one of the generic flux drivers
//
  GFluxI * flux_driver = 0;

  if(gOptNuEnergyRange<0) flux_driver = MonoEnergeticFluxDriver();
  else flux_driver = TH1FluxDriver();

  return flux_driver;
}

//____________________________________________________________________________
GFluxI * MonoEnergeticFluxDriver(void)
{
//
//
  flux::GMonoEnergeticFlux * flux = 
              new flux::GMonoEnergeticFlux(gOptNuEnergy, gOptNuPdgCode);
  GFluxI * flux_driver = dynamic_cast<GFluxI *>(flux);
  return flux_driver;
}

//____________________________________________________________________________
GFluxI * TH1FluxDriver(void)
{
// 
//
  flux::GCylindTH1Flux * flux = new flux::GCylindTH1Flux;
  TH1D * spectrum = 0;

  int flux_entries = 100000;

  double emin = gOptNuEnergy;
  double emax = gOptNuEnergy+gOptNuEnergyRange;
  double de   = gOptNuEnergyRange;

  // check whether the input flux is a file or a functional form
  //
  bool input_is_text_file = ! gSystem->AccessPathName(gOptFlux.c_str());
  bool input_is_root_file = gOptFlux.find(".root") != string::npos &&
                            gOptFlux.find(",") != string::npos;
  if (input_is_text_file) {
    //
    // ** generate the flux histogram from the x,y pairs in the input text file
    //
    Spline * input_flux = new Spline(gOptFlux.c_str());
    int  n = 100;
    double estep = (emax-emin)/(n-1);
    double ymax  = -1, ry = -1, gy = -1, e = -1;
    for(int i=0; i<n; i++) {
      e = emin + i*estep;
      ymax = TMath::Max(ymax, input_flux->Evaluate(e));
    }
    ymax *= 1.3;

    RandomGen * r = RandomGen::Instance();
    spectrum  = new TH1D("spectrum","neutrino flux", 300, emin, emax);
    spectrum->SetDirectory(0);

    for(int ientry=0; ientry<flux_entries; ientry++) {
      bool accept = false;
      unsigned int iter=0;
      while(!accept) {
        iter++;
        if(iter > kRjMaxIterations) {
           LOG("gevgen_capt", pFATAL) << "Couldn't generate a flux histogram";
           exit(1);
        }
        e = emin + de * r->RndGen().Rndm();
        gy = ymax * r->RndGen().Rndm();
        ry = input_flux->Evaluate(e);
        accept = gy < ry;
        if(accept) spectrum->Fill(e);
      }
    }
    delete input_flux;
  } 
  else if (input_is_root_file) {
    //
    // ** extract specified flux histogram from the input root file
    //
    vector<string> fv = utils::str::Split(gOptFlux,",");
    assert(fv.size()==2); 
    assert( !gSystem->AccessPathName(fv[0].c_str()) );

    LOG("gevgen_capt", pNOTICE) << "Getting input flux from root file: " << fv[0];
    TFile * flux_file = new TFile(fv[0].c_str(),"read");

    LOG("gevgen_capt", pNOTICE) << "Flux name: " << fv[1];
    TH1D * hst = (TH1D *)flux_file->Get(fv[1].c_str());
    assert(hst);

    LOG("gevgen_capt", pNOTICE) << hst->GetEntries();

    // copy all bins between emin,emax
    spectrum  = new TH1D("spectrum","neutrino flux",
                         hst->GetNbinsX(),
                         hst->GetXaxis()->GetXmin(),
                         hst->GetXaxis()->GetXmax());
    spectrum->SetDirectory(0);
    for(int ibin = 1; ibin <= hst->GetNbinsX(); ibin++) {
       if(hst->GetBinLowEdge(ibin) < emax &&
         hst->GetBinLowEdge(ibin) + hst->GetBinWidth(ibin) > emin) {
           spectrum->SetBinContent(ibin, hst->GetBinContent(ibin));
            LOG("gevgen_capt", pNOTICE) 
              << "adding => " << ibin << ": " << hst->GetBinContent(ibin);
         }
    }

    LOG("gevgen_capt", pNOTICE) << spectrum->GetEntries();

    flux_file->Close();        
    delete flux_file;

    LOG("gevgen_capt", pNOTICE) << spectrum->GetEntries();

  } 
  else {
    //
    // ** generate the flux histogram from the input functional form
    //
    TF1 *  input_func = new TF1("input_func", gOptFlux.c_str(), emin, emax);
    spectrum  = new TH1D("spectrum","neutrino flux", 300, emin, emax);
    spectrum->SetDirectory(0);
    spectrum->FillRandom("input_func", flux_entries);
    delete input_func;
  }

  // save input flux into a standard location.
  TFile f("./input-flux.root","recreate");
  spectrum->Write();
  f.Close();

  TVector3 bspot(0,0,0);

  flux->SetNuDirection      (gOptNuDir);
  flux->SetBeamSpot         (bspot);
  flux->SetTransverseRadius (-1);
  flux->AddEnergySpectrum   (gOptNuPdgCode, spectrum);

  GFluxI * flux_driver = dynamic_cast<GFluxI *>(flux);
  return flux_driver;
}
//............................................................................

//____________________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
  LOG("gevgen_capt", pINFO) << "Parsing command line arguments";

  // Common run options. Set defaults and read.
  RunOpt::Instance()->EnableBareXSecPreCalc(true);
  RunOpt::Instance()->ReadFromCommandLine(argc,argv);

  // Parse run options for this app

  CmdLnArgParser parser(argc,argv);

  // help?
  bool help = parser.OptionExists('h');
  if(help) {
      PrintSyntax();
      exit(0);
  }

  // number of events
  if( parser.OptionExists('n') ) {
    LOG("gevgen_capt", pINFO) << "Reading number of events to generate";
    gOptNevents = parser.ArgAsInt('n');
  } else {
    LOG("gevgen_capt", pINFO)
       << "Unspecified number of events to generate - Using default";
    gOptNevents = kDefOptNevents;
  }

  // run number
  if( parser.OptionExists('r') ) {
    LOG("gevgen_capt", pINFO) << "Reading MC run number";
    gOptRunNu = parser.ArgAsLong('r');
  } else {
    LOG("gevgen_capt", pINFO) << "Unspecified run number - Using default";
    gOptRunNu = kDefOptRunNu;
  }

  // event file prefix
  if( parser.OptionExists('o') ) {
    LOG("gevgen_capt", pINFO) << "Setting the event filename prefix";
    gOptEvFilePrefix = parser.ArgAsString('o');
  } else {
    LOG("gevgen_capt", pINFO)
      << "Will set the default event filename prefix";
    gOptEvFilePrefix = kDefOptEvFilePrefix;
  } //-o

  // flux functional form
  bool using_flux = false;
  if( parser.OptionExists('f') ) {
    LOG("gevgen_capt", pINFO) << "Reading flux function";
    gOptFlux = parser.ArgAsString('f');
    using_flux = true;
  }

  if(parser.OptionExists('s')) {
    LOG("gevgen_capt", pWARN) 
      << "-s option no longer available. Please read the revised code documentation";
    gAbortingInErr = true;
    exit(1);
  }


  // generate weighted events option (only relevant if using a flux)
  gOptWeighted = parser.OptionExists('w');

  // neutrino energy
  if( parser.OptionExists('e') ) {
    LOG("gevgen_capt", pINFO) << "Reading neutrino energy";
    string nue = parser.ArgAsString('e');

    // is it just a value or a range (comma separated set of values)
    if(nue.find(",") != string::npos) {
       // split the comma separated list
       vector<string> nurange = utils::str::Split(nue, ",");
       assert(nurange.size() == 2);   
       double emin = atof(nurange[0].c_str());
       double emax = atof(nurange[1].c_str());
       assert(emax>emin && emin>=0);
       gOptNuEnergy      = emin;
       gOptNuEnergyRange = emax-emin;
       if(!using_flux) {
          LOG("gevgen_capt", pWARN) 
             << "No flux was specified but an energy range was input!";
          LOG("gevgen_capt", pWARN) 
             << "Events will be generated at fixed E = " 
             << gOptNuEnergy << " GeV";
          gOptNuEnergyRange = -1;
       }
    } else {
       gOptNuEnergy       = atof(nue.c_str());
       gOptNuEnergyRange = -1;
    }
  } else {
    LOG("gevgen_capt", pFATAL) << "Unspecified neutrino energy - Exiting";
    PrintSyntax();
    exit(1);
  }

  // neutrino PDG code
  if( parser.OptionExists('p') ) {
    LOG("gevgen_capt", pINFO) << "Reading neutrino PDG code";
    gOptNuPdgCode = parser.ArgAsInt('p');
  } else {
    LOG("gevgen_capt", pFATAL) << "Unspecified neutrino PDG code - Exiting";
    PrintSyntax();
    exit(1);
  }

  // target mix (their PDG codes with their corresponding weights)
  bool using_tgtmix = false;
  if( parser.OptionExists('t') ) {
    LOG("gevgen_capt", pINFO) << "Reading target mix";
    string stgtmix = parser.ArgAsString('t');
    gOptTgtMix.clear();
    vector<string> tgtmix = utils::str::Split(stgtmix,",");
    if (tgtmix.size()==1) {
         int    pdg = atoi(tgtmix[0].c_str());
         double wgt = 1.0;
         gOptTgtMix.insert(map<int, double>::value_type(pdg, wgt));
    } 
    else {
      using_tgtmix = true;
      vector<string>::const_iterator tgtmix_iter = tgtmix.begin();
      for( ; tgtmix_iter != tgtmix.end(); ++tgtmix_iter) {
        string tgt_with_wgt = *tgtmix_iter;
        string::size_type open_bracket  = tgt_with_wgt.find("[");
        string::size_type close_bracket = tgt_with_wgt.find("]");
        string::size_type ibeg = 0;
        string::size_type iend = open_bracket;
        string::size_type jbeg = open_bracket+1;
        string::size_type jend = close_bracket-1;
        int    pdg = atoi(tgt_with_wgt.substr(ibeg,iend).c_str());
        double wgt = atof(tgt_with_wgt.substr(jbeg,jend).c_str());
        LOG("Main", pNOTICE)
          << "Adding to target mix: pdg = " << pdg << ", wgt = " << wgt;
        gOptTgtMix.insert(map<int, double>::value_type(pdg, wgt));
      }//tgtmix_iter
    }
  } else {
    LOG("gevgen_capt", pINFO) 
      << "Unspecified target PDG code - Assuming natural Argon";
    gOptTgtMix[1000180400] = 0.9960;
    gOptTgtMix[1000180360] = 0.0034;
    gOptTgtMix[1000180380] = 0.0006;
  }

  gOptUsingFluxOrTgtMix = using_flux || using_tgtmix;

  // Set the neutrino direction.
  if (parser.OptionExists('d')) {
    LOG("gevgen_capt", pINFO) << "Reading neutrino direction";
    std::vector<double> dcos = parser.ArgAsDoubleTokens('d',",");
    if (dcos.size() != 3) {
      LOG("gevgen_capt", pFATAL) << "Wrong number of direction cosines";
      exit(1);
    }
    gOptNuDir = TVector3(dcos[0], dcos[1], dcos[2]).Unit();
  }

  // random number seed
  if( parser.OptionExists("seed") ) {
    gOptRanSeed = parser.ArgAsLong("seed");
    LOG("gevgen_capt", pINFO) 
      << "Setting random number seed to " 
      << gOptRanSeed;
  } else {
    gOptRanSeed = 0;
    LOG("gevgen_capt", pINFO) 
      << "Unspecified random number seed - Using time based seed.";
  }

  // input cross-section file
  if( parser.OptionExists("cross-sections") ) {
    gOptInpXSecFile = parser.ArgAsString("cross-sections");
    LOG("gevgen_capt", pINFO) << "Reading cross-section file"
                              << gOptInpXSecFile;
  } 
  else {
    const char* genie = gSystem->Getenv("GENIE");
    if (!genie) {
      LOG("gevgen_capt", pFATAL) << "GENIE source location not set.";
      std::exit(1);
    }
    gOptInpXSecFile = std::string(genie) 
      + "/data/cross-section/"
      + "gxspl-t2k-v" + __GENIE_RELEASE__ + ".xml";
    LOG("gevgen_capt", pINFO) << "Reading cross-section file"
                              << gOptInpXSecFile;
  }

  //
  // print-out the command line options
  //
  LOG("gevgen_capt", pNOTICE) 
     << "\n" 
     << utils::print::PrintFramedMesg("gevgen_capt job configuration");
  LOG("gevgen_capt", pNOTICE) 
     << "MC Run Number: " << gOptRunNu;
  if (gOptRanSeed != -1) {
     if (gOptRanSeed == 0) {
        gRandom->SetSeed(0);
        gOptRanSeed = gRandom->GetSeed();
     }
     LOG("gevgen_capt", pNOTICE) 
       << "Random number seed: " << gOptRanSeed;
  } 
  else {
     LOG("gevgen_capt", pNOTICE) 
       << "Random number set to default value";
  }
  LOG("gevgen_capt", pNOTICE) 
       << "Number of events requested: " << gOptNevents;
  if(gOptInpXSecFile.size() > 0) {
     LOG("gevgen_capt", pNOTICE) 
       << "Using cross-section splines read from: " << gOptInpXSecFile;
  } else {
     LOG("gevgen_capt", pNOTICE) 
       << "No input cross-section spline file";
  }
  LOG("gevgen_capt", pNOTICE) 
       << "Flux: " << gOptFlux;
  LOG("gevgen_capt", pNOTICE) 
       << "Generate weighted events? " << gOptWeighted;
  if(gOptNuEnergyRange>0) {
     LOG("gevgen_capt", pNOTICE) 
        << "Neutrino energy: [" 
        << gOptNuEnergy << ", " << gOptNuEnergy+gOptNuEnergyRange << "]";
  } else {
     LOG("gevgen_capt", pNOTICE) 
        << "Neutrino energy: " << gOptNuEnergy;
  }
  LOG("gevgen_capt", pNOTICE) 
      << "Neutrino code (PDG): " << gOptNuPdgCode;
  LOG("gevgen_capt", pNOTICE) 
      << "Target code (PDG) & weight fraction (in case of multiple targets): ";
  map<int,double>::const_iterator iter;
  for(iter = gOptTgtMix.begin(); iter != gOptTgtMix.end(); ++iter) {
      int    tgtpdgc = iter->first;
      double wgt     = iter->second;
      LOG("gevgen_capt", pNOTICE) 
          << " >> " <<  tgtpdgc << " (weight fraction = " << wgt << ")";
  }
  LOG("gevgen_capt", pNOTICE) << "\n";

  LOG("gevgen_capt", pNOTICE) << *RunOpt::Instance();

}
//____________________________________________________________________________
void PrintSyntax(void)
{
  LOG("gevgen_capt", pNOTICE)
    << "\n\n" << "Syntax:" << "\n"
    << "\n      gevgen_capt [-h]"
    << "\n              [-r run#]"
    << "\n               -n nev"
    << "\n               -e energy or lowEnergy,highEnergy "
    << "\n               -p neutrino_pdg" 
    << "\n              [-t target_pdg] "
    << "\n              [-d dCosX,dCosY,dCosZ] (default: 1,0,0)"
    << "\n              [-o file_prefix]"
    << "\n              [-f flux_description]"
    << "\n              [-w]"
    << "\n              [--seed random_number_seed]"
    << "\n              [--cross-sections xml_file]"
    << "\n              [--event-generator-list list_name]"
    << "\n              [--message-thresholds xml_file]"
    << "\n              [--unphysical-event-mask mask]"
    << "\n              [--event-record-print-level level]"
    << "\n              [--mc-job-status-refresh-rate  rate]"
    << "\n              [--cache-file root_file]"
    << "\n";
}
//____________________________________________________________________________

// Local Variables:
// mode:c++
// c-basic-offset:2
// End:
