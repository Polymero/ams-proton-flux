// Written by Sebastiaan Venendaal (University of Groningen, the Netherlands)
// C++ class for generating histograms of proton-like AMS-02 data, used for flux analysis
// Created          16-05-23
// Last modified    16-05-23
//
// Usage ::
// ...

//-----------------------------------------------------------------------------------
// HEADER FILES
//-----------------------------------------------------------------------------------

// Native C headers
#include <algorithm>
#include <iostream>
#include <string>
// Native ROOT headers
#include "TChain.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2.h"
#include "TCanvas.h"
#include "Tobject.h"
#include "TString.h"
// Local headers
#include "Header Files/Ntp.h"


//-----------------------------------------------------------------------------------
// CLASS DEFINITION
//-----------------------------------------------------------------------------------

class MIRJA {
    // Access Specifier
    public:


    //-------------------------------------------------------------------------------
    // CLASS ATTRIBUTES
    //-------------------------------------------------------------------------------

    // Rigidity bins (based on equal logarithmic widths)
    const int binNumber = 32;
    double binEdges[32 + 1] = {
        1.00, 1.16, 1.33, 1.51, 1.71, 1.92, 2.15, 2.40, 2.67, 2.97, 3.29, 3.64, 4.02,
        4.43, 4.88, 5.37, 5.90, 6.47, 7.09, 7.76, 8.48, 9.26, 10.1, 11.0, 12.0, 13.0,
        14.1, 15.3, 16.6, 18.0, 19.5, 21.1, 22.8
    };
    double binErrors[32];
    double binCentres[32];

    // Rigidity cut-off level (based on ?)
    double rigidityCutOff = 1.2;

    // New file object
    TFile *f = new TFile("/afs/cern.ch/usr/s/svenenda/public/ams-proton-flux/ProtonHistogramsAMS02.root", "recreate");

    // RTI map
    map<int, std::pair<float, float>> RTIMap = map<int, std::pair<float, float>>();

    // List of histograms
    // Proton compact data
    TH1F *exposureTime          = new TH1F("exposureTime", "Exposure Time per Rigidity Bin", 32, binEdges);
    TH1F *eventsDetected        = new TH1F("eventsDetected", "Detected Events per Rigidity Bin", 32, binEdges);
    TH1F *eventsSelected        = new TH1F("eventsSelected", "Selected Proton Events per Rigidity Bin", 32, binEdges);
    TH1F *triggersPhysical      = new TH1F("triggersPhysical", "Proton Physical Triggers per Rigidity Bin", 32, binEdges);
    TH1F *triggersBias          = new TH1F("triggersBias", "Proton Bias Triggers per Rigidity Bin", 32, binEdges);
    TH1F *baseTracker           = new TH1F("baseTracker", "Proton Tracker Base", 32, binEdges);
    TH1F *baseTOF               = new TH1F("baseTOF", "Proton TOF Base", 32, binEdges);
    TH1F *cutParticle           = new TH1F("cutParticle", "Proton Particle Cut", 32, binEdges);
    TH1F *cutBeta               = new TH1F("cutBeta", "Proton Beta Cut", 32, binEdges);
    TH1F *cutChiSquared         = new TH1F("cutChiSquared", "Proton Chi Squared Cut", 32, binEdges);
    TH1F *cutInnerLayer         = new TH1F("cutInnerLayer", "Proton Inner Layer Cut", 32, binEdges);
    // Monte-Carlo proton compact data
    TH1F *montecarloDetected    = new TH1F("montecarloDetected", "MC Proton Detected Events per Rigidity Bin", 32, binEdges);
    TH1F *montecarloSelected    = new TH1F("montecarloSelected", "MC Proton Selected Events per Rigidity Bin", 32, binEdges);
    TH1F *montecarloPhysical    = new TH1F("montecarloPhysical", "MC Proton Physical Triggers per Rigidity Bin", 32, binEdges);
    TH1F *montecarloBias        = new TH1F("montecarloBias", "MC Proton Bias Triggers per Rigidity Bin", 32, binEdges);
    TH1F *montecarloTracker     = new TH1F("montecarloTracker", "MC Proton Tracker Base", 32, binEdges);
    TH1F *montecarloTOF         = new TH1F("montecarloTOF", "MC Proton TOF Base", 32, binEdges);
    TH1F *montecarloParticle    = new TH1F("montecarloParticle", "MC Proton Particle Cut", 32, binEdges);
    TH1F *montecarloBeta        = new TH1F("montecarloBeta", "MC Proton Beta Cut", 32, binEdges);
    TH1F *montecarloChiSquared  = new TH1F("montecarloChiSquared", "MC Proton Chi Squared Cut", 32, binEdges);
    TH1F *montecarloInnerLayer  = new TH1F("montecarloInnerLayer", "MC Proton Inner Layer Cut", 32, binEdges);
    // Monte-Carlo proton FileMCInfo
    TH1F *montecarloGenerated   = new TH1F("montecarloGenerated", "MC Proton Generated Events per Rigidity Bin", 32, binEdges);

    // List of data objects
    // Chains
    TChain *chainCompact        = new TChain("Compact");
    TChain *chainRTI            = new TChain("RTI");
    TChain *chainMCCompact      = new TChain("Compact");
    TChain *chainMCInfo         = new TChain("File");
    // Classes
    NtpCompact *classCompact    = new NtpCompact();
    NtpSHeader *classSHeader    = new class NtpSHeader();
    RTIInfo *classRTI           = new class RTIInfo();
    NtpCompact *classMCCompact  = new NtpCompact();
    FileMCInfo *classMCInfo     = new FileMCInfo();


    //-------------------------------------------------------------------------------
    // CLASS CONSTRUCTORS
    //-------------------------------------------------------------------------------

    MIRJA(string rootFiles = "13300*.root") { // Default constructor

        // ROOT gStyle configuration
        gStyle->SetOptTitle(0);
        gStyle->SetOptStat(0);
        gStyle->SetOptLogx(1);

        // Bin properties
        for (int i=0; i < binNumber; i++) {
            binErrors[i]  = (binEdges[i + 1] - binEdges[i]) / 2;
            binCentres[i] = (binEdges[i + 1] + binEdges[i]) / 2;
        }

        // Read the data trees
        chainCompact->Add("/eos/ams/group/dbar/release_v7/e1_vdev_200421/neg/ISS.B1130/pass7/" + rootFiles);
        chainRTI->Add("/eos/ams/group/dbar/release_v7/e1_vdev_200421/neg/ISS.B1130/pass7/" + rootFiles);
        chainMCCompact->Add("/eos/ams/group/dbar/release_v7/e1_vdev_200421/full/Pr.B1200/pr.pl1.05100.4_00/*.root");
        chainMCInfo->Add("/eos/ams/group/dbar/release_v7/e1_vdev_200421/full/Pr.B1200/pr.pl1.05100.4_00/*.root");

        // Set branch addresses
        chainCompact->SetBranchAddress("Compact", &classCompact);
        chainCompact->SetBranchAddress("SHeader", &classSHeader);
        chainRTI->SetBranchAddress("RTIInfo", &classRTI);
        chainMCCompact->SetBranchAddress("Compact", &classMCCompact);
        chainMCInfo->SetBranchAddress("FileMCInfo", &classMCInfo);

        cout << "\nClass succesfully constructed!\n" << endl;

    };


    //-------------------------------------------------------------------------------
    // CLASS METHODS
    //-------------------------------------------------------------------------------

    void runAnalysis();

};


//-----------------------------------------------------------------------------------
// SUPPORT FUNCTIONS
//-----------------------------------------------------------------------------------

// No support functions


//-----------------------------------------------------------------------------------
// METHOD FUNCTIONS
//-----------------------------------------------------------------------------------

void MIRJA::runAnalysis() {

    cout << "Starting MIRJA.runAnalysis()..." << endl;


    //-------------------------------------------------------------------------------
    // (1/6)
    //-------------------------------------------------------------------------------
    cout << "Looping over RTIInfo data... (1/6)" << endl;

    int chainRTINumber = chainRTI->GetEntries();
    cout << "Number of RTIInfo entries: " << chainRTINumber << endl;

    // Looping over RTI files
    for (int i=0; i < chainRTINumber; i++) {

        // Get entry
        chainRTI->GetEntry(i);

        // Fill RTI map
        RTIMap.insert({classRTI->utime, std::pair<float, float>(classRTI->lf, classRTI->cf[0][3][1])});

        // Loop over rigidity bins
        for (int j=0; j < binNumber; j++) {

            // Exposure() --> Get total livetime as a function of rigidity
            // If bin centre is above geo-matgnetic cut-off, include the livetime
            if (binCentres[j] > rigidityCutOff * classRTI->cf[0][3][1]) {
                exposureTime->SetBinContent(j + 1, exposureTime->GetBinContent(j + 1) + classRTI->lf);
            }

        }

        // Progress tracker
        int progress = (chainRTINumber / 100);
        if (i % progress == 0) {
            cout << "#" << flush;
        }

    }


    //-------------------------------------------------------------------------------
    // (2/6)
    //-------------------------------------------------------------------------------
    cout << "\nLooping over Compact data... (2/6)" << endl;

    int chainCompactNumber = chainCompact->GetEntries();
    cout << "Number of Compact entries: " << chainCompactNumber << endl;

    // Loop over Compact data entries
    for (int i=0; i < chainCompactNumber; i++){

        // Get entry
        chainCompact->GetEntry(i);

        // List of boolean cuts
        // Geomagnetic cut-off
        bool boolCutOff     = classCompact->trk_rig[0] > rigidityCutOff * RTIMap[classSHeader->utime].second;
        // Within our rigidity range
        bool boolRigidity   = (classCompact->trk_rig[0] > binEdges[0]) && (classCompact->trk_rig[0] <= binEdges[binNumber]);
        // Correct trigger pattern
        bool boolTriggers   = ((classCompact->sublvl1 & 0x3E) != 0) && ((classCompact->trigpatt & 0x02) != 0);
        // Particle-like events
        bool boolParticle   = classCompact->status % 10 == 1;
        // TOF Beta selection
        bool boolBeta       = classCompact->tof_beta > 0.3;
        // Chi-Squared selection
        bool boolChiSquared = (classCompact->trk_chisqn[0][0] < 10) && (classCompact->trk_chisqn[0][1] < 10) && (classCompact->trk_chisqn[0][0] > 0) && (classCompact->trk_chisqn[0][1] > 0);
        // Inner Layer selection
        bool boolInnerLayer = (classCompact->trk_q_inn > 0.80) && (classCompact->trk_q_inn < 1.30);
        
        int boolBit = boolCutOff + (boolRigidity << 1) + (boolTriggers << 2) + (boolParticle << 3) + 
                      (boolBeta << 4) + (boolChiSquared << 5) + (boolInnerLayer << 6);

        // RigBinner() --> Bin events as a function of rigidity
        eventsDetected->Fill(classCompact->trk_rig[0]);
        if ((boolBit & 0x7F) == 0x7F) { // 0x7F = 0b01111111 (All)
            eventsSelected->Fill(classCompact->trk_rig[0]);
        }

        // TrigEff(): Data --> Trigger efficiency as a function fo rigidity
        // List of trigger booleans
        bool boolPhysical   = ((classCompact->sublvl1 & 0x3E) != 0) && ((classCompact->trigpatt & 0x02) != 0);
        bool boolUnphysical = ((classCompact->sublvl1 & 0x3E) == 0) && ((classCompact->trigpatt & 0x02) != 0);

        // Trigger histograms
        if ((boolBit & 0x7B) == 0x7B) { // 0x7B = 0b01111011 (All but Triggers)
            if (boolPhysical) {
                triggersPhysical->Fill(classCompact->trk_rig[0]);
            }
            if (boolUnphysical) {
                triggersBias->Fill(classCompact->trk_rig[0]);
            }
        }

        // SelEff(): Data --> Selection efficiency of applied cuts as a function of rigidity
        // Additional TOF charge cuts (to replace TRK charge cuts)
        bool boolTOFCharge = (classCompact->tof_q_lay[0] > 0.8) && (classCompact->tof_q_lay[0] < 1.5);

        // TRK base histogram
        if ((boolBit & 0x17) == 0x17) { // 0x17 = 0b00010111 (Beta, Triggers, Rigidity, CutOff)
            if (boolTOFCharge) {
                baseTracker->Fill(classCompact->trk_rig[0]);
            }
        }

        // TOF base histogram
        if ((boolBit & 0x6F) == 0x6F) { // 0x6F = 0b01101111 (All but Beta)
            baseTOF->Fill(classCompact->trk_rig[0]);
        }

        // Particle-like selection (TRK base)
        if ((boolBit & 0x3E) == 0x3E) { // 0x3E = 0b00011111 (All but InnerLayer, ChiSquared)
            if (boolTOFCharge) {
                cutParticle->Fill(classCompact->trk_rig[0]);
            }
        }

        // Beta selection (TOF base)
        if ((boolBit & 0x7F) == 0x7F) { // 0x7F = 0b01111111 (All)
            cutBeta->Fill(classCompact->trk_rig[0]);
        }

        // Chi-Squared selection (TRK base)
        if ((boolBit & 0x37) == 0x37) { // 0x37 = 0b00110111 (All but Innerlayer, Particle)
            if (boolTOFCharge) {
                cutChiSquared->Fill(classCompact->trk_rig[0]);
            }
        }

        // Inner Layer selection (TRK base w/o TOFCharge cut)
        if ((boolBit & 0x57) == 0x57) { // 0x57 = 0b01010111 (All but Particle, ChiSquared)
            cutInnerLayer->Fill(classCompact->trk_rig[0]);
        }

        // Progress tracker
        int progress = (chainCompactNumber / 100);
        if (i % progress == 0) {
            cout << "#" << flush;
        }

    }


    //-------------------------------------------------------------------------------
    // (3/6)
    //-------------------------------------------------------------------------------
    cout << "\nLooping over Proton Monte-Carlo Compact data... (3/6)" << endl;

    int chainMCCompactNumber = chainMCCompact->GetEntries();
    cout << "Number of Proton Monte-Carlo Compact entries: " << chainMCCompactNumber << endl;

    // Loop over Proton MC Compact entries
    for (int i=0; i < chainMCCompactNumber; i++) {

        // Get entry
        chainMCCompact->GetEntry(i);

        // List of boolean cuts
        // Within our rigidity range
        bool boolRigidity   = (classMCCompact->trk_rig[0] > binEdges[0]) && (classMCCompact->trk_rig[0] <= binEdges[binNumber]);
        // Correct trigger pattern
        bool boolTriggers   = ((classMCCompact->sublvl1 & 0x3E) != 0) && ((classMCCompact->trigpatt & 0x02) != 0);
        // Particle-like events
        bool boolParticle   = classMCCompact->status % 10 == 1;
        // TOF Beta selection
        bool boolBeta       = classMCCompact->tof_beta > 0.3;
        // Chi-Squared selection
        bool boolChiSquared = (classMCCompact->trk_chisqn[0][0] < 10) && (classMCCompact->trk_chisqn[0][1] < 10) && (classMCCompact->trk_chisqn[0][0] > 0) && (classMCCompact->trk_chisqn[0][1] > 0);
        // Inner Layer selection
        bool boolInnerLayer = (classMCCompact->trk_q_inn > 0.80) && (classMCCompact->trk_q_inn < 1.30);
        
        int boolBit = 1 + (boolRigidity << 1) + (boolTriggers << 2) + (boolParticle << 3) + 
                      (boolBeta << 4) + (boolChiSquared << 5) + (boolInnerLayer << 6);
        
        // Acceptance() --> Geometric Aperature Acceptance as a function of rigidity
        montecarloDetected->Fill(chainMCCompact->trk_rig[0]);
        if ((boolBit & 0x7F) == 0x7F) { // 0x7F = 0b01111111 (All)
            montecarloSelected->Fill(chainMCCompact->trk_rig[0]);
        }

        // TrigEff(): MC --> Trigger efficiency as a function of rigidity
        // List of trigger booleans
        bool boolPhysical   = ((classMCCompact->sublvl1 & 0x3E) != 0) && ((classMCCompact->trigpatt & 0x02) != 0);
        bool boolUnphysical = ((classMCCompact->sublvl1 & 0x3E) == 0) && ((classMCCompact->trigpatt & 0x02) != 0);

        // Trigger histograms
        if ((boolBit & 0x7B) == 0x7B) { // 0x7B = 0b01111011 (All but Triggers)
            if (boolPhysical) {
                montecarloPhysical->Fill(classMCCompact->trk_rig[0]);
            }
            if (boolUnphysical) {
                montecarloBias->Fill(classMCCompact->trk_rig[0]);
            }
        }

        // SelEff(): MC --> Selection efficiency of applied cuts as a function of rigidity
        // Additional TOF charge cuts (to replace TRK charge cuts)
        bool boolTOFCharge = (classMCCompact->tof_q_lay[0] > 0.8) && (classMCCompact->tof_q_lay[0] < 1.5);

        // TRK base histogram
        if ((boolBit & 0x17) == 0x17) { // 0x17 = 0b00010111 (Beta, Triggers, Rigidity)
            if (boolTOFCharge) {
                montecarloTracker->Fill(classMCCompact->trk_rig[0]);
            }
        }

        // TOF base histogram
        if ((boolBit & 0x6F) == 0x6F) { // 0x6F = 0b01101111 (All but Beta)
            montecarloTOF->Fill(classMCCompact->trk_rig[0]);
        }

        // Particle-like selection (TRK base)
        if ((boolBit & 0x3E) == 0x3E) { // 0x3E = 0b00011111 (All but InnerLayer, ChiSquared)
            if (boolTOFCharge) {
                montecarloParticle->Fill(classMCCompact->trk_rig[0]);
            }
        }

        // Beta selection (TOF base)
        if ((boolBit & 0x7F) == 0x7F) { // 0x7F = 0b01111111 (All)
            montecarloBeta->Fill(classMCCompact->trk_rig[0]);
        }

        // Chi-Squared selection (TRK base)
        if ((boolBit & 0x37) == 0x37) { // 0x37 = 0b00110111 (All but Innerlayer, Particle)
            if (boolTOFCharge) {
                montecarloChiSquared->Fill(classMCCompact->trk_rig[0]);
            }
        }

        // Inner Layer selection (TRK base w/o TOFCharge cut)
        if ((boolBit & 0x57) == 0x57) { // 0x57 = 0b01010111 (All but Particle, ChiSquared)
            montecarloInnerLayer->Fill(classMCCompact->trk_rig[0]);
        }

        // Progress tracker
        int progress = (chainCompactNumber / 100);
        if (i % progress == 0) {
            cout << "#" << flush;
        }

    }


    //-------------------------------------------------------------------------------
    // (4/6)
    //-------------------------------------------------------------------------------
    cout << "\nLooping over Proton Monte-Carlo File Info data... (4/6)" << endl;

    int chainMCInfoNumber = chainMCInfo->GetEntries();
    cout << "Number of Proton Monte-Carlo FileMCInfo entries: " << chainMCInfoNumber << endl;

    // Loop over MC Proton FileMCInfo entries
    for (int i=0; i < chainMCInfoNumber; i++) {

        // Get entry
        chainMCInfo->GetEntry(i);

        // Get MC generation parameters
        double generatedNumber = chainMCInfo->ngen_datacard;
        double rigidityMinimum = chainMCInfo->momentum[0];
        double rigidityMaximum = chainMCInfo->momentum[1];

        // 1/R generation spectrum (normalised)
        TF1 *generatedFlux = new TF1("generatedFlux", "[0]/(x)", rigidityMinimum, rigidityMaximum);
        generatedFlux->SetParameter(0, log(rigidityMaximum / rigidityMinimum));

        // Loop over rigidity bins
        for (int j=0; j < binNumber; j++) {

            // Fraction of spectrum within the bin
            double fraction = generatedFlux->Integral(binEdges[j], binEdges[j + 1]) / generatedFlux->Integral(rigidityMinimum, rigidityMaximum);

            // Number of events within said fraction
            double eventsNumber = fraction * generatedNumber;

            // Set bin value to current bin value + number of events
            montecarloGenerated->SetBinContent(j + 1, montecarloGenerated->GetBinContent(j + 1) + eventsNumber);

        }

        // Progress tracker
        int progress = (chainMCInfoNumber / 100);
        if (i % progress = 0) {
            cout << "#" << flush;
        }

    }


    //-------------------------------------------------------------------------------
    // (45/6)
    //-------------------------------------------------------------------------------
    cout << "\nSaving all my hard work..." << endl;

    f->Write();
    f->Close();

    cout << "|nAll done! :)\n" << endl;

}


//-----------------------------------------------------------------------------------
// MAIN
//-----------------------------------------------------------------------------------

void HistMaker() {

    MIRJA *classMirja = new class MIRJA();

    classMirja->runAnalysis();

}
