// Written by Sebastiaan Venendaal (University of Groningen, the Netherlands)
// C++ class for generating histograms of proton-like AMS-02 data, used for flux analysis
// Created          16-05-23
// Last modified    16-05-23
//
// Usage ::
// Condor job submission through SSH connection to CERN

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
#include "TObject.h"
#include "TString.h"
#include "TSystem.h"
// Local headers
#include "../Header Files/Ntp.h"


//-----------------------------------------------------------------------------------
// CLASS DEFINITION
//-----------------------------------------------------------------------------------

class MIRJA {
    // Access specifier
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

    // List of data objects
    // Chains
    TChain *chainCompact        = new TChain("Compact");
    TChain *chainRTI            = new TChain("RTI");
    // Classes
    NtpCompact *classCompact    = new NtpCompact();
    NtpSHeader *classSHeader    = new class NtpSHeader();
    RTIInfo *classRTI           = new class RTIInfo();
    // Files
    TFile *f = new TFile();


    //-------------------------------------------------------------------------------
    // CLASS CONSTRUCTORS
    //-------------------------------------------------------------------------------

    MIRJA(int zoneIndex) { // Default constructor

        // ROOT gStyle configuration
        gStyle->SetOptTitle(0);
        gStyle->SetOptStat(0);
        gStyle->SetOptLogx(1);

        // Bin properties
        for (int i=0; i < binNumber; i++) {
            binErrors[i]  = (binEdges[i + 1] - binEdges[i]) / 2;
            binCentres[i] = (binEdges[i + 1] + binEdges[i]) / 2;
        }

        // New file object
        f = (TFile*)TFile::Open(Form("/afs/cern.ch/user/s/svenenda/public/ams-proton-flux/ZoneLoader/AMS02Zone%d.root", zoneIndex), "recreate");

        // Get correct root files according to the zone index
        // int utcint[129] = {
        //     1307499168, 1309717509, 1311935851, 1314154192, 1316372533, 1318590875, 1320809216, 
        //     1323027558, 1325245899, 1327464240, 1329682582, 1331900923, 1334119264, 1336337606, 
        //     1338555947, 1340774288, 1342992630, 1345210971, 1347429312, 1349647654, 1351865995, 
        //     1354084337, 1356302678, 1358521019, 1360739361, 1362957702, 1365176043, 1367394385, 
        //     1369612726, 1371831067, 1374049409, 1376267750, 1378486092, 1380704433, 1382922774, 
        //     1385141116, 1387359457, 1389577798, 1391796140, 1394014481, 1396232822, 1398451164, 
        //     1400669505, 1402887846, 1405106188, 1407324529, 1409542871, 1411761212, 1413979553, 
        //     1416197895, 1418416236, 1420634577, 1422852919, 1425071260, 1427289601, 1429507943, 
        //     1431726284, 1433944625, 1436162967, 1438381308, 1440599650, 1442817991, 1445036332, 
        //     1447254674, 1449473015, 1451691356, 1453909698, 1456128039, 1458346380, 1460564722, 
        //     1462783063, 1465001405, 1467219746, 1469438087, 1471656429, 1473874770, 1476093111, 
        //     1478311453, 1480529794, 1482748135, 1484966477, 1487184818, 1489403159, 1491621501, 
        //     1493839842, 1496058184, 1498276525, 1500494866, 1502713208, 1504931549, 1507149890, 
        //     1509368232, 1511586573, 1513804914, 1516023256, 1518241597, 1520459938, 1522678280, 
        //     1524896621, 1527114963, 1529333304, 1531551645, 1533769987, 1535988328, 1538206669, 
        //     1540425011, 1542643352, 1544861693, 1547080035, 1549298376, 1551516718, 1553735059, 
        //     1555953400, 1558171742, 1560390083, 1562608424, 1564826766, 1567045107, 1569263448, 
        //     1571481790, 1573700131, 1575918472, 1578136814, 1580355155, 1582573497, 1584791838, 
        //     1587010179, 1589228521, 1591446862
        // };

        int utcint[3] = {
            1309704319, 1309705723, 1309707106
        };

        // Read the data trees
        for (int i = utcint[zoneIndex]; i < utcint[zoneIndex + 1]; i++) {

            if (gSystem->AccessPathName(Form("/eos/ams/group/dbar/release_v7/e1_vdev_200421/neg/ISS.B1130/pass7/%d.root", i))) {

                ;

            } else {

                chainCompact->Add(Form("/eos/ams/group/dbar/release_v7/e1_vdev_200421/neg/ISS.B1130/pass7/%d.root", i));
                chainRTI->Add(Form("/eos/ams/group/dbar/release_v7/e1_vdev_200421/neg/ISS.B1130/pass7/%d.root", i));

            }

        }

        // Set branch addresses
        chainCompact->SetBranchAddress("Compact", &classCompact);
        chainCompact->SetBranchAddress("SHeader", &classSHeader);
        chainRTI->SetBranchAddress("RTIInfo", &classRTI);

        cout << "\nClass succesfully constructed!\n" << endl;

    };


    //-------------------------------------------------------------------------------
    // CLASS METHODS
    //-------------------------------------------------------------------------------

    void run();

};


//-----------------------------------------------------------------------------------
// SUPPORT FUNCTIONS
//-----------------------------------------------------------------------------------

// No support functions


//-----------------------------------------------------------------------------------
// METHOD FUNCTIONS
//-----------------------------------------------------------------------------------
void MIRJA::run() {

    cout << "Starting MIRJA.run()..." << endl;


    //-------------------------------------------------------------------------------
    // (1/2)
    //-------------------------------------------------------------------------------
    cout << "Looping over RTIInfo data... (1/2)" << endl;

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

    }


    //-------------------------------------------------------------------------------
    // (2/2)
    //-------------------------------------------------------------------------------
    cout << "\nLooping over Compact data... (2/2)" << endl;

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
        
        // Selection parameter
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

    }
 

    //-------------------------------------------------------------------------------
    // SAVE
    //-------------------------------------------------------------------------------
    cout << "\nSaving all my hard work..." << endl;

    f->Write();
    f->Close();

    cout << "\nAll done! :)\n" << endl;

}


//-----------------------------------------------------------------------------------
// MAIN
//-----------------------------------------------------------------------------------

void ZoneLooper(int zoneIndex) {

    MIRJA *classMirja = new class MIRJA(zoneIndex);

    classMirja->run();

}
