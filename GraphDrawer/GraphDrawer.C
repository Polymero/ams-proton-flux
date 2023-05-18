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
#include "TObject.h"
#include "TString.h"
// Local headers
#include "../Header Files/Ntp.h"


//-----------------------------------------------------------------------------------
// CLASS DEFINITION
//-----------------------------------------------------------------------------------

// Class
class LIMI {
    public: // Access specifier


    //-------------------------------------------------------------------------------
    // ATTRIBUTES
    //-------------------------------------------------------------------------------

    // Rigidity bins
    const int Bin_num = 32;
    double binEdges[32 + 1] = {
        1.00, 1.16, 1.33, 1.51, 1.71, 1.92, 2.15, 2.40, 2.67, 2.97, 3.29, 3.64, 4.02,
        4.43, 4.88, 5.37, 5.90, 6.47, 7.09, 7.76, 8.48, 9.26, 10.1, 11.0, 12.0, 13.0,
        14.1, 15.3, 16.6, 18.0, 19.5, 21.1, 22.8
    };
    double binErrors[32];
    double binCentres[32];

    // Rigidity cut-off level
    double Rig_Cut_Level = 1.2;

    //-------------------------------------------------------------------------------
    // CONSTRUCTORS
    //-------------------------------------------------------------------------------


    LIMI() { // Default constructor

        // gStyle
        gStyle->SetOptTitle(0);
        gStyle->SetOptStat(0);
        gStyle->SetOptLogx(1);

        // Bin properties
        for (int i=0; i<Bin_num; i++) {
            binErrors[i] = (binEdges[i + 1] - binEdges[i]) / 2;
            binCentres[i] = (binEdges[i + 1] + binEdges[i]) / 2;
        }

        // Retrieve histogram ROOT file $ Hischaajat.C
        cout << "Trying to load Histogram file..." << flush
        TFile *histFile = new TFile("../ZoneLoader/ZoneAMS02.root");
        cout << "   ...File loaded!" << endl;

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
void LIMI::run() {

    cout << "Starting LIMI.run()..." << endl;


    //-------------------------------------------------------------------------------
    // (1/?)
    //-------------------------------------------------------------------------------
    cout << "Creating Events graph... (1/?)" << endl;

    // Get relevant histograms
    hEventsDetected = (TH1F*) histFile->Get("eventsDetected");
    hEventsSelected = (TH1F*) histFile->Get("eventsSelected");

    // Set arrays
    double eventsDetected[binNumber]; double eventsDetectedErrors[binNumber];
    double eventsSelected[binNumber]; double eventsSelectedErrors[binNumber];

    // Fill arrays
    for (int i=0; i < binNumber; i++) {

        eventsDetected[i] = hEventsDetected->GetBinContent(i + 1);
        eventsSelected[i] = hEventsSelected->GetBinContent(i + 1);

        eventsDetectedErrors[i] = TMath::Sqrt(eventsDetected[i]);
        eventsSelectedErrors[i] = TMath::Sqrt(eventsSelected[i]);

    }

    // Create TGraphs
    TGraphErrors *gEventsDetected = new TGraphErrors(binNumber, binCentres, eventsDetected, binErrors, eventsDetectedErrors);
    TGraphErrors *gEventsSelected = new TGraphErrors(binNumber, binCentres, eventsSelected, binErrors, eventsSelectedErrors);

    // Create TCanvas
    TCanvas *cEvents = new TCanvas("cEvents", "Events per Rigidity Bin");
    gEventsDetected->Draw("AP");
    gEventsSelected->Draw("P");

    // Styling
    gEventsDetected->SetMarkerStyle(20);
    gEventsDetected->SetMarkerSize(1);
    gEventsDetected->SetMarkerColor(kBlack);
    gEventsSelected->SetMarkerStyle(20);
    gEventsSelected->SetMarkerSize(1);
    gEventsSelected->SetMarkerColor(kBlue);

    // Axes
    cEvents->SetLogy();
    gEventsDetected->SetMinimum(1);
    gEventsDetected->GetXaxis()->SetTitle("R [GV]");
    gEventsDetected->GetYaxis()->SetTitle("Events");

    // Print
    cEvents->Draw();
    cEvents->print(("Events.png").c_str());


    //-------------------------------------------------------------------------------
    // (2/?)
    //-------------------------------------------------------------------------------
    cout << "Creating ExposureTime graph... (2/?)" << endl;

    // Get relevant histograms
    hExposureTime = (TH1F*) histFile->Get("exposureTime");

    // Set arrays
    double exposureTime[binNumber];

    // Fill arrays
    for (int i=0; i < binNumber; i++) {
        exposureTime[i] = hExposureTime->GetBinContent(i + 1);
    }

    // Create TGraphs
    TGraphErrors *gExposureTime = new TGraphErrors(binNumber, binCentres, exposureTime, binErrors, 0);

    // Create TCanvas
    TCanvas *cExposureTime = new TCanvas("cExposureTime", "Exposure Time per Rigidity Bin");
    gExposureTime->Draw("AP");

    // Styling
    gExposureTime->SetMarkerSize(20);
    gExposureTime->SetMarkerSize(1);
    gExposureTime->SetMarkerColor(kRed);

    // Axes
    cExposureTime->SetLogy();
    gExposureTime->GetXaxis()->SetTitle("R [GV]");
    gExposureTime->GetYaxis()->SetTitle("Exposure Time [s]");

    // Print
    cExposureTime->Draw();
    cExposureTime->Print(("ExposureTime.png").c_str());


    //-------------------------------------------------------------------------------
    // (3/?)
    //-------------------------------------------------------------------------------
    // cout << "Creating Acceptance graph... (3/?)" << endl;




    cout << "\nAll done! :)\n" << endl;

}


//-----------------------------------------------------------------------------------
// MAIN
//-----------------------------------------------------------------------------------

void GraphDrawer() {

    LIMI *classLimi = new class LIMI();

    classLimi->run();

}
