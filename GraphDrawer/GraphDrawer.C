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

    //-------------------------------------------------------------------------------
    // CONSTRUCTORS
    //-------------------------------------------------------------------------------


    LIMI() { // Default constructor

        // gStyle
        gStyle->SetOptTitle(0);
        gStyle->SetOptStat(0);
        gStyle->SetOptLogx(1);

        // Bin properties
        for (int i=0; i < binNumber; i++) {
            binErrors[i] = (binEdges[i + 1] - binEdges[i]) / 2;
            binCentres[i] = (binEdges[i + 1] + binEdges[i]) / 2;
        }

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
    // (0/?)
    //-------------------------------------------------------------------------------
    cout << "Trying to load Histogram file..." << endl;

    // Retrieve histogram ROOT file
    TFile *histFile = new TFile("../ZoneLoader/ZoneAMS02.root");
    TFile *mcFile   = new TFile("../HistMaker/ProtonHistogramsAMS02.root");

    cout << "   ...File loaded!" << endl;


    //-------------------------------------------------------------------------------
    // (1/?)
    //-------------------------------------------------------------------------------
    cout << "Creating Events graph... (1/?)" << endl;

    // Create histograms
    TH1F *hEventsDetected = new TH1F();
    TH1F *hEventsSelected = new TH1F();

    // Get relevant histograms
    hEventsDetected = (TH1F*)histFile->Get("eventsDetected");
    hEventsSelected = (TH1F*)histFile->Get("eventsSelected");

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
    cEvents->Print("Events.png");


    //-------------------------------------------------------------------------------
    // (2/?)
    //-------------------------------------------------------------------------------
    cout << "Creating ExposureTime graph... (2/?)" << endl;

    // Create histograms
    TH1F *hExposureTime = new TH1F();

    // Get relevant histograms
    hExposureTime = (TH1F*)histFile->Get("exposureTime");

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
    gExposureTime->SetMarkerStyle(20);
    gExposureTime->SetMarkerSize(1);
    gExposureTime->SetMarkerColor(kRed);

    // Axes
    cExposureTime->SetLogy();
    gExposureTime->GetXaxis()->SetTitle("R [GV]");
    gExposureTime->GetYaxis()->SetTitle("Exposure Time [s]");

    // Print
    cExposureTime->Draw();
    cExposureTime->Print("ExposureTime.png");


    //-------------------------------------------------------------------------------
    // (3/?)
    //-------------------------------------------------------------------------------
    cout << "Creating Acceptance graph... (3/?)" << endl;

    // Create histograms
    TH1F *hMCGenerated = new TH1F();
    TH1F *hMCDetected  = new TH1F();
    TH1F *hMCSelected  = new TH1F();

    // Get relevant histograms
    hMCGenerated = (TH1F*)mcFile->Get("montecarloGenerated");
    hMCDetected  = (TH1F*)mcFile->Get("montecarloDetected");
    hMCSelected  = (TH1F*)mcFile->Get("montecarloSelected");

    // Set arrays
    double acceptanceDetected[binNumber];
    double acceptanceSelected[binNumber];

    // Fill arrays
    for (int i=0; i < binNumber; i++) {

        acceptanceDetected[i] = TMath::Pi() * 3.9 * 3.9 * hMCDetected->GetBinContent(i + 1) / hMCGenerated->GetBinContent(i + 1);
        acceptanceSelected[i] = TMath::Pi() * 3.9 * 3.9 * hMCSelected->GetBinContent(i + 1) / hMCGenerated->GetBinContent(i + 1);

    }

    // Create TGraphs
    TGraphErrors *gAcceptanceDetected = new TGraphErrors(binNumber, binCentres, acceptanceDetected, binErrors, 0);
    TGraphErrors *gAcceptanceSelected = new TGraphErrors(binNumber, binCentres, acceptanceSelected, binErrors, 0);

    // Create TCanvas
    TCanvas *cAcceptance = new TCanvas("cAcceptance", "Acceptance per Rigidity Bin");
    gAcceptanceDetected->Draw("AP");
    gAcceptanceSelected->Draw("P");

    // Styling
    gAcceptanceDetected->SetMarkerStyle(20);
    gAcceptanceDetected->SetMarkerSize(1);
    gAcceptanceDetected->SetMarkerColor(kBlack);
    gAcceptanceSelected->SetMarkerStyle(20);
    gAcceptanceSelected->SetMarkerSize(1);
    gAcceptanceSelected->SetMarkerColor(kBlue);

    // Axes
    gAcceptanceDetected->SetMinimum(0);
    gAcceptanceDetected->GetXaxis()->SetTitle("R [GV]");
    gAcceptanceDetected->GetYaxis()->SetTitle("Acceptance [m^2 sr]");

    // Print
    cAcceptance->Draw();
    cAcceptance->Print("Acceptance.png");


    //-------------------------------------------------------------------------------
    // (4/?)
    //-------------------------------------------------------------------------------
    cout << "Creating TriggerEfficiency graph... (4/?)" << endl;

    // Create histograms
    TH1F *hTriggersPhysical = new TH1F();
    TH1F *hTriggersBias     = new TH1F();
    TH1F *hMCPhysical       = new TH1F();
    TH1F *hMCBias           = new TH1F();

    // Get relevant histograms
    hTriggersPhysical = (TH1F*)histFile->Get("triggersPhysical");
    hTriggersBias     = (TH1F*)histFile->Get("triggersBias");
    hMCPhysical       = (TH1F*)mcFile->Get("montecarloPhysical");
    hMCBias           = (TH1F*)mcFile->Get("montecarloBias");

    // Set arrays
    double triggerEfficiency[binNumber]; double triggerEfficiencyErrors[binNumber];
    double mcTriggerEfficiency[binNumber]; double mcTriggerEfficiencyErrors[binNumber];

    // Fill arrays
    for (int i=0; i < binNumber; i++) {

        triggerEfficiency[i]      = hTriggersPhysical->GetBinContent(i + 1) / (hTriggersPhysical->GetBinContent(i + 1) + 100 * hTriggersBias->GetBinContent(i + 1));
        mcTriggerEfficiency[i]    = hMCPhysical->GetBinContent(i + 1) / (hMCPhysical->GetBinContent(i + 1) + hMCBias->GetBinContent(i + 1));

        triggerEfficiencyErrors[i]      = 100 * TMath::Sqrt( hTriggersPhysical->GetBinContent(i + 1) * pow(hTriggersBias->GetBinContent(i + 1), 2)
                                          + hTriggersBias->GetBinContent(i + 1) * pow(hTriggersPhysical->GetBinContent(i + 1), 2) )
                                          / pow( hTriggersPhysical->GetBinContent(i + 1) + 100 * hTriggersBias->GetBinContent(i + 1) , 2);
        mcTriggerEfficiencyErrors[i]    = TMath::Sqrt( hTriggersPhysical->GetBinContent(i + 1) * pow(hTriggersBias->GetBinContent(i + 1), 2)
                                          + hTriggersBias->GetBinContent(i + 1) * pow(hTriggersPhysical->GetBinContent(i + 1), 2) )
                                          / pow( hTriggersPhysical->GetBinContent(i + 1) + hTriggersBias->GetBinContent(i + 1) , 2);

    }

    // Create TGraphs
    TGraphErrors *gTriggerEfficiency      = new TGraphErrors(binNumber, binCentres, triggerEfficiency, binErrors, triggerEfficiencyErrors);
    TGraphErrors *gMCTriggerEfficiency    = new TGraphErrors(binNumber, binCentres, mcTriggerEfficiency, binErrors, mcTriggerEfficiencyErrors);

    // Create TCanvas
    TCanvas *cTriggerEfficiency = new TCanvas("cTriggerEfficiency", "Trigger Efficiency per Rigidity Bin");
    gTriggerEfficiency->Draw("AP");
    gMCTriggerEfficiency->Draw("P");

    // Styling
    gTriggerEfficiency->SetMarkerStyle(20);
    gTriggerEfficiency->SetMarkerSize(1);
    gTriggerEfficiency->SetMarkerColor(kBlue);
    gMCTriggerEfficiency->SetMarkerStyle(20);
    gMCTriggerEfficiency->SetMarkerSize(1);
    gMCTriggerEfficiency->SetMarkerColor(kRed);

    // Axes
    gTriggerEfficiency->SetMaximum(1.0);
    gTriggerEfficiency->SetMinimum(0.0);
    gTriggerEfficiency->GetXaxis()->SetTitle("R [GV]");
    gTriggerEfficiency->GetYaxis()->SetTitle("Trigger Efficiency");

    // Print
    cTriggerEfficiency->Draw();
    cTriggerEfficiency->Print("Trigger Efficiency.png");


    //-------------------------------------------------------------------------------
    // (5/?)
    //-------------------------------------------------------------------------------
    cout << "Creating SelectionEfficiency graph... (5/?)" << endl;

    // Create histograms
    TH1F *hTracker      = new TH1F();
    TH1F *hTOF          = new TH1F();
    TH1F *hParticle     = new TH1F();
    TH1F *hBeta         = new TH1F();
    TH1F *hChiSquared   = new TH1F();
    TH1F *hInnerLayer   = new TH1F();
    TH1F *hMCTracker    = new TH1F();
    TH1F *hMCTOF        = new TH1F();
    TH1F *hMCParticle   = new TH1F();
    TH1F *hMCBeta       = new TH1F();
    TH1F *hMCChiSquared = new TH1F();
    TH1F *hMCInnerLayer = new TH1F();

    // Get relevant histograms
    hTracker      = (TH1F*)histFile->Get("triggersPhysical");
    hTOF          = (TH1F*)histFile->Get("baseTOF");
    hParticle     = (TH1F*)histFile->Get("cutParticle");
    hBeta         = (TH1F*)histFile->Get("cutBeta");
    hChiSquared   = (TH1F*)histFile->Get("cutChiSquared");
    hInnerLayer   = (TH1F*)histFile->Get("cutInnerLayer");
    hMCTracker    = (TH1F*)mcFile->Get("montecarloTracker");
    hMCTOF        = (TH1F*)mcFile->Get("montecarloTOF");
    hMCParticle   = (TH1F*)mcFile->Get("montecarloParticle");
    hMCBeta       = (TH1F*)mcFile->Get("montecarloBeta");
    hMCChiSquared = (TH1F*)mcFile->Get("montecarloChiSquared");
    hMCInnerLayer = (TH1F*)mcFile->Get("montecarloInnerLayer");

    // Normalise by corresponding instrument base
    hParticle->Divide(hTOF);
    hBeta->Divide(hTOF);
    hChiSquared->Divide(hTracker);
    hInnerLayer->Divide(hTracker);
    hMCParticle->Divide(hMCTOF);
    hMCBeta->Divide(hMCTOF);
    hMCChiSquared->Divide(hMCTracker);
    hMCInnerLayer->Divide(hMCTracker);

    // Set arrays
    double selectionEfficiency[binNumber];   // double selectionEfficiencyErrors[binNumber];
    double mcSelectionEfficiency[binNumber]; // double mcSelectionEfficiencyErrors[binNumber];

    // Fill arrays
    for (int i=0; i < binNumber; i++) {

        selectionEfficiency[i]   = hParticle->GetBinContent(i + 1) * hBeta->GetBinContent(i + 1) * hChiSquared->GetBinContent(i + 1) * hInnerLayer->GetBinContent(i + 1);
        mcSelectionEfficiency[i] = hMCParticle->GetBinContent(i + 1) *hMCBeta->GetBinContent(i + 1) * hMCChiSquared->GetBinContent(i + 1) * hMCInnerLayer->GetBinContent(i + 1);

    }

    // Create TGraphs
    TGraphErrors *gSelectionEfficiency = new TGraphErrors(binNumber, binCentres, selectionEfficiency, binErrors, 0);
    TGraphErrors *gMCSelectionEfficiency = new TGraphErrors(binNumber, binCentres, mcSelectionEfficiency, binErrors, 0);

    // Create TCanvas
    TCanvas *cSelectionEfficiency = new TCanvas("cSelectionEfficiency", "Selection Efficiency per Rigidity Bin");
    gSelectionEfficiency->Draw("AP");
    gMCSelectionEfficiency->Draw("P");

    // Styling
    gSelectionEfficiency->SetMarkerStyle(20);
    gSelectionEfficiency->SetMarkerSize(1);
    gSelectionEfficiency->SetMarkerColor(kBlue);
    gMCSelectionEfficiency->SetMarkerStyle(20);
    gMCSelectionEfficiency->SetMarkerSize(1);
    gMCSelectionEfficiency->SetMarkerColor(kRed);

    // Axes
    gSelectionEfficiency->SetMaximum(1.0);
    gSelectionEfficiency->SetMinimum(0.0);
    gSelectionEfficiency->GetXaxis()->SetTitle("R [GV]");
    gSelectionEfficiency->GetYaxis()->SetTitle("Selection Efficiency");

    // Print
    cSelectionEfficiency->Draw();
    cSelectionEfficiency->Print("Selection Efficiency.png");


    //-------------------------------------------------------------------------------
    // (6/?)
    //-------------------------------------------------------------------------------
    cout << "Creating ProtonRate graph... (6/?)" << endl;

    // Set arrays
    double rate[binNumber]; double rateErrors[binNumber];

    // Fill arrays
    for (int i=0; i < binNumber; i++) {

        if (exposureTime[i] == 0) {

            rate[i] = 0;

            rateErrors[i] = 0;

        } else {

            rate[i] = eventsSelected[i] / exposureTime[i] / (2 * binErrors[i]);

            rateErrors[i] = rate[i] / TMath::Sqrt(eventsSelected[i]);
            
        }

    }

    // Create TGraphs
    TGraphErrors *gRate = new TGraphErrors(binNumber, binCentres, rate, binErrors, rateErrors);

    // Create TCanvas
    TCanvas *cRate = new TCanvas("cRate", "Proton Rate as a function of Rigidity");
    gRate->Draw("AP");

    // Styling
    gRate->SetMarkerStyle(20);
    gRate->SetMarkerSize(1);
    gRate->SetMarkerColor(kRed);

    // Axes
    cRate->SetLogy();
    gRate->GetXaxis()->SetTitle("R [GV]");
    gRate->GetYaxis()->SetTitle("Rate [s^-1]");

    // Print
    cRate->Draw();
    cRate->Print("Proton Rate.png");


    //-------------------------------------------------------------------------------
    // (7/?)
    //-------------------------------------------------------------------------------
    cout << "Creating ProtonFlux graph... (7/?)" << endl;

    // Set arrays
    double flux[binNumber]; double fluxErrors[binNumber];

    // Fill arrays
    for (int i=0; i < binNumber; i++) {

        if (exposureTime[i] == 0) {

            flux[i] = 0;

            fluxErrors[i] = 0;

        } else {

            flux[i] = rate[i] / acceptanceSelected[i] 
                      / triggerEfficiency[i] * mcTriggerEfficiency[i]
                      / selectionEfficiency[i] * mcSelectionEfficiency[i];

            fluxErrors[i] = flux[i] * TMath::Sqrt( 1 / eventsSelected[i]
                            + pow(triggerEfficiencyErrors[i] / triggerEfficiency[i], 2) 
                            + pow(mcTriggerEfficiencyErrors[i] / mcTriggerEfficiency[i], 2) );

        }

    }

    // Create TGraphs
    TGraphErrors *gFlux = new TGraphErrors(binNumber, binCentres, flux, binErrors, fluxErrors);

    // Create TCanvas
    TCanvas *cFlux = new TCanvas("cFlux", "Proton Flux as a function of Rigidity");
    gFlux->Draw("AP");

    // Styling
    gFlux->SetMarkerStyle(20);
    gFlux->SetMarkerSize(1);
    gFlux->SetMarkerColor(kRed);

    // Axes
    cFlux->SetLogy();
    gFlux->GetXaxis()->SetTitle("R [GV]");
    gFlux->GetYaxis()->SetTitle("Flux [m^-2 sr^-1 s^-1 GV^-1]");

    // Print
    cFlux->Draw();
    cFlux->Print("Proton Flux.png");


    //-------------------------------------------------------------------------------
    // (8/?)
    //-------------------------------------------------------------------------------
    cout << "Creating ProtonScaledFlux graph... (8/?)" << endl;

    // Set arrays
    double scaledFlux[binNumber]; double scaledFluxErrors[binNumber];

    // Fill arrays
    for (int i=0; i < binNumber; i++) {

        scaledFlux[i] = flux[i] * pow(binCentres[i], 2.7);

        scaledFluxErrors[i] = fluxErrors[i] * pow(binCentres[i], 2.7);

    }

    // Create TGraphs
    TGraphErrors *gScaledFlux = new TGraphErrors(binNumber, binCentres, scaledFlux, binErrors, scaledFluxErrors);

    // Create TCanvas
    TCanvas *cScaledFlux = new TCanvas("cScaledFlux", "Scaled Proton Flux as a function of Rigidity");
    gScaledFlux->Draw("AP");

    // Styling
    gScaledFlux->SetMarkerStyle(20);
    gScaledFlux->SetMarkerSize(1);
    gScaledFlux->SetMarkerColor(kRed);

    // Axes
    cScaledFlux->SetLogy();
    gScaledFlux->GetXaxis()->SetTitle("R [GV]");
    gScaledFlux->GetYaxis()->SetTitle("Flux R^2.7 [m^-2 sr^-1 s^-1 GV^1.7]");

    // Print
    cScaledFlux->Draw();
    cScaledFlux->Print("Scaled Proton Flux.png");


    cout << "\nAll done! :)\n" << endl;

}


//-----------------------------------------------------------------------------------
// MAIN
//-----------------------------------------------------------------------------------

void GraphDrawer() {

    LIMI *classLimi = new class LIMI();

    classLimi->run();

}
