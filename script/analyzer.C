//////////  Reads a RAT-PAC output file and does a quick analysis ////////
/////////// Author: Ayse Bat  ///////

#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TLeaf.h>
#include <Rtypes.h>
#include <TROOT.h>
#include <TRandom3.h>
#include <TH2.h>
#include <TH3.h>
#include <TPad.h>
#include <TVector3.h>
#include <TString.h>
#include <TPRegexp.h>
#include <TGraph.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TClonesArray.h>

#include <sys/stat.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <string>
#include <vector>
#include <ctime>

//#include "rootstart.h"

// Header file for the classes stored in the TTree if any.
#include <RAT/DS/Root.hh>
#include <RAT/DS/MC.hh>
#include <RAT/DS/Calib.hh>
#include <RAT/DS/EV.hh>
#include <RAT/DS/PMT.hh>
#include <RAT/DS/PMTInfo.hh>
#include <RAT/DS/RunStore.hh>
#include <RAT/DS/Run.hh>
#include <RAT/DSReader.hh>
#include <RAT/TrackNav.hh>
#include <RAT/TrackCursor.hh>
#include <RAT/TrackNode.hh>
#include <RAT/DB.hh>

using namespace std;
using namespace TMath;

void analyzer(const char* filename_ratpac) {
	
	//Conversion Factor
	double CFactor = 16;

	// ------------------------------------------------------------------------------------------- //
	// Load RAT libraries (for dsReader)
	gSystem->Load("$(RATPAC_PATH)/lib/libRATEvent.so");
	
	// Initialization
	RAT::DSReader *dsReader;
	RAT::DS::Root *ds;
	RAT::TrackNav *nav;
	RAT::TrackCursor *cursor;
	RAT::TrackNode *node;
	TChain* tri;
	TChain* runtri;
	
	RAT::DS::Run* run;
	RAT::DS::PMTInfo* pmtInfo;
	
	std::clock_t start;
	double duration;
	
	ULong64_t NbEntries;
	
	TVector3 muon_momentum;
	TVector3 unit_vect(0.,0.,1.);
	
	TVector3 init_pos;
	TVector3 fin_pos;
	Double_t init_time, fin_time;
	TString nucl_cap_pdg_code;
	
	TH1::SetDefaultSumw2();
	//SetMyStyle();
	
	// ------------------------------------------------------------------------------------------- //
	// Starts the timer
	start = clock();
	
	gRandom = new TRandom3();
	
	// Output file
	TFile f_output("Scintillation.root","RECREATE");
	
	// Load the files
	dsReader = new RAT::DSReader(filename_ratpac);
	NbEntries = dsReader->GetTotal();
	
	// Load the trees 
	// load the ratpac trees and DS
	tri = new TChain("T");
	runtri = new TChain("runT");
	
	if (TString(filename_ratpac).MaybeWildcard()) {
		// Assume there is a runT in all files
		runtri->Add(filename_ratpac);
		RAT::DS::RunStore::SetReadTree(runtri);
	} else {
		// In single file case, we can check
		TFile *ftemp = TFile::Open(filename_ratpac);
		if (ftemp->Get("runT")) {
			runtri->Add(filename_ratpac);
			RAT::DS::RunStore::SetReadTree(runtri);
		} // else, no runT, so don't register runtri with RunStore
		
		delete ftemp;
	}
	
	RAT::DS::Root *branchDS = new RAT::DS::Root();
	tri->SetBranchAddress("ds", &branchDS);
	RAT::DS::RunStore::GetRun(branchDS);
	
	TH1::SetDefaultSumw2(kTRUE);
	
	// Create some histograms
        // mc 
	TH1F* h_numPE = new TH1F("h_numPE","Number of Photoelectron(p.e)",500,0, 500);
	// mc.summary
        TH1F* h_numScintPhoton = new TH1F("h_numScintPhoton","Number of Scintilation Photon",10000,0, 10000);	     
        TH1F* h_numCerenkovPhoton = new TH1F("h_numCerenkovPhoton","Number of Cherenkov",10000,0, 10000);
	
	// position(x,y,z) mc.particle
        TH1F* h_particleX = new TH1F("h_particleX","X (mm)",4000,-2000, 2000);
	TH1F* h_particleY = new TH1F("h_particleY","Y (mm)",4000,-2000, 2000);
        TH1F* h_particleZ = new TH1F("h_particleZ","Z (mm)",4000,-2000, 2000);
	
	TH1F* h_hitTime = new TH1F("h_hitTime","Hit Time (ns)",1000,0, 1000);
	TH1F* h_charge = new TH1F("h_Charge","Charge Distribution",1000,0, 1000);
	
	


	vector<Double_t> v_prompt_hits_times, v_delayed_hits_times;
        vector<Double_t> v_photon_charge, v_photon_hitTime;
        
	std::fstream output;
        output.open("Scintillation.csv",ios::app);
        output<<"Event"<<","<<"numPE"<<","<<"ScintPhoton"<<","<<"CerenkovPhoton"<<","<<"HitTime"<<","<<"Charge"<<","<<"x"<<","<<"y"<<","<<"z"<<","<<std::endl;
	// ------------------------------------------------------------------------------------------- //
	// Analysis loop over all the events
	for (ULong64_t entry=0; entry<NbEntries; ++entry) {
		
		// 		cout << "Entry: " << entry << endl;
		ds = dsReader->GetEvent(entry);    
		run = RAT::DS::RunStore::Get()->GetRun(ds);
		
		// Some initializations
                v_photon_charge.clear(); v_photon_hitTime.clear();
		
		// ---------- PMT loop ---------------- //
		for(long iPMT = 0; iPMT < ds->GetMC()->GetMCPMTCount(); iPMT++ ){
			for(long iPhot = 0; iPhot < ds->GetMC()->GetMCPMT(iPMT)->GetMCPhotonCount(); iPhot++){
                        	v_photon_charge.push_back(ds->GetMC()->GetMCPMT(iPMT)->GetMCPhoton(iPhot)->GetCharge());
                                v_photon_hitTime.push_back(ds->GetMC()->GetMCPMT(iPMT)->GetMCPhoton(iPhot)->GetHitTime());


				
					//std::cout<<"Pos Z: "<<ds->GetMC()->GetMCParticle(0)->GetPosition().Z()<<std::endl;
					//std::cout<<"####################: "<<std::endl;
			}
		}
		// -------------------------------------- //
                output << entry<<","<<ds->GetMC()->GetNumPE()<<","<<ds->GetMC()->GetMCSummary()->GetNumScintPhoton()<<","<<ds->GetMC()->GetMCSummary()->GetNumCerenkovPhoton()<<","<<v_photon_hitTime.size()<<","<<v_photon_charge.size()<<","<<ds->GetMC()->GetMCParticle(0)->GetPosition().x()<<","<<ds->GetMC()->GetMCParticle(0)->GetPosition().y()<<","<<ds->GetMC()->GetMCParticle(0)->GetPosition().z()<<std::endl;

		// Fills NHits histos
	        h_hitTime->Fill(v_photon_hitTime.size());
                h_charge->Fill(v_photon_charge.size());	
	
                //Fill num PE
                h_numPE->Fill(ds->GetMC()->GetNumPE());
                
               // fill Scintialation and Cheremkov Light
                h_numScintPhoton->Fill(ds->GetMC()->GetMCSummary()->GetNumScintPhoton());
                h_numCerenkovPhoton->Fill(ds->GetMC()->GetMCSummary()->GetNumCerenkovPhoton());
	
		
		h_particleX->Fill(ds->GetMC()->GetMCParticle(0)->GetPosition().x());
                h_particleY->Fill(ds->GetMC()->GetMCParticle(0)->GetPosition().y());
                h_particleZ->Fill(ds->GetMC()->GetMCParticle(0)->GetPosition().z());

	}
        output.close();
	// ------------------------------------------------------------------------------------------- //
	
	delete run;
	delete tri, runtri, branchDS;
	delete dsReader;
	
	
	f_output.cd();
       
        h_numPE->SetXTitle("Photoelectrons");
        h_numScintPhoton->SetXTitle("Scitillation");
        h_numCerenkovPhoton->SetXTitle("Cherenkov");
        h_particleX->SetXTitle("X (mm)");
        h_particleY->SetXTitle("Y(mm)");
        h_particleZ->SetXTitle("Z(mm)");
        h_hitTime->SetXTitle("Hit Time (ns)");
        h_charge->SetXTitle("Charge");
   
 	h_numPE->SetYTitle("Events");
        h_numScintPhoton->SetYTitle("Events");
        h_numCerenkovPhoton->SetYTitle("Events");
        h_particleX->SetYTitle("Events");
        h_particleY->SetYTitle("Events");
        h_particleZ->SetYTitle("Events");
        h_hitTime->SetYTitle("Events");
        h_charge->SetYTitle("Events");


 
	h_numPE->SetOption("hist");
        h_numScintPhoton->SetOption("hist");
        h_numCerenkovPhoton->SetOption("hist");
        h_particleX->SetOption("hist");
        h_particleY->SetOption("hist");
        h_particleZ->SetOption("hist");
        h_hitTime->SetOption("hist");
        h_charge->SetOption("hist");

	
        h_numPE->Write();	
        h_numScintPhoton->Write();
        h_numCerenkovPhoton->Write();
        h_particleX->Write();
        h_particleY->Write();
        h_particleZ->Write();
        h_hitTime->Write();
        h_charge->Write();


	f_output.Close();
	
	// Ends the timer
	duration = ( clock() - start ) / (double) CLOCKS_PER_SEC;
	cout << "Execution time: " << duration << " seconds\n";
	
}
