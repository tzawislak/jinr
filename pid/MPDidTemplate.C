/* Macro reads DST file produced by macro reco.C */
/* Standard readDST.root macro modified by Tomasz Zawislak */

#include <Rtypes.h>
#include <TChain.h>
#include <TClonesArray.h>
#include <TStopwatch.h>

R__ADD_INCLUDE_PATH($VMCWORKDIR)
#include "macro/mpd/mpdloadlibs.C"

TChain* readAllRootFiles( TString &in, Int_t &nof){
	TChain *dstTree = new TChain("cbmsim");
	//Add all exsisting mpddst.root files to chain
	ifstream files(in);
	if(!files.good()){
		cout << "Error while loading files" << endl;
	}
	Int_t job_name=0;
	for( Int_t iter=0; iter<=nof; ++iter){ 
		files>>job_name;
		TString fileName = "/eos/nica/mpd/data/Urqmd.4GeV/" + to_string(job_name) + "/TestEvUrqmd/mpddst.root";
		dstTree->Add(fileName.Data());
	}
	files.close();
	return dstTree;
}

void MPDidTemplate(TString in = "/eos/nica/mpd/users/zawislak/pid/4gev/4gev.txt", Int_t nOfFiles=300 ) {
        
    TStopwatch timer;    timer.Start();

    TChain *dstTree = readAllRootFiles(in, nOfFiles);
   
    // Activate branches
    MpdEvent *event = nullptr;
    dstTree->SetBranchAddress("MPDEvent.", &event);
    TClonesArray *fMCTracks = nullptr;
    dstTree->SetBranchAddress("MCTrack", &fMCTracks);

    Int_t events = dstTree->GetEntries();
    cout << " Number of events in DST file = " << events << endl;
	
    //Declare your histograms here:
	const Int_t NpBins=6;
	const Int_t Nparts=4;
	const Float_t ptb[7] = {0, 0.2, 0.5, 0.9, 1.5, 2.2, 3.};
	const string parts[Nparts] = {"pions", "protons", "kaons", "electrons"};
	const Int_t pid[Nparts] = {211, 2212, 321, 11};
	
    TH2F* dEdx[2][Nparts];
    TH2F* tpctof[2][Nparts][NpBins];
    TH2F* tof[2][Nparts];
    TH2F* beta[2][Nparts];
    
    for( Int_t i=0; i<Nparts; ++i){
		dEdx[0][i] = new TH2F( "dedxp", "Positive particles;log(p) [log(GeV/c)];dE/dx [a.u]", 600, -2, 10, 1000, 0, 100000); 
		dEdx[1][i] = new TH2F( "dedxn", "Negative particles;log(p) [log(GeV/c)];dE/dx [a.u]", 600, -2, 10, 1000, 0, 100000); 
		tof[0][i] = new TH2F( "tofp", "Positive particles; p [GeV/c]; m^{2};", 500, 0, 5, 600, -1, 2); 
		tof[1][i] = new TH2F( "tofn", "Negative particles; p [GeV/c]; m^{2};", 500, 0, 5, 600, -1, 2); 
		beta[0][i] = new TH2F( "betap", "Positive particles; p [GeV/c]; #beta;",  500, 0, 5, 500, -0.1, 1.9); 
		beta[1][i] = new TH2F( "betan", "Negative particles; p [GeV/c]; #beta;",  500, 0, 5, 500, -0.1, 1.9); 
		for( Int_t i0=0; i0<NpBins; ++i0){
			tpctof[0][i][i0] = new TH2F( Form("tpctofp_%d", i0), Form("Pos p #in (%f: %f) GeV; dE/dx [a.u]; m^{2}", ptb[i0], ptb[i0+1]),  1000, 0, 100000, 600, -1, 2); 
			tpctof[1][i][i0] = new TH2F( Form("tpctofn_%d", i0), Form("Neg p #in (%f: %f) GeV; dE/dx [a.u]; m^{2}", ptb[i0], ptb[i0+1]),  1000, 0, 100000, 600, -1, 2); 
		}	
	}
	TH1I* err = new TH1I( "nans", "nans", 12, -0.5, 11.5);
	
    //Event loop

    for (Int_t i = 0; i < events; i++) {
        dstTree->GetEntry(i);
        for (Int_t iTrack = 0; iTrack < event->GetGlobalTracks()->GetEntriesFast(); iTrack++) {
            MpdTrack* track = (MpdTrack*) event->GetGlobalTracks()->UncheckedAt(iTrack);
            auto mctrack = (FairMCTrack*) fMCTracks->UncheckedAt(track->GetID());
	        auto pdgID =  mctrack->GetPdgCode();
	    //for templates one needs only pid[] IDs
			Int_t ID=-1;		//to iterate through arrays
			for( Int_t p=0; p<Nparts; ++p){
				if(abs(pdgID)==pid[p]){
					ID=p;
					break;
				}
			}
			if(ID==-1) continue;
		//ToF existance flag 
			if( track->GetTofFlag() != 6 ) continue;  
		//Set variables
			Int_t CH=0, PTOT=0;
			if(track->GetCharge()==0 ) continue;
			if(track->GetCharge()<0 )  CH=1;
			Float_t pt=track->GetPt(); 
			Float_t pz=track->GetPz();
			Float_t ptotal = sqrt(pz*pz+pt*pt);
			Float_t dedx=track->GetdEdXTPC();
			Float_t m2=track->GetTofMass2();
			Float_t vbeta=track->GetTofBeta();
		
		//  Search for appropriate bin in p_total
			for( Int_t i1=0; i1<NpBins; ++i1){
				if( ptotal > 3 ){ 
					PTOT=5;
					break;
				}
				if( ptb[i1]<=ptotal && ptotal<ptb[i1+1]){
					PTOT=i1;
					break;
				}
			}
	
				tof[CH][ID]->Fill( ptotal , m2 );
				tpctof[CH][ID][PTOT]->Fill( dedx, m2 );
				dEdx[CH][ID]->Fill(  ptotal , dedx );
				beta[CH][ID]->Fill(  ptotal , vbeta);
			/* See mpddata/MpdTrack.h for more methods */
        } // track loop
    } // event loop

	//Create TFile to save histograms:
	TFile *tfhistos;
	string templPath="/eos/nica/mpd/users/zawislak/pid/4gev/raw/";
	for( Int_t i=0; i<Nparts; ++i){
		string tpath = templPath + parts[i] + ".root";
		tfhistos = new TFile(tpath.c_str(), "recreate");
		
		for( Int_t k=0; k<2; ++k){
			dEdx[k][i]->Write();
			tof[k][i]->Write();
			beta[k][i]->Write();
			for( Int_t i2=0; i2<NpBins; ++i2)
				tpctof[k][i][i2]->Write();
		}
		tfhistos->Close();
	}
		
	
    timer.Print();
    exit(0);
}
