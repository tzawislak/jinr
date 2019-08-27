/* Macro reads DST file produced by macro reco.C */
/* Standard readDST.root macro modified by Tomasz Zawislak */

#include <Rtypes.h>
#include <TChain.h>
#include <TClonesArray.h>
#include <TStopwatch.h>
#include "MPDidentify.h"

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
	for( Int_t iter=0; iter<nof; ++iter){ 
		files>>job_name;
		TString fileName = "/eos/nica/mpd/data/Urqmd.11GeV/" + to_string(job_name) + "/TestEvUrqmd/mpddst.root";
		dstTree->Add(fileName.Data());
	}
	files.close();
	return dstTree;
}

void idDST(TString in = "/eos/nica/mpd/users/zawislak/pid/11gev/11gev.txt", Int_t nOfFiles=2 ) {
        
    TStopwatch timer;    timer.Start();

    TChain *dstTree = readAllRootFiles(in, nOfFiles);
   
    // Activate branches
    MpdEvent *event = nullptr;
    dstTree->SetBranchAddress("MPDEvent.", &event);
    TClonesArray *fMCTracks = nullptr;
    dstTree->SetBranchAddress("MCTrack", &fMCTracks);

    Int_t events = dstTree->GetEntries();
    cout << " Number of events in DST file = " << events << endl;
	
    //Event loop
    for (Int_t i = 0; i < events; i++) {
        dstTree->GetEntry(i);
        for (Int_t iTrack = 0; iTrack < event->GetGlobalTracks()->GetEntriesFast(); iTrack++) {
            MpdTrack* track = (MpdTrack*) event->GetGlobalTracks()->UncheckedAt(iTrack);
            auto mctrack = (FairMCTrack*) fMCTracks->UncheckedAt(track->GetID());
	        auto pdgID =  mctrack->GetPdgCode();
			Int_t CH=1;
			if(track->GetCharge()==0 ) continue;
			if(track->GetCharge()<0 ) CH=-1;
			if( track->GetTofFlag() != 6) continue;

			Float_t pt=track->GetPt();
			Float_t pz=track->GetPz();
			Float_t dedx=track->GetdEdXTPC();
			Float_t m2=track->GetTofMass2();
			Float_t beta=track->GetTofBeta();
			Float_t ptotal = sqrt(pz*pz+pt*pt);
			
		//Identification

			Int_t myID = identify( ptotal, dedx, m2, beta, CH);
	
		/* See mpddata/MpdTrack.h for more methods */

        } // track loop
    } // event loop

    timer.Print();
    exit(0);
}
