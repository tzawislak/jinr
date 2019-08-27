#ifndef _CreateIDTemplates_h_
#define _CreateIDTemplates_h_

void crophist(Float_t bias=0.01) {
	const Int_t Ntplt=18;
    const Int_t Npar=4;
    
    TFile* files[Npar];
    string fnam[Npar]={"pions.root", "kaons.root", "protons.root", "electrons.root"};

    files[0] = new TFile(fnam[0].c_str());
    files[1] = new TFile(fnam[1].c_str());
    files[2] = new TFile(fnam[2].c_str());
    files[3] = new TFile(fnam[3].c_str());
    string names[Ntplt]={"dedxp","betap","tpctofp_0","tpctofp_1","tpctofp_2","tpctofp_3","tpctofp_4","tpctofp_5",
			   "tofp","dedxn","betan","tpctofn_0","tpctofn_1","tpctofn_2","tpctofn_3","tpctofn_4","tpctofn_5",
					  "tofn"};
    TH2F *his;

    for(Int_t t=0; t<Npar; ++t){
		string title = "t"+fnam[t];
		TFile* tf = new TFile(title.c_str(), "recreate");
		
		for(Int_t k=0; k<Ntplt; ++k){
			his = (TH2F*) files[t]->Get(names[k].c_str());
		
			Int_t nx = his->GetNbinsX();
			Int_t ny = his->GetNbinsY();
			Float_t max = his->GetMaximum();
		
			for( Int_t i=0; i<nx; ++i){
				for( Int_t j=0; j<ny; ++j){
					if( his->GetBinContent(i, j) == 0) continue;
					if( his->GetBinContent(i, j) < bias*max)//a moze procent z binu o najwiekszej l wejsc?
						his->SetBinContent(i, j, 0);
				}
			}
			tf->cd();
			his->Write();
		}
		tf->Close();
				
    }	
}



void checkCropResult(){
	const Int_t Ntplt=18;
    const Int_t Npar=4;
    
    TFile* files[Npar];
    TFile* tfiles[Npar];
    string fnam[Npar]={"pions.root", "kaons.root", "protons.root", "electrons.root"};
    
    files[0] = new TFile(fnam[0].c_str());
    files[1] = new TFile(fnam[1].c_str());
    files[2] = new TFile(fnam[2].c_str());
    files[3] = new TFile(fnam[3].c_str());
    for(Int_t i=0; i<Npar; ++i)
		fnam[i] = "t" + fnam[i];
    tfiles[0] = new TFile(fnam[0].c_str());
    tfiles[1] = new TFile(fnam[1].c_str());
    tfiles[2] = new TFile(fnam[2].c_str());
    tfiles[3] = new TFile(fnam[3].c_str());
    string names[Ntplt]={"dedxp","betap","tpctofp_0","tpctofp_1","tpctofp_2","tpctofp_3","tpctofp_4","tpctofp_5",
			   "tofp","dedxn","betan","tpctofn_0","tpctofn_1","tpctofn_2","tpctofn_3","tpctofn_4","tpctofn_5",
					  "tofn"};
    TH2F *his;
    TH2F *hhis;
    
    TCanvas* canv[4];
    for(Int_t p=0; p<4; ++p){
		canv[p] = new TCanvas();
		canv[p]->Divide( 4, 5);
	}
   
    for(Int_t t=0; t<Npar; ++t){
		for(Int_t k=0; k<Ntplt; ++k){
			his = (TH2F*) files[t]->Get(names[k].c_str());
			hhis = (TH2F*) tfiles[t]->Get(names[k].c_str());
			canv[t]->cd(k+1);
			Double_t ratio = hhis->Integral(0, hhis->GetNbinsX()-1, 0, hhis->GetNbinsY()-1 )
							/his->Integral(0, his->GetNbinsX()-1, 0, his->GetNbinsY()-1 );
			his->SetTitle(Form("%f",ratio));
			his->Draw();
			hhis->SetMarkerColor(kRed);
			hhis->Draw("same");
			
		}			
    }	
    /*
    for(Int_t i=0; i<Npar; ++i){
		files[i]->Close();
		tfiles[i]->Close();
	}
    */
}
 
 
#endif
