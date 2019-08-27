#ifndef _identify_h_
#define _identify_h_

/*
 * 0. 		dedx
 * 1. 		tof
 * 2.-7. 	tpctof
 * 8. 		beta
 */

const Int_t Ntplt=18;
const Int_t Npar=4;
const Int_t NPBINS=6;
const Float_t bias=0.0006;
const Float_t PTB[7] = {0, 0.2, 0.5, 0.9, 1.5, 2.2, 3.};
static const string fnam[Npar]={"tpions.root", "tkaons.root", "tprotons.root", "telectrons.root"};
static const Int_t pid[Npar]={211, 321, 2212, 11};
static TFile* files[Npar] = {
	new TFile(fnam[0].c_str()),
    new TFile(fnam[1].c_str()),
    new TFile(fnam[2].c_str()),
    new TFile(fnam[3].c_str())};

static const string names[Ntplt]={"dedxp","betap","tpctofp_0","tpctofp_1","tpctofp_2","tpctofp_3","tpctofp_4","tpctofp_5",
			   "tofp","dedxn","betan","tpctofn_0","tpctofn_1","tpctofn_2","tpctofn_3","tpctofn_4","tpctofn_5",
					  "tofn"};
					  
TH2F* GT(Int_t part, Int_t name){
	return (TH2F*) files[part]->Get(names[name].c_str());
}
//All templates
const TH2F* tplts[Npar][Ntplt] = 
{ 
  { GT(0,0), GT(0,1), GT(0,2), GT(0,3), GT(0,4), GT(0,5), GT(0,6), GT(0,7), GT(0,8),
   GT(0,9), GT(0,10), GT(0,11), GT(0,12), GT(0,13), GT(0,14), GT(0,15), GT(0,16), GT(0,17)},
   
  { GT(1,0), GT(1,1), GT(1,2), GT(1,3), GT(1,4), GT(1,5), GT(1,6), GT(1,7), GT(1,8),
   GT(1,9), GT(1,10), GT(1,11), GT(1,12), GT(1,13), GT(1,14), GT(1,15), GT(1,16), GT(1,17)},
   
  { GT(2,0), GT(2,1), GT(2,2), GT(2,3), GT(2,4), GT(2,5), GT(2,6), GT(2,7), GT(2,8),
   GT(2,9), GT(2,10), GT(2,11), GT(2,12), GT(2,13), GT(2,14), GT(2,15), GT(2,16), GT(2,17)},
   
  { GT(3,0), GT(3,1), GT(3,2), GT(3,3), GT(3,4), GT(3,5), GT(3,6), GT(3,7), GT(3,8),
   GT(3,9), GT(3,10), GT(3,11), GT(3,12), GT(3,13), GT(3,14), GT(3,15), GT(3,16), GT(3,17)},
};

Int_t checkTemplate(Int_t id, Int_t &charge, Float_t x, Float_t y, Int_t tID ){
	/*
	 * return -1	no bin found (should not appear)
	 * return 0		should never appear
	 * return 1 	bin found, possibly identified	
	 * return 2		bin found, rather not a seeked particle
	 */
	Int_t ch=0;
	if(charge<0) ch=Ntplt/2;	//IF ORDER OF HISTOGRAMS CHANGES THIS MUST CHANGE!!!!!

	//searching appropriate bin in a template
	
	//X bin
	Int_t xBin=-1;
	for( Int_t j=0; j<tplts[id][tID+ch]->GetNbinsX(); ++j){
		Float_t x0 = tplts[id][tID+ch]->GetXaxis()->GetBinCenter(j) - tplts[id][tID+ch]->GetXaxis()->GetBinWidth(j)/2;
		Float_t x1 = tplts[id][tID+ch]->GetXaxis()->GetBinCenter(j) + tplts[id][tID+ch]->GetXaxis()->GetBinWidth(j)/2;
		if( x0<=x && x<x1 ){
			xBin=j;
			break;
		}
	}
	
	//Y bin
	Int_t yBin=-1;
	for( Int_t j=0; j<tplts[id][tID+ch]->GetNbinsY(); ++j){
		Float_t y0 = tplts[id][tID+ch]->GetYaxis()->GetBinCenter(j) - tplts[id][tID+ch]->GetYaxis()->GetBinWidth(j)/2;
		Float_t y1 = tplts[id][tID+ch]->GetYaxis()->GetBinCenter(j) + tplts[id][tID+ch]->GetYaxis()->GetBinWidth(j)/2;
		if( y0<=y && y<y1 ){
			yBin=j;
			break;
		}
	}
	if( xBin==-1 || yBin==-1 ) return -1;	//No bin found. Beyond histograms
	if( tplts[id][tID+ch]->GetBinContent(xBin, yBin)!=0 ) return 1;
	if( tplts[id][tID+ch]->GetBinContent(xBin, yBin)==0 ) return 2;
	return 0;
}

Int_t identify( Float_t _ptotal, Float_t _dedx, Float_t _m2, Float_t _beta, Int_t _charge) {
	Int_t c[4][4]={0};
	Int_t FLAG=1;
    Int_t &_pbin;
	//  Search for appropriate bin in p_total
	for( Int_t i1=0; i1<NPBINS; ++i1){
		if( _ptotal > 3 ){ 
			_pbin=5;
			break;
		}
		if( PTB[i1]<=_ptotal && _ptotal<PTB[i1+1]){
			_pbin=i1;
			break;
		}
	}
	
	for( Int_t i=0; i<Npar; ++i){
		//dedx
		if(checkTemplate(i, _charge, log10(_ptotal), _dedx, 0 )==FLAG) c[i][0]+=10;
		//tof
		if(checkTemplate(i, _charge, _ptotal, _m2, 1 )==FLAG) c[i][1]+=1;
		//tpctof
		if(checkTemplate(i, _charge, _dedx, _m2, _pbin+2  )==FLAG) c[i][2]+=1000;
		//beta
		if(checkTemplate(i, _charge, _ptotal, _beta, 8 )==FLAG) c[i][3]+=100;
	}
	//c []-particle []-method
	//cout<<"\t"<<c[0][0]+c[0][1]+c[0][2]+c[0][3]<<"\t"<<c[1][0]+c[1][1]+c[1][2]+c[1][3]<<"\t"
	//	<<c[2][0]+c[2][1]+c[2][2]+c[2][3]<<"\t"<<c[3][0]+c[3][1]+c[3][2]+c[3][3]<<endl;

	//Criterium 0: method A identifies only particle B 
	Int_t record[2]={0};
	Int_t sum[4]={0};
	Bool_t ambig=false;
	for( Int_t i=0; i<4; ++i){
		//iterate through particles and sum up the id points
		Int_t s=c[i][0]+c[i][1]+c[i][2]+c[i][3];
		sum[i]=s;
		if(s==record[0])	ambig=true;
		if(s>record[0]){
			record[0]=s;
			record[1]=i;
			ambig=false;
		}
		
	}
	if(record[0]!=0){
		if(!ambig){
			 return pid[record[1]];
		}else{
			//score is the same, so let's guess (or maybe not)
		//	if(sum[0]!=0) return pid[0]; //majority of tracks are pions
		//	if(sum[2]!=0) return pid[2]; //then protons
		//	if(sum[1]!=0) return pid[1]; //then kaons
		//	if(sum[3]!=0) return pid[3]; //it may be electron
		}
	}
	//cout<<c[0][0]<<"\t"<<c[0][1]<<"\t"<<c[0][2]<<"\t"<<c[0][3]<<endl<<c[1][0]<<"\t"<<c[1][1]<<"\t"<<c[1][2]<<"\t"<<c[1][3]<<endl
	//	<<c[2][0]<<"\t"<<c[2][1]<<"\t"<<c[2][2]<<"\t"<<c[2][3]<<endl<<c[3][0]<<"\t"<<c[3][1]<<"\t"<<c[3][2]<<"\t"<<c[3][3]<<endl<<endl<<endl;
	
	
	return -1;	//pile-up
}
 
 
#endif
