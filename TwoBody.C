#define M 1019.445
#define m 493.667
void TwoBody(){
	auto canv=new TCanvas("Two Body Decay","Two Body Decay",1000,800);
	
	canv->Divide(1,3);
	canv->cd(1);
	gPad->SetTitle("energy");

	//TH1F* mother_rpd;
	//TH1F* mother_pT;
	TH1F* par1_E,par1_pT,par1_pp;
	TH1F* par2_E,par2_pT,par2_pp;
	TH1F* mother_E,mother_pT,mother_pp;
	//mother_rpd=new TH1F("pseudorapidity","pseudorapidity",100,-7.,7.);
	//mother_pT=new TH1F("transverse momentum","transverse momentum",100,0.,5000.);
	
	par1_E=new TH1F("particle_1","particle 1",200,350.,1400.);
	par2_E=new TH1F("particle_2","particle 2",200,350.,1400.);
	mother_E=new TH1F("mother","mother",200,350.,1400.);
	TRandom3* rndgen=new TRandom3();
	float ps_rapid,pT,E1,E2,pT1,pT2,pp1,pp2;
	
	for(int i=0;i<5000;i++){
		
		ps_rapid=rndgen->Gaus(0.,3);
		pT=rndgen->Exp(2000.);
		float theta=2*atan(exp(-2*ps_rapid));
		float p=pT*tan(theta);
		float E=sqrt(p*p+M*M);
		float beta=p/E;
		float gamma=1./sqrt(1-beta*beta);
		E1=gamma*(M/2+beta*.5*sqrt(M*M-4*m*m));
		E2=gamma*(M/2-beta*.5*sqrt(M*M-4*m*m));
		par1_E->Fill(E1);
		par2_E->Fill(E2);
		mother_E->Fill(E);
	}
	
	par1_E->SetLineColor(kBlue);
	par2_E->SetLineColor(kBlue+3);
	mother_E->SetLineColor(kRed);
	
	THStack* En=new THStack("energy","energy of mother and daughter particles");
	En->Add(par1_E);
	En->Add(par2_E);
	En->Add(mother_E);
	En->Draw();
	En->GetXaxis()->SetTitle("energy (in MeV)");
	En->GetYaxis()->SetTitle("cross-section");
	gPad->Modified();
}
