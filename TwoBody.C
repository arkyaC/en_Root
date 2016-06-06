#define M 1019.445
#define m 493.667
void TwoBody(){
	TH1F* par1_E,par1_pT,par1_pp;
	TH1F* par2_E,par2_pT,par2_pp;
	TH1F* mother_E,mother_pT,mother_pp;

	auto canv=new TCanvas("Two Body Decay","Two Body Decay",1000,800);
	canv->Divide(1,3);

	canv->cd(1);
	gPad->SetTitle("energy");
	canv->cd(2);
	gPad->SetTitle("p_{T}");
	canv->cd(3);
	gPad->SetTitle("p_{P}");
	canv->cd(1);
	//TH1F* mother_rpd=new TH1F("pseudorapidity","pseudorapidity",100,-7.,7.);
	//TH1F* mother_pT=new TH1F("transverse momentum","transverse momentum",100,0.,5000.);
	
	par1_E=new TH1F("particle_1_E","particle 1",200,350.,1300.);
	par2_E=new TH1F("particle_2_E","particle 2",200,350.,1300.);
	mother_E=new TH1F("mother_E","mother particle",200,350.,1300.);
	/*par1_pT=new TH1F("mother","particle 1 pT",200,350.,1400.);
	par2_pT=new TH1F("mother","particle 2 pT",200,350.,1400.);
	mother_pT=new TH1F("mother","mother pT",200,0.,4000.);
	par1_pp=new TH1F("mother","particle 1 pP",200,350.,1400.);
	par2_pp=new TH1F("mother","particle 2 pP",200,350.,1400.);
	mother_pp=new TH1F("mother","mother pP",200,350.,1400.);*/
	TRandom3* rndgen=new TRandom3();
	float ps_rapid,pT,E1,E2,pT1,pT2,pp1,pp2;
	
	for(int i=0;i<5000;i++){
		
		ps_rapid=rndgen->Gaus(0.,3);
		pT=rndgen->Exp(2000.);
		float theta=2*atan(exp(-1*ps_rapid));
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
	par2_E->SetLineColor(kBlue+10);
	mother_E->SetLineColor(kRed);
	
	THStack* En=new THStack("energy","energy of mother and daughter particles");
	En->Add(par1_E);
	En->Add(par2_E);
	En->Add(mother_E);
	En->Draw();
	En->GetXaxis()->SetTitle("energy (in MeV)");
	En->GetYaxis()->SetTitle("cross-section");
	gPad->BuildLegend();
	gPad->Modified();

}
