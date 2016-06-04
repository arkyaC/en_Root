void TwoBody(){
	auto canv=new TCanvas("Two Body Decay","Mother Particle",1000,800);
	canv->Divide(2,1);
	canv->cd(1);
	TH1F* mother_rpd;
	TH1F* mother_p;
	TH1F* par1;
	TH1F* par2;
	mother_rpd=new TH1F("pseudorapidity","pseudorapidity",100,-7.,7.);
	mother_p=new TH1F("transverse momentum","transverse momentum",100,0.,5000.);
	TRandom3 rndgen;
	float ps_rapid,pT;
	for(int i=0;i<1000;i++){
		ps_rapid=rndgen.Gaus(0.,3);
		pT=rndgen.Exp(2000.);
		mother_rpd->Fill(ps_rapid);
		mother_p->Fill(pT);
	}
	mother_p->DrawClone();
	canv->cd(2);
	mother_rpd->DrawClone();
}
