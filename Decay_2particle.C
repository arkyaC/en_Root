#define M 1019.445
#define m 493.667
#define pi 3.141592653589
#define n_bins 200

void init(TCanvas* canv){
	canv->Divide(1,3);
	canv->cd(1);
	gPad->SetTitle("particle 1");
	canv->cd(2);
	gPad->SetTitle("particle 2");
	canv->cd(3);
	gPad->SetTitle("mother");
}
void Decay_2particle(){

	auto canv1=new TCanvas("p_x","p_x",900,900);
	auto canv2=new TCanvas("p_y","p_y",900,900);
	auto canv3=new TCanvas("p_z","p_z",900,900);
	auto canv4=new TCanvas("E","E",900,900);
	init(canv1);
	init(canv2);
	init(canv3);
	init(canv4);
	
	TH1F* par1_E=new TH1F("particle_1_E","particle 1",n_bins,350.,1600.);
	TH1F* par2_E=new TH1F("particle_2_E","particle 2",n_bins,350.,1600.);
	TH1F* mother_E=new TH1F("mother_E","mother particle",n_bins,350.,1600.);
	TH1F* par1_px=new TH1F("particle_1_px","particle 1",n_bins,-700.,700.);
	TH1F* par2_px=new TH1F("particle_2_px","particle 2",n_bins,-700.,700.);
	TH1F* mother_px=new TH1F("mother_px","mother particle",n_bins,-800.,800.);
	TH1F* par1_py=new TH1F("particle_1_py","particle 1",n_bins,-700.,700.);
	TH1F* par2_py=new TH1F("particle_2_py","particle 2",n_bins,-700.,700.);
	TH1F* mother_py=new TH1F("mother_py","mother particle",n_bins,-800.,800.);
	TH1F* par1_pz=new TH1F("particle_1_pz","particle 1",n_bins,-700.,700.);
	TH1F* par2_pz=new TH1F("particle_2_pz","particle 2",n_bins,-700.,700.);
	TH1F* mother_pz=new TH1F("mother_pz","mother particle",n_bins,-2000.,2000.);
	
	TRandom3* rndgen=new TRandom3();
	float ps_rapid,pT;
	
	for(int i=0;i<1000000;i++){
		//if(i%100==0)cout<<i<<endl;
		ps_rapid=rndgen->Gaus(0.,3);
		pT=0.01*(rndgen->Exp(2000.));
		float theta_mother=2*atan(exp(-1*ps_rapid));
		float p_mother=pT/sin(theta_mother);
		float E=sqrt(p_mother*p_mother+M*M);

		float beta=p_mother/E;
		float p=0.5*sqrt(M*M - 4*m*m);
		float gamma=1./sqrt(1-beta*beta);
		float theta=acos(rndgen->Uniform(-1,1));
		float phi=2*pi*rndgen->Uniform(0,1);
		float phi_mother=2*pi*rndgen->Uniform(0,1);

		float px_m=p_mother*sin(theta_mother)*cos(phi_mother);
		float py_m=p_mother*sin(theta_mother)*sin(phi_mother);
		float pz_m=p_mother*cos(theta_mother);

	    float beta_x = beta*sin(theta_mother)*cos(phi_mother);
	    float beta_y = beta*sin(theta_mother)*sin(phi_mother);
	    float beta_z = beta*cos(theta_mother);

	    float px1=p*sin(theta)*cos(phi);
	    float py1=p*sin(theta)*sin(phi);
	    float pz1=p*cos(theta);

	    float Px1=(gamma-1)*((px1*beta_x*beta_x)/(beta*beta) + (py1*beta_x*beta_y)/(beta*beta) + (pz1*beta_x*beta_z)/(beta*beta))  +  px1  +  gamma*beta_x*(M/2);
	    float Py1=(gamma-1)*((px1*beta_y*beta_x)/(beta*beta) + (py1*beta_y*beta_y)/(beta*beta) + (pz1*beta_y*beta_z)/(beta*beta))  +  py1  +  gamma*beta_y*(M/2);
	    float Pz1=(gamma-1)*((px1*beta_z*beta_x)/(beta*beta) + (py1*beta_z*beta_y)/(beta*beta) + (pz1*beta_z*beta_z)/(beta*beta))  +  pz1  +  gamma*beta_z*(M/2);
	    float E1=(gamma*(M/2))  +  gamma*((px1*beta_x) + (py1*beta_y) + (pz1*beta_z));
	    px1=-px1;
	    py1=-py1;
	    pz1=-pz1;
	    float Px2=(gamma-1)*((px1*beta_x*beta_x)/(beta*beta) + (py1*beta_x*beta_y)/(beta*beta) + (pz1*beta_x*beta_z)/(beta*beta))  +  px1  +  gamma*beta_x*(M/2);
	    float Py2=(gamma-1)*((px1*beta_y*beta_x)/(beta*beta) + (py1*beta_y*beta_y)/(beta*beta) + (pz1*beta_y*beta_z)/(beta*beta))  +  py1  +  gamma*beta_y*(M/2);
	    float Pz2=(gamma-1)*((px1*beta_z*beta_x)/(beta*beta) + (py1*beta_z*beta_y)/(beta*beta) + (pz1*beta_z*beta_z)/(beta*beta))  +  pz1  +  gamma*beta_z*(M/2);
	    float E2=(gamma*(M/2))  +  gamma*((px1*beta_x) + (py1*beta_y) + (pz1*beta_z));

	    //testing momentum and energy conservation at an interval
	    if(i%10000==0) cout<<Px1+Px2-px_m<<" "<<Py1+Py2-py_m<<" "<<Pz1+Pz2-pz_m<<" "<<E1+E2-E<<endl;

	    par1_px->Fill(Px1);
	    par1_py->Fill(Py1);
	    par1_pz->Fill(Pz1);
	    par1_E->Fill(E1);
		par2_px->Fill(Px2);
	    par2_py->Fill(Py2);
	    par2_pz->Fill(Pz2);
	    par2_E->Fill(E2);
	    mother_px->Fill(px_m);
		mother_py->Fill(py_m);
		mother_pz->Fill(pz_m);
		mother_E->Fill(E);
	}
	canv1->cd(1);
	par1_px->Draw();
	canv1->cd(2);
	par2_px->Draw();
	canv1->cd(3);
	mother_px->Draw();
	canv2->cd(1);
	par1_py->Draw();
	canv2->cd(2);
	par2_py->Draw();
	canv2->cd(3);
	mother_py->Draw();
	canv3->cd(1);
	par1_pz->Draw();
	canv3->cd(2);
	par2_pz->Draw();
	canv3->cd(3);
	mother_pz->Draw();
	canv4->cd(1);
	par1_E->Draw();
	canv4->cd(2);
	par2_E->Draw();
	canv4->cd(3);
	mother_E->Draw();
	/*canv1->Modified();
	canv2->Modified();
	canv3->Modified();
	canv4->Modified();*/
}
