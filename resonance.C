#define M 1.019445
#define m .493667
#define m3 .23033
#define pi 3.141592653589
#define n_bins 300
/*void init(TCanvas* canv){
	canv->Divide(1,3);
	canv->cd(1);
	gPad->SetTitle("particle 1");
	canv->cd(2);
	gPad->SetTitle("particle 2");
	canv->cd(3);
	gPad->SetTitle("mother");
}
*/
void resonance(){

	/*auto canv1=new TCanvas("p_x","p_x",900,900);
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
	*/
	auto canv1=new TCanvas("inv_mass","inv_mass",900,900);
	TH1F* inv_mass=new TH1F("inv_mass","inv_mass",n_bins,0,10);
	TRandom3* rndgen=new TRandom3();
	float ps_rapid,pT;
	float** p_1;
	float** p_2;
	p_1=new float*[1000];
	p_2=new float*[1000];
	for(int i=0;i<1000;i++){
		p_1[i]=new float[4];
		p_2[i]=new float[4];
	}
	for(int i=0;i<1000;i++){
		ps_rapid=rndgen->Gaus(0.,3);
		pT=0.1*rndgen->Exp(2.);
		float theta_mother=2*atan(exp(-1*ps_rapid));
		float p_mother=pT/sin(theta_mother);
		float E=sqrt(p_mother*p_mother+M*M);

		float beta=p_mother/E;
		float p=0.5*sqrt(M*M - 4*m*m);
		float gamma=1./sqrt(1-beta*beta);
		float theta=acos(rndgen->Uniform(-1,1));
		float phi=2*pi*rndgen->Uniform(0,1);
		float phi_mother=2*pi*rndgen->Uniform(0,1);
		float px_mother=p_mother*sin(theta_mother)*cos(phi_mother);
		float py_mother=p_mother*sin(theta_mother)*sin(phi_mother);
		float pz_mother=p_mother*cos(theta_mother);

		//mother_px->Fill(p_mother*sin(theta_mother)*cos(phi_mother));
		//mother_py->Fill(p_mother*sin(theta_mother)*sin(phi_mother));
		//mother_pz->Fill(p_mother*cos(theta_mother));

	    float beta_x = beta*sin(theta_mother)*cos(phi_mother);
	    float beta_y = beta*sin(theta_mother)*sin(phi_mother);
	    float beta_z = beta*cos(theta_mother);
	    float px1=p*sin(theta)*cos(phi);
	    float py1=p*sin(theta)*sin(phi);
	    float pz1=p*cos(theta);
	    float px1_lab=(gamma-1)*((px1*beta_x*beta_x)/(beta*beta) + (py1*beta_x*beta_y)/(beta*beta) + (pz1*beta_x*beta_z)/(beta*beta))  +  px1  +  gamma*beta_x*(M/2);
	    float py1_lab=(gamma-1)*((px1*beta_y*beta_x)/(beta*beta) + (py1*beta_y*beta_y)/(beta*beta) + (pz1*beta_y*beta_z)/(beta*beta))  +  py1  +  gamma*beta_y*(M/2);
	    float pz1_lab=(gamma-1)*((px1*beta_z*beta_x)/(beta*beta) + (py1*beta_z*beta_y)/(beta*beta) + (pz1*beta_z*beta_z)/(beta*beta))  +  pz1  +  gamma*beta_z*(M/2);
	    float E1_lab=(gamma*(M/2))  +  gamma*((px1*beta_x) + (py1*beta_y) + (pz1*beta_z));
	    //par1_px->Fill((gamma-1)*((px1*beta_x*beta_x)/(beta*beta) + (py1*beta_x*beta_y)/(beta*beta) + (pz1*beta_x*beta_z)/(beta*beta))  +  px1  +  gamma*beta_x*(M/2));
	    //par1_py->Fill((gamma-1)*((px1*beta_y*beta_x)/(beta*beta) + (py1*beta_y*beta_y)/(beta*beta) + (pz1*beta_y*beta_z)/(beta*beta))  +  py1  +  gamma*beta_y*(M/2));
	    //par1_pz->Fill((gamma-1)*((px1*beta_z*beta_x)/(beta*beta) + (py1*beta_z*beta_y)/(beta*beta) + (pz1*beta_z*beta_z)/(beta*beta))  +  pz1  +  gamma*beta_z*(M/2));
	    //par1_E->Fill((gamma*(M/2))  +  gamma*((px1*beta_x) + (py1*beta_y) + (pz1*beta_z)));
	    px1=-px1;
	    py1=-py1;
	    pz1=-pz1;
	    float px2_lab=(gamma-1)*((px1*beta_x*beta_x)/(beta*beta) + (py1*beta_x*beta_y)/(beta*beta) + (pz1*beta_x*beta_z)/(beta*beta))  +  px1  +  gamma*beta_x*(M/2);
	    float py2_lab=(gamma-1)*((px1*beta_y*beta_x)/(beta*beta) + (py1*beta_y*beta_y)/(beta*beta) + (pz1*beta_y*beta_z)/(beta*beta))  +  py1  +  gamma*beta_y*(M/2);
	    float pz2_lab=(gamma-1)*((px1*beta_z*beta_x)/(beta*beta) + (py1*beta_z*beta_y)/(beta*beta) + (pz1*beta_z*beta_z)/(beta*beta))  +  pz1  +  gamma*beta_z*(M/2);
	    float E2_lab=(gamma*(M/2))  +  gamma*((px1*beta_x) + (py1*beta_y) + (pz1*beta_z));
		//par2_px->Fill((gamma-1)*((px1*beta_x*beta_x)/(beta*beta) + (py1*beta_x*beta_y)/(beta*beta) + (pz1*beta_x*beta_z)/(beta*beta))  +  px1  +  gamma*beta_x*(M/2));
	    //par2_py->Fill((gamma-1)*((px1*beta_y*beta_x)/(beta*beta) + (py1*beta_y*beta_y)/(beta*beta) + (pz1*beta_y*beta_z)/(beta*beta))  +  py1  +  gamma*beta_y*(M/2));
	    //par2_pz->Fill((gamma-1)*((px1*beta_z*beta_x)/(beta*beta) + (py1*beta_z*beta_y)/(beta*beta) + (pz1*beta_z*beta_z)/(beta*beta))  +  pz1  +  gamma*beta_z*(M/2));
	    //par2_E->Fill((gamma*(M/2))  +  gamma*((px1*beta_x) + (py1*beta_y) + (pz1*beta_z)));
		//mother_E->Fill(E);
		p_1[i][0]=px1_lab;
		p_1[i][1]=py1_lab;
		p_1[i][2]=pz1_lab;
		p_1[i][3]=E1_lab;
		p_2[i][0]=px2_lab;
		p_2[i][1]=py2_lab;
		p_2[i][2]=pz2_lab;
		p_2[i][3]=E2_lab;
	}
	for(int i=0; i<1000; i++){
		for(int j=0; j<5; j++){
			float px=p_1[i][0] + p_2[j][0];
			float py=p_1[i][1] + p_2[j][1];
			float pz=p_1[i][2] + p_2[j][2];
			float E=p_1[i][3] + p_2[j][3];
			cout<<px<<"  "<<py<<"  "<<pz<<"  "<<E<<endl<<endl;
			inv_mass->Fill(sqrt((E*E) - (px*px) - (py*py) - (pz*pz)));
		}	
	}
	//canv1->cd(0);
	inv_mass->DrawNormalized();
	canv1->Modified();
	canv1->update();
	/*canv1->cd(1);
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
	canv1->Modified();
	canv2->Modified();
	canv3->Modified();
	canv4->Modified();*/
}
