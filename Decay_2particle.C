#define M 1019.445
#define m 493.667
#define pi 3.141592653589
#include <cmath> 
void TwoBody(){
	TH1F* par1_E,par1_px,par1_py,par1_pz;
	TH1F* par2_E,par2_px,par2_py,par2_pz;
	TH1F* mother_E,mother_px,mother_py,nother_pz;

	auto canv1=new TCanvas("p_x","p_x",900,900);
	auto canv2=new TCanvas("p_y","p_y",900,900);
	auto canv3=new TCanvas("p_z","p_z",900,900);
	auto canv4=new TCanvas("E","E",900,900);
	canv1->Divide(1,3);
	canv2->Divide(1,3);
	canv3->Divide(1,3);
	canv4->Divide(1,3);
	canv1->cd(1);
	gPad->SetTitle("Particle 1");
	canv1->cd(2);
	gPad->SetTitle("Particle 2");
	canv1->cd(3);
	gPad->SetTitle("Mother");
	canv2->cd(1);
	gPad->SetTitle("Particle 1");
	canv2->cd(2);
	gPad->SetTitle("Particle 2");
	canv2->cd(3);
	gPad->SetTitle("Mother");
	canv3->cd(1);
	gPad->SetTitle("Particle 1");
	canv3->cd(2);
	gPad->SetTitle("Particle 2");
	canv3->cd(3);
	gPad->SetTitle("Mother");
	canv4->cd(1);
	gPad->SetTitle("Particle 1");
	canv4->cd(2);
	gPad->SetTitle("Particle 2");
	canv4->cd(3);
	gPad->SetTitle("Mother");
	
	par1_E=new TH1F("particle_1_E","particle 1",200,350.,1300.);
	par2_E=new TH1F("particle_2_E","particle 2",200,350.,1300.);
	mother_E=new TH1F("mother_E","mother particle",200,350.,1300.);
	par1_px=new TH1F("particle_1_px","particle 1",200,350.,1300.);
	par2_px=new TH1F("particle_2_px","particle 2",200,350.,1300.);
	mother_px=new TH1F("mother_px","mother particle",200,350.,1300.);
	par1_py=new TH1F("particle_1_py","particle 1",200,350.,1300.);
	par2_py=new TH1F("particle_2_py","particle 2",200,350.,1300.);
	mother_py=new TH1F("mother_py","mother particle",200,350.,1300.);
	par1_pz=new TH1F("particle_1_pz","particle 1",200,350.,1300.);
	par2_pz=new TH1F("particle_2_pz","particle 2",200,350.,1300.);
	mother_pz=new TH1F("mother_pz","mother particle",200,350.,1300.);
	
	TRandom3* rndgen=new TRandom3();
	float ps_rapid,pT;
	
	for(int i=0;i<5000;i++){
		
		ps_rapid=rndgen->Gaus(0.,3);
		pT=100*rndgen->Exp(2000.);
		float theta_mother=2*atan(exp(-1*ps_rapid));
		float p_mother=pT*tan(theta_mother);
		float E=sqrt(p_mother*p_mother+M*M);
		float beta=p_mother/E;
		float p=0.5*sqrt(M*M - 4*m*m);
		float gamma=1./sqrt(1-beta*beta);
		float theta=acos(rndgen->Uniform(-1,1));
		float phi=2*pi*rndgen->Uniform(0,1);
		float phi_mother=2*pi*rndgen->Uniform(0,1);
	    float beta_x = beta*sin(theta_mother)*cos(phi_mother);
	    float beta_y = beta*sin(theta_mother)*sin(phi_mother);
	    float beta_z = beta*cos(theta_mother);
	    float px1=p*sin(theta)*cos(phi);
	    float py1=p*sin(theta)*sin(phi);
	    float pz1=p*cos(theta);
	    par1_px->Fill((gamma-1)*((px1*beta_x*beta_x)/(beta*beta) + (py1*beta_x*beta_y)/(beta*beta) + (pz1*beta_x*beta_z)/(beta*beta))  +  px1  +  gamma*beta_x*(M/2));
	    par1_py->Fill((gamma-1)*((px1*beta_y*beta_x)/(beta*beta) + (py1*beta_y*beta_y)/(beta*beta) + (pz1*beta_y*beta_z)/(beta*beta))  +  py1  +  gamma*beta_y*(M/2));
	    par1_pz->Fill((gamma-1)*((px1*beta_z*beta_x)/(beta*beta) + (py1*beta_z*beta_y)/(beta*beta) + (pz1*beta_z*beta_z)/(beta*beta))  +  pz1  +  gamma*beta_z*(M/2));
	    par1_E->Fill((gamma*(M/2))  +  gamma*((px1*beta_x) + (py1*beta_y) + (pz1*beta_z)));
	    px1=-px1;
	    py1=-py1;
	    pz1=-pz1;
		par2_px->Fill((gamma-1)*((px1*beta_x*beta_x)/(beta*beta) + (py1*beta_x*beta_y)/(beta*beta) + (pz1*beta_x*beta_z)/(beta*beta))  +  px1  +  gamma*beta_x*(M/2));
	    par2_py->Fill((gamma-1)*((px1*beta_y*beta_x)/(beta*beta) + (py1*beta_y*beta_y)/(beta*beta) + (pz1*beta_y*beta_z)/(beta*beta))  +  py1  +  gamma*beta_y*(M/2));
	    par2_pz->Fill((gamma-1)*((px1*beta_z*beta_x)/(beta*beta) + (py1*beta_z*beta_y)/(beta*beta) + (pz1*beta_z*beta_z)/(beta*beta))  +  pz1  +  gamma*beta_z*(M/2));
	    par2_E->Fill((gamma*(M/2))  +  gamma*((px1*beta_x) + (py1*beta_y) + (pz1*beta_z)));
		mother_E->Fill(E);
	}

}
