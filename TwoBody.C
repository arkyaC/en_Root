#define M 1019.445
#define m 493.667
#define pi 3.141592653589
#include <cmath>
void init(TCanvas* canv){
	canv->Divide(1,3);
	canv->cd(1);
	gPad->SetTitle("particle 1");
	canv->cd(2);
	gPad->SetTitle("particle 2");
	canv->cd(3);
	gPad->SetTitle("mother");
}
void TwoBody(){
	TH1F* par1_E,par1_px,par1_py,par1_pz;
	TH1F* par2_E,par2_px,par2_py,par2_pz;
	TH1F* mother_E,mother_px,mother_py,mother_pz;

	auto canv1=new TCanvas("Two Body Decay","Energy",900,900);
	auto canv2=new TCanvas("Two Body Decay","Px",900,900);
	auto canv3=new TCanvas("Two Body Decay","Py",900,900);
	auto canv4=new TCanvas("Two Body Decay","Pz",900,900);
//initialising canvases
	init(canv1);
	init(canv2);
	init(canv3);
	init(canv4);
	//TH1F* mother_rpd=new TH1F("pseudorapidity","pseudorapidity",100,-7.,7.);
	//TH1F* mother_pT=new TH1F("transverse momentum","transverse momentum",100,0.,5000.);
	
	par1_E=new TH1F("particle_1_E","particle 1",200,350.,1300.);
	par2_E=new TH1F("particle_2_E","particle 2",200,350.,1300.);
	mother_E=new TH1F("mother_E","mother particle",200,350.,1300.);
	par1_px=new TH1F("particle_1_pX","particle 1",200,350.,1300.);
	par2_px=new TH1F("particle_2_pX","particle 2",200,350.,1300.);
	mother_px=new TH1F("mother_pX","mother particle",200,350.,1300.);
	par1_py=new TH1F("particle_1_pY","particle 1",200,350.,1300.);
	par2_py=new TH1F("particle_2_pY","particle 2",200,350.,1300.);
	mother_py=new TH1F("mother_pY","mother particle",200,350.,1300.);
	par1_pZ=new TH1F("particle_1_pZ","particle 1",200,350.,1300.);
	par2_pZ=new TH1F("particle_2_pZ","particle 2",200,350.,1300.);
	mother_pZ=new TH1F("mother_pZ","mother particle",200,350.,1300.);
	
	TRandom3* rndgen=new TRandom3();
	float ps_rapid,pT;
	
	for(int i=0;i<5000;i++){
		ps_rapid=rndgen->Gaus(0.,3);
		pT=100*rndgen->Exp(2000.);
		float theta_mother=2*atan(exp(-1*ps_rapid));
		float phi_mother=2*pi*rndgen->Uniform(0,1);
		float p_mother=pT*tan(theta_mother);
		float E=sqrt(p_mother*p_mother+M*M);
		float beta=p_mother/E;
		float gamma=1./sqrt(1-beta*beta);
		float beta_x=beta*sin(theta_mother)*cos(phi_mother);
		float beta_y=beta*sin(theta_mother)*sin(phi_mother);
		float beta_z=beta*cos(theta_mother);

		float theta=acos(rndgen->Uniform(-1,1));
		float phi=2*pi*rndgen->Uniform(0,1);
		float p=0.5*sqrt(M*M-4m*m);
		float px1=p*sin(theta)*cos(phi);
		float py1=p*sin(theta)*sin(phi);
		float pz1=p*cos(theta);
		par1_px->Fill((gamma-1)*(px1*beta_x*beta_x/(beta*beta) + py1*beta_y*beta_x/(beta*beta) + pz1*beta_z*beta_x/(beta*beta)) + px1 + gamma*beta_x*M/2);
		par1_py->Fill((gamma-1)*(px1*beta_x*beta_y/(beta*beta) + py1*beta_y*beta_y/(beta*beta) + pz1*beta_z*beta_y/(beta*beta)) + px1 + gamma*beta_y*M/2);
		par1_pz->Fill((gamma-1)*(px1*beta_x*beta_z/(beta*beta) + py1*beta_y*beta_z/(beta*beta) + pz1*beta_z*beta_z/(beta*beta)) + px1 + gamma*beta_z*M/2);
		par1_E->Fill(gamma*(M/2+beta_x*px1+beta_y*py1+beta_z*pz1));
		px1=-px1;
		py1=-py1;
		pz1=-pz1;
		par2_px->Fill((gamma-1)*(px1*beta_x*beta_x/(beta*beta)+py1*beta_y*beta_x/(beta*beta)+pz1*beta_z*beta_x/(beta*beta))+px1-gamma*beta_x*M/2);
		par2_py->Fill((gamma-1)*(px1*beta_x*beta_y/(beta*beta) + py1*beta_y*beta_y/(beta*beta) + pz1*beta_z*beta_y/(beta*beta)) + px1 + gamma*beta_y*M/2);
		par2_pz->Fill((gamma-1)*(px1*beta_x*beta_z/(beta*beta) + py1*beta_y*beta_z/(beta*beta) + pz1*beta_z*beta_z/(beta*beta)) + px1 + gamma*beta_z*M/2);
		par2_E->Fill(gamma*(M/2+beta_x*px1+beta_y*py1+beta_z*pz1));
		mother_E->Fill(E);
	}

}
