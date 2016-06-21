#define M 1019.445
#define m 493.667
#define m3 230.33
#define pi 3.141592653589
#define n_bins 200
#define n 100000
/*void init(TCanvas* canv){
	canv->Divide(1,3);
	canv->cd(1);
	gPad->SetTitle("particle 1");
	canv->cd(2);
	gPad->SetTitle("particle 2");
	canv->cd(3);
	gPad->SetTitle("mother");
}*/
void Threebodyviatwobody(){

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
	TRandom3* rndgen=new TRandom3();
	float ps_rapid,pT;
	float x[n],y[n],z[n];
	
	for(int i=0;i<n;i++){
		ps_rapid=rndgen->Gaus(0.,3);
		pT=200+rndgen->Exp(500);
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
		float p1=sqrt((px1_lab*px1_lab) + (py1_lab*py1_lab) + (pz1_lab*pz1_lab));
		float E1=sqrt((p1*p1) + (m*m));
		float theta_p3cm=acos(rndgen->Uniform(-1,1));
		float phi_p3cm=2*pi*rndgen->Uniform(0,1);
		float p3_cm=0.5*sqrt(m*m - 4*m3*m3);
		float p3x=p3_cm*sin(theta_p3cm)*cos(phi_p3cm);
		float p3y=p3_cm*sin(theta_p3cm)*sin(phi_p3cm);
		float p3z=p3_cm*cos(theta_p3cm);
		float beta3=p1/E1;
		float beta3_x=px1_lab/E1;
		float beta3_y=py1_lab/E1;
		float beta3_z=pz1_lab/E1;
		float gamma3=1./sqrt(1-(beta3*beta3));
		float px3_lab=(gamma3-1)*((p3x*beta3_x*beta3_x)/(beta3*beta3) + (p3y*beta3_x*beta3_y)/(beta3*beta3) + (p3z*beta3_x*beta3_z)/(beta3*beta3))  +  p3x  +  gamma3*beta3_x*(m/2);
		float py3_lab=(gamma3-1)*((p3x*beta3_y*beta3_x)/(beta3*beta3) + (p3y*beta3_y*beta3_y)/(beta3*beta3) + (p3z*beta3_y*beta3_z)/(beta3*beta3))  +  p3y  +  gamma3*beta3_y*(m/2);
		float pz3_lab=(gamma3-1)*((p3x*beta3_z*beta3_x)/(beta3*beta3) + (p3y*beta3_z*beta3_y)/(beta3*beta3) + (p3z*beta3_z*beta3_z)/(beta3*beta3))  +  p3z  +  gamma3*beta3_z*(m/2);
		float E3_lab=(gamma3*(m/2))  +  gamma3*((p3x*beta3_x) + (p3y*beta3_y) + (p3z*beta3_z));
		p3x=-p3x;
		p3y=-p3y;
		p3z=-p3z;
		float px4_lab=(gamma3-1)*((p3x*beta3_x*beta3_x)/(beta3*beta3) + (p3y*beta3_x*beta3_y)/(beta3*beta3) + (p3z*beta3_x*beta3_z)/(beta3*beta3))  +  p3x  +  gamma3*beta3_x*(m/2);
		float py4_lab=(gamma3-1)*((p3x*beta3_y*beta3_x)/(beta3*beta3) + (p3y*beta3_y*beta3_y)/(beta3*beta3) + (p3z*beta3_y*beta3_z)/(beta3*beta3))  +  p3y  +  gamma3*beta3_y*(m/2);
		float pz4_lab=(gamma3-1)*((p3x*beta3_z*beta3_x)/(beta3*beta3) + (p3y*beta3_z*beta3_y)/(beta3*beta3) + (p3z*beta3_z*beta3_z)/(beta3*beta3))  +  p3z  +  gamma3*beta3_z*(m/2);
		float E4_lab=(gamma3*(m/2))  +  gamma3*((p3x*beta3_x) + (p3y*beta3_y) + (p3z*beta3_z));
		float px_sum=px2_lab + px3_lab + px4_lab;
		float py_sum=py2_lab + py3_lab + py4_lab;
		float pz_sum=pz2_lab + pz3_lab + pz4_lab;
		float p23x=px2_lab + px3_lab;
		float p23y=py2_lab + py3_lab;
		float p23z=pz2_lab + pz3_lab;
		float E23=E2_lab + E3_lab;
		float p24x=px4_lab + px2_lab;
		float p24y=py4_lab + py2_lab;
		float p24z=pz4_lab + pz2_lab;
		float E24=E2_lab + E4_lab;
		x[i]=(E23)*(E23) - (p23x)*(p23x) - (p23y)*(p23y) - (p23z)*(p23z);
		y[i]=(E24)*(E24) - (p24x)*(p24x) - (p24y)*(p24y) - (p24z)*(p24z);
		z[i]=((E3_lab+E4_lab)*(E3_lab+E4_lab)) - ((px3_lab+px4_lab)*(px3_lab+px4_lab)) - ((py3_lab+py4_lab)*(py3_lab+py4_lab)) - ((pz3_lab+pz4_lab)*(pz3_lab+pz4_lab));
		//cout<<px_sum<<"="<<px_mother<<"   ,   "<<py_sum<<"="<<py_mother<<"   ,   "<<pz_sum<<"="<<pz_mother<<"   ,   "<<E1_lab+E2_lab+E3_lab<<endl<<endl;
		if(i%10000==0){
			cout<<i<<endl;
		}
	}
	TGraph* dalitz1=new TGraph(n,x,y);
	TGraph* dalitz2=new TGraph(n,y,z);
	TGraph* dalitz3=new TGraph(n,x,z);
	TCanvas* can=new TCanvas("dalitz","dalitz",1500,700);
	can->Divide(3,1);
	can->cd(1);
	dalitz1->Draw("A*");
	can->cd(2);
	dalitz2->Draw("A*");
	can->cd(3);
	dalitz3->Draw("A*");
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
