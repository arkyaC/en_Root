#define M 1.019445
#define m 0.000510999
#define m3 .23033
#define pi 3.141592653589
#define n_bins 200
void resonance_correct(){
	auto canv1=new TCanvas("inv_mass","inv_mass",900,900);
	TH1F* inv_mass=new TH1F("inv_mass","inv_mass",n_bins,0.,2.5);
	TRandom3* rndgen=new TRandom3();
	float ps_rapid,pT;
	float** p_1;
	float** p_2;
	p_1=new float*[50];
	p_2=new float*[50];
	for(int i=0;i<50;i++){
		p_1[i]=new float[4];
		p_2[i]=new float[4];
	}
	for(int i=0; i<100; i++){
		for(int j=0; j<50; j++){
		ps_rapid=rndgen->Gaus(0.,3);
		pT=0.2+rndgen->Exp(0.5);
		//inv_mass->Fill(pT);
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
	    px1=-px1;
	    py1=-py1;
	    pz1=-pz1;
	    float px2_lab=(gamma-1)*((px1*beta_x*beta_x)/(beta*beta) + (py1*beta_x*beta_y)/(beta*beta) + (pz1*beta_x*beta_z)/(beta*beta))  +  px1  +  gamma*beta_x*(M/2);
	    float py2_lab=(gamma-1)*((px1*beta_y*beta_x)/(beta*beta) + (py1*beta_y*beta_y)/(beta*beta) + (pz1*beta_y*beta_z)/(beta*beta))  +  py1  +  gamma*beta_y*(M/2);
	    float pz2_lab=(gamma-1)*((px1*beta_z*beta_x)/(beta*beta) + (py1*beta_z*beta_y)/(beta*beta) + (pz1*beta_z*beta_z)/(beta*beta))  +  pz1  +  gamma*beta_z*(M/2);
	    float E2_lab=(gamma*(M/2))  +  gamma*((px1*beta_x) + (py1*beta_y) + (pz1*beta_z));
	    p_1[j][0]=px1_lab;
		p_1[j][1]=py1_lab;
		p_1[j][2]=pz1_lab;
		p_1[j][3]=E1_lab;
		p_2[j][0]=px2_lab;
		p_2[j][1]=py2_lab;
		p_2[j][2]=pz2_lab;
		p_2[j][3]=E2_lab;
		cout<<px1_lab<<" "<<py1_lab<<" "<<pz1_lab<<endl;
		cout<<px2_lab<<" "<<py2_lab<<" "<<pz2_lab<<endl;
		cout<<px_mother<<" "<<py_mother<<" "<<pz_mother<<endl;
		}
		for(int j=0; j<50; j++){
			for(int k=0; k<50; k++){
				float px=p_1[j][0] + p_2[k][0];
				float py=p_1[j][1] + p_2[k][1];
				float pz=p_1[j][2] + p_2[k][2];
				float E=p_1[j][3] + p_2[k][3];
				cout<<px<<"  "<<py<<"  "<<pz<<"  "<<E<<endl<<endl;
				inv_mass->Fill(sqrt((E*E) - (px*px) - (py*py) - (pz*pz)));
				//cout<<sqrt((E*E) - (px*px) - (py*py) - (pz*pz))<<endl;
			}	
		}
	}
	inv_mass->Draw();
	canv1->Modified();
}