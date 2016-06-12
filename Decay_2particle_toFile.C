
#define M 1019.445
#define m 493.667
#define pi 3.141592653589

void Decay_2particle(){
	TFile* ofile=TFile::Open("decay_data_phi.root","RECREATE");
	TNtuple* decay_data=new TNtuple("decay_data","phi decay data","px:py:pz:En");
	TNtuple* mother_data=new TNtuple("mother_data","phi mother data","px:py:pz:En");
	TRandom3* rndgen=new TRandom3();
	float ps_rapid,pT;
	
	for(int i=0;i<10000;i++){
		//if(i%1000==0)cout<<i<<endl;
		int N=(int)(rndgen->Gaus(5,2));//number of phi particles
		if(i%1000==0) cout<<N<<endl;

		for(int j=1;j<=N;j++){
			ps_rapid=rndgen->Gaus(0.,3);//pseudorapidity (gaussian) of mother
			pT=100*(rndgen->Exp(2000.));//transverse momentum of mother
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

		    //calculating momentum and energy in COM frame
		    float px1=p*sin(theta)*cos(phi);
		    float py1=p*sin(theta)*sin(phi);
		    float pz1=p*cos(theta);

		    //converting to LAB frame
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
		    //if(i%10000==0) cout<<Px1+Px2-px_m<<" "<<Py1+Py2-py_m<<" "<<Pz1+Pz2-pz_m<<" "<<E1+E2-E<<endl;

			decay_data->Fill(Px1,Py1,Pz1,E1);
			decay_data->Fill(Px2,Py2,Pz2,E2);
			mother_data->Fill(px_m,py_m,pz_m,E);
		}
	}
	decay_data->Write();
	mother_data->Write();
	delete ofile;
	delete rndgen;
}