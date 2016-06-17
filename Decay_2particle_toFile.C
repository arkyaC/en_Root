
#define M 1.019445
#define m .493667
#define pi 3.141592653589

void Decay_2particle_toFile(){
	TFile* ofile=TFile::Open("decay_data_phi.root","RECREATE");
	TNtuple* decay_data=new TNtuple("decay_data","phi decay data","px:py:pz:En");
	TNtuple* number_of_phi=new TNtuple("number_of_phi","number of phi particles in each run","event_num:number_of_phi");
	TRandom3* rndgen=new TRandom3();
	float ps_rapid,pT;
	
	for(int i=0;i<1000;i++){
		//if(i%1000==0)cout<<i<<endl;
		int N=(int)(rndgen->Gaus(250.,5.));//number of phi particles
		if(i%1000==0) cout<<N<<endl;
		number_of_phi->Fill(i+1,N);

			
		for(int j=1;j<=N;j++){
			ps_rapid=rndgen->Gaus(0.,3);//pseudorapidity (gaussian) of mother
			pT=.2+(rndgen->Exp(.5));//transverse momentum of mother	
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
			
		}
	}
	decay_data->Write();
	number_of_phi->Write();
	delete ofile;
	delete rndgen;
}
