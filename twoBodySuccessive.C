
#define M 1019.445 //neutron mass
#define m 493.667 //proton mass
#define m1 230.33 //electron mass
#define pi 3.141592653589

void twoBodySuccessive(){

	TRandom3* rndgen=new TRandom3();
	TNtuple* dalitzPlot41=new TNtuple("dalitz","s4 versus s1","s4:s1");
	TNtuple* dalitzPlot34=new TNtuple("dalitz","s3 versus s4","s3:s4");
	TNtuple* dalitzPlot13=new TNtuple("dalitz","s1 versus s3","s1:s3");
	TCanvas* canv=new TCanvas("dalitz_plot","Dalitz plots",900,320);
	canv->Divide(3,1);

	canv->cd(1);
	gPad->SetTitle("s4 versus s1");
	canv->cd(2);
	gPad->SetTitle("s3 versus s4");
	canv->cd(3);
	gPad->SetTitle("s1 versus s3");	

	float ps_rapid,pT;
	
	for(int i=0;i<1000;i++){
		//if(i%1000==0)cout<<i<<endl;
		int N=1;//(int)(rndgen->Gaus(250.,5.));//number of phi particles
		//if(i%1000==0) cout<<N<<endl;
			
		for(int j=1;j<=N;j++){
			ps_rapid=rndgen->Gaus(0.,3);//pseudorapidity (gaussian) of mother
			pT=200+(rndgen->Exp(500));//transverse momentum of mother	
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

		    //converting to LAB frame...

		    //particle 1 in LAB
		    float Px1=(gamma-1)*((px1*beta_x*beta_x)/(beta*beta) + (py1*beta_x*beta_y)/(beta*beta) + (pz1*beta_x*beta_z)/(beta*beta))  +  px1  +  gamma*beta_x*(M/2);
		    float Py1=(gamma-1)*((px1*beta_y*beta_x)/(beta*beta) + (py1*beta_y*beta_y)/(beta*beta) + (pz1*beta_y*beta_z)/(beta*beta))  +  py1  +  gamma*beta_y*(M/2);
		    float Pz1=(gamma-1)*((px1*beta_z*beta_x)/(beta*beta) + (py1*beta_z*beta_y)/(beta*beta) + (pz1*beta_z*beta_z)/(beta*beta))  +  pz1  +  gamma*beta_z*(M/2);
		    float E1=(gamma*(M/2))  +  gamma*((px1*beta_x) + (py1*beta_y) + (pz1*beta_z));
		    px1=-px1;
		    py1=-py1;
		    pz1=-pz1;
		    
		    //particle 2 in LAB
		    float Px2=(gamma-1)*((px1*beta_x*beta_x)/(beta*beta) + (py1*beta_x*beta_y)/(beta*beta) + (pz1*beta_x*beta_z)/(beta*beta))  +  px1  +  gamma*beta_x*(M/2);
		    float Py2=(gamma-1)*((px1*beta_y*beta_x)/(beta*beta) + (py1*beta_y*beta_y)/(beta*beta) + (pz1*beta_y*beta_z)/(beta*beta))  +  py1  +  gamma*beta_y*(M/2);
		    float Pz2=(gamma-1)*((px1*beta_z*beta_x)/(beta*beta) + (py1*beta_z*beta_y)/(beta*beta) + (pz1*beta_z*beta_z)/(beta*beta))  +  pz1  +  gamma*beta_z*(M/2);
		    float E2=(gamma*(M/2))  +  gamma*((px1*beta_x) + (py1*beta_y) + (pz1*beta_z));

		    //testing momentum and energy conservation at an interval
		    //if(i%10000==0) cout<<Px1+Px2-px_m<<" "<<Py1+Py2-py_m<<" "<<Pz1+Pz2-pz_m<<" "<<E1+E2-E<<endl;

//new decay
		    beta_x=Px2/E2;
		    beta_y=Py2/E2;
		    beta_z=Pz2/E2;
		    beta=sqrt(beta_x*beta_x+beta_y*beta_y+beta_z*beta_z);
		    gamma=1./sqrt(1-beta*beta);
		    p=0.5*sqrt(m*m-4*m1*m1);
		    theta=acos(rndgen->Uniform(-1,1));
			phi=2*pi*rndgen->Uniform(0,1);

		    px1=p*sin(theta)*cos(phi);
		    py1=p*sin(theta)*sin(phi);
		    pz1=p*cos(theta);

		    float Px3=(gamma-1)*((px1*beta_x*beta_x)/(beta*beta) + (py1*beta_x*beta_y)/(beta*beta) + (pz1*beta_x*beta_z)/(beta*beta))  +  px1  +  gamma*beta_x*(m/2);
		    float Py3=(gamma-1)*((px1*beta_y*beta_x)/(beta*beta) + (py1*beta_y*beta_y)/(beta*beta) + (pz1*beta_y*beta_z)/(beta*beta))  +  py1  +  gamma*beta_y*(m/2);
		    float Pz3=(gamma-1)*((px1*beta_z*beta_x)/(beta*beta) + (py1*beta_z*beta_y)/(beta*beta) + (pz1*beta_z*beta_z)/(beta*beta))  +  pz1  +  gamma*beta_z*(m/2);
		    float E3=(gamma*(m/2))  +  gamma*((px1*beta_x) + (py1*beta_y) + (pz1*beta_z));
		    px1=-px1;
		    py1=-py1;
		    pz1=-pz1;
		    
		    //particle 2 in LAB
		    float Px4=(gamma-1)*((px1*beta_x*beta_x)/(beta*beta) + (py1*beta_x*beta_y)/(beta*beta) + (pz1*beta_x*beta_z)/(beta*beta))  +  px1  +  gamma*beta_x*(m/2);
		    float Py4=(gamma-1)*((px1*beta_y*beta_x)/(beta*beta) + (py1*beta_y*beta_y)/(beta*beta) + (pz1*beta_y*beta_z)/(beta*beta))  +  py1  +  gamma*beta_y*(m/2);
		    float Pz4=(gamma-1)*((px1*beta_z*beta_x)/(beta*beta) + (py1*beta_z*beta_y)/(beta*beta) + (pz1*beta_z*beta_z)/(beta*beta))  +  pz1  +  gamma*beta_z*(m/2);
		    float E4=(gamma*(m/2))  +  gamma*((px1*beta_x) + (py1*beta_y) + (pz1*beta_z));

		    // s4 aka s13
		    float s13=(pow(E1+E3,2.)-pow(Px1+Px3,2.)-pow(Py1+Py3,2.)-pow(Pz1+Pz3,2.));
		    // s3 aka s14
		    float s14=(pow(E1+E4,2.)-pow(Px1+Px4,2.)-pow(Py1+Py4,2.)-pow(Pz1+Pz4,2.));
		    // s1 aka s43
		    float s43=(pow(E4+E3,2.)-pow(Px4+Px3,2.)-pow(Py4+Py3,2.)-pow(Pz4+Pz3,2.));

		    cout<<Px1+Px3+Px4-px_m<<endl<<endl;

		    dalitzPlot41->Fill(s13,s43);
		    dalitzPlot34->Fill(s14,s13);
		    dalitzPlot13->Fill(s43,s14);
		}
	}

	canv->cd(1);
	dalitzPlot41->Draw("s4:s1");
	TH2F* htemp1=(TH2F*)gPad->GetPrimitive("htemp");
	htemp1->SetTitle("dalitz plot of s13 versus s34");
	gPad->Modified();gPad->Update();

	canv->cd(2);
	dalitzPlot34->Draw("s3:s4");
	TH2F* htemp2=(TH2F*)gPad->GetPrimitive("htemp");
	htemp2->SetTitle("dalitz plot of s14 versus s13");
	gPad->Modified();gPad->Update();	

	canv->cd(3);
	dalitzPlot13->Draw("s1:s3");
	TH2F* htemp3=(TH2F*)gPad->GetPrimitive("htemp");
	htemp3->SetTitle("dalitz plot of s34 versus s14");
	gPad->Modified();gPad->Update();

	delete rndgen;
}
