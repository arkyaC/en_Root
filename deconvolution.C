#define M0 1019.445
#define m 0.510999
#define pi 3.141592653589
#define n_bins 200
#define n 100000
#define e 2.71828183
#define width 4.43
Double_t lorrentzian(Double_t *x, Double_t *par){
	return (0.5*par[0]*par[1]/TMath::Pi()) / TMath::Max(1.e-10,(x[0]-par[2])*(x[0]-par[2])+ .25*par[1]*par[1]);
}
void deconvolution(){
	
	TH1F* inv_mass=new TH1F("inv_mass","inv_mass",n_bins,950,1090);
	TRandom3* rndgen=new TRandom3();
	float ps_rapid,pT;
	TF1* f= new TF1("lorrentz",lorrentzian,0,2000,3);
	f->SetParameter(0,1);
	f->SetParameter(1,width);
	f->SetParameter(2,M0);
	for(int i=0; i<n; i++)
		Double_t M=f->GetRandom();
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
	    //float E1_lab=(gamma*(M/2))  +  gamma*((px1*beta_x) + (py1*beta_y) + (pz1*beta_z));
	    float pp1=sqrt((px1_lab)*(px1_lab) + (py1_lab)*(py1_lab) + (pz1_lab)*(pz1_lab));
	    float p1_lab=rndgen->Gaus(pp1,0.05*pp1);
	    float E1_lab=sqrt(m*m + p1_lab*p1_lab);
	    px1=-px1;
	    py1=-py1;
	    pz1=-pz1;
	    float px2_lab=(gamma-1)*((px1*beta_x*beta_x)/(beta*beta) + (py1*beta_x*beta_y)/(beta*beta) + (pz1*beta_x*beta_z)/(beta*beta))  +  px1  +  gamma*beta_x*(M/2);
	    float py2_lab=(gamma-1)*((px1*beta_y*beta_x)/(beta*beta) + (py1*beta_y*beta_y)/(beta*beta) + (pz1*beta_y*beta_z)/(beta*beta))  +  py1  +  gamma*beta_y*(M/2);
	    float pz2_lab=(gamma-1)*((px1*beta_z*beta_x)/(beta*beta) + (py1*beta_z*beta_y)/(beta*beta) + (pz1*beta_z*beta_z)/(beta*beta))  +  pz1  +  gamma*beta_z*(M/2);
	    //float E2_lab=(gamma*(M/2))  +  gamma*((px1*beta_x) + (py1*beta_y) + (pz1*beta_z));
	    float pp2=sqrt((px2_lab)*(px2_lab) + (py2_lab)*(py2_lab) + (pz2_lab)*(pz2_lab));
	    float costheta=(px1_lab*px2_lab+py1_lab*py2_lab+pz1_lab*pz2_lab)/(pp1*pp2);
	    float p2_lab=rndgen->Gaus(pp2,0.05*pp2);
	    float E2_lab=sqrt(m*m + p2_lab*p2_lab);
		inv_mass->Fill(sqrt(pow(E1_lab+E2_lab,2)-p1_lab*p1_lab-p2_lab*p2_lab-(2*p1_lab*p2_lab*costheta)));

		
	}
	auto canv1=new TCanvas("inv_ass","inv_mss",1500,900);
	canv1->Divide(2,1);
	canv1->cd(1);
	inv_mass->Draw();
	TFitResultPtr r= inv_mass->Fit("gaus","S");
	double par0=r->Parameter(0);
	double par1=r->Parameter(1);
	double par2=r->Parameter(2);
	TF1* g=new TF1("achha","[0]*exp(-0.5*(x-[1])*(x-[1])/([2]*[2]))",950,1090);
	g->SetParameter(0,par0);
	g->SetParameter(1,par1);
	g->SetParameter(2,par2);
	int size = inv_mass->GetXaxis()->GetNbins();
	float *source,*resp;
	source= new float[size];
	resp= new float[size];
	for(int i=0;i<size;i++){
		source[i]=inv_mass->GetBinContent(i+1);
		resp[i]=g->Eval(950.35+0.7*i);
	}
	TSpectrum decon;
	cout<<"chut"<<endl;
	decon.Deconvolution(source,resp,size,3,9,1);
	cout<<endl<<par0<<endl<<par1<<endl<<par2<<endl;
	TH1F *hout = new TH1F("hout","hout",size,1115,1030);
	for(int i=0;i<size;i++){
		hout->Fill(950.35+0.7*i,source[i]);
	}
	canv1->cd(2);
	hout->Draw();
	hout->Fit("lorrentz","SR");
	cout<<"gfdgsd"<<endl;
	canv1->Modified();

}