#define pi 3.141592653589
#define M 493.7
#define m1 139.6
#define m2 139.6
#define m3 135.0
#define n 100000
#define n_bins 150

void P3E3plots(){
	TRandom3* rndgen=new TRandom3();
	float p1max,p2max;
	p1max=sqrt((M*M - ((m1+m2+m3)*(m1+m2+m3)))*(M*M - ((-m1+m2+m3)*(-m1+m2+m3))))/(2*M);
	p2max=sqrt((M*M - ((m1+m2+m3)*(m1+m2+m3)))*(M*M - ((-m2+m1+m3)*(-m2+m1+m3))))/(2*M);
	float ps_rapid=rndgen->Gaus(0.,3);
	float pT=0.2+rndgen->Exp(0.5);
	float theta_mother=2*atan(exp(-1*ps_rapid));
	float phi_mother=2*pi*rndgen->Uniform(0,1);
	float p_mother=pT/sin(theta_mother);
	float E=sqrt(p_mother*p_mother+M*M);
	float beta=p_mother/E;
	float gamma=1./sqrt(1-beta*beta);
	float px_mother=p_mother*sin(theta_mother)*cos(phi_mother);
	float py_mother=p_mother*sin(theta_mother)*sin(phi_mother);
	float pz_mother=p_mother*cos(theta_mother);
	float beta_x=px_mother/E;
	float beta_y=py_mother/E;
	float beta_z=pz_mother/E;
	int r = 0;
	TH1F* P3=new TH1F("P3","P3",n_bins,0.,150.);
	TH1F* e3=new TH1F("e3","e3",n_bins,100.,200.);
	for(int i=0; i<n; i++){
		float p1=rndgen->Gaus(0,p1max);
		float phi1=2*pi*rndgen->Uniform(0,1);
		float theta1=acos(rndgen->Uniform(-1,1));
		float p2=rndgen->Gaus(0,p2max);
		float phi2=2*pi*rndgen->Uniform(0,1);
		float theta2=acos(rndgen->Uniform(-1,1));
		float p1xcm=p1*sin(theta1)*sin(phi1);
		float p1ycm=p1*sin(theta1)*cos(phi1);
		float p1zcm=p1*cos(theta1);
		float p2xcm=p2*sin(theta2)*sin(phi2);
		float p2ycm=p2*sin(theta2)*cos(phi2);
		float p2zcm=p2*cos(theta2);
		float p3xcm=-1.*(p1xcm + p2xcm);
		float p3ycm=-1.*(p1ycm + p2ycm);
		float p3zcm=-1.*(p1zcm + p2zcm);
		float p3cm=sqrt((p3xcm*p3xcm) + (p3ycm*p3ycm) + (p3zcm*p3zcm));
		float E1cm=sqrt(p1*p1 + m1*m1);
		float E2cm=sqrt(p2*p2 + m2*m2);
		float E3cm=sqrt(p3cm*p3cm + m3*m3);
		float Ediff=M-(E1cm+E2cm+E3cm);
		float p3x=(gamma-1)*(p3xcm*beta_x*beta_x/(beta*beta) + p3ycm*beta_y*beta_x/(beta*beta) + p3zcm*beta_z*beta_x/(beta*beta)) + p3xcm + gamma*beta_x*E3cm;
		float p3y=(gamma-1)*(p3xcm*beta_x*beta_y/(beta*beta) + p3ycm*beta_y*beta_y/(beta*beta) + p3zcm*beta_z*beta_y/(beta*beta)) + p3ycm + gamma*beta_x*E3cm;
		float p3z=(gamma-1)*(p3xcm*beta_x*beta_z/(beta*beta) + p3ycm*beta_y*beta_z/(beta*beta) + p3zcm*beta_z*beta_z/(beta*beta)) + p3zcm + gamma*beta_x*E3cm;
		float E3=gamma*(E3cm + beta_x*p3xcm + beta_y*p3ycm + beta_z*p3zcm) ;
		float p3=sqrt((p3x*p3x) + (p3y*p3y) + (p3z*p3z));
		if(Ediff*Ediff < 50){
			P3->Fill(p3);
			e3->Fill(E3);
			r++;
		}
	}
	TCanvas* can2=new TCanvas("3","3",1500,700);
	can2->Divide(2,1);
	can2->cd(1);
	P3->Draw();
	can2->cd(2);
	e3->Draw();
}