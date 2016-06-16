#define pi 3.141592653589
#define M 493.7
#define m1 139.6
#define m2 139.6
#define m3 135.0
#define n 100000
#define n_bins 150

void ThreeBody(){
	
	TH1F* inv_mass=new TH1F("inv_mass","inv_mass",n_bins,75000.,140000);
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
	float x[n],y[n];
	for(int i=0; i<n; i++){
		float p1cm=rndgen->Uniform(0,p1max);
		float phi1=2*pi*rndgen->Uniform(0,1);
		float theta1=acos(rndgen->Uniform(-1,1));
		float p2cm=rndgen->Uniform(0,p2max);
		float phi2=2*pi*rndgen->Uniform(0,1);
		float theta2=acos(rndgen->Uniform(-1,1));
		float p1xcm=p1cm*sin(theta1)*sin(phi1);
		float p1ycm=p1cm*sin(theta1)*cos(phi1);
		float p1zcm=p1cm*cos(theta1);
		float p2xcm=p2cm*sin(theta2)*sin(phi2);
		float p2ycm=p2cm*sin(theta2)*cos(phi2);
		float p2zcm=p2cm*cos(theta2);
		float p3xcm=-1.*(p1xcm + p2xcm);
		float p3ycm=-1.*(p1ycm + p2ycm);
		float p3zcm=-1.*(p1zcm + p2zcm);
		float p3cm=sqrt((p3xcm*p3xcm) + (p3ycm*p3ycm) + (p3zcm*p3zcm));
		float E1cm=sqrt(p1cm*p1cm + m1*m1);
		float E2cm=sqrt(p2cm*p2cm + m2*m2);
		float E3cm=sqrt(p3cm*p3cm + m3*m3);
		//Particle 1
		float p1x=(gamma-1)*(p1xcm*beta_x*beta_x/(beta*beta) + p1ycm*beta_y*beta_x/(beta*beta) + p1zcm*beta_z*beta_x/(beta*beta)) + p1xcm + gamma*beta_x*E1cm;
		float p1y=(gamma-1)*(p1xcm*beta_x*beta_y/(beta*beta) + p1ycm*beta_y*beta_y/(beta*beta) + p1zcm*beta_z*beta_y/(beta*beta)) + p1ycm + gamma*beta_x*E1cm;
		float p1z=(gamma-1)*(p1xcm*beta_x*beta_z/(beta*beta) + p1ycm*beta_y*beta_z/(beta*beta) + p1zcm*beta_z*beta_z/(beta*beta)) + p1zcm + gamma*beta_x*E1cm;
		float E1=gamma*(E1cm + beta_x*p1xcm + beta_y*p1ycm + beta_z*p1zcm);
		//Particle 2
		float p2x=(gamma-1)*(p2xcm*beta_x*beta_x/(beta*beta) + p2ycm*beta_y*beta_x/(beta*beta) + p2zcm*beta_z*beta_x/(beta*beta)) + p2xcm + gamma*beta_x*E2cm;
		float p2y=(gamma-1)*(p2xcm*beta_x*beta_y/(beta*beta) + p2ycm*beta_y*beta_y/(beta*beta) + p2zcm*beta_z*beta_y/(beta*beta)) + p2ycm + gamma*beta_x*E2cm;
		float p2z=(gamma-1)*(p2xcm*beta_x*beta_z/(beta*beta) + p2ycm*beta_y*beta_z/(beta*beta) + p2zcm*beta_z*beta_z/(beta*beta)) + p2zcm + gamma*beta_x*E2cm;
		float E2=gamma*(E2cm + beta_x*p2xcm + beta_y*p2ycm + beta_z*p2zcm);
		//Particle 3
		float p3x=(gamma-1)*(p3xcm*beta_x*beta_x/(beta*beta) + p3ycm*beta_y*beta_x/(beta*beta) + p3zcm*beta_z*beta_x/(beta*beta)) + p3xcm + gamma*beta_x*E3cm;
		float p3y=(gamma-1)*(p3xcm*beta_x*beta_y/(beta*beta) + p3ycm*beta_y*beta_y/(beta*beta) + p3zcm*beta_z*beta_y/(beta*beta)) + p3ycm + gamma*beta_x*E3cm;
		float p3z=(gamma-1)*(p3xcm*beta_x*beta_z/(beta*beta) + p3ycm*beta_y*beta_z/(beta*beta) + p3zcm*beta_z*beta_z/(beta*beta)) + p3zcm + gamma*beta_x*E3cm;
		float E3=gamma*(E3cm + beta_x*p3xcm + beta_y*p3ycm + beta_z*p3zcm);
		//tests
		float su=px_mother - (p1x + p2x + p3x);
		float p3=sqrt((p3x*p3x) + (p3y*p3y) + (p3z*p3z));
		float p2=sqrt((p2x*p2x) + (p2y*p2y) + (p2z*p2z));
		float p1=sqrt((p1x*p1x) + (p1y*p1y) + (p1z*p1z));
		float Ediff=E - (E1 + E2 + E3);
		cout<<Ediff<<"  "<<su<<endl;
		float p12x=p1x+p2x;
		float p12y=p1y+p2y;
		float p12z=p1z+p2z;
		float E12=E1+E2;
		float p23x=p3x+p2x;
		float p23y=p3y+p2y;
		float p23z=p3z+p2z;
		float E23=E3+E2;
		if(Ediff*Ediff < 50){
			x[r]=(E12*E12) - (p12x*p12x) - (p12y*p12y) - (p12z*p12z);
			y[r]=(E23*E23) - (p23x*p23x) - (p23y*p23y) - (p23z*p23z);
			r++;
			inv_mass->Fill(y[r-1]);
		}
	}
	TGraph* dalitz=new TGraph(r,x,y);
	TCanvas* can=new TCanvas("dalitz","dalitz",1500,700);
	can->Divide(2,1);
	can->cd(1);
	dalitz->Draw("A*");
	can->cd(2);
	inv_mass->Draw();
	cout<<r<<endl;
}