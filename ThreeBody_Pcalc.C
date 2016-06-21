#define pi 3.141592653589
#define M 939.57
#define m1 938.28
#define m2 0.000000320
#define m3 0.511
#define n 100000
#define n_bins 150

void ThreeBody_Pcalc(){
	TRandom3* rndgen=new TRandom3();
	TH1F* Mom3=new TH1F("P3","P3",n_bins,0.,1.2);
	TH1F* Energy3=new TH1F("e3","e3",n_bins,0.3,1.4);
	float p1max,p2max;
	p1max=sqrt((M*M - ((m1+m2+m3)*(m1+m2+m3)))*(M*M - ((-m1+m2+m3)*(-m1+m2+m3))))/(2*M);
	p2max=sqrt((M*M - ((m1+m2+m3)*(m1+m2+m3)))*(M*M - ((-m2+m1+m3)*(-m2+m1+m3))))/(2*M);
	int r = 0;
	float x[n],y[n],z[n];
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
		float E1=sqrt(p1cm*p1cm + m1*m1);
		float E2=sqrt(p2cm*p2cm + m2*m2);
		float E3=sqrt(p3cm*p3cm + m3*m3);
		/*//Particle 1
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
		//tests*/
		float Ediff=M-E1-E2-E3;
		cout<<Ediff<<endl;
		if(Ediff*Ediff < 0.00001){
			x[r]=((E1+E2)*(E1+E2))-((p1xcm+p2xcm)*(p1xcm+p2xcm))-((p1ycm+p2ycm)*(p1ycm+p2ycm))-((p1zcm+p2zcm)*(p1zcm+p2zcm));
			y[r]=((E3+E2)*(E3+E2))-((p3xcm+p3xcm)*(p3xcm+p2xcm))-((p3ycm+p2ycm)*(p3ycm+p2ycm))-((p3zcm+p2zcm)*(p3zcm+p2zcm));
			z[r]=((E1+E3)*(E1+E3))-((p1xcm+p3xcm)*(p1xcm+p3xcm))-((p1ycm+p3ycm)*(p1ycm+p3ycm))-((p1zcm+p3zcm)*(p1zcm+p3zcm));
			Mom3->Fill(p3cm);
			Energy3->Fill(E3);
			r++;
		}
	}
	TGraph* dalitz1=new TGraph(r,x,y);
	TGraph* dalitz2=new TGraph(r,y,z);
	TGraph* dalitz3=new TGraph(r,x,z);
	TCanvas* can=new TCanvas("dalitz","dalitz",1500,700);
	TCanvas* can2=new TCanvas("3","3",1500,700);
	can->Divide(3,1);
	can->cd(1);
	dalitz1->Draw("A*");
	can->cd(2);
	dalitz2->Draw("A*");
	can->cd(3);
	dalitz3->Draw("A*");
	can2->Divide(2,1);
	can2->cd(1);
	Mom3->Draw();
	can2->cd(2);
	Energy3->Draw();
	cout<<r<<endl;
}