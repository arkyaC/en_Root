#define pi 3.141592653589
#define M 939.57
#define m1 938.28
#define m2 0.000000320
#define m3 0.511
#define n 1000000
#define n_bins 150

void ThreeBodyinCM(){
	TRandom3* rndgen=new TRandom3();
	TH1F* m12=new TH1F("m12","m12",n_bins,75000.,140000.);
	TH1F* m23=new TH1F("m23","m23",n_bins,75000.,140000.);
	TH1F* m13=new TH1F("m13","m13",n_bins,75000.,140000.);
	//TH1F* Mom3=new TH1F("P3","P3",n_bins,0.,150.);
	//TH1F* Energy3=new TH1F("e3","e3",n_bins,100.,200.);
	float p1max,p2max;
	p1max=sqrt((M*M - ((m1+m2+m3)*(m1+m2+m3)))*(M*M - ((-m1+m2+m3)*(-m1+m2+m3))))/(2*M);
	p2max=sqrt((M*M - ((m1+m2+m3)*(m1+m2+m3)))*(M*M - ((-m2+m1+m3)*(-m2+m1+m3))))/(2*M);
	int r=0;
	float x12[n], y23[n], z13[n];
	for(int i=0; i<n; i++){
		float p1=rndgen->Uniform(0,p1max);
		float phi1=2*pi*rndgen->Uniform(0,1);
		float theta1=acos(rndgen->Uniform(-1,1));
		float p2=rndgen->Uniform(0,p2max);
		float phi2=2*pi*rndgen->Uniform(0,1);
		float theta2=acos(rndgen->Uniform(-1,1));
		float p1x=p1*sin(theta1)*sin(phi1);
		float p1y=p1*sin(theta1)*cos(phi1);
		float p1z=p1*cos(theta1);
		float p2x=p2*sin(theta2)*sin(phi2);
		float p2y=p2*sin(theta2)*cos(phi2);
		float p2z=p2*cos(theta2);
		float p3x=-1.*(p1x + p2x);
		float p3y=-1.*(p1y + p2y);
		float p3z=-1.*(p1z + p2z);
		float p3=sqrt((p3x*p3x) + (p3y*p3y) + (p3z*p3z));
		float E1=sqrt(p1*p1 + m1*m1);
		float E2=sqrt(p2*p2 + m2*m2);
		float E3=sqrt(p3*p3 + m3*m3);
		float Ediff=M-(E1+E2+E3);
		cout<<Ediff<<endl;
		if(Ediff*Ediff < 0.00001){
			x12[r]=((E1+E2)*(E1+E2)) - ((p1x+p2x)*(p1x+p2x)) - ((p1y+p2y)*(p1y+p2y)) - ((p1z+p2z)*(p1z+p2z));
			y23[r]=((E3+E2)*(E3+E2)) - ((p3x+p2x)*(p3x+p2x)) - ((p3y+p2y)*(p3y+p2y)) - ((p3z+p2z)*(p3z+p2z));
			z13[r]=((E1+E3)*(E1+E3)) - ((p1x+p3x)*(p1x+p3x)) - ((p1y+p3y)*(p1y+p3y)) - ((p1z+p3z)*(p1z+p3z));
			m12->Fill(x12[r]);
			m23->Fill(y23[r]);
			m13->Fill(z13[r]);
			//Mom3->Fill(p3);
			//Energy3->Fill(E3);
			r++;
		}
	}
	TGraph* dalitz1=new TGraph(r,x12,y23);
	TGraph* dalitz2=new TGraph(r,x12,z13);
	TGraph* dalitz3=new TGraph(r,y23,z13);
	TCanvas* can=new TCanvas("dalitz","dalitz",1500,700);
	//TCanvas* can1=new TCanvas("inv_mass","inv_mass",1500,700);
	//TCanvas* can2=new TCanvas("3","3",1500,700);
	//can2->Divide(2,1);
	//can2->cd(1);
	//Mom3->Draw();
	//can2->cd(2);
	//Energy3->Draw();
	/*can1->Divide(3,1);
	can1->cd(1);
	m12->Draw();
	can1->cd(2);
	m23->Draw();
	can1->cd(3);
	m13->Draw();
	*/
	can->Divide(3,1);
	can->cd(1);
	dalitz1->Draw("A*");
	can->cd(2);
	dalitz2->Draw("A*");
	can->cd(3);
	dalitz3->Draw("A*");
	cout<<r;
}