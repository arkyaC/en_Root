#define pi 3.141592653589
#define M 939.57
#define m1 938.28
#define m2 0.000000320
#define m3 0.511
#define n 10000000
#define n_bins 150

void ThreeBody_Allrand(){
	TH1F* Mom3=new TH1F("P3","P3",n_bins,0.,1.2);
	TH1F* Energy3=new TH1F("e3","e3",n_bins,0.3,1.4);
	float p1max,p2max,p3max;
	p1max=sqrt((M*M - ((m1+m2+m3)*(m1+m2+m3)))*(M*M - ((-m1+m2+m3)*(-m1+m2+m3))))/(2*M);
	p2max=sqrt((M*M - ((m1+m2+m3)*(m1+m2+m3)))*(M*M - ((-m2+m1+m3)*(-m2+m1+m3))))/(2*M);
	p3max=sqrt((M*M - ((m1+m2+m3)*(m1+m2+m3)))*(M*M - ((-m3+m2+m1)*(-m3+m2+m1))))/(2*M);
	int r=0;
	float x12[n], y23[n], z13[n];
	for(int i=0; i<n; i++){
		float p1=rndgen->Uniform(0,p1max);
		float p2=rndgen->Uniform(0,p2max);
		float p3=rndgen->Uniform(0,p3max);
		float Ediff=M-(sqrt(p1*p1 + m1*m1)+sqrt(p2*p2 + m2*m2)+sqrt(p3*p3 + m3*m3));
		cout<<Ediff<<endl;
		if(p1+p2>p3 && p2+p3>p1 && p3+p1>p2 && Ediff*Ediff < 0.0001){
			x12[r]=(M-sqrt(p3*p3+m3*m3))*(M-sqrt(p3*p3+m3*m3)) - p3*p3;
			y23[r]=(M-sqrt(p1*p1+m1*m1))*(M-sqrt(p1*p1+m1*m1)) - p1*p1;
			z13[r]=(M-sqrt(p2*p2+m2*m2))*(M-sqrt(p2*p2+m2*m2)) - p2*p2;;
			Mom3->Fill(p3);
			Energy3->Fill(sqrt(p3*p3 + m3*m3));
			r++;
		}
	}
	TGraph* dalitz1=new TGraph(r,x12,y23);
	TGraph* dalitz2=new TGraph(r,x12,z13);
	TGraph* dalitz3=new TGraph(r,y23,z13);
	TCanvas* can=new TCanvas("dalitz","dalitz",1500,700);
	TCanvas* can2=new TCanvas("3","3",1500,700);
	can2->Divide(2,1);
	can2->cd(1);
	Mom3->Draw();
	can2->cd(2);
	Energy3->Draw();
	can->Divide(3,1);
	can->cd(1);
	dalitz1->Draw("A*");
	can->cd(2);
	dalitz2->Draw("A*");
	can->cd(3);
	dalitz3->Draw("A*");
	float p=r;
	cout<<r<<"   "<<p/n<<endl;
}