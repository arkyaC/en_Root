#define M 939.57
#define m1 938.28
#define m2 0.000000320
#define m3  0.511
#define pi 3.141592653589
void tbody(){
	TRandom3 rndgen;
	const int runs=10000;
	float arem12[runs],arem23[runs],arem31[runs];
	int count;
	float p1max,p2max,p3max;
	p1max=sqrt((M*M-pow(m1+m2+m3,2))*(M*M-pow(-m1+m2+m3,2)))/(2*M);
	p2max=sqrt((M*M-pow(m1+m2+m3,2))*(M*M-pow(m1-m2+m3,2)))/(2*M);
	p3max=sqrt((M*M-pow(m1+m2+m3,2))*(M*M-pow(m1+m2-m3,2)))/(2*M);
	TH1F* E3 = new TH1F("E3","E3",150,0.5,1.4);
	TH1F* P3 = new TH1F("P3","P3",150,0,1.2);
	for(int j=0;j<runs;j++){
		float p1=rndgen.Uniform(0,p1max);
		float p2=rndgen.Uniform(0,p2max);
		float p3=rndgen.Uniform(0,p3max);
		float cur1=M-sqrt(p1*p1+m1*m1)-sqrt(p2*p2+m2*m2)-sqrt(p3*p3+m3*m3);
		cout<<cur1<<endl;
		if (cur1*cur1<0.00001 && p1+p2>p3 && p2+p3>p1 && p3+p1>p2){
				arem12[count]=pow(M-sqrt(p3*p3+m3*m3),2)-p3*p3;
				arem23[count]=pow(M-sqrt(p1*p1+m1*m1),2)-p1*p1;
				arem31[count]=pow(M-sqrt(p2*p2+m2*m2),2)-p2*p2;
				E3->Fill(sqrt(p3*p3+m3*m3));
				P3->Fill(p3);
				count++;
			
		}
	}
	TCanvas* c= new TCanvas("sagr","three particle decay",1800,1000);
	c->Divide(3,1);
	TGraph* dalitz1=new TGraph(count,arem12,arem31);
	TGraph* dalitz2=new TGraph(count,arem23,arem12);
	TGraph* dalitz3=new TGraph(count,arem31,arem23);
	c->cd(1);
	dalitz1->Draw("A*");
	c->cd(2);
	dalitz2->Draw("A*");
	c->cd(3);
	dalitz3->Draw("A*");
	TCanvas* c1= new TCanvas("sagar","three particle decay",1800,1000);
	c1->Divide(2,1);
	c1->cd(1);
	E3->Draw();
	c1->cd(2);
	P3->Draw();
	cout<<count<<endl;
}