/*
1 represents pi+
2 represents pi+
3 represents pi-
*/
#define M 493.6 //Kaon mass
#define m1 139.6 //pi+
#define m2 139.6 //pi+
#define m3 139.6 //pi-

#define pi 3.141592653589
#define nRuns 40000

void KaonThreeBody(){
	TRandom3* rndgen=new TRandom3();
	TNtuple* dalitzPlot21=new TNtuple("dalitz","s2 versus s1","s2:s1");
	TNtuple* dalitzPlot32=new TNtuple("dalitz","s3 versus s2","s3:s2");
	TNtuple* dalitzPlot13=new TNtuple("dalitz","s1 versus s3","s1:s3");
	TCanvas* canv=new TCanvas("dalitz_plot","Dalitz plots",900,300);
	canv->Divide(3,1);

	canv->cd(1);
	gPad->SetTitle("s2 versus s1");
	canv->cd(2);
	gPad->SetTitle("s3 versus s2");
	canv->cd(3);
	gPad->SetTitle("s1 versus s3");

// s2 versus s1
	canv->cd(1);
	for(int i=0;i<nRuns;i++){
		float s1=rndgen->Uniform((m2+m3)*(m2+m3),(M-m1)*(M-m1));
		float s2,s2Min,s2Max,E1o,E3o,p1o,p3o;
		E1o=abs(s1- M*M + m1*m1)/(2*sqrt(s1));
		E3o=abs(s1- m2*m2 + m3*m3)/(2*sqrt(s1));
		p1o=sqrt(E1o*E1o - m1*m1);
		p3o=sqrt(E3o*E3o - m3*m3);
		s2Min=m1*m1 + m3*m3 + 2*(E1o*E3o-p1o*p3o);
		s2Max=m1*m1 + m3*m3 + 2*(E1o*E3o+p1o*p3o);
		s2=rndgen->Uniform(s2Min,s2Max);
		dalitzPlot21->Fill(s2,s1);
		
	}
	dalitzPlot21->Draw("s2:s1");
	TH2F* htemp=(TH2F*)gPad->GetPrimitive("htemp");
	htemp->SetTitle("dalitz plot of s2 versus s1");
	gPad->Modified();gPad->Update();

// s3 versus s2 
	canv->cd(2);
	for(int i=0;i<nRuns;i++){
		float s2=rndgen->Uniform((m3+m1)*(m3+m1),(M-m2)*(M-m2));
		float s3,s3Min,s3Max,E2o,E1o,p2o,p1o;
		E2o=abs(s2- M*M + m2*m2)/(2*sqrt(s2));
		E1o=abs(s2- m3*m3 + m1*m1)/(2*sqrt(s2));
		p2o=sqrt(E2o*E2o - m2*m2);
		p1o=sqrt(E1o*E1o - m1*m1);
		s3Min=m1*m1 + m2*m2 + 2*(E2o*E1o-p2o*p1o);
		s3Max=m1*m1 + m2*m2 + 2*(E2o*E1o+p2o*p1o);
		s3=rndgen->Uniform(s3Min,s3Max);
		dalitzPlot32->Fill(s3,s2);
		
	}
	dalitzPlot32->Draw("s3:s2");
	TH2F* htemp1=(TH2F*)gPad->GetPrimitive("htemp");
	htemp1->SetTitle("dalitz plot of s3 versus s2");
	gPad->Modified();gPad->Update();

// s1 versus s3 1->2 2->3 3->1
	canv->cd(3);
	for(int i=0;i<nRuns;i++){
		float s3=rndgen->Uniform((m1+m2)*(m1+m2),(M-m3)*(M-m3));
		float s1,s1Min,s1Max,E3o,E2o,p3o,p2o;
		E3o=abs(s3- M*M + m3*m3)/(2*sqrt(s3));
		E2o=abs(s3- m1*m1 + m2*m2)/(2*sqrt(s3));
		p3o=sqrt(E3o*E3o - m3*m3);
		p2o=sqrt(E2o*E2o - m2*m2);
		s1Min=m2*m2 + m3*m3 + 2*(E3o*E2o-p3o*p2o);
		s1Max=m2*m2 + m3*m3 + 2*(E3o*E2o+p3o*p2o);
		s1=rndgen->Uniform(s1Min,s1Max);
		dalitzPlot13->Fill(s1,s3);
		
	}
	dalitzPlot13->Draw("s1:s3");
	TH2F* htemp2=(TH2F*)gPad->GetPrimitive("htemp");
	htemp2->SetTitle("dalitz plot of s1 versus s3");
	gPad->Modified();gPad->Update();

	delete rndgen;
}
