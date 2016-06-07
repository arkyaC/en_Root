void mac(){
	TRandom3 rnd;
	TH1F* pX1=new TH1F("part1x","px of particle 1",500,0,1500);
	//TH1F* pY1=new TH1F("part1y","py of particle 1",500,0,1500);
	TH1F* pZ1=new TH1F("part1z","pz of particle 1",500,0,1500);
	TH1F* e1=new TH1F("Particle1e","energy of particle1",5000,0,1500);
	pX1->SetBit(TH1::kCanRebin);
//	pY1->SetBit(TH1::kCanRebin);
	pZ1->SetBit(TH1::kCanRebin);
	e1->SetBit(TH1::kCanRebin);
	TCanvas* c = new TCanvas("sagar","sagar",900,900);
	c->Divide(2,2);
	double pi=3.14;
	float prapidity,alpha,theta,phi,px1,py1,pz1,energy1,energy2,px2,py2,pz2,gamma,beta,energy,pt,modp,energy;
	float M=1020,m=493.6;//in MeV/c^2
	float pval=0.5*sqrt(((M*M)-(4*m*m)));//Mev/c
	for (int i=0;i<100000;i++){
		pt=100*rnd.Exp(2);//in MeV/c
		prapidity=100*rnd.Gaus(0,0.5);
		alpha=2*(atan(exp(-prapidity)));
		modp=abs(pt/sin(alpha));
		energy=sqrt((modp*modp)+(M*M));
		beta = modp/energy;
		gamma = energy / M;
		phi=2*pi*rnd.Uniform(0,1);
		theta=acos(rnd.Uniform(0,1));
		px1=(pval*sin(theta)*cos(phi)*((cos(alpha)*cos(alpha))+(gamma*(sin(alpha)*sin(alpha)))))+(pval*cos(theta)*cos(alpha)*sin(alpha)*(gamma-1))+(gamma*beta*M*0.5*sin(alpha));
		py1=pval*sin(theta)*sin(phi);
		pz1=(pval*cos(theta)*((cos(alpha)*cos(alpha))+(gamma*(sin(alpha)*sin(alpha)))))+(pval*(gamma-1)*sin(theta)*cos(phi)*sin(alpha)*cos(alpha))+(gamma*beta*0.5*M*cos(alpha));
		energy1=gamma*((0.5*M)+(beta*pval*(sin(theta)*cos(phi)*sin(alpha)+cos(theta)*cos(alpha))));
		cout<<energy1<<endl;
		pX1->Fill(sqrt((px1*px1)+(py1*py1)));
		//pY1->Fill(py1);
		pZ1->Fill(pz1);
		e1->Fill(energy1);
	}
c->cd(1);
pX1->Draw();
c->cd(2);
//pY1->Draw();
c->cd(3);
pZ1->Draw();
c->cd(4);
e1->Draw();
}
