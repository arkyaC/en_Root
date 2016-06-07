void mac(){
	TRandom3 rnd;
	TH1F* pX1=new TH1F("part1x","px of particle 1",500,-700,700);
	TH1F* pY1=new TH1F("part1y","py of particle 1",500,-700,700);
	TH1F* pZ1=new TH1F("part1z","pz of particle 1",500,-500,500);
	TH1F* e1=new TH1F("Particle1e","energy of particle1",500,0,1500);
	TH1F* pX2=new TH1F("part2x","px of particle 2",500,-700,700);
	TH1F* pY2=new TH1F("part2y","py of particle 2",500,-700,700);
	TH1F* pZ2=new TH1F("part2z","pz of particle 2",500,-500,500);
	TH1F* e2=new TH1F("Particle2e","energy of particle2",500,0,1500);
	//pX1->SetBit(TH1::kCanRebin);
	//pY1->SetBit(TH1::kCanRebin);
	//pZ1->SetBit(TH1::kCanRebin);
	//e1->SetBit(TH1::kCanRebin);
	TCanvas* c = new TCanvas("particle1","particle1",900,900);
	TCanvas* c1= new TCanvas("particle2","particle2",900,900);
	c->Divide(2,2);
	c1->Divide(2,2);
	double pi=3.14;
	float prapidity,alpha,theta,phi,px1,py1,pz1,energy1,energy2,px2,py2,pz2,gamma,beta,energy,pt,modp,energy;
	float M=1020,m=493.6;//in MeV/c^2
	float pval=0.5*sqrt(((M*M)-(4*m*m)));//Mev/c
	for (int i=0;i<10000;i++){
		pt=100*rnd.Exp(2);//in MeV/c
		prapidity=100*rnd.Gaus(0,0.5);
		alpha=2*(atan(exp(-prapidity)));
		modp=abs(pt/sin(alpha));
		energy=sqrt((modp*modp)+(M*M));
		beta = modp/energy;
		gamma = energy / M;
		phi=2*pi*rnd.Uniform(0,1);
		theta=acos(rnd.Uniform(-1,1));
		px1=(pval*sin(theta)*cos(phi)*((cos(alpha)*cos(alpha))+(gamma*(sin(alpha)*sin(alpha)))))+(pval*cos(theta)*cos(alpha)*sin(alpha)*(gamma-1))+(gamma*beta*M*0.5*sin(alpha));
		py1=pval*sin(theta)*sin(phi);
		pz1=(pval*cos(theta)*((cos(alpha)*cos(alpha))+(gamma*(sin(alpha)*sin(alpha)))))+(pval*(gamma-1)*sin(theta)*cos(phi)*sin(alpha)*cos(alpha))+(gamma*beta*0.5*M*cos(alpha));
		energy1=gamma*((0.5*M)+(beta*pval*(sin(theta)*cos(phi)*sin(alpha)+cos(theta)*cos(alpha))));
		px2=(-pval*sin(theta)*cos(phi)*((cos(alpha)*cos(alpha))+(gamma*(sin(alpha)*sin(alpha)))))+(-pval*cos(theta)*cos(alpha)*sin(alpha)*(gamma-1))+(gamma*beta*M*0.5*sin(alpha));
		py2=-py1;
		pz2=(-pval*cos(theta)*((cos(alpha)*cos(alpha))+(gamma*(sin(alpha)*sin(alpha)))))+(-pval*(gamma-1)*sin(theta)*cos(phi)*sin(alpha)*cos(alpha))+(gamma*beta*0.5*M*cos(alpha));
		energy2	= gamma*((0.5*M)-(beta*pval*(sin(theta)*cos(phi)*sin(alpha)+cos(theta)*cos(alpha))));
		pX1->Fill(px1);
		pY1->Fill(py1);
		pZ1->Fill(pz1);
		e1->Fill(energy1);

		pX2->Fill(px2);
		pY2->Fill(py2);
		pZ2->Fill(pz2);
		e2->Fill(energy2);
		bool b=energy-energy2-energy1<5;
		cout<<b<<endl;
	}
	c->cd(1);
	pX1->Draw();
	c->cd(2);
	pY1->Draw();
	c->cd(3);
	pZ1->Draw();
	c->cd(4);
	e1->Draw();
	c1->cd(1);
	pX2->Draw();
	c1->cd(2);
	pY2->Draw();
	c1->cd(3);
	pZ2->Draw();
	c1->cd(4);
	e2->Draw();
}
