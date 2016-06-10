#define M 1019.445
#define m 493.667
#define m1 230.33
#define pi 3.141592653589
class vector{
	double x,y,z;
	public:
		vector(){
			x=0;y=0;z=0;
		}
		vector(double a,double b, double c)	{
			x=a;y=b;z=c;
		}
		void set(double a=0, double b=0, double c=0){
			x=a;y=b;z=c;
		}
		double value(){
			return sqrt(x*x+y*y+z*z);
		}
		double x(){
			return x;
		}
		double y(){
			return y;
		}
		double z(){
			return z;
		} 
		vector operator+ (const vector &b) const{
			return vector(x+b.x,y+b.y,z+b.z);
		}
		vector operator* (double t) const{
			return vector(x*t,y*t,z*t);
		}
		void print(){
			cout<<x<<"  "<<y<<"   "<<z<<endl;
		}

};
double dot(vector a, vector b){
	return a.x()*b.x()+a.y()*b.y()+a.z()*b.z();
}
vector boost(vector p,vector beta,double energy){
	double gamma= 1/sqrt(1-beta.value()*beta.value());
	double f=gamma*((gamma/(gamma+1)*dot(beta,p))-energy);
	return p+(beta*f);
}
void init(TCanvas* canv){
	canv->Divide(1,3);
	canv->cd(1);
	gPad->SetTitle("particle 1");
	canv->cd(2);
	gPad->SetTitle("particle 2");
	canv->cd(3);
	gPad->SetTitle("mother");
}
void via2(){
	/*
	auto canv1=new TCanvas("p_x","p_x",900,900);
	auto canv2=new TCanvas("p_y","p_y",900,900);
	auto canv3=new TCanvas("p_z","p_z",900,900);
	auto canv4=new TCanvas("E","E",900,900);
	init(canv1);
	init(canv2);
	init(canv3);
	init(canv4);
	TH1F* par1_E=new TH1F("particle_1_E","particle 1",200,350.,1600.);
	TH1F* par2_E=new TH1F("particle_2_E","particle 2",200,350.,1600.);
	TH1F* mother_E=new TH1F("mother_E","mother particle",200,350.,1600.);
	TH1F* par1_px=new TH1F("particle_1_px","particle 1",200,-700.,700.);
	TH1F* par2_px=new TH1F("particle_2_px","particle 2",200,-700.,700.);
	TH1F* mother_px=new TH1F("mother_px","mother particle",200,0.,800.);
	TH1F* par1_py=new TH1F("particle_1_py","particle 1",200,-700.,700.);
	TH1F* par2_py=new TH1F("particle_2_py","particle 2",200,-700.,700.);
	TH1F* mother_py=new TH1F("mother_py","mother particle",200,0.,800.);
	TH1F* par1_pz=new TH1F("particle_1_pz","particle 1",200,-700.,700.);
	TH1F* par2_pz=new TH1F("particle_2_pz","particle 2",200,-700.,700.);
	TH1F* mother_pz=new TH1F("mother_pz","mother particle",200,0.,1400.);
*/
	float ps_rapid,pT;
	
	TRandom3* rndgen=new TRandom3();
	for(int i=0;i<10;i++){
		if(i%100==0)cout<<i<<endl;
		ps_rapid=rndgen->Gaus(0.,3);
		pT=rndgen->Exp(2000.);
		float theta_mother=2*atan(exp(-1*ps_rapid));
		float phi_mother=2*pi*rndgen->Uniform(0,1);
		float p_mother=pT/sin(theta_mother);
		vector pMother(p_mother*sin(theta_mother)*cos(phi_mother),p_mother*sin(theta_mother)*sin(phi_mother),p_mother*cos(theta_mother));
		float E=sqrt(p_mother*p_mother+M*M);
		float beta=p_mother/E;
		vector beta1(beta*sin(theta_mother)*cos(phi_mother),beta*sin(theta_mother)*sin(phi_mother),beta*cos(theta_mother));
		float p=0.5*sqrt(M*M - 4*m*m);
		float theta=acos(rndgen->Uniform(-1,1));
		float phi=2*pi*rndgen->Uniform(0,1);
		vector p1cm(p*sin(theta)*cos(phi),p*sin(theta)*sin(phi),p*cos(theta));
		vector p1lab= boost(p1cm,beta1,-M/2);
		vector p2cm=p1cm*(-1);
		vector p2lab=boost(p2cm,beta1,-M/2);
		/*
		mother_px->Fill(pMother.x());
		mother_py->Fill(pMother.y());
		mother_pz->Fill(pMother.z());
		par1_px->Fill(p1lab.x());
		par1_py->Fill(p1lab.y());
		par1_pz->Fill(p1lab.z());
		par1_E->Fill(sqrt(pow(p1lab.value(),2)+m*m));
		*/
		p1=0.5*sqrt(m*m-4*m1*m1);
		float energyp1=sqrt(p1lab.value()*p1lab.value()+m*m);
		float theta1=acos(rndgen->Uniform(-1,1));
		float phi1=2*pi*rndgen->Uniform(0,1);
		vector p3cm(p1*sin(theta1)*cos(phi1),p1*sin(theta1)*sin(phi1),p1*cos(theta1));
		vector p4cm=p3cm*(-1);
		vector beta2=p1lab*(1/energyp1);
		vector p3lab=boost(p3cm,beta2,-m/2);
		vector p4lab=boost(p4cm,beta2,-m/2);
		
		pMother.print();
		vector su=p4lab+p3lab+p2lab;
		su.print();
		cout<<E<<endl;
		cout<<sqrt(pow(p4lab.value(),2)+m1*m1)+sqrt(pow(p3lab.value(),2)+m1*m1)+sqrt(pow(p2lab.value(),2)+m*m)<<endl;
		cout<<endl;
	}
	/*
	canv1->cd(1);
	par1_px->Draw();
	canv1->cd(2);
	par2_px->Draw();
	canv1->cd(3);
	mother_px->Draw();
	canv2->cd(1);
	par1_py->Draw();
	canv2->cd(2);
	par2_py->Draw();
	canv2->cd(3);
	mother_py->Draw();
	canv3->cd(1);
	par1_pz->Draw();
	canv3->cd(2);
	par2_pz->Draw();
	canv3->cd(3);
	mother_pz->Draw();
	canv4->cd(1);
	par1_E->Draw();
	canv4->cd(2);
	par2_E->Draw();
	canv4->cd(3);
	mother_E->Draw();
*/
}	