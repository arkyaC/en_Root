#define M 10.19445
#define m 4.93667
#define m1 2.3033
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
		double energ(double mass){
			return sqrt(pow(value(),2)+mass*mass);
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

void via21(){
	TCanvas* c =new TCanvas("sagar","sagar",1500,700);
	long double ps_rapid,pT;
	TH1F* inmass= new TH1F("inmass","inmass",100,20,10000);
	//inmass.SetBit(TH1::kCanRebin);
	TRandom3* rndgen=new TRandom3();
	for(int i=0;i<200;i++){
		if(i%100==0)cout<<i<<endl;
		
		
		vector are[20];
		vector are1[20];
		int r =(int)(rndgen->Gaus(5,2))+1;
		if(r<0) continue;
		cout<< r<<endl;
		for (int j=0;j<10 	;j++){
			ps_rapid=rndgen->Gaus(0.,3);
			pT=rndgen->Exp(2000.);
			long double theta_mother=2*atan(exp(-1*ps_rapid));
			long double phi_mother=2*pi*rndgen->Uniform(0,1);
			long double p_mother=pT/sin(theta_mother);
			vector pMother(p_mother*sin(theta_mother)*cos(phi_mother),p_mother*sin(theta_mother)*sin(phi_mother),p_mother*cos(theta_mother));
			long double E=sqrt(p_mother*p_mother+M*M);
			long double beta=p_mother/E;
			vector beta1(beta*sin(theta_mother)*cos(phi_mother),beta*sin(theta_mother)*sin(phi_mother),beta*cos(theta_mother));
			long double p=0.5*sqrt(M*M - 4*m*m);
			long double theta=acos(rndgen->Uniform(-1,1));
			long double phi=2*pi*rndgen->Uniform(0,1);
			vector p1cm(p*sin(theta)*cos(phi),p*sin(theta)*sin(phi),p*cos(theta));
			vector p1lab= boost(p1cm,beta1,-M/2);
			vector p2cm=p1cm*(-1);
			vector p2lab=boost(p2cm,beta1,-M/2);
			/*
			p1=0.5*sqrt(m*m-4*m1*m1);
			long double energyp1=sqrt(p1lab.value()*p1lab.value()+m*m);
			long double theta1=acos(rndgen->Uniform(-1,1));
			long double phi1=2*pi*rndgen->Uniform(0,1);
			vector p3cm(p1*sin(theta1)*cos(phi1),p1*sin(theta1)*sin(phi1),p1*cos(theta1));
			vector p4cm=p3cm*(-1);
			vector beta2=p1lab*(1/energyp1);
			vector p3lab=boost(p3cm,beta2,-m/2);
			vector p4lab=boost(p4cm,beta2,-m/2);
			*//*
			pMother.print();
			vector su=p1lab+p2lab;
			su.print();
			cout<<E<<endl;
			cout<<p1lab.energ(m)+p2lab.energ(m)<<endl;
			cout<<endl;*/
			are[j]=p1lab;
			are1[j]=p2lab;
			if(j==r) break;
		}
		
		for (int k=0;k<r; k++){
			for(int l=0;l<r ; l++){
				vector pdash= are[k]+are1[l];
				double edash =are[k].energ(m)+are1[l].energ(m);

				inmass.Fill(sqrt(-pow(pdash.value(),2)+edash*edash));
			}
		}
		

	}
	inmass.Draw();
	
}	