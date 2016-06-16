#define M 1.019445
#define m 0.493667
#define m0 0.493667
#define m1 0.23033
#define m2 0.23033
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

void via2(){
	
	long double ps_rapid,pT;
	double x[1000],y[1000];
	TRandom3* rndgen=new TRandom3();
	for(int i=0;i<1000;i++){
			if(i%100==0)cout<<i<<endl;
			ps_rapid=rndgen->Gaus(0.,3);
			pT=0.2+rndgen->Exp(0.5);
			long double theta_mother=2*atan(exp(-1*ps_rapid));
			long double phi_mother=2*pi*rndgen->Uniform(0,1);
			long double p_mother=pT/sin(theta_mother);
			vector pMother(p_mother*sin(theta_mother)*cos(phi_mother),p_mother*sin(theta_mother)*sin(phi_mother),p_mother*cos(theta_mother));
			long double E=sqrt(p_mother*p_mother+M*M);
			long double beta=p_mother/E;
			vector beta1(beta*sin(theta_mother)*cos(phi_mother),beta*sin(theta_mother)*sin(phi_mother),beta*cos(theta_mother));
			long double p=sqrt(pow((M*M-m*m+m0*m0)/(2*M),2)-m0*m0);
			long double theta=acos(rndgen->Uniform(-1,1));
			long double phi=2*pi*rndgen->Uniform(0,1);
			vector p1cm(p*sin(theta)*cos(phi),p*sin(theta)*sin(phi),p*cos(theta));
			vector p1lab= boost(p1cm,beta1,-(M*M+m*m-m0*m0)/(2*M));
			vector p2cm=p1cm*(-1);
			vector p2lab=boost(p2cm,beta1,-(M*M-m*m+m0*m0)/(2*M));
			
			p1=sqrt(pow((m*m-m2*m2+m1*m1)/(2*m),2)-m1*m1);
			long double energyp1=sqrt(p1lab.value()*p1lab.value()+m*m);
			long double theta1=acos(rndgen->Uniform(-1,1));
			long double phi1=2*pi*rndgen->Uniform(0,1);
			vector p3cm(p1*sin(theta1)*cos(phi1),p1*sin(theta1)*sin(phi1),p1*cos(theta1));
			vector p4cm=p3cm*(-1);
			vector beta2=p1lab*(1/energyp1);
			vector p3lab=boost(p3cm,beta2,-(m*m-m2*m2+m1*m1)/(2*m));
			vector p4lab=boost(p4cm,beta2,-(m*m+m2*m2-m1*m1)/(2*m));
			pMother.print();
			vector su=p3lab+p2lab+ p4lab;
			su.print();
			(p1lab+p2lab).print();
			cout<<E<<endl;
			cout<<p3lab.energ(m1)+p4lab.energ(m2)+p2lab.energ(m0)<<endl;
			cout<<endl;
			x[i]=pow(p4lab.energ(m2)+p2lab.energ(m0),2)-pow((p4lab+p2lab).value(),2);
			y[i]=pow(p3lab.energ(m1)+p2lab.energ(m0),2)-pow((p3lab+p2lab).value(),2);
			
	}
	TGraph* dalitz=new TGraph(1000,x,y);
	TCanvas* c =new TCanvas("sagar","sagar",1500,700);
	dalitz->Draw("A*");

	
}	
