#define pi 3.141592653589
#define M
#define m1
#define m2
#define m3
void ThreebodyDecay(){
	double px1,py1,pz1,py2,px2,pz2,p1max;
	Trandom3* rndgen=new TRandom3();
	p1max=sqrt((M*M - pow(m1+m2+m3,2))*(M*M - pow(-m1+m2+m3,2)))/(2*M);
	for(int i=0; i<100; i++){
		double p1cm=rndgen->Uniform(0,p1max);
		double phi1=2*pi*rndgen->Uniform(0,1);
	}
}