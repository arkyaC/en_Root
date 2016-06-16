
#define M 493.7
#define m1 139.6
#define m2 139.6
#define m3  135.0
#define pi 3.141592653589
double* subtract(double* are,double* are1){
	double* are2 = new double[4];
	for (int i=0;i<4;i++) are2[i]=are[i]-are1[i];
	return are2;
}
double * multiply(double* are,double fac){
	double* are2=new double[4];
	for(int i=0;i<4;i++){
		are2[i]=are[i]*fac;
	}
	return are2;
}
void prin(double are[4]){
	for(int i=0;i<4;i++) cout<<are[i]<<" ";
	cout<<endl;
}
void printer(double **are){
	for (int i=0;i<3;i++){
		for (int j=0;j<4;j++)
			cout<<are[i][j]<<"  ";
		cout<<endl;
	}
}
void solver(double **are,double &x,double &y, double &z){
	int c;
	for (int i=0;i<3;i++){
		if (are[i][0]!=0){
			double fac=are[i][0];
			for (int j=0; j<4;j++)	are[i][j]=are[i][j]/fac;
			for (int k=0; k<3; k++){ 
				if(k!=i) are[k]=subtract(are[k],multiply(are[i],are[k][0]));
			}
			for (int j=0;j<3;j++){
				if (are[j][1]!=0 && j!=i ){
					c=3-i-j;
					double fac1=are[j][1];
					for (int k=1; k<4;k++)	are[j][k]=are[j][k]/fac1;
					are[c]=subtract(are[c],multiply(are[j],are[c][1]));
					if (are[c][2]!=0){
						z=are[c][3]/are[c][2];
						y=are[j][3]-z*are[j][2];
						x=are[i][3]-z*are[i][2]-y*are[i][1];
					}
				}
				
			}
			return;
		}
		
	}	
}
double f1(double alpha,double beta,double gamma,double theta, double phi){
	return m1*tan(alpha)+m2*cos(theta)*tan(beta)+m3*cos(phi)*tan(gamma);
}
double f2(double alpha,double beta,double gamma,double theta, double phi){
	return m2*sin(theta)*tan(beta)-m3*sin(phi)*tan(gamma);
}
double f3(double alpha,double beta,double gamma){
	return (m1/cos(alpha))+(m1/cos(beta))+(m1/cos(gamma))-M;
}
double f1x(double alpha){
	return m1/pow(cos(alpha),2);
}
double f1y(double beta,double theta){
	return m2*cos(theta)/pow(cos(beta),2);
}
double f1z(double gamma,double phi){
	return m3*cos(phi)/pow(cos(gamma),2);
}
double f2y(double beta,double theta){
	return m2*sin(theta)/pow(cos(beta),2);
}
double f2z(double gamma,double phi){
	return -m3*sin(phi)/pow(cos(gamma),2);
}
double f3x(double alpha){
	return m1*tan(alpha)/cos(alpha);
}
double f3y(double alpha){
	return m2*tan(alpha)/cos(alpha);
}
double f3z(double alpha){
	return m3*tan(alpha)/cos(alpha);
}
void dec(){
	
	TH1F* m23=new TH1F("sagar","m12",100,78400,129600);
	TRandom3 rndgen;
	double are55[1800],are66[1800];
	int count=0;
	for(int u=0; u<1800;u++){
		double x=pi/4,y=pi/4,z=pi/4;
		cout<<u<<endl;
		while(1){
			double theta=pi*rndgen.Uniform(0,1);
			double phi=pi*rndgen.Uniform(0,1);
			//cout<<theta<<" "<<phi<<endl;	
			for (int po=0;po<20;po++){
				double dx=0,dy=0,dz=0;
				double are1[]={f1x(x),f1y(y,theta),-f1z(z,phi),-f1(x,y,z,theta,phi)};
				double are2[]={0,f2y(y,theta),f2z(z,phi),-f2(x,y,z,theta,phi)};
				double are3[]={f3x(x),f3y(y),f3z(z),-f3(x,y,z)};
				double* are[3]={are1,are2,are3};
				solver(are,dx,dy,dz);
				x=x+dx;y=y+dy;z=z+dz;
				double p=f1(x,y,z,theta,phi),q=f2(x,y,z,theta,phi),r=f3(x,y,z);
				if(p*p+q*q+r*r<0.01) {break;}
			}
			if(tan(x)<0 | tan(y)<0 | tan(y)<0 | sin(x)<0 | sin(y)<0 | sin(z)<0   ) {cout<<"chutiyapa"<<endl;continue;}
			//cout<<endl<<m1*tan(x)<<" "<<m2*tan(y)<<" "<<m3*tan(z)<<endl;
			//cout<<endl<<m1/cos(x)+m2/cos(y)+m3/cos(z)<<endl;
			
			double mas23=sqrt(pow(M-(m1/cos(x)),2)-pow(m1*tan(x),2));
			double mas13=sqrt(pow(M-(m2/cos(y)),2)-pow(m2*tan(y),2));
			
			bool b= (mas23>m2+m3) && (mas23<M-m1);
			cout<<b<<endl;
			double p=f1(x,y,z,theta,phi),q=f2(x,y,z,theta,phi),r=f3(x,y,z);
			if(p*p+q*q+r*r<0.01) {
				m23->Fill(mas23*mas23);
				count++;
				are55[count-1]=mas23*mas23;
				are66[count-1]=mas13*mas13;
			}

			break;
		}
	}
	TGraph* dalitz=new TGraph(count,are66,are55);
	TCanvas* c3 =new TCanvas("sagar","sagar",1500,700);
	c3->Divide(2,1);
	c3->cd(1);
	dalitz->Draw("AF");
	c3->cd(2);
	m23->Draw();
	
}
