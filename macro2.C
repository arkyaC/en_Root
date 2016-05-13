#include "TMath.h"
#include "TGraphPolar.h"
void macro2(){
	auto canv=new TCanvas("Polar plot","Polar plot",600,600);
	Double_t tmin=0.;
	Double_t tmax=TMath::Pi()*6.;
	const Int_t npoints=1000;
	Double_t t[npoints];
	Double_t r[npoints];
	for (Int_t i=0;i<npoints;i++){
		t[i]=tmin+i*(tmax-tmin)/npoints;
		r[i]=TMath::Sin(t[i]);
	}
	TGraphPolar grp1(npoints,t,r);
	grp1.SetTitle("A fan");
	grp1.SetLineWidth(3);
	grp1.SetLineColor(kRed);
	grp1.DrawClone("AL");//A specifies the angle markers

	canv->Print("Polar_Graph.jpg");
}