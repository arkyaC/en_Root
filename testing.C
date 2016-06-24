Double_t lorrentzian(Double_t *x, Double_t *par){
	return (0.5*par[0]*par[1]/TMath::Pi()) / TMath::Max(1.e-10,
(x[0]-par[2])*(x[0]-1019.455)+ .25*par[1]*par[1]);
}
Double_t gauss(Double_t *x, Double_t *par){
	return exp(-0.5*pow((x[0]-par[2])/par[1],2));
}

Double_t combination(Double_t *x,Double_t *par){
	TF1* f1= new TF1("lorrentz",lorrentzian,0,2000,3);
	f1->SetParameter(0,par[0]);
	f1->SetParameter(1,par[1]);
	f1->SetParameter(2,par[2]);
	TF1* f= new TF1("gaussian",gauss,0,2000,3);
	f->SetParameter(0,par[3]);
	f->SetParameter(1,par[4]);
	f->SetParameter(2,f1->GetRandom());
	return f->Eval(x[0]);
}
void testing(){
	TF1* d= new TF1("combined",combination,0,2000,5);
	d->SetParameter(0,1);
	d->SetParameter(1,4.43);
	d->SetParameter(2,1019.455);
	d->SetParameter(3,10);
	d->SetParameter(4,100);
	TCanvas* c= new TCanvas("Sagar","Sagar",900,900);
	d->Draw();
}