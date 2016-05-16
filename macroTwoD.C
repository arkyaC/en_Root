int macroTwoD(){
	gStyle->SetPalette(kBird);//sets color scheme for 2D fitted plot
	const double e=0.3;//arbitrary error value
	const int nd=500;

	TRandom3 rand_gen;
	TF2 f2("f2","1000*(([0]*sin(x)/x)*([1]*sin(y)/y))",-5,5,-5,5);
	f2.SetParameters(1,1);
	TGraph2DErrors grp(nd);
	double rnd,z,x,y,ex,ey,ez;
	for(Int_t i=0;i<nd;i++){
		f2.GetRandom2(x,y);//choose a random point in xy plane
		rnd=rand_gen.Uniform(-e,e);//choose a random error in z
		z=f2.Eval(x,y)*(1+rnd);
		grp.SetPoint(i,x,y,z);
		ex=.05*rand_gen.Uniform();//generates random number in [0,1]
		ey=.05*rand_gen.Uniform();
		ez=fabs(z*rnd);
		grp.SetPointError(i,ex,ey,ez);
	}

	f2.SetParameters(0.7,1.5);//initial params for fitting
	f2.SetTitle("Fitted 2D function");
	grp.Fit(&f2);//fitting f2 to graph data

	auto canv=new TCanvas("2D graph with fit","2D graph with fit",600,600);
	canv->SetGrid();
	f2.SetLineWidth(1);
	f2.SetLineColor(kBlue-5);//sets colour for the mesh produced due to surf1

	TF2 *f2c = (TF2*)f2.DrawClone("Surf1");
	TAxis *Xaxis = f2c->GetXaxis();
	TAxis *Yaxis = f2c->GetYaxis();
	TAxis *Zaxis = f2c->GetZaxis();
	Xaxis->SetTitle("X Title"); Xaxis->SetTitleOffset(1.5);
	Yaxis->SetTitle("Y Title"); Yaxis->SetTitleOffset(1.5);
	Zaxis->SetTitle("Z Title"); Zaxis->SetTitleOffset(1.5);
	grp.DrawClone("P0 Same");
	//new window to display x and y projections of fit
	auto new_canv=new TCanvas("ProjCan","The Projections",1000,400);
	new_canv->Divide(2,1);
	new_canv->cd(1);
	grp.Project("x")->Draw();
	new_canv->cd(2);
	grp.Project("y")->Draw();

	//canv->Print("with projections.pdf");
	//new_canv->Print("the projections.pdf");
	return 0;
}