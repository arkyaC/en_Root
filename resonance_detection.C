//not working
#define n_bins 200

void resonance_detection(){
	TFile in_file("decay_data_phi.root");
	TNtuple* phi_data;
	in_file.GetObject("decay_data",phi_data);

	auto canv=new TCanvas("reso_mass","mass distribution of mother particle;cross-section;mass (in MeV)",900,900);

	TH1F* mass=new TH1F("mother_py","mother particle",n_bins,-800.,800.);
	int n=(int)(phi_data->GetEntries()/2);
	float* k_plus[n];
	float* k_minus[n];
	float px1,px2,py1,py2,pz1,pz2,E1,E2,px_m,py_m,pz_m,E_m;

	for(int i=0;i<phi_data->GetEntries();i++){
		cout<<i<<endl;
		phi_data->GetEntry(i);
		if(i%2==0) cout<<(k_plus[i/2]=phi_data->GetArgs());
		else k_minus[i/2]=phi_data->GetArgs();
	}
	for(int i=0;i<n;i++){
		px1=k_plus[i][0];
		py1=k_plus[i][1];
		pz1=k_plus[i][2];
		E1=k_plus[i][3];
		for(int j=0;j<n;j++){
			cout<<i<<"\t"<<j<<endl;
			px2=k_minus[i][0];
			py2=k_minus[i][1];
			pz2=k_minus[i][2];
			E2=k_minus[i][3];
			E_m=E1+E2;
			px_m=px1+px2;
			py_m=py1+py2;
			pz_m=pz1+pz2;
			mass->Fill(sqrt(E_m*E_m - px_m*px_m - py_m*py_m - pz_m*pz_m));
		}
	}
	canv->Draw();
	//canv->Modified();
}
