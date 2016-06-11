//not working
#define n_bins 200

void resonance_detection(){
	TFile in_file("decay_data_phi.root");
	TNtuple* phi_data;
	in_file.GetObject("decay_data",phi_data);
	gROOT->cd();

	TCanvas* canv=new TCanvas("reso_mass","mass distribution of mother particle;cross-section;mass (in MeV)",900,900);

	TH1F* mass=new TH1F("mother_py","mother particle",n_bins,-800.,800.);
	int n=(int)(phi_data->GetEntries()/2);
	//cout<<n<<endl;
	float** k_plus;
	float** k_minus;
	k_plus=new float*[n];
	k_minus=new float*[n];
	for(int i=0;i<n;i++){
		k_plus[i]=new float[4];
		k_minus[i]=new float[4];
	}
	float px1,px2,py1,py2,pz1,pz2,E1,E2,px_m,py_m,pz_m,E_m;
	//cout<<n<<endl;
	for(int i=0;i<phi_data->GetEntries();i++){
		//cout<<i<<endl;
		phi_data->GetEntry(i);
		if(i%2==0) k_plus[i/2]=phi_data->GetArgs();
		else k_minus[i/2]=phi_data->GetArgs();
	}
	for(int i=0;i<n;i++){
		//if (i%10000==0) cout<<i<<endl;
		px1=k_plus[i][0];
		py1=k_plus[i][1];
		pz1=k_plus[i][2];
		E1=k_plus[i][3];
		if(i%100000==0) cout<<i<<endl;
		for(int j=0;j<n;j++){
			
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

	mass->Draw();
	canv->cd(0);canv->Modified();canv->Update();
	//canv->Modified();
	for(int i=0;i<n;i++){
		delete k_plus[i];
		delete k_minus[i];
	}
	delete[] k_plus;
	delete[] k_minus;
}
