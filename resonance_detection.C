
#define n_bins 300

void resonance_detection(){
	TFile in_file("decay_data_phi.root");
	TNtuple* phi_data;
	in_file.GetObject("decay_data",phi_data);
	gROOT->cd();

	TCanvas* canv=new TCanvas("reso_mass","mass distribution of mother particle;cross-section;mass (in MeV)",900,900);

	TH1F* mass=new TH1F("mother_py","mother particle",n_bins,0.,5000.);
	cout<<"number of entries: "<<phi_data->GetEntries()<<endl;
	int n=(int)((phi_data->GetEntries())/2);
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
	cout<<n<<endl;
	float* row_content;
	for(int i=0;i<phi_data->GetEntries();i++){
		//cout<<i<<endl;
		phi_data->GetEntry(i);
		row_content=phi_data->GetArgs();
		if(i%2==0) {
			
			for(int temp=0;temp<4;temp++){
				k_plus[i/2][temp]=row_content[temp];
			}
			
		}
		else{ 
			
			for(int temp=0;temp<4;temp++){
				k_minus[i/2][temp]=row_content[temp];
			}
			
		}
	}

	//delete[] row_content;
	cout<<"stored values"<<endl<<endl;
	
	for(int i=0;i<n;i++){
		//if (i%10000==0) cout<<i<<endl;
		px1=k_plus[i][0];
		py1=k_plus[i][1];
		pz1=k_plus[i][2];
		E1=k_plus[i][3];
		
		for(int j=0;j<n;j++){
			if(j%1000==0) cout<<i<<"\t"<<j<<endl;
			px2=k_minus[j][0];
			py2=k_minus[j][1];
			pz2=k_minus[j][2];
			E2=k_minus[j][3];
			E_m=E1+E2;
			px_m=px1+px2;
			py_m=py1+py2;
			pz_m=pz1+pz2;
			mass->Fill(sqrt(E_m*E_m - px_m*px_m - py_m*py_m - pz_m*pz_m));
			/*cout<<px1<<"\t "<<px2<<"\t "<<endl;
			cout<<py1<<" \t"<<py2<<"\t "<<endl;
			cout<<pz1<<" \t"<<pz2<<" \t"<<endl;
			cout<<E1<<"\t "<<E2<<"\t "<<endl<<endl;*/
		}
	}

	mass->Draw();
	/*canv->cd(0);*/canv->Modified();canv->Update();
	//canv->Modified();
	/*for(int i=0;i<n;i++){
		delete k_plus[i];
		delete k_minus[i];
	}
	delete[] k_plus;
	delete[] k_minus;*/
}
