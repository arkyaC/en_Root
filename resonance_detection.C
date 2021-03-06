
#define n_bins 300

void resonance_detection(){
	TFile in_file("decay_data_phi.root");
	TNtuple* phi_data;
	TNtuple* event_data;
	in_file.GetObject("decay_data",phi_data);
	in_file.GetObject("number_of_phi",event_data);
	gROOT->cd();

	TCanvas* canv=new TCanvas("inv_mass","Distribution of invariant mass",900,900);

	TH1F* mass=new TH1F("inv_mass","Invariant Mass",n_bins,0.98,1.35);
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

	cout<<"stored values"<<endl<<endl;
	int event_ctr=0;
	for(int i=0;i<event_data->GetEntries();i++){
		event_data->GetEntry(i);
		row_content=event_data->GetArgs();
		int n_events=(int)row_content[1];
		for(int j=event_ctr;j<event_ctr+n_events;j++){
			px1=k_plus[j][0];
			py1=k_plus[j][1];
			pz1=k_plus[j][2];
			E1=k_plus[j][3];
		
			for(int k=event_ctr;k<event_ctr+n_events;k++){
				if(k%1000==0) cout<<j<<"\t"<<k<<endl;
				px2=k_minus[k][0];
				py2=k_minus[k][1];
				pz2=k_minus[k][2];
				E2=k_minus[k][3];
				E_m=E1+E2;
				px_m=px1+px2;
				py_m=py1+py2;
				pz_m=pz1+pz2;
				mass->Fill(sqrt(E_m*E_m - px_m*px_m - py_m*py_m - pz_m*pz_m));
				//cout<<px1<<"\t "<<px2<<"\t "<<endl;
				//cout<<py1<<" \t"<<py2<<"\t "<<endl;
				//cout<<pz1<<" \t"<<pz2<<" \t"<<endl;
				//cout<<E1<<"\t "<<E2<<"\t "<<endl<<endl;
			}
		}
		event_ctr+=n_events;
	}

	mass->Draw();
	mass->GetXaxis()->SetTitle("Invariant Mass (in GeV)");
	mass->GetYaxis()->SetTitle("Cross Section");
	canv->Modified();canv->Update();
}
