
#define n_bins 300

void read_ConvDecay(){
	TFile in_file("conv_decay.root");
	TNtuple* conv_decay;
	TNtuple* event_data;
	in_file.GetObject("decay_data",conv_decay);
	in_file.GetObject("number_of_phi",event_data);
	gROOT->cd();

	TCanvas* canv=new TCanvas("inv_mass","Distribution of convoluted invariant mass",900,1200);
	canv->Divide(1,2);
	canv->cd(1);
	gPad->SetTitle("actual data");
	canv->cd(2);
	gPad->SetTitle("Background removed");

	TH1F* mass=new TH1F("inv_mass","Invariant Mass (convoluted)",n_bins,0.,1350.);
	cout<<"number of entries: "<<conv_decay->GetEntries()<<endl;
	int n=(int)((conv_decay->GetEntries())/2);
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
	for(int i=0;i<conv_decay->GetEntries();i++){
		//cout<<i<<endl;
		conv_decay->GetEntry(i);
		row_content=conv_decay->GetArgs();
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

	canv->cd(1);
	mass->Draw();
	mass->GetXaxis()->SetTitle("Invariant Mass (in GeV)");
	mass->GetYaxis()->SetTitle("Cross Section");
	gPad->Modified();gPad->Update();
	canv->cd(2);
	TH1* bg=mass->ShowBackground();
	cout<<bg->Add(mass,-1)<<endl;
	bg->Scale(-1.,"");
	bg->Draw();
	gPad->Modified();gPad->Update();
}
