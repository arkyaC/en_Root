//deconvolution not working
#define pi 3.141592653589
#define n_bins 1000
#define rangeMin 942.
#define rangeMax 1101.
#define massMin 0.
#define massMax 1500.

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

	TH1F* mass=new TH1F("inv_mass","Invariant Mass (convoluted)",n_bins,massMin,massMax);
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
	mass->GetXaxis()->SetTitle("Invariant Mass (in MeV)");
	mass->GetYaxis()->SetTitle("Cross Section");
	gPad->Modified();gPad->Update();
	canv->cd(2);
	TH1* bg=mass->ShowBackground();
	cout<<bg->Add(mass,-1)<<endl;
	bg->Scale(-1.,"");
	bg->Draw("Same");
	bg->SetTitle("Background removed");
	bg->GetXaxis()->SetRangeUser(rangeMin,rangeMax);

	TFitResultPtr r= bg->Fit("gaus","S");
	double par0=r->Parameter(0);
	double par1=r->Parameter(1);
	double par2=r->Parameter(2);
	TF1* g=new TF1("fitfn","[0]*exp(-0.5*(x-[1])*(x-[1])/([2]*[2]))",rangeMin,rangeMax); //fitted function -to be deconvoluted
	g->SetParameter(0,par0);
	g->SetParameter(1,par1);
	g->SetParameter(2,par2);
	g->Draw("Same");
	gPad->Modified();gPad->Update();

	//deconvolution code (with FT)
	float binWidth = (massMax - massMin)/n_bins;
	int N=rangeMax-rangeMin;
	float* obsFT_real=new float[N];
	float* obsFT_imag=new float[N];
	float* noiseFT_real= new float[N];
	float* noiseFT_imag= new float[N];
	float* realmassFT_real= new float[N];
	float* realmassFT_imag= new float[N];
	float* observation= new float[N];
	float* noise= new float[N];
	float* realmass= new float[N];

	for(int i=0;i<N;i++){ //initialising obs and noise data points
		observation[i]=bg->GetBinContent(i+1+(int)(rangeMin/binWidth));
		noise[i]=g->Eval(rangeMin+binWidth*i);
	}
	for(int i=0;i<N;i++){ //evaluating DFT of obs
		obsFT_real[i]=0.;
		obsFT_imag[i]=0.;
		for(int j=0;j<N;j++){
			obsFT_real[i]+=observation[j] * cos(2*pi*i*j/N);
			obsFT_imag[i]+=-observation[j] * sin(2*pi*i*j/N);
		}
	}
	for(int i=0;i<N;i++){ //evaluating DFT of noise
		noiseFT_real[i]=0.;
		noiseFT_imag[i]=0.;
		for(int j=0;j<N;j++){
			noiseFT_real[i]+=noise[j] * cos(2*pi*i*j/N);
			noiseFT_imag[i]+=-noise[j] * sin(2*pi*i*j/N);
		}
	}
	for(int i=0;i<N;i++){
		realmassFT_real[i] = (obsFT_real[i]*noiseFT_real[i] + obsFT_imag[i]*noiseFT_imag[i])/(pow(noiseFT_real[i],2)+pow(noiseFT_imag[i],2));
		realmassFT_imag[i] = (obsFT_imag[i]*noiseFT_real[i] - obsFT_real[i]*noiseFT_imag[i])/(pow(noiseFT_real[i],2)+pow(noiseFT_imag[i],2));
	}
	//time for inverse DFT!!!...
	for(int i=0;i<N;i++){
		realmass[i]=0.;
		for(int j=0;j<N;j++){
			realmass[i]+=realmassFT_real[j] * cos(2*pi*j*i/N) - realmassFT_imag[j] * sin(2*pi*j*i/N);
		}
		realmass[i]=realmass[i]/N;
	}

	TH1F *deconv_mass = new TH1F("deconv_mass","deconvoluted mass",N,rangeMin,rangeMax);
	for(int i=0;i<N;i++){
		deconv_mass->Fill(rangeMin+i*binWidth,realmass[i]);
	}
	TCanvas* canv1=new TCanvas("deconvoluted_inv_mass","Distribution of deconvoluted invariant mass",900,450);
	canv1->cd();
	deconv_mass->Draw();
}
