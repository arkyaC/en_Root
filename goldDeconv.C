#define pi 3.141592653589
#define n_bins 1500
#define rangeMin 900.
#define rangeMax 1140.
#define massMin 0.
#define massMax 1500.
#define nIter 2000
#include <cmath>

void goldDeconv(){
	TFile in_file("conv_decay.root");
	TNtuple* conv_decay;
	TNtuple* event_data;
	in_file.GetObject("decay_data",conv_decay);
	in_file.GetObject("number_of_phi",event_data);
	gROOT->cd();

	TCanvas* canv=new TCanvas("inv_mass","Distribution of convolved invariant mass",900,1200);
	canv->Divide(1,2);
	canv->cd(1);
	gPad->SetTitle("actual data");
	canv->cd(2);
	gPad->SetTitle("Background removed");
//pre-processing
	TH1F* mass=new TH1F("inv_mass","Invariant Mass (convolved)",n_bins,massMin,massMax);
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
	//cout<<n<<endl;
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

	//cout<<"stored values"<<endl<<endl;
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
				//if(k%1000==0) cout<<j<<"\t"<<k<<endl;
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
//background removal
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
	TF1* g=new TF1("fitfn","[0]*exp(-0.5*(x-[1])*(x-[1])/([2]*[2]))",massMin,massMax); //fitted function -to be deconvolved
	g->SetParameter(0,par0);
	g->SetParameter(1,par1);
	g->SetParameter(2,par2);
	g->Draw("Same");
	gPad->Modified();gPad->Update();

//deconv starts
	float binWidth = (massMax - massMin)/n_bins;
	int N=(rangeMax-rangeMin)/binWidth - 100;
	float* observation= new float[2*N-1];
	float* noise= new float[N];
	float* invmass= new float[N];
	float** conv=new float*[2*N-1];
	float** convT=new float*[N];
//initialising obs and noise data points
	for(int i=0;i<N;i++){
		noise[i]=g->Eval(rangeMin+binWidth*(i+50));
		//cout<<rangeMin+binWidth*(i+50)<<endl;
		convT[i]=new float[2*N-1];
	}
	for(int i=0;i<2*N-1;i++){
		observation[i]=bg->GetBinContent(2*(bg->GetXaxis()->GetFirst()+49)+i-par1);//changed
		//cout<<bg->GetXaxis()->GetFirst()+49+i<<endl;
		conv[i]=new float[N];
	}
	/*cout<<noise[0]<<endl;
	cout<<observation[0]<<endl;
	cout<<observation[N/2]<<endl;
	cout<<observation[N]<<endl;
	cout<<observation[2*N-2]<<endl;*/
//initialising the convolution matrix conv
	for(int i=0;i<2*N-1;i++){
		if(i<N){
			int j=0;
			for(;j<i+1;j++)
				conv[i][j]=noise[i-j];
			for(;j<N;j++)
				conv[i][j]=0;
		}
		else{
			int j=0;
			for(;j<i-N+1;j++)
				conv[i][j]=0;
			for(;j<N;j++)
				conv[i][j]=noise[i-j];
		}
	}

//initialising conv transpose (convT)
	for(int i=0;i<N;i++){
		for(int j=0;j<2*N-1;j++){
			convT[i][j]=conv[j][i];
		}
	}

//C'*y=y1
	float* y1=new float[N];
	for(int i=0;i<N;i++){
		y1[i]=0;
		for(int j=0;j<(2*N-1);j++){
			y1[i]+=convT[i][j]*observation[j];
		}
	}

//C'*C=H
	float** H=new float*[N];
	for(int i=0;i<N;i++)
		H[i]=new float[N];
	for(int i=0;i<N;i++){
		for(int j=0;j<N;j++){
			H[i][j]=0;
			for(int k=0;k<2*N-1;k++)
				H[i][j]+=convT[i][k]*conv[k][j];
		}
	}

//H^2=C1
	float** C1=new float*[N];
	for(int i=0;i<N;i++)
		C1[i]=new float[N];
	for(int i=0;i<N;i++){
		for(int j=0;j<N;j++){
			C1[i][j]=0;
			for(int k=0;k<N;k++)
				C1[i][j]+=H[i][k]*H[k][j];
		}
	} 
//y2=H*y1, and set invmass initially to y2
	float* y2=new float[N];
	for(int i=0;i<N;i++){
		y2[i]=0;
		for(int j=0;j<N;j++){
			y2[i]+=H[i][j]*y1[j];
		}
		invmass[i]=y2[i];
	}
//gold deconvolution iterations
	float* invmassCopy=new float[N];
	for(int i=1;i<=nIter;i++){
		for(int j=0;j<N;j++){
			invmassCopy[j]=invmass[j];
		}
		for(int j=0;j<N;j++){
			float C1Xj=0;
			for(int k=0;k<N;k++)
				C1Xj+=C1[j][k]*invmassCopy[k];
			invmass[j]=invmass[j]*(y2[j]/C1Xj);
		}
	}
//drawing the histogram after deconvolution
	TH1F *deconv_mass = new TH1F("deconv_mass","deconvolved mass",N,rangeMin,rangeMax);
	
	for(int i=0;i<N;i++){
		deconv_mass->Fill(rangeMin+(i+50)*binWidth,max(invmass[i],1.e-4f));
		//cout<<i<<" "<<rangeMin+(i+50)*binWidth<<" "<<invmass[i]<<endl;
		//deconv_mass->Fill(rangeMin+i*binWidth,abs(invmass_real[i]));
	}
	TCanvas* canv1=new TCanvas("deconvolved_inv_mass","Distribution of deconvolved invariant mass",900,450);
	canv1->cd();
	deconv_mass->DrawNormalized();//don't know what option to set
	gPad->Modified();gPad->Update();
}