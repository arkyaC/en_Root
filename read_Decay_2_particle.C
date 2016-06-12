//not working
#define n_bins 200

void init(TCanvas* canv){
	canv->Divide(1,3);
	canv->cd(1);
	gPad->SetTitle("particle 1");
	canv->cd(2);
	gPad->SetTitle("particle 2");
	canv->cd(3);
	gPad->SetTitle("mother");
}
void read_Decay_2_particle(){
	TFile in_file("decay_data_phi.root");
	TNtuple* phi_data;
	TNtuple* mother_data;
	in_file.GetObject("decay_data",phi_data);
	in_file.GetObject("mother_data",mother_data);
	gROOT->cd();

	TCanvas* canv1=new TCanvas("p_x","p_x",900,900);
	TCanvas* canv2=new TCanvas("p_y","p_y",900,900);
	TCanvas* canv3=new TCanvas("p_z","p_z",900,900);
	TCanvas* canv4=new TCanvas("E","E",900,900);
	init(canv1);
	init(canv2);
	init(canv3);
	init(canv4);
	
	TH1F* par1_E=new TH1F("particle_1_E","particle 1",n_bins,350.,1600.);
	TH1F* par2_E=new TH1F("particle_2_E","particle 2",n_bins,350.,1600.);
	TH1F* mother_E=new TH1F("mother_E","mother particle",n_bins,350.,1600.);
	TH1F* par1_px=new TH1F("particle_1_px","particle 1",n_bins,-700.,700.);
	TH1F* par2_px=new TH1F("particle_2_px","particle 2",n_bins,-700.,700.);
	TH1F* mother_px=new TH1F("mother_px","mother particle",n_bins,-800.,800.);
	TH1F* par1_py=new TH1F("particle_1_py","particle 1",n_bins,-700.,700.);
	TH1F* par2_py=new TH1F("particle_2_py","particle 2",n_bins,-700.,700.);
	TH1F* mother_py=new TH1F("mother_py","mother particle",n_bins,-800.,800.);
	TH1F* par1_pz=new TH1F("particle_1_pz","particle 1",n_bins,-700.,700.);
	TH1F* par2_pz=new TH1F("particle_2_pz","particle 2",n_bins,-700.,700.);
	TH1F* mother_pz=new TH1F("mother_pz","mother particle",n_bins,-2000.,2000.);
	
	float* row_content;
	float* mother_content;
	for(int i=0;i<phi_data->GetEntries();i++){
		phi_data->GetEntry(i);
		row_content=phi_data->GetArgs();
		if(i%2==0){
			//cout<<row_content[0]<<endl;
		    par1_px->Fill(row_content[0]);
		    par1_py->Fill(row_content[1]);
		    par1_pz->Fill(row_content[2]);
		    par1_E->Fill(row_content[3]);
		}
		else{
			par2_px->Fill(row_content[0]);
		    par2_py->Fill(row_content[1]);
		    par2_pz->Fill(row_content[2]);
		    par2_E->Fill(row_content[3]);
		}
	}
	for(int i=0;i<mother_data->GetEntries();i++){
		mother_data->GetEntry(i);
		mother_content=phi_data->GetArgs();
		mother_px->Fill(mother_content[0]);
		mother_py->Fill(mother_content[1]);
		mother_pz->Fill(mother_content[2]);
		mother_E->Fill(mother_content[3]);
	}

	canv1->cd(1);
	par1_px->Draw();
	canv1->cd(2);
	par2_px->Draw();
	canv1->cd(3);
	mother_px->Draw();
	canv2->cd(1);
	par1_py->Draw();
	canv2->cd(2);
	par2_py->Draw();
	canv2->cd(3);
	mother_py->Draw();
	canv3->cd(1);
	par1_pz->Draw();
	canv3->cd(2);
	par2_pz->Draw();
	canv3->cd(3);
	mother_pz->Draw();
	canv4->cd(1);
	par1_E->Draw();
	canv4->cd(2);
	par2_E->Draw();
	canv4->cd(3);
	mother_E->Draw();

	canv1->cd(0);canv1->Modified();canv1->Update();
	canv2->cd(0);canv2->Modified();canv2->Update();
	canv3->cd(0);canv3->Modified();canv3->Update();
	canv4->cd(0);canv4->Modified();canv4->Update();
}
