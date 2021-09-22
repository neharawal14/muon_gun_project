#include "TFile.h"
#include <stdio.h>
#include <iostream>
#include <vector>
#include "TChain.h"
#include "TTree.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TROOT.h"
#include "TMath.h"
#include "TStyle.h"
#include "TColor.h"
#include "TLegend.h"
#include "TPad.h"
#include "TLine.h"
#include "TH1F.h"
#include "TPaveStats.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooBreitWigner.h"
#include "RooCBShape.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooBernstein.h"
#include "RooExponential.h"
#include "RooFitResult.h"
#include "RooSimultaneous.h"
#include "RooCategory.h"
#include "RooAddPdf.h"
#include "RooPolynomial.h"
#include "RooVoigtian.h"
#include <fstream>
#include "RooLandau.h"
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

using namespace RooFit;

float Max(float a, float b);
void Draw(TH1F *h1, TH1F *h2, TH1F *h3, TString nome_canvas, TString save, TString x_axis, TLegend *legend, int fit, bool LogY = 0);
void Draw(TH1F *h1, TH1F *h2, TString nome_canvas, TString save, TString x_axis, TLegend *legend, int fit, bool LogY = 0);
void Draw(TH1F *h1, TH1F *h2, TH1F *h3, RooRealVar* rv_Mass_1, RooRealVar* rv_Mass_2, RooRealVar* rv_Mass_3, RooDataSet* Mass_1, RooDataSet* Mass_2, RooDataSet* Mass_3, TString nome_canvas, TString save, TString x_name, TLegend *legend, bool LogY);
std::vector<float> FitMass(TH1F *h1, TString nome_canvas, TString save);
std::vector<float> FitMass(RooRealVar* rv_Mass_1, RooDataSet* Mass_1, TString title, TString save);
std::vector<int> PickMuon(std::vector<float> pt_b, std::vector<float> eta_b, float pt, float eta);
void dxyBS(TString physics){

	gROOT->Reset();
	gROOT->SetBatch();
	
	TLegend *legend = new TLegend(0.75,0.75,0.9,0.9);
	TLegend *legend_noGen = new TLegend(0.75,0.75,0.9,0.9);
	TLegend *legend_reco = new TLegend(0.75,0.75,0.9,0.9);

	float BS_x, BS_y, BS_z;
	float BS_xErr, BS_yErr, BS_zErr, BeamWidth_x, BeamWidth_y, BeamWidth_xErr, BeamWidth_yErr;
	int nVtx;

	std::vector<float>* pt_GEN=0; std::vector<float>* eta_GEN=0; std::vector<float>* phi_GEN=0; std::vector<float>* mass_GEN=0;
	std::vector<float>* pt_FLO=0; std::vector<float>* eta_FLO=0; std::vector<float>* phi_FLO=0; std::vector<float>* mass_FLO=0; std::vector<float>* dxy_BS_FLO=0; std::vector<float>* dxy_PV_FLO=0; 
	std::vector<float>* pt_glb=0; std::vector<float>* eta_glb=0; std::vector<float>* phi_glb=0; std::vector<float>* mass_glb=0; std::vector<float>* dxy_BS_glb=0; std::vector<float>* dxy_PV_glb=0; 
	std::vector<float>* pt_stdAlone=0; std::vector<float>* eta_stdAlone=0; std::vector<float>* phi_stdAlone=0; std::vector<float>* dxy_BS_stdAlone=0; 
	std::vector<float>* pt_trk=0; std::vector<float>* eta_trk=0; std::vector<float>* phi_trk=0; std::vector<float>* dxy_BS_trk=0; 

	std::vector<float>* pt_FLO_vtx=0; std::vector<float>* eta_FLO_vtx=0; std::vector<float>* phi_FLO_vtx=0; std::vector<float>* dxy_BS_FLO_vtx=0; 
	std::vector<float>* pt_FLO_single=0; std::vector<float>* eta_FLO_single=0; std::vector<float>* phi_FLO_single=0; std::vector<float>* dxy_BS_FLO_single=0; 
	std::vector<float>* pt_glb_vtx=0; std::vector<float>* eta_glb_vtx=0; std::vector<float>* phi_glb_vtx=0; std::vector<float>* dxy_BS_glb_vtx=0; 
	std::vector<float>* pt_glb_single=0; std::vector<float>* eta_glb_single=0; std::vector<float>* phi_glb_single=0; std::vector<float>* dxy_BS_glb_single=0; 
	std::vector<float>* pt_stdAlone_vtx=0; std::vector<float>* eta_stdAlone_vtx=0; std::vector<float>* phi_stdAlone_vtx=0; std::vector<float>* dxy_BS_stdAlone_vtx=0; 
	std::vector<float>* pt_trk_vtx=0; std::vector<float>* eta_trk_vtx=0; std::vector<float>* phi_trk_vtx=0; std::vector<float>* dxy_BS_trk_vtx=0; 

// 	TString filename = "DUMMYFILENAME_2.root";
	TString filename = "DUMMYFILENAME_" + physics + ".root";
// 	TString filename = physics + ".root";
// 	TString filename = "root://cmsio5.rc.ufl.edu:1094///store/user/t2/users/ferrico/ZMuMu";// + ".root";
	TFile* _file0 = new TFile(filename );
	
	TTree* tree;
	
	if(_file0)
		tree = (TTree*)_file0->Get("passedEvents");		
	else std::cout<<"ERROR could not find the file"<<std::endl;


	std::vector<TString> lepton_type;
	lepton_type.push_back("GEN");
	lepton_type.push_back("FLO");
	lepton_type.push_back("glb");
	lepton_type.push_back("FLO_vtx");
	lepton_type.push_back("glb_vtx");
	lepton_type.push_back("FLO_single");
	lepton_type.push_back("glb_single");
// 	lepton_type.push_back("stdAlone");
// 	lepton_type.push_back("stdAlone_vtx");
	lepton_type.push_back("stdAlone");
	lepton_type.push_back("stdAlone_vtx");
	
	TH1F* lepton_pt[9];
	TH1F* lepton_deltaR[8];
	TH1F* lepton_d0BS[8];
	TH1F* lepton_d0PV[8];
	TH1F* lepton_ptRes[8];
	TH1F* dilepton_mass[9];


	std::vector<float> eta_bins;
	eta_bins.push_back(0);
	eta_bins.push_back(0.9);
	eta_bins.push_back(1.4);
	eta_bins.push_back(2.4);
	
	std::vector<float> pt_bins;
	pt_bins.push_back(5);
	pt_bins.push_back(20);
	pt_bins.push_back(30);
	pt_bins.push_back(40);
	pt_bins.push_back(50);
	pt_bins.push_back(60);
	pt_bins.push_back(100);
	pt_bins.push_back(200);
	
	
	std::vector<Double_t> vtx_bins;
	vtx_bins.push_back(0);
	vtx_bins.push_back(17);
	vtx_bins.push_back(20);
	vtx_bins.push_back(22);
	vtx_bins.push_back(24);
	vtx_bins.push_back(25);
	vtx_bins.push_back(26);
	vtx_bins.push_back(27);
	vtx_bins.push_back(28);
	vtx_bins.push_back(30);
	vtx_bins.push_back(32);
	vtx_bins.push_back(34);
	vtx_bins.push_back(36);
	vtx_bins.push_back(38);
	vtx_bins.push_back(40);
	vtx_bins.push_back(42);
	vtx_bins.push_back(44);
	vtx_bins.push_back(46);
	vtx_bins.push_back(48);
	vtx_bins.push_back(50);
	vtx_bins.push_back(60);
	
	RooRealVar* rv_ps[9][eta_bins.size()][pt_bins.size()];
	RooArgSet*  rastmp_ps[9][eta_bins.size()][pt_bins.size()];
	RooDataSet* Data_ps[9][eta_bins.size()][pt_bins.size()];

	RooRealVar* rv_nVtx[9][vtx_bins.size()];
	RooArgSet*  rastmp_nVtx[9][vtx_bins.size()];
	RooDataSet* Data_nVtx[9][vtx_bins.size()];
	
	TH1F* sigma_FLO = new TH1F("sigma_FLO", "sigma_FLO", 25, -0.5, 24.5);
	TH1F* sigma_glb = new TH1F("sigma_glb", "sigma_glb", 25, -0.5, 24.5);
	TH1F* sigma_FLO_vtx = new TH1F("sigma_FLO_vtx", "sigma_FLO_vtx", 25, -0.5, 24.5);
	TH1F* sigma_glb_vtx = new TH1F("sigma_glb_vtx", "sigma_glb_vtx", 25, -0.5, 24.5);
	TH1F* sigma_FLO_single = new TH1F("sigma_FLO_single", "sigma_FLO_single", 25, -0.5, 24.5);
	TH1F* sigma_glb_single = new TH1F("sigma_glb_single", "sigma_glb_single", 25, -0.5, 24.5);
	TH1F* sigma_stdAlone = new TH1F("sigma_stdAlone", "sigma_stdAlone", 25, -0.5, 24.5);
	TH1F* sigma_stdAlone_vtx = new TH1F("sigma_stdAlone_vtx", "sigma_stdAlone_vtx", 25, -0.5, 24.5);


	TH1F* sigma_FLO_nVtx = new TH1F("sigma_FLO nVtx", "sigma_FLO nVtx", vtx_bins.size()-1, &vtx_bins[0]);
	TH1F* sigma_glb_nVtx = new TH1F("sigma_glb nVtx", "sigma_glb nVtx", vtx_bins.size()-1, &vtx_bins[0]);
	TH1F* sigma_FLO_vtx_nVtx = new TH1F("sigma_FLO_vtx nVtx", "sigma_FLO_vtx nVtx", vtx_bins.size()-1, &vtx_bins[0]);
	TH1F* sigma_glb_vtx_nVtx = new TH1F("sigma_glb_vtx nVtx", "sigma_glb_vtx nVtx", vtx_bins.size()-1, &vtx_bins[0]);
	TH1F* sigma_FLO_single_nVtx = new TH1F("sigma_FLO_single nVtx", "sigma_FLO_single nVtx", vtx_bins.size()-1, &vtx_bins[0]);
	TH1F* sigma_glb_single_nVtx = new TH1F("sigma_glb_single nVtx", "sigma_glb_single nVtx", vtx_bins.size()-1, &vtx_bins[0]);

	
	RooRealVar* rv_Mass[9];
	RooArgSet*  rastmp_Mass[9];
	RooDataSet* Data_Mass[9];

	
	TH2F* BS_position = new TH2F("BS position", "BS position", 200, 0.0, 0.02, 200, 0.0, 0.05);
	TH2F* BS_position_error = new TH2F("BS position error", "BS position error", 200, -0.1, 0.1, 200, -0.1, 0.1);
	TH1F* BS_position_Z = new TH1F("BS position_Z", "BS position_Z", 200, 0.015, 0.035);
	TH1F* BS_position_error_Z = new TH1F("BS position error Z", "BS position error", 200, -0.1, 0.1);
	TH2F* BS_width = new TH2F("BS width", "BS width", 200, 0.0, 0.001, 200, 0.0, 0.001);
	
	for(int i = 0; i < lepton_type.size(); i++){

		TString name = "pt_" + lepton_type.at(i);
		lepton_pt[i] = new TH1F(name, name, 100, 0, 100);
		for(int	eta = 0; eta < eta_bins.size()-1; eta++){
			for(int pt = 0; pt < pt_bins.size()-1; pt++){
				name = Form("pt_%s_%.0f_%.0f_%.1f_%.1f", lepton_type[i].Data(), pt_bins.at(pt), pt_bins.at(pt+1),  eta_bins.at(eta), eta_bins.at(eta+1));
// 				std::cout<<name<<std::endl;
				rv_ps[i][eta][pt] = new RooRealVar(name, name, -0.25, 0.25);
				if(i > 6)
					rv_ps[i][eta][pt] = new RooRealVar(name, name, -1, 1);
				rastmp_ps[i][eta][pt] = new RooArgSet(*rv_ps[i][eta][pt]);
				Data_ps[i][eta][pt] = new RooDataSet("Data" + name, "Data" + name, *rastmp_ps[i][eta][pt]);
			}
		}

		name = "mass_" + lepton_type.at(i);
		dilepton_mass[i] = new TH1F(name, name, 100, 60, 120);		

		rv_Mass[i] = new RooRealVar(name, name, 105, 140);
		rastmp_Mass[i] = new RooArgSet(*rv_Mass[i]);
		Data_Mass[i] = new RooDataSet("Data" + name, "Data" + name, *rastmp_Mass[i]);

		if(i !=0){
			name = "deltaR_" + lepton_type.at(i);
			lepton_deltaR[i-1] = new TH1F(name, name, 100, 0, 0.005);		
// 			lepton_deltaR[i-1] = new TH1F(name, name, 100, 0, 0.1);		

			name = "d0BS" + lepton_type.at(i);
			lepton_d0BS[i-1] = new TH1F(name, name, 200, -0.01, 0.01);		
// 			lepton_d0BS[i-1] = new TH1F(name, name, 200, -0.05, 0.05);		

			name = "d0PV" + lepton_type.at(i);
			lepton_d0PV[i-1] = new TH1F(name, name, 200, -0.01, 0.01);		

			name = "pTRes" + lepton_type.at(i);
			lepton_ptRes[i-1] = new TH1F(name, name, 200, -0.1, 0.1);		
// 			lepton_ptRes[i-1] = new TH1F(name, name, 500, -0.025, 0.025);	

			for(int vtx = 1; vtx < vtx_bins.size(); vtx++){
				name = Form("RecoVtx_%s_%.0f_%.0f", lepton_type[i].Data(), vtx_bins.at(vtx-1), vtx_bins.at(vtx));	
				rv_nVtx[i-1][vtx-1] = new RooRealVar(name, name, -0.25, 0.25);
				rastmp_nVtx[i-1][vtx-1] = new RooArgSet(*rv_nVtx[i-1][vtx-1]);
				Data_nVtx[i-1][vtx-1] = new RooDataSet("Data" + name, "Data" + name, *rastmp_nVtx[i-1][vtx-1]);				
			}
	
		}
	
	}
	
	tree->SetBranchStatus("*",1);				

	tree->SetBranchAddress("BS_x", &BS_x);
	tree->SetBranchAddress("BS_y", &BS_y);
	tree->SetBranchAddress("BS_z", &BS_z);
	tree->SetBranchAddress("BS_xErr", &BS_xErr);
	tree->SetBranchAddress("BS_yErr", &BS_yErr);
	tree->SetBranchAddress("BS_zErr", &BS_zErr);
	tree->SetBranchAddress("BeamWidth_x", &BeamWidth_x);
	tree->SetBranchAddress("BeamWidth_y", &BeamWidth_y);
	tree->SetBranchAddress("BeamWidth_xErr", &BeamWidth_xErr);
	tree->SetBranchAddress("BeamWidth_yErr", &BeamWidth_yErr);
	tree->SetBranchAddress("nVtx", &nVtx);
				
	tree->SetBranchAddress("pt_GEN", &pt_GEN);
	tree->SetBranchAddress("eta_GEN", &eta_GEN);
	tree->SetBranchAddress("phi_GEN", &phi_GEN);
	tree->SetBranchAddress("mass_GEN", &mass_GEN);
	tree->SetBranchAddress("pt_FLO", &pt_FLO);
	tree->SetBranchAddress("eta_FLO", &eta_FLO);
	tree->SetBranchAddress("phi_FLO", &phi_FLO);
	tree->SetBranchAddress("mass_FLO", &mass_FLO);
	tree->SetBranchAddress("dxy_BS_FLO", &dxy_BS_FLO);
	tree->SetBranchAddress("dxy_PV_FLO", &dxy_PV_FLO);
	tree->SetBranchAddress("pt_glb", &pt_glb);
	tree->SetBranchAddress("eta_glb", &eta_glb);
	tree->SetBranchAddress("phi_glb", &phi_glb);
	tree->SetBranchAddress("mass_glb", &mass_glb);
	tree->SetBranchAddress("dxy_BS_glb", &dxy_BS_glb);
	tree->SetBranchAddress("dxy_PV_glb", &dxy_PV_glb);
	tree->SetBranchAddress("pt_trk", &pt_stdAlone);
	tree->SetBranchAddress("eta_trk", &eta_stdAlone);
	tree->SetBranchAddress("phi_trk", &phi_stdAlone);
	tree->SetBranchAddress("dxy_BS_trk", &dxy_BS_stdAlone);
// 	tree->SetBranchAddress("pt_stdAlone", &pt_stdAlone);
// 	tree->SetBranchAddress("eta_stdAlone", &eta_stdAlone);
// 	tree->SetBranchAddress("phi_stdAlone", &phi_stdAlone);
// 	tree->SetBranchAddress("dxy_BS_stdAlone", &dxy_BS_stdAlone);

	tree->SetBranchAddress("pt_FLO_vtx", &pt_FLO_vtx);
	tree->SetBranchAddress("eta_FLO_vtx", &eta_FLO_vtx);
	tree->SetBranchAddress("phi_FLO_vtx", &phi_FLO_vtx);
	tree->SetBranchAddress("dxy_BS_FLO_vtx", &dxy_BS_FLO_vtx);
	tree->SetBranchAddress("pt_FLO_single", &pt_FLO_single);
	tree->SetBranchAddress("eta_FLO_single", &eta_FLO_single);
	tree->SetBranchAddress("phi_FLO_single", &phi_FLO_single);
	tree->SetBranchAddress("dxy_BS_FLO_single", &dxy_BS_FLO_single);

	tree->SetBranchAddress("pt_glb_vtx", &pt_glb_vtx);
	tree->SetBranchAddress("eta_glb_vtx", &eta_glb_vtx);
	tree->SetBranchAddress("phi_glb_vtx", &phi_glb_vtx);
	tree->SetBranchAddress("dxy_BS_glb_vtx", &dxy_BS_glb_vtx);
	tree->SetBranchAddress("pt_glb_single", &pt_glb_single);
	tree->SetBranchAddress("eta_glb_single", &eta_glb_single);
	tree->SetBranchAddress("phi_glb_single", &phi_glb_single);
	tree->SetBranchAddress("dxy_BS_glb_single", &dxy_BS_glb_single);

	tree->SetBranchAddress("pt_trk_vtx", &pt_stdAlone_vtx);
	tree->SetBranchAddress("eta_trk_vtx", &eta_stdAlone_vtx);
	tree->SetBranchAddress("phi_trk_vtx", &phi_stdAlone_vtx);
	tree->SetBranchAddress("dxy_BS_trk_vtx", &dxy_BS_stdAlone_vtx);
// 	tree->SetBranchAddress("pt_stdAlone_vtx", &pt_stdAlone_vtx);
// 	tree->SetBranchAddress("eta_stdAlone_vtx", &eta_stdAlone_vtx);
// 	tree->SetBranchAddress("phi_stdAlone_vtx", &phi_stdAlone_vtx);
// 	tree->SetBranchAddress("dxy_BS_stdAlone_vtx", &dxy_BS_stdAlone_vtx);
	
	float min_dR;			
	float deltaEta;
	float deltaPhi;
	float DR2;
	int match;
	int kvtx;
	

	Long64_t nentries = tree->GetEntries();

	
	for(int entry = 0; entry < nentries; entry++){
// 	for(int entry = 0; entry < 10000; entry++){
	
		tree->GetEntry(entry);

		if(entry % 100000 == 0)       
			std::cout<<entry<<" su "<<nentries<<std::endl;      	
		
		BS_position->Fill(BS_x, BS_y);
		BS_position_error->Fill(BS_xErr, BS_yErr);
		BS_position_Z->Fill(BS_z);
		BS_position_error_Z->Fill(BS_zErr);
		BS_width->Fill(BeamWidth_x, BeamWidth_y);
		
		for(int vtx = 1; vtx < vtx_bins.size(); vtx++){
			if(nVtx < vtx_bins.at(vtx)){
				kvtx = vtx-1;
				break;
			}
		}

		TLorentzVector mu1, mu2, Zboson;

		for(int gen = 0; gen < pt_GEN->size(); gen++){
			lepton_pt[0]->Fill(pt_GEN->at(gen));
			if(gen == 0)			
				mu1.SetPtEtaPhiM(pt_GEN->at(gen), eta_GEN->at(gen), phi_GEN->at(gen), mass_GEN->at(gen));
			else
				mu2.SetPtEtaPhiM(pt_GEN->at(gen), eta_GEN->at(gen), phi_GEN->at(gen), mass_GEN->at(gen));
		}
		Zboson = mu1+mu2;
		dilepton_mass[0]->Fill(Zboson.M());
		rv_Mass[0]->setVal(Zboson.M());
		Data_Mass[0]->add(*rastmp_Mass[0]);
		
		std::vector<int> pick_muon;

		for(int flo = 0; flo < pt_FLO->size(); flo++){
			lepton_pt[1]->Fill(pt_FLO->at(flo));
			lepton_d0BS[0]->Fill(dxy_BS_FLO->at(flo));
// 			lepton_d0PV[0]->Fill(dxy_PV_FLO->at(flo));

			min_dR = 999;

			for(int gen = 0; gen < pt_GEN->size(); gen++){
				deltaEta = eta_FLO->at(flo) - eta_GEN->at(gen);
				deltaPhi = phi_FLO->at(flo) - phi_GEN->at(gen);
				DR2 = sqrt(pow(deltaEta, 2) + pow(deltaPhi, 2));	
				if(DR2 < min_dR){
					min_dR = DR2;
					if(gen == 1) match = 1;
					else match = 0;
				}				
			}	
			
			
			if(pt_FLO->size() == 2){
				if(flo == 0)
					mu1.SetPtEtaPhiM(pt_FLO->at(flo), eta_FLO->at(flo), phi_FLO->at(flo), mass_FLO->at(flo));
				else{
					mu2.SetPtEtaPhiM(pt_FLO->at(flo), eta_FLO->at(flo), phi_FLO->at(flo), mass_FLO->at(flo));
					Zboson = mu1+mu2;
					dilepton_mass[1]->Fill(Zboson.M());
					rv_Mass[1]->setVal(Zboson.M());
					Data_Mass[1]->add(*rastmp_Mass[1]);
				}
			}
					
			lepton_deltaR[0]->Fill(min_dR);
			lepton_ptRes[0]->Fill((pt_FLO->at(flo) - pt_GEN->at(match))/pt_GEN->at(match));
			pick_muon = PickMuon(pt_bins, eta_bins, pt_FLO->at(flo), fabs(eta_FLO->at(flo)));
			if(pick_muon.size() > 1){
// 				std::cout<<pick_muon.at(0)<<"\t"<<pt_FLO->at(flo)<<"\t"<<pick_muon.at(1)<<"\t"<<fabs(eta_FLO->at(flo))<<std::endl;
				rv_ps[1][pick_muon.at(0)-1][pick_muon.at(1)-1]->setVal((pt_FLO->at(flo) - pt_GEN->at(match))/pt_GEN->at(match));
				Data_ps[1][pick_muon.at(0)-1][pick_muon.at(1)-1]->add(*rastmp_ps[1][pick_muon.at(0)-1][pick_muon.at(1)-1]);

				rv_nVtx[0][kvtx]->setVal((pt_FLO->at(flo) - pt_GEN->at(match))/pt_GEN->at(match));
				Data_nVtx[0][kvtx]->add(*rastmp_nVtx[1][kvtx]);

			}			
		}

		for(int glb = 0; glb < pt_glb->size(); glb++){
			lepton_pt[2]->Fill(pt_glb->at(glb));
			lepton_d0BS[1]->Fill(dxy_BS_glb->at(glb));
// 			lepton_d0PV[1]->Fill(dxy_PV_glb->at(glb));
	
			min_dR = 999;

			for(int gen = 0; gen < pt_GEN->size(); gen++){

				deltaEta = eta_glb->at(glb) - eta_GEN->at(gen);
				deltaPhi = phi_glb->at(glb) - phi_GEN->at(gen);
				DR2 = sqrt(pow(deltaEta, 2) + pow(deltaPhi, 2));	
				if(DR2 < min_dR){
					min_dR = DR2;
					if(gen == 1) match = 1;
					else match = 0;
				}

			}
			
			if(pt_glb->size() == 2){
				if(glb == 0)
					mu1.SetPtEtaPhiM(pt_glb->at(glb), eta_glb->at(glb), phi_glb->at(glb), mass_glb->at(glb));
				else{
					mu2.SetPtEtaPhiM(pt_glb->at(glb), eta_glb->at(glb), phi_glb->at(glb), mass_glb->at(glb));
					Zboson = mu1+mu2;
					dilepton_mass[2]->Fill(Zboson.M());
					rv_Mass[2]->setVal(Zboson.M());
					Data_Mass[2]->add(*rastmp_Mass[2]);
				}
			}
			
			lepton_deltaR[1]->Fill(min_dR);
			lepton_ptRes[1]->Fill((pt_glb->at(glb) - pt_GEN->at(match))/pt_GEN->at(match));
			pick_muon = PickMuon(pt_bins, eta_bins, pt_glb->at(glb), fabs(eta_glb->at(glb)));
			if(pick_muon.size() > 1){
// 				std::cout<<pick_muon.at(0)<<"\t"<<pt_FLO->at(glb)<<"\t"<<pick_muon.at(1)<<"\t"<<fabs(eta_FLO->at(glb))<<std::endl;
				rv_ps[2][pick_muon.at(0)-1][pick_muon.at(1)-1]->setVal((pt_glb->at(glb) - pt_GEN->at(match))/pt_GEN->at(match));
				Data_ps[2][pick_muon.at(0)-1][pick_muon.at(1)-1]->add(*rastmp_ps[2][pick_muon.at(0)-1][pick_muon.at(1)-1]);

				rv_nVtx[1][kvtx]->setVal((pt_glb->at(glb) - pt_GEN->at(match))/pt_GEN->at(match));
				Data_nVtx[1][kvtx]->add(*rastmp_nVtx[1][kvtx]);
			}			
		}

		for(int flo = 0; flo < pt_FLO_vtx->size(); flo++){
			lepton_pt[3]->Fill(pt_FLO_vtx->at(flo));
			lepton_d0BS[2]->Fill(dxy_BS_FLO_vtx->at(flo));

			min_dR = 999;

			for(int gen = 0; gen < pt_GEN->size(); gen++){
				deltaEta = eta_FLO_vtx->at(flo) - eta_GEN->at(gen);
				deltaPhi = phi_FLO_vtx->at(flo) - phi_GEN->at(gen);
				DR2 = sqrt(pow(deltaEta, 2) + pow(deltaPhi, 2));	
				if(DR2 < min_dR){
					min_dR = DR2;
					if(gen == 1) match = 1;
					else match = 0;
				}				
			}	
			
			if(pt_FLO_vtx->size() == 2){
				if(flo == 0)
					mu1.SetPtEtaPhiM(pt_FLO_vtx->at(flo), eta_FLO_vtx->at(flo), phi_FLO_vtx->at(flo), mass_FLO->at(flo));
				else{
					mu2.SetPtEtaPhiM(pt_FLO_vtx->at(flo), eta_FLO_vtx->at(flo), phi_FLO_vtx->at(flo), mass_FLO->at(flo));
					Zboson = mu1+mu2;
					dilepton_mass[3]->Fill(Zboson.M());
					rv_Mass[3]->setVal(Zboson.M());
					Data_Mass[3]->add(*rastmp_Mass[3]);
				}
			}		
			
			lepton_deltaR[2]->Fill(min_dR);
			lepton_ptRes[2]->Fill((pt_FLO_vtx->at(flo) - pt_GEN->at(match))/pt_GEN->at(match));
			pick_muon = PickMuon(pt_bins, eta_bins, pt_FLO_vtx->at(flo), fabs(eta_FLO_vtx->at(flo)));
			if(pick_muon.size() > 1){
// 				std::cout<<pick_muon.at(0)<<"\t"<<pt_FLO->at(flo)<<"\t"<<pick_muon.at(1)<<"\t"<<fabs(eta_FLO->at(flo))<<std::endl;
				rv_ps[3][pick_muon.at(0)-1][pick_muon.at(1)-1]->setVal((pt_FLO_vtx->at(flo) - pt_GEN->at(match))/pt_GEN->at(match));
				Data_ps[3][pick_muon.at(0)-1][pick_muon.at(1)-1]->add(*rastmp_ps[3][pick_muon.at(0)-1][pick_muon.at(1)-1]);

				rv_nVtx[2][kvtx]->setVal((pt_FLO_vtx->at(flo) - pt_GEN->at(match))/pt_GEN->at(match));
				Data_nVtx[2][kvtx]->add(*rastmp_nVtx[1][kvtx]);
			}			
		}

		for(int glb_vtx = 0; glb_vtx < pt_glb_vtx->size(); glb_vtx++){
			lepton_pt[4]->Fill(pt_glb_vtx->at(glb_vtx));
			lepton_d0BS[3]->Fill(dxy_BS_glb_vtx->at(glb_vtx));
	
			min_dR = 999;

			for(int gen = 0; gen < pt_GEN->size(); gen++){

				deltaEta = eta_glb_vtx->at(glb_vtx) - eta_GEN->at(gen);
				deltaPhi = phi_glb_vtx->at(glb_vtx) - phi_GEN->at(gen);
				DR2 = sqrt(pow(deltaEta, 2) + pow(deltaPhi, 2));	
				if(DR2 < min_dR){
					min_dR = DR2;
					if(gen == 1) match = 1;
					else match = 0;
				}

			}
			
			if(pt_glb_vtx->size() == 2){
				if(glb_vtx == 0)
					mu1.SetPtEtaPhiM(pt_glb_vtx->at(glb_vtx), eta_glb_vtx->at(glb_vtx), phi_glb_vtx->at(glb_vtx), mass_glb->at(glb_vtx));
				else{
					mu2.SetPtEtaPhiM(pt_glb_vtx->at(glb_vtx), eta_glb_vtx->at(glb_vtx), phi_glb_vtx->at(glb_vtx), mass_glb->at(glb_vtx));
					Zboson = mu1+mu2;
					dilepton_mass[4]->Fill(Zboson.M());
					rv_Mass[4]->setVal(Zboson.M());
					Data_Mass[4]->add(*rastmp_Mass[4]);
				}
			}
			
			lepton_deltaR[3]->Fill(min_dR);
			lepton_ptRes[3]->Fill((pt_glb_vtx->at(glb_vtx) - pt_GEN->at(match))/pt_GEN->at(match));
			pick_muon = PickMuon(pt_bins, eta_bins, pt_glb_vtx->at(glb_vtx), fabs(eta_glb_vtx->at(glb_vtx)));
			if(pick_muon.size() > 1){
// 				std::cout<<pick_muon.at(0)<<"\t"<<pt_FLO->at(glb_vtx)<<"\t"<<pick_muon.at(1)<<"\t"<<fabs(eta_FLO->at(glb_vtx))<<std::endl;
				rv_ps[4][pick_muon.at(0)-1][pick_muon.at(1)-1]->setVal((pt_glb_vtx->at(glb_vtx) - pt_GEN->at(match))/pt_GEN->at(match));
				Data_ps[4][pick_muon.at(0)-1][pick_muon.at(1)-1]->add(*rastmp_ps[4][pick_muon.at(0)-1][pick_muon.at(1)-1]);

				rv_nVtx[3][kvtx]->setVal((pt_glb_vtx->at(glb_vtx) - pt_GEN->at(match))/pt_GEN->at(match));
				Data_nVtx[3][kvtx]->add(*rastmp_nVtx[1][kvtx]);
			}			
		}

		for(int flo = 0; flo < pt_FLO_single->size(); flo++){
			lepton_pt[5]->Fill(pt_FLO_single->at(flo));
			lepton_d0BS[4]->Fill(dxy_BS_FLO_single->at(flo));

			min_dR = 999;

			for(int gen = 0; gen < pt_GEN->size(); gen++){
				deltaEta = eta_FLO_single->at(flo) - eta_GEN->at(gen);
				deltaPhi = phi_FLO_single->at(flo) - phi_GEN->at(gen);
				DR2 = sqrt(pow(deltaEta, 2) + pow(deltaPhi, 2));	
				if(DR2 < min_dR){
					min_dR = DR2;
					if(gen == 1) match = 1;
					else match = 0;
				}				
			}	
			
			if(pt_FLO_single->size() == 2){
				if(flo == 0)
					mu1.SetPtEtaPhiM(pt_FLO_single->at(flo), eta_FLO_single->at(flo), phi_FLO_single->at(flo), mass_FLO->at(flo));
				else{
					mu2.SetPtEtaPhiM(pt_FLO_single->at(flo), eta_FLO_single->at(flo), phi_FLO_single->at(flo), mass_FLO->at(flo));
					Zboson = mu1+mu2;
					dilepton_mass[5]->Fill(Zboson.M());
					rv_Mass[5]->setVal(Zboson.M());
					Data_Mass[5]->add(*rastmp_Mass[5]);
				}
			}
		
			
			lepton_deltaR[4]->Fill(min_dR);
			lepton_ptRes[4]->Fill((pt_FLO_single->at(flo) - pt_GEN->at(match))/pt_GEN->at(match));
			pick_muon = PickMuon(pt_bins, eta_bins, pt_FLO_single->at(flo), fabs(eta_FLO_single->at(flo)));
			if(pick_muon.size() > 1){
// 				std::cout<<pick_muon.at(0)<<"\t"<<pt_FLO->at(flo)<<"\t"<<pick_muon.at(1)<<"\t"<<fabs(eta_FLO->at(flo))<<std::endl;
				rv_ps[5][pick_muon.at(0)-1][pick_muon.at(1)-1]->setVal((pt_FLO_single->at(flo) - pt_GEN->at(match))/pt_GEN->at(match));
				Data_ps[5][pick_muon.at(0)-1][pick_muon.at(1)-1]->add(*rastmp_ps[5][pick_muon.at(0)-1][pick_muon.at(1)-1]);

				rv_nVtx[4][kvtx]->setVal((pt_FLO_single->at(flo) - pt_GEN->at(match))/pt_GEN->at(match));
				Data_nVtx[4][kvtx]->add(*rastmp_nVtx[1][kvtx]);
			}			
		}

		for(int glb_single = 0; glb_single < pt_glb_single->size(); glb_single++){
			lepton_pt[6]->Fill(pt_glb_single->at(glb_single));
			lepton_d0BS[5]->Fill(dxy_BS_glb_single->at(glb_single));
	
			min_dR = 999;

			for(int gen = 0; gen < pt_GEN->size(); gen++){

				deltaEta = eta_glb_single->at(glb_single) - eta_GEN->at(gen);
				deltaPhi = phi_glb_single->at(glb_single) - phi_GEN->at(gen);
				DR2 = sqrt(pow(deltaEta, 2) + pow(deltaPhi, 2));	
				if(DR2 < min_dR){
					min_dR = DR2;
					if(gen == 1) match = 1;
					else match = 0;
				}

			}
			
			if(pt_glb_single->size() == 2){
				if(glb_single == 0)
					mu1.SetPtEtaPhiM(pt_glb_single->at(glb_single), eta_glb_single->at(glb_single), phi_glb_single->at(glb_single), mass_glb->at(glb_single));
				else{
					mu2.SetPtEtaPhiM(pt_glb_single->at(glb_single), eta_glb_single->at(glb_single), phi_glb_single->at(glb_single), mass_glb->at(glb_single));
					Zboson = mu1+mu2;
					dilepton_mass[6]->Fill(Zboson.M());
					rv_Mass[6]->setVal(Zboson.M());
					Data_Mass[6]->add(*rastmp_Mass[6]);
				}
			}

			
			lepton_deltaR[5]->Fill(min_dR);
			lepton_ptRes[5]->Fill((pt_glb_single->at(glb_single) - pt_GEN->at(match))/pt_GEN->at(match));
			pick_muon = PickMuon(pt_bins, eta_bins, pt_glb_single->at(glb_single), fabs(eta_glb_single->at(glb_single)));
			if(pick_muon.size() > 1){
// 				std::cout<<pick_muon.at(0)<<"\t"<<pt_FLO->at(glb_single)<<"\t"<<pick_muon.at(1)<<"\t"<<fabs(eta_FLO->at(glb_single))<<std::endl;
				rv_ps[6][pick_muon.at(0)-1][pick_muon.at(1)-1]->setVal((pt_glb_single->at(glb_single) - pt_GEN->at(match))/pt_GEN->at(match));
				Data_ps[6][pick_muon.at(0)-1][pick_muon.at(1)-1]->add(*rastmp_ps[6][pick_muon.at(0)-1][pick_muon.at(1)-1]);

				rv_nVtx[5][kvtx]->setVal((pt_glb_single->at(glb_single) - pt_GEN->at(match))/pt_GEN->at(match));
				Data_nVtx[5][kvtx]->add(*rastmp_nVtx[1][kvtx]);
			}			
		}
		
		for(int stdAlone = 0; stdAlone < pt_stdAlone->size(); stdAlone++){
			lepton_pt[7]->Fill(pt_stdAlone->at(stdAlone));
			lepton_d0BS[6]->Fill(dxy_BS_stdAlone->at(stdAlone));
// 			lepton_d0PV[1]->Fill(dxy_PV_stdAlone->at(stdAlone));
	
			min_dR = 999;

			for(int gen = 0; gen < pt_GEN->size(); gen++){

				deltaEta = eta_stdAlone->at(stdAlone) - eta_GEN->at(gen);
				deltaPhi = phi_stdAlone->at(stdAlone) - phi_GEN->at(gen);
				DR2 = sqrt(pow(deltaEta, 2) + pow(deltaPhi, 2));	
				if(DR2 < min_dR){
					min_dR = DR2;
					if(gen == 1) match = 1;
					else match = 0;
				}

			}
			
			if(pt_stdAlone->size() == 2){
				if(stdAlone == 0)
					mu1.SetPtEtaPhiM(pt_stdAlone->at(stdAlone), eta_stdAlone->at(stdAlone), phi_stdAlone->at(stdAlone), 0.10565837);
				else{
					mu2.SetPtEtaPhiM(pt_stdAlone->at(stdAlone), eta_stdAlone->at(stdAlone), phi_stdAlone->at(stdAlone), 0.10565837);
					Zboson = mu1+mu2;
					dilepton_mass[7]->Fill(Zboson.M());
					rv_Mass[7]->setVal(Zboson.M());
					Data_Mass[7]->add(*rastmp_Mass[7]);
				}
			}
			
			lepton_deltaR[6]->Fill(min_dR);
			lepton_ptRes[6]->Fill((pt_stdAlone->at(stdAlone) - pt_GEN->at(match))/pt_GEN->at(match));
			pick_muon = PickMuon(pt_bins, eta_bins, pt_stdAlone->at(stdAlone), fabs(eta_stdAlone->at(stdAlone)));
			if(pick_muon.size() > 1){
// 				std::cout<<"\t\t\tCIOLA\t"<<(pt_stdAlone->at(stdAlone) - pt_GEN->at(match))/pt_GEN->at(match)<<std::endl;
// 				std::cout<<pt_stdAlone->at(stdAlone)<<"\t"<<fabs(eta_stdAlone->at(stdAlone))<<"\t"<<pt_GEN->at(match)<<std::endl;
// 				std::cout<<pick_muon.at(0)<<"\t"<<pt_FLO->at(stdAlone)<<"\t"<<pick_muon.at(1)<<"\t"<<fabs(eta_FLO->at(stdAlone))<<std::endl;
				rv_ps[7][pick_muon.at(0)-1][pick_muon.at(1)-1]->setVal((pt_stdAlone->at(stdAlone) - pt_GEN->at(match))/pt_GEN->at(match));
				Data_ps[7][pick_muon.at(0)-1][pick_muon.at(1)-1]->add(*rastmp_ps[7][pick_muon.at(0)-1][pick_muon.at(1)-1]);
			}	
// 			else{
// 				std::cout<<pt_stdAlone->at(stdAlone)<<"\t"<<fabs(eta_stdAlone->at(stdAlone))<<std::endl;
// 				std::cout<<"CAZZO"<<std::endl;	
// 			}	

			rv_nVtx[6][kvtx]->setVal((pt_stdAlone->at(stdAlone) - pt_GEN->at(match))/pt_GEN->at(match));
			Data_nVtx[6][kvtx]->add(*rastmp_nVtx[1][kvtx]);
		}	
		
		for(int stdAlone_vtx = 0; stdAlone_vtx < pt_stdAlone_vtx->size(); stdAlone_vtx++){
			lepton_pt[8]->Fill(pt_stdAlone_vtx->at(stdAlone_vtx));
			lepton_d0BS[7]->Fill(dxy_BS_stdAlone_vtx->at(stdAlone_vtx));
	
			min_dR = 999;

			for(int gen = 0; gen < pt_GEN->size(); gen++){

				deltaEta = eta_stdAlone_vtx->at(stdAlone_vtx) - eta_GEN->at(gen);
				deltaPhi = phi_stdAlone_vtx->at(stdAlone_vtx) - phi_GEN->at(gen);
				DR2 = sqrt(pow(deltaEta, 2) + pow(deltaPhi, 2));	
				if(DR2 < min_dR){
					min_dR = DR2;
					if(gen == 1) match = 1;
					else match = 0;
				}

			}
			
			if(pt_stdAlone_vtx->size() == 2){
				if(stdAlone_vtx == 0)
					mu1.SetPtEtaPhiM(pt_stdAlone_vtx->at(stdAlone_vtx), eta_stdAlone_vtx->at(stdAlone_vtx), phi_stdAlone_vtx->at(stdAlone_vtx), 0.10565837);
				else{
					mu2.SetPtEtaPhiM(pt_stdAlone_vtx->at(stdAlone_vtx), eta_stdAlone_vtx->at(stdAlone_vtx), phi_stdAlone_vtx->at(stdAlone_vtx), 0.10565837);
					Zboson = mu1+mu2;
					dilepton_mass[8]->Fill(Zboson.M());
					rv_Mass[8]->setVal(Zboson.M());
					Data_Mass[8]->add(*rastmp_Mass[8]);
				}
			}
			
			lepton_deltaR[7]->Fill(min_dR);
			lepton_ptRes[7]->Fill((pt_stdAlone_vtx->at(stdAlone_vtx) - pt_GEN->at(match))/pt_GEN->at(match));
			if(pick_muon.size() > 1){
// 				std::cout<<pick_muon.at(0)<<"\t"<<pt_FLO->at(stdAlone_vtx)<<"\t"<<pick_muon.at(1)<<"\t"<<fabs(eta_FLO->at(stdAlone_vtx))<<std::endl;
				rv_ps[8][pick_muon.at(0)-1][pick_muon.at(1)-1]->setVal((pt_stdAlone_vtx->at(stdAlone_vtx) - pt_GEN->at(match))/pt_GEN->at(match));
				Data_ps[8][pick_muon.at(0)-1][pick_muon.at(1)-1]->add(*rastmp_ps[8][pick_muon.at(0)-1][pick_muon.at(1)-1]);
			}			

			rv_nVtx[7][kvtx]->setVal((pt_stdAlone_vtx->at(stdAlone_vtx) - pt_GEN->at(match))/pt_GEN->at(match));
			Data_nVtx[7][kvtx]->add(*rastmp_nVtx[1][kvtx]);
		}		
	
	}
		
	legend->AddEntry(lepton_pt[0], "GEN");
	lepton_pt[1]->SetLineColor(kRed+1);
	legend->AddEntry(lepton_pt[1], "FLO");
	lepton_pt[2]->SetLineColor(kGreen+1);
	legend->AddEntry(lepton_pt[2], "glb");

	legend_noGen->AddEntry(lepton_pt[1], "FLO");
	legend_noGen->AddEntry(lepton_pt[2], "glb");

	legend_reco->AddEntry(lepton_pt[1], "base");
	legend_reco->AddEntry(lepton_pt[3], "VX+BS");
	legend_reco->AddEntry(lepton_pt[5], "Single");
	
	gStyle->SetOptFit(1);	


	TCanvas* c1 = new TCanvas("pT", "pT", 750, 500);
	c1->Update();
	BS_position->Draw("COLZ");
	c1->Print(physics + ".pdf[");
	c1->Print(physics + ".pdf");

	BS_position_error->Draw("COLZ");
	c1->Print(physics + ".pdf");

	BS_position_Z->Draw();
	c1->Print(physics + ".pdf");

	BS_position_error_Z->Draw();
	c1->Print(physics + ".pdf");

	BS_width->Draw("COLZ");
	c1->Print(physics + ".pdf");
	
// 	Draw(dilepton_mass[0], dilepton_mass[1], dilepton_mass[0], "Mass", physics + ".pdf[", "m_{ll} [GeV]", legend, 0);
	Draw(dilepton_mass[0], dilepton_mass[1], dilepton_mass[2], "Mass", physics + ".pdf", "m_{ll} [GeV]", legend, 0, 0);
// 	Draw(dilepton_mass[1], dilepton_mass[3], dilepton_mass[5], rv_Mass[1], rv_Mass[3], rv_Mass[5], Data_Mass[1], Data_Mass[3], Data_Mass[5], "Mass", physics + ".pdf", "m_{ll} [GeV]", legend_reco, 0); 
// 	Draw(dilepton_mass[2], dilepton_mass[4], dilepton_mass[6], rv_Mass[2], rv_Mass[4], rv_Mass[6], Data_Mass[2], Data_Mass[4], Data_Mass[6], "Mass", physics + ".pdf", "m_{ll} [GeV]", legend_reco, 0); 
	Draw(dilepton_mass[1], dilepton_mass[3], dilepton_mass[5], "Mass", physics + ".pdf", "m_{ll} [GeV]", legend_reco, 2, 0);
	Draw(dilepton_mass[2], dilepton_mass[4], dilepton_mass[6], "Mass", physics + ".pdf", "m_{ll} [GeV]", legend_reco, 2, 0);

	Draw(lepton_deltaR[0], lepton_deltaR[1], "#DeltaR", physics + ".pdf", "#DeltaR", legend_noGen, 0, 0);
	Draw(lepton_deltaR[0], lepton_deltaR[2], lepton_deltaR[4], "#DeltaR", physics + ".pdf", "#DeltaR", legend_reco, 0, 0);
	Draw(lepton_deltaR[1], lepton_deltaR[3], lepton_deltaR[5], "#DeltaR", physics + ".pdf", "#DeltaR", legend_reco, 0, 0);

	Draw(lepton_ptRes[0], lepton_ptRes[1], "#Deltap_{T}/p_{T}", physics + ".pdf", "#Deltap_{T}/p_{T}", legend_noGen, 1, 0);

	Draw(lepton_ptRes[0], lepton_ptRes[2], "#Deltap_{T}/p_{T}", physics + ".pdf", "#Deltap_{T}/p_{T}", legend_noGen, 1, 0);
	Draw(lepton_ptRes[1], lepton_ptRes[3], "#Deltap_{T}/p_{T}", physics + ".pdf", "#Deltap_{T}/p_{T}", legend_noGen, 1, 0);

	Draw(lepton_ptRes[6], lepton_ptRes[7], "#Deltap_{T}/p_{T}", physics + ".pdf", "#Deltap_{T}/p_{T}", legend_noGen, 1, 0);

	Draw(lepton_ptRes[0], lepton_ptRes[2], lepton_ptRes[4], "#Deltap_{T}/p_{T}", physics + ".pdf", "#Deltap_{T}/p_{T}", legend_reco, 1, 0);
	Draw(lepton_ptRes[1], lepton_ptRes[3], lepton_ptRes[5], "#Deltap_{T}/p_{T}", physics + ".pdf", "#Deltap_{T}/p_{T}", legend_reco, 1, 0);

	Draw(lepton_d0BS[0], lepton_d0BS[1], "d0 (wrt BS)", physics + ".pdf", "d0 (cm)", legend_noGen, 1, 0);

	Draw(lepton_d0BS[0], lepton_d0BS[2], "d0 (wrt BS)", physics + ".pdf", "d0 (cm)", legend_noGen, 1, 0);
	Draw(lepton_d0BS[1], lepton_d0BS[3], "d0 (wrt BS)", physics + ".pdf", "d0 (cm)", legend_noGen, 1, 0);

	Draw(lepton_d0PV[0], lepton_d0PV[1], "d0 (wrt PV)" , physics + ".pdf", "d0 (cm)", legend_noGen, 1, 0);
	Draw(lepton_d0PV[0], lepton_d0PV[1], "d0 (wrt PV)", physics + ".pdf]", "d0 (cm)", legend_noGen, 1, 0);


	for(int i = 1; i < lepton_type.size(); i++){
// 	for(int i = 1; i < 2; i++){
		int bin = 1;
		for(int	eta = 0; eta < eta_bins.size()-1; eta++){
			for(int pt = 0; pt < pt_bins.size()-1; pt++){
				TString name = Form("pt_%s_%.0f_%.0f_%.1f_%.1f", lepton_type[i].Data(), pt_bins.at(pt), pt_bins.at(pt+1),  eta_bins.at(eta), eta_bins.at(eta+1));
				RooDataSet h_Higgs = RooDataSet(Data_ps[i][eta][pt]->GetName(), Data_ps[i][eta][pt]->GetTitle(), Data_ps[i][eta][pt], *Data_ps[i][eta][pt]->get());

				RooRealVar Mean("Mean", "Mean", 0, -0.05, 0.05);
				RooRealVar Sigma("Sigma", "Sigma", 0.005, 0.001, 0.1);

				RooGaussian Gauss("Gauss", "Gauss", *rv_ps[i][eta][pt], Mean, Sigma);

				TCanvas *c_MC = new TCanvas(name, name, 750, 500);
				c_MC->SetFrameFillColor(0);
				c_MC->cd(1)->SetBottomMargin(0.2);
				RooPlot* xframe = rv_ps[i][eta][pt]->frame(Title(name));
				h_Higgs.plotOn(xframe);
				Gauss.fitTo(h_Higgs, Range(-0.1, 0.1));
// 				if(i > 6)
// 					Gauss.fitTo(h_Higgs, Range(-0.05, 0.05));
				Gauss.plotOn(xframe,RooFit::LineColor(kBlue));
				Gauss.paramOn(xframe, RooFit::Layout(0.13, 0.5, 0.80));

				xframe->Draw();
				
				if(i == 1){
					sigma_FLO->SetBinContent(bin, Sigma.getVal());
					sigma_FLO->SetBinError(bin, Sigma.getError());
				}
				if(i == 2){
					sigma_glb->SetBinContent(bin, Sigma.getVal());
					sigma_glb->SetBinError(bin, Sigma.getError());
				}
				if(i == 3){
					sigma_FLO_vtx->SetBinContent(bin, Sigma.getVal());
					sigma_FLO_vtx->SetBinError(bin, Sigma.getError());
				}
				if(i == 4){
					sigma_glb_vtx->SetBinContent(bin, Sigma.getVal());
					sigma_glb_vtx->SetBinError(bin, Sigma.getError());
				}
				if(i == 5){
					sigma_FLO_single->SetBinContent(bin, Sigma.getVal());
					sigma_FLO_single->SetBinError(bin, Sigma.getError());
				}
				if(i == 6){
					sigma_glb_single->SetBinContent(bin, Sigma.getVal());
					sigma_glb_single->SetBinError(bin, Sigma.getError());
				}
				if(i == 7){
					sigma_stdAlone->SetBinContent(bin, Sigma.getVal());
					sigma_stdAlone->SetBinError(bin, Sigma.getError());
				}
				if(i == 8){
					sigma_stdAlone_vtx->SetBinContent(bin, Sigma.getVal());
					sigma_stdAlone_vtx->SetBinError(bin, Sigma.getError());
				}

				bin++;				
	
				if(pt == 0 && eta == 0)
					c_MC->Print("pt_Fit/" + lepton_type.at(i) + ".pdf[");
				c_MC->Print("pt_Fit/" + lepton_type.at(i) + ".pdf");
				if(pt == pt_bins.size()-2 && eta == eta_bins.size()-2)
					c_MC->Print("pt_Fit/" + lepton_type.at(i) + ".pdf]");
									
			}
		}
	}
	
	TCanvas* b1 = new TCanvas("p_{T} resolution", "p_{T} resolution", 750, 500);
   	b1->cd();
   	TPad* pad11 = new TPad("pad1", "pad1", 0, 0.25, 1, 1);
   	pad11->SetGrid();
   	pad11->SetBottomMargin(0.1);
   	pad11->Draw();
   	pad11->cd();
   	pad11->SetTicks();
   	b1->Update();	
   	sigma_FLO->SetMarkerColor(kRed+1);
	sigma_FLO->SetStats(0);
	sigma_FLO->SetLineColor(kRed+1);
	sigma_glb->SetMarkerColor(kGreen+1);
	sigma_glb->SetLineColor(kGreen+1);
	sigma_FLO->Draw();
	sigma_glb->SetStats(0);
	sigma_glb->Draw("same");
	b1->Update();		
	b1->cd();		
	TPad* pad22 = new TPad("pad2", "pad2", 0, 0.01, 1, 0.285);
	pad22->SetGrid();
	pad22->Draw();
	pad22->cd();
	pad22->SetTicks();
	TH1F* ratio = (TH1F*) sigma_FLO->Clone();
	ratio->Divide(sigma_glb);
	ratio->GetYaxis()->SetRangeUser(0.6, 1.3);
	ratio->SetTitle("");
	ratio->SetStats(0);
	ratio->GetYaxis()->SetTitleSize(0.1);
	ratio->GetYaxis()->SetLabelSize(0.14);
	ratio->GetYaxis()->SetTitle("FLO / glb");
	ratio->GetYaxis()->SetTitleOffset(0.40);
	ratio->GetYaxis()->SetNdivisions(506); 
	ratio->Draw();		
	b1->Update();
	pad22->Update();
	b1->Print("pt_Fit/FLO_glb.pdf");







	b1 = new TCanvas("p_{T} resolution", "p_{T} resolution", 750, 500);
   	b1->cd();
   	pad11 = new TPad("pad1", "pad1", 0, 0.25, 1, 1);
   	pad11->SetGrid();
   	pad11->SetBottomMargin(0.1);
   	pad11->Draw();
   	pad11->cd();
   	pad11->SetTicks();
   	b1->Update();	
   	sigma_FLO->SetMarkerColor(kRed+1);
	sigma_FLO->SetStats(0);
	sigma_FLO->SetLineColor(kRed+1);
	sigma_FLO_vtx->SetMarkerColor(kGreen+1);
	sigma_FLO_vtx->SetLineColor(kGreen+1);
	sigma_FLO_vtx->SetStats(0);
	sigma_FLO_single->SetMarkerColor(kBlue-5);
	sigma_FLO_single->SetLineColor(kBlue-5);
	sigma_FLO_single->SetStats(0);
	sigma_FLO->Draw();
	sigma_FLO_vtx->Draw("same");
	sigma_FLO_single->Draw("same");
	b1->Update();		
	b1->cd();		
	pad22 = new TPad("pad2", "pad2", 0, 0.01, 1, 0.285);
	pad22->SetGrid();
	pad22->Draw();
	pad22->cd();
	pad22->SetTicks();
	TH1F* ratio_1 = (TH1F*) sigma_FLO->Clone();
	TH1F* ratio_2 = (TH1F*) sigma_FLO->Clone();
	ratio_1->Divide(sigma_FLO_vtx);
	ratio_1->GetYaxis()->SetRangeUser(0.6, 1.3);
	ratio_1->SetTitle("");
	ratio_1->SetStats(0);
	ratio_1->GetYaxis()->SetTitleSize(0.1);
	ratio_1->GetYaxis()->SetLabelSize(0.14);
	ratio_1->GetYaxis()->SetTitle("X / Base");
	ratio_1->GetYaxis()->SetTitleOffset(0.40);
	ratio_1->GetYaxis()->SetNdivisions(506); 
	ratio_1->SetLineColor(kGreen+1);
	ratio_1->SetMarkerColor(kGreen+1);
	ratio_1->Draw();		
	ratio_2->Divide(sigma_FLO_single);
	ratio_2->GetYaxis()->SetRangeUser(0.6, 1.3);
	ratio_2->SetTitle("");
	ratio_2->SetStats(0);
	ratio_2->GetYaxis()->SetTitleSize(0.1);
	ratio_2->GetYaxis()->SetLabelSize(0.14);
	ratio_2->GetYaxis()->SetTitle("X / Base");
	ratio_2->GetYaxis()->SetTitleOffset(0.40);
	ratio_2->GetYaxis()->SetNdivisions(506); 
	ratio_2->SetLineColor(kBlue-5);
	ratio_2->SetMarkerColor(kBlue-5);
	ratio_2->Draw();		
	b1->Update();
	pad22->Update();
	b1->Print("pt_Fit/diff_FLO.pdf");

	b1 = new TCanvas("p_{T} resolution", "p_{T} resolution", 750, 500);
   	b1->cd();
   	pad11 = new TPad("pad1", "pad1", 0, 0.25, 1, 1);
   	pad11->SetGrid();
   	pad11->SetBottomMargin(0.1);
   	pad11->Draw();
   	pad11->cd();
   	pad11->SetTicks();
   	b1->Update();	
   	sigma_glb->SetMarkerColor(kRed+1);
	sigma_glb->SetStats(0);
	sigma_glb->SetLineColor(kRed+1);
	sigma_glb_vtx->SetMarkerColor(kGreen+1);
	sigma_glb_vtx->SetLineColor(kGreen+1);
	sigma_glb_vtx->SetStats(0);
	sigma_glb_single->SetMarkerColor(kBlue-5);
	sigma_glb_single->SetLineColor(kBlue-5);
	sigma_glb_single->SetStats(0);
	sigma_glb->Draw();
	sigma_glb_vtx->Draw("same");
	sigma_glb_single->Draw("same");
	b1->Update();		
	b1->cd();		
	pad22 = new TPad("pad2", "pad2", 0, 0.01, 1, 0.285);
	pad22->SetGrid();
	pad22->Draw();
	pad22->cd();
	pad22->SetTicks();
	ratio_1 = (TH1F*) sigma_glb->Clone();
	ratio_2 = (TH1F*) sigma_glb->Clone();
	ratio_1->Divide(sigma_glb_vtx);
	ratio_1->GetYaxis()->SetRangeUser(0.6, 1.3);
	ratio_1->SetTitle("");
	ratio_1->SetStats(0);
	ratio_1->GetYaxis()->SetTitleSize(0.1);
	ratio_1->GetYaxis()->SetLabelSize(0.14);
	ratio_1->GetYaxis()->SetTitle("X / Base");
	ratio_1->GetYaxis()->SetTitleOffset(0.40);
	ratio_1->GetYaxis()->SetNdivisions(506); 
	ratio_1->SetLineColor(kGreen+1);
	ratio_1->SetMarkerColor(kGreen+1);
	ratio_1->Draw();		
	ratio_2->Divide(sigma_glb_single);
	ratio_2->GetYaxis()->SetRangeUser(0.6, 1.3);
	ratio_2->SetTitle("");
	ratio_2->SetStats(0);
	ratio_2->GetYaxis()->SetTitleSize(0.1);
	ratio_2->GetYaxis()->SetLabelSize(0.14);
	ratio_2->GetYaxis()->SetTitle("X / Base");
	ratio_2->GetYaxis()->SetTitleOffset(0.40);
	ratio_2->GetYaxis()->SetNdivisions(506); 
	ratio_2->SetLineColor(kBlue-5);
	ratio_2->SetMarkerColor(kBlue-5);
	ratio_2->Draw();		
	b1->Update();
	pad22->Update();
	b1->Print("pt_Fit/diff_glb.pdf");









	for(int i = 1; i < lepton_type.size()-2; i++){
// 	for(int i = 1; i < 2; i++){
		int bin = 1;
		for(int	vtx = 0; vtx < vtx_bins.size()-1; vtx++){
			TString name = Form("RecoVtx_%s_%.0f_%.0f", lepton_type[i].Data(), vtx_bins.at(vtx), vtx_bins.at(vtx+1));	
			RooDataSet h_Higgs = RooDataSet(Data_nVtx[i][vtx]->GetName(), Data_nVtx[i][vtx]->GetTitle(), Data_nVtx[i][vtx], *Data_nVtx[i][vtx]->get());

			RooRealVar Mean("Mean", "Mean", 0, -0.05, 0.05);
			RooRealVar Sigma("Sigma", "Sigma", 0.05, -0.5, 0.5);

			RooGaussian Gauss("Gauss", "Gauss", *rv_nVtx[i][vtx], Mean, Sigma);

			TCanvas *c_MC = new TCanvas(name, name, 750, 500);
			c_MC->SetFrameFillColor(0);
			c_MC->cd(1)->SetBottomMargin(0.2);
			RooPlot* xframe = rv_nVtx[i][vtx]->frame(Title(name));
			h_Higgs.plotOn(xframe);
// 			Gauss.fitTo(h_Higgs, Range(-0.1, 0.1));
// 			Gauss.plotOn(xframe,RooFit::LineColor(kBlue));
// 			Gauss.paramOn(xframe, RooFit::Layout(0.13, 0.5, 0.80));

			xframe->Draw();
				
			if(i == 1){
				sigma_FLO_nVtx->SetBinContent(bin, Sigma.getVal());
				sigma_FLO_nVtx->SetBinError(bin, Sigma.getError());
			}
			if(i == 2){
				sigma_glb_nVtx->SetBinContent(bin, Sigma.getVal());
				sigma_glb_nVtx->SetBinError(bin, Sigma.getError());
			}
			if(i == 3){
				sigma_FLO_vtx_nVtx->SetBinContent(bin, Sigma.getVal());
				sigma_FLO_vtx_nVtx->SetBinError(bin, Sigma.getError());
			}
			if(i == 4){
				sigma_glb_vtx_nVtx->SetBinContent(bin, Sigma.getVal());
				sigma_glb_vtx_nVtx->SetBinError(bin, Sigma.getError());
			}
			if(i == 5){
				sigma_FLO_single_nVtx->SetBinContent(bin, Sigma.getVal());
				sigma_FLO_single_nVtx->SetBinError(bin, Sigma.getError());
			}
			if(i == 6){
				sigma_glb_single_nVtx->SetBinContent(bin, Sigma.getVal());
				sigma_glb_single_nVtx->SetBinError(bin, Sigma.getError());
			}

			bin++;				
	
			if(vtx == 0)
				c_MC->Print("pt_Fit/" + lepton_type.at(i) + "_nVtx.pdf[");
			c_MC->Print("pt_Fit/" + lepton_type.at(i) + "_nVtx.pdf");
			if(vtx == vtx_bins.size()-2)
				c_MC->Print("pt_Fit/" + lepton_type.at(i) + "_nVtx.pdf]");
									
		}
	}
	
	b1 = new TCanvas("p_{T} resolution", "p_{T} resolution", 750, 500);
   	b1->cd();
   	pad11 = new TPad("pad1", "pad1", 0, 0.25, 1, 1);
   	pad11->SetGrid();
   	pad11->SetBottomMargin(0.1);
   	pad11->Draw();
   	pad11->cd();
   	pad11->SetTicks();
   	b1->Update();	
   	sigma_FLO_nVtx->SetMarkerColor(kRed+1);
	sigma_FLO_nVtx->SetStats(0);
	sigma_FLO_nVtx->SetLineColor(kRed+1);
	sigma_glb->SetMarkerColor(kGreen+1);
	sigma_glb->SetLineColor(kGreen+1);
	sigma_FLO_nVtx->Draw();
	sigma_glb_nVtx->SetStats(0);
	sigma_glb_nVtx->Draw("same");
	b1->Update();		
	b1->cd();		
	pad22 = new TPad("pad2", "pad2", 0, 0.01, 1, 0.285);
	pad22->SetGrid();
	pad22->Draw();
	pad22->cd();
	pad22->SetTicks();
	ratio = (TH1F*) sigma_FLO_nVtx->Clone();
	ratio->Divide(sigma_glb_nVtx);
	ratio->GetYaxis()->SetRangeUser(0.6, 1.3);
	ratio->SetTitle("");
	ratio->SetStats(0);
	ratio->GetYaxis()->SetTitleSize(0.1);
	ratio->GetYaxis()->SetLabelSize(0.14);
	ratio->GetYaxis()->SetTitle("FLO / glb");
	ratio->GetYaxis()->SetTitleOffset(0.40);
	ratio->GetYaxis()->SetNdivisions(506); 
	ratio->Draw();		
	b1->Update();
	pad22->Update();
	b1->Print("pt_Fit/FLO_glb_nVtx.pdf");







	b1 = new TCanvas("p_{T} resolution", "p_{T} resolution", 750, 500);
   	b1->cd();
   	pad11 = new TPad("pad1", "pad1", 0, 0.25, 1, 1);
   	pad11->SetGrid();
   	pad11->SetBottomMargin(0.1);
   	pad11->Draw();
   	pad11->cd();
   	pad11->SetTicks();
   	b1->Update();	
   	sigma_FLO_nVtx->SetMarkerColor(kRed+1);
	sigma_FLO_nVtx->SetStats(0);
	sigma_FLO_nVtx->SetLineColor(kRed+1);
	sigma_FLO_vtx_nVtx->SetMarkerColor(kGreen+1);
	sigma_FLO_vtx_nVtx->SetLineColor(kGreen+1);
	sigma_FLO_vtx_nVtx->SetStats(0);
	sigma_FLO_single_nVtx->SetMarkerColor(kBlue-5);
	sigma_FLO_single_nVtx->SetLineColor(kBlue-5);
	sigma_FLO_single_nVtx->SetStats(0);
	sigma_FLO_nVtx->Draw();
	sigma_FLO_vtx_nVtx->Draw("same");
	sigma_FLO_single_nVtx->Draw("same");
	b1->Update();		
	b1->cd();		
	pad22 = new TPad("pad2", "pad2", 0, 0.01, 1, 0.285);
	pad22->SetGrid();
	pad22->Draw();
	pad22->cd();
	pad22->SetTicks();
	ratio_1 = (TH1F*) sigma_FLO_nVtx->Clone();
	ratio_2 = (TH1F*) sigma_FLO_nVtx->Clone();
	ratio_1->Divide(sigma_FLO_vtx_nVtx);
	ratio_1->GetYaxis()->SetRangeUser(0.6, 1.3);
	ratio_1->SetTitle("");
	ratio_1->SetStats(0);
	ratio_1->GetYaxis()->SetTitleSize(0.1);
	ratio_1->GetYaxis()->SetLabelSize(0.14);
	ratio_1->GetYaxis()->SetTitle("Baseline / X");
	ratio_1->GetYaxis()->SetTitleOffset(0.40);
	ratio_1->GetYaxis()->SetNdivisions(506); 
	ratio_1->SetLineColor(kGreen+1);
	ratio_1->SetMarkerColor(kGreen+1);
	ratio_1->Draw();		
	ratio_2->Divide(sigma_FLO_single_nVtx);
	ratio_2->GetYaxis()->SetRangeUser(0.6, 1.3);
	ratio_2->SetTitle("");
	ratio_2->SetStats(0);
	ratio_2->GetYaxis()->SetTitleSize(0.1);
	ratio_2->GetYaxis()->SetLabelSize(0.14);
	ratio_2->GetYaxis()->SetTitle("Baseline / X");
	ratio_2->GetYaxis()->SetTitleOffset(0.40);
	ratio_2->GetYaxis()->SetNdivisions(506); 
	ratio_2->SetLineColor(kBlue-5);
	ratio_2->SetMarkerColor(kBlue-5);
	ratio_2->Draw();		
	b1->Update();
	pad22->Update();
	b1->Print("pt_Fit/diff_FLO_nVtx.pdf");

	b1 = new TCanvas("p_{T} resolution", "p_{T} resolution", 750, 500);
   	b1->cd();
   	pad11 = new TPad("pad1", "pad1", 0, 0.25, 1, 1);
   	pad11->SetGrid();
   	pad11->SetBottomMargin(0.1);
   	pad11->Draw();
   	pad11->cd();
   	pad11->SetTicks();
   	b1->Update();	
   	sigma_glb_nVtx->SetMarkerColor(kRed+1);
	sigma_glb_nVtx->SetStats(0);
	sigma_glb_nVtx->SetLineColor(kRed+1);
	sigma_glb_vtx_nVtx->SetMarkerColor(kGreen+1);
	sigma_glb_vtx_nVtx->SetLineColor(kGreen+1);
	sigma_glb_vtx_nVtx->SetStats(0);
	sigma_glb_single_nVtx->SetMarkerColor(kBlue-5);
	sigma_glb_single_nVtx->SetLineColor(kBlue-5);
	sigma_glb_single_nVtx->SetStats(0);
	sigma_glb_nVtx->Draw();
	sigma_glb_vtx_nVtx->Draw("same");
	sigma_glb_single_nVtx->Draw("same");
	b1->Update();		
	b1->cd();		
	pad22 = new TPad("pad2", "pad2", 0, 0.01, 1, 0.285);
	pad22->SetGrid();
	pad22->Draw();
	pad22->cd();
	pad22->SetTicks();
	ratio_1 = (TH1F*) sigma_glb_nVtx->Clone();
	ratio_2 = (TH1F*) sigma_glb_nVtx->Clone();
	ratio_1->Divide(sigma_glb_vtx_nVtx);
	ratio_1->GetYaxis()->SetRangeUser(0.6, 1.3);
	ratio_1->SetTitle("");
	ratio_1->SetStats(0);
	ratio_1->GetYaxis()->SetTitleSize(0.1);
	ratio_1->GetYaxis()->SetLabelSize(0.14);
	ratio_1->GetYaxis()->SetTitle("Baseline / X");
	ratio_1->GetYaxis()->SetTitleOffset(0.40);
	ratio_1->GetYaxis()->SetNdivisions(506); 
	ratio_1->SetLineColor(kGreen+1);
	ratio_1->SetMarkerColor(kGreen+1);
	ratio_1->Draw();		
	ratio_2->Divide(sigma_glb_single_nVtx);
	ratio_2->GetYaxis()->SetRangeUser(0.6, 1.3);
	ratio_2->SetTitle("");
	ratio_2->SetStats(0);
	ratio_2->GetYaxis()->SetTitleSize(0.1);
	ratio_2->GetYaxis()->SetLabelSize(0.14);
	ratio_2->GetYaxis()->SetTitle("Baseline / X");
	ratio_2->GetYaxis()->SetTitleOffset(0.40);
	ratio_2->GetYaxis()->SetNdivisions(506); 
	ratio_2->SetLineColor(kBlue-5);
	ratio_2->SetMarkerColor(kBlue-5);
	ratio_2->Draw();		
	b1->Update();
	pad22->Update();
	b1->Print("pt_Fit/diff_glb_nVtx.pdf");	







	b1 = new TCanvas("p_{T} resolution", "p_{T} resolution", 750, 500);
   	b1->cd();
   	pad11 = new TPad("pad1", "pad1", 0, 0.25, 1, 1);
   	pad11->SetGrid();
   	pad11->SetBottomMargin(0.1);
   	pad11->Draw();
   	pad11->cd();
   	pad11->SetTicks();
   	b1->Update();	
   	sigma_stdAlone_vtx->SetMarkerColor(kRed+1);
	sigma_stdAlone_vtx->SetStats(0);
	sigma_stdAlone_vtx->SetLineColor(kRed+1);
	sigma_stdAlone->SetMarkerColor(kGreen+1);
	sigma_stdAlone->SetLineColor(kGreen+1);
	sigma_stdAlone_vtx->Draw();
	sigma_stdAlone->SetStats(0);
	sigma_stdAlone->Draw("same");
	b1->Update();		
	b1->cd();		
	pad22 = new TPad("pad2", "pad2", 0, 0.01, 1, 0.285);
	pad22->SetGrid();
	pad22->Draw();
	pad22->cd();
	pad22->SetTicks();
	ratio = (TH1F*) sigma_stdAlone_vtx->Clone();
	ratio->Divide(sigma_stdAlone);
	ratio->GetYaxis()->SetRangeUser(0.6, 1.3);
	ratio->SetTitle("");
	ratio->SetStats(0);
	ratio->GetYaxis()->SetTitleSize(0.1);
	ratio->GetYaxis()->SetLabelSize(0.14);
	ratio->GetYaxis()->SetTitle("stdAlone_vtx / stdAlone");
	ratio->GetYaxis()->SetTitleOffset(0.40);
	ratio->GetYaxis()->SetNdivisions(506); 
	ratio->Draw();		
	b1->Update();
	pad22->Update();
	b1->Print("pt_Fit/stdAlone_vtx_stdAlone.pdf");




}




std::vector<int> PickMuon(std::vector<float> pt_b, std::vector<float> eta_b, float pt, float eta){

	std::vector<int> k;
	for(int e = 1; e < eta_b.size(); e++){
		if(eta < eta_b.at(e)){
			k.push_back(e);
			break;
		}	
	}
	for(int p = 1; p < pt_b.size(); p++){
		if(pt < pt_b.at(p)){
			k.push_back(p);
			break;
		}	
	}

	return k;

}



float Max(float a, float b){

	if(a > b) return a;
	else return b;

}

void Draw(TH1F *h1, TH1F *h2, TH1F *h3, TString nome_canvas, TString save, TString x_name, TLegend *legend, int fit, bool LogY){

	std::vector<float> param_1;
	std::vector<float> param_2;
	std::vector<float> param_3;
	
	if(fit == 2){
		param_1 = FitMass(h1, "Mass", save);
		param_2 = FitMass(h2, "Mass", save);
		param_3 = FitMass(h3, "Mass", save);
	}
	
	TCanvas *canvas = new TCanvas(nome_canvas, nome_canvas, 750, 750);
   	canvas->cd();
   	TPad* pad11 = new TPad("pad1", "pad1", 0, 0.25, 1, 1);
   	pad11->SetGrid();
   	pad11->SetBottomMargin(0.1);
   	pad11->Draw();
   	pad11->cd();
   	pad11->SetTicks();
   	
	Double_t min, max;

   	if(LogY){
   		pad11->SetLogy();
   		min = 100;
   	}
	else
		min = 0;
	
	float tmp_max;
	tmp_max = Max(h1->GetMaximum(), h2->GetMaximum());
	max = Max(h3->GetMaximum(), tmp_max);
		
	max = 1.1 * max;
	
	if(max == 0) max = 1;	
	
	h1->SetTitle(nome_canvas);
	h2->SetTitle(nome_canvas);
	h3->SetTitle(nome_canvas);
	h1->SetLineColor(kBlack);
	h2->SetLineColor(kRed+1);
	h3->SetLineColor(kGreen+1);

// 	h1->SetStats(0);
// 	h2->SetStats(0);
// 	h3->SetStats(0);

	h1->GetYaxis()->SetRangeUser(min, max);//SetMaximum(max);
	h2->GetYaxis()->SetRangeUser(min, max);//SetMaximum(max);
	h3->GetYaxis()->SetRangeUser(min, max);//SetMaximum(max);
	
	h1->GetXaxis()->SetTitle(x_name);
	h2->GetXaxis()->SetTitle(x_name);
	h3->GetXaxis()->SetTitle(x_name);
	
	canvas->Update();
	
	h1->Draw("E");
	h2->Draw("Same E");
	h3->Draw("Same E");

	TH1F *ratio_2 = (TH1F*) h2->Clone();
	ratio_2->Divide(h1);
	TH1F *ratio_3 = (TH1F*) h3->Clone();
	ratio_3->Divide(h1);

	if(fit == 1){
	
		h1->SetLineColor(kRed+1);
		h2->SetLineColor(kGreen+1);
		h3->SetLineColor(kBlue-5);

		h1->Draw();	
		h1->Fit("gaus", "","", -0.02, 0.02);
		canvas->Update();
		TPaveStats* stats = (TPaveStats*)h1->GetListOfFunctions()->FindObject("stats");
		if(stats){
			stats->SetName("Florida");
			stats->SetY1NDC(.6);
			stats->SetY2NDC(.9);
			stats->SetTextColor(kRed+1);
		}
		else
			std::cout<<"NON 1"<<std::endl;

		h2->Draw();	
		h2->Fit("gaus", "","", -0.02, 0.02);
		canvas->Update();
		TPaveStats* stats_glb = (TPaveStats*)h2->GetListOfFunctions()->FindObject("stats");
		if(stats_glb){
			stats_glb->SetName("glb");
			stats_glb->SetY1NDC(.3);
			stats_glb->SetY2NDC(.6);
			stats_glb->SetTextColor(kGreen+1);
		}
		else
			std::cout<<"NON 2"<<std::endl;

		h3->Draw();	
		h3->Fit("gaus", "","", -0.02, 0.02);
		canvas->Update();
		TPaveStats* stats_single= (TPaveStats*)h3->GetListOfFunctions()->FindObject("stats");
		if(stats_single){
			stats_single->SetName("single");
			stats_single->SetY1NDC(.0);
			stats_single->SetY2NDC(.3);
			stats_single->SetTextColor(kBlue-5);
		}
		else
			std::cout<<"NON 3"<<std::endl;

		h1->Draw();	
		h2->Draw("same");	
		TF1* f = (TF1*) h2->GetListOfFunctions()->FindObject("gaus");
		if(f){
			f->SetLineColor(kGreen+1);
			std::cout<<"OOOK"<<std::endl;
		}
		else
			std::cout<<"\t\t\t\tMMMERDA"<<std::endl;
			
		h3->Draw("same");	
		TF1* f3 = (TF1*) h3->GetListOfFunctions()->FindObject("gaus");
		if(f3){
			f3->SetLineColor(kBlue-5);
			std::cout<<"OOOK"<<std::endl;
		}
		else
			std::cout<<"\t\t\t\tMMMERDA"<<std::endl;

		h1->SetLineColor(kRed+1);
		h2->SetLineColor(kGreen+1);
		h3->SetLineColor(kBlue-5);
		if(stats_glb) stats_glb->Draw();
		if(stats) stats->Draw();	
		if(stats_single) stats_single->Draw();
	}

	if(fit == 2){
	
		TString ciao = Form("Base = %.3f", param_1.at(1));
		TLatex* tex1 = new TLatex(0.65,0.57, ciao);
		tex1->SetNDC();
		tex1->SetTextColor(kRed+1);
		tex1->Draw();

		ciao = Form("VX+BS = %.3f", param_2.at(1));
		tex1 = new TLatex(0.65,0.52, ciao);
		tex1->SetNDC();
		tex1->SetTextColor(kGreen+1);
		tex1->Draw();

		ciao = Form("Single = %.3f", param_3.at(1));
		tex1 = new TLatex(0.65,0.47, ciao);
		tex1->SetNDC();
		tex1->SetTextColor(kBlue-5);
		tex1->Draw();	
	}


/*
// 	legend->Draw();	
	canvas->Update();

	canvas->cd();

	if(fit == 0){
		h1->SetStats(0);
		h2->SetStats(0);	
		h3->SetStats(0);	
		legend->Draw();		

	TPad* pad22 = new TPad("pad2", "pad2", 0, 0.01, 1, 0.25);
	pad22->SetGrid();
	pad22->Draw();
	pad22->cd();
	pad22->SetTicks();
	
	ratio_2->GetYaxis()->SetRangeUser(0.0, 2.0);
	ratio_2->SetTitle("");
	ratio_2->SetStats(0);
	ratio_2->GetYaxis()->SetTitleSize(0.1);
	ratio_2->GetYaxis()->SetLabelSize(0.14);
	ratio_2->GetYaxis()->SetTitle("X / GEN");
	ratio_2->GetYaxis()->SetTitleOffset(0.50);
	ratio_2->GetYaxis()->SetNdivisions(506); 
	ratio_2->SetLineColor(kRed+1);
	
	ratio_3->GetYaxis()->SetRangeUser(0.0, 2.0);
	ratio_3->SetTitle("");
	ratio_3->SetStats(0);
	ratio_3->GetYaxis()->SetTitleSize(0.2);
	ratio_3->GetYaxis()->SetLabelSize(0.14);
	ratio_3->GetYaxis()->SetTitle("X / GEN");
	ratio_3->GetYaxis()->SetTitleOffset(0.50);
	ratio_3->GetYaxis()->SetNdivisions(506); 
	ratio_3->SetLineColor(kGreen+1);

	ratio_2->Draw();
	ratio_3->Draw("same");
	
	canvas->Update();
	pad22->Update();
	TLine *line = new TLine(pad22->GetUxmin(), 1, pad22->GetUxmax(), 1);
	line->SetLineColor(kBlack);
	line->SetLineWidth(1);
	line->Draw();
	TLegend *legend_comparison = new TLegend(0.8,0.6,0.9,0.9);
	legend_comparison->AddEntry(ratio_2, "FLO");
	legend_comparison->AddEntry(ratio_3, "glb");
	legend_comparison->Draw();
	}
*/
	canvas->Print(save);

	std::cout<<"Esco"<<std::endl;
	
}
void Draw(TH1F *h1, TH1F *h2, TString nome_canvas, TString save, TString x_name, TLegend *legend, int fit, bool LogY){

	TCanvas *canvas = new TCanvas(nome_canvas, nome_canvas, 750, 750);
   	canvas->cd();
   	TPad* pad11 = new TPad("pad1", "pad1", 0, 0.25, 1, 1);
   	pad11->SetGrid();
   	pad11->SetBottomMargin(0.1);
   	pad11->Draw();
   	pad11->cd();
   	pad11->SetTicks();
   	
	Double_t min, max;

   	if(LogY){
   		pad11->SetLogy();
   		min = 1;
   	}
	else
		min = 0;
	
	float tmp_max;
	max = Max(h1->GetMaximum(), h2->GetMaximum());
		
	max = 1.1 * max;
	
	if(max == 0) max = 1;	
	
	h1->SetTitle(nome_canvas);
	h2->SetTitle(nome_canvas);
	h1->SetLineColor(kRed+1);
	h2->SetLineColor(kGreen+1);

	h1->GetYaxis()->SetRangeUser(min, max);//SetMaximum(max);
	h2->GetYaxis()->SetRangeUser(min, max);//SetMaximum(max);
	
	h1->GetXaxis()->SetTitle(x_name);
	h2->GetXaxis()->SetTitle(x_name);
	
	canvas->Update();
	
	h1->Draw("E");
	h2->Draw("Same E");

	TH1F *ratio_2 = (TH1F*) h2->Clone();	
	ratio_2->Divide(h1);
	
	if(fit == 1){
		h1->Draw();	
		h1->Fit("gaus", "","", -0.02, 0.02);
		canvas->Update();
		TPaveStats* stats = (TPaveStats*)h1->GetListOfFunctions()->FindObject("stats");
		stats->SetName("Florida");
		stats->SetY1NDC(.6);
		stats->SetY2NDC(.9);
		stats->SetTextColor(kRed+1);

		h2->Draw();	
// 		h2->SetLineColor(kGreen+1);
// 		TF1* gaus = new TF1("gaus","gaus", -0.025, 0.025);
		h2->Fit("gaus", "","", -0.02, 0.02);
		canvas->Update();
// 		h2->GetFunction("gaus")->SetLineColor(kGreen+1);
		TPaveStats* stats_glb = (TPaveStats*)h2->GetListOfFunctions()->FindObject("stats");
		stats_glb->SetName("glb");
		stats_glb->SetY1NDC(.3);
		stats_glb->SetY2NDC(.6);
		stats_glb->SetTextColor(kGreen+1);

		h1->Draw();	
		h2->Draw("same");	
// 		h2->GetFunction("gaus")->SetLineColor(kGreen+1);
		TF1* f = (TF1*) h2->GetListOfFunctions()->FindObject("gaus");
		if(f){
			f->SetLineColor(kGreen+1);
			std::cout<<"OK"<<std::endl;
		}
		else
			std::cout<<"\t\t\t\tMERDA"<<std::endl;
		h1->SetLineColor(kRed+1);
		h2->SetLineColor(kGreen+1);
		stats_glb->Draw();
		stats->Draw();	
	}

	if(fit == 0){
		h1->SetStats(0);
		h2->SetStats(0);	
		legend->Draw();	
	}
			
	canvas->Update();

// 	TH1F *ratio_2 = (TH1F*) h2->Clone();

	canvas->Update();
	canvas->cd();

	TPad* pad22 = new TPad("pad2", "pad2", 0, 0.01, 1, 0.25);
	pad22->SetGrid();
	pad22->Draw();
	pad22->cd();
	pad22->SetTicks();
	
// 	ratio_2->Divide(h1);
	ratio_2->GetYaxis()->SetRangeUser(0.0, 2.0);
	ratio_2->SetTitle("");
	ratio_2->SetStats(0);
	ratio_2->GetYaxis()->SetTitleSize(0.1);
	ratio_2->GetYaxis()->SetLabelSize(0.14);
	ratio_2->GetYaxis()->SetTitle("FLO / glb");
	ratio_2->GetYaxis()->SetTitleOffset(0.50);
	ratio_2->GetYaxis()->SetNdivisions(506); 
	ratio_2->SetLineColor(kRed+1);
	
	ratio_2->Draw();
	
	canvas->Update();
	pad22->Update();
	TLine *line = new TLine(pad22->GetUxmin(), 1, pad22->GetUxmax(), 1);
	line->SetLineColor(kBlack);
	line->SetLineWidth(1);
	line->Draw();
// 	TLegend *legend_comparison = new TLegend(0.8,0.75,0.9,0.9);
// 	legend_comparison->AddEntry(ratio_2, "FLO");
// 	if(fit != 1)
// 		legend_comparison->Draw();

	canvas->Print(save);
	
}
void Draw(TH1F *h1, TH1F *h2, TH1F *h3, RooRealVar* rv_Mass_1, RooRealVar* rv_Mass_2, RooRealVar* rv_Mass_3, RooDataSet* Mass_1, RooDataSet* Mass_2, RooDataSet* Mass_3, TString nome_canvas, TString save, TString x_name, TLegend *legend, bool LogY){
	
	std::vector<float> baseline;
	std::vector<float> VX_BS;
	std::vector<float> Single;

	baseline = FitMass(rv_Mass_1, Mass_1, "Baseline", save);
	VX_BS = FitMass(rv_Mass_2, Mass_2, "VX+BS", save);
	Single = FitMass(rv_Mass_3, Mass_3, "Single", save);

	TCanvas *canvas = new TCanvas(nome_canvas, nome_canvas, 750, 750);
   	canvas->cd();
   	TPad* pad11 = new TPad("pad1", "pad1", 0, 0.25, 1, 1);
   	pad11->SetGrid();
   	pad11->SetBottomMargin(0.1);
   	pad11->Draw();
   	pad11->cd();
   	pad11->SetTicks();
   	
	Double_t min, max;

   	if(LogY){
   		pad11->SetLogy();
   		min = 0.000001;
   	}
	else
		min = 0;
	
// 	h1->Scale(1/h1->Integral());
// 	h2->Scale(1/h2->Integral());
// 	h3->Scale(1/h3->Integral());

	float tmp_max;
	tmp_max = Max(h1->GetMaximum(), h2->GetMaximum());
	max = Max(h3->GetMaximum(), tmp_max);
		
	max = 1.1 * max;
	
	if(max == 0) max = 1;	
	
	h1->SetTitle(nome_canvas);
	h2->SetTitle(nome_canvas);
	h3->SetTitle(nome_canvas);
	h1->SetLineColor(kRed+1);
	h2->SetLineColor(kGreen+1);
	h3->SetLineColor(kBlue-5);

	h1->GetYaxis()->SetRangeUser(min, max);//SetMaximum(max);
	h2->GetYaxis()->SetRangeUser(min, max);//SetMaximum(max);
	h3->GetYaxis()->SetRangeUser(min, max);//SetMaximum(max);
	
	h1->GetXaxis()->SetTitle(x_name);
	h2->GetXaxis()->SetTitle(x_name);
	h3->GetXaxis()->SetTitle(x_name);
	
	canvas->Update();
	
	h1->Draw("E");
	h2->Draw("Same E");
	h3->Draw("Same E");

	legend->Draw();	

	TLatex* tex1 = new TLatex(0.65,0.62, "#sigma:");
	tex1->SetNDC();
	tex1->Draw();
	
	TString ciao = Form("Base = %.3f", baseline.at(1));
	tex1 = new TLatex(0.65,0.57, ciao);
	tex1->SetNDC();
	tex1->SetTextColor(kRed+1);
	tex1->Draw();

	ciao = Form("VX+BS = %.3f", VX_BS.at(1));
	tex1 = new TLatex(0.65,0.52, ciao);
	tex1->SetNDC();
	tex1->SetTextColor(kGreen+1);
	tex1->Draw();

	ciao = Form("Single = %.3f", Single.at(1));
	tex1 = new TLatex(0.65,0.47, ciao);
	tex1->SetNDC();
	tex1->SetTextColor(kBlue-5);
	tex1->Draw();	
		
	canvas->Update();

	TH1F *ratio_2 = (TH1F*) h2->Clone();
	TH1F *ratio_3 = (TH1F*) h3->Clone();

	canvas->Update();
	canvas->cd();

	TPad* pad22 = new TPad("pad2", "pad2", 0, 0.01, 1, 0.25);
	pad22->SetGrid();
	pad22->Draw();
	pad22->cd();
	pad22->SetTicks();
	
	ratio_2->Divide(h1);
	ratio_2->GetYaxis()->SetRangeUser(0.0, 2.0);
	ratio_2->SetTitle("");
	ratio_2->SetStats(0);
	ratio_2->SetLineColor(kGreen+1);
	ratio_2->GetYaxis()->SetTitleSize(0.1);
	ratio_2->GetYaxis()->SetLabelSize(0.14);
	ratio_2->GetYaxis()->SetTitle("X / Baseline");//Reco / Gen");
	ratio_2->GetYaxis()->SetTitleOffset(0.50);
	ratio_2->GetYaxis()->SetNdivisions(506); 
	
	ratio_3->Divide(h1);
	ratio_3->GetYaxis()->SetRangeUser(0.0, 2.0);
	ratio_3->SetTitle("");
	ratio_3->SetStats(0);
	ratio_3->SetLineColor(kBlue-5);
	ratio_3->GetYaxis()->SetTitleSize(0.2);
	ratio_3->GetYaxis()->SetLabelSize(0.14);
	ratio_3->GetYaxis()->SetTitle("X / Baseline");//Reco / Gen");
	ratio_3->GetYaxis()->SetTitleOffset(0.50);
	ratio_3->GetYaxis()->SetNdivisions(506); 

// 	ratio_2->Draw();
// 	ratio_3->Draw("same");
	ratio_3->Draw();
	
	canvas->Update();
	pad22->Update();
	TLine *line = new TLine(pad22->GetUxmin(), 1, pad22->GetUxmax(), 1);
	line->SetLineColor(kBlack);
	line->SetLineWidth(1);
	line->Draw();
	TLegend *legend_comparison = new TLegend(0.75,0.55,0.9,0.9);
	legend_comparison->AddEntry(ratio_2, "VX+BS");
	legend_comparison->AddEntry(ratio_3, "Single");
	legend_comparison->Draw();

	canvas->Print(save);
	
}
std::vector<float> FitMass(TH1F *h1, TString nome_canvas, TString save){

	std::vector<float> param;
	RooRealVar x("x", "Mass (GeV/c^{2})", 60, 120);
	RooDataHist h_Higgs("h_Zboson", "h_Zboson", x, h1);

	//// BW
	RooRealVar BW_mean_DY("BW_mean_DY", "BW_mean_DY", 91.91);
	RooRealVar BW_sigma_DY("BW_sigma_DY", "BW_sigma_DY", 2.44);
	RooBreitWigner BW_DY("BW_DY", "BW_DY", x, BW_mean_DY, BW_sigma_DY);
	//// BW

	//// CB			
	RooRealVar Mean("Mean", "Mean", 0, -2, 2);
	RooRealVar Sigma("Sigma", "Sigma", 1, 0., 10);
	RooRealVar CB_alpha("CB_alpha", "CB_alpha", 5, 0., 30);
	RooRealVar CB_exp("CB_exp", "CB_exp", 5, 0., 30);
	RooCBShape CB("CB", "CB", x, Mean, Sigma, CB_alpha, CB_exp);
	//// CB		

	//// expo	
	RooRealVar tau("tau", "tau", 0, 0, 1);
	RooExponential bkg("bkg","bkg", x, tau);
	//// expo	
		
	RooRealVar fsig("fsig","signal fraction", 0.7, 0.5, 1);

	RooFFTConvPdf BWxCB("BWxCB","BWxCB", x, BW_DY, CB);
	RooAddPdf Final_DY("Final_DY","Final_DY", RooArgList(BWxCB, bkg), fsig);

	TCanvas *c_MC = new TCanvas(nome_canvas, nome_canvas, 900, 700);
	c_MC->SetFrameFillColor(0);
	c_MC->cd(1)->SetBottomMargin(0.2);
	RooPlot* xframe = x.frame(Title(nome_canvas));
	h_Higgs.plotOn(xframe);
	Final_DY.fitTo(h_Higgs, Range(60, 120));
	Final_DY.plotOn(xframe,RooFit::LineColor(kBlue));
// 		RooPlot *framePull_DATA = x.frame("");
// 		framePull_DATA->addObject((TObject*)xframe->pullHist(), "p");
// 		Final_DY.plotOn(xframe, Components("bkg"), LineColor(kBlue), LineStyle(kDashed));
	Final_DY.paramOn(xframe, RooFit::Layout(0.13, 0.5, 0.80));
	xframe->Draw();	
// 		TPad* padPull_DATA =  new  TPad("padPull","padPull",0.,0.,1.,0.2);
// 		padPull_DATA->Draw();
// 		padPull_DATA->cd(0);
// 		framePull_DATA->GetYaxis()->SetLabelSize(0.1);
// 		framePull_DATA->GetXaxis()->SetLabelSize(0.1);
// 		framePull_DATA->SetMinimum(-5.);
// 		framePull_DATA->SetMaximum(5.);
// 		framePull_DATA->Draw();
	c_MC->Print(save);
	
	param.push_back(Mean.getVal());
	param.push_back(Sigma.getVal());
	
	return param;

}
std::vector<float> FitMass(RooRealVar* rv_Mass_1, RooDataSet* Mass_1, TString title, TString save){
			
	std::vector<float> param;
	
	RooDataSet h_Higgs = RooDataSet(Mass_1->GetName(), Mass_1->GetTitle(), Mass_1, *Mass_1->get());//, "1", "weight");

	//// BW
	RooRealVar BW_mean_DY("BW_mean_DY", "BW_mean_DY", 91.91);
	RooRealVar BW_sigma_DY("BW_sigma_DY", "BW_sigma_DY", 2.44);
	RooBreitWigner BW_DY("BW_DY", "BW_DY", *rv_Mass_1, BW_mean_DY, BW_sigma_DY);
	//// BW

	//// CB			
	RooRealVar Mean("Mean", "Mean", 0, -2, 2);
	RooRealVar Sigma("Sigma", "Sigma", 1, 0., 10);
	RooRealVar CB_alpha("CB_alpha", "CB_alpha", 5, 0., 30);
	RooRealVar CB_exp("CB_exp", "CB_exp", 5, 0., 30);
	RooCBShape CB("CB", "CB", *rv_Mass_1, Mean, Sigma, CB_alpha, CB_exp);
	//// CB		

	//// expo	
	RooRealVar tau("tau", "tau", 0, 0, 5);
	RooExponential bkg("bkg","bkg", *rv_Mass_1, tau);
	//// expo	
		
	RooRealVar fsig("fsig","signal fraction", 0.7, 0.5, 1);

	RooFFTConvPdf BWxCB("BWxCB","BWxCB", *rv_Mass_1, BW_DY, CB);
	RooAddPdf Final_DY("Final_DY","Final_DY", RooArgList(BWxCB, bkg), fsig);
		
	TCanvas *c_MC = new TCanvas(title, title, 900, 700);
	c_MC->SetFrameFillColor(0);
	c_MC->cd(1)->SetBottomMargin(0.2);
	RooPlot* xframe = rv_Mass_1->frame(Title(title));
	h_Higgs.plotOn(xframe);
	Final_DY.fitTo(h_Higgs, Range(60, 120));
	Final_DY.plotOn(xframe,RooFit::LineColor(kBlue));
// 	RooPlot *framePull_DATA = rv_Mass_1->frame("");
// 	framePull_DATA->addObject((TObject*)xframe->pullHist(), "p");
	Final_DY.paramOn(xframe, RooFit::Layout(0.13, 0.5, 0.80));

	xframe->Draw();
	
	c_MC->Print(save);
	
	param.push_back(Mean.getVal());
	param.push_back(Sigma.getVal());
	
	return param;

}

