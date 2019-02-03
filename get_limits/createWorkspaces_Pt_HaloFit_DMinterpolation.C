#include <vector>
#include <utility>
#include "TFile.h"
#include "TMath.h"

void addTemplate(string procname, RooArgList& varlist, RooWorkspace& ws, TH1F* hist) {
  RooDataHist rhist(procname.c_str(), "", varlist, hist);
  ws.import(rhist);
}

void makeBinList(string procname, RooRealVar& var, RooWorkspace& ws, TH1F* hist, RooArgList& binlist, bool setConst = false) {
  for (int i = 1; i <= hist->GetNbinsX(); i++) {
    stringstream binss;
    binss << procname << "_bin" << i;
    RooRealVar* binvar;
    if (!setConst) binvar = new RooRealVar(binss.str().c_str(), "", hist->GetBinContent(i), 0.0, 2000.0);
    else           binvar = new RooRealVar(binss.str().c_str(), "", hist->GetBinContent(i));
    binlist.add(*binvar);
  }

  stringstream normss;
  normss << procname << "_norm";

  RooParametricHist phist(procname.c_str(), "", var, binlist, *hist);
  RooAddition norm(normss.str().c_str(), "", binlist);

  // ws.import(phist,RooFit::RecycleConflictNodes());
  ws.import(phist);
  ws.import(norm, RooFit::RecycleConflictNodes());

}


void makeBinList_halo(string procname, RooRealVar& var, RooWorkspace& ws, TH1F* hist, vector<pair<RooRealVar*, TH1*>> syst, RooArgList& binlist) {
  RooRealVar* mu_halo_monophHighPhi_scale = new RooRealVar("mu_halo_monophHighPhi_scale", "mu_halo_monophHighPhi_scale", 1.0, 0.0, 5.0);
  RooRealVar* shape_err = new RooRealVar("Halo_SR_above0p5_MIPTotalEnergy", "Halo_SR_above0p5_MIPTotalEnergy", 0.0, -5.0, 5.0);
  for (int i = 1; i <= hist->GetNbinsX(); i++) {
    stringstream rawmuss;
    rawmuss << procname << "_rawmu_bin" << i;
    RooRealVar* rawmu_bin = new RooRealVar(rawmuss.str().c_str(), "", hist->GetBinContent(i)); // Constant
    RooArgList fobinlist;
    fobinlist.add(*mu_halo_monophHighPhi_scale);
    fobinlist.add(*rawmu_bin);
    fobinlist.add(*shape_err);
    
    stringstream binss;
    binss << procname << "_bin" << i;
    stringstream formss;
    formss << "@0*@1";
    formss << "*(TMath::Exp(" << hist->GetBinError(i)/hist->GetBinContent(i) << "*@2))";
    for (int j = 0; j < syst.size(); j++) {
      stringstream systbinss;
      if (syst[j].first == NULL) {
        systbinss << procname << "_" << syst[j].second->GetName() << "_bin" << i;
        RooRealVar* systbinvar = new RooRealVar(systbinss.str().c_str(), "", 0., -5., 5.);
        fobinlist.add(*systbinvar);
      }
      else {
        fobinlist.add(*syst[j].first);
      }
      formss << "*(TMath::Exp(" << syst[j].second->GetBinContent(i) << "*@" << j+3 << "))";
    }
    RooFormulaVar* binvar = new RooFormulaVar(binss.str().c_str(), "", formss.str().c_str(), RooArgList(fobinlist));
    binlist.add(*binvar);
  }

  stringstream normss;
  normss << procname << "_norm";

  RooParametricHist phist(procname.c_str(), "", var, binlist, *hist);
  RooAddition norm(normss.str().c_str(), "", binlist);

  // ws.import(phist,RooFit::RecycleConflictNodes());
  ws.import(phist);
  ws.import(norm, RooFit::RecycleConflictNodes());

}
// ///////////////////////
// // RooProduct version
// ///////////////////////
// void makeBinList_halo(string procname, RooRealVar& var, RooWorkspace& ws, TH1F* hist, RooArgList& binlist) {
//   RooRealVar* mu_halo_monophHighPhi_scale = new RooRealVar("mu_halo_monophHighPhi_scale", "mu_halo_monophHighPhi_scale", 1.0, 0.0, 5.0);
//   for (int i = 1; i <= hist->GetNbinsX(); i++) {
//     stringstream rawmuss;
//     rawmuss << procname << "_rawmu_bin" << i;
//     stringstream prodss;
//     prodss << procname << "_prod_bin" << i;
//     RooRealVar* rawmu_bin = new RooRealVar(rawmuss.str().c_str(), "", hist->GetBinContent(i)); // Constant
//     RooProduct* binvar = new RooProduct(prodss.str().c_str(), "", RooArgSet(*mu_halo_monophHighPhi_scale,*rawmu_bin));
//     binlist.add(*binvar);
//   }

//   stringstream normss;
//   normss << procname << "_norm";

//   RooParametricHist phist(procname.c_str(), "", var, binlist, *hist);
//   RooAddition norm(normss.str().c_str(), "", binlist);

//   // ws.import(phist,RooFit::RecycleConflictNodes());
//   ws.import(phist);
//   ws.import(norm, RooFit::RecycleConflictNodes());

// }

void makeConnectedBinList(string procname, RooRealVar& var, RooWorkspace& ws, TH1F* rhist, vector<pair<RooRealVar*, TH1*>> syst, const RooArgList& srbinlist, RooArgList* crbinlist=NULL) {
  if (crbinlist == NULL) crbinlist = new RooArgList();

  for (int i = 1; i <= rhist->GetNbinsX(); i++) {
    stringstream rbinss;
    rbinss << "r_" << procname << "_bin" << i;
    RooRealVar* rbinvar = new RooRealVar(rbinss.str().c_str(), "", rhist->GetBinContent(i));

    stringstream rerrbinss;
    rerrbinss << procname << "_bin" << i << "_Runc";
    RooRealVar* rerrbinvar = new RooRealVar(rerrbinss.str().c_str(), "", 0., -5., 5.);

    stringstream binss;
    binss << procname << "_bin" << i;

    RooArgList fobinlist;
    fobinlist.add(srbinlist[i-1]);
    fobinlist.add(*rbinvar);
    fobinlist.add(*rerrbinvar);

    stringstream formss;
    formss << "@0/";
    formss << "(";
    formss << "@1";
    formss << "*(TMath::Exp(" << rhist->GetBinError(i)/rhist->GetBinContent(i) << "*@2))";
    for (int j = 0; j < syst.size(); j++) {
      stringstream systbinss;
      if (syst[j].first == NULL) {
  systbinss << procname << "_" << syst[j].second->GetName() << "_bin" << i;
  RooRealVar* systbinvar = new RooRealVar(systbinss.str().c_str(), "", 0., -5., 5.);
  fobinlist.add(*systbinvar);
      }
      else {
  fobinlist.add(*syst[j].first);
      }
      formss << "*(TMath::Exp(" << syst[j].second->GetBinContent(i) << "*@" << j+3 << "))";
    }
    formss << ")";

    RooFormulaVar* binvar = new RooFormulaVar(binss.str().c_str(), "", formss.str().c_str(), RooArgList(fobinlist));
    crbinlist->add(*binvar);
  }

  stringstream normss;
  normss << procname << "_norm";

  RooParametricHist phist(procname.c_str(), "", var, *crbinlist, *rhist);
  RooAddition norm(normss.str().c_str(),"", *crbinlist);

  ws.import(phist, RooFit::RecycleConflictNodes());
  ws.import(norm, RooFit::RecycleConflictNodes());
}

void makeConnectedBinList_noStatErr(string procname, RooRealVar& var, RooWorkspace& ws, TH1F* rhist, vector<pair<RooRealVar*, TH1*>> syst, const RooArgList& srbinlist, RooArgList* crbinlist=NULL) {
  if (crbinlist == NULL) crbinlist = new RooArgList();

  for (int i = 1; i <= rhist->GetNbinsX(); i++) {
    stringstream rbinss;
    rbinss << "r_" << procname << "_bin" << i;
    RooRealVar* rbinvar = new RooRealVar(rbinss.str().c_str(), "", rhist->GetBinContent(i));

    stringstream binss;
    binss << procname << "_bin" << i;

    RooArgList fobinlist;
    fobinlist.add(srbinlist[i-1]);
    fobinlist.add(*rbinvar);

    stringstream formss;
    formss << "@0/";
    formss << "(";
    formss << "@1";
    for (int j = 0; j < syst.size(); j++) {
      stringstream systbinss;
      if (syst[j].first == NULL) {
  systbinss << procname << "_" << syst[j].second->GetName() << "_bin" << i;
  RooRealVar* systbinvar = new RooRealVar(systbinss.str().c_str(), "", 0., -5., 5.);
  fobinlist.add(*systbinvar);
      }
      else {
  fobinlist.add(*syst[j].first);
      }
      formss << "*(TMath::Exp(" << syst[j].second->GetBinContent(i) << "*@" << j+3 << "))";
    }
    formss << ")";

    RooFormulaVar* binvar = new RooFormulaVar(binss.str().c_str(), "", formss.str().c_str(), RooArgList(fobinlist));
    crbinlist->add(*binvar);
  }

  stringstream normss;
  normss << procname << "_norm";

  RooParametricHist phist(procname.c_str(), "", var, *crbinlist, *rhist);
  RooAddition norm(normss.str().c_str(),"", *crbinlist);

  ws.import(phist, RooFit::RecycleConflictNodes());
  ws.import(norm, RooFit::RecycleConflictNodes());
}

void do_createWorkspace(string signal_samplename, bool connectWZ, std::vector<string> &signal_samplenames, std::vector<float> &signal_multipliers, string signal_histos_filename_stub){
  gSystem->Load("libHiggsAnalysisCombinedLimit.so");
  
  string outfile_stub = "/nfs_scratch/jjbuchanan/CMSSW_8_1_0/src/HiggsAnalysis/CombinedLimit/SimultaneousFit/workspaces_DMinterpolation/";
  outfile_stub += "workspace_Pt_";
  if(connectWZ) outfile_stub += "yesZW_";
  else outfile_stub += "noZW_";
  
  TFile *outfile = new TFile(TString(outfile_stub+signal_samplename+".root"),"RECREATE");
  RooWorkspace wspace("w","w");

  // RooRealVar mt("mt","M_{T}",0,1200);
  // RooRealVar pfmet("pfmet","pfMET",170,1000);
  RooRealVar phopt("phopt","E_{T}",175,1000);
  // RooArgList vars(mt);
  // RooArgList vars(pfmet);
  RooArgList vars(phopt);

  // Templates
  TFile* transfer_factors_file = new TFile("transfer_factors_Pt.root");
    
  TFile* signalabove0p5file = new TFile(TString(signal_histos_filename_stub+"_histos_above0p5_Pt.root"));
  TFile* signalbelow0p5file = new TFile(TString(signal_histos_filename_stub+"_histos_below0p5_Pt.root"));
  
  TFile* ZnnGabove0p5file = new TFile("ZnnG_histos_above0p5_Pt.root");
  TFile* ZnnGbelow0p5file = new TFile("ZnnG_histos_below0p5_Pt.root");

  TFile* WenGfile = new TFile("WenG_histos_Pt.root");
  
  TFile* WmnGfile = new TFile("WmnG_histos_Pt.root");
  
  TFile* ZeeGfile = new TFile("ZeeG_histos_Pt.root");
  
  TFile* ZmmGfile = new TFile("ZmmG_histos_Pt.root");
  
  // ---------------------------- SIGNAL REGION (above0p5) -------------------------------------------------------------------//
  string procname_ZG_SA = "ZnunuG_SR_above0p5";
  string procname_WG_SA = "WG_SR_above0p5";
  string procname_halo_SA = "Halo_SR_above0p5";
  
  TH1F* data_obs_SR_above0p5 = (TH1F*)ZnnGabove0p5file->Get("data_obs");
  const int nBins = data_obs_SR_above0p5->GetNbinsX();
  
  // Data
  addTemplate("data_obs_SR_above0p5", vars, wspace, data_obs_SR_above0p5);

  // Signal shape
  TH1F* signal_SR_above0p5_hist = (TH1F*)signalabove0p5file->Get(TString("histo_"+signal_samplename));
  // For printout at the end
  signal_samplenames.push_back(signal_samplename);
  float signal_sampleyield_above0p5 = signal_SR_above0p5_hist->Integral();
  float signal_multiplier = 55.0 / signal_sampleyield_above0p5; // Normalize so that combine limits are close to 1
  signal_multipliers.push_back(signal_multiplier);
    
  TH1F* signal_SR_above0p5_DMgenPhoETinterpStatUp_hist = (TH1F*)signalabove0p5file->Get(TString("histo_"+signal_samplename+"_DMgenPhoETinterpStatUp"));
  TH1F* signal_SR_above0p5_DMgenPhoETinterpStatDown_hist = (TH1F*)signalabove0p5file->Get(TString("histo_"+signal_samplename+"_DMgenPhoETinterpStatDown"));
  TH1F* signal_SR_above0p5_DMNLOoverLOStatUp_hist = (TH1F*)signalabove0p5file->Get(TString("histo_"+signal_samplename+"_DMNLOoverLOStatUp"));
  TH1F* signal_SR_above0p5_DMNLOoverLOStatDown_hist = (TH1F*)signalabove0p5file->Get(TString("histo_"+signal_samplename+"_DMNLOoverLOStatDown"));
  TH1F* signal_SR_above0p5_JESUp_hist = (TH1F*)signalabove0p5file->Get(TString("histo_"+signal_samplename+"_JESUp"));
  TH1F* signal_SR_above0p5_JESDown_hist = (TH1F*)signalabove0p5file->Get(TString("histo_"+signal_samplename+"_JESDown"));
  TH1F* signal_SR_above0p5_PESUp_hist = (TH1F*)signalabove0p5file->Get(TString("histo_"+signal_samplename+"_PESUp"));
  TH1F* signal_SR_above0p5_PESDown_hist = (TH1F*)signalabove0p5file->Get(TString("histo_"+signal_samplename+"_PESDown"));
  signal_SR_above0p5_hist->Scale(signal_multiplier);
  signal_SR_above0p5_DMgenPhoETinterpStatUp_hist->Scale(signal_multiplier);
  signal_SR_above0p5_DMgenPhoETinterpStatDown_hist->Scale(signal_multiplier);
  signal_SR_above0p5_DMNLOoverLOStatUp_hist->Scale(signal_multiplier);
  signal_SR_above0p5_DMNLOoverLOStatDown_hist->Scale(signal_multiplier);
  signal_SR_above0p5_JESUp_hist->Scale(signal_multiplier);
  signal_SR_above0p5_JESDown_hist->Scale(signal_multiplier);
  signal_SR_above0p5_PESUp_hist->Scale(signal_multiplier);
  signal_SR_above0p5_PESDown_hist->Scale(signal_multiplier);
  addTemplate("Signal_SR_above0p5", vars, wspace, signal_SR_above0p5_hist);
  addTemplate("Signal_SR_above0p5_DMgenPhoETinterpStatUp", vars, wspace, signal_SR_above0p5_DMgenPhoETinterpStatUp_hist);
  addTemplate("Signal_SR_above0p5_DMgenPhoETinterpStatDown", vars, wspace, signal_SR_above0p5_DMgenPhoETinterpStatDown_hist);
  addTemplate("Signal_SR_above0p5_DMNLOoverLOStatUp", vars, wspace, signal_SR_above0p5_DMNLOoverLOStatUp_hist);
  addTemplate("Signal_SR_above0p5_DMNLOoverLOStatDown", vars, wspace, signal_SR_above0p5_DMNLOoverLOStatDown_hist);
  addTemplate("Signal_SR_above0p5_JESUp", vars, wspace, signal_SR_above0p5_JESUp_hist);
  addTemplate("Signal_SR_above0p5_JESDown", vars, wspace, signal_SR_above0p5_JESDown_hist);
  addTemplate("Signal_SR_above0p5_PESUp", vars, wspace, signal_SR_above0p5_PESUp_hist);
  addTemplate("Signal_SR_above0p5_PESDown", vars, wspace, signal_SR_above0p5_PESDown_hist);
  
  TH1F* znng_SR_above0p5_hist = (TH1F*)ZnnGabove0p5file->Get("histo_ZNuNuG");
  RooArgList znng_SR_above0p5_bins;
  makeBinList(procname_ZG_SA, phopt, wspace, znng_SR_above0p5_hist, znng_SR_above0p5_bins);
  
  TH1F* halo_SR_above0p5_hist = (TH1F*)ZnnGabove0p5file->Get("histo_bhalo");
  float halo_above0p5_int = halo_SR_above0p5_hist->Integral();
  TH1F* halo_SR_below0p5_hist = (TH1F*)((TH1F*)ZnnGbelow0p5file->Get("histo_bhalo"))->Clone("histo_bhalo_below");
  halo_SR_above0p5_hist->Add(halo_SR_below0p5_hist);
  halo_SR_above0p5_hist->Scale(halo_above0p5_int/(halo_SR_above0p5_hist->Integral()));
  RooArgList halo_SR_above0p5_bins;
  TH1F* halo_SR_above0p5_MIPTotEnergy_shiftUp = (TH1F*)ZnnGabove0p5file->Get("histo_bhalo_MIPTotEnergyUp");
  float halo_above0p5_MIPTotEnergy_shiftUp_int = halo_SR_above0p5_MIPTotEnergy_shiftUp->Integral();
  TH1F* halo_SR_below0p5_MIPTotEnergy_shiftUp = (TH1F*)((TH1F*)ZnnGbelow0p5file->Get("histo_bhalo_MIPTotEnergyUp"))->Clone("histo_bhalo_MIPTotEnergyUp_below");
  halo_SR_above0p5_MIPTotEnergy_shiftUp->Add(halo_SR_below0p5_MIPTotEnergy_shiftUp);
  halo_SR_above0p5_MIPTotEnergy_shiftUp->Scale(halo_above0p5_MIPTotEnergy_shiftUp_int/(halo_SR_above0p5_MIPTotEnergy_shiftUp->Integral()));
  TH1F* halo_SR_above0p5_MIPTotEnergy_shiftDown = (TH1F*)ZnnGabove0p5file->Get("histo_bhalo_MIPTotEnergyUp");
  float halo_above0p5_MIPTotEnergy_shiftDown_int = halo_SR_above0p5_MIPTotEnergy_shiftDown->Integral();
  TH1F* halo_SR_below0p5_MIPTotEnergy_shiftDown = (TH1F*)((TH1F*)ZnnGbelow0p5file->Get("histo_bhalo_MIPTotEnergyUp"))->Clone("histo_bhalo_MIPTotEnergyUp_below");
  halo_SR_above0p5_MIPTotEnergy_shiftDown->Add(halo_SR_below0p5_MIPTotEnergy_shiftDown);
  halo_SR_above0p5_MIPTotEnergy_shiftDown->Scale(halo_above0p5_MIPTotEnergy_shiftDown_int/(halo_SR_above0p5_MIPTotEnergy_shiftDown->Integral()));
  TH1F* halo_SR_above0p5_MIPTotEnergy_fractionalShifts = (TH1F*)halo_SR_above0p5_MIPTotEnergy_shiftUp->Clone("halo_SR_above0p5_MIPTotEnergy_fractionalShifts");
  for (int i = 1; i <= halo_SR_above0p5_MIPTotEnergy_fractionalShifts->GetNbinsX(); i++) {
    Float_t upshift = halo_SR_above0p5_MIPTotEnergy_shiftUp->GetBinContent(i)/halo_SR_above0p5_hist->GetBinContent(i) - 1.0;
    Float_t downshift = halo_SR_above0p5_MIPTotEnergy_shiftDown->GetBinContent(i)/halo_SR_above0p5_hist->GetBinContent(i) - 1.0;
    Float_t shiftEnvelope = TMath::Max(fabs(upshift), fabs(downshift));
    halo_SR_above0p5_MIPTotEnergy_fractionalShifts->SetBinContent(i, shiftEnvelope);
  }
  vector<pair<RooRealVar*, TH1*>> halo_SR_above0p5_syst;
  RooRealVar* halo_SR_above0p5_MIPTotEnergy = new RooRealVar("Halo_SR_above0p5_MIPTotEnergy", "", 0., -5., 5.);
  halo_SR_above0p5_syst.push_back(pair<RooRealVar*, TH1*>(halo_SR_above0p5_MIPTotEnergy, halo_SR_above0p5_MIPTotEnergy_fractionalShifts));
  makeBinList_halo(procname_halo_SA, phopt, wspace, halo_SR_above0p5_hist, halo_SR_above0p5_syst, halo_SR_above0p5_bins);

  // Without Z/W link
  TH1F* wlng_SR_above0p5_hist = (TH1F*)ZnnGabove0p5file->Get("histo_WG");
  RooArgList wlng_SR_above0p5_bins;
  // With Z/W link
  TH1F* zoverw_SRr_above0p5_hist = (TH1F*)transfer_factors_file->Get("transfer_factor_WGinZnnG_to_ZnnGinZnnG_above0p5_noShift");
  
  TH1F* zoverw_SRr_above0p5_ZNuNuGJets_qcdscale_shiftUp = (TH1F*)transfer_factors_file->Get("transfer_factor_WGinZnnG_to_ZnnGinZnnG_above0p5_ZNuNuGJets_qcdscale_shiftUp");
  TH1F* zoverw_SRr_above0p5_ZNuNuGJets_qcdscale_shiftDown = (TH1F*)transfer_factors_file->Get("transfer_factor_WGinZnnG_to_ZnnGinZnnG_above0p5_ZNuNuGJets_qcdscale_shiftDown");
  TH1F* zoverw_SRr_above0p5_ZNuNuGJets_qcdscale_fractionalShifts = (TH1F*)zoverw_SRr_above0p5_ZNuNuGJets_qcdscale_shiftUp->Clone(TString(procname_WG_SA+"_ZNuNuGJets_qcdscale"));
  for (int i = 1; i <= zoverw_SRr_above0p5_ZNuNuGJets_qcdscale_fractionalShifts->GetNbinsX(); i++) {
    Float_t upshift = zoverw_SRr_above0p5_ZNuNuGJets_qcdscale_shiftUp->GetBinContent(i)/zoverw_SRr_above0p5_hist->GetBinContent(i) - 1.0;
    Float_t downshift = zoverw_SRr_above0p5_ZNuNuGJets_qcdscale_shiftDown->GetBinContent(i)/zoverw_SRr_above0p5_hist->GetBinContent(i) - 1.0;
    Float_t shiftEnvelope = TMath::Max(fabs(upshift), fabs(downshift));
    zoverw_SRr_above0p5_ZNuNuGJets_qcdscale_fractionalShifts->SetBinContent(i, shiftEnvelope);
  }
  TH1F* zoverw_SRr_above0p5_ZNuNuGJets_ewkscale_shiftUp = (TH1F*)transfer_factors_file->Get("transfer_factor_WGinZnnG_to_ZnnGinZnnG_above0p5_ZNuNuGJets_ewkscale_shiftUp");
  TH1F* zoverw_SRr_above0p5_ZNuNuGJets_ewkscale_shiftDown = (TH1F*)transfer_factors_file->Get("transfer_factor_WGinZnnG_to_ZnnGinZnnG_above0p5_ZNuNuGJets_ewkscale_shiftDown");
  TH1F* zoverw_SRr_above0p5_ZNuNuGJets_ewkscale_fractionalShifts = (TH1F*)zoverw_SRr_above0p5_ZNuNuGJets_ewkscale_shiftUp->Clone(TString(procname_WG_SA+"_ZNuNuGJets_ewkscale"));
  for (int i = 1; i <= zoverw_SRr_above0p5_ZNuNuGJets_ewkscale_fractionalShifts->GetNbinsX(); i++) {
    Float_t upshift = zoverw_SRr_above0p5_ZNuNuGJets_ewkscale_shiftUp->GetBinContent(i)/zoverw_SRr_above0p5_hist->GetBinContent(i) - 1.0;
    Float_t downshift = zoverw_SRr_above0p5_ZNuNuGJets_ewkscale_shiftDown->GetBinContent(i)/zoverw_SRr_above0p5_hist->GetBinContent(i) - 1.0;
    Float_t shiftEnvelope = TMath::Max(fabs(upshift), fabs(downshift));
    zoverw_SRr_above0p5_ZNuNuGJets_ewkscale_fractionalShifts->SetBinContent(i, shiftEnvelope);
  }
  TH1F* zoverw_SRr_above0p5_ZNuNuGJets_ewkshape_shiftUp = (TH1F*)transfer_factors_file->Get("transfer_factor_WGinZnnG_to_ZnnGinZnnG_above0p5_ZNuNuGJets_ewkshape_shiftUp");
  TH1F* zoverw_SRr_above0p5_ZNuNuGJets_ewkshape_shiftDown = (TH1F*)transfer_factors_file->Get("transfer_factor_WGinZnnG_to_ZnnGinZnnG_above0p5_ZNuNuGJets_ewkshape_shiftDown");
  TH1F* zoverw_SRr_above0p5_ZNuNuGJets_ewkshape_fractionalShifts = (TH1F*)zoverw_SRr_above0p5_ZNuNuGJets_ewkshape_shiftUp->Clone(TString(procname_WG_SA+"_ZNuNuGJets_ewkshape"));
  for (int i = 1; i <= zoverw_SRr_above0p5_ZNuNuGJets_ewkshape_fractionalShifts->GetNbinsX(); i++) {
    Float_t upshift = zoverw_SRr_above0p5_ZNuNuGJets_ewkshape_shiftUp->GetBinContent(i)/zoverw_SRr_above0p5_hist->GetBinContent(i) - 1.0;
    Float_t downshift = zoverw_SRr_above0p5_ZNuNuGJets_ewkshape_shiftDown->GetBinContent(i)/zoverw_SRr_above0p5_hist->GetBinContent(i) - 1.0;
    Float_t shiftEnvelope = TMath::Max(fabs(upshift), fabs(downshift));
    zoverw_SRr_above0p5_ZNuNuGJets_ewkshape_fractionalShifts->SetBinContent(i, shiftEnvelope);
  }
  TH1F* zoverw_SRr_above0p5_ZNuNuGJets_gamma_shiftUp = (TH1F*)transfer_factors_file->Get("transfer_factor_WGinZnnG_to_ZnnGinZnnG_above0p5_ZNuNuGJets_gamma_shiftUp");
  TH1F* zoverw_SRr_above0p5_ZNuNuGJets_gamma_shiftDown = (TH1F*)transfer_factors_file->Get("transfer_factor_WGinZnnG_to_ZnnGinZnnG_above0p5_ZNuNuGJets_gamma_shiftDown");
  TH1F* zoverw_SRr_above0p5_ZNuNuGJets_gamma_fractionalShifts = (TH1F*)zoverw_SRr_above0p5_ZNuNuGJets_gamma_shiftUp->Clone(TString(procname_WG_SA+"_ZNuNuGJets_gamma"));
  for (int i = 1; i <= zoverw_SRr_above0p5_ZNuNuGJets_gamma_fractionalShifts->GetNbinsX(); i++) {
    Float_t upshift = zoverw_SRr_above0p5_ZNuNuGJets_gamma_shiftUp->GetBinContent(i)/zoverw_SRr_above0p5_hist->GetBinContent(i) - 1.0;
    Float_t downshift = zoverw_SRr_above0p5_ZNuNuGJets_gamma_shiftDown->GetBinContent(i)/zoverw_SRr_above0p5_hist->GetBinContent(i) - 1.0;
    Float_t shiftEnvelope = TMath::Max(fabs(upshift), fabs(downshift));
    zoverw_SRr_above0p5_ZNuNuGJets_gamma_fractionalShifts->SetBinContent(i, shiftEnvelope);
  }
  TH1F* zoverw_SRr_above0p5_WGJets_qcdscale_shiftUp = (TH1F*)transfer_factors_file->Get("transfer_factor_WGinZnnG_to_ZnnGinZnnG_above0p5_WGJets_qcdscale_shiftUp");
  TH1F* zoverw_SRr_above0p5_WGJets_qcdscale_shiftDown = (TH1F*)transfer_factors_file->Get("transfer_factor_WGinZnnG_to_ZnnGinZnnG_above0p5_WGJets_qcdscale_shiftDown");
  TH1F* zoverw_SRr_above0p5_WGJets_qcdscale_fractionalShifts = (TH1F*)zoverw_SRr_above0p5_WGJets_qcdscale_shiftUp->Clone(TString(procname_WG_SA+"_WGJets_qcdscale"));
  for (int i = 1; i <= zoverw_SRr_above0p5_WGJets_qcdscale_fractionalShifts->GetNbinsX(); i++) {
    Float_t upshift = zoverw_SRr_above0p5_WGJets_qcdscale_shiftUp->GetBinContent(i)/zoverw_SRr_above0p5_hist->GetBinContent(i) - 1.0;
    Float_t downshift = zoverw_SRr_above0p5_WGJets_qcdscale_shiftDown->GetBinContent(i)/zoverw_SRr_above0p5_hist->GetBinContent(i) - 1.0;
    Float_t shiftEnvelope = TMath::Max(fabs(upshift), fabs(downshift));
    zoverw_SRr_above0p5_WGJets_qcdscale_fractionalShifts->SetBinContent(i, shiftEnvelope);
  }
  TH1F* zoverw_SRr_above0p5_WGJets_ewkscale_shiftUp = (TH1F*)transfer_factors_file->Get("transfer_factor_WGinZnnG_to_ZnnGinZnnG_above0p5_WGJets_ewkscale_shiftUp");
  TH1F* zoverw_SRr_above0p5_WGJets_ewkscale_shiftDown = (TH1F*)transfer_factors_file->Get("transfer_factor_WGinZnnG_to_ZnnGinZnnG_above0p5_WGJets_ewkscale_shiftDown");
  TH1F* zoverw_SRr_above0p5_WGJets_ewkscale_fractionalShifts = (TH1F*)zoverw_SRr_above0p5_WGJets_ewkscale_shiftUp->Clone(TString(procname_WG_SA+"_WGJets_ewkscale"));
  for (int i = 1; i <= zoverw_SRr_above0p5_WGJets_ewkscale_fractionalShifts->GetNbinsX(); i++) {
    Float_t upshift = zoverw_SRr_above0p5_WGJets_ewkscale_shiftUp->GetBinContent(i)/zoverw_SRr_above0p5_hist->GetBinContent(i) - 1.0;
    Float_t downshift = zoverw_SRr_above0p5_WGJets_ewkscale_shiftDown->GetBinContent(i)/zoverw_SRr_above0p5_hist->GetBinContent(i) - 1.0;
    Float_t shiftEnvelope = TMath::Max(fabs(upshift), fabs(downshift));
    zoverw_SRr_above0p5_WGJets_ewkscale_fractionalShifts->SetBinContent(i, shiftEnvelope);
  }
  TH1F* zoverw_SRr_above0p5_WGJets_ewkshape_shiftUp = (TH1F*)transfer_factors_file->Get("transfer_factor_WGinZnnG_to_ZnnGinZnnG_above0p5_WGJets_ewkshape_shiftUp");
  TH1F* zoverw_SRr_above0p5_WGJets_ewkshape_shiftDown = (TH1F*)transfer_factors_file->Get("transfer_factor_WGinZnnG_to_ZnnGinZnnG_above0p5_WGJets_ewkshape_shiftDown");
  TH1F* zoverw_SRr_above0p5_WGJets_ewkshape_fractionalShifts = (TH1F*)zoverw_SRr_above0p5_WGJets_ewkshape_shiftUp->Clone(TString(procname_WG_SA+"_WGJets_ewkshape"));
  for (int i = 1; i <= zoverw_SRr_above0p5_WGJets_ewkshape_fractionalShifts->GetNbinsX(); i++) {
    Float_t upshift = zoverw_SRr_above0p5_WGJets_ewkshape_shiftUp->GetBinContent(i)/zoverw_SRr_above0p5_hist->GetBinContent(i) - 1.0;
    Float_t downshift = zoverw_SRr_above0p5_WGJets_ewkshape_shiftDown->GetBinContent(i)/zoverw_SRr_above0p5_hist->GetBinContent(i) - 1.0;
    Float_t shiftEnvelope = TMath::Max(fabs(upshift), fabs(downshift));
    zoverw_SRr_above0p5_WGJets_ewkshape_fractionalShifts->SetBinContent(i, shiftEnvelope);
  }
  TH1F* zoverw_SRr_above0p5_WGJets_gamma_shiftUp = (TH1F*)transfer_factors_file->Get("transfer_factor_WGinZnnG_to_ZnnGinZnnG_above0p5_WGJets_gamma_shiftUp");
  TH1F* zoverw_SRr_above0p5_WGJets_gamma_shiftDown = (TH1F*)transfer_factors_file->Get("transfer_factor_WGinZnnG_to_ZnnGinZnnG_above0p5_WGJets_gamma_shiftDown");
  TH1F* zoverw_SRr_above0p5_WGJets_gamma_fractionalShifts = (TH1F*)zoverw_SRr_above0p5_WGJets_gamma_shiftUp->Clone(TString(procname_WG_SA+"_WGJets_gamma"));
  for (int i = 1; i <= zoverw_SRr_above0p5_WGJets_gamma_fractionalShifts->GetNbinsX(); i++) {
    Float_t upshift = zoverw_SRr_above0p5_WGJets_gamma_shiftUp->GetBinContent(i)/zoverw_SRr_above0p5_hist->GetBinContent(i) - 1.0;
    Float_t downshift = zoverw_SRr_above0p5_WGJets_gamma_shiftDown->GetBinContent(i)/zoverw_SRr_above0p5_hist->GetBinContent(i) - 1.0;
    Float_t shiftEnvelope = TMath::Max(fabs(upshift), fabs(downshift));
    zoverw_SRr_above0p5_WGJets_gamma_fractionalShifts->SetBinContent(i, shiftEnvelope);
  }
  // TH1F* zoverw_SRr_above0p5_anticorrelated_ewkscale_shiftUp = (TH1F*)transfer_factors_file->Get("transfer_factor_WGinZnnG_to_ZnnGinZnnG_above0p5_anticorrelated_ewkscale_shiftUp");
  // TH1F* zoverw_SRr_above0p5_anticorrelated_ewkscale_shiftDown = (TH1F*)transfer_factors_file->Get("transfer_factor_WGinZnnG_to_ZnnGinZnnG_above0p5_anticorrelated_ewkscale_shiftDown");
  // TH1F* zoverw_SRr_above0p5_anticorrelated_ewkscale_fractionalShifts = (TH1F*)zoverw_SRr_above0p5_anticorrelated_ewkscale_shiftUp->Clone(TString(procname_WG_SA+"_anticorrelated_ewkscale"));
  // for (int i = 1; i <= zoverw_SRr_above0p5_anticorrelated_ewkscale_fractionalShifts->GetNbinsX(); i++) {
  //   Float_t upshift = zoverw_SRr_above0p5_anticorrelated_ewkscale_shiftUp->GetBinContent(i)/zoverw_SRr_above0p5_hist->GetBinContent(i) - 1.0;
  //   Float_t downshift = zoverw_SRr_above0p5_anticorrelated_ewkscale_shiftDown->GetBinContent(i)/zoverw_SRr_above0p5_hist->GetBinContent(i) - 1.0;
  //   Float_t shiftEnvelope = TMath::Max(fabs(upshift), fabs(downshift));
  //   zoverw_SRr_above0p5_anticorrelated_ewkscale_fractionalShifts->SetBinContent(i, shiftEnvelope);
  // }
  // TH1F* zoverw_SRr_above0p5_anticorrelated_ewkshape_shiftUp = (TH1F*)transfer_factors_file->Get("transfer_factor_WGinZnnG_to_ZnnGinZnnG_above0p5_anticorrelated_ewkshape_shiftUp");
  // TH1F* zoverw_SRr_above0p5_anticorrelated_ewkshape_shiftDown = (TH1F*)transfer_factors_file->Get("transfer_factor_WGinZnnG_to_ZnnGinZnnG_above0p5_anticorrelated_ewkshape_shiftDown");
  // TH1F* zoverw_SRr_above0p5_anticorrelated_ewkshape_fractionalShifts = (TH1F*)zoverw_SRr_above0p5_anticorrelated_ewkshape_shiftUp->Clone(TString(procname_WG_SA+"_anticorrelated_ewkshape"));
  // for (int i = 1; i <= zoverw_SRr_above0p5_anticorrelated_ewkshape_fractionalShifts->GetNbinsX(); i++) {
  //   Float_t upshift = zoverw_SRr_above0p5_anticorrelated_ewkshape_shiftUp->GetBinContent(i)/zoverw_SRr_above0p5_hist->GetBinContent(i) - 1.0;
  //   Float_t downshift = zoverw_SRr_above0p5_anticorrelated_ewkshape_shiftDown->GetBinContent(i)/zoverw_SRr_above0p5_hist->GetBinContent(i) - 1.0;
  //   Float_t shiftEnvelope = TMath::Max(fabs(upshift), fabs(downshift));
  //   zoverw_SRr_above0p5_anticorrelated_ewkshape_fractionalShifts->SetBinContent(i, shiftEnvelope);
  // }
  // TH1F* zoverw_SRr_above0p5_anticorrelated_gamma_shiftUp = (TH1F*)transfer_factors_file->Get("transfer_factor_WGinZnnG_to_ZnnGinZnnG_above0p5_anticorrelated_gamma_shiftUp");
  // TH1F* zoverw_SRr_above0p5_anticorrelated_gamma_shiftDown = (TH1F*)transfer_factors_file->Get("transfer_factor_WGinZnnG_to_ZnnGinZnnG_above0p5_anticorrelated_gamma_shiftDown");
  // TH1F* zoverw_SRr_above0p5_anticorrelated_gamma_fractionalShifts = (TH1F*)zoverw_SRr_above0p5_anticorrelated_gamma_shiftUp->Clone(TString(procname_WG_SA+"_anticorrelated_gamma"));
  // for (int i = 1; i <= zoverw_SRr_above0p5_anticorrelated_gamma_fractionalShifts->GetNbinsX(); i++) {
  //   Float_t upshift = zoverw_SRr_above0p5_anticorrelated_gamma_shiftUp->GetBinContent(i)/zoverw_SRr_above0p5_hist->GetBinContent(i) - 1.0;
  //   Float_t downshift = zoverw_SRr_above0p5_anticorrelated_gamma_shiftDown->GetBinContent(i)/zoverw_SRr_above0p5_hist->GetBinContent(i) - 1.0;
  //   Float_t shiftEnvelope = TMath::Max(fabs(upshift), fabs(downshift));
  //   zoverw_SRr_above0p5_anticorrelated_gamma_fractionalShifts->SetBinContent(i, shiftEnvelope);
  // }
  TH1F* zoverw_SRr_above0p5_pdf_shiftUp = (TH1F*)transfer_factors_file->Get("transfer_factor_WGinZnnG_to_ZnnGinZnnG_above0p5_pdf_shiftUp");
  TH1F* zoverw_SRr_above0p5_pdf_shiftDown = (TH1F*)transfer_factors_file->Get("transfer_factor_WGinZnnG_to_ZnnGinZnnG_above0p5_pdf_shiftDown");
  TH1F* zoverw_SRr_above0p5_pdf_fractionalShifts = (TH1F*)zoverw_SRr_above0p5_pdf_shiftUp->Clone(TString(procname_WG_SA+"_pdf"));
  for (int i = 1; i <= zoverw_SRr_above0p5_pdf_fractionalShifts->GetNbinsX(); i++) {
    Float_t upshift = zoverw_SRr_above0p5_pdf_shiftUp->GetBinContent(i)/zoverw_SRr_above0p5_hist->GetBinContent(i) - 1.0;
    Float_t downshift = zoverw_SRr_above0p5_pdf_shiftDown->GetBinContent(i)/zoverw_SRr_above0p5_hist->GetBinContent(i) - 1.0;
    Float_t shiftEnvelope = TMath::Max(fabs(upshift), fabs(downshift));
    zoverw_SRr_above0p5_pdf_fractionalShifts->SetBinContent(i, shiftEnvelope);
  }
  vector<pair<RooRealVar*, TH1*>> zoverw_SR_above0p5_syst;
  RooRealVar* zoverw_SR_above0p5_pdf = new RooRealVar("ZNuNuGoverWG_SR_above0p5_pdf", "", 0., -5., 5.);
  zoverw_SR_above0p5_syst.push_back(pair<RooRealVar*, TH1*>(zoverw_SR_above0p5_pdf, zoverw_SRr_above0p5_pdf_fractionalShifts));
  RooRealVar* zoverw_SR_above0p5_ZNuNuGJets_qcdscale = new RooRealVar("ZNuNuGoverWG_SR_above0p5_ZNuNuGJets_qcdscale", "", 0., -5., 5.);
  zoverw_SR_above0p5_syst.push_back(pair<RooRealVar*, TH1*>(zoverw_SR_above0p5_ZNuNuGJets_qcdscale, zoverw_SRr_above0p5_ZNuNuGJets_qcdscale_fractionalShifts));
  RooRealVar* zoverw_SR_above0p5_ZNuNuGJets_ewkscale = new RooRealVar("ZNuNuGoverWG_SR_above0p5_ZNuNuGJets_ewkscale", "", 0., -5., 5.);
  zoverw_SR_above0p5_syst.push_back(pair<RooRealVar*, TH1*>(zoverw_SR_above0p5_ZNuNuGJets_ewkscale, zoverw_SRr_above0p5_ZNuNuGJets_ewkscale_fractionalShifts));
  RooRealVar* zoverw_SR_above0p5_ZNuNuGJets_ewkshape = new RooRealVar("ZNuNuGoverWG_SR_above0p5_ZNuNuGJets_ewkshape", "", 0., -5., 5.);
  zoverw_SR_above0p5_syst.push_back(pair<RooRealVar*, TH1*>(zoverw_SR_above0p5_ZNuNuGJets_ewkshape, zoverw_SRr_above0p5_ZNuNuGJets_ewkshape_fractionalShifts));
  RooRealVar* zoverw_SR_above0p5_ZNuNuGJets_gamma = new RooRealVar("ZNuNuGoverWG_SR_above0p5_ZNuNuGJets_gamma", "", 0., -5., 5.);
  zoverw_SR_above0p5_syst.push_back(pair<RooRealVar*, TH1*>(zoverw_SR_above0p5_ZNuNuGJets_gamma, zoverw_SRr_above0p5_ZNuNuGJets_gamma_fractionalShifts));
  RooRealVar* zoverw_SR_above0p5_WGJets_qcdscale = new RooRealVar("ZNuNuGoverWG_SR_above0p5_WGJets_qcdscale", "", 0., -5., 5.);
  zoverw_SR_above0p5_syst.push_back(pair<RooRealVar*, TH1*>(zoverw_SR_above0p5_WGJets_qcdscale, zoverw_SRr_above0p5_WGJets_qcdscale_fractionalShifts));
  RooRealVar* zoverw_SR_above0p5_WGJets_ewkscale = new RooRealVar("ZNuNuGoverWG_SR_above0p5_WGJets_ewkscale", "", 0., -5., 5.);
  zoverw_SR_above0p5_syst.push_back(pair<RooRealVar*, TH1*>(zoverw_SR_above0p5_WGJets_ewkscale, zoverw_SRr_above0p5_WGJets_ewkscale_fractionalShifts));
  RooRealVar* zoverw_SR_above0p5_WGJets_ewkshape = new RooRealVar("ZNuNuGoverWG_SR_above0p5_WGJets_ewkshape", "", 0., -5., 5.);
  zoverw_SR_above0p5_syst.push_back(pair<RooRealVar*, TH1*>(zoverw_SR_above0p5_WGJets_ewkshape, zoverw_SRr_above0p5_WGJets_ewkshape_fractionalShifts));
  RooRealVar* zoverw_SR_above0p5_WGJets_gamma = new RooRealVar("ZNuNuGoverWG_SR_above0p5_WGJets_gamma", "", 0., -5., 5.);
  zoverw_SR_above0p5_syst.push_back(pair<RooRealVar*, TH1*>(zoverw_SR_above0p5_WGJets_gamma, zoverw_SRr_above0p5_WGJets_gamma_fractionalShifts));
  // RooRealVar* zoverw_SR_above0p5_anticorrelated_qcdscale = new RooRealVar("ZNuNuGoverWG_SR_above0p5_anticorrelated_qcdscale", "", 0., -5., 5.);
  // zoverw_SR_above0p5_syst.push_back(pair<RooRealVar*, TH1*>(zoverw_SR_above0p5_anticorrelated_qcdscale, zoverw_SRr_above0p5_anticorrelated_qcdscale_fractionalShifts));
  // RooRealVar* zoverw_SR_above0p5_anticorrelated_ewkscale = new RooRealVar("ZNuNuGoverWG_SR_above0p5_anticorrelated_ewkscale", "", 0., -5., 5.);
  // zoverw_SR_above0p5_syst.push_back(pair<RooRealVar*, TH1*>(zoverw_SR_above0p5_anticorrelated_ewkscale, zoverw_SRr_above0p5_anticorrelated_ewkscale_fractionalShifts));
  // RooRealVar* zoverw_SR_above0p5_anticorrelated_ewkshape = new RooRealVar("ZNuNuGoverWG_SR_above0p5_anticorrelated_ewkshape", "", 0., -5., 5.);
  // zoverw_SR_above0p5_syst.push_back(pair<RooRealVar*, TH1*>(zoverw_SR_above0p5_anticorrelated_ewkshape, zoverw_SRr_above0p5_anticorrelated_ewkshape_fractionalShifts));
  // RooRealVar* zoverw_SR_above0p5_anticorrelated_gamma = new RooRealVar("ZNuNuGoverWG_SR_above0p5_anticorrelated_gamma", "", 0., -5., 5.);
  // zoverw_SR_above0p5_syst.push_back(pair<RooRealVar*, TH1*>(zoverw_SR_above0p5_anticorrelated_gamma, zoverw_SRr_above0p5_anticorrelated_gamma_fractionalShifts));

  if (!connectWZ) makeBinList(procname_WG_SA, phopt, wspace, wlng_SR_above0p5_hist, wlng_SR_above0p5_bins);
  else   makeConnectedBinList(procname_WG_SA, phopt, wspace, zoverw_SRr_above0p5_hist, zoverw_SR_above0p5_syst, znng_SR_above0p5_bins, &wlng_SR_above0p5_bins);
  
  // TH1F* hist_ZNuNuG_SR_above0p5_corrected = (TH1F*)ZnnGabove0p5file->Get("histo_ZNuNuG"); // Only if not fitting
  // TH1F* hist_ZNuNuG_SR_above0p5_uncorrected = (TH1F*)ZnnGabove0p5file->Get("histo_ZNuNuG_uncorrected"); // Only if not fitting
  // TH1F* hist_ZNuNuG_SR_above0p5_qcdscale = (TH1F*)ZnnGabove0p5file->Get("histo_ZNuNuG_qcdscale"); // Only if not fitting
  // TH1F* hist_ZNuNuG_SR_above0p5_EWKUp = (TH1F*)hist_ZNuNuG_SR_above0p5_corrected->Clone("hist_ZNuNuG_SR_above0p5_EWKUp"); // Only if not fitting
  // TH1F* hist_ZNuNuG_SR_above0p5_EWKDown = (TH1F*)hist_ZNuNuG_SR_above0p5_corrected->Clone("hist_ZNuNuG_SR_above0p5_EWKDown"); // Only if not fitting
  // TH1F* hist_ZNuNuG_SR_above0p5_qcdscaleUp = (TH1F*)hist_ZNuNuG_SR_above0p5_corrected->Clone("hist_ZNuNuG_SR_above0p5_qcdscaleUp"); // Only if not fitting
  // TH1F* hist_ZNuNuG_SR_above0p5_qcdscaleDown = (TH1F*)hist_ZNuNuG_SR_above0p5_corrected->Clone("hist_ZNuNuG_SR_above0p5_qcdscaleDown"); // Only if not fitting
  // TH1F* hist_WG_SR_above0p5_corrected = (TH1F*)ZnnGabove0p5file->Get("histo_WG"); // Only if not fitting
  // TH1F* hist_WG_SR_above0p5_uncorrected = (TH1F*)ZnnGabove0p5file->Get("histo_WG_uncorrected"); // Only if not fitting
  // TH1F* hist_WG_SR_above0p5_qcdscale = (TH1F*)ZnnGabove0p5file->Get("histo_WG_qcdscale"); // Only if not fitting
  // TH1F* hist_WG_SR_above0p5_EWKUp = (TH1F*)hist_WG_SR_above0p5_corrected->Clone("hist_WG_SR_above0p5_EWKUp"); // Only if not fitting
  // TH1F* hist_WG_SR_above0p5_EWKDown = (TH1F*)hist_WG_SR_above0p5_corrected->Clone("hist_WG_SR_above0p5_EWKDown"); // Only if not fitting
  // TH1F* hist_WG_SR_above0p5_qcdscaleUp = (TH1F*)hist_WG_SR_above0p5_corrected->Clone("hist_WG_SR_above0p5_qcdscaleUp"); // Only if not fitting
  // TH1F* hist_WG_SR_above0p5_qcdscaleDown = (TH1F*)hist_WG_SR_above0p5_corrected->Clone("hist_WG_SR_above0p5_qcdscaleDown"); // Only if not fitting
  // for(int i = 1; i <= nBins; i++){ // Only if not fitting
  //   Float_t diff = fabs(hist_ZNuNuG_SR_above0p5_corrected->GetBinContent(i) - hist_ZNuNuG_SR_above0p5_uncorrected->GetBinContent(i)); // Only if not fitting
  //   hist_ZNuNuG_SR_above0p5_EWKUp->SetBinContent(i, hist_ZNuNuG_SR_above0p5_corrected->GetBinContent(i) + diff); // Only if not fitting
  //   hist_ZNuNuG_SR_above0p5_EWKDown->SetBinContent(i, hist_ZNuNuG_SR_above0p5_corrected->GetBinContent(i) - diff); // Only if not fitting
  //   diff = fabs(hist_ZNuNuG_SR_above0p5_corrected->GetBinContent(i) - hist_ZNuNuG_SR_above0p5_qcdscale->GetBinContent(i)); // Only if not fitting
  //   hist_ZNuNuG_SR_above0p5_qcdscaleUp->SetBinContent(i, hist_ZNuNuG_SR_above0p5_corrected->GetBinContent(i) + diff); // Only if not fitting
  //   hist_ZNuNuG_SR_above0p5_qcdscaleDown->SetBinContent(i, hist_ZNuNuG_SR_above0p5_corrected->GetBinContent(i) - diff); // Only if not fitting
  //   diff = fabs(hist_WG_SR_above0p5_corrected->GetBinContent(i) - hist_WG_SR_above0p5_uncorrected->GetBinContent(i)); // Only if not fitting
  //   hist_WG_SR_above0p5_EWKUp->SetBinContent(i, hist_WG_SR_above0p5_corrected->GetBinContent(i) + diff); // Only if not fitting
  //   hist_WG_SR_above0p5_EWKDown->SetBinContent(i, hist_WG_SR_above0p5_corrected->GetBinContent(i) - diff); // Only if not fitting
  //   diff = fabs(hist_WG_SR_above0p5_corrected->GetBinContent(i) - hist_WG_SR_above0p5_qcdscale->GetBinContent(i)); // Only if not fitting
  //   hist_WG_SR_above0p5_qcdscaleUp->SetBinContent(i, hist_WG_SR_above0p5_corrected->GetBinContent(i) + diff); // Only if not fitting
  //   hist_WG_SR_above0p5_qcdscaleDown->SetBinContent(i, hist_WG_SR_above0p5_corrected->GetBinContent(i) - diff); // Only if not fitting
  // } // Only if not fitting

  
  // Data driven backgrounds
  addTemplate("QCD_SR_above0p5"                , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_jetfake"));
  addTemplate("QCD_SR_above0p5_QCDrUp"         , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_jetfake_errUp"));
  addTemplate("QCD_SR_above0p5_QCDrDown"       , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_jetfake_errDown"));
  addTemplate("Elefake_SR_above0p5"            , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_elefake"));
  // addTemplate("BHalo_SR_above0p5"              , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_bhalo"));
  addTemplate("Spikes_SR_above0p5"             , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_spikes"));
  // MC backgrounds
  // addTemplate("ZNuNuG_SR_above0p5"              , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_ZNuNuG")); // Only if not fitting
  // addTemplate("ZNuNuG_SR_above0p5_JESUp"        , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_ZNuNuG_JESUp")); // Only if not fitting
  // addTemplate("ZNuNuG_SR_above0p5_JESDown"      , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_ZNuNuG_JESDown")); // Only if not fitting
  // addTemplate("ZNuNuG_SR_above0p5_PESUp"        , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_ZNuNuG_PESUp")); // Only if not fitting
  // addTemplate("ZNuNuG_SR_above0p5_PESDown"      , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_ZNuNuG_PESDown")); // Only if not fitting
  // addTemplate("ZNuNuG_SR_above0p5_phoSFUp"        , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_ZNuNuG_phoSFUp")); // Only if not fitting
  // addTemplate("ZNuNuG_SR_above0p5_phoSFDown"      , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_ZNuNuG_phoSFDown")); // Only if not fitting
  // addTemplate("ZNuNuG_SR_above0p5_EWKUp"      , vars, wspace, hist_ZNuNuG_SR_above0p5_EWKUp); // Only if not fitting
  // addTemplate("ZNuNuG_SR_above0p5_EWKDown"      , vars, wspace, hist_ZNuNuG_SR_above0p5_EWKDown); // Only if not fitting
  // addTemplate("ZNuNuG_SR_above0p5_qcdscaleUp"      , vars, wspace, hist_ZNuNuG_SR_above0p5_qcdscaleUp); // Only if not fitting
  // addTemplate("ZNuNuG_SR_above0p5_qcdscaleDown"      , vars, wspace, hist_ZNuNuG_SR_above0p5_qcdscaleDown); // Only if not fitting
  // addTemplate("WG_SR_above0p5"              , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_WG")); // Only if not fitting
  // addTemplate("WG_SR_above0p5_JESUp"        , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_WG_JESUp")); // Only if not fitting
  // addTemplate("WG_SR_above0p5_JESDown"      , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_WG_JESDown")); // Only if not fitting
  // addTemplate("WG_SR_above0p5_PESUp"        , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_WG_PESUp")); // Only if not fitting
  // addTemplate("WG_SR_above0p5_PESDown"      , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_WG_PESDown")); // Only if not fitting
  // addTemplate("WG_SR_above0p5_phoSFUp"        , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_WG_phoSFUp")); // Only if not fitting
  // addTemplate("WG_SR_above0p5_phoSFDown"      , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_WG_phoSFDown")); // Only if not fitting
  // addTemplate("WG_SR_above0p5_EWKUp"      , vars, wspace, hist_WG_SR_above0p5_EWKUp); // Only if not fitting
  // addTemplate("WG_SR_above0p5_EWKDown"      , vars, wspace, hist_WG_SR_above0p5_EWKDown); // Only if not fitting
  // addTemplate("WG_SR_above0p5_qcdscaleUp"      , vars, wspace, hist_WG_SR_above0p5_qcdscaleUp); // Only if not fitting
  // addTemplate("WG_SR_above0p5_qcdscaleDown"      , vars, wspace, hist_WG_SR_above0p5_qcdscaleDown); // Only if not fitting
  addTemplate("GJets_SR_above0p5"              , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_GJets"));
  addTemplate("GJets_SR_above0p5_JESUp"        , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_GJets_JESUp"));
  addTemplate("GJets_SR_above0p5_JESDown"      , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_GJets_JESDown"));
  addTemplate("GJets_SR_above0p5_PESUp"        , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_GJets_PESUp"));
  addTemplate("GJets_SR_above0p5_PESDown"      , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_GJets_PESDown"));
  addTemplate("ZllG_SR_above0p5"               , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_ZllG_combined"));
  addTemplate("ZllG_SR_above0p5_JESUp"         , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_ZllG_JESUp_combined"));
  addTemplate("ZllG_SR_above0p5_JESDown"       , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_ZllG_JESDown_combined"));
  addTemplate("ZllG_SR_above0p5_PESUp"         , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_ZllG_PESUp_combined"));
  addTemplate("ZllG_SR_above0p5_PESDown"       , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_ZllG_PESDown_combined"));
  addTemplate("TTG_SR_above0p5"                , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_TTG"));
  addTemplate("TTG_SR_above0p5_JESUp"          , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_TTG_JESUp"));
  addTemplate("TTG_SR_above0p5_JESDown"        , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_TTG_JESDown"));
  addTemplate("TTG_SR_above0p5_PESUp"          , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_TTG_PESUp"));
  addTemplate("TTG_SR_above0p5_PESDown"        , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_TTG_PESDown"));
  addTemplate("TG_SR_above0p5"                 , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_TG"));
  addTemplate("TG_SR_above0p5_JESUp"           , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_TG_JESUp"));
  addTemplate("TG_SR_above0p5_JESDown"         , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_TG_JESDown"));
  addTemplate("TG_SR_above0p5_PESUp"           , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_TG_PESUp"));
  addTemplate("TG_SR_above0p5_PESDown"         , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_TG_PESDown"));
  addTemplate("Diphoton_SR_above0p5"           , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_diphoton"));
  addTemplate("Diphoton_SR_above0p5_JESUp"     , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_diphoton_JESUp"));
  addTemplate("Diphoton_SR_above0p5_JESDown"   , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_diphoton_JESDown"));
  addTemplate("Diphoton_SR_above0p5_PESUp"     , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_diphoton_PESUp"));
  addTemplate("Diphoton_SR_above0p5_PESDown"   , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_diphoton_PESDown"));
  addTemplate("WZ_SR_above0p5"                 , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_WZ"));
  addTemplate("WZ_SR_above0p5_JESUp"           , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_WZ_JESUp"));
  addTemplate("WZ_SR_above0p5_JESDown"         , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_WZ_JESDown"));
  addTemplate("WZ_SR_above0p5_PESUp"           , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_WZ_PESUp"));
  addTemplate("WZ_SR_above0p5_PESDown"         , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_WZ_PESDown"));
  addTemplate("ZZ_SR_above0p5"                 , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_ZZ"));
  addTemplate("ZZ_SR_above0p5_JESUp"           , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_ZZ_JESUp"));
  addTemplate("ZZ_SR_above0p5_JESDown"         , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_ZZ_JESDown"));
  addTemplate("ZZ_SR_above0p5_PESUp"           , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_ZZ_PESUp"));
  addTemplate("ZZ_SR_above0p5_PESDown"         , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_ZZ_PESDown"));
  addTemplate("WMuNu_SR_above0p5"              , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_WMuNu"));
  addTemplate("WMuNu_SR_above0p5_JESUp"        , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_WMuNu_JESUp"));
  addTemplate("WMuNu_SR_above0p5_JESDown"      , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_WMuNu_JESDown"));
  addTemplate("WMuNu_SR_above0p5_PESUp"        , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_WMuNu_PESUp"));
  addTemplate("WMuNu_SR_above0p5_PESDown"      , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_WMuNu_PESDown"));
  addTemplate("WTauNu_SR_above0p5"             , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_WTauNu"));
  addTemplate("WTauNu_SR_above0p5_JESUp"       , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_WTauNu_JESUp"));
  addTemplate("WTauNu_SR_above0p5_JESDown"     , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_WTauNu_JESDown"));
  addTemplate("WTauNu_SR_above0p5_PESUp"       , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_WTauNu_PESUp"));
  addTemplate("WTauNu_SR_above0p5_PESDown"     , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_WTauNu_PESDown"));
  addTemplate("WW_SR_above0p5"                 , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_WW"));
  addTemplate("WW_SR_above0p5_JESUp"           , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_WW_JESUp"));
  addTemplate("WW_SR_above0p5_JESDown"         , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_WW_JESDown"));
  addTemplate("WW_SR_above0p5_PESUp"           , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_WW_PESUp"));
  addTemplate("WW_SR_above0p5_PESDown"         , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_WW_PESDown"));
  
  // ---------------------------- SIGNAL REGION (below0p5) -------------------------------------------------------------------//
  string procname_ZG_SB = "ZnunuG_SR_below0p5";
  string procname_WG_SB = "WG_SR_below0p5";
  string procname_halo_SB = "Halo_SR_below0p5";
  
  // Data
  addTemplate("data_obs_SR_below0p5", vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("data_obs"));
  
  // Signal shape
  TH1F* signal_SR_below0p5_hist = (TH1F*)signalbelow0p5file->Get(TString("histo_"+signal_samplename));
  TH1F* signal_SR_below0p5_DMgenPhoETinterpStatUp_hist = (TH1F*)signalbelow0p5file->Get(TString("histo_"+signal_samplename+"_DMgenPhoETinterpStatUp"));
  TH1F* signal_SR_below0p5_DMgenPhoETinterpStatDown_hist = (TH1F*)signalbelow0p5file->Get(TString("histo_"+signal_samplename+"_DMgenPhoETinterpStatDown"));
  TH1F* signal_SR_below0p5_DMNLOoverLOStatUp_hist = (TH1F*)signalbelow0p5file->Get(TString("histo_"+signal_samplename+"_DMNLOoverLOStatUp"));
  TH1F* signal_SR_below0p5_DMNLOoverLOStatDown_hist = (TH1F*)signalbelow0p5file->Get(TString("histo_"+signal_samplename+"_DMNLOoverLOStatDown"));
  TH1F* signal_SR_below0p5_JESUp_hist = (TH1F*)signalbelow0p5file->Get(TString("histo_"+signal_samplename+"_JESUp"));
  TH1F* signal_SR_below0p5_JESDown_hist = (TH1F*)signalbelow0p5file->Get(TString("histo_"+signal_samplename+"_JESDown"));
  TH1F* signal_SR_below0p5_PESUp_hist = (TH1F*)signalbelow0p5file->Get(TString("histo_"+signal_samplename+"_PESUp"));
  TH1F* signal_SR_below0p5_PESDown_hist = (TH1F*)signalbelow0p5file->Get(TString("histo_"+signal_samplename+"_PESDown"));
  // For printout at the end
  signal_SR_below0p5_hist->Scale(signal_multiplier);
  signal_SR_below0p5_DMgenPhoETinterpStatUp_hist->Scale(signal_multiplier);
  signal_SR_below0p5_DMgenPhoETinterpStatDown_hist->Scale(signal_multiplier);
  signal_SR_below0p5_DMNLOoverLOStatUp_hist->Scale(signal_multiplier);
  signal_SR_below0p5_DMNLOoverLOStatDown_hist->Scale(signal_multiplier);
  signal_SR_below0p5_JESUp_hist->Scale(signal_multiplier);
  signal_SR_below0p5_JESDown_hist->Scale(signal_multiplier);
  signal_SR_below0p5_PESUp_hist->Scale(signal_multiplier);
  signal_SR_below0p5_PESDown_hist->Scale(signal_multiplier);
  addTemplate("Signal_SR_below0p5", vars, wspace, signal_SR_below0p5_hist);
  addTemplate("Signal_SR_below0p5_DMgenPhoETinterpStatUp", vars, wspace, signal_SR_below0p5_DMgenPhoETinterpStatUp_hist);
  addTemplate("Signal_SR_below0p5_DMgenPhoETinterpStatDown", vars, wspace, signal_SR_below0p5_DMgenPhoETinterpStatDown_hist);
  addTemplate("Signal_SR_below0p5_DMNLOoverLOStatUp", vars, wspace, signal_SR_below0p5_DMNLOoverLOStatUp_hist);
  addTemplate("Signal_SR_below0p5_DMNLOoverLOStatDown", vars, wspace, signal_SR_below0p5_DMNLOoverLOStatDown_hist);
  addTemplate("Signal_SR_below0p5_JESUp", vars, wspace, signal_SR_below0p5_JESUp_hist);
  addTemplate("Signal_SR_below0p5_JESDown", vars, wspace, signal_SR_below0p5_JESDown_hist);
  addTemplate("Signal_SR_below0p5_PESUp", vars, wspace, signal_SR_below0p5_PESUp_hist);
  addTemplate("Signal_SR_below0p5_PESDown", vars, wspace, signal_SR_below0p5_PESDown_hist);
    
  TH1F* znng_SRr_aboveOverBelow_hist = (TH1F*)transfer_factors_file->Get("transfer_factor_SBtoSA");
  vector<pair<RooRealVar*, TH1*>>   znng_SR_below0p5_syst;
  makeConnectedBinList_noStatErr(procname_ZG_SB, phopt, wspace, znng_SRr_aboveOverBelow_hist, znng_SR_below0p5_syst, znng_SR_above0p5_bins);
  
  TH1F* wlng_SRr_aboveOverBelow_hist = (TH1F*)znng_SRr_aboveOverBelow_hist->Clone("wlng_SRr_aboveOverBelow_hist");
  vector<pair<RooRealVar*, TH1*>>   wlng_SR_below0p5_syst;
  makeConnectedBinList_noStatErr(procname_WG_SB, phopt, wspace, wlng_SRr_aboveOverBelow_hist, wlng_SR_below0p5_syst, wlng_SR_above0p5_bins);
  
  TH1F* halo_SRr_aboveOverBelow_hist = (TH1F*)transfer_factors_file->Get("transfer_factor_SBtoSA_halo");
  vector<pair<RooRealVar*, TH1*>>   halo_SR_below0p5_syst;
  makeConnectedBinList(procname_halo_SB, phopt, wspace, halo_SRr_aboveOverBelow_hist, halo_SR_below0p5_syst, halo_SR_above0p5_bins);
  
  // TH1F* hist_ZNuNuG_SR_below0p5_corrected = (TH1F*)ZnnGbelow0p5file->Get("histo_ZNuNuG"); // Only if not fitting
  // TH1F* hist_ZNuNuG_SR_below0p5_uncorrected = (TH1F*)ZnnGbelow0p5file->Get("histo_ZNuNuG_uncorrected"); // Only if not fitting
  // TH1F* hist_ZNuNuG_SR_below0p5_qcdscale = (TH1F*)ZnnGbelow0p5file->Get("histo_ZNuNuG_qcdscale"); // Only if not fitting
  // TH1F* hist_ZNuNuG_SR_below0p5_EWKUp = (TH1F*)hist_ZNuNuG_SR_below0p5_corrected->Clone("hist_ZNuNuG_SR_below0p5_EWKUp"); // Only if not fitting
  // TH1F* hist_ZNuNuG_SR_below0p5_EWKDown = (TH1F*)hist_ZNuNuG_SR_below0p5_corrected->Clone("hist_ZNuNuG_SR_below0p5_EWKDown"); // Only if not fitting
  // TH1F* hist_ZNuNuG_SR_below0p5_qcdscaleUp = (TH1F*)hist_ZNuNuG_SR_below0p5_corrected->Clone("hist_ZNuNuG_SR_below0p5_qcdscaleUp"); // Only if not fitting
  // TH1F* hist_ZNuNuG_SR_below0p5_qcdscaleDown = (TH1F*)hist_ZNuNuG_SR_below0p5_corrected->Clone("hist_ZNuNuG_SR_below0p5_qcdscaleDown"); // Only if not fitting
  // TH1F* hist_WG_SR_below0p5_corrected = (TH1F*)ZnnGbelow0p5file->Get("histo_WG"); // Only if not fitting
  // TH1F* hist_WG_SR_below0p5_uncorrected = (TH1F*)ZnnGbelow0p5file->Get("histo_WG_uncorrected"); // Only if not fitting
  // TH1F* hist_WG_SR_below0p5_qcdscale = (TH1F*)ZnnGbelow0p5file->Get("histo_WG_qcdscale"); // Only if not fitting
  // TH1F* hist_WG_SR_below0p5_EWKUp = (TH1F*)hist_WG_SR_below0p5_corrected->Clone("hist_WG_SR_below0p5_EWKUp"); // Only if not fitting
  // TH1F* hist_WG_SR_below0p5_EWKDown = (TH1F*)hist_WG_SR_below0p5_corrected->Clone("hist_WG_SR_below0p5_EWKDown"); // Only if not fitting
  // TH1F* hist_WG_SR_below0p5_qcdscaleUp = (TH1F*)hist_WG_SR_below0p5_corrected->Clone("hist_WG_SR_below0p5_qcdscaleUp"); // Only if not fitting
  // TH1F* hist_WG_SR_below0p5_qcdscaleDown = (TH1F*)hist_WG_SR_below0p5_corrected->Clone("hist_WG_SR_below0p5_qcdscaleDown"); // Only if not fitting
  // for(int i = 1; i <= nBins; i++){ // Only if not fitting
  //   Float_t diff = fabs(hist_ZNuNuG_SR_below0p5_corrected->GetBinContent(i) - hist_ZNuNuG_SR_below0p5_uncorrected->GetBinContent(i)); // Only if not fitting
  //   hist_ZNuNuG_SR_below0p5_EWKUp->SetBinContent(i, hist_ZNuNuG_SR_below0p5_corrected->GetBinContent(i) + diff); // Only if not fitting
  //   hist_ZNuNuG_SR_below0p5_EWKDown->SetBinContent(i, hist_ZNuNuG_SR_below0p5_corrected->GetBinContent(i) - diff); // Only if not fitting
  //   diff = fabs(hist_ZNuNuG_SR_below0p5_corrected->GetBinContent(i) - hist_ZNuNuG_SR_below0p5_qcdscale->GetBinContent(i)); // Only if not fitting
  //   hist_ZNuNuG_SR_below0p5_qcdscaleUp->SetBinContent(i, hist_ZNuNuG_SR_below0p5_corrected->GetBinContent(i) + diff); // Only if not fitting
  //   hist_ZNuNuG_SR_below0p5_qcdscaleDown->SetBinContent(i, hist_ZNuNuG_SR_below0p5_corrected->GetBinContent(i) - diff); // Only if not fitting
  //   diff = fabs(hist_WG_SR_below0p5_corrected->GetBinContent(i) - hist_WG_SR_below0p5_uncorrected->GetBinContent(i)); // Only if not fitting
  //   hist_WG_SR_below0p5_EWKUp->SetBinContent(i, hist_WG_SR_below0p5_corrected->GetBinContent(i) + diff); // Only if not fitting
  //   hist_WG_SR_below0p5_EWKDown->SetBinContent(i, hist_WG_SR_below0p5_corrected->GetBinContent(i) - diff); // Only if not fitting
  //   diff = fabs(hist_WG_SR_below0p5_corrected->GetBinContent(i) - hist_WG_SR_below0p5_qcdscale->GetBinContent(i)); // Only if not fitting
  //   hist_WG_SR_below0p5_qcdscaleUp->SetBinContent(i, hist_WG_SR_below0p5_corrected->GetBinContent(i) + diff); // Only if not fitting
  //   hist_WG_SR_below0p5_qcdscaleDown->SetBinContent(i, hist_WG_SR_below0p5_corrected->GetBinContent(i) - diff); // Only if not fitting
  // } // Only if not fitting
  
  // Data driven backgrounds
  addTemplate("QCD_SR_below0p5"                , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_jetfake"));
  addTemplate("QCD_SR_below0p5_QCDrUp"         , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_jetfake_errUp"));
  addTemplate("QCD_SR_below0p5_QCDrDown"       , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_jetfake_errDown"));
  addTemplate("Elefake_SR_below0p5"            , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_elefake"));
  // addTemplate("BHalo_SR_below0p5"              , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_bhalo"));
  addTemplate("Spikes_SR_below0p5"             , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_spikes"));
  // MC backgrounds
  // addTemplate("ZNuNuG_SR_below0p5"              , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_ZNuNuG")); // Only if not fitting
  // addTemplate("ZNuNuG_SR_below0p5_JESUp"        , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_ZNuNuG_JESUp")); // Only if not fitting
  // addTemplate("ZNuNuG_SR_below0p5_JESDown"      , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_ZNuNuG_JESDown")); // Only if not fitting
  // addTemplate("ZNuNuG_SR_below0p5_PESUp"        , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_ZNuNuG_PESUp")); // Only if not fitting
  // addTemplate("ZNuNuG_SR_below0p5_PESDown"      , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_ZNuNuG_PESDown")); // Only if not fitting
  // addTemplate("ZNuNuG_SR_below0p5_phoSFUp"        , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_ZNuNuG_phoSFUp")); // Only if not fitting
  // addTemplate("ZNuNuG_SR_below0p5_phoSFDown"      , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_ZNuNuG_phoSFDown")); // Only if not fitting
  // addTemplate("ZNuNuG_SR_below0p5_EWKUp"      , vars, wspace, hist_ZNuNuG_SR_below0p5_EWKUp); // Only if not fitting
  // addTemplate("ZNuNuG_SR_below0p5_EWKDown"      , vars, wspace, hist_ZNuNuG_SR_below0p5_EWKDown); // Only if not fitting
  // addTemplate("ZNuNuG_SR_below0p5_qcdscaleUp"      , vars, wspace, hist_ZNuNuG_SR_below0p5_qcdscaleUp); // Only if not fitting
  // addTemplate("ZNuNuG_SR_below0p5_qcdscaleDown"      , vars, wspace, hist_ZNuNuG_SR_below0p5_qcdscaleDown); // Only if not fitting
  // addTemplate("WG_SR_below0p5"              , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_WG")); // Only if not fitting
  // addTemplate("WG_SR_below0p5_JESUp"        , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_WG_JESUp")); // Only if not fitting
  // addTemplate("WG_SR_below0p5_JESDown"      , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_WG_JESDown")); // Only if not fitting
  // addTemplate("WG_SR_below0p5_PESUp"        , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_WG_PESUp")); // Only if not fitting
  // addTemplate("WG_SR_below0p5_PESDown"      , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_WG_PESDown")); // Only if not fitting
  // addTemplate("WG_SR_below0p5_phoSFUp"        , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_WG_phoSFUp")); // Only if not fitting
  // addTemplate("WG_SR_below0p5_phoSFDown"      , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_WG_phoSFDown")); // Only if not fitting
  // addTemplate("WG_SR_below0p5_EWKUp"      , vars, wspace, hist_WG_SR_below0p5_EWKUp); // Only if not fitting
  // addTemplate("WG_SR_below0p5_EWKDown"      , vars, wspace, hist_WG_SR_below0p5_EWKDown); // Only if not fitting
  // addTemplate("WG_SR_below0p5_qcdscaleUp"      , vars, wspace, hist_WG_SR_below0p5_qcdscaleUp); // Only if not fitting
  // addTemplate("WG_SR_below0p5_qcdscaleDown"      , vars, wspace, hist_WG_SR_below0p5_qcdscaleDown); // Only if not fitting
  addTemplate("GJets_SR_below0p5"              , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_GJets"));
  addTemplate("GJets_SR_below0p5_JESUp"        , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_GJets_JESUp"));
  addTemplate("GJets_SR_below0p5_JESDown"      , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_GJets_JESDown"));
  addTemplate("GJets_SR_below0p5_PESUp"        , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_GJets_PESUp"));
  addTemplate("GJets_SR_below0p5_PESDown"      , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_GJets_PESDown"));
  addTemplate("ZllG_SR_below0p5"               , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_ZllG_combined"));
  addTemplate("ZllG_SR_below0p5_JESUp"         , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_ZllG_JESUp_combined"));
  addTemplate("ZllG_SR_below0p5_JESDown"       , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_ZllG_JESDown_combined"));
  addTemplate("ZllG_SR_below0p5_PESUp"         , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_ZllG_PESUp_combined"));
  addTemplate("ZllG_SR_below0p5_PESDown"       , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_ZllG_PESDown_combined"));
  addTemplate("TTG_SR_below0p5"                , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_TTG"));
  addTemplate("TTG_SR_below0p5_JESUp"          , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_TTG_JESUp"));
  addTemplate("TTG_SR_below0p5_JESDown"        , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_TTG_JESDown"));
  addTemplate("TTG_SR_below0p5_PESUp"          , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_TTG_PESUp"));
  addTemplate("TTG_SR_below0p5_PESDown"        , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_TTG_PESDown"));
  addTemplate("TG_SR_below0p5"                 , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_TG"));
  addTemplate("TG_SR_below0p5_JESUp"           , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_TG_JESUp"));
  addTemplate("TG_SR_below0p5_JESDown"         , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_TG_JESDown"));
  addTemplate("TG_SR_below0p5_PESUp"           , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_TG_PESUp"));
  addTemplate("TG_SR_below0p5_PESDown"         , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_TG_PESDown"));
  addTemplate("Diphoton_SR_below0p5"           , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_diphoton"));
  addTemplate("Diphoton_SR_below0p5_JESUp"     , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_diphoton_JESUp"));
  addTemplate("Diphoton_SR_below0p5_JESDown"   , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_diphoton_JESDown"));
  addTemplate("Diphoton_SR_below0p5_PESUp"     , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_diphoton_PESUp"));
  addTemplate("Diphoton_SR_below0p5_PESDown"   , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_diphoton_PESDown"));
  addTemplate("WZ_SR_below0p5"                 , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_WZ"));
  addTemplate("WZ_SR_below0p5_JESUp"           , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_WZ_JESUp"));
  addTemplate("WZ_SR_below0p5_JESDown"         , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_WZ_JESDown"));
  addTemplate("WZ_SR_below0p5_PESUp"           , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_WZ_PESUp"));
  addTemplate("WZ_SR_below0p5_PESDown"         , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_WZ_PESDown"));
  addTemplate("ZZ_SR_below0p5"                 , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_ZZ"));
  addTemplate("ZZ_SR_below0p5_JESUp"           , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_ZZ_JESUp"));
  addTemplate("ZZ_SR_below0p5_JESDown"         , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_ZZ_JESDown"));
  addTemplate("ZZ_SR_below0p5_PESUp"           , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_ZZ_PESUp"));
  addTemplate("ZZ_SR_below0p5_PESDown"         , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_ZZ_PESDown"));
  addTemplate("WMuNu_SR_below0p5"              , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_WMuNu"));
  addTemplate("WMuNu_SR_below0p5_JESUp"        , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_WMuNu_JESUp"));
  addTemplate("WMuNu_SR_below0p5_JESDown"      , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_WMuNu_JESDown"));
  addTemplate("WMuNu_SR_below0p5_PESUp"        , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_WMuNu_PESUp"));
  addTemplate("WMuNu_SR_below0p5_PESDown"      , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_WMuNu_PESDown"));
  addTemplate("WTauNu_SR_below0p5"             , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_WTauNu"));
  addTemplate("WTauNu_SR_below0p5_JESUp"       , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_WTauNu_JESUp"));
  addTemplate("WTauNu_SR_below0p5_JESDown"     , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_WTauNu_JESDown"));
  addTemplate("WTauNu_SR_below0p5_PESUp"       , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_WTauNu_PESUp"));
  addTemplate("WTauNu_SR_below0p5_PESDown"     , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_WTauNu_PESDown"));
  addTemplate("WW_SR_below0p5"                 , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_WW"));
  addTemplate("WW_SR_below0p5_JESUp"           , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_WW_JESUp"));
  addTemplate("WW_SR_below0p5_JESDown"         , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_WW_JESDown"));
  addTemplate("WW_SR_below0p5_PESUp"           , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_WW_PESUp"));
  addTemplate("WW_SR_below0p5_PESDown"         , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_WW_PESDown"));

  // ---------------------------- CONTROL REGION (Dimuon) -----------------------------------------------------------------//
  string procname_ZM = "ZllG_ZM";
  
  addTemplate("data_obs_ZM", vars, wspace, (TH1F*)ZmmGfile->Get("data_obs"));
  
  TH1F* znng_ZMr_above0p5_hist = (TH1F*)transfer_factors_file->Get("transfer_factor_ZllGinZmmG_to_ZnnGinZnnG_above0p5_noShift");
  TH1F* znng_ZMr_above0p5_muEff_shiftUp = (TH1F*)transfer_factors_file->Get("transfer_factor_ZllGinZmmG_to_ZnnGinZnnG_above0p5_muEff_shiftUp");
  TH1F* znng_ZMr_above0p5_muEff_shiftDown = (TH1F*)transfer_factors_file->Get("transfer_factor_ZllGinZmmG_to_ZnnGinZnnG_above0p5_muEff_shiftDown");
  TH1F* znng_ZMr_above0p5_muEff_fractionalShifts = (TH1F*)znng_ZMr_above0p5_muEff_shiftUp->Clone(TString(procname_ZM+"_muEff"));
  for (int i = 1; i <= znng_ZMr_above0p5_muEff_fractionalShifts->GetNbinsX(); i++) {
    Float_t upshift = znng_ZMr_above0p5_muEff_shiftUp->GetBinContent(i)/znng_ZMr_above0p5_hist->GetBinContent(i) - 1.0;
    Float_t downshift = znng_ZMr_above0p5_muEff_shiftDown->GetBinContent(i)/znng_ZMr_above0p5_hist->GetBinContent(i) - 1.0;
    Float_t shiftEnvelope = TMath::Max(fabs(upshift), fabs(downshift));
    znng_ZMr_above0p5_muEff_fractionalShifts->SetBinContent(i, shiftEnvelope);
  }
  vector<pair<RooRealVar*, TH1*>> znng_ZM_above0p5_syst;
  RooRealVar* znng_ZM_above0p5_muEff = new RooRealVar("ZNuNuGoverZLLG_ZM_above0p5_muEff", "", 0., -5., 5.);
  // znng_ZM_above0p5_syst.push_back(pair<RooRealVar*, TH1*>(znng_ZM_above0p5_muEff, znng_ZMr_above0p5_muEff_fractionalShifts));
  if (!connectWZ || connectWZ) makeConnectedBinList(procname_ZM, phopt, wspace, znng_ZMr_above0p5_hist, znng_ZM_above0p5_syst, znng_SR_above0p5_bins);
  
  // Data driven backgrounds
  addTemplate("QCD_ZM"              , vars, wspace, (TH1F*)ZmmGfile->Get("histo_jetfake"));
  addTemplate("QCD_ZM_QCDrUp"       , vars, wspace, (TH1F*)ZmmGfile->Get("histo_jetfake_errUp"));
  addTemplate("QCD_ZM_QCDrDown"     , vars, wspace, (TH1F*)ZmmGfile->Get("histo_jetfake_errDown"));
  // MC backgrounds
  addTemplate("TTG_ZM"              , vars, wspace, (TH1F*)ZmmGfile->Get("histo_TTG"));
  addTemplate("TTG_ZM_JESUp"        , vars, wspace, (TH1F*)ZmmGfile->Get("histo_TTG_JESUp"));
  addTemplate("TTG_ZM_JESDown"      , vars, wspace, (TH1F*)ZmmGfile->Get("histo_TTG_JESDown"));
  addTemplate("TTG_ZM_PESUp"        , vars, wspace, (TH1F*)ZmmGfile->Get("histo_TTG_PESUp"));
  addTemplate("TTG_ZM_PESDown"      , vars, wspace, (TH1F*)ZmmGfile->Get("histo_TTG_PESDown"));
  addTemplate("WZ_ZM"               , vars, wspace, (TH1F*)ZmmGfile->Get("histo_WZ"));
  addTemplate("WZ_ZM_JESUp"         , vars, wspace, (TH1F*)ZmmGfile->Get("histo_WZ_JESUp"));
  addTemplate("WZ_ZM_JESDown"       , vars, wspace, (TH1F*)ZmmGfile->Get("histo_WZ_JESDown"));
  addTemplate("WZ_ZM_PESUp"         , vars, wspace, (TH1F*)ZmmGfile->Get("histo_WZ_PESUp"));
  addTemplate("WZ_ZM_PESDown"       , vars, wspace, (TH1F*)ZmmGfile->Get("histo_WZ_PESDown"));

  // ---------------------------- CONTROL REGION (Dielectron) -----------------------------------------------------------------//
  string procname_ZE = "ZllG_ZE";
  
  addTemplate("data_obs_ZE", vars, wspace, (TH1F*)ZeeGfile->Get("data_obs"));
  
  TH1F* znng_ZEr_above0p5_hist = (TH1F*)transfer_factors_file->Get("transfer_factor_ZllGinZeeG_to_ZnnGinZnnG_above0p5_noShift");
  TH1F* znng_ZEr_above0p5_eleEff_shiftUp = (TH1F*)transfer_factors_file->Get("transfer_factor_ZllGinZeeG_to_ZnnGinZnnG_above0p5_eleEff_shiftUp");
  TH1F* znng_ZEr_above0p5_eleEff_shiftDown = (TH1F*)transfer_factors_file->Get("transfer_factor_ZllGinZeeG_to_ZnnGinZnnG_above0p5_eleEff_shiftDown");
  TH1F* znng_ZEr_above0p5_eleEff_fractionalShifts = (TH1F*)znng_ZEr_above0p5_eleEff_shiftUp->Clone(TString(procname_ZE+"_eleEff"));
  for (int i = 1; i <= znng_ZEr_above0p5_eleEff_fractionalShifts->GetNbinsX(); i++) {
    Float_t upshift = znng_ZEr_above0p5_eleEff_shiftUp->GetBinContent(i)/znng_ZEr_above0p5_hist->GetBinContent(i) - 1.0;
    Float_t downshift = znng_ZEr_above0p5_eleEff_shiftDown->GetBinContent(i)/znng_ZEr_above0p5_hist->GetBinContent(i) - 1.0;
    Float_t shiftEnvelope = TMath::Max(fabs(upshift), fabs(downshift));
    znng_ZEr_above0p5_eleEff_fractionalShifts->SetBinContent(i, shiftEnvelope);
  }
  vector<pair<RooRealVar*, TH1*>> znng_ZE_above0p5_syst;
  RooRealVar* znng_ZE_above0p5_eleEff = new RooRealVar("ZNuNuGoverZLLG_ZE_above0p5_eleEff", "", 0., -5., 5.);
  // znng_ZE_above0p5_syst.push_back(pair<RooRealVar*, TH1*>(znng_ZE_above0p5_eleEff, znng_ZEr_above0p5_eleEff_fractionalShifts));
  if (!connectWZ || connectWZ) makeConnectedBinList(procname_ZE, phopt, wspace, znng_ZEr_above0p5_hist, znng_ZE_above0p5_syst, znng_SR_above0p5_bins);
  
  // Data driven backgrounds
  addTemplate("QCD_ZE"              , vars, wspace, (TH1F*)ZeeGfile->Get("histo_jetfake"));
  addTemplate("QCD_ZE_QCDrUp"       , vars, wspace, (TH1F*)ZeeGfile->Get("histo_jetfake_errUp"));
  addTemplate("QCD_ZE_QCDrDown"     , vars, wspace, (TH1F*)ZeeGfile->Get("histo_jetfake_errDown"));
  // MC backgrounds
  addTemplate("TTG_ZE"              , vars, wspace, (TH1F*)ZeeGfile->Get("histo_TTG"));
  addTemplate("TTG_ZE_JESUp"        , vars, wspace, (TH1F*)ZeeGfile->Get("histo_TTG_JESUp"));
  addTemplate("TTG_ZE_JESDown"      , vars, wspace, (TH1F*)ZeeGfile->Get("histo_TTG_JESDown"));
  addTemplate("TTG_ZE_PESUp"        , vars, wspace, (TH1F*)ZeeGfile->Get("histo_TTG_PESUp"));
  addTemplate("TTG_ZE_PESDown"      , vars, wspace, (TH1F*)ZeeGfile->Get("histo_TTG_PESDown"));
  addTemplate("WZ_ZE"               , vars, wspace, (TH1F*)ZeeGfile->Get("histo_WZ"));
  addTemplate("WZ_ZE_JESUp"         , vars, wspace, (TH1F*)ZeeGfile->Get("histo_WZ_JESUp"));
  addTemplate("WZ_ZE_JESDown"       , vars, wspace, (TH1F*)ZeeGfile->Get("histo_WZ_JESDown"));
  addTemplate("WZ_ZE_PESUp"         , vars, wspace, (TH1F*)ZeeGfile->Get("histo_WZ_PESUp"));
  addTemplate("WZ_ZE_PESDown"       , vars, wspace, (TH1F*)ZeeGfile->Get("histo_WZ_PESDown"));

  // ---------------------------- CONTROL REGION (Single muon) -----------------------------------------------------------------//
  string procname_WM = "WG_above0p5_WM";
  
  addTemplate("data_obs_WM"  , vars, wspace, (TH1F*)WmnGfile->Get("data_obs"));

  // Without Z/W link
  TH1F* wlng_WMr_above0p5_hist = (TH1F*)transfer_factors_file->Get("transfer_factor_WGinWmnG_to_WGinZnnG_above0p5_noShift");
  TH1F* wlng_WMr_above0p5_muEff_shiftUp = (TH1F*)transfer_factors_file->Get("transfer_factor_WGinWmnG_to_WGinZnnG_above0p5_muEff_shiftUp");
  TH1F* wlng_WMr_above0p5_muEff_shiftDown = (TH1F*)transfer_factors_file->Get("transfer_factor_WGinWmnG_to_WGinZnnG_above0p5_muEff_shiftDown");
  TH1F* wlng_WMr_above0p5_muEff_fractionalShifts = (TH1F*)wlng_WMr_above0p5_muEff_shiftUp->Clone(TString(procname_WM+"_muEff"));
  for (int i = 1; i <= wlng_WMr_above0p5_muEff_fractionalShifts->GetNbinsX(); i++) {
    Float_t upshift = wlng_WMr_above0p5_muEff_shiftUp->GetBinContent(i)/wlng_WMr_above0p5_hist->GetBinContent(i) - 1.0;
    Float_t downshift = wlng_WMr_above0p5_muEff_shiftDown->GetBinContent(i)/wlng_WMr_above0p5_hist->GetBinContent(i) - 1.0;
    Float_t shiftEnvelope = TMath::Max(fabs(upshift), fabs(downshift));
    wlng_WMr_above0p5_muEff_fractionalShifts->SetBinContent(i, shiftEnvelope);
  }
  vector<pair<RooRealVar*, TH1*>> wlng_WM_above0p5_syst;
  RooRealVar* wlng_WM_above0p5_muEff = new RooRealVar("WLNuG_WM_above0p5_muEff", "", 0., -5., 5.);
  // wlng_WM_above0p5_syst.push_back(pair<RooRealVar*, TH1*>(wlng_WM_above0p5_muEff, wlng_WMr_above0p5_muEff_fractionalShifts));
  
  makeConnectedBinList(procname_WM, phopt, wspace, wlng_WMr_above0p5_hist, wlng_WM_above0p5_syst, wlng_SR_above0p5_bins);
  
  // Data driven backgrounds
  addTemplate("QCD_WM"              , vars, wspace, (TH1F*)WmnGfile->Get("histo_jetfake"));
  addTemplate("QCD_WM_QCDrUp"       , vars, wspace, (TH1F*)WmnGfile->Get("histo_jetfake_errUp"));
  addTemplate("QCD_WM_QCDrDown"     , vars, wspace, (TH1F*)WmnGfile->Get("histo_jetfake_errDown"));
  // MC backgrounds
  addTemplate("ZllG_WM"             , vars, wspace, (TH1F*)WmnGfile->Get("histo_ZllG_combined"));
  addTemplate("ZllG_WM_JESUp"       , vars, wspace, (TH1F*)WmnGfile->Get("histo_ZllG_JESUp_combined"));
  addTemplate("ZllG_WM_JESDown"     , vars, wspace, (TH1F*)WmnGfile->Get("histo_ZllG_JESDown_combined"));
  addTemplate("ZllG_WM_PESUp"       , vars, wspace, (TH1F*)WmnGfile->Get("histo_ZllG_PESUp_combined"));
  addTemplate("ZllG_WM_PESDown"     , vars, wspace, (TH1F*)WmnGfile->Get("histo_ZllG_PESDown_combined"));
  addTemplate("TTG_WM"              , vars, wspace, (TH1F*)WmnGfile->Get("histo_TTG"));
  addTemplate("TTG_WM_JESUp"        , vars, wspace, (TH1F*)WmnGfile->Get("histo_TTG_JESUp"));
  addTemplate("TTG_WM_JESDown"      , vars, wspace, (TH1F*)WmnGfile->Get("histo_TTG_JESDown"));
  addTemplate("TTG_WM_PESUp"        , vars, wspace, (TH1F*)WmnGfile->Get("histo_TTG_PESUp"));
  addTemplate("TTG_WM_PESDown"      , vars, wspace, (TH1F*)WmnGfile->Get("histo_TTG_PESDown"));
  addTemplate("TG_WM"               , vars, wspace, (TH1F*)WmnGfile->Get("histo_TG"));
  addTemplate("TG_WM_JESUp"         , vars, wspace, (TH1F*)WmnGfile->Get("histo_TG_JESUp"));
  addTemplate("TG_WM_JESDown"       , vars, wspace, (TH1F*)WmnGfile->Get("histo_TG_JESDown"));
  addTemplate("TG_WM_PESUp"         , vars, wspace, (TH1F*)WmnGfile->Get("histo_TG_PESUp"));
  addTemplate("TG_WM_PESDown"       , vars, wspace, (TH1F*)WmnGfile->Get("histo_TG_PESDown"));
  addTemplate("Diphoton_WM"         , vars, wspace, (TH1F*)WmnGfile->Get("histo_diphoton"));
  addTemplate("Diphoton_WM_JESUp"   , vars, wspace, (TH1F*)WmnGfile->Get("histo_diphoton_JESUp"));
  addTemplate("Diphoton_WM_JESDown" , vars, wspace, (TH1F*)WmnGfile->Get("histo_diphoton_JESDown"));
  addTemplate("Diphoton_WM_PESUp"   , vars, wspace, (TH1F*)WmnGfile->Get("histo_diphoton_PESUp"));
  addTemplate("Diphoton_WM_PESDown" , vars, wspace, (TH1F*)WmnGfile->Get("histo_diphoton_PESDown"));
  addTemplate("WZ_WM"               , vars, wspace, (TH1F*)WmnGfile->Get("histo_WZ"));
  addTemplate("WZ_WM_JESUp"         , vars, wspace, (TH1F*)WmnGfile->Get("histo_WZ_JESUp"));
  addTemplate("WZ_WM_JESDown"       , vars, wspace, (TH1F*)WmnGfile->Get("histo_WZ_JESDown"));
  addTemplate("WZ_WM_PESUp"         , vars, wspace, (TH1F*)WmnGfile->Get("histo_WZ_PESUp"));
  addTemplate("WZ_WM_PESDown"       , vars, wspace, (TH1F*)WmnGfile->Get("histo_WZ_PESDown"));
  addTemplate("WW_WM"               , vars, wspace, (TH1F*)WmnGfile->Get("histo_WW"));
  addTemplate("WW_WM_JESUp"         , vars, wspace, (TH1F*)WmnGfile->Get("histo_WW_JESUp"));
  addTemplate("WW_WM_JESDown"       , vars, wspace, (TH1F*)WmnGfile->Get("histo_WW_JESDown"));
  addTemplate("WW_WM_PESUp"         , vars, wspace, (TH1F*)WmnGfile->Get("histo_WW_PESUp"));
  addTemplate("WW_WM_PESDown"       , vars, wspace, (TH1F*)WmnGfile->Get("histo_WW_PESDown"));
  
  // ---------------------------- CONTROL REGION (Single electron) -----------------------------------------------------------------//
  string procname_WE = "WG_above0p5_WE";
  
  addTemplate("data_obs_WE"  , vars, wspace, (TH1F*)WenGfile->Get("data_obs"));

  // Without Z/W link
  TH1F* wlng_WEr_above0p5_hist = (TH1F*)transfer_factors_file->Get("transfer_factor_WGinWenG_to_WGinZnnG_above0p5_noShift");
  TH1F* wlng_WEr_above0p5_eleEff_shiftUp = (TH1F*)transfer_factors_file->Get("transfer_factor_WGinWenG_to_WGinZnnG_above0p5_eleEff_shiftUp");
  TH1F* wlng_WEr_above0p5_eleEff_shiftDown = (TH1F*)transfer_factors_file->Get("transfer_factor_WGinWenG_to_WGinZnnG_above0p5_eleEff_shiftDown");
  TH1F* wlng_WEr_above0p5_eleEff_fractionalShifts = (TH1F*)wlng_WEr_above0p5_eleEff_shiftUp->Clone(TString(procname_WE+"_eleEff"));
  for (int i = 1; i <= wlng_WEr_above0p5_eleEff_fractionalShifts->GetNbinsX(); i++) {
    Float_t upshift = wlng_WEr_above0p5_eleEff_shiftUp->GetBinContent(i)/wlng_WEr_above0p5_hist->GetBinContent(i) - 1.0;
    Float_t downshift = wlng_WEr_above0p5_eleEff_shiftDown->GetBinContent(i)/wlng_WEr_above0p5_hist->GetBinContent(i) - 1.0;
    Float_t shiftEnvelope = TMath::Max(fabs(upshift), fabs(downshift));
    wlng_WEr_above0p5_eleEff_fractionalShifts->SetBinContent(i, shiftEnvelope);
  }
  vector<pair<RooRealVar*, TH1*>> wlng_WE_above0p5_syst;
  RooRealVar* wlng_WE_above0p5_eleEff = new RooRealVar("WLNuG_WE_above0p5_eleEff", "", 0., -5., 5.);
  // wlng_WE_above0p5_syst.push_back(pair<RooRealVar*, TH1*>(wlng_WE_above0p5_eleEff, wlng_WEr_above0p5_eleEff_fractionalShifts));

  makeConnectedBinList(procname_WE, phopt, wspace, wlng_WEr_above0p5_hist, wlng_WE_above0p5_syst, wlng_SR_above0p5_bins);
  
  // Data driven backgrounds
  addTemplate("QCD_WE"              , vars, wspace, (TH1F*)WenGfile->Get("histo_jetfake"));
  addTemplate("QCD_WE_QCDrUp"       , vars, wspace, (TH1F*)WenGfile->Get("histo_jetfake_errUp"));
  addTemplate("QCD_WE_QCDrDown"     , vars, wspace, (TH1F*)WenGfile->Get("histo_jetfake_errDown"));
  addTemplate("Elefake_WE"          , vars, wspace, (TH1F*)WenGfile->Get("histo_elefake"));
  // MC backgrounds
  addTemplate("ZllG_WE"             , vars, wspace, (TH1F*)WenGfile->Get("histo_ZllG_combined"));
  addTemplate("ZllG_WE_JESUp"       , vars, wspace, (TH1F*)WenGfile->Get("histo_ZllG_JESUp_combined"));
  addTemplate("ZllG_WE_JESDown"     , vars, wspace, (TH1F*)WenGfile->Get("histo_ZllG_JESDown_combined"));
  addTemplate("ZllG_WE_PESUp"       , vars, wspace, (TH1F*)WenGfile->Get("histo_ZllG_PESUp_combined"));
  addTemplate("ZllG_WE_PESDown"     , vars, wspace, (TH1F*)WenGfile->Get("histo_ZllG_PESDown_combined"));
  addTemplate("TTG_WE"              , vars, wspace, (TH1F*)WenGfile->Get("histo_TTG"));
  addTemplate("TTG_WE_JESUp"        , vars, wspace, (TH1F*)WenGfile->Get("histo_TTG_JESUp"));
  addTemplate("TTG_WE_JESDown"      , vars, wspace, (TH1F*)WenGfile->Get("histo_TTG_JESDown"));
  addTemplate("TTG_WE_PESUp"        , vars, wspace, (TH1F*)WenGfile->Get("histo_TTG_PESUp"));
  addTemplate("TTG_WE_PESDown"      , vars, wspace, (TH1F*)WenGfile->Get("histo_TTG_PESDown"));
  addTemplate("TG_WE"               , vars, wspace, (TH1F*)WenGfile->Get("histo_TG"));
  addTemplate("TG_WE_JESUp"         , vars, wspace, (TH1F*)WenGfile->Get("histo_TG_JESUp"));
  addTemplate("TG_WE_JESDown"       , vars, wspace, (TH1F*)WenGfile->Get("histo_TG_JESDown"));
  addTemplate("TG_WE_PESUp"         , vars, wspace, (TH1F*)WenGfile->Get("histo_TG_PESUp"));
  addTemplate("TG_WE_PESDown"       , vars, wspace, (TH1F*)WenGfile->Get("histo_TG_PESDown"));
  addTemplate("Diphoton_WE"         , vars, wspace, (TH1F*)WenGfile->Get("histo_diphoton"));
  addTemplate("Diphoton_WE_JESUp"   , vars, wspace, (TH1F*)WenGfile->Get("histo_diphoton_JESUp"));
  addTemplate("Diphoton_WE_JESDown" , vars, wspace, (TH1F*)WenGfile->Get("histo_diphoton_JESDown"));
  addTemplate("Diphoton_WE_PESUp"   , vars, wspace, (TH1F*)WenGfile->Get("histo_diphoton_PESUp"));
  addTemplate("Diphoton_WE_PESDown" , vars, wspace, (TH1F*)WenGfile->Get("histo_diphoton_PESDown"));
  addTemplate("WZ_WE"               , vars, wspace, (TH1F*)WenGfile->Get("histo_WZ"));
  addTemplate("WZ_WE_JESUp"         , vars, wspace, (TH1F*)WenGfile->Get("histo_WZ_JESUp"));
  addTemplate("WZ_WE_JESDown"       , vars, wspace, (TH1F*)WenGfile->Get("histo_WZ_JESDown"));
  addTemplate("WZ_WE_PESUp"         , vars, wspace, (TH1F*)WenGfile->Get("histo_WZ_PESUp"));
  addTemplate("WZ_WE_PESDown"       , vars, wspace, (TH1F*)WenGfile->Get("histo_WZ_PESDown"));
  addTemplate("WW_WE"               , vars, wspace, (TH1F*)WenGfile->Get("histo_WW"));
  addTemplate("WW_WE_JESUp"         , vars, wspace, (TH1F*)WenGfile->Get("histo_WW_JESUp"));
  addTemplate("WW_WE_JESDown"       , vars, wspace, (TH1F*)WenGfile->Get("histo_WW_JESDown"));
  addTemplate("WW_WE_PESUp"         , vars, wspace, (TH1F*)WenGfile->Get("histo_WW_PESUp"));
  addTemplate("WW_WE_PESDown"       , vars, wspace, (TH1F*)WenGfile->Get("histo_WW_PESDown"));
  
  //Statistical errors
  for(int i = 1; i <= nBins; i++){
    char binChar[10];
    sprintf(binChar, "%d", i);
    std::string binNumber(binChar);
    
    
// Primary
    // TH1F* ZNuNuG_SR_above0p5_histBinUp = (TH1F*)((TH1F*)ZnnGabove0p5file->Get("histo_ZNuNuG"))->Clone("ZNuNuG_SR_above0p5_histBinUp"); // Only if not fitting
    // TH1F* ZNuNuG_SR_above0p5_histBinDown = (TH1F*)((TH1F*)ZnnGabove0p5file->Get("histo_ZNuNuG"))->Clone("ZNuNuG_SR_above0p5_histBinDown"); // Only if not fitting
    // ZNuNuG_SR_above0p5_histBinUp->SetBinContent(i, ZNuNuG_SR_above0p5_histBinUp->GetBinContent(i)+ZNuNuG_SR_above0p5_histBinUp->GetBinError(i)); // Only if not fitting
    // ZNuNuG_SR_above0p5_histBinDown->SetBinContent(i, TMath::Max(0.01*ZNuNuG_SR_above0p5_histBinDown->GetBinContent(i), ZNuNuG_SR_above0p5_histBinDown->GetBinContent(i)-ZNuNuG_SR_above0p5_histBinDown->GetBinError(i))); // Only if not fitting
    // addTemplate("ZNuNuG_SR_above0p5_ZNuNuGSignalSBin"+binNumber+"Up", vars, wspace, ZNuNuG_SR_above0p5_histBinUp); // Only if not fitting
    // addTemplate("ZNuNuG_SR_above0p5_ZNuNuGSignalSBin"+binNumber+"Down", vars, wspace, ZNuNuG_SR_above0p5_histBinDown); // Only if not fitting
    // TH1F* WG_SR_above0p5_histBinUp = (TH1F*)((TH1F*)ZnnGabove0p5file->Get("histo_WG"))->Clone("WG_SR_above0p5_histBinUp"); // Only if not fitting
    // TH1F* WG_SR_above0p5_histBinDown = (TH1F*)((TH1F*)ZnnGabove0p5file->Get("histo_WG"))->Clone("WG_SR_above0p5_histBinDown"); // Only if not fitting
    // WG_SR_above0p5_histBinUp->SetBinContent(i, WG_SR_above0p5_histBinUp->GetBinContent(i)+WG_SR_above0p5_histBinUp->GetBinError(i)); // Only if not fitting
    // WG_SR_above0p5_histBinDown->SetBinContent(i, TMath::Max(0.01*WG_SR_above0p5_histBinDown->GetBinContent(i), WG_SR_above0p5_histBinDown->GetBinContent(i)-WG_SR_above0p5_histBinDown->GetBinError(i))); // Only if not fitting
    // addTemplate("WG_SR_above0p5_WGSignalSBin"+binNumber+"Up", vars, wspace, WG_SR_above0p5_histBinUp); // Only if not fitting
    // addTemplate("WG_SR_above0p5_WGSignalSBin"+binNumber+"Down", vars, wspace, WG_SR_above0p5_histBinDown); // Only if not fitting
    // TH1F* ZNuNuG_SR_below0p5_histBinUp = (TH1F*)((TH1F*)ZnnGbelow0p5file->Get("histo_ZNuNuG"))->Clone("ZNuNuG_SR_below0p5_histBinUp"); // Only if not fitting
    // TH1F* ZNuNuG_SR_below0p5_histBinDown = (TH1F*)((TH1F*)ZnnGbelow0p5file->Get("histo_ZNuNuG"))->Clone("ZNuNuG_SR_below0p5_histBinDown"); // Only if not fitting
    // ZNuNuG_SR_below0p5_histBinUp->SetBinContent(i, ZNuNuG_SR_below0p5_histBinUp->GetBinContent(i)+ZNuNuG_SR_below0p5_histBinUp->GetBinError(i)); // Only if not fitting
    // ZNuNuG_SR_below0p5_histBinDown->SetBinContent(i, TMath::Max(0.01*ZNuNuG_SR_below0p5_histBinDown->GetBinContent(i), ZNuNuG_SR_below0p5_histBinDown->GetBinContent(i)-ZNuNuG_SR_below0p5_histBinDown->GetBinError(i))); // Only if not fitting
    // addTemplate("ZNuNuG_SR_below0p5_ZNuNuGSignalSBin"+binNumber+"Up", vars, wspace, ZNuNuG_SR_below0p5_histBinUp); // Only if not fitting
    // addTemplate("ZNuNuG_SR_below0p5_ZNuNuGSignalSBin"+binNumber+"Down", vars, wspace, ZNuNuG_SR_below0p5_histBinDown); // Only if not fitting
    // TH1F* WG_SR_below0p5_histBinUp = (TH1F*)((TH1F*)ZnnGbelow0p5file->Get("histo_WG"))->Clone("WG_SR_below0p5_histBinUp"); // Only if not fitting
    // TH1F* WG_SR_below0p5_histBinDown = (TH1F*)((TH1F*)ZnnGbelow0p5file->Get("histo_WG"))->Clone("WG_SR_below0p5_histBinDown"); // Only if not fitting
    // WG_SR_below0p5_histBinUp->SetBinContent(i, WG_SR_below0p5_histBinUp->GetBinContent(i)+WG_SR_below0p5_histBinUp->GetBinError(i)); // Only if not fitting
    // WG_SR_below0p5_histBinDown->SetBinContent(i, TMath::Max(0.01*WG_SR_below0p5_histBinDown->GetBinContent(i), WG_SR_below0p5_histBinDown->GetBinContent(i)-WG_SR_below0p5_histBinDown->GetBinError(i))); // Only if not fitting
    // addTemplate("WG_SR_below0p5_WGSignalSBin"+binNumber+"Up", vars, wspace, WG_SR_below0p5_histBinUp); // Only if not fitting
    // addTemplate("WG_SR_below0p5_WGSignalSBin"+binNumber+"Down", vars, wspace, WG_SR_below0p5_histBinDown); // Only if not fitting
    
    TH1F* GJets_SR_above0p5_histBinUp = (TH1F*)((TH1F*)ZnnGabove0p5file->Get("histo_GJets"))->Clone("GJets_SR_above0p5_histBinUp");
    TH1F* GJets_SR_above0p5_histBinDown = (TH1F*)((TH1F*)ZnnGabove0p5file->Get("histo_GJets"))->Clone("GJets_SR_above0p5_histBinDown");
    GJets_SR_above0p5_histBinUp->SetBinContent(i, GJets_SR_above0p5_histBinUp->GetBinContent(i)+GJets_SR_above0p5_histBinUp->GetBinError(i));
    GJets_SR_above0p5_histBinDown->SetBinContent(i, TMath::Max(0.01*GJets_SR_above0p5_histBinDown->GetBinContent(i), GJets_SR_above0p5_histBinDown->GetBinContent(i)-GJets_SR_above0p5_histBinDown->GetBinError(i)));
    addTemplate("GJets_SR_above0p5_GJetsSignalSBin"+binNumber+"Up", vars, wspace, GJets_SR_above0p5_histBinUp);
    addTemplate("GJets_SR_above0p5_GJetsSignalSBin"+binNumber+"Down", vars, wspace, GJets_SR_above0p5_histBinDown);
    TH1F* WZ_SR_above0p5_histBinUp = (TH1F*)((TH1F*)ZnnGabove0p5file->Get("histo_WZ"))->Clone("WZ_SR_above0p5_histBinUp");
    TH1F* WZ_SR_above0p5_histBinDown = (TH1F*)((TH1F*)ZnnGabove0p5file->Get("histo_WZ"))->Clone("WZ_SR_above0p5_histBinDown");
    WZ_SR_above0p5_histBinUp->SetBinContent(i, WZ_SR_above0p5_histBinUp->GetBinContent(i)+WZ_SR_above0p5_histBinUp->GetBinError(i));
    WZ_SR_above0p5_histBinDown->SetBinContent(i, TMath::Max(0.01*WZ_SR_above0p5_histBinDown->GetBinContent(i), WZ_SR_above0p5_histBinDown->GetBinContent(i)-WZ_SR_above0p5_histBinDown->GetBinError(i)));
    addTemplate("WZ_SR_above0p5_WZSignalSBin"+binNumber+"Up", vars, wspace, WZ_SR_above0p5_histBinUp);
    addTemplate("WZ_SR_above0p5_WZSignalSBin"+binNumber+"Down", vars, wspace, WZ_SR_above0p5_histBinDown);
    TH1F* WMuNu_SR_above0p5_histBinUp = (TH1F*)((TH1F*)ZnnGabove0p5file->Get("histo_WMuNu"))->Clone("WMuNu_SR_above0p5_histBinUp");
    TH1F* WMuNu_SR_above0p5_histBinDown = (TH1F*)((TH1F*)ZnnGabove0p5file->Get("histo_WMuNu"))->Clone("WMuNu_SR_above0p5_histBinDown");
    WMuNu_SR_above0p5_histBinUp->SetBinContent(i, WMuNu_SR_above0p5_histBinUp->GetBinContent(i)+WMuNu_SR_above0p5_histBinUp->GetBinError(i));
    WMuNu_SR_above0p5_histBinDown->SetBinContent(i, TMath::Max(0.01*WMuNu_SR_above0p5_histBinDown->GetBinContent(i), WMuNu_SR_above0p5_histBinDown->GetBinContent(i)-WMuNu_SR_above0p5_histBinDown->GetBinError(i)));
    addTemplate("WMuNu_SR_above0p5_WMuNuSignalSBin"+binNumber+"Up", vars, wspace, WMuNu_SR_above0p5_histBinUp);
    addTemplate("WMuNu_SR_above0p5_WMuNuSignalSBin"+binNumber+"Down", vars, wspace, WMuNu_SR_above0p5_histBinDown);
    TH1F* WTauNu_SR_above0p5_histBinUp = (TH1F*)((TH1F*)ZnnGabove0p5file->Get("histo_WTauNu"))->Clone("WTauNu_SR_above0p5_histBinUp");
    TH1F* WTauNu_SR_above0p5_histBinDown = (TH1F*)((TH1F*)ZnnGabove0p5file->Get("histo_WTauNu"))->Clone("WTauNu_SR_above0p5_histBinDown");
    WTauNu_SR_above0p5_histBinUp->SetBinContent(i, WTauNu_SR_above0p5_histBinUp->GetBinContent(i)+WTauNu_SR_above0p5_histBinUp->GetBinError(i));
    WTauNu_SR_above0p5_histBinDown->SetBinContent(i, TMath::Max(0.01*WTauNu_SR_above0p5_histBinDown->GetBinContent(i), WTauNu_SR_above0p5_histBinDown->GetBinContent(i)-WTauNu_SR_above0p5_histBinDown->GetBinError(i)));
    addTemplate("WTauNu_SR_above0p5_WTauNuSignalSBin"+binNumber+"Up", vars, wspace, WTauNu_SR_above0p5_histBinUp);
    addTemplate("WTauNu_SR_above0p5_WTauNuSignalSBin"+binNumber+"Down", vars, wspace, WTauNu_SR_above0p5_histBinDown);
    TH1F* GJets_SR_below0p5_histBinUp = (TH1F*)((TH1F*)ZnnGbelow0p5file->Get("histo_GJets"))->Clone("GJets_SR_below0p5_histBinUp");
    TH1F* GJets_SR_below0p5_histBinDown = (TH1F*)((TH1F*)ZnnGbelow0p5file->Get("histo_GJets"))->Clone("GJets_SR_below0p5_histBinDown");
    GJets_SR_below0p5_histBinUp->SetBinContent(i, GJets_SR_below0p5_histBinUp->GetBinContent(i)+GJets_SR_below0p5_histBinUp->GetBinError(i));
    GJets_SR_below0p5_histBinDown->SetBinContent(i, TMath::Max(0.01*GJets_SR_below0p5_histBinDown->GetBinContent(i), GJets_SR_below0p5_histBinDown->GetBinContent(i)-GJets_SR_below0p5_histBinDown->GetBinError(i)));
    addTemplate("GJets_SR_below0p5_GJetsSignalSBin"+binNumber+"Up", vars, wspace, GJets_SR_below0p5_histBinUp);
    addTemplate("GJets_SR_below0p5_GJetsSignalSBin"+binNumber+"Down", vars, wspace, GJets_SR_below0p5_histBinDown);
    TH1F* WZ_SR_below0p5_histBinUp = (TH1F*)((TH1F*)ZnnGbelow0p5file->Get("histo_WZ"))->Clone("WZ_SR_below0p5_histBinUp");
    TH1F* WZ_SR_below0p5_histBinDown = (TH1F*)((TH1F*)ZnnGbelow0p5file->Get("histo_WZ"))->Clone("WZ_SR_below0p5_histBinDown");
    WZ_SR_below0p5_histBinUp->SetBinContent(i, WZ_SR_below0p5_histBinUp->GetBinContent(i)+WZ_SR_below0p5_histBinUp->GetBinError(i));
    WZ_SR_below0p5_histBinDown->SetBinContent(i, TMath::Max(0.01*WZ_SR_below0p5_histBinDown->GetBinContent(i), WZ_SR_below0p5_histBinDown->GetBinContent(i)-WZ_SR_below0p5_histBinDown->GetBinError(i)));
    addTemplate("WZ_SR_below0p5_WZSignalSBin"+binNumber+"Up", vars, wspace, WZ_SR_below0p5_histBinUp);
    addTemplate("WZ_SR_below0p5_WZSignalSBin"+binNumber+"Down", vars, wspace, WZ_SR_below0p5_histBinDown);
    TH1F* WMuNu_SR_below0p5_histBinUp = (TH1F*)((TH1F*)ZnnGbelow0p5file->Get("histo_WMuNu"))->Clone("WMuNu_SR_below0p5_histBinUp");
    TH1F* WMuNu_SR_below0p5_histBinDown = (TH1F*)((TH1F*)ZnnGbelow0p5file->Get("histo_WMuNu"))->Clone("WMuNu_SR_below0p5_histBinDown");
    WMuNu_SR_below0p5_histBinUp->SetBinContent(i, WMuNu_SR_below0p5_histBinUp->GetBinContent(i)+WMuNu_SR_below0p5_histBinUp->GetBinError(i));
    WMuNu_SR_below0p5_histBinDown->SetBinContent(i, TMath::Max(0.01*WMuNu_SR_below0p5_histBinDown->GetBinContent(i), WMuNu_SR_below0p5_histBinDown->GetBinContent(i)-WMuNu_SR_below0p5_histBinDown->GetBinError(i)));
    addTemplate("WMuNu_SR_below0p5_WMuNuSignalSBin"+binNumber+"Up", vars, wspace, WMuNu_SR_below0p5_histBinUp);
    addTemplate("WMuNu_SR_below0p5_WMuNuSignalSBin"+binNumber+"Down", vars, wspace, WMuNu_SR_below0p5_histBinDown);
    TH1F* WTauNu_SR_below0p5_histBinUp = (TH1F*)((TH1F*)ZnnGbelow0p5file->Get("histo_WTauNu"))->Clone("WTauNu_SR_below0p5_histBinUp");
    TH1F* WTauNu_SR_below0p5_histBinDown = (TH1F*)((TH1F*)ZnnGbelow0p5file->Get("histo_WTauNu"))->Clone("WTauNu_SR_below0p5_histBinDown");
    WTauNu_SR_below0p5_histBinUp->SetBinContent(i, WTauNu_SR_below0p5_histBinUp->GetBinContent(i)+WTauNu_SR_below0p5_histBinUp->GetBinError(i));
    WTauNu_SR_below0p5_histBinDown->SetBinContent(i, TMath::Max(0.01*WTauNu_SR_below0p5_histBinDown->GetBinContent(i), WTauNu_SR_below0p5_histBinDown->GetBinContent(i)-WTauNu_SR_below0p5_histBinDown->GetBinError(i)));
    addTemplate("WTauNu_SR_below0p5_WTauNuSignalSBin"+binNumber+"Up", vars, wspace, WTauNu_SR_below0p5_histBinUp);
    addTemplate("WTauNu_SR_below0p5_WTauNuSignalSBin"+binNumber+"Down", vars, wspace, WTauNu_SR_below0p5_histBinDown);
    TH1F* WZ_WM_histBinUp = (TH1F*)((TH1F*)WmnGfile->Get("histo_WZ"))->Clone("WZ_WM_histBinUp");
    TH1F* WZ_WM_histBinDown = (TH1F*)((TH1F*)WmnGfile->Get("histo_WZ"))->Clone("WZ_WM_histBinDown");
    WZ_WM_histBinUp->SetBinContent(i, WZ_WM_histBinUp->GetBinContent(i)+WZ_WM_histBinUp->GetBinError(i));
    WZ_WM_histBinDown->SetBinContent(i, TMath::Max(0.01*WZ_WM_histBinDown->GetBinContent(i), WZ_WM_histBinDown->GetBinContent(i)-WZ_WM_histBinDown->GetBinError(i)));
    addTemplate("WZ_WM_WZMonomuSBin"+binNumber+"Up", vars, wspace, WZ_WM_histBinUp);
    addTemplate("WZ_WM_WZMonomuSBin"+binNumber+"Down", vars, wspace, WZ_WM_histBinDown);
    TH1F* WW_WM_histBinUp = (TH1F*)((TH1F*)WmnGfile->Get("histo_WW"))->Clone("WW_WM_histBinUp");
    TH1F* WW_WM_histBinDown = (TH1F*)((TH1F*)WmnGfile->Get("histo_WW"))->Clone("WW_WM_histBinDown");
    WW_WM_histBinUp->SetBinContent(i, WW_WM_histBinUp->GetBinContent(i)+WW_WM_histBinUp->GetBinError(i));
    WW_WM_histBinDown->SetBinContent(i, TMath::Max(0.01*WW_WM_histBinDown->GetBinContent(i), WW_WM_histBinDown->GetBinContent(i)-WW_WM_histBinDown->GetBinError(i)));
    addTemplate("WW_WM_WWMonomuSBin"+binNumber+"Up", vars, wspace, WW_WM_histBinUp);
    addTemplate("WW_WM_WWMonomuSBin"+binNumber+"Down", vars, wspace, WW_WM_histBinDown);
    TH1F* WZ_WE_histBinUp = (TH1F*)((TH1F*)WenGfile->Get("histo_WZ"))->Clone("WZ_WE_histBinUp");
    TH1F* WZ_WE_histBinDown = (TH1F*)((TH1F*)WenGfile->Get("histo_WZ"))->Clone("WZ_WE_histBinDown");
    WZ_WE_histBinUp->SetBinContent(i, WZ_WE_histBinUp->GetBinContent(i)+WZ_WE_histBinUp->GetBinError(i));
    WZ_WE_histBinDown->SetBinContent(i, TMath::Max(0.01*WZ_WE_histBinDown->GetBinContent(i), WZ_WE_histBinDown->GetBinContent(i)-WZ_WE_histBinDown->GetBinError(i)));
    addTemplate("WZ_WE_WZMonoeleSBin"+binNumber+"Up", vars, wspace, WZ_WE_histBinUp);
    addTemplate("WZ_WE_WZMonoeleSBin"+binNumber+"Down", vars, wspace, WZ_WE_histBinDown);
    TH1F* WW_WE_histBinUp = (TH1F*)((TH1F*)WenGfile->Get("histo_WW"))->Clone("WW_WE_histBinUp");
    TH1F* WW_WE_histBinDown = (TH1F*)((TH1F*)WenGfile->Get("histo_WW"))->Clone("WW_WE_histBinDown");
    WW_WE_histBinUp->SetBinContent(i, WW_WE_histBinUp->GetBinContent(i)+WW_WE_histBinUp->GetBinError(i));
    WW_WE_histBinDown->SetBinContent(i, TMath::Max(0.01*WW_WE_histBinDown->GetBinContent(i), WW_WE_histBinDown->GetBinContent(i)-WW_WE_histBinDown->GetBinError(i)));
    addTemplate("WW_WE_WWMonoeleSBin"+binNumber+"Up", vars, wspace, WW_WE_histBinUp);
    addTemplate("WW_WE_WWMonoeleSBin"+binNumber+"Down", vars, wspace, WW_WE_histBinDown);

// Others    
    // TH1F* signal_SR_above0p5_histBinUp = (TH1F*)signal_SR_above0p5_hist->Clone("signal_SR_above0p5_histBinUp");
    // TH1F* signal_SR_above0p5_histBinDown = (TH1F*)signal_SR_above0p5_hist->Clone("signal_SR_above0p5_histBinDown");
    // signal_SR_above0p5_histBinUp->SetBinContent(i, signal_SR_above0p5_hist->GetBinContent(i)+signal_SR_above0p5_hist->GetBinError(i));
    // signal_SR_above0p5_histBinDown->SetBinContent(i, TMath::Max(0.01*signal_SR_above0p5_hist->GetBinContent(i), signal_SR_above0p5_hist->GetBinContent(i)-signal_SR_above0p5_hist->GetBinError(i)));
    // addTemplate("Signal_SR_above0p5_DMSignalSBin"+binNumber+"Up", vars, wspace, signal_SR_above0p5_histBinUp);
    // addTemplate("Signal_SR_above0p5_DMSignalSBin"+binNumber+"Down", vars, wspace, signal_SR_above0p5_histBinDown);
    
    TH1F* QCD_SR_above0p5_histBinUp = (TH1F*)((TH1F*)ZnnGabove0p5file->Get("histo_jetfake"))->Clone("QCD_SR_above0p5_histBinUp");
    TH1F* QCD_SR_above0p5_histBinDown = (TH1F*)((TH1F*)ZnnGabove0p5file->Get("histo_jetfake"))->Clone("QCD_SR_above0p5_histBinDown");
    QCD_SR_above0p5_histBinUp->SetBinContent(i, QCD_SR_above0p5_histBinUp->GetBinContent(i)+QCD_SR_above0p5_histBinUp->GetBinError(i));
    QCD_SR_above0p5_histBinDown->SetBinContent(i, TMath::Max(0.01*QCD_SR_above0p5_histBinDown->GetBinContent(i), QCD_SR_above0p5_histBinDown->GetBinContent(i)-QCD_SR_above0p5_histBinDown->GetBinError(i)));
    addTemplate("QCD_SR_above0p5_QCDHiPhiSignalSBin"+binNumber+"Up", vars, wspace, QCD_SR_above0p5_histBinUp);
    addTemplate("QCD_SR_above0p5_QCDHiPhiSignalSBin"+binNumber+"Down", vars, wspace, QCD_SR_above0p5_histBinDown);
    // TH1F* Elefake_SR_above0p5_histBinUp = (TH1F*)((TH1F*)ZnnGabove0p5file->Get("histo_elefake"))->Clone("Elefake_SR_above0p5_histBinUp");
    // TH1F* Elefake_SR_above0p5_histBinDown = (TH1F*)((TH1F*)ZnnGabove0p5file->Get("histo_elefake"))->Clone("Elefake_SR_above0p5_histBinDown");
    // Elefake_SR_above0p5_histBinUp->SetBinContent(i, Elefake_SR_above0p5_histBinUp->GetBinContent(i)+Elefake_SR_above0p5_histBinUp->GetBinError(i));
    // Elefake_SR_above0p5_histBinDown->SetBinContent(i, TMath::Max(0.01*Elefake_SR_above0p5_histBinDown->GetBinContent(i), Elefake_SR_above0p5_histBinDown->GetBinContent(i)-Elefake_SR_above0p5_histBinDown->GetBinError(i)));
    // addTemplate("Elefake_SR_above0p5_EleHiPhiSignalSBin"+binNumber+"Up", vars, wspace, Elefake_SR_above0p5_histBinUp);
    // addTemplate("Elefake_SR_above0p5_EleHiPhiSignalSBin"+binNumber+"Down", vars, wspace, Elefake_SR_above0p5_histBinDown);
    // TH1F* ZllG_SR_above0p5_histBinUp = (TH1F*)((TH1F*)ZnnGabove0p5file->Get("histo_ZllG_combined"))->Clone("ZllG_SR_above0p5_histBinUp");
    // TH1F* ZllG_SR_above0p5_histBinDown = (TH1F*)((TH1F*)ZnnGabove0p5file->Get("histo_ZllG_combined"))->Clone("ZllG_SR_above0p5_histBinDown");
    // ZllG_SR_above0p5_histBinUp->SetBinContent(i, ZllG_SR_above0p5_histBinUp->GetBinContent(i)+ZllG_SR_above0p5_histBinUp->GetBinError(i));
    // ZllG_SR_above0p5_histBinDown->SetBinContent(i, TMath::Max(0.01*ZllG_SR_above0p5_histBinDown->GetBinContent(i), ZllG_SR_above0p5_histBinDown->GetBinContent(i)-ZllG_SR_above0p5_histBinDown->GetBinError(i)));
    // addTemplate("ZllG_SR_above0p5_ZllGSignalSBin"+binNumber+"Up", vars, wspace, ZllG_SR_above0p5_histBinUp);
    // addTemplate("ZllG_SR_above0p5_ZllGSignalSBin"+binNumber+"Down", vars, wspace, ZllG_SR_above0p5_histBinDown);
    // TH1F* TTG_SR_above0p5_histBinUp = (TH1F*)((TH1F*)ZnnGabove0p5file->Get("histo_TTG"))->Clone("TTG_SR_above0p5_histBinUp");
    // TH1F* TTG_SR_above0p5_histBinDown = (TH1F*)((TH1F*)ZnnGabove0p5file->Get("histo_TTG"))->Clone("TTG_SR_above0p5_histBinDown");
    // TTG_SR_above0p5_histBinUp->SetBinContent(i, TTG_SR_above0p5_histBinUp->GetBinContent(i)+TTG_SR_above0p5_histBinUp->GetBinError(i));
    // TTG_SR_above0p5_histBinDown->SetBinContent(i, TMath::Max(0.01*TTG_SR_above0p5_histBinDown->GetBinContent(i), TTG_SR_above0p5_histBinDown->GetBinContent(i)-TTG_SR_above0p5_histBinDown->GetBinError(i)));
    // addTemplate("TTG_SR_above0p5_TTGSignalSBin"+binNumber+"Up", vars, wspace, TTG_SR_above0p5_histBinUp);
    // addTemplate("TTG_SR_above0p5_TTGSignalSBin"+binNumber+"Down", vars, wspace, TTG_SR_above0p5_histBinDown);
    // TH1F* TG_SR_above0p5_histBinUp = (TH1F*)((TH1F*)ZnnGabove0p5file->Get("histo_TG"))->Clone("TG_SR_above0p5_histBinUp");
    // TH1F* TG_SR_above0p5_histBinDown = (TH1F*)((TH1F*)ZnnGabove0p5file->Get("histo_TG"))->Clone("TG_SR_above0p5_histBinDown");
    // TG_SR_above0p5_histBinUp->SetBinContent(i, TG_SR_above0p5_histBinUp->GetBinContent(i)+TG_SR_above0p5_histBinUp->GetBinError(i));
    // TG_SR_above0p5_histBinDown->SetBinContent(i, TMath::Max(0.01*TG_SR_above0p5_histBinDown->GetBinContent(i), TG_SR_above0p5_histBinDown->GetBinContent(i)-TG_SR_above0p5_histBinDown->GetBinError(i)));
    // addTemplate("TG_SR_above0p5_TGSignalSBin"+binNumber+"Up", vars, wspace, TG_SR_above0p5_histBinUp);
    // addTemplate("TG_SR_above0p5_TGSignalSBin"+binNumber+"Down", vars, wspace, TG_SR_above0p5_histBinDown);
    // TH1F* Diphoton_SR_above0p5_histBinUp = (TH1F*)((TH1F*)ZnnGabove0p5file->Get("histo_diphoton"))->Clone("Diphoton_SR_above0p5_histBinUp");
    // TH1F* Diphoton_SR_above0p5_histBinDown = (TH1F*)((TH1F*)ZnnGabove0p5file->Get("histo_diphoton"))->Clone("Diphoton_SR_above0p5_histBinDown");
    // Diphoton_SR_above0p5_histBinUp->SetBinContent(i, Diphoton_SR_above0p5_histBinUp->GetBinContent(i)+Diphoton_SR_above0p5_histBinUp->GetBinError(i));
    // Diphoton_SR_above0p5_histBinDown->SetBinContent(i, TMath::Max(0.01*Diphoton_SR_above0p5_histBinDown->GetBinContent(i), Diphoton_SR_above0p5_histBinDown->GetBinContent(i)-Diphoton_SR_above0p5_histBinDown->GetBinError(i)));
    // addTemplate("Diphoton_SR_above0p5_DiphotonSignalSBin"+binNumber+"Up", vars, wspace, Diphoton_SR_above0p5_histBinUp);
    // addTemplate("Diphoton_SR_above0p5_DiphotonSignalSBin"+binNumber+"Down", vars, wspace, Diphoton_SR_above0p5_histBinDown);
    TH1F* ZZ_SR_above0p5_histBinUp = (TH1F*)((TH1F*)ZnnGabove0p5file->Get("histo_ZZ"))->Clone("ZZ_SR_above0p5_histBinUp");
    TH1F* ZZ_SR_above0p5_histBinDown = (TH1F*)((TH1F*)ZnnGabove0p5file->Get("histo_ZZ"))->Clone("ZZ_SR_above0p5_histBinDown");
    ZZ_SR_above0p5_histBinUp->SetBinContent(i, ZZ_SR_above0p5_histBinUp->GetBinContent(i)+ZZ_SR_above0p5_histBinUp->GetBinError(i));
    ZZ_SR_above0p5_histBinDown->SetBinContent(i, TMath::Max(0.01*ZZ_SR_above0p5_histBinDown->GetBinContent(i), ZZ_SR_above0p5_histBinDown->GetBinContent(i)-ZZ_SR_above0p5_histBinDown->GetBinError(i)));
    addTemplate("ZZ_SR_above0p5_ZZSignalSBin"+binNumber+"Up", vars, wspace, ZZ_SR_above0p5_histBinUp);
    addTemplate("ZZ_SR_above0p5_ZZSignalSBin"+binNumber+"Down", vars, wspace, ZZ_SR_above0p5_histBinDown);
    TH1F* WW_SR_above0p5_histBinUp = (TH1F*)((TH1F*)ZnnGabove0p5file->Get("histo_WW"))->Clone("WW_SR_above0p5_histBinUp");
    TH1F* WW_SR_above0p5_histBinDown = (TH1F*)((TH1F*)ZnnGabove0p5file->Get("histo_WW"))->Clone("WW_SR_above0p5_histBinDown");
    WW_SR_above0p5_histBinUp->SetBinContent(i, WW_SR_above0p5_histBinUp->GetBinContent(i)+WW_SR_above0p5_histBinUp->GetBinError(i));
    WW_SR_above0p5_histBinDown->SetBinContent(i, TMath::Max(0.01*WW_SR_above0p5_histBinDown->GetBinContent(i), WW_SR_above0p5_histBinDown->GetBinContent(i)-WW_SR_above0p5_histBinDown->GetBinError(i)));
    addTemplate("WW_SR_above0p5_WWSignalSBin"+binNumber+"Up", vars, wspace, WW_SR_above0p5_histBinUp);
    addTemplate("WW_SR_above0p5_WWSignalSBin"+binNumber+"Down", vars, wspace, WW_SR_above0p5_histBinDown);
    
    // TH1F* signal_SR_below0p5_histBinUp = (TH1F*)signal_SR_below0p5_hist->Clone("signal_SR_below0p5_histBinUp");
    // TH1F* signal_SR_below0p5_histBinDown = (TH1F*)signal_SR_below0p5_hist->Clone("signal_SR_below0p5_histBinDown");
    // signal_SR_below0p5_histBinUp->SetBinContent(i, signal_SR_below0p5_hist->GetBinContent(i)+signal_SR_below0p5_hist->GetBinError(i));
    // signal_SR_below0p5_histBinDown->SetBinContent(i, TMath::Max(0.01*signal_SR_below0p5_hist->GetBinContent(i), signal_SR_below0p5_hist->GetBinContent(i)-signal_SR_below0p5_hist->GetBinError(i)));
    // addTemplate("Signal_SR_below0p5_DMSignalSBin"+binNumber+"Up", vars, wspace, signal_SR_below0p5_histBinUp);
    // addTemplate("Signal_SR_below0p5_DMSignalSBin"+binNumber+"Down", vars, wspace, signal_SR_below0p5_histBinDown);
    
    TH1F* QCD_SR_below0p5_histBinUp = (TH1F*)((TH1F*)ZnnGbelow0p5file->Get("histo_jetfake"))->Clone("QCD_SR_below0p5_histBinUp");
    TH1F* QCD_SR_below0p5_histBinDown = (TH1F*)((TH1F*)ZnnGbelow0p5file->Get("histo_jetfake"))->Clone("QCD_SR_below0p5_histBinDown");
    QCD_SR_below0p5_histBinUp->SetBinContent(i, QCD_SR_below0p5_histBinUp->GetBinContent(i)+QCD_SR_below0p5_histBinUp->GetBinError(i));
    QCD_SR_below0p5_histBinDown->SetBinContent(i, TMath::Max(0.01*QCD_SR_below0p5_histBinDown->GetBinContent(i), QCD_SR_below0p5_histBinDown->GetBinContent(i)-QCD_SR_below0p5_histBinDown->GetBinError(i)));
    addTemplate("QCD_SR_below0p5_QCDLoPhiSignalSBin"+binNumber+"Up", vars, wspace, QCD_SR_below0p5_histBinUp);
    addTemplate("QCD_SR_below0p5_QCDLoPhiSignalSBin"+binNumber+"Down", vars, wspace, QCD_SR_below0p5_histBinDown);
    // TH1F* Elefake_SR_below0p5_histBinUp = (TH1F*)((TH1F*)ZnnGbelow0p5file->Get("histo_elefake"))->Clone("Elefake_SR_below0p5_histBinUp");
    // TH1F* Elefake_SR_below0p5_histBinDown = (TH1F*)((TH1F*)ZnnGbelow0p5file->Get("histo_elefake"))->Clone("Elefake_SR_below0p5_histBinDown");
    // Elefake_SR_below0p5_histBinUp->SetBinContent(i, Elefake_SR_below0p5_histBinUp->GetBinContent(i)+Elefake_SR_below0p5_histBinUp->GetBinError(i));
    // Elefake_SR_below0p5_histBinDown->SetBinContent(i, TMath::Max(0.01*Elefake_SR_below0p5_histBinDown->GetBinContent(i), Elefake_SR_below0p5_histBinDown->GetBinContent(i)-Elefake_SR_below0p5_histBinDown->GetBinError(i)));
    // addTemplate("Elefake_SR_below0p5_EleLoPhiSignalSBin"+binNumber+"Up", vars, wspace, Elefake_SR_below0p5_histBinUp);
    // addTemplate("Elefake_SR_below0p5_EleLoPhiSignalSBin"+binNumber+"Down", vars, wspace, Elefake_SR_below0p5_histBinDown);
    // TH1F* ZllG_SR_below0p5_histBinUp = (TH1F*)((TH1F*)ZnnGbelow0p5file->Get("histo_ZllG_combined"))->Clone("ZllG_SR_below0p5_histBinUp");
    // TH1F* ZllG_SR_below0p5_histBinDown = (TH1F*)((TH1F*)ZnnGbelow0p5file->Get("histo_ZllG_combined"))->Clone("ZllG_SR_below0p5_histBinDown");
    // ZllG_SR_below0p5_histBinUp->SetBinContent(i, ZllG_SR_below0p5_histBinUp->GetBinContent(i)+ZllG_SR_below0p5_histBinUp->GetBinError(i));
    // ZllG_SR_below0p5_histBinDown->SetBinContent(i, TMath::Max(0.01*ZllG_SR_below0p5_histBinDown->GetBinContent(i), ZllG_SR_below0p5_histBinDown->GetBinContent(i)-ZllG_SR_below0p5_histBinDown->GetBinError(i)));
    // addTemplate("ZllG_SR_below0p5_ZllGSignalSBin"+binNumber+"Up", vars, wspace, ZllG_SR_below0p5_histBinUp);
    // addTemplate("ZllG_SR_below0p5_ZllGSignalSBin"+binNumber+"Down", vars, wspace, ZllG_SR_below0p5_histBinDown);
    // TH1F* TTG_SR_below0p5_histBinUp = (TH1F*)((TH1F*)ZnnGbelow0p5file->Get("histo_TTG"))->Clone("TTG_SR_below0p5_histBinUp");
    // TH1F* TTG_SR_below0p5_histBinDown = (TH1F*)((TH1F*)ZnnGbelow0p5file->Get("histo_TTG"))->Clone("TTG_SR_below0p5_histBinDown");
    // TTG_SR_below0p5_histBinUp->SetBinContent(i, TTG_SR_below0p5_histBinUp->GetBinContent(i)+TTG_SR_below0p5_histBinUp->GetBinError(i));
    // TTG_SR_below0p5_histBinDown->SetBinContent(i, TMath::Max(0.01*TTG_SR_below0p5_histBinDown->GetBinContent(i), TTG_SR_below0p5_histBinDown->GetBinContent(i)-TTG_SR_below0p5_histBinDown->GetBinError(i)));
    // addTemplate("TTG_SR_below0p5_TTGSignalSBin"+binNumber+"Up", vars, wspace, TTG_SR_below0p5_histBinUp);
    // addTemplate("TTG_SR_below0p5_TTGSignalSBin"+binNumber+"Down", vars, wspace, TTG_SR_below0p5_histBinDown);
    // TH1F* TG_SR_below0p5_histBinUp = (TH1F*)((TH1F*)ZnnGbelow0p5file->Get("histo_TG"))->Clone("TG_SR_below0p5_histBinUp");
    // TH1F* TG_SR_below0p5_histBinDown = (TH1F*)((TH1F*)ZnnGbelow0p5file->Get("histo_TG"))->Clone("TG_SR_below0p5_histBinDown");
    // TG_SR_below0p5_histBinUp->SetBinContent(i, TG_SR_below0p5_histBinUp->GetBinContent(i)+TG_SR_below0p5_histBinUp->GetBinError(i));
    // TG_SR_below0p5_histBinDown->SetBinContent(i, TMath::Max(0.01*TG_SR_below0p5_histBinDown->GetBinContent(i), TG_SR_below0p5_histBinDown->GetBinContent(i)-TG_SR_below0p5_histBinDown->GetBinError(i)));
    // addTemplate("TG_SR_below0p5_TGSignalSBin"+binNumber+"Up", vars, wspace, TG_SR_below0p5_histBinUp);
    // addTemplate("TG_SR_below0p5_TGSignalSBin"+binNumber+"Down", vars, wspace, TG_SR_below0p5_histBinDown);
    // TH1F* Diphoton_SR_below0p5_histBinUp = (TH1F*)((TH1F*)ZnnGbelow0p5file->Get("histo_diphoton"))->Clone("Diphoton_SR_below0p5_histBinUp");
    // TH1F* Diphoton_SR_below0p5_histBinDown = (TH1F*)((TH1F*)ZnnGbelow0p5file->Get("histo_diphoton"))->Clone("Diphoton_SR_below0p5_histBinDown");
    // Diphoton_SR_below0p5_histBinUp->SetBinContent(i, Diphoton_SR_below0p5_histBinUp->GetBinContent(i)+Diphoton_SR_below0p5_histBinUp->GetBinError(i));
    // Diphoton_SR_below0p5_histBinDown->SetBinContent(i, TMath::Max(0.01*Diphoton_SR_below0p5_histBinDown->GetBinContent(i), Diphoton_SR_below0p5_histBinDown->GetBinContent(i)-Diphoton_SR_below0p5_histBinDown->GetBinError(i)));
    // addTemplate("Diphoton_SR_below0p5_DiphotonSignalSBin"+binNumber+"Up", vars, wspace, Diphoton_SR_below0p5_histBinUp);
    // addTemplate("Diphoton_SR_below0p5_DiphotonSignalSBin"+binNumber+"Down", vars, wspace, Diphoton_SR_below0p5_histBinDown);
    TH1F* ZZ_SR_below0p5_histBinUp = (TH1F*)((TH1F*)ZnnGbelow0p5file->Get("histo_ZZ"))->Clone("ZZ_SR_below0p5_histBinUp");
    TH1F* ZZ_SR_below0p5_histBinDown = (TH1F*)((TH1F*)ZnnGbelow0p5file->Get("histo_ZZ"))->Clone("ZZ_SR_below0p5_histBinDown");
    ZZ_SR_below0p5_histBinUp->SetBinContent(i, ZZ_SR_below0p5_histBinUp->GetBinContent(i)+ZZ_SR_below0p5_histBinUp->GetBinError(i));
    ZZ_SR_below0p5_histBinDown->SetBinContent(i, TMath::Max(0.01*ZZ_SR_below0p5_histBinDown->GetBinContent(i), ZZ_SR_below0p5_histBinDown->GetBinContent(i)-ZZ_SR_below0p5_histBinDown->GetBinError(i)));
    addTemplate("ZZ_SR_below0p5_ZZSignalSBin"+binNumber+"Up", vars, wspace, ZZ_SR_below0p5_histBinUp);
    addTemplate("ZZ_SR_below0p5_ZZSignalSBin"+binNumber+"Down", vars, wspace, ZZ_SR_below0p5_histBinDown);
    TH1F* WW_SR_below0p5_histBinUp = (TH1F*)((TH1F*)ZnnGbelow0p5file->Get("histo_WW"))->Clone("WW_SR_below0p5_histBinUp");
    TH1F* WW_SR_below0p5_histBinDown = (TH1F*)((TH1F*)ZnnGbelow0p5file->Get("histo_WW"))->Clone("WW_SR_below0p5_histBinDown");
    WW_SR_below0p5_histBinUp->SetBinContent(i, WW_SR_below0p5_histBinUp->GetBinContent(i)+WW_SR_below0p5_histBinUp->GetBinError(i));
    WW_SR_below0p5_histBinDown->SetBinContent(i, TMath::Max(0.01*WW_SR_below0p5_histBinDown->GetBinContent(i), WW_SR_below0p5_histBinDown->GetBinContent(i)-WW_SR_below0p5_histBinDown->GetBinError(i)));
    addTemplate("WW_SR_below0p5_WWSignalSBin"+binNumber+"Up", vars, wspace, WW_SR_below0p5_histBinUp);
    addTemplate("WW_SR_below0p5_WWSignalSBin"+binNumber+"Down", vars, wspace, WW_SR_below0p5_histBinDown);
       
    TH1F* QCD_ZM_histBinUp = (TH1F*)((TH1F*)ZmmGfile->Get("histo_jetfake"))->Clone("QCD_ZM_histBinUp");
    TH1F* QCD_ZM_histBinDown = (TH1F*)((TH1F*)ZmmGfile->Get("histo_jetfake"))->Clone("QCD_ZM_histBinDown");
    QCD_ZM_histBinUp->SetBinContent(i, QCD_ZM_histBinUp->GetBinContent(i)+QCD_ZM_histBinUp->GetBinError(i));
    QCD_ZM_histBinDown->SetBinContent(i, TMath::Max(0.01*QCD_ZM_histBinDown->GetBinContent(i), QCD_ZM_histBinDown->GetBinContent(i)-QCD_ZM_histBinDown->GetBinError(i)));
    addTemplate("QCD_ZM_QCDDimuSBin"+binNumber+"Up", vars, wspace, QCD_ZM_histBinUp);
    addTemplate("QCD_ZM_QCDDimuSBin"+binNumber+"Down", vars, wspace, QCD_ZM_histBinDown);
    // TH1F* TTG_ZM_histBinUp = (TH1F*)((TH1F*)ZmmGfile->Get("histo_TTG"))->Clone("TTG_ZM_histBinUp");
    // TH1F* TTG_ZM_histBinDown = (TH1F*)((TH1F*)ZmmGfile->Get("histo_TTG"))->Clone("TTG_ZM_histBinDown");
    // TTG_ZM_histBinUp->SetBinContent(i, TTG_ZM_histBinUp->GetBinContent(i)+TTG_ZM_histBinUp->GetBinError(i));
    // TTG_ZM_histBinDown->SetBinContent(i, TMath::Max(0.01*TTG_ZM_histBinDown->GetBinContent(i), TTG_ZM_histBinDown->GetBinContent(i)-TTG_ZM_histBinDown->GetBinError(i)));
    // addTemplate("TTG_ZM_TTGDimuSBin"+binNumber+"Up", vars, wspace, TTG_ZM_histBinUp);
    // addTemplate("TTG_ZM_TTGDimuSBin"+binNumber+"Down", vars, wspace, TTG_ZM_histBinDown);
    // TH1F* WZ_ZM_histBinUp = (TH1F*)((TH1F*)ZmmGfile->Get("histo_WZ"))->Clone("WZ_ZM_histBinUp");
    // TH1F* WZ_ZM_histBinDown = (TH1F*)((TH1F*)ZmmGfile->Get("histo_WZ"))->Clone("WZ_ZM_histBinDown");
    // WZ_ZM_histBinUp->SetBinContent(i, WZ_ZM_histBinUp->GetBinContent(i)+WZ_ZM_histBinUp->GetBinError(i));
    // WZ_ZM_histBinDown->SetBinContent(i, TMath::Max(0.01*WZ_ZM_histBinDown->GetBinContent(i), WZ_ZM_histBinDown->GetBinContent(i)-WZ_ZM_histBinDown->GetBinError(i)));
    // addTemplate("WZ_ZM_WZDimuSBin"+binNumber+"Up", vars, wspace, WZ_ZM_histBinUp);
    // addTemplate("WZ_ZM_WZDimuSBin"+binNumber+"Down", vars, wspace, WZ_ZM_histBinDown);

    TH1F* QCD_ZE_histBinUp = (TH1F*)((TH1F*)ZeeGfile->Get("histo_jetfake"))->Clone("QCD_ZE_histBinUp");
    TH1F* QCD_ZE_histBinDown = (TH1F*)((TH1F*)ZeeGfile->Get("histo_jetfake"))->Clone("QCD_ZE_histBinDown");
    QCD_ZE_histBinUp->SetBinContent(i, QCD_ZE_histBinUp->GetBinContent(i)+QCD_ZE_histBinUp->GetBinError(i));
    QCD_ZE_histBinDown->SetBinContent(i, TMath::Max(0.01*QCD_ZE_histBinDown->GetBinContent(i), QCD_ZE_histBinDown->GetBinContent(i)-QCD_ZE_histBinDown->GetBinError(i)));
    addTemplate("QCD_ZE_QCDDieleSBin"+binNumber+"Up", vars, wspace, QCD_ZE_histBinUp);
    addTemplate("QCD_ZE_QCDDieleSBin"+binNumber+"Down", vars, wspace, QCD_ZE_histBinDown);
    // TH1F* TTG_ZE_histBinUp = (TH1F*)((TH1F*)ZeeGfile->Get("histo_TTG"))->Clone("TTG_ZE_histBinUp");
    // TH1F* TTG_ZE_histBinDown = (TH1F*)((TH1F*)ZeeGfile->Get("histo_TTG"))->Clone("TTG_ZE_histBinDown");
    // TTG_ZE_histBinUp->SetBinContent(i, TTG_ZE_histBinUp->GetBinContent(i)+TTG_ZE_histBinUp->GetBinError(i));
    // TTG_ZE_histBinDown->SetBinContent(i, TMath::Max(0.01*TTG_ZE_histBinDown->GetBinContent(i), TTG_ZE_histBinDown->GetBinContent(i)-TTG_ZE_histBinDown->GetBinError(i)));
    // addTemplate("TTG_ZE_TTGDieleSBin"+binNumber+"Up", vars, wspace, TTG_ZE_histBinUp);
    // addTemplate("TTG_ZE_TTGDieleSBin"+binNumber+"Down", vars, wspace, TTG_ZE_histBinDown);
    // TH1F* WZ_ZE_histBinUp = (TH1F*)((TH1F*)ZeeGfile->Get("histo_WZ"))->Clone("WZ_ZE_histBinUp");
    // TH1F* WZ_ZE_histBinDown = (TH1F*)((TH1F*)ZeeGfile->Get("histo_WZ"))->Clone("WZ_ZE_histBinDown");
    // WZ_ZE_histBinUp->SetBinContent(i, WZ_ZE_histBinUp->GetBinContent(i)+WZ_ZE_histBinUp->GetBinError(i));
    // WZ_ZE_histBinDown->SetBinContent(i, TMath::Max(0.01*WZ_ZE_histBinDown->GetBinContent(i), WZ_ZE_histBinDown->GetBinContent(i)-WZ_ZE_histBinDown->GetBinError(i)));
    // addTemplate("WZ_ZE_WZDieleSBin"+binNumber+"Up", vars, wspace, WZ_ZE_histBinUp);
    // addTemplate("WZ_ZE_WZDieleSBin"+binNumber+"Down", vars, wspace, WZ_ZE_histBinDown);

    TH1F* QCD_WM_histBinUp = (TH1F*)((TH1F*)WmnGfile->Get("histo_jetfake"))->Clone("QCD_WM_histBinUp");
    TH1F* QCD_WM_histBinDown = (TH1F*)((TH1F*)WmnGfile->Get("histo_jetfake"))->Clone("QCD_WM_histBinDown");
    QCD_WM_histBinUp->SetBinContent(i, QCD_WM_histBinUp->GetBinContent(i)+QCD_WM_histBinUp->GetBinError(i));
    QCD_WM_histBinDown->SetBinContent(i, TMath::Max(0.01*QCD_WM_histBinDown->GetBinContent(i), QCD_WM_histBinDown->GetBinContent(i)-QCD_WM_histBinDown->GetBinError(i)));
    addTemplate("QCD_WM_QCDMonomuSBin"+binNumber+"Up", vars, wspace, QCD_WM_histBinUp);
    addTemplate("QCD_WM_QCDMonomuSBin"+binNumber+"Down", vars, wspace, QCD_WM_histBinDown);
    // TH1F* ZllG_WM_histBinUp = (TH1F*)((TH1F*)WmnGfile->Get("histo_ZllG_combined"))->Clone("ZllG_WM_histBinUp");
    // TH1F* ZllG_WM_histBinDown = (TH1F*)((TH1F*)WmnGfile->Get("histo_ZllG_combined"))->Clone("ZllG_WM_histBinDown");
    // ZllG_WM_histBinUp->SetBinContent(i, ZllG_WM_histBinUp->GetBinContent(i)+ZllG_WM_histBinUp->GetBinError(i));
    // ZllG_WM_histBinDown->SetBinContent(i, TMath::Max(0.01*ZllG_WM_histBinDown->GetBinContent(i), ZllG_WM_histBinDown->GetBinContent(i)-ZllG_WM_histBinDown->GetBinError(i)));
    // addTemplate("ZllG_WM_ZllGMonomuSBin"+binNumber+"Up", vars, wspace, ZllG_WM_histBinUp);
    // addTemplate("ZllG_WM_ZllGMonomuSBin"+binNumber+"Down", vars, wspace, ZllG_WM_histBinDown);
    // TH1F* TTG_WM_histBinUp = (TH1F*)((TH1F*)WmnGfile->Get("histo_TTG"))->Clone("TTG_WM_histBinUp");
    // TH1F* TTG_WM_histBinDown = (TH1F*)((TH1F*)WmnGfile->Get("histo_TTG"))->Clone("TTG_WM_histBinDown");
    // TTG_WM_histBinUp->SetBinContent(i, TTG_WM_histBinUp->GetBinContent(i)+TTG_WM_histBinUp->GetBinError(i));
    // TTG_WM_histBinDown->SetBinContent(i, TMath::Max(0.01*TTG_WM_histBinDown->GetBinContent(i), TTG_WM_histBinDown->GetBinContent(i)-TTG_WM_histBinDown->GetBinError(i)));
    // addTemplate("TTG_WM_TTGMonomuSBin"+binNumber+"Up", vars, wspace, TTG_WM_histBinUp);
    // addTemplate("TTG_WM_TTGMonomuSBin"+binNumber+"Down", vars, wspace, TTG_WM_histBinDown);
    // TH1F* TG_WM_histBinUp = (TH1F*)((TH1F*)WmnGfile->Get("histo_TG"))->Clone("TG_WM_histBinUp");
    // TH1F* TG_WM_histBinDown = (TH1F*)((TH1F*)WmnGfile->Get("histo_TG"))->Clone("TG_WM_histBinDown");
    // TG_WM_histBinUp->SetBinContent(i, TG_WM_histBinUp->GetBinContent(i)+TG_WM_histBinUp->GetBinError(i));
    // TG_WM_histBinDown->SetBinContent(i, TMath::Max(0.01*TG_WM_histBinDown->GetBinContent(i), TG_WM_histBinDown->GetBinContent(i)-TG_WM_histBinDown->GetBinError(i)));
    // addTemplate("TG_WM_TGMonomuSBin"+binNumber+"Up", vars, wspace, TG_WM_histBinUp);
    // addTemplate("TG_WM_TGMonomuSBin"+binNumber+"Down", vars, wspace, TG_WM_histBinDown);
    // TH1F* Diphoton_WM_histBinUp = (TH1F*)((TH1F*)WmnGfile->Get("histo_diphoton"))->Clone("Diphoton_WM_histBinUp");
    // TH1F* Diphoton_WM_histBinDown = (TH1F*)((TH1F*)WmnGfile->Get("histo_diphoton"))->Clone("Diphoton_WM_histBinDown");
    // Diphoton_WM_histBinUp->SetBinContent(i, Diphoton_WM_histBinUp->GetBinContent(i)+Diphoton_WM_histBinUp->GetBinError(i));
    // Diphoton_WM_histBinDown->SetBinContent(i, TMath::Max(0.01*Diphoton_WM_histBinDown->GetBinContent(i), Diphoton_WM_histBinDown->GetBinContent(i)-Diphoton_WM_histBinDown->GetBinError(i)));
    // addTemplate("Diphoton_WM_DiphotonMonomuSBin"+binNumber+"Up", vars, wspace, Diphoton_WM_histBinUp);
    // addTemplate("Diphoton_WM_DiphotonMonomuSBin"+binNumber+"Down", vars, wspace, Diphoton_WM_histBinDown);

    TH1F* QCD_WE_histBinUp = (TH1F*)((TH1F*)WenGfile->Get("histo_jetfake"))->Clone("QCD_WE_histBinUp");
    TH1F* QCD_WE_histBinDown = (TH1F*)((TH1F*)WenGfile->Get("histo_jetfake"))->Clone("QCD_WE_histBinDown");
    QCD_WE_histBinUp->SetBinContent(i, QCD_WE_histBinUp->GetBinContent(i)+QCD_WE_histBinUp->GetBinError(i));
    QCD_WE_histBinDown->SetBinContent(i, TMath::Max(0.01*QCD_WE_histBinDown->GetBinContent(i), QCD_WE_histBinDown->GetBinContent(i)-QCD_WE_histBinDown->GetBinError(i)));
    addTemplate("QCD_WE_QCDMonoeleSBin"+binNumber+"Up", vars, wspace, QCD_WE_histBinUp);
    addTemplate("QCD_WE_QCDMonoeleSBin"+binNumber+"Down", vars, wspace, QCD_WE_histBinDown);
    // TH1F* Elefake_WE_histBinUp = (TH1F*)((TH1F*)WenGfile->Get("histo_elefake"))->Clone("Elefake_WE_histBinUp");
    // TH1F* Elefake_WE_histBinDown = (TH1F*)((TH1F*)WenGfile->Get("histo_elefake"))->Clone("Elefake_WE_histBinDown");
    // Elefake_WE_histBinUp->SetBinContent(i, Elefake_WE_histBinUp->GetBinContent(i)+Elefake_WE_histBinUp->GetBinError(i));
    // Elefake_WE_histBinDown->SetBinContent(i, TMath::Max(0.01*Elefake_WE_histBinDown->GetBinContent(i), Elefake_WE_histBinDown->GetBinContent(i)-Elefake_WE_histBinDown->GetBinError(i)));
    // addTemplate("Elefake_WE_EleMonoeleSBin"+binNumber+"Up", vars, wspace, Elefake_WE_histBinUp);
    // addTemplate("Elefake_WE_EleMonoeleSBin"+binNumber+"Down", vars, wspace, Elefake_WE_histBinDown);
    // TH1F* ZllG_WE_histBinUp = (TH1F*)((TH1F*)WenGfile->Get("histo_ZllG_combined"))->Clone("ZllG_WE_histBinUp");
    // TH1F* ZllG_WE_histBinDown = (TH1F*)((TH1F*)WenGfile->Get("histo_ZllG_combined"))->Clone("ZllG_WE_histBinDown");
    // ZllG_WE_histBinUp->SetBinContent(i, ZllG_WE_histBinUp->GetBinContent(i)+ZllG_WE_histBinUp->GetBinError(i));
    // ZllG_WE_histBinDown->SetBinContent(i, TMath::Max(0.01*ZllG_WE_histBinDown->GetBinContent(i), ZllG_WE_histBinDown->GetBinContent(i)-ZllG_WE_histBinDown->GetBinError(i)));
    // addTemplate("ZllG_WE_ZllGMonoeleSBin"+binNumber+"Up", vars, wspace, ZllG_WE_histBinUp);
    // addTemplate("ZllG_WE_ZllGMonoeleSBin"+binNumber+"Down", vars, wspace, ZllG_WE_histBinDown);
    // TH1F* TTG_WE_histBinUp = (TH1F*)((TH1F*)WenGfile->Get("histo_TTG"))->Clone("TTG_WE_histBinUp");
    // TH1F* TTG_WE_histBinDown = (TH1F*)((TH1F*)WenGfile->Get("histo_TTG"))->Clone("TTG_WE_histBinDown");
    // TTG_WE_histBinUp->SetBinContent(i, TTG_WE_histBinUp->GetBinContent(i)+TTG_WE_histBinUp->GetBinError(i));
    // TTG_WE_histBinDown->SetBinContent(i, TMath::Max(0.01*TTG_WE_histBinDown->GetBinContent(i), TTG_WE_histBinDown->GetBinContent(i)-TTG_WE_histBinDown->GetBinError(i)));
    // addTemplate("TTG_WE_TTGMonoeleSBin"+binNumber+"Up", vars, wspace, TTG_WE_histBinUp);
    // addTemplate("TTG_WE_TTGMonoeleSBin"+binNumber+"Down", vars, wspace, TTG_WE_histBinDown);
    // TH1F* TG_WE_histBinUp = (TH1F*)((TH1F*)WenGfile->Get("histo_TG"))->Clone("TG_WE_histBinUp");
    // TH1F* TG_WE_histBinDown = (TH1F*)((TH1F*)WenGfile->Get("histo_TG"))->Clone("TG_WE_histBinDown");
    // TG_WE_histBinUp->SetBinContent(i, TG_WE_histBinUp->GetBinContent(i)+TG_WE_histBinUp->GetBinError(i));
    // TG_WE_histBinDown->SetBinContent(i, TMath::Max(0.01*TG_WE_histBinDown->GetBinContent(i), TG_WE_histBinDown->GetBinContent(i)-TG_WE_histBinDown->GetBinError(i)));
    // addTemplate("TG_WE_TGMonoeleSBin"+binNumber+"Up", vars, wspace, TG_WE_histBinUp);
    // addTemplate("TG_WE_TGMonoeleSBin"+binNumber+"Down", vars, wspace, TG_WE_histBinDown);
    // TH1F* Diphoton_WE_histBinUp = (TH1F*)((TH1F*)WenGfile->Get("histo_diphoton"))->Clone("Diphoton_WE_histBinUp");
    // TH1F* Diphoton_WE_histBinDown = (TH1F*)((TH1F*)WenGfile->Get("histo_diphoton"))->Clone("Diphoton_WE_histBinDown");
    // Diphoton_WE_histBinUp->SetBinContent(i, Diphoton_WE_histBinUp->GetBinContent(i)+Diphoton_WE_histBinUp->GetBinError(i));
    // Diphoton_WE_histBinDown->SetBinContent(i, TMath::Max(0.01*Diphoton_WE_histBinDown->GetBinContent(i), Diphoton_WE_histBinDown->GetBinContent(i)-Diphoton_WE_histBinDown->GetBinError(i)));
    // addTemplate("Diphoton_WE_DiphotonMonoeleSBin"+binNumber+"Up", vars, wspace, Diphoton_WE_histBinUp);
    // addTemplate("Diphoton_WE_DiphotonMonoeleSBin"+binNumber+"Down", vars, wspace, Diphoton_WE_histBinDown);
  }

  // ---------------------------- Write out the workspace -----------------------------------------------------------------//
  outfile->cd();
  wspace.Write();
  outfile->Close();
  
  transfer_factors_file->Close();
  signalabove0p5file->Close();
  signalbelow0p5file->Close();
  ZnnGabove0p5file->Close();
  ZnnGbelow0p5file->Close();
  WenGfile->Close();
  WmnGfile->Close();
  ZeeGfile->Close();
  ZmmGfile->Close();
}
 

void createWorkspaces_Pt_HaloFit_DMinterpolation(){
  bool connectWZ = true;
  
  std::vector<string> sample_names_DM_NLOfromInterpolation;
  sample_names_DM_NLOfromInterpolation.clear();
  std::vector<float> signal_multipliers_DM_NLOfromInterpolation;
  signal_multipliers_DM_NLOfromInterpolation.clear();
  
  do_createWorkspace("V_MxTarget400_MvTarget850", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget325_MvTarget600", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget125_MvTarget1400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget200_MvTarget1200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget1_MvTarget450", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget400_MvTarget1175", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget25_MvTarget450", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget5_MvTarget625", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget10_MvTarget1250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget300_MvTarget800", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget100_MvTarget1150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget10_MvTarget1125", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget5_MvTarget1150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget225_MvTarget150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget650_MvTarget1100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget100_MvTarget325", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget275_MvTarget150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget375_MvTarget400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget900_MvTarget1725", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget475_MvTarget400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget50_MvTarget1400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget500_MvTarget1075", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget300_MvTarget500", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget900_MvTarget1525", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget400_MvTarget200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget50_MvTarget950", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget425_MvTarget600", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget475_MvTarget600", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget25_MvTarget350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget50_MvTarget575", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget175_MvTarget650", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget1_MvTarget350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget25_MvTarget300", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget1_MvTarget400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget350_MvTarget1100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget250_MvTarget1200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget50_MvTarget1200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget275_MvTarget450", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget25_MvTarget1200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget125_MvTarget350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget5_MvTarget150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget350_MvTarget350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget250_MvTarget300", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget375_MvTarget800", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget525_MvTarget400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget400_MvTarget1400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget600_MvTarget500", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget350_MvTarget850", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget100_MvTarget1400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget275_MvTarget10", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget475_MvTarget1100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget225_MvTarget650", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget650_MvTarget950", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget225_MvTarget950", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget50_MvTarget225", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget200_MvTarget1100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget225_MvTarget550", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget475_MvTarget1150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget575_MvTarget400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget275_MvTarget550", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget200_MvTarget1175", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget1_MvTarget625", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget1_MvTarget1350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget600_MvTarget1000", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget50_MvTarget1150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget200_MvTarget350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget25_MvTarget175", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget125_MvTarget200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget275_MvTarget850", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget5_MvTarget30", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget450_MvTarget950", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget275_MvTarget950", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget575_MvTarget1025", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget325_MvTarget250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget125_MvTarget10", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget325_MvTarget700", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget400_MvTarget800", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget25_MvTarget500", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget50_MvTarget850", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget125_MvTarget600", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget700_MvTarget1925", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget50_MvTarget150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget250_MvTarget1400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget300_MvTarget1175", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget150_MvTarget1125", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget25_MvTarget60", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget175_MvTarget1250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget275_MvTarget1175", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget5_MvTarget350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget100_MvTarget1000", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget625_MvTarget1100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget525_MvTarget800", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget100_MvTarget725", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget200_MvTarget1050", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget1000_MvTarget1800", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget300_MvTarget750", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget10_MvTarget650", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget650_MvTarget1250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget25_MvTarget800", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget75_MvTarget100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget650_MvTarget1400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget1_MvTarget800", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget10_MvTarget675", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget625_MvTarget1125", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget500_MvTarget1150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget25_MvTarget1175", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget75_MvTarget1125", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget175_MvTarget575", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget250_MvTarget800", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget5_MvTarget1100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget1_MvTarget850", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget400_MvTarget1350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget400_MvTarget650", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget625_MvTarget1000", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget575_MvTarget1200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget500_MvTarget450", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget50_MvTarget500", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget1_MvTarget1125", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget125_MvTarget1100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget525_MvTarget700", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget75_MvTarget1000", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget250_MvTarget1150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget325_MvTarget550", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget200_MvTarget500", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget300_MvTarget700", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget50_MvTarget900", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget250_MvTarget200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget350_MvTarget750", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget10_MvTarget225", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget700_MvTarget1800", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget450_MvTarget1050", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget50_MvTarget100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget175_MvTarget950", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget300_MvTarget100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget1_MvTarget1400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget250_MvTarget500", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget500_MvTarget500", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget400_MvTarget700", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget375_MvTarget100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget150_MvTarget1050", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget5_MvTarget850", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget325_MvTarget650", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget50_MvTarget1250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget5_MvTarget175", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget200_MvTarget1400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget150_MvTarget1250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget10_MvTarget400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget275_MvTarget1100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget450_MvTarget1100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget350_MvTarget1250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget200_MvTarget800", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget225_MvTarget1000", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget625_MvTarget1175", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget250_MvTarget1350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget325_MvTarget1175", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget500_MvTarget1200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget425_MvTarget550", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget75_MvTarget450", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget250_MvTarget1125", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget475_MvTarget550", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget300_MvTarget600", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget125_MvTarget1125", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget5_MvTarget325", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget225_MvTarget700", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget625_MvTarget450", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget650_MvTarget650", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget450_MvTarget1025", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget425_MvTarget950", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget50_MvTarget425", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget200_MvTarget850", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget300_MvTarget300", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget325_MvTarget1050", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget1_MvTarget250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget10_MvTarget375", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget350_MvTarget200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget75_MvTarget600", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget575_MvTarget1250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget5_MvTarget100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget375_MvTarget850", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget450_MvTarget550", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget550_MvTarget1125", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget225_MvTarget250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget275_MvTarget1250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget325_MvTarget1350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget400_MvTarget1250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget175_MvTarget1100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget25_MvTarget950", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget1_MvTarget90", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget75_MvTarget300", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget1_MvTarget1200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget25_MvTarget900", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget75_MvTarget50", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget100_MvTarget150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget700_MvTarget1125", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget175_MvTarget1050", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget75_MvTarget250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget100_MvTarget100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget375_MvTarget1250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget10_MvTarget1300", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget1250_MvTarget2000", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget475_MvTarget1050", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget150_MvTarget950", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget300_MvTarget525", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget200_MvTarget400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget5_MvTarget1400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget50_MvTarget650", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget900_MvTarget1800", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget150_MvTarget750", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget125_MvTarget100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget475_MvTarget650", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget150_MvTarget1300", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget350_MvTarget1175", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget10_MvTarget950", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget1_MvTarget700", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget1_MvTarget40", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget5_MvTarget1300", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget550_MvTarget1075", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget650_MvTarget1125", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget575_MvTarget500", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget350_MvTarget450", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget325_MvTarget400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget325_MvTarget950", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget250_MvTarget350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget1_MvTarget125", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget200_MvTarget950", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget10_MvTarget500", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget275_MvTarget50", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget200_MvTarget900", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget75_MvTarget525", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget650_MvTarget1175", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget300_MvTarget1150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget100_MvTarget250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget200_MvTarget550", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget5_MvTarget1000", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget600_MvTarget1125", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget1_MvTarget375", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget525_MvTarget650", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget5_MvTarget725", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget625_MvTarget650", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget375_MvTarget1200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget300_MvTarget250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget75_MvTarget200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget650_MvTarget750", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget525_MvTarget450", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget50_MvTarget400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget600_MvTarget600", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget125_MvTarget1050", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget625_MvTarget1075", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget525_MvTarget1050", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget10_MvTarget350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget75_MvTarget850", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget475_MvTarget950", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget200_MvTarget725", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget625_MvTarget1150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget650_MvTarget700", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget5_MvTarget950", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget175_MvTarget850", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget150_MvTarget900", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget350_MvTarget600", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget100_MvTarget800", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget700_MvTarget1725", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget425_MvTarget1250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget275_MvTarget1400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget200_MvTarget250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget50_MvTarget90", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget450_MvTarget1000", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget500_MvTarget1350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget50_MvTarget475", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget375_MvTarget1150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget625_MvTarget700", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget175_MvTarget1200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget350_MvTarget1125", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget525_MvTarget600", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget50_MvTarget350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget75_MvTarget700", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget25_MvTarget550", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget275_MvTarget400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget250_MvTarget750", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget225_MvTarget900", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget175_MvTarget300", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget100_MvTarget850", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget325_MvTarget150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget525_MvTarget1125", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget200_MvTarget150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget550_MvTarget1350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget400_MvTarget600", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget275_MvTarget100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget175_MvTarget525", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget275_MvTarget800", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget200_MvTarget100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget300_MvTarget1100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget100_MvTarget1250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget225_MvTarget400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget525_MvTarget1250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget425_MvTarget1150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget225_MvTarget1200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget575_MvTarget1175", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget450_MvTarget600", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget275_MvTarget1200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget500_MvTarget800", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget325_MvTarget500", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget25_MvTarget70", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget25_MvTarget250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget100_MvTarget750", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget1_MvTarget70", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget175_MvTarget900", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget175_MvTarget450", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget300_MvTarget725", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget625_MvTarget900", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget400_MvTarget50", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget250_MvTarget700", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget5_MvTarget80", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget325_MvTarget1000", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget400_MvTarget400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget1_MvTarget950", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget300_MvTarget1250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget5_MvTarget300", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget550_MvTarget800", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget275_MvTarget1300", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget400_MvTarget500", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget600_MvTarget1100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget10_MvTarget30", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget10_MvTarget1050", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget25_MvTarget150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget625_MvTarget600", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget150_MvTarget1175", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget475_MvTarget500", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget50_MvTarget1175", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget1_MvTarget1175", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget25_MvTarget100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget500_MvTarget950", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget50_MvTarget175", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget100_MvTarget950", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget225_MvTarget450", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget1_MvTarget80", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget375_MvTarget10", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget575_MvTarget950", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget150_MvTarget1000", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget325_MvTarget850", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget600_MvTarget1350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget5_MvTarget125", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget400_MvTarget450", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget225_MvTarget1250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget5_MvTarget10", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget100_MvTarget1300", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget550_MvTarget450", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget400_MvTarget10", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget50_MvTarget675", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget175_MvTarget1150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget325_MvTarget1200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget650_MvTarget1025", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget10_MvTarget200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget350_MvTarget1000", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget75_MvTarget1050", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget25_MvTarget1250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget175_MvTarget200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget600_MvTarget1075", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget275_MvTarget350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget400_MvTarget350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget125_MvTarget300", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget625_MvTarget500", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget525_MvTarget900", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget900_MvTarget1400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget225_MvTarget50", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget1000_MvTarget1600", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget5_MvTarget250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget400_MvTarget1000", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget425_MvTarget1350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget25_MvTarget50", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget150_MvTarget1400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget25_MvTarget575", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget175_MvTarget1175", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget650_MvTarget1350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget525_MvTarget1100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget375_MvTarget950", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget425_MvTarget1200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget475_MvTarget1200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget550_MvTarget1100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget50_MvTarget10", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget650_MvTarget900", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget5_MvTarget90", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget650_MvTarget1150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget150_MvTarget700", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget500_MvTarget750", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget100_MvTarget1050", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget450_MvTarget700", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget550_MvTarget600", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget10_MvTarget175", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget250_MvTarget600", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget225_MvTarget500", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget10_MvTarget625", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget375_MvTarget50", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget1_MvTarget275", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget75_MvTarget1350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget5_MvTarget800", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget225_MvTarget1300", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget200_MvTarget750", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget200_MvTarget1250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget625_MvTarget1350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget150_MvTarget1150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget250_MvTarget1100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget10_MvTarget525", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget250_MvTarget1000", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget600_MvTarget1050", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget325_MvTarget10", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget400_MvTarget550", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget1_MvTarget550", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget275_MvTarget1125", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget600_MvTarget450", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget450_MvTarget750", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget525_MvTarget1400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget25_MvTarget650", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget500_MvTarget1025", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget1_MvTarget1250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget250_MvTarget1050", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget75_MvTarget575", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget125_MvTarget1250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget625_MvTarget400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget475_MvTarget1350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget225_MvTarget850", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget10_MvTarget1175", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget1_MvTarget30", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget500_MvTarget1175", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget5_MvTarget1250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget350_MvTarget650", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget1_MvTarget425", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget50_MvTarget375", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget125_MvTarget750", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget5_MvTarget40", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget75_MvTarget1100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget300_MvTarget350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget650_MvTarget800", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget125_MvTarget150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget500_MvTarget900", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget275_MvTarget300", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget100_MvTarget500", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget400_MvTarget1050", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget375_MvTarget450", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget625_MvTarget550", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget425_MvTarget1025", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget125_MvTarget700", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget575_MvTarget900", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget700_MvTarget1325", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget75_MvTarget550", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget200_MvTarget650", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget25_MvTarget90", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget250_MvTarget10", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget25_MvTarget400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget275_MvTarget650", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget450_MvTarget1250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget1_MvTarget900", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget575_MvTarget1125", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget325_MvTarget350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget1_MvTarget150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget500_MvTarget1000", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget600_MvTarget800", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget375_MvTarget550", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget100_MvTarget1200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget550_MvTarget1050", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget600_MvTarget750", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget10_MvTarget125", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget550_MvTarget700", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget350_MvTarget550", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget25_MvTarget700", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget575_MvTarget600", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget250_MvTarget450", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget25_MvTarget80", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget400_MvTarget950", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget350_MvTarget950", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget1_MvTarget225", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget350_MvTarget1400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget10_MvTarget60", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget900_MvTarget1325", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget300_MvTarget550", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget550_MvTarget1150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget100_MvTarget10", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget50_MvTarget1100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget375_MvTarget350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget525_MvTarget1000", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget125_MvTarget1000", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget250_MvTarget1300", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget50_MvTarget725", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget125_MvTarget525", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget175_MvTarget550", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget475_MvTarget1125", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget125_MvTarget550", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget300_MvTarget900", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget10_MvTarget900", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget25_MvTarget125", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget50_MvTarget1125", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget1_MvTarget1300", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget10_MvTarget80", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget175_MvTarget750", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget25_MvTarget1300", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget175_MvTarget700", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget325_MvTarget750", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget525_MvTarget750", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget275_MvTarget200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget50_MvTarget1000", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget575_MvTarget450", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget5_MvTarget50", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget450_MvTarget1125", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget550_MvTarget550", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget375_MvTarget1350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget275_MvTarget900", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget75_MvTarget800", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget550_MvTarget1400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget525_MvTarget1175", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget625_MvTarget1250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget550_MvTarget950", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget200_MvTarget1350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget5_MvTarget1200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget400_MvTarget1300", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget75_MvTarget1250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget100_MvTarget400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget375_MvTarget150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget600_MvTarget700", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget225_MvTarget100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget700_MvTarget1000", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget350_MvTarget900", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget250_MvTarget900", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget125_MvTarget650", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget300_MvTarget450", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget375_MvTarget200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget325_MvTarget1150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget350_MvTarget150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget10_MvTarget1350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget25_MvTarget10", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget500_MvTarget1125", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget575_MvTarget1050", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget450_MvTarget1150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget425_MvTarget1400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget500_MvTarget700", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget5_MvTarget275", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget75_MvTarget1400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget600_MvTarget1150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget75_MvTarget1300", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget325_MvTarget100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget100_MvTarget525", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget1_MvTarget725", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget1_MvTarget525", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget375_MvTarget1300", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget300_MvTarget400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget200_MvTarget600", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget400_MvTarget150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget25_MvTarget1350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget225_MvTarget1400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget500_MvTarget600", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget175_MvTarget100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget100_MvTarget650", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget600_MvTarget550", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget500_MvTarget1100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget225_MvTarget1175", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget275_MvTarget1000", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget225_MvTarget1350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget600_MvTarget1250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget25_MvTarget1050", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget5_MvTarget1050", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget200_MvTarget10", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget5_MvTarget60", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget1_MvTarget60", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget225_MvTarget300", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget1_MvTarget1150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget125_MvTarget400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget5_MvTarget200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget5_MvTarget500", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget150_MvTarget850", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget1_MvTarget1100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget25_MvTarget675", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget350_MvTarget100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget325_MvTarget50", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget600_MvTarget1175", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget300_MvTarget1300", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget575_MvTarget1075", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget5_MvTarget550", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget75_MvTarget1150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget150_MvTarget800", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget375_MvTarget1400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget175_MvTarget250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget300_MvTarget650", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget350_MvTarget500", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget575_MvTarget1150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget250_MvTarget1175", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget25_MvTarget600", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget10_MvTarget450", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget525_MvTarget1150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget350_MvTarget1200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget650_MvTarget1000", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget75_MvTarget150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget75_MvTarget500", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget225_MvTarget200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget550_MvTarget750", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget25_MvTarget225", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget50_MvTarget1300", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget10_MvTarget1400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget75_MvTarget950", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget125_MvTarget50", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget575_MvTarget700", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget175_MvTarget600", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget525_MvTarget1200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget250_MvTarget50", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget10_MvTarget850", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget125_MvTarget1350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget125_MvTarget575", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget25_MvTarget325", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget50_MvTarget250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget525_MvTarget1350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget1_MvTarget325", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget175_MvTarget1000", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget10_MvTarget600", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget125_MvTarget1300", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget100_MvTarget1350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget600_MvTarget650", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget225_MvTarget1125", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget400_MvTarget1200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget350_MvTarget400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget350_MvTarget1350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget550_MvTarget1250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget10_MvTarget275", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget300_MvTarget50", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget325_MvTarget800", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget200_MvTarget700", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget5_MvTarget525", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget25_MvTarget475", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget100_MvTarget300", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget225_MvTarget600", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget5_MvTarget375", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget275_MvTarget600", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget25_MvTarget375", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget325_MvTarget200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget1_MvTarget575", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget100_MvTarget50", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget50_MvTarget125", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget575_MvTarget1100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget125_MvTarget1150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget5_MvTarget70", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget100_MvTarget575", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget275_MvTarget250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget10_MvTarget325", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget225_MvTarget1150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget500_MvTarget1250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget50_MvTarget325", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget1_MvTarget475", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget600_MvTarget1400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget650_MvTarget1075", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget375_MvTarget500", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget600_MvTarget950", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget5_MvTarget400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget250_MvTarget150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget200_MvTarget325", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget250_MvTarget400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget225_MvTarget1050", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget525_MvTarget1025", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget10_MvTarget1200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget350_MvTarget300", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget100_MvTarget900", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget475_MvTarget1025", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget625_MvTarget1200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget200_MvTarget50", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget5_MvTarget600", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget350_MvTarget800", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget600_MvTarget400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget175_MvTarget1300", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget5_MvTarget225", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget75_MvTarget1200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget10_MvTarget800", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget175_MvTarget10", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget10_MvTarget40", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget125_MvTarget1175", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget375_MvTarget1000", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget425_MvTarget1175", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget150_MvTarget1100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget100_MvTarget550", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget325_MvTarget1125", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget550_MvTarget650", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget25_MvTarget275", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget5_MvTarget475", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget10_MvTarget1100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget450_MvTarget900", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget200_MvTarget1300", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget125_MvTarget950", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget10_MvTarget475", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget225_MvTarget750", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget625_MvTarget1050", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget425_MvTarget750", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget125_MvTarget1200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget475_MvTarget750", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget125_MvTarget850", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget25_MvTarget625", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget5_MvTarget450", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget475_MvTarget1400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget275_MvTarget700", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget5_MvTarget20", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget325_MvTarget1250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget300_MvTarget1200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget5_MvTarget650", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget10_MvTarget700", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget75_MvTarget1175", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget400_MvTarget1100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget250_MvTarget650", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget50_MvTarget600", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget500_MvTarget650", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget575_MvTarget800", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget50_MvTarget525", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget10_MvTarget150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget25_MvTarget850", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget600_MvTarget1025", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget425_MvTarget1125", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget75_MvTarget750", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget250_MvTarget550", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget400_MvTarget1150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget275_MvTarget1150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget25_MvTarget1000", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget1000_MvTarget1925", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget25_MvTarget1100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget275_MvTarget1050", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget100_MvTarget200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget1_MvTarget675", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget450_MvTarget1075", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget325_MvTarget300", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget625_MvTarget750", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget5_MvTarget675", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget1_MvTarget1050", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget50_MvTarget1050", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget375_MvTarget700", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget100_MvTarget700", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget1000_MvTarget1725", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget550_MvTarget400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget550_MvTarget1175", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget550_MvTarget1200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget300_MvTarget200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget450_MvTarget650", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget450_MvTarget1350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget425_MvTarget1000", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget100_MvTarget1175", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget575_MvTarget1400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget375_MvTarget1125", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget275_MvTarget1350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget400_MvTarget1125", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget600_MvTarget900", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget450_MvTarget1400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget10_MvTarget90", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget575_MvTarget1350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget500_MvTarget1400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget525_MvTarget950", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget650_MvTarget1200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget1_MvTarget650", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget700_MvTarget1200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget10_MvTarget550", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget125_MvTarget900", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget475_MvTarget900", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget425_MvTarget900", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget5_MvTarget1175", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget10_MvTarget1000", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget375_MvTarget1050", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget200_MvTarget1000", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget75_MvTarget400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget175_MvTarget1400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget325_MvTarget1100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget250_MvTarget1250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget5_MvTarget1350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget375_MvTarget1100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget75_MvTarget650", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget375_MvTarget250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget450_MvTarget800", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget150_MvTarget1350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget650_MvTarget1050", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget550_MvTarget1025", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget50_MvTarget1350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget425_MvTarget1075", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget300_MvTarget150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget525_MvTarget1075", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget400_MvTarget300", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget1_MvTarget1000", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget325_MvTarget900", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget10_MvTarget1150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget275_MvTarget500", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget50_MvTarget550", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget350_MvTarget250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget475_MvTarget700", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget350_MvTarget50", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget5_MvTarget575", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget575_MvTarget650", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget75_MvTarget350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget50_MvTarget700", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget25_MvTarget1400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget175_MvTarget800", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget250_MvTarget100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget200_MvTarget1125", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget50_MvTarget450", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget375_MvTarget600", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget600_MvTarget1200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget575_MvTarget550", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget50_MvTarget80", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget10_MvTarget300", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget10_MvTarget250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget1_MvTarget600", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget475_MvTarget1000", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget325_MvTarget1300", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget550_MvTarget900", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget700_MvTarget1400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget5_MvTarget900", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget525_MvTarget500", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget250_MvTarget950", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget50_MvTarget625", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget100_MvTarget600", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget175_MvTarget500", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget125_MvTarget250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget425_MvTarget800", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget550_MvTarget500", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget10_MvTarget70", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget200_MvTarget200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget375_MvTarget300", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget300_MvTarget1125", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget250_MvTarget250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget300_MvTarget10", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget200_MvTarget1150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget475_MvTarget1250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget225_MvTarget1100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget25_MvTarget525", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget475_MvTarget1175", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget425_MvTarget1100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget375_MvTarget1175", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget175_MvTarget1350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget350_MvTarget10", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget175_MvTarget1125", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget550_MvTarget1000", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget375_MvTarget650", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget125_MvTarget725", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget350_MvTarget1300", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget225_MvTarget350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget500_MvTarget1050", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget10_MvTarget725", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget25_MvTarget1125", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget225_MvTarget800", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget25_MvTarget40", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget625_MvTarget1025", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget400_MvTarget900", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget300_MvTarget950", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget100_MvTarget450", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget5_MvTarget425", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget1_MvTarget175", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget25_MvTarget425", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget525_MvTarget550", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget175_MvTarget50", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget700_MvTarget1600", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget50_MvTarget275", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget300_MvTarget1050", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget450_MvTarget1200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget450_MvTarget1175", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget500_MvTarget400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget300_MvTarget1000", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget5_MvTarget700", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget625_MvTarget800", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget100_MvTarget1125", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget25_MvTarget200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget350_MvTarget1150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget75_MvTarget900", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget175_MvTarget350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget25_MvTarget1150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget325_MvTarget450", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget200_MvTarget525", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget275_MvTarget750", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget125_MvTarget500", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget200_MvTarget450", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget150_MvTarget725", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget350_MvTarget1050", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget250_MvTarget850", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget10_MvTarget425", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget425_MvTarget1050", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget575_MvTarget750", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget425_MvTarget700", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget500_MvTarget550", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget575_MvTarget1000", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget225_MvTarget10", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget125_MvTarget800", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget150_MvTarget1200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget300_MvTarget1350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget100_MvTarget1100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget425_MvTarget650", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget1000_MvTarget1525", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget10_MvTarget575", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget10_MvTarget20", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget325_MvTarget1400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget625_MvTarget1400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget300_MvTarget1400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget1250_MvTarget1925", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget300_MvTarget850", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget5_MvTarget1125", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget625_MvTarget950", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget175_MvTarget150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget175_MvTarget400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget475_MvTarget1075", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget475_MvTarget450", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget125_MvTarget450", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget50_MvTarget800", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget375_MvTarget900", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget200_MvTarget300", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget25_MvTarget725", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget400_MvTarget100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget475_MvTarget800", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("V_MxTarget75_MvTarget10", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget400_MvTarget850", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget325_MvTarget600", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget125_MvTarget1400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget200_MvTarget1200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget1_MvTarget450", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget400_MvTarget1175", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget25_MvTarget450", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget5_MvTarget625", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget10_MvTarget1250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget300_MvTarget800", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget100_MvTarget1150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget10_MvTarget1125", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget5_MvTarget1150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget225_MvTarget150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget650_MvTarget1100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget100_MvTarget325", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget275_MvTarget150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget375_MvTarget400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget900_MvTarget1725", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget475_MvTarget400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget350_MvTarget300", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget500_MvTarget1075", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget300_MvTarget500", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget900_MvTarget1525", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget400_MvTarget250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget400_MvTarget200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget50_MvTarget950", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget425_MvTarget600", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget475_MvTarget600", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget25_MvTarget350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget50_MvTarget575", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget175_MvTarget650", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget1_MvTarget350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget25_MvTarget300", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget1_MvTarget400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget350_MvTarget1100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget250_MvTarget1200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget50_MvTarget1200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget275_MvTarget450", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget25_MvTarget1200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget125_MvTarget350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget5_MvTarget150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget350_MvTarget350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget250_MvTarget300", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget375_MvTarget800", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget525_MvTarget400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget400_MvTarget1400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget600_MvTarget500", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget350_MvTarget850", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget100_MvTarget1400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget275_MvTarget10", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget475_MvTarget1100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget225_MvTarget650", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget650_MvTarget950", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget225_MvTarget950", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget50_MvTarget225", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget200_MvTarget1100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget225_MvTarget550", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget475_MvTarget1150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget575_MvTarget400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget275_MvTarget550", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget200_MvTarget1175", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget1_MvTarget625", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget1_MvTarget1350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget600_MvTarget1000", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget50_MvTarget1150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget200_MvTarget350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget25_MvTarget175", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget125_MvTarget200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget175_MvTarget600", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget5_MvTarget30", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget450_MvTarget950", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget275_MvTarget950", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget575_MvTarget1025", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget325_MvTarget250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget125_MvTarget10", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget325_MvTarget700", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget400_MvTarget800", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget25_MvTarget500", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget50_MvTarget850", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget125_MvTarget600", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget700_MvTarget1925", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget50_MvTarget150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget250_MvTarget1400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget300_MvTarget1175", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget150_MvTarget1125", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget25_MvTarget60", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget175_MvTarget1250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget275_MvTarget1175", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget5_MvTarget350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget100_MvTarget1000", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget625_MvTarget1100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget525_MvTarget800", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget100_MvTarget725", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget200_MvTarget1050", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget1000_MvTarget1800", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget300_MvTarget750", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget650_MvTarget1250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget25_MvTarget800", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget75_MvTarget100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget650_MvTarget1400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget1_MvTarget800", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget10_MvTarget675", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget625_MvTarget1125", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget500_MvTarget1150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget25_MvTarget1175", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget75_MvTarget1125", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget175_MvTarget575", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget250_MvTarget800", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget5_MvTarget1100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget1_MvTarget850", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget400_MvTarget1350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget400_MvTarget650", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget625_MvTarget1000", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget575_MvTarget1200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget500_MvTarget450", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget50_MvTarget500", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget1_MvTarget1125", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget125_MvTarget1100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget525_MvTarget700", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget75_MvTarget1000", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget250_MvTarget1150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget325_MvTarget550", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget200_MvTarget500", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget300_MvTarget700", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget50_MvTarget900", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget250_MvTarget200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget350_MvTarget750", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget10_MvTarget225", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget700_MvTarget1800", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget450_MvTarget1050", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget50_MvTarget100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget175_MvTarget950", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget300_MvTarget100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget1_MvTarget1400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget250_MvTarget500", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget500_MvTarget500", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget400_MvTarget700", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget375_MvTarget100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget150_MvTarget1050", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget5_MvTarget850", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget325_MvTarget650", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget50_MvTarget1250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget5_MvTarget175", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget200_MvTarget1400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget150_MvTarget1250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget10_MvTarget400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget275_MvTarget1100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget450_MvTarget1100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget350_MvTarget1250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget200_MvTarget800", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget225_MvTarget1000", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget625_MvTarget1175", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget250_MvTarget1350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget325_MvTarget1175", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget500_MvTarget1200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget425_MvTarget550", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget75_MvTarget450", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget250_MvTarget1125", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget475_MvTarget550", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget300_MvTarget600", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget125_MvTarget1125", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget5_MvTarget325", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget225_MvTarget700", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget625_MvTarget450", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget650_MvTarget650", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget450_MvTarget1025", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget425_MvTarget950", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget150_MvTarget525", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget200_MvTarget850", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget300_MvTarget300", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget325_MvTarget1050", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget1_MvTarget250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget10_MvTarget375", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget350_MvTarget200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget75_MvTarget600", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget575_MvTarget1250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget5_MvTarget100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget375_MvTarget850", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget450_MvTarget550", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget550_MvTarget1125", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget225_MvTarget250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget275_MvTarget1250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget325_MvTarget1350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget400_MvTarget1250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget175_MvTarget1100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget25_MvTarget950", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget1_MvTarget90", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget75_MvTarget300", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget1_MvTarget1200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget25_MvTarget900", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget75_MvTarget50", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget100_MvTarget150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget700_MvTarget1125", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget175_MvTarget1050", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget75_MvTarget250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget100_MvTarget100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget375_MvTarget1250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget10_MvTarget1300", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget1250_MvTarget2000", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget475_MvTarget1050", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget150_MvTarget950", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget300_MvTarget525", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget200_MvTarget400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget5_MvTarget1400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget50_MvTarget650", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget900_MvTarget1800", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget150_MvTarget750", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget125_MvTarget100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget475_MvTarget650", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget150_MvTarget1300", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget375_MvTarget10", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget10_MvTarget950", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget1_MvTarget700", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget1_MvTarget40", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget5_MvTarget1300", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget550_MvTarget1075", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget350_MvTarget50", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget650_MvTarget1125", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget450_MvTarget700", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget350_MvTarget450", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget325_MvTarget400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget325_MvTarget950", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget250_MvTarget350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget1_MvTarget125", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget200_MvTarget950", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget10_MvTarget500", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget275_MvTarget50", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget200_MvTarget900", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget75_MvTarget525", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget650_MvTarget1175", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget300_MvTarget1150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget100_MvTarget250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget200_MvTarget550", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget5_MvTarget1000", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget600_MvTarget1125", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget1_MvTarget375", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget525_MvTarget650", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget5_MvTarget725", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget625_MvTarget650", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget375_MvTarget1200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget300_MvTarget250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget75_MvTarget200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget650_MvTarget750", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget525_MvTarget450", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget50_MvTarget400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget600_MvTarget600", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget125_MvTarget1050", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget150_MvTarget400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget625_MvTarget1075", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget525_MvTarget1050", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget10_MvTarget350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget75_MvTarget850", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget475_MvTarget950", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget200_MvTarget725", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget625_MvTarget1150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget650_MvTarget700", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget5_MvTarget950", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget150_MvTarget150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget175_MvTarget850", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget150_MvTarget900", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget350_MvTarget600", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget100_MvTarget800", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget700_MvTarget1725", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget425_MvTarget1250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget275_MvTarget1400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget200_MvTarget250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget50_MvTarget90", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget450_MvTarget1000", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget500_MvTarget1350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget50_MvTarget475", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget375_MvTarget1150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget625_MvTarget700", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget175_MvTarget1200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget350_MvTarget1125", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget525_MvTarget600", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget50_MvTarget350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget75_MvTarget700", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget25_MvTarget550", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget275_MvTarget400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget250_MvTarget750", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget225_MvTarget900", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget175_MvTarget300", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget100_MvTarget850", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget325_MvTarget150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget525_MvTarget1125", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget200_MvTarget150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget550_MvTarget1350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget400_MvTarget600", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget275_MvTarget100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget175_MvTarget525", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget275_MvTarget800", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget200_MvTarget100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget300_MvTarget1100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget100_MvTarget1250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget225_MvTarget400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget525_MvTarget1250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget425_MvTarget1150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget225_MvTarget1200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget575_MvTarget1175", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget450_MvTarget600", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget275_MvTarget1200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget500_MvTarget800", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget325_MvTarget500", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget25_MvTarget70", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget25_MvTarget250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget100_MvTarget750", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget1_MvTarget70", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget175_MvTarget900", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget175_MvTarget450", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget300_MvTarget725", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget625_MvTarget900", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget400_MvTarget50", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget250_MvTarget700", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget5_MvTarget80", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget325_MvTarget1000", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget400_MvTarget400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget150_MvTarget100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget300_MvTarget1250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget5_MvTarget300", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget550_MvTarget800", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget275_MvTarget1300", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget400_MvTarget500", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget600_MvTarget1100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget10_MvTarget30", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget10_MvTarget1050", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget25_MvTarget150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget625_MvTarget600", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget150_MvTarget1175", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget475_MvTarget500", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget50_MvTarget1175", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget1_MvTarget1175", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget25_MvTarget100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget500_MvTarget950", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget50_MvTarget175", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget100_MvTarget950", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget225_MvTarget450", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget1_MvTarget80", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget575_MvTarget950", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget150_MvTarget1000", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget325_MvTarget850", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget650_MvTarget400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget600_MvTarget1350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget5_MvTarget125", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget400_MvTarget450", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget225_MvTarget1250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget5_MvTarget10", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget100_MvTarget1300", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget550_MvTarget450", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget400_MvTarget10", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget50_MvTarget675", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget175_MvTarget1150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget325_MvTarget1200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget650_MvTarget1025", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget10_MvTarget200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget350_MvTarget1000", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget75_MvTarget1050", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget25_MvTarget1250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget175_MvTarget200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget600_MvTarget1075", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget275_MvTarget350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget400_MvTarget350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget125_MvTarget300", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget625_MvTarget500", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget525_MvTarget900", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget900_MvTarget1400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget225_MvTarget50", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget1000_MvTarget1600", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget5_MvTarget250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget400_MvTarget1000", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget425_MvTarget1350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget25_MvTarget50", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget150_MvTarget1400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget25_MvTarget575", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget175_MvTarget1175", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget650_MvTarget1350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget525_MvTarget1100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget375_MvTarget950", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget425_MvTarget1200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget475_MvTarget1200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget550_MvTarget1100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget50_MvTarget10", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget650_MvTarget900", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget5_MvTarget90", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget650_MvTarget1150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget150_MvTarget700", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget500_MvTarget750", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget100_MvTarget1050", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget575_MvTarget500", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget550_MvTarget600", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget10_MvTarget175", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget700_MvTarget1000", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget225_MvTarget500", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget10_MvTarget625", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget375_MvTarget50", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget1_MvTarget275", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget75_MvTarget1350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget5_MvTarget800", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget225_MvTarget1300", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget200_MvTarget750", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget200_MvTarget1250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget625_MvTarget1350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget150_MvTarget1150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget250_MvTarget1100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget10_MvTarget525", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget650_MvTarget600", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget250_MvTarget1000", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget600_MvTarget1050", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget325_MvTarget10", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget400_MvTarget550", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget1_MvTarget550", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget275_MvTarget1125", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget600_MvTarget450", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget450_MvTarget750", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget525_MvTarget1400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget25_MvTarget650", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget500_MvTarget1025", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget1_MvTarget1250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget250_MvTarget1050", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget75_MvTarget575", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget125_MvTarget1250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget625_MvTarget400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget475_MvTarget1350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget225_MvTarget850", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget10_MvTarget1175", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget1_MvTarget30", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget500_MvTarget1175", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget150_MvTarget650", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget5_MvTarget1250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget350_MvTarget650", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget1_MvTarget425", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget125_MvTarget750", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget125_MvTarget650", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget5_MvTarget40", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget75_MvTarget1100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget300_MvTarget350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget650_MvTarget800", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget125_MvTarget150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget500_MvTarget900", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget275_MvTarget300", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget100_MvTarget500", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget400_MvTarget1050", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget375_MvTarget450", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget625_MvTarget550", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget425_MvTarget1025", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget125_MvTarget700", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget575_MvTarget900", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget700_MvTarget1325", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget75_MvTarget550", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget200_MvTarget650", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget25_MvTarget90", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget250_MvTarget10", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget25_MvTarget400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget275_MvTarget650", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget450_MvTarget1250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget1_MvTarget950", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget1_MvTarget900", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget575_MvTarget1125", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget325_MvTarget350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget1_MvTarget150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget500_MvTarget1000", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget600_MvTarget800", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget375_MvTarget550", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget100_MvTarget1200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget550_MvTarget1050", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget600_MvTarget750", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget10_MvTarget125", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget550_MvTarget700", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget350_MvTarget550", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget25_MvTarget700", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget575_MvTarget600", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget250_MvTarget450", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget25_MvTarget80", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget400_MvTarget950", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget350_MvTarget950", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget150_MvTarget300", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget1_MvTarget225", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget350_MvTarget1400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget900_MvTarget1325", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget300_MvTarget550", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget550_MvTarget1150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget100_MvTarget10", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget50_MvTarget1100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget375_MvTarget350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget525_MvTarget1000", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget125_MvTarget1000", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget250_MvTarget1300", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget50_MvTarget725", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget125_MvTarget525", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget175_MvTarget550", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget475_MvTarget1125", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget125_MvTarget550", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget300_MvTarget900", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget10_MvTarget900", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget25_MvTarget125", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget50_MvTarget1125", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget1_MvTarget1300", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget10_MvTarget80", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget175_MvTarget750", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget25_MvTarget1300", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget175_MvTarget700", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget325_MvTarget750", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget525_MvTarget750", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget275_MvTarget200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget50_MvTarget1000", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget575_MvTarget450", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget5_MvTarget50", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget450_MvTarget1125", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget550_MvTarget550", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget375_MvTarget1350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget275_MvTarget900", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget75_MvTarget800", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget550_MvTarget1400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget525_MvTarget1175", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget625_MvTarget1250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget550_MvTarget950", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget200_MvTarget1350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget5_MvTarget1200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget400_MvTarget1300", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget75_MvTarget1250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget650_MvTarget550", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget100_MvTarget400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget375_MvTarget150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget600_MvTarget700", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget225_MvTarget100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget250_MvTarget600", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget250_MvTarget900", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget150_MvTarget600", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget300_MvTarget450", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget375_MvTarget200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget325_MvTarget1150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget350_MvTarget150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget10_MvTarget1350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget25_MvTarget10", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget500_MvTarget1125", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget1_MvTarget175", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget450_MvTarget1150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget425_MvTarget1400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget500_MvTarget700", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget5_MvTarget275", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget75_MvTarget1400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget600_MvTarget1150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget75_MvTarget1300", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget325_MvTarget100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget100_MvTarget525", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget1_MvTarget725", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget1_MvTarget525", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget375_MvTarget1300", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget300_MvTarget400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget200_MvTarget600", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget400_MvTarget150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget25_MvTarget1350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget225_MvTarget1400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget500_MvTarget600", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget175_MvTarget100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget100_MvTarget650", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget600_MvTarget550", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget500_MvTarget1100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget225_MvTarget1175", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget275_MvTarget1000", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget225_MvTarget1350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget650_MvTarget450", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget600_MvTarget1250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget25_MvTarget1050", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget5_MvTarget1050", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget200_MvTarget10", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget5_MvTarget60", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget1_MvTarget60", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget225_MvTarget300", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget1_MvTarget1150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget125_MvTarget400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget5_MvTarget200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget5_MvTarget500", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget150_MvTarget850", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget1_MvTarget1100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget25_MvTarget675", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget350_MvTarget100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget325_MvTarget50", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget600_MvTarget1175", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget300_MvTarget1300", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget575_MvTarget1075", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget5_MvTarget550", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget75_MvTarget1150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget150_MvTarget800", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget375_MvTarget1400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget175_MvTarget250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget300_MvTarget650", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget350_MvTarget500", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget575_MvTarget1150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget250_MvTarget1175", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget25_MvTarget600", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget10_MvTarget450", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget525_MvTarget1150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget350_MvTarget1200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget650_MvTarget1000", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget75_MvTarget150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget75_MvTarget500", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget225_MvTarget200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget550_MvTarget750", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget25_MvTarget225", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget50_MvTarget1300", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget10_MvTarget1400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget75_MvTarget950", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget125_MvTarget50", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget575_MvTarget700", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget275_MvTarget850", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget525_MvTarget1200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget250_MvTarget50", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget10_MvTarget850", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget125_MvTarget1350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget125_MvTarget575", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget25_MvTarget325", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget50_MvTarget250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget525_MvTarget1350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget1_MvTarget325", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget175_MvTarget1000", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget10_MvTarget600", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget125_MvTarget1300", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget100_MvTarget1350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget600_MvTarget650", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget225_MvTarget1125", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget400_MvTarget1200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget350_MvTarget400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget550_MvTarget1250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget300_MvTarget50", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget325_MvTarget800", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget200_MvTarget700", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget150_MvTarget325", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget5_MvTarget525", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget25_MvTarget475", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget100_MvTarget300", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget225_MvTarget600", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget5_MvTarget375", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget275_MvTarget600", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget25_MvTarget375", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget325_MvTarget200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget1_MvTarget575", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget100_MvTarget50", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget50_MvTarget125", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget575_MvTarget1100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget125_MvTarget1150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget5_MvTarget70", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget100_MvTarget575", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget275_MvTarget250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget10_MvTarget325", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget225_MvTarget1150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget150_MvTarget550", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget125_MvTarget900", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget500_MvTarget1250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget50_MvTarget325", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget1_MvTarget475", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget600_MvTarget1400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget650_MvTarget1075", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget375_MvTarget500", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget600_MvTarget950", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget5_MvTarget400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget250_MvTarget150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget200_MvTarget325", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget250_MvTarget400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget225_MvTarget1050", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget525_MvTarget1025", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget10_MvTarget1200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget100_MvTarget900", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget475_MvTarget1025", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget625_MvTarget1200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget200_MvTarget50", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget5_MvTarget600", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget350_MvTarget800", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget600_MvTarget400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget175_MvTarget1300", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget5_MvTarget225", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget75_MvTarget1200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget10_MvTarget800", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget175_MvTarget10", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget150_MvTarget250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget10_MvTarget40", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget125_MvTarget1175", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget375_MvTarget1000", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget425_MvTarget1175", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget150_MvTarget1100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget100_MvTarget550", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget325_MvTarget1125", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget550_MvTarget650", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget25_MvTarget275", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget5_MvTarget475", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget10_MvTarget1100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget450_MvTarget900", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget200_MvTarget1300", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget125_MvTarget950", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget10_MvTarget475", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget225_MvTarget750", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget625_MvTarget1050", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget425_MvTarget750", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget125_MvTarget1200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget475_MvTarget750", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget125_MvTarget850", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget25_MvTarget625", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget5_MvTarget450", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget475_MvTarget1400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget275_MvTarget700", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget5_MvTarget20", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget325_MvTarget1250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget300_MvTarget1200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget5_MvTarget650", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget10_MvTarget700", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget10_MvTarget275", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget400_MvTarget1100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget250_MvTarget650", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget50_MvTarget600", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget500_MvTarget650", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget575_MvTarget800", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget50_MvTarget525", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget10_MvTarget150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget25_MvTarget850", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget600_MvTarget1025", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget425_MvTarget1125", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget75_MvTarget750", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget250_MvTarget550", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget400_MvTarget1150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget275_MvTarget1150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget25_MvTarget1000", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget1000_MvTarget1925", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget25_MvTarget1100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget275_MvTarget1050", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget100_MvTarget200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget1_MvTarget675", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget450_MvTarget1075", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget325_MvTarget300", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget625_MvTarget750", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget5_MvTarget675", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget1_MvTarget1050", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget50_MvTarget1050", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget375_MvTarget700", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget100_MvTarget700", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget1000_MvTarget1725", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget550_MvTarget400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget550_MvTarget1175", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget550_MvTarget1200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget300_MvTarget200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget450_MvTarget650", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget450_MvTarget1350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget425_MvTarget1000", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget100_MvTarget1175", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget575_MvTarget1400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget375_MvTarget1125", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget275_MvTarget1350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget400_MvTarget1125", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget600_MvTarget900", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget450_MvTarget1400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget10_MvTarget90", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget575_MvTarget1350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget500_MvTarget1400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget525_MvTarget950", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget650_MvTarget1200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget1_MvTarget650", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget700_MvTarget1200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget10_MvTarget550", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget150_MvTarget575", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget475_MvTarget900", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget425_MvTarget900", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget5_MvTarget1175", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget10_MvTarget1000", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget375_MvTarget1050", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget200_MvTarget1000", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget75_MvTarget400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget175_MvTarget1400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget325_MvTarget1100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget250_MvTarget1250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget5_MvTarget1350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget75_MvTarget650", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget375_MvTarget250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget450_MvTarget800", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget150_MvTarget1350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget650_MvTarget1050", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget550_MvTarget1025", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget50_MvTarget1350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget425_MvTarget1075", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget300_MvTarget150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget525_MvTarget1075", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget400_MvTarget300", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget1_MvTarget1000", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget325_MvTarget900", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget10_MvTarget1150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget275_MvTarget500", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget50_MvTarget550", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget350_MvTarget250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget475_MvTarget700", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget650_MvTarget500", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget5_MvTarget575", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget575_MvTarget650", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget75_MvTarget350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget50_MvTarget700", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget25_MvTarget1400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget175_MvTarget800", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget250_MvTarget100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget200_MvTarget1125", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget50_MvTarget450", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget375_MvTarget600", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget600_MvTarget1200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget575_MvTarget550", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget50_MvTarget80", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget10_MvTarget300", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget10_MvTarget250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget1_MvTarget600", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget475_MvTarget1000", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget325_MvTarget1300", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget550_MvTarget900", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget700_MvTarget1400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget5_MvTarget900", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget525_MvTarget500", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget250_MvTarget950", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget100_MvTarget600", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget175_MvTarget500", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget125_MvTarget250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget425_MvTarget800", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget550_MvTarget500", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget10_MvTarget70", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget200_MvTarget200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget375_MvTarget300", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget300_MvTarget1125", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget250_MvTarget250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget300_MvTarget10", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget200_MvTarget1150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget475_MvTarget1250", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget225_MvTarget1100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget25_MvTarget525", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget475_MvTarget1175", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget425_MvTarget1100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget375_MvTarget1175", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget175_MvTarget1350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget350_MvTarget10", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget175_MvTarget1125", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget550_MvTarget1000", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget375_MvTarget650", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget125_MvTarget725", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget350_MvTarget1300", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget225_MvTarget350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget500_MvTarget1050", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget10_MvTarget725", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget25_MvTarget1125", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget225_MvTarget800", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget25_MvTarget40", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget625_MvTarget1025", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget400_MvTarget900", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget300_MvTarget950", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget100_MvTarget450", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget150_MvTarget450", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget5_MvTarget425", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget575_MvTarget1050", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget25_MvTarget425", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget150_MvTarget50", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget525_MvTarget550", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget175_MvTarget50", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget700_MvTarget1600", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget50_MvTarget275", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget300_MvTarget1050", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget450_MvTarget1200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget450_MvTarget1175", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget500_MvTarget400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget300_MvTarget1000", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget5_MvTarget700", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget625_MvTarget800", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget100_MvTarget1125", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget25_MvTarget200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget350_MvTarget1150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget75_MvTarget900", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget175_MvTarget350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget25_MvTarget1150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget325_MvTarget450", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget275_MvTarget750", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget125_MvTarget500", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget200_MvTarget450", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget150_MvTarget725", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget350_MvTarget1050", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget250_MvTarget850", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget10_MvTarget425", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget425_MvTarget1050", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget575_MvTarget750", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget425_MvTarget700", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget500_MvTarget550", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget575_MvTarget1000", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget225_MvTarget10", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget125_MvTarget800", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget150_MvTarget1200", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget300_MvTarget1350", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget100_MvTarget1100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget425_MvTarget650", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget1000_MvTarget1525", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget10_MvTarget575", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget10_MvTarget20", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget325_MvTarget1400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget625_MvTarget1400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget300_MvTarget1400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget1250_MvTarget1925", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget300_MvTarget850", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget5_MvTarget1125", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget625_MvTarget950", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget175_MvTarget150", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget175_MvTarget400", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget475_MvTarget1075", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget475_MvTarget450", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget125_MvTarget450", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget50_MvTarget800", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget375_MvTarget900", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget200_MvTarget300", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget25_MvTarget725", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget400_MvTarget100", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget475_MvTarget800", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  do_createWorkspace("AV_MxTarget75_MvTarget10", connectWZ, sample_names_DM_NLOfromInterpolation, signal_multipliers_DM_NLOfromInterpolation, "DM_NLOfromInterpolation");
  
  // Printout for createDatacards.py
  
  cout<<"\n\n\n"<<endl;
  
  for(int i = 0; i < sample_names_DM_NLOfromInterpolation.size(); i++){
    if(i == 0)
      cout<<"DM_NLOfromInterpolation_samples_and_multipliers = {"<<"'"<<sample_names_DM_NLOfromInterpolation[i]<<"': '"<<signal_multipliers_DM_NLOfromInterpolation[i]<<"',"<<endl;
    else if(i == sample_names_DM_NLOfromInterpolation.size()-1)
      cout<<"'"<<sample_names_DM_NLOfromInterpolation[i]<<"': '"<<signal_multipliers_DM_NLOfromInterpolation[i]<<"'}"<<endl;
    else
      cout<<"'"<<sample_names_DM_NLOfromInterpolation[i]<<"': '"<<signal_multipliers_DM_NLOfromInterpolation[i]<<"',"<<endl;
  }
}