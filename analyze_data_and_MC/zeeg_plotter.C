#include <fstream>
#include <vector>
#include <iomanip>
#include "TFile.h"
#include "TH2.h"
#include "TH2F.h"
#include "TGraph2D.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TString.h"
#include "THStack.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TColor.h"
#include "boost/format.hpp"

void plot(string histname_string, Double_t leg_xoffset, Double_t leg_yoffset, TString xaxis_title, TString plotname)
{
  TString histname = TString(histname_string+"_4");
  TString histname_JESUp = TString(histname_string+"_9");
  TString histname_JESDown = TString(histname_string+"_14");
  TString histname_PESUp = TString(histname_string+"_19");
  TString histname_PESDown = TString(histname_string+"_24");
  TString histname_straightUp = TString(histname_string+"_26");
  TString histname_straightDown = TString(histname_string+"_27");
  TString histname_twistedUp = TString(histname_string+"_29");
  TString histname_twistedDown = TString(histname_string+"_30");
  TString histname_gammaUp = TString(histname_string+"_31");
  TString histname_gammaDown = TString(histname_string+"_32");
  TString histname_uncorrected = TString(histname_string+"_25");
  TString histname_qcdscale = TString(histname_string+"_28");
  TString histname_qcd_sidebandUp = TString(histname_string+"_5");
  TString histname_qcd_sidebandDown = TString(histname_string+"_6");
  TString histname_qcd_METUp = TString(histname_string+"_7");
  TString histname_qcd_METDown = TString(histname_string+"_8");
  TString histname_qcd_binningUp = TString(histname_string+"_9");
  TString histname_qcd_binningDown = TString(histname_string+"_10");
  TString histname_qcd_sieieLeft = TString(histname_string+"_11");
  TString histname_qcd_sieieRight = TString(histname_string+"_12");
  TString histname_qcd_templateUp = TString(histname_string+"_13");
  TString histname_qcd_templateDown = TString(histname_string+"_14");
  TString histname_qcd_unweighted = TString(histname_string+"_15");

  std::vector<TH1F*> histo_vector;
  histo_vector.clear();
  
  Float_t int_lumi = 35900.0;
  Float_t scale_factor = 0.984*1.002; // Flat pixel seed veto SF times flat pho ID SF
  
  double photon_scale_factor_unc = sqrt(pow(0.009/0.984,2)+pow(0.007/1.002,2)); // combining pix seed and pho ID components
  double electron_scale_factor_unc = 0.04; // uncertainty on the square of the electron scale factor
  
  double total_background = 0.0;
  double background_unc_sumsquares = 0.0;
  
  TCanvas *c = new TCanvas("c", "canvas",700,640);
  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  c->SetFillColor(0);
  c->SetBorderMode(0);
  c->SetBorderSize(2);
  float t_m = 0.08; //top margin
  float b_m = 0.4; //botton margin
  float l_m = 0.09; //left margin
  float r_m = 0.05; //right margin
  c->SetTopMargin(t_m);
  c->SetBottomMargin(b_m);
  c->SetLeftMargin(l_m);
  c->SetRightMargin(r_m);
  c->SetFrameFillStyle(0);
  c->SetFrameBorderMode(0);
  c->SetFrameFillStyle(0);
  c->SetFrameBorderMode(0);
  c->cd();
  
  TFile *f_data = new TFile("ZeeG_data_all.root");
  TH1F* histo_data = (TH1F*)((TH1F*)f_data->Get(histname))->Clone("data_obs");
  const int nBins = histo_data->GetXaxis()->GetNbins();
  histo_data->SetBinContent(nBins, histo_data->GetBinContent(nBins)+histo_data->GetBinContent(nBins+1));
  histo_data->ClearUnderflowAndOverflow();
  Float_t int_data = histo_data->Integral();
  histo_data->SetLineWidth(3);
  histo_data->SetLineColor(kWhite);
  histo_data->SetMarkerStyle(kFullSquare);
  histo_data->SetMarkerColor(kWhite);
  histo_vector.push_back(histo_data);
  
  //Now that nBins has been specified, initialize binned MET shift repositories to the appropriate length
  std::vector<double> jesup_shift;
  jesup_shift.clear();
  std::vector<double> jesdown_shift;
  jesdown_shift.clear();
  std::vector<double> pesup_shift;
  pesup_shift.clear();
  std::vector<double> pesdown_shift;
  pesdown_shift.clear();
  std::vector<double> renup_shift;
  renup_shift.clear();
  std::vector<double> rendown_shift;
  rendown_shift.clear();
  std::vector<double> facup_shift;
  facup_shift.clear();
  std::vector<double> facdown_shift;
  facdown_shift.clear();
  std::vector<double> pdfup_shift;
  pdfup_shift.clear();
  std::vector<double> pdfdown_shift;
  pdfdown_shift.clear();
  std::vector<double> straightup_shift_ZllG;
  straightup_shift_ZllG.clear();
  std::vector<double> straightdown_shift_ZllG;
  straightdown_shift_ZllG.clear();
  std::vector<double> twistedup_shift_ZllG;
  twistedup_shift_ZllG.clear();
  std::vector<double> twisteddown_shift_ZllG;
  twisteddown_shift_ZllG.clear();
  std::vector<double> gammaup_shift_ZllG;
  gammaup_shift_ZllG.clear();
  std::vector<double> gammadown_shift_ZllG;
  gammadown_shift_ZllG.clear();
  std::vector<double> qcdscale_shift_ZllG;
  qcdscale_shift_ZllG.clear();
  std::vector<double> syst_shiftUp_jetfake;
  syst_shiftUp_jetfake.clear();
  std::vector<double> syst_shiftDown_jetfake;
  syst_shiftDown_jetfake.clear();
  for(int i = 1; i <= nBins; i++){
    jesup_shift.push_back(0);
    jesdown_shift.push_back(0);
    pesup_shift.push_back(0);
    pesdown_shift.push_back(0);
    renup_shift.push_back(0);
    rendown_shift.push_back(0);
    facup_shift.push_back(0);
    facdown_shift.push_back(0);
    pdfup_shift.push_back(0);
    pdfdown_shift.push_back(0);
    straightup_shift_ZllG.push_back(0);
    straightdown_shift_ZllG.push_back(0);
    twistedup_shift_ZllG.push_back(0);
    twisteddown_shift_ZllG.push_back(0);
    gammaup_shift_ZllG.push_back(0);
    gammadown_shift_ZllG.push_back(0);
    qcdscale_shift_ZllG.push_back(0);
    syst_shiftUp_jetfake.push_back(0);
    syst_shiftDown_jetfake.push_back(0);
  }
  
  TFile *f_jetfake = new TFile("ZeeG_qcd_all.root");
  TH1F* histo_jetfake = (TH1F*)((TH1F*)f_jetfake->Get(histname))->Clone("histo_jetfake");
  histo_jetfake->SetBinContent(nBins, histo_jetfake->GetBinContent(nBins)+histo_jetfake->GetBinContent(nBins+1));
  histo_jetfake->ClearUnderflowAndOverflow();
  TH1F* histo_jetfake_sidebandUp = (TH1F*)((TH1F*)f_jetfake->Get(histname_qcd_sidebandUp))->Clone("histo_jetfake_sidebandUp");
  histo_jetfake_sidebandUp->SetBinContent(nBins, histo_jetfake_sidebandUp->GetBinContent(nBins)+histo_jetfake_sidebandUp->GetBinContent(nBins+1));
  histo_jetfake_sidebandUp->ClearUnderflowAndOverflow();
  TH1F* histo_jetfake_sidebandDown = (TH1F*)((TH1F*)f_jetfake->Get(histname_qcd_sidebandDown))->Clone("histo_jetfake_sidebandDown");
  histo_jetfake_sidebandDown->SetBinContent(nBins, histo_jetfake_sidebandDown->GetBinContent(nBins)+histo_jetfake_sidebandDown->GetBinContent(nBins+1));
  histo_jetfake_sidebandDown->ClearUnderflowAndOverflow();
  TH1F* histo_jetfake_METUp = (TH1F*)((TH1F*)f_jetfake->Get(histname_qcd_METUp))->Clone("histo_jetfake_METUp");
  histo_jetfake_METUp->SetBinContent(nBins, histo_jetfake_METUp->GetBinContent(nBins)+histo_jetfake_METUp->GetBinContent(nBins+1));
  histo_jetfake_METUp->ClearUnderflowAndOverflow();
  TH1F* histo_jetfake_METDown = (TH1F*)((TH1F*)f_jetfake->Get(histname_qcd_METDown))->Clone("histo_jetfake_METDown");
  histo_jetfake_METDown->SetBinContent(nBins, histo_jetfake_METDown->GetBinContent(nBins)+histo_jetfake_METDown->GetBinContent(nBins+1));
  histo_jetfake_METDown->ClearUnderflowAndOverflow();
  TH1F* histo_jetfake_binningUp = (TH1F*)((TH1F*)f_jetfake->Get(histname_qcd_binningUp))->Clone("histo_jetfake_binningUp");
  histo_jetfake_binningUp->SetBinContent(nBins, histo_jetfake_binningUp->GetBinContent(nBins)+histo_jetfake_binningUp->GetBinContent(nBins+1));
  histo_jetfake_binningUp->ClearUnderflowAndOverflow();
  TH1F* histo_jetfake_binningDown = (TH1F*)((TH1F*)f_jetfake->Get(histname_qcd_binningDown))->Clone("histo_jetfake_binningDown");
  histo_jetfake_binningDown->SetBinContent(nBins, histo_jetfake_binningDown->GetBinContent(nBins)+histo_jetfake_binningDown->GetBinContent(nBins+1));
  histo_jetfake_binningDown->ClearUnderflowAndOverflow();
  TH1F* histo_jetfake_sieieLeft = (TH1F*)((TH1F*)f_jetfake->Get(histname_qcd_sieieLeft))->Clone("histo_jetfake_sieieLeft");
  histo_jetfake_sieieLeft->SetBinContent(nBins, histo_jetfake_sieieLeft->GetBinContent(nBins)+histo_jetfake_sieieLeft->GetBinContent(nBins+1));
  histo_jetfake_sieieLeft->ClearUnderflowAndOverflow();
  TH1F* histo_jetfake_sieieRight = (TH1F*)((TH1F*)f_jetfake->Get(histname_qcd_sieieRight))->Clone("histo_jetfake_sieieRight");
  histo_jetfake_sieieRight->SetBinContent(nBins, histo_jetfake_sieieRight->GetBinContent(nBins)+histo_jetfake_sieieRight->GetBinContent(nBins+1));
  histo_jetfake_sieieRight->ClearUnderflowAndOverflow();
  TH1F* histo_jetfake_templateUp = (TH1F*)((TH1F*)f_jetfake->Get(histname_qcd_templateUp))->Clone("histo_jetfake_templateUp");
  histo_jetfake_templateUp->SetBinContent(nBins, histo_jetfake_templateUp->GetBinContent(nBins)+histo_jetfake_templateUp->GetBinContent(nBins+1));
  histo_jetfake_templateUp->ClearUnderflowAndOverflow();
  TH1F* histo_jetfake_templateDown = (TH1F*)((TH1F*)f_jetfake->Get(histname_qcd_templateDown))->Clone("histo_jetfake_templateDown");
  histo_jetfake_templateDown->SetBinContent(nBins, histo_jetfake_templateDown->GetBinContent(nBins)+histo_jetfake_templateDown->GetBinContent(nBins+1));
  histo_jetfake_templateDown->ClearUnderflowAndOverflow();
  TH1F* histo_jetfake_unweighted = (TH1F*)((TH1F*)f_jetfake->Get(histname_qcd_unweighted))->Clone("histo_jetfake_unweighted");
  histo_jetfake_unweighted->SetBinContent(nBins, histo_jetfake_unweighted->GetBinContent(nBins)+histo_jetfake_unweighted->GetBinContent(nBins+1));
  histo_jetfake_unweighted->ClearUnderflowAndOverflow();
  TH1F* histo_jetfake_errUp = (TH1F*)histo_jetfake->Clone("histo_jetfake_errUp");
  TH1F* histo_jetfake_errDown = (TH1F*)histo_jetfake->Clone("histo_jetfake_errDown");
  Float_t int_jetfake = histo_jetfake->Integral();
  Float_t max_int_jetfake = 0.0;
  Float_t min_int_jetfake = 0.0;
  Float_t stat_jetfake = 0.0;
  for(int i = 1; i <= nBins; i++){
    double int_bin_jetfake = histo_jetfake->GetBinContent(i);
    double int_bin_jetfake_sidebandUp = histo_jetfake_sidebandUp->GetBinContent(i);
    double int_bin_jetfake_sidebandDown = histo_jetfake_sidebandDown->GetBinContent(i);
    double int_bin_jetfake_METUp = histo_jetfake_METUp->GetBinContent(i);
    double int_bin_jetfake_METDown = histo_jetfake_METDown->GetBinContent(i);
    double int_bin_jetfake_binningUp = histo_jetfake_binningUp->GetBinContent(i);
    double int_bin_jetfake_binningDown = histo_jetfake_binningDown->GetBinContent(i);
    double int_bin_jetfake_sieieLeft = histo_jetfake_sieieLeft->GetBinContent(i);
    double int_bin_jetfake_sieieRight = histo_jetfake_sieieRight->GetBinContent(i);
    double int_bin_jetfake_templateUp = histo_jetfake_templateUp->GetBinContent(i);
    double int_bin_jetfake_templateDown = histo_jetfake_templateDown->GetBinContent(i);
    double int_bin_jetfake_unweighted = histo_jetfake_unweighted->GetBinContent(i);
    double ints_bin[] = {int_bin_jetfake, int_bin_jetfake_sidebandUp, int_bin_jetfake_sidebandDown, int_bin_jetfake_METUp, int_bin_jetfake_METDown, int_bin_jetfake_binningUp, int_bin_jetfake_binningDown, int_bin_jetfake_sieieLeft, int_bin_jetfake_sieieRight, int_bin_jetfake_templateUp, int_bin_jetfake_templateDown};
    double max_int_bin = *max_element(ints_bin, ints_bin+11);
    double min_int_bin = *min_element(ints_bin, ints_bin+11);
    histo_jetfake_errUp->SetBinContent(i, max_int_bin);
    histo_jetfake_errDown->SetBinContent(i, min_int_bin);
    max_int_jetfake += max_int_bin;
    min_int_jetfake += min_int_bin;
    syst_shiftUp_jetfake[i-1] = max_int_bin-int_bin_jetfake;
    syst_shiftDown_jetfake[i-1] = min_int_bin-int_bin_jetfake;
    double stat_bin_jetfake = 0.0;
    if (int_bin_jetfake_unweighted > 0)
      stat_bin_jetfake = int_bin_jetfake/sqrt(int_bin_jetfake_unweighted);
    histo_jetfake->SetBinError(i, stat_bin_jetfake);
    stat_jetfake += stat_bin_jetfake*stat_bin_jetfake;
  }
  Float_t syst_jetfake = TMath::Max(max_int_jetfake-int_jetfake, int_jetfake-min_int_jetfake);
  stat_jetfake = sqrt(stat_jetfake);
  Float_t err_jetfake = sqrt(syst_jetfake*syst_jetfake + stat_jetfake*stat_jetfake);
  total_background += int_jetfake;
  background_unc_sumsquares += err_jetfake*err_jetfake;
  histo_jetfake->SetFillColor(kBlue-4);
  histo_vector.push_back(histo_jetfake);
    
  TFile *f_ZllG_130_under300 = new TFile("ZeeG_JESPES_ZLLGJets_130_under300.root");
  // TFile *f_ZllG_130_over300 = new TFile("ZeeG_JESPES_ZLLGJets_130_over300.root");
  TFile *f_ZllG_300 = new TFile("ZeeG_JESPES_ZLLGJets_300.root");
  TH1F* histo_ZllG_130_under300 = (TH1F*)((TH1F*)f_ZllG_130_under300->Get(histname))->Clone("histo_ZllG_130_under300");
  histo_ZllG_130_under300->SetBinContent(nBins, histo_ZllG_130_under300->GetBinContent(nBins)+histo_ZllG_130_under300->GetBinContent(nBins+1));
  histo_ZllG_130_under300->ClearUnderflowAndOverflow();
  TH1F* histo_ZllG_JESUp_130_under300 = (TH1F*)((TH1F*)f_ZllG_130_under300->Get(histname_JESUp))->Clone("histo_ZllG_JESUp_130_under300");
  histo_ZllG_JESUp_130_under300->SetBinContent(nBins, histo_ZllG_JESUp_130_under300->GetBinContent(nBins)+histo_ZllG_JESUp_130_under300->GetBinContent(nBins+1));
  histo_ZllG_JESUp_130_under300->ClearUnderflowAndOverflow();
  TH1F* histo_ZllG_JESDown_130_under300 = (TH1F*)((TH1F*)f_ZllG_130_under300->Get(histname_JESDown))->Clone("histo_ZllG_JESDown_130_under300");
  histo_ZllG_JESDown_130_under300->SetBinContent(nBins, histo_ZllG_JESDown_130_under300->GetBinContent(nBins)+histo_ZllG_JESDown_130_under300->GetBinContent(nBins+1));
  histo_ZllG_JESDown_130_under300->ClearUnderflowAndOverflow();
  TH1F* histo_ZllG_PESUp_130_under300 = (TH1F*)((TH1F*)f_ZllG_130_under300->Get(histname_PESUp))->Clone("histo_ZllG_PESUp_130_under300");
  histo_ZllG_PESUp_130_under300->SetBinContent(nBins, histo_ZllG_PESUp_130_under300->GetBinContent(nBins)+histo_ZllG_PESUp_130_under300->GetBinContent(nBins+1));
  histo_ZllG_PESUp_130_under300->ClearUnderflowAndOverflow();
  TH1F* histo_ZllG_PESDown_130_under300 = (TH1F*)((TH1F*)f_ZllG_130_under300->Get(histname_PESDown))->Clone("histo_ZllG_PESDown_130_under300");
  histo_ZllG_PESDown_130_under300->SetBinContent(nBins, histo_ZllG_PESDown_130_under300->GetBinContent(nBins)+histo_ZllG_PESDown_130_under300->GetBinContent(nBins+1));
  histo_ZllG_PESDown_130_under300->ClearUnderflowAndOverflow();
  TH1F* histo_ZllG_straightUp_130_under300 = (TH1F*)((TH1F*)f_ZllG_130_under300->Get(histname_straightUp))->Clone("histo_ZllG_straightUp_130_under300");
  histo_ZllG_straightUp_130_under300->SetBinContent(nBins, histo_ZllG_straightUp_130_under300->GetBinContent(nBins)+histo_ZllG_straightUp_130_under300->GetBinContent(nBins+1));
  histo_ZllG_straightUp_130_under300->ClearUnderflowAndOverflow();
  TH1F* histo_ZllG_straightDown_130_under300 = (TH1F*)((TH1F*)f_ZllG_130_under300->Get(histname_straightDown))->Clone("histo_ZllG_straightDown_130_under300");
  histo_ZllG_straightDown_130_under300->SetBinContent(nBins, histo_ZllG_straightDown_130_under300->GetBinContent(nBins)+histo_ZllG_straightDown_130_under300->GetBinContent(nBins+1));
  histo_ZllG_straightDown_130_under300->ClearUnderflowAndOverflow();
  TH1F* histo_ZllG_twistedUp_130_under300 = (TH1F*)((TH1F*)f_ZllG_130_under300->Get(histname_twistedUp))->Clone("histo_ZllG_twistedUp_130_under300");
  histo_ZllG_twistedUp_130_under300->SetBinContent(nBins, histo_ZllG_twistedUp_130_under300->GetBinContent(nBins)+histo_ZllG_twistedUp_130_under300->GetBinContent(nBins+1));
  histo_ZllG_twistedUp_130_under300->ClearUnderflowAndOverflow();
  TH1F* histo_ZllG_twistedDown_130_under300 = (TH1F*)((TH1F*)f_ZllG_130_under300->Get(histname_twistedDown))->Clone("histo_ZllG_twistedDown_130_under300");
  histo_ZllG_twistedDown_130_under300->SetBinContent(nBins, histo_ZllG_twistedDown_130_under300->GetBinContent(nBins)+histo_ZllG_twistedDown_130_under300->GetBinContent(nBins+1));
  histo_ZllG_twistedDown_130_under300->ClearUnderflowAndOverflow();
  TH1F* histo_ZllG_gammaUp_130_under300 = (TH1F*)((TH1F*)f_ZllG_130_under300->Get(histname_gammaUp))->Clone("histo_ZllG_gammaUp_130_under300");
  histo_ZllG_gammaUp_130_under300->SetBinContent(nBins, histo_ZllG_gammaUp_130_under300->GetBinContent(nBins)+histo_ZllG_gammaUp_130_under300->GetBinContent(nBins+1));
  histo_ZllG_gammaUp_130_under300->ClearUnderflowAndOverflow();
  TH1F* histo_ZllG_gammaDown_130_under300 = (TH1F*)((TH1F*)f_ZllG_130_under300->Get(histname_gammaDown))->Clone("histo_ZllG_gammaDown_130_under300");
  histo_ZllG_gammaDown_130_under300->SetBinContent(nBins, histo_ZllG_gammaDown_130_under300->GetBinContent(nBins)+histo_ZllG_gammaDown_130_under300->GetBinContent(nBins+1));
  histo_ZllG_gammaDown_130_under300->ClearUnderflowAndOverflow();
  TH1F* histo_ZllG_qcdscale_130_under300 = (TH1F*)((TH1F*)f_ZllG_130_under300->Get(histname_qcdscale))->Clone("histo_ZllG_qcdscale_130_under300");
  histo_ZllG_qcdscale_130_under300->SetBinContent(nBins, histo_ZllG_qcdscale_130_under300->GetBinContent(nBins)+histo_ZllG_qcdscale_130_under300->GetBinContent(nBins+1));
  histo_ZllG_qcdscale_130_under300->ClearUnderflowAndOverflow();
  TH1F* histo_ZllG_300 = (TH1F*)((TH1F*)f_ZllG_300->Get(histname))->Clone("histo_ZllG_300");
  histo_ZllG_300->SetBinContent(nBins, histo_ZllG_300->GetBinContent(nBins)+histo_ZllG_300->GetBinContent(nBins+1));
  histo_ZllG_300->ClearUnderflowAndOverflow();
  TH1F* histo_ZllG_JESUp_300 = (TH1F*)((TH1F*)f_ZllG_300->Get(histname_JESUp))->Clone("histo_ZllG_JESUp_300");
  histo_ZllG_JESUp_300->SetBinContent(nBins, histo_ZllG_JESUp_300->GetBinContent(nBins)+histo_ZllG_JESUp_300->GetBinContent(nBins+1));
  histo_ZllG_JESUp_300->ClearUnderflowAndOverflow();
  TH1F* histo_ZllG_JESDown_300 = (TH1F*)((TH1F*)f_ZllG_300->Get(histname_JESDown))->Clone("histo_ZllG_JESDown_300");
  histo_ZllG_JESDown_300->SetBinContent(nBins, histo_ZllG_JESDown_300->GetBinContent(nBins)+histo_ZllG_JESDown_300->GetBinContent(nBins+1));
  histo_ZllG_JESDown_300->ClearUnderflowAndOverflow();
  TH1F* histo_ZllG_PESUp_300 = (TH1F*)((TH1F*)f_ZllG_300->Get(histname_PESUp))->Clone("histo_ZllG_PESUp_300");
  histo_ZllG_PESUp_300->SetBinContent(nBins, histo_ZllG_PESUp_300->GetBinContent(nBins)+histo_ZllG_PESUp_300->GetBinContent(nBins+1));
  histo_ZllG_PESUp_300->ClearUnderflowAndOverflow();
  TH1F* histo_ZllG_PESDown_300 = (TH1F*)((TH1F*)f_ZllG_300->Get(histname_PESDown))->Clone("histo_ZllG_PESDown_300");
  histo_ZllG_PESDown_300->SetBinContent(nBins, histo_ZllG_PESDown_300->GetBinContent(nBins)+histo_ZllG_PESDown_300->GetBinContent(nBins+1));
  histo_ZllG_PESDown_300->ClearUnderflowAndOverflow();
  TH1F* histo_ZllG_straightUp_300 = (TH1F*)((TH1F*)f_ZllG_300->Get(histname_straightUp))->Clone("histo_ZllG_straightUp_300");
  histo_ZllG_straightUp_300->SetBinContent(nBins, histo_ZllG_straightUp_300->GetBinContent(nBins)+histo_ZllG_straightUp_300->GetBinContent(nBins+1));
  histo_ZllG_straightUp_300->ClearUnderflowAndOverflow();
  TH1F* histo_ZllG_straightDown_300 = (TH1F*)((TH1F*)f_ZllG_300->Get(histname_straightDown))->Clone("histo_ZllG_straightDown_300");
  histo_ZllG_straightDown_300->SetBinContent(nBins, histo_ZllG_straightDown_300->GetBinContent(nBins)+histo_ZllG_straightDown_300->GetBinContent(nBins+1));
  histo_ZllG_straightDown_300->ClearUnderflowAndOverflow();
  TH1F* histo_ZllG_twistedUp_300 = (TH1F*)((TH1F*)f_ZllG_300->Get(histname_twistedUp))->Clone("histo_ZllG_twistedUp_300");
  histo_ZllG_twistedUp_300->SetBinContent(nBins, histo_ZllG_twistedUp_300->GetBinContent(nBins)+histo_ZllG_twistedUp_300->GetBinContent(nBins+1));
  histo_ZllG_twistedUp_300->ClearUnderflowAndOverflow();
  TH1F* histo_ZllG_twistedDown_300 = (TH1F*)((TH1F*)f_ZllG_300->Get(histname_twistedDown))->Clone("histo_ZllG_twistedDown_300");
  histo_ZllG_twistedDown_300->SetBinContent(nBins, histo_ZllG_twistedDown_300->GetBinContent(nBins)+histo_ZllG_twistedDown_300->GetBinContent(nBins+1));
  histo_ZllG_twistedDown_300->ClearUnderflowAndOverflow();
  TH1F* histo_ZllG_gammaUp_300 = (TH1F*)((TH1F*)f_ZllG_300->Get(histname_gammaUp))->Clone("histo_ZllG_gammaUp_300");
  histo_ZllG_gammaUp_300->SetBinContent(nBins, histo_ZllG_gammaUp_300->GetBinContent(nBins)+histo_ZllG_gammaUp_300->GetBinContent(nBins+1));
  histo_ZllG_gammaUp_300->ClearUnderflowAndOverflow();
  TH1F* histo_ZllG_gammaDown_300 = (TH1F*)((TH1F*)f_ZllG_300->Get(histname_gammaDown))->Clone("histo_ZllG_gammaDown_300");
  histo_ZllG_gammaDown_300->SetBinContent(nBins, histo_ZllG_gammaDown_300->GetBinContent(nBins)+histo_ZllG_gammaDown_300->GetBinContent(nBins+1));
  histo_ZllG_gammaDown_300->ClearUnderflowAndOverflow();
  TH1F* histo_ZllG_qcdscale_300 = (TH1F*)((TH1F*)f_ZllG_300->Get(histname_qcdscale))->Clone("histo_ZllG_qcdscale_300");
  histo_ZllG_qcdscale_300->SetBinContent(nBins, histo_ZllG_qcdscale_300->GetBinContent(nBins)+histo_ZllG_qcdscale_300->GetBinContent(nBins+1));
  histo_ZllG_qcdscale_300->ClearUnderflowAndOverflow();
  histo_ZllG_130_under300->Scale(int_lumi*scale_factor*0.148/489423.0);
  histo_ZllG_JESUp_130_under300->Scale(int_lumi*scale_factor*0.148/489423.0);
  histo_ZllG_JESDown_130_under300->Scale(int_lumi*scale_factor*0.148/489423.0);
  histo_ZllG_PESUp_130_under300->Scale(int_lumi*scale_factor*0.148/489423.0);
  histo_ZllG_PESDown_130_under300->Scale(int_lumi*scale_factor*0.148/489423.0);
  histo_ZllG_straightUp_130_under300->Scale(int_lumi*scale_factor*0.148/489423.0);
  histo_ZllG_straightDown_130_under300->Scale(int_lumi*scale_factor*0.148/489423.0);
  histo_ZllG_twistedUp_130_under300->Scale(int_lumi*scale_factor*0.148/489423.0);
  histo_ZllG_twistedDown_130_under300->Scale(int_lumi*scale_factor*0.148/489423.0);
  histo_ZllG_gammaUp_130_under300->Scale(int_lumi*scale_factor*0.148/489423.0);
  histo_ZllG_gammaDown_130_under300->Scale(int_lumi*scale_factor*0.148/489423.0);
  histo_ZllG_qcdscale_130_under300->Scale(int_lumi*scale_factor*0.148/489423.0);
  histo_ZllG_300->Scale(int_lumi*scale_factor*0.0092/1707773.0);
  histo_ZllG_JESUp_300->Scale(int_lumi*scale_factor*0.0092/1707773.0);
  histo_ZllG_JESDown_300->Scale(int_lumi*scale_factor*0.0092/1707773.0);
  histo_ZllG_PESUp_300->Scale(int_lumi*scale_factor*0.0092/1707773.0);
  histo_ZllG_PESDown_300->Scale(int_lumi*scale_factor*0.0092/1707773.0);
  histo_ZllG_straightUp_300->Scale(int_lumi*scale_factor*0.0092/1707773.0);
  histo_ZllG_straightDown_300->Scale(int_lumi*scale_factor*0.0092/1707773.0);
  histo_ZllG_twistedUp_300->Scale(int_lumi*scale_factor*0.0092/1707773.0);
  histo_ZllG_twistedDown_300->Scale(int_lumi*scale_factor*0.0092/1707773.0);
  histo_ZllG_gammaUp_300->Scale(int_lumi*scale_factor*0.0092/1707773.0);
  histo_ZllG_gammaDown_300->Scale(int_lumi*scale_factor*0.0092/1707773.0);
  histo_ZllG_qcdscale_300->Scale(int_lumi*scale_factor*0.0092/1707773.0);
  TH1F* histo_ZllG_combined = (TH1F*)histo_ZllG_300->Clone("histo_ZllG_combined");
  TH1F* histo_ZllG_JESUp_combined = (TH1F*)histo_ZllG_JESUp_300->Clone("histo_ZllG_JESUp_combined");
  TH1F* histo_ZllG_JESDown_combined = (TH1F*)histo_ZllG_JESDown_300->Clone("histo_ZllG_JESDown_combined");
  TH1F* histo_ZllG_PESUp_combined = (TH1F*)histo_ZllG_PESUp_300->Clone("histo_ZllG_PESUp_combined");
  TH1F* histo_ZllG_PESDown_combined = (TH1F*)histo_ZllG_PESDown_300->Clone("histo_ZllG_PESDown_combined");
  TH1F* histo_ZllG_straightUp_combined = (TH1F*)histo_ZllG_straightUp_300->Clone("histo_ZllG_straightUp_combined");
  TH1F* histo_ZllG_straightDown_combined = (TH1F*)histo_ZllG_straightDown_300->Clone("histo_ZllG_straightDown_combined");
  TH1F* histo_ZllG_twistedUp_combined = (TH1F*)histo_ZllG_twistedUp_300->Clone("histo_ZllG_twistedUp_combined");
  TH1F* histo_ZllG_twistedDown_combined = (TH1F*)histo_ZllG_twistedDown_300->Clone("histo_ZllG_twistedDown_combined");
  TH1F* histo_ZllG_gammaUp_combined = (TH1F*)histo_ZllG_gammaUp_300->Clone("histo_ZllG_gammaUp_combined");
  TH1F* histo_ZllG_gammaDown_combined = (TH1F*)histo_ZllG_gammaDown_300->Clone("histo_ZllG_gammaDown_combined");
  TH1F* histo_ZllG_qcdscale_combined = (TH1F*)histo_ZllG_qcdscale_300->Clone("histo_ZllG_qcdscale_combined");
  histo_ZllG_combined->SetStats(0);
  histo_ZllG_combined->Add(histo_ZllG_130_under300);
  histo_ZllG_JESUp_combined->Add(histo_ZllG_JESUp_130_under300);
  histo_ZllG_JESDown_combined->Add(histo_ZllG_JESDown_130_under300);
  histo_ZllG_PESUp_combined->Add(histo_ZllG_PESUp_130_under300);
  histo_ZllG_PESDown_combined->Add(histo_ZllG_PESDown_130_under300);
  histo_ZllG_straightUp_combined->Add(histo_ZllG_straightUp_130_under300);
  histo_ZllG_straightDown_combined->Add(histo_ZllG_straightDown_130_under300);
  histo_ZllG_twistedUp_combined->Add(histo_ZllG_twistedUp_130_under300);
  histo_ZllG_twistedDown_combined->Add(histo_ZllG_twistedDown_130_under300);
  histo_ZllG_gammaUp_combined->Add(histo_ZllG_gammaUp_130_under300);
  histo_ZllG_gammaDown_combined->Add(histo_ZllG_gammaDown_130_under300);
  histo_ZllG_qcdscale_combined->Add(histo_ZllG_qcdscale_130_under300);
  Float_t int_ZllG_130_under300 = histo_ZllG_130_under300->Integral();
  Float_t int_ZllG_300 = histo_ZllG_300->Integral();
  Float_t int_ZllG = histo_ZllG_combined->Integral();
  Float_t int_ZllG_JESUp = histo_ZllG_JESUp_combined->Integral();
  Float_t int_ZllG_JESDown = histo_ZllG_JESDown_combined->Integral();
  Float_t int_ZllG_PESUp = histo_ZllG_PESUp_combined->Integral();
  Float_t int_ZllG_PESDown = histo_ZllG_PESDown_combined->Integral();
  Float_t int_ZllG_straightUp = histo_ZllG_straightUp_combined->Integral();
  Float_t int_ZllG_straightDown = histo_ZllG_straightDown_combined->Integral();
  Float_t int_ZllG_twistedUp = histo_ZllG_twistedUp_combined->Integral();
  Float_t int_ZllG_twistedDown = histo_ZllG_twistedDown_combined->Integral();
  Float_t int_ZllG_gammaUp = histo_ZllG_gammaUp_combined->Integral();
  Float_t int_ZllG_gammaDown = histo_ZllG_gammaDown_combined->Integral();
  Float_t int_ZllG_qcdscale = histo_ZllG_qcdscale_combined->Integral();
  double jeserr_ZllG = (fabs(int_ZllG_JESUp-int_ZllG)+fabs(int_ZllG_JESDown-int_ZllG))/2.0;
  double peserr_ZllG = (fabs(int_ZllG_PESUp-int_ZllG)+fabs(int_ZllG_PESDown-int_ZllG))/2.0;
  double straighterr_ZllG = (fabs(int_ZllG_straightUp-int_ZllG)+fabs(int_ZllG_straightDown-int_ZllG))/2.0;
  double twistederr_ZllG = (fabs(int_ZllG_twistedUp-int_ZllG)+fabs(int_ZllG_twistedDown-int_ZllG))/2.0;
  double gammaerr_ZllG = (fabs(int_ZllG_gammaUp-int_ZllG)+fabs(int_ZllG_gammaDown-int_ZllG))/2.0;
  double qcdscaleerr_ZllG = fabs(int_ZllG_qcdscale-int_ZllG);
  // TH1F* histo_ZllG_RenUp_130_under300;
  // TH1F* histo_ZllG_RenDown_130_under300;
  // TH1F* histo_ZllG_FacUp_130_under300;
  // TH1F* histo_ZllG_FacDown_130_under300;
  TH1F* histo_ZllG_PDFUp_130_under300;
  TH1F* histo_ZllG_PDFDown_130_under300;
  // TH1F* histo_ZllG_RenUp_130_over300;
  // TH1F* histo_ZllG_RenDown_130_over300;
  // TH1F* histo_ZllG_FacUp_130_over300;
  // TH1F* histo_ZllG_FacDown_130_over300;
  // TH1F* histo_ZllG_PDFUp_130_over300;
  // TH1F* histo_ZllG_PDFDown_130_over300;
  // TH1F* histo_ZllG_RenUp_300;
  // TH1F* histo_ZllG_RenDown_300;
  // TH1F* histo_ZllG_FacUp_300;
  // TH1F* histo_ZllG_FacDown_300;
  TH1F* histo_ZllG_PDFUp_300;
  TH1F* histo_ZllG_PDFDown_300;
  // double renerr_ZllG = 0.0;
  // double facerr_ZllG = 0.0;
  double pdferr_ZllG = 0.0;
  if(histname == "Photon_Et_range_4" || histname == "h_photonic_recoil_4" || histname == "h_phoRecoilMt_4"){
    // These should each already be scaled to the appropriate luminosity, scale factor, cross section, and total number of events
    string xaxis_variable = "XAXIS_VARIABLE_NOT_SET";
    if (histname == "Photon_Et_range_4") xaxis_variable = "Pt";
    else if (histname == "h_photonic_recoil_4") xaxis_variable = "MET";
    else if (histname == "h_phoRecoilMt_4") xaxis_variable = "Mt";
    TString filename = TString("histos_"+xaxis_variable+"_ZeeG_pdfscale_ZLLGJets_130_under300.root");
    TString filename_ext = TString("histos_"+xaxis_variable+"_ZeeG_pdfscale_ZLLGJets_300.root");
    TFile* f_pdfscale_ZllG_130_under300 = new TFile(filename);
    // histo_ZllG_RenUp_130_under300 = (TH1F*)f_pdfscale_ZllG_130_under300->Get("h_ZeeG_pdfscale_ZLLGJets_130_under300_renUp");
    // histo_ZllG_RenDown_130_under300 = (TH1F*)f_pdfscale_ZllG_130_under300->Get("h_ZeeG_pdfscale_ZLLGJets_130_under300_renDown");
    // histo_ZllG_FacUp_130_under300 = (TH1F*)f_pdfscale_ZllG_130_under300->Get("h_ZeeG_pdfscale_ZLLGJets_130_under300_facUp");
    // histo_ZllG_FacDown_130_under300 = (TH1F*)f_pdfscale_ZllG_130_under300->Get("h_ZeeG_pdfscale_ZLLGJets_130_under300_facDown");
    histo_ZllG_PDFUp_130_under300 = (TH1F*)f_pdfscale_ZllG_130_under300->Get("h_ZeeG_pdfscale_ZLLGJets_130_under300_pdfUp");
    histo_ZllG_PDFDown_130_under300 = (TH1F*)f_pdfscale_ZllG_130_under300->Get("h_ZeeG_pdfscale_ZLLGJets_130_under300_pdfDown");
    TFile* f_pdfscale_ZllG_300 = new TFile(filename_ext);
    // histo_ZllG_RenUp_300 = (TH1F*)f_pdfscale_ZllG_300->Get("h_ZeeG_pdfscale_ZLLGJets_300_renUp");
    // histo_ZllG_RenDown_300 = (TH1F*)f_pdfscale_ZllG_300->Get("h_ZeeG_pdfscale_ZLLGJets_300_renDown");
    // histo_ZllG_FacUp_300 = (TH1F*)f_pdfscale_ZllG_300->Get("h_ZeeG_pdfscale_ZLLGJets_300_facUp");
    // histo_ZllG_FacDown_300 = (TH1F*)f_pdfscale_ZllG_300->Get("h_ZeeG_pdfscale_ZLLGJets_300_facDown");
    histo_ZllG_PDFUp_300 = (TH1F*)f_pdfscale_ZllG_300->Get("h_ZeeG_pdfscale_ZLLGJets_300_pdfUp");
    histo_ZllG_PDFDown_300 = (TH1F*)f_pdfscale_ZllG_300->Get("h_ZeeG_pdfscale_ZLLGJets_300_pdfDown");
    // TFile* f_pdfscale_ZllG_130_over300 = new TFile("histos_ZeeG_pdfscale_ZLLGJets_130_over300.root");
    // histo_ZllG_RenUp_130_over300 = (TH1F*)f_pdfscale_ZllG_130_over300->Get("h_ZeeG_pdfscale_ZLLGJets_130_over300_renUp");
    // histo_ZllG_RenDown_130_over300 = (TH1F*)f_pdfscale_ZllG_130_over300->Get("h_ZeeG_pdfscale_ZLLGJets_130_over300_renDown");
    // histo_ZllG_FacUp_130_over300 = (TH1F*)f_pdfscale_ZllG_130_over300->Get("h_ZeeG_pdfscale_ZLLGJets_130_over300_facUp");
    // histo_ZllG_FacDown_130_over300 = (TH1F*)f_pdfscale_ZllG_130_over300->Get("h_ZeeG_pdfscale_ZLLGJets_130_over300_facDown");
    // histo_ZllG_PDFUp_130_over300 = (TH1F*)f_pdfscale_ZllG_130_over300->Get("h_ZeeG_pdfscale_ZLLGJets_130_over300_pdfUp");
    // histo_ZllG_PDFDown_130_over300 = (TH1F*)f_pdfscale_ZllG_130_over300->Get("h_ZeeG_pdfscale_ZLLGJets_130_over300_pdfDown");
    // histo_ZllG_RenUp_300->Add(histo_ZllG_RenUp_130_under300);
    // histo_ZllG_RenDown_300->Add(histo_ZllG_RenDown_130_under300);
    // histo_ZllG_FacUp_300->Add(histo_ZllG_FacUp_130_under300);
    // histo_ZllG_FacDown_300->Add(histo_ZllG_FacDown_130_under300);
    histo_ZllG_PDFUp_300->Add(histo_ZllG_PDFUp_130_under300);
    histo_ZllG_PDFDown_300->Add(histo_ZllG_PDFDown_130_under300);
    // Float_t int_ZllG_RenUp = histo_ZllG_RenUp_300->Integral();
    // Float_t int_ZllG_RenDown = histo_ZllG_RenDown_300->Integral();
    // Float_t int_ZllG_FacUp = histo_ZllG_FacUp_300->Integral();
    // Float_t int_ZllG_FacDown = histo_ZllG_FacDown_300->Integral();
    Float_t int_ZllG_PDFUp = histo_ZllG_PDFUp_300->Integral();
    Float_t int_ZllG_PDFDown = histo_ZllG_PDFDown_300->Integral();
    // renerr_ZllG = (fabs(int_ZllG_RenUp-int_ZllG)+fabs(int_ZllG_RenDown-int_ZllG))/2.0;
    // facerr_ZllG = (fabs(int_ZllG_FacUp-int_ZllG)+fabs(int_ZllG_FacDown-int_ZllG))/2.0;
    pdferr_ZllG = (fabs(int_ZllG_PDFUp-int_ZllG)+fabs(int_ZllG_PDFDown-int_ZllG))/2.0;
  }
  for(int i = 1; i <= nBins; i++){
    double int_bin_130_under300 = histo_ZllG_130_under300->GetBinContent(i);
    double int_bin_300 = histo_ZllG_300->GetBinContent(i);
    double int_bin = histo_ZllG_combined->GetBinContent(i);
    double jesup = histo_ZllG_JESUp_combined->GetBinContent(i);
    double jesdown = histo_ZllG_JESDown_combined->GetBinContent(i);
    double pesup = histo_ZllG_PESUp_combined->GetBinContent(i);
    double pesdown = histo_ZllG_PESDown_combined->GetBinContent(i);
    double straightup = histo_ZllG_straightUp_combined->GetBinContent(i);
    double straightdown = histo_ZllG_straightDown_combined->GetBinContent(i);
    double twistedup = histo_ZllG_twistedUp_combined->GetBinContent(i);
    double twisteddown = histo_ZllG_twistedDown_combined->GetBinContent(i);
    double gammaup = histo_ZllG_gammaUp_combined->GetBinContent(i);
    double gammadown = histo_ZllG_gammaDown_combined->GetBinContent(i);
    jesup_shift[i-1] += jesup-int_bin;
    jesdown_shift[i-1] += jesdown-int_bin;
    pesup_shift[i-1] += pesup-int_bin;
    pesdown_shift[i-1] += pesdown-int_bin;
    straightup_shift_ZllG[i-1] += straightup-int_bin;
    straightdown_shift_ZllG[i-1] += straightdown-int_bin;
    twistedup_shift_ZllG[i-1] += twistedup-int_bin;
    twisteddown_shift_ZllG[i-1] += twisteddown-int_bin;
    gammaup_shift_ZllG[i-1] += gammaup-int_bin;
    gammadown_shift_ZllG[i-1] += gammadown-int_bin;
    double qcdscale = histo_ZllG_qcdscale_combined->GetBinContent(i);
    qcdscale_shift_ZllG[i-1] = fabs(qcdscale-int_bin);
    if(histname == "Photon_Et_range_4" || histname == "h_photonic_recoil_4" || histname == "h_phoRecoilMt_4"){
      // double renup = histo_ZllG_RenUp_300->GetBinContent(i);
      // double rendown = histo_ZllG_RenDown_300->GetBinContent(i);
      // double facup = histo_ZllG_FacUp_300->GetBinContent(i);
      // double facdown = histo_ZllG_FacDown_300->GetBinContent(i);
      double pdfup = histo_ZllG_PDFUp_300->GetBinContent(i);
      double pdfdown = histo_ZllG_PDFDown_300->GetBinContent(i);
      // renup_shift[i-1] += renup-int_bin;
      // rendown_shift[i-1] += rendown-int_bin;
      // facup_shift[i-1] += facup-int_bin;
      // facdown_shift[i-1] += facdown-int_bin;
      pdfup_shift[i-1] += pdfup-int_bin;
      pdfdown_shift[i-1] += pdfdown-int_bin;
    }
    // cout<<"ZllG: bin="<<i<<",int_bin="<<int_bin<<",jesup="<<jesup<<",jesdown="<<jesdown<<",pesup="<<pesup<<",pesdown="<<pesdown<<endl;
    
    double err_bin = 0.0;
    if(int_bin > 0)
      err_bin = sqrt(int_bin_130_under300*(0.148*int_lumi*scale_factor-int_bin_130_under300)/489423.0 + int_bin_300*(0.0092*int_lumi*scale_factor-int_bin_300)/1707773.0);
    histo_ZllG_combined->SetBinError(i,err_bin);
  }
  Float_t err_ZllG = 0.0;
  if(histname == "Photon_Et_range_4"){
    cout<<"ZllG fractional errors"<<endl;
    cout<<"ZllG_130_under300 stat err: "<<sqrt((0.148*int_lumi*scale_factor-int_ZllG_130_under300)/(489423.0*int_ZllG_130_under300))/int_ZllG_130_under300<<endl;
    cout<<"ZllG_300 stat err: "<<sqrt((0.0092*int_lumi*scale_factor-int_ZllG_300)/(1707773.0*int_ZllG_300))/int_ZllG_300<<endl;
    cout<<"phoSF err: "<<photon_scale_factor_unc<<endl;
    cout<<"qcdscaleerr_ZllG: "<<qcdscaleerr_ZllG/int_ZllG<<endl;
    cout<<"jeserr_ZllG: "<<jeserr_ZllG/int_ZllG<<endl;
    cout<<"peserr_ZllG: "<<peserr_ZllG/int_ZllG<<endl;
    cout<<"straighterr_ZllG: "<<straighterr_ZllG/int_ZllG<<endl;
    cout<<"twistederr_ZllG: "<<twistederr_ZllG/int_ZllG<<endl;
    cout<<"gammaerr_ZllG: "<<gammaerr_ZllG/int_ZllG<<endl;
    cout<<"pdferr_ZllG: "<<pdferr_ZllG/int_ZllG<<endl;
  }
  if(int_ZllG > 0.0)
    err_ZllG = sqrt(int_ZllG_130_under300*int_ZllG_130_under300*((photon_scale_factor_unc*photon_scale_factor_unc)+(electron_scale_factor_unc*electron_scale_factor_unc)+(0.148*int_lumi*scale_factor-int_ZllG_130_under300)/(489423.0*int_ZllG_130_under300))+int_ZllG_300*int_ZllG_300*((photon_scale_factor_unc*photon_scale_factor_unc)+(electron_scale_factor_unc*electron_scale_factor_unc)+(0.0092*int_lumi*scale_factor-int_ZllG_300)/(1707773.0*int_ZllG_300))+(qcdscaleerr_ZllG*qcdscaleerr_ZllG)+(jeserr_ZllG*jeserr_ZllG)+(peserr_ZllG*peserr_ZllG)+(straighterr_ZllG*straighterr_ZllG)+(twistederr_ZllG*twistederr_ZllG)+(gammaerr_ZllG*gammaerr_ZllG)+(pdferr_ZllG*pdferr_ZllG));
  total_background += int_ZllG;
  background_unc_sumsquares += err_ZllG*err_ZllG;
  histo_ZllG_combined->SetFillColor(kSpring-9);
  histo_vector.push_back(histo_ZllG_combined);
    
  //DEBUG
  // cout<<"ZllG got"<<endl;
  
  TFile *f_TTG = new TFile("ZeeG_JESPES_TTGJets.root");
  TH1F* histo_TTG = (TH1F*)((TH1F*)f_TTG->Get(histname))->Clone("histo_TTG");
  histo_TTG->SetBinContent(nBins, histo_TTG->GetBinContent(nBins)+histo_TTG->GetBinContent(nBins+1));
  histo_TTG->ClearUnderflowAndOverflow();
  TH1F* histo_TTG_JESUp = (TH1F*)((TH1F*)f_TTG->Get(histname_JESUp))->Clone("histo_TTG_JESUp");
  histo_TTG_JESUp->SetBinContent(nBins, histo_TTG_JESUp->GetBinContent(nBins)+histo_TTG_JESUp->GetBinContent(nBins+1));
  histo_TTG_JESUp->ClearUnderflowAndOverflow();
  TH1F* histo_TTG_JESDown = (TH1F*)((TH1F*)f_TTG->Get(histname_JESDown))->Clone("histo_TTG_JESDown");
  histo_TTG_JESDown->SetBinContent(nBins, histo_TTG_JESDown->GetBinContent(nBins)+histo_TTG_JESDown->GetBinContent(nBins+1));
  histo_TTG_JESDown->ClearUnderflowAndOverflow();
  TH1F* histo_TTG_PESUp = (TH1F*)((TH1F*)f_TTG->Get(histname_PESUp))->Clone("histo_TTG_PESUp");
  histo_TTG_PESUp->SetBinContent(nBins, histo_TTG_PESUp->GetBinContent(nBins)+histo_TTG_PESUp->GetBinContent(nBins+1));
  histo_TTG_PESUp->ClearUnderflowAndOverflow();
  TH1F* histo_TTG_PESDown = (TH1F*)((TH1F*)f_TTG->Get(histname_PESDown))->Clone("histo_TTG_PESDown");
  histo_TTG_PESDown->SetBinContent(nBins, histo_TTG_PESDown->GetBinContent(nBins)+histo_TTG_PESDown->GetBinContent(nBins+1));
  histo_TTG_PESDown->ClearUnderflowAndOverflow();
  histo_TTG->SetStats(0);
  histo_TTG->Scale(int_lumi*scale_factor*3.697/3170400.0);
  histo_TTG_JESUp->Scale(int_lumi*scale_factor*3.697/3170400.0);
  histo_TTG_JESDown->Scale(int_lumi*scale_factor*3.697/3170400.0);
  histo_TTG_PESUp->Scale(int_lumi*scale_factor*3.697/3170400.0);
  histo_TTG_PESDown->Scale(int_lumi*scale_factor*3.697/3170400.0);
  Float_t int_TTG = histo_TTG->Integral();
  Float_t int_TTG_JESUp = histo_TTG_JESUp->Integral();
  Float_t int_TTG_JESDown = histo_TTG_JESDown->Integral();
  Float_t int_TTG_PESUp = histo_TTG_PESUp->Integral();
  Float_t int_TTG_PESDown = histo_TTG_PESDown->Integral();
  double jeserr_TTG = (fabs(int_TTG_JESUp-int_TTG)+fabs(int_TTG_JESDown-int_TTG))/2.0;
  double peserr_TTG = (fabs(int_TTG_PESUp-int_TTG)+fabs(int_TTG_PESDown-int_TTG))/2.0;
  Float_t err_TTG = 0.0;
  if(int_TTG > 0.0)
    err_TTG = sqrt(int_TTG*int_TTG*((photon_scale_factor_unc*photon_scale_factor_unc)+(electron_scale_factor_unc*electron_scale_factor_unc)+(3.697*int_lumi*scale_factor-int_TTG)/(3170400.0*int_TTG))+(jeserr_TTG*jeserr_TTG)+(peserr_TTG*peserr_TTG));
  total_background += int_TTG;
  background_unc_sumsquares += err_TTG*err_TTG;
  for(int i = 1; i <= nBins; i++){
    double int_bin = histo_TTG->GetBinContent(i);
    double jesup = histo_TTG_JESUp->GetBinContent(i);
    double jesdown = histo_TTG_JESDown->GetBinContent(i);
    double pesup = histo_TTG_PESUp->GetBinContent(i);
    double pesdown = histo_TTG_PESDown->GetBinContent(i);
    jesup_shift[i-1] += jesup-int_bin;
    jesdown_shift[i-1] += jesdown-int_bin;
    pesup_shift[i-1] += pesup-int_bin;
    pesdown_shift[i-1] += pesdown-int_bin;
    // cout<<"TTG: bin="<<i<<",int_bin="<<int_bin<<",jesup="<<jesup<<",jesdown="<<jesdown<<",pesup="<<pesup<<",pesdown="<<pesdown<<endl;
    
    double err_bin = 0.0;
    if(int_bin > 0)
      err_bin = int_bin*sqrt((3.697*int_lumi*scale_factor-int_bin)/(3170400.0*int_bin));
    histo_TTG->SetBinError(i,err_bin);
  }
  histo_TTG->SetFillColor(kOrange-3);
  histo_vector.push_back(histo_TTG);
      
  // TFile* f_WWG = new TFile("ZeeG_JESPES_WWG.root");
  // TH1F* histo_WWG = (TH1F*)f_WWG->Get(histname);
  // TH1F* histo_WWG_JESUp = (TH1F*)f_WWG->Get(histname_JESUp);
  // TH1F* histo_WWG_JESDown = (TH1F*)f_WWG->Get(histname_JESDown);
  // TH1F* histo_WWG_PESUp = (TH1F*)f_WWG->Get(histname_PESUp);
  // TH1F* histo_WWG_PESDown = (TH1F*)f_WWG->Get(histname_PESDown);
  // histo_WWG->SetStats(0);
  // histo_WWG->Scale(int_lumi*scale_factor*0.2147/827604.0);
  // histo_WWG_JESUp->Scale(int_lumi*scale_factor*0.2147/827604.0);
  // histo_WWG_JESDown->Scale(int_lumi*scale_factor*0.2147/827604.0);
  // histo_WWG_PESUp->Scale(int_lumi*scale_factor*0.2147/827604.0);
  // histo_WWG_PESDown->Scale(int_lumi*scale_factor*0.2147/827604.0);
  // Float_t int_WWG = histo_WWG->Integral()+histo_WWG->GetBinContent(0)+histo_WWG->GetBinContent(nBins+1);
  // Float_t int_WWG_JESUp = histo_WWG_JESUp->Integral()+histo_WWG_JESUp->GetBinContent(0)+histo_WWG_JESUp->GetBinContent(nBins+1);
  // Float_t int_WWG_JESDown = histo_WWG_JESDown->Integral()+histo_WWG_JESDown->GetBinContent(0)+histo_WWG_JESDown->GetBinContent(nBins+1);
  // Float_t int_WWG_PESUp = histo_WWG_PESUp->Integral()+histo_WWG_PESUp->GetBinContent(0)+histo_WWG_PESUp->GetBinContent(nBins+1);
  // Float_t int_WWG_PESDown = histo_WWG_PESDown->Integral()+histo_WWG_PESDown->GetBinContent(0)+histo_WWG_PESDown->GetBinContent(nBins+1);
  // double jeserr_WWG = (fabs(int_WWG_JESUp-int_WWG)+fabs(int_WWG_JESDown-int_WWG))/2.0;
  // double peserr_WWG = (fabs(int_WWG_PESUp-int_WWG)+fabs(int_WWG_PESDown-int_WWG))/2.0;
  // Float_t err_WWG = 0.0;
  // if(int_WWG > 0.0)
  //   err_WWG = sqrt(int_WWG*int_WWG*((photon_scale_factor_unc+photon_scale_factor_unc)+(electron_scale_factor_unc+electron_scale_factor_unc)+(0.2147*int_lumi*scale_factor-int_WWG)/(827604.0*int_WWG))+(jeserr_WWG*jeserr_WWG)+(peserr_WWG*peserr_WWG));
  // total_background += int_WWG;
  // background_unc_sumsquares += err_WWG*err_WWG;
  // for(int i = 1; i <= nBins; i++){
  //   double int_bin = histo_WWG->GetBinContent(i);
  //   double jesup = histo_WWG_JESUp->GetBinContent(i);
  //   double jesdown = histo_WWG_JESDown->GetBinContent(i);
  //   double pesup = histo_WWG_PESUp->GetBinContent(i);
  //   double pesdown = histo_WWG_PESDown->GetBinContent(i);
  //   jesup_shift[i-1] += jesup-int_bin;
  //   jesdown_shift[i-1] += jesdown-int_bin;
  //   pesup_shift[i-1] += pesup-int_bin;
  //   pesdown_shift[i-1] += pesdown-int_bin;
  //   // cout<<"WWG: bin="<<i<<",int_bin="<<int_bin<<",jesup="<<jesup<<",jesdown="<<jesdown<<",pesup="<<pesup<<",pesdown="<<pesdown<<endl;
    
  //   double err_bin = 0.0;
  //   if(int_bin > 0)
  //     err_bin = int_bin*sqrt((0.2147*int_lumi*scale_factor-int_bin)/(827604.0*int_bin));
  //   histo_WWG->SetBinError(i,err_bin);
  // }
  // // histo_WWG->SetFillColor(kBlue-4);
  // histo_WWG->SetFillColor(kTeal+3);
  // histo_vector.push_back(histo_WWG);
    
  TFile *f_WZ = new TFile("ZeeG_JESPES_WZ.root");
  TH1F* histo_WZ = (TH1F*)((TH1F*)f_WZ->Get(histname))->Clone("histo_WZ");
  histo_WZ->SetBinContent(nBins, histo_WZ->GetBinContent(nBins)+histo_WZ->GetBinContent(nBins+1));
  histo_WZ->ClearUnderflowAndOverflow();
  TH1F* histo_WZ_JESUp = (TH1F*)((TH1F*)f_WZ->Get(histname_JESUp))->Clone("histo_WZ_JESUp");
  histo_WZ_JESUp->SetBinContent(nBins, histo_WZ_JESUp->GetBinContent(nBins)+histo_WZ_JESUp->GetBinContent(nBins+1));
  histo_WZ_JESUp->ClearUnderflowAndOverflow();
  TH1F* histo_WZ_JESDown = (TH1F*)((TH1F*)f_WZ->Get(histname_JESDown))->Clone("histo_WZ_JESDown");
  histo_WZ_JESDown->SetBinContent(nBins, histo_WZ_JESDown->GetBinContent(nBins)+histo_WZ_JESDown->GetBinContent(nBins+1));
  histo_WZ_JESDown->ClearUnderflowAndOverflow();
  TH1F* histo_WZ_PESUp = (TH1F*)((TH1F*)f_WZ->Get(histname_PESUp))->Clone("histo_WZ_PESUp");
  histo_WZ_PESUp->SetBinContent(nBins, histo_WZ_PESUp->GetBinContent(nBins)+histo_WZ_PESUp->GetBinContent(nBins+1));
  histo_WZ_PESUp->ClearUnderflowAndOverflow();
  TH1F* histo_WZ_PESDown = (TH1F*)((TH1F*)f_WZ->Get(histname_PESDown))->Clone("histo_WZ_PESDown");
  histo_WZ_PESDown->SetBinContent(nBins, histo_WZ_PESDown->GetBinContent(nBins)+histo_WZ_PESDown->GetBinContent(nBins+1));
  histo_WZ_PESDown->ClearUnderflowAndOverflow();
  histo_WZ->SetStats(0);
  histo_WZ->Scale(int_lumi*scale_factor*66.1/2995783.0);
  histo_WZ_JESUp->Scale(int_lumi*scale_factor*66.1/2995783.0);
  histo_WZ_JESDown->Scale(int_lumi*scale_factor*66.1/2995783.0);
  histo_WZ_PESUp->Scale(int_lumi*scale_factor*66.1/2995783.0);
  histo_WZ_PESDown->Scale(int_lumi*scale_factor*66.1/2995783.0);
  Float_t int_WZ = histo_WZ->Integral();
  Float_t int_WZ_JESUp = histo_WZ_JESUp->Integral();
  Float_t int_WZ_JESDown = histo_WZ_JESDown->Integral();
  Float_t int_WZ_PESUp = histo_WZ_PESUp->Integral();
  Float_t int_WZ_PESDown = histo_WZ_PESDown->Integral();
  double jeserr_WZ = (fabs(int_WZ_JESUp-int_WZ)+fabs(int_WZ_JESDown-int_WZ))/2.0;
  double peserr_WZ = (fabs(int_WZ_PESUp-int_WZ)+fabs(int_WZ_PESDown-int_WZ))/2.0;
  Float_t err_WZ = 0.0;
  if(int_WZ > 0.0)
    err_WZ = sqrt(int_WZ*int_WZ*((photon_scale_factor_unc*photon_scale_factor_unc)+(electron_scale_factor_unc*electron_scale_factor_unc)+(66.1*int_lumi*scale_factor-int_WZ)/(2995783.0*int_WZ))+(jeserr_WZ*jeserr_WZ)+(peserr_WZ*peserr_WZ));
  total_background += int_WZ;
  background_unc_sumsquares += err_WZ*err_WZ;
  for(int i = 1; i <= nBins; i++){
    double int_bin = histo_WZ->GetBinContent(i);
    double jesup = histo_WZ_JESUp->GetBinContent(i);
    double jesdown = histo_WZ_JESDown->GetBinContent(i);
    double pesup = histo_WZ_PESUp->GetBinContent(i);
    double pesdown = histo_WZ_PESDown->GetBinContent(i);
    jesup_shift[i-1] += jesup-int_bin;
    jesdown_shift[i-1] += jesdown-int_bin;
    pesup_shift[i-1] += pesup-int_bin;
    pesdown_shift[i-1] += pesdown-int_bin;
    // cout<<"WZ: bin="<<i<<",int_bin="<<int_bin<<",jesup="<<jesup<<",jesdown="<<jesdown<<",pesup="<<pesup<<",pesdown="<<pesdown<<endl;
    
    double err_bin = 0.0;
    if(int_bin > 0)
      err_bin = int_bin*sqrt((66.1*int_lumi*scale_factor-int_bin)/(2995783.0*int_bin));
    histo_WZ->SetBinError(i,err_bin);
  }
  // histo_WZ->SetFillColor(kYellow);
  histo_WZ->SetFillColor(kRed-5);
  histo_vector.push_back(histo_WZ);
  
  //DEBUG
  // cout<<"WZ got"<<endl;
  
  // TFile *f_ZGGToLLGG = new TFile("ZeeG_JESPES_ZGGToLLGG.root");
  // TH1F* histo_ZGGToLLGG = (TH1F*)((TH1F*)f_ZGGToLLGG->Get(histname))->Clone("histo_ZGGToLLGG");
  // histo_ZGGToLLGG->SetBinContent(nBins, histo_ZGGToLLGG->GetBinContent(nBins)+histo_ZGGToLLGG->GetBinContent(nBins+1));
  // histo_ZGGToLLGG->ClearUnderflowAndOverflow();
  // TH1F* histo_ZGGToLLGG_JESUp = (TH1F*)((TH1F*)f_ZGGToLLGG->Get(histname_JESUp))->Clone("histo_ZGGToLLGG_JESUp");
  // histo_ZGGToLLGG_JESUp->SetBinContent(nBins, histo_ZGGToLLGG_JESUp->GetBinContent(nBins)+histo_ZGGToLLGG_JESUp->GetBinContent(nBins+1));
  // histo_ZGGToLLGG_JESUp->ClearUnderflowAndOverflow();
  // TH1F* histo_ZGGToLLGG_JESDown = (TH1F*)((TH1F*)f_ZGGToLLGG->Get(histname_JESDown))->Clone("histo_ZGGToLLGG_JESDown");
  // histo_ZGGToLLGG_JESDown->SetBinContent(nBins, histo_ZGGToLLGG_JESDown->GetBinContent(nBins)+histo_ZGGToLLGG_JESDown->GetBinContent(nBins+1));
  // histo_ZGGToLLGG_JESDown->ClearUnderflowAndOverflow();
  // TH1F* histo_ZGGToLLGG_PESUp = (TH1F*)((TH1F*)f_ZGGToLLGG->Get(histname_PESUp))->Clone("histo_ZGGToLLGG_PESUp");
  // histo_ZGGToLLGG_PESUp->SetBinContent(nBins, histo_ZGGToLLGG_PESUp->GetBinContent(nBins)+histo_ZGGToLLGG_PESUp->GetBinContent(nBins+1));
  // histo_ZGGToLLGG_PESUp->ClearUnderflowAndOverflow();
  // TH1F* histo_ZGGToLLGG_PESDown = (TH1F*)((TH1F*)f_ZGGToLLGG->Get(histname_PESDown))->Clone("histo_ZGGToLLGG_PESDown");
  // histo_ZGGToLLGG_PESDown->SetBinContent(nBins, histo_ZGGToLLGG_PESDown->GetBinContent(nBins)+histo_ZGGToLLGG_PESDown->GetBinContent(nBins+1));
  // histo_ZGGToLLGG_PESDown->ClearUnderflowAndOverflow();
  // TH1F* histo_ZGGToLLGG_phoSFUp = (TH1F*)((TH1F*)f_ZGGToLLGG->Get(histname_phoSFUp))->Clone("histo_ZGGToLLGG_phoSFUp");
  // histo_ZGGToLLGG_phoSFUp->SetBinContent(nBins, histo_ZGGToLLGG_phoSFUp->GetBinContent(nBins)+histo_ZGGToLLGG_phoSFUp->GetBinContent(nBins+1));
  // histo_ZGGToLLGG_phoSFUp->ClearUnderflowAndOverflow();
  // TH1F* histo_ZGGToLLGG_phoSFDown = (TH1F*)((TH1F*)f_ZGGToLLGG->Get(histname_phoSFDown))->Clone("histo_ZGGToLLGG_phoSFDown");
  // histo_ZGGToLLGG_phoSFDown->SetBinContent(nBins, histo_ZGGToLLGG_phoSFDown->GetBinContent(nBins)+histo_ZGGToLLGG_phoSFDown->GetBinContent(nBins+1));
  // histo_ZGGToLLGG_phoSFDown->ClearUnderflowAndOverflow();
  // histo_ZGGToLLGG->SetStats(0);
  // histo_ZGGToLLGG->Scale(int_lumi*scale_factor*0.6289/731966.0);
  // histo_ZGGToLLGG_JESUp->Scale(int_lumi*scale_factor*0.6289/731966.0);
  // histo_ZGGToLLGG_JESDown->Scale(int_lumi*scale_factor*0.6289/731966.0);
  // histo_ZGGToLLGG_PESUp->Scale(int_lumi*scale_factor*0.6289/731966.0);
  // histo_ZGGToLLGG_PESDown->Scale(int_lumi*scale_factor*0.6289/731966.0);
  // histo_ZGGToLLGG_phoSFUp->Scale(int_lumi*scale_factor*0.6289/731966.0);
  // histo_ZGGToLLGG_phoSFDown->Scale(int_lumi*scale_factor*0.6289/731966.0);
  // Float_t int_ZGGToLLGG = histo_ZGGToLLGG->Integral();
  // Float_t int_ZGGToLLGG_JESUp = histo_ZGGToLLGG_JESUp->Integral();
  // Float_t int_ZGGToLLGG_JESDown = histo_ZGGToLLGG_JESDown->Integral();
  // Float_t int_ZGGToLLGG_PESUp = histo_ZGGToLLGG_PESUp->Integral();
  // Float_t int_ZGGToLLGG_PESDown = histo_ZGGToLLGG_PESDown->Integral();
  // Float_t int_ZGGToLLGG_phoSFUp = histo_ZGGToLLGG_phoSFUp->Integral();
  // Float_t int_ZGGToLLGG_phoSFDown = histo_ZGGToLLGG_phoSFDown->Integral();
  // double jeserr_ZGGToLLGG = (fabs(int_ZGGToLLGG_JESUp-int_ZGGToLLGG)+fabs(int_ZGGToLLGG_JESDown-int_ZGGToLLGG))/2.0;
  // double peserr_ZGGToLLGG = (fabs(int_ZGGToLLGG_PESUp-int_ZGGToLLGG)+fabs(int_ZGGToLLGG_PESDown-int_ZGGToLLGG))/2.0;
  // double phosferr_ZGGToLLGG = (fabs(int_ZGGToLLGG_phoSFUp-int_ZGGToLLGG)+fabs(int_ZGGToLLGG_phoSFDown-int_ZGGToLLGG))/2.0;
  // Float_t err_ZGGToLLGG = 0.0;
  // if(int_ZGGToLLGG > 0.0)
  //   err_ZGGToLLGG = sqrt(int_ZGGToLLGG*int_ZGGToLLGG*((photon_scale_factor_unc*photon_scale_factor_unc)+(electron_scale_factor_unc*electron_scale_factor_unc)+(0.6289*int_lumi*scale_factor-int_ZGGToLLGG)/(731966.0*int_ZGGToLLGG))+(jeserr_ZGGToLLGG*jeserr_ZGGToLLGG)+(peserr_ZGGToLLGG*peserr_ZGGToLLGG)+(phosferr_ZGGToLLGG*phosferr_ZGGToLLGG));
  // total_background += int_ZGGToLLGG;
  // background_unc_sumsquares += err_ZGGToLLGG*err_ZGGToLLGG;
  // for(int i = 1; i <= nBins; i++){
  //   double int_bin = histo_ZGGToLLGG->GetBinContent(i);
  //   double jesup = histo_ZGGToLLGG_JESUp->GetBinContent(i);
  //   double jesdown = histo_ZGGToLLGG_JESDown->GetBinContent(i);
  //   double pesup = histo_ZGGToLLGG_PESUp->GetBinContent(i);
  //   double pesdown = histo_ZGGToLLGG_PESDown->GetBinContent(i);
  //   double phosfup = histo_ZGGToLLGG_phoSFUp->GetBinContent(i);
  //   double phosfdown = histo_ZGGToLLGG_phoSFDown->GetBinContent(i);
  //   jesup_shift[i-1] += jesup-int_bin;
  //   jesdown_shift[i-1] += jesdown-int_bin;
  //   pesup_shift[i-1] += pesup-int_bin;
  //   pesdown_shift[i-1] += pesdown-int_bin;
  //   phosfup_shift[i-1] += phosfup-int_bin;
  //   phosfdown_shift[i-1] += phosfdown-int_bin;
  //   // cout<<"ZGGToLLGG: bin="<<i<<",int_bin="<<int_bin<<",jesup="<<jesup<<",jesdown="<<jesdown<<",pesup="<<pesup<<",pesdown="<<pesdown<<endl;
    
  //   double err_bin = 0.0;
  //   if(int_bin > 0)
  //     err_bin = int_bin*sqrt((0.6289*int_lumi*scale_factor-int_bin)/(731966.0*int_bin));
  //   histo_ZGGToLLGG->SetBinError(i,err_bin);
  // }
  // histo_ZGGToLLGG->SetFillColor(kRed-5);
  // histo_vector.push_back(histo_ZGGToLLGG);
  
  // TFile *f_Zll_100to200 = new TFile(TString("ZeeG_JESPES_DYJetsToLL_HT-100to200.root"));
  // TFile *f_Zll_200to400 = new TFile(TString("ZeeG_JESPES_DYJetsToLL_HT-200to400.root"));
  // TFile *f_Zll_400to600 = new TFile(TString("ZeeG_JESPES_DYJetsToLL_HT-400to600.root"));
  // TFile *f_Zll_600to800 = new TFile(TString("ZeeG_JESPES_DYJetsToLL_HT-600to800.root"));
  // TFile *f_Zll_800to1200 = new TFile(TString("ZeeG_JESPES_DYJetsToLL_HT-800to1200.root"));
  // TFile *f_Zll_1200to2500 = new TFile(TString("ZeeG_JESPES_DYJetsToLL_HT-1200to2500.root"));
  // TFile *f_Zll_2500toInf = new TFile(TString("ZeeG_JESPES_DYJetsToLL_HT-2500toInf.root"));
  // TH1F* histo_Zll_100to200 = (TH1F*)f_Zll_100to200->Get(histname);
  // TH1F* histo_Zll_100to200_JESUp = (TH1F*)f_Zll_100to200->Get(histname_JESUp);
  // TH1F* histo_Zll_100to200_JESDown = (TH1F*)f_Zll_100to200->Get(histname_JESDown);
  // TH1F* histo_Zll_100to200_PESUp = (TH1F*)f_Zll_100to200->Get(histname_PESUp);
  // TH1F* histo_Zll_100to200_PESDown = (TH1F*)f_Zll_100to200->Get(histname_PESDown);
  // TH1F* histo_Zll_200to400 = (TH1F*)f_Zll_200to400->Get(histname);
  // TH1F* histo_Zll_200to400_JESUp = (TH1F*)f_Zll_200to400->Get(histname_JESUp);
  // TH1F* histo_Zll_200to400_JESDown = (TH1F*)f_Zll_200to400->Get(histname_JESDown);
  // TH1F* histo_Zll_200to400_PESUp = (TH1F*)f_Zll_200to400->Get(histname_PESUp);
  // TH1F* histo_Zll_200to400_PESDown = (TH1F*)f_Zll_200to400->Get(histname_PESDown);
  // TH1F* histo_Zll_400to600 = (TH1F*)f_Zll_400to600->Get(histname);
  // TH1F* histo_Zll_400to600_JESUp = (TH1F*)f_Zll_400to600->Get(histname_JESUp);
  // TH1F* histo_Zll_400to600_JESDown = (TH1F*)f_Zll_400to600->Get(histname_JESDown);
  // TH1F* histo_Zll_400to600_PESUp = (TH1F*)f_Zll_400to600->Get(histname_PESUp);
  // TH1F* histo_Zll_400to600_PESDown = (TH1F*)f_Zll_400to600->Get(histname_PESDown);
  // TH1F* histo_Zll_600to800 = (TH1F*)f_Zll_600to800->Get(histname);
  // TH1F* histo_Zll_600to800_JESUp = (TH1F*)f_Zll_600to800->Get(histname_JESUp);
  // TH1F* histo_Zll_600to800_JESDown = (TH1F*)f_Zll_600to800->Get(histname_JESDown);
  // TH1F* histo_Zll_600to800_PESUp = (TH1F*)f_Zll_600to800->Get(histname_PESUp);
  // TH1F* histo_Zll_600to800_PESDown = (TH1F*)f_Zll_600to800->Get(histname_PESDown);
  // TH1F* histo_Zll_800to1200 = (TH1F*)f_Zll_800to1200->Get(histname);
  // TH1F* histo_Zll_800to1200_JESUp = (TH1F*)f_Zll_800to1200->Get(histname_JESUp);
  // TH1F* histo_Zll_800to1200_JESDown = (TH1F*)f_Zll_800to1200->Get(histname_JESDown);
  // TH1F* histo_Zll_800to1200_PESUp = (TH1F*)f_Zll_800to1200->Get(histname_PESUp);
  // TH1F* histo_Zll_800to1200_PESDown = (TH1F*)f_Zll_800to1200->Get(histname_PESDown);
  // TH1F* histo_Zll_1200to2500 = (TH1F*)f_Zll_1200to2500->Get(histname);
  // TH1F* histo_Zll_1200to2500_JESUp = (TH1F*)f_Zll_1200to2500->Get(histname_JESUp);
  // TH1F* histo_Zll_1200to2500_JESDown = (TH1F*)f_Zll_1200to2500->Get(histname_JESDown);
  // TH1F* histo_Zll_1200to2500_PESUp = (TH1F*)f_Zll_1200to2500->Get(histname_PESUp);
  // TH1F* histo_Zll_1200to2500_PESDown = (TH1F*)f_Zll_1200to2500->Get(histname_PESDown);
  // TH1F* histo_Zll_2500toInf = (TH1F*)f_Zll_2500toInf->Get(histname);
  // TH1F* histo_Zll_2500toInf_JESUp = (TH1F*)f_Zll_2500toInf->Get(histname_JESUp);
  // TH1F* histo_Zll_2500toInf_JESDown = (TH1F*)f_Zll_2500toInf->Get(histname_JESDown);
  // TH1F* histo_Zll_2500toInf_PESUp = (TH1F*)f_Zll_2500toInf->Get(histname_PESUp);
  // TH1F* histo_Zll_2500toInf_PESDown = (TH1F*)f_Zll_2500toInf->Get(histname_PESDown);
  // histo_Zll_100to200->Scale(int_lumi*scale_factor*147.4/2489448.0);
  // histo_Zll_100to200_JESUp->Scale(int_lumi*scale_factor*147.4/2489448.0);
  // histo_Zll_100to200_JESDown->Scale(int_lumi*scale_factor*147.4/2489448.0);
  // histo_Zll_100to200_PESUp->Scale(int_lumi*scale_factor*147.4/2489448.0);
  // histo_Zll_100to200_PESDown->Scale(int_lumi*scale_factor*147.4/2489448.0);
  // histo_Zll_200to400->Scale(int_lumi*scale_factor*40.99/962154.0);
  // histo_Zll_200to400_JESUp->Scale(int_lumi*scale_factor*40.99/962154.0);
  // histo_Zll_200to400_JESDown->Scale(int_lumi*scale_factor*40.99/962154.0);
  // histo_Zll_200to400_PESUp->Scale(int_lumi*scale_factor*40.99/962154.0);
  // histo_Zll_200to400_PESDown->Scale(int_lumi*scale_factor*40.99/962154.0);
  // histo_Zll_400to600->Scale(int_lumi*scale_factor*5.678/1069553.0);
  // histo_Zll_400to600_JESUp->Scale(int_lumi*scale_factor*5.678/1069553.0);
  // histo_Zll_400to600_JESDown->Scale(int_lumi*scale_factor*5.678/1069553.0);
  // histo_Zll_400to600_PESUp->Scale(int_lumi*scale_factor*5.678/1069553.0);
  // histo_Zll_400to600_PESDown->Scale(int_lumi*scale_factor*5.678/1069553.0);
  // histo_Zll_600to800->Scale(int_lumi*scale_factor*1.367/8292160.0);
  // histo_Zll_600to800_JESUp->Scale(int_lumi*scale_factor*1.367/8292160.0);
  // histo_Zll_600to800_JESDown->Scale(int_lumi*scale_factor*1.367/8292160.0);
  // histo_Zll_600to800_PESUp->Scale(int_lumi*scale_factor*1.367/8292160.0);
  // histo_Zll_600to800_PESDown->Scale(int_lumi*scale_factor*1.367/8292160.0);
  // histo_Zll_800to1200->Scale(int_lumi*scale_factor*0.6304/2668311.0);
  // histo_Zll_800to1200_JESUp->Scale(int_lumi*scale_factor*0.6304/2668311.0);
  // histo_Zll_800to1200_JESDown->Scale(int_lumi*scale_factor*0.6304/2668311.0);
  // histo_Zll_800to1200_PESUp->Scale(int_lumi*scale_factor*0.6304/2668311.0);
  // histo_Zll_800to1200_PESDown->Scale(int_lumi*scale_factor*0.6304/2668311.0);
  // histo_Zll_1200to2500->Scale(int_lumi*scale_factor*0.1514/575630.0);
  // histo_Zll_1200to2500_JESUp->Scale(int_lumi*scale_factor*0.1514/575630.0);
  // histo_Zll_1200to2500_JESDown->Scale(int_lumi*scale_factor*0.1514/575630.0);
  // histo_Zll_1200to2500_PESUp->Scale(int_lumi*scale_factor*0.1514/575630.0);
  // histo_Zll_1200to2500_PESDown->Scale(int_lumi*scale_factor*0.1514/575630.0);
  // histo_Zll_2500toInf->Scale(int_lumi*scale_factor*0.003565/350990.0);
  // histo_Zll_2500toInf_JESUp->Scale(int_lumi*scale_factor*0.003565/350990.0);
  // histo_Zll_2500toInf_JESDown->Scale(int_lumi*scale_factor*0.003565/350990.0);
  // histo_Zll_2500toInf_PESUp->Scale(int_lumi*scale_factor*0.003565/350990.0);
  // histo_Zll_2500toInf_PESDown->Scale(int_lumi*scale_factor*0.003565/350990.0);
  // Float_t int_Zll_100to200 = histo_Zll_100to200->Integral()+histo_Zll_100to200->GetBinContent(0)+histo_Zll_100to200->GetBinContent(nBins+1);
  // Float_t int_Zll_200to400 = histo_Zll_200to400->Integral()+histo_Zll_200to400->GetBinContent(0)+histo_Zll_200to400->GetBinContent(nBins+1);
  // Float_t int_Zll_400to600 = histo_Zll_400to600->Integral()+histo_Zll_400to600->GetBinContent(0)+histo_Zll_400to600->GetBinContent(nBins+1);
  // Float_t int_Zll_600to800 = histo_Zll_600to800->Integral()+histo_Zll_600to800->GetBinContent(0)+histo_Zll_600to800->GetBinContent(nBins+1);
  // Float_t int_Zll_800to1200 = histo_Zll_800to1200->Integral()+histo_Zll_800to1200->GetBinContent(0)+histo_Zll_800to1200->GetBinContent(nBins+1);
  // Float_t int_Zll_1200to2500 = histo_Zll_1200to2500->Integral()+histo_Zll_1200to2500->GetBinContent(0)+histo_Zll_1200to2500->GetBinContent(nBins+1);
  // Float_t int_Zll_2500toInf = histo_Zll_2500toInf->Integral()+histo_Zll_2500toInf->GetBinContent(0)+histo_Zll_2500toInf->GetBinContent(nBins+1);
  // Float_t err_Zll_100to200 = 0.0;
  // Float_t err_Zll_200to400 = 0.0;
  // Float_t err_Zll_400to600 = 0.0;
  // Float_t err_Zll_600to800 = 0.0;
  // Float_t err_Zll_800to1200 = 0.0;
  // Float_t err_Zll_1200to2500 = 0.0;
  // Float_t err_Zll_2500toInf = 0.0;
  // if(int_Zll_100to200>0.0)
  //   err_Zll_100to200 = int_Zll_100to200*sqrt((147.4*int_lumi*scale_factor-int_Zll_100to200)/(2489448.0*int_Zll_100to200));
  // if(int_Zll_200to400>0.0)
  //   err_Zll_200to400 = int_Zll_200to400*sqrt((40.99*int_lumi*scale_factor-int_Zll_200to400)/(962154.0*int_Zll_200to400));
  // if(int_Zll_400to600>0.0)
  //   err_Zll_400to600 = int_Zll_400to600*sqrt((5.678*int_lumi*scale_factor-int_Zll_400to600)/(1069553.0*int_Zll_400to600));
  // if(int_Zll_600to800>0.0)
  //   err_Zll_600to800 = int_Zll_600to800*sqrt((1.367*int_lumi*scale_factor-int_Zll_600to800)/(8292160.0*int_Zll_600to800));
  // if(int_Zll_800to1200>0.0)
  //   err_Zll_800to1200 = int_Zll_800to1200*sqrt((0.6304*int_lumi*scale_factor-int_Zll_800to1200)/(2668311.0*int_Zll_800to1200));
  // if(int_Zll_1200to2500>0.0)
  //   err_Zll_1200to2500 = int_Zll_1200to2500*sqrt((0.1514*int_lumi*scale_factor-int_Zll_1200to2500)/(575630.0*int_Zll_1200to2500));
  // if(int_Zll_2500toInf>0.0)
  //   err_Zll_2500toInf = int_Zll_2500toInf*sqrt((0.003565*int_lumi*scale_factor-int_Zll_2500toInf)/(350990.0*int_Zll_2500toInf));
  // Float_t int_sum_Zll = int_Zll_100to200+int_Zll_200to400+int_Zll_400to600+int_Zll_600to800+int_Zll_800to1200+int_Zll_1200to2500+int_Zll_2500toInf;
  // Float_t err_sum_Zll = sqrt(err_Zll_100to200*err_Zll_100to200+err_Zll_200to400*err_Zll_200to400+err_Zll_400to600*err_Zll_400to600+err_Zll_600to800*err_Zll_600to800+err_Zll_800to1200*err_Zll_800to1200+err_Zll_1200to2500*err_Zll_1200to2500+err_Zll_2500toInf*err_Zll_2500toInf);
  // for(int i; i <= nBins; i++){
  //   //Skipping jes and pes for now
  //   double int_bin_100to200 = histo_Zll_100to200->GetBinContent(i);
  //   double int_bin_200to400 = histo_Zll_200to400->GetBinContent(i);
  //   double int_bin_400to600 = histo_Zll_400to600->GetBinContent(i);
  //   double int_bin_600to800 = histo_Zll_600to800->GetBinContent(i);
  //   double int_bin_800to1200 = histo_Zll_800to1200->GetBinContent(i);
  //   double int_bin_1200to2500 = histo_Zll_1200to2500->GetBinContent(i);
  //   double int_bin_2500toInf = histo_Zll_2500toInf->GetBinContent(i);
  //   double err_bin_100to200 = 0.0;
  //   double err_bin_200to400 = 0.0;
  //   double err_bin_400to600 = 0.0;
  //   double err_bin_600to800 = 0.0;
  //   double err_bin_800to1200 = 0.0;
  //   double err_bin_1200to2500 = 0.0;
  //   double err_bin_2500toInf = 0.0;
  //   if(int_bin_100to200>0.0)
  //     err_bin_100to200 = int_bin_100to200*sqrt((147.4*int_lumi*scale_factor-int_bin_100to200)/(2489448.0*int_bin_100to200));
  //   if(int_bin_200to400>0.0)
  //     err_bin_200to400 = int_bin_200to400*sqrt((40.99*int_lumi*scale_factor-int_bin_200to400)/(962154.0*int_bin_200to400));
  //   if(int_bin_400to600>0.0)
  //     err_bin_400to600 = int_bin_400to600*sqrt((5.678*int_lumi*scale_factor-int_bin_400to600)/(1069553.0*int_bin_400to600));
  //   if(int_bin_600to800>0.0)
  //     err_bin_600to800 = int_bin_600to800*sqrt((1.367*int_lumi*scale_factor-int_bin_600to800)/(8292160.0*int_bin_600to800));
  //   if(int_bin_800to1200>0.0)
  //     err_bin_800to1200 = int_bin_800to1200*sqrt((0.6304*int_lumi*scale_factor-int_bin_800to1200)/(2668311.0*int_bin_800to1200));
  //   if(int_bin_1200to2500>0.0)
  //     err_bin_1200to2500 = int_bin_1200to2500*sqrt((0.1514*int_lumi*scale_factor-int_bin_1200to2500)/(575630.0*int_bin_1200to2500));
  //   if(int_bin_2500toInf>0.0)
  //     err_bin_2500toInf = int_bin_2500toInf*sqrt((0.003565*int_lumi*scale_factor-int_bin_2500toInf)/(350990.0*int_bin_2500toInf));
  //   histo_Zll_100to200->SetBinError(i,err_bin_100to200);
  //   histo_Zll_200to400->SetBinError(i,err_bin_200to400);
  //   histo_Zll_400to600->SetBinError(i,err_bin_400to600);
  //   histo_Zll_600to800->SetBinError(i,err_bin_600to800);
  //   histo_Zll_800to1200->SetBinError(i,err_bin_800to1200);
  //   histo_Zll_1200to2500->SetBinError(i,err_bin_1200to2500);
  //   histo_Zll_2500toInf->SetBinError(i,err_bin_2500toInf);
  // }
  // TH1F* histo_Zll_100toInf = (TH1F*)histo_Zll_100to200->Clone("");
  // histo_Zll_100toInf->Add(histo_Zll_200to400);
  // histo_Zll_100toInf->Add(histo_Zll_400to600);
  // histo_Zll_100toInf->Add(histo_Zll_600to800);
  // histo_Zll_100toInf->Add(histo_Zll_800to1200);
  // histo_Zll_100toInf->Add(histo_Zll_1200to2500);
  // histo_Zll_100toInf->Add(histo_Zll_2500toInf);
  // histo_Zll_100toInf->SetFillColor(kPink+1);
  // histo_vector.push_back(histo_Zll_100toInf);
  
  // Print bin contents
  cout<<endl;
  if (histname=="Photon_Et_range_4"){
    vector<float> total_background_binned;
    total_background_binned.clear();
    for(int i = 1; i <= nBins; i ++){
      total_background_binned.push_back(0.0);
    }
    
    cout<<"$E_{T}^{\\gamma}$ &        [ 175,  200] &        [ 200,  250] &        [ 250,  300] &        [ 300,  400] &        [ 400,  600] &        [ 600, 1000] \\\\"<<endl;
    
    cout<<"\\hline"<<endl;
    
    cout<<"        jetfake ";
    for(int i = 1; i <= nBins; i++){
      cout<<"& $ "<<boost::format("%.2f")%histo_jetfake->GetBinContent(i)<<" $ ";
      total_background_binned[i-1] += histo_jetfake->GetBinContent(i);
    }
    cout<<"\\\\"<<endl;
    // cout<<"         spikes ";
    // for(int i = 1; i <= nBins; i++){
    //   cout<<"& $ "<<boost::format("%.2f")%histo_spikes->GetBinContent(i)<<" $ ";
    //   total_background_binned[i-1] += histo_spikes->GetBinContent(i);
    // }
    // cout<<"\\\\"<<endl;
    // cout<<"        elefake ";
    // for(int i = 1; i <= nBins; i++){
    //   cout<<"& $ "<<boost::format("%.2f")%histo_elefake->GetBinContent(i)<<" $ ";
    //   total_background_binned[i-1] += histo_elefake->GetBinContent(i);
    // }
    // cout<<"\\\\"<<endl;
    // cout<<"          bhalo ";
    // for(int i = 1; i <= nBins; i++){
    //   cout<<"& $ "<<boost::format("%.2f")%histo_bhalo->GetBinContent(i)<<" $ ";
    //   total_background_binned[i-1] += histo_bhalo->GetBinContent(i);
    // }
    // cout<<"\\\\"<<endl;
    // cout<<"         ZNuNuG ";
    // for(int i = 1; i <= nBins; i++){
    //   cout<<"& $ "<<boost::format("%.2f")%histo_ZNuNuG->GetBinContent(i)<<" $ ";
    //   total_background_binned[i-1] += histo_ZNuNuG->GetBinContent(i);
    // }
    // cout<<"\\\\"<<endl;
    // cout<<"             WG ";
    // for(int i = 1; i <= nBins; i++){
    //   cout<<"& $ "<<boost::format("%.2f")%histo_WG->GetBinContent(i)<<" $ ";
    //   total_background_binned[i-1] += histo_WG->GetBinContent(i);
    // }
    // cout<<"\\\\"<<endl;
    // cout<<"          GJets ";
    // for(int i = 1; i <= nBins; i++){
    //   cout<<"& $ "<<boost::format("%.2f")%histo_GJets_40toInf->GetBinContent(i)<<" $ ";
    //   total_background_binned[i-1] += histo_GJets_40toInf->GetBinContent(i);
    // }
    // cout<<"\\\\"<<endl;
    cout<<"           ZllG ";
    for(int i = 1; i <= nBins; i++){
      cout<<"& $ "<<boost::format("%.2f")%histo_ZllG_combined->GetBinContent(i)<<" $ ";
      total_background_binned[i-1] += histo_ZllG_combined->GetBinContent(i);
    }
    cout<<"\\\\"<<endl;
    cout<<"            TTG ";
    for(int i = 1; i <= nBins; i++){
      if(histo_TTG->GetBinContent(i) < 0.0){
        int_TTG -= histo_TTG->GetBinContent(i);
        histo_TTG->SetBinContent(i, 0.0);
      }
      cout<<"& $ "<<boost::format("%.2f")%histo_TTG->GetBinContent(i)<<" $ ";
      total_background_binned[i-1] += histo_TTG->GetBinContent(i);
    }
    cout<<"\\\\"<<endl;
    // cout<<"             TG ";
    // for(int i = 1; i <= nBins; i++){
    //   if(histo_TG->GetBinContent(i) < 0.0){
    //     int_TG -= histo_TG->GetBinContent(i);
    //     histo_TG->SetBinContent(i, 0.0);
    //   }
    //   cout<<"& $ "<<boost::format("%.2f")%histo_TG->GetBinContent(i)<<" $ ";
    //   total_background_binned[i-1] += histo_TG->GetBinContent(i);
    // }
    // cout<<"\\\\"<<endl;
    // cout<<"        diphton ";
    // for(int i = 1; i <= nBins; i++){
    //   if(histo_diphoton->GetBinContent(i) < 0.0){
    //     int_diphoton -= histo_diphoton->GetBinContent(i);
    //     histo_diphoton->SetBinContent(i, 0.0);
    //   }
    //   cout<<"& $ "<<boost::format("%.2f")%histo_diphoton->GetBinContent(i)<<" $ ";
    //   total_background_binned[i-1] += histo_diphoton->GetBinContent(i);
    // }
    // cout<<"\\\\"<<endl;
    cout<<"             WZ ";
    for(int i = 1; i <= nBins; i++){
      cout<<"& $ "<<boost::format("%.2f")%histo_WZ->GetBinContent(i)<<" $ ";
      total_background_binned[i-1] += histo_WZ->GetBinContent(i);
    }
    cout<<"\\\\"<<endl;
    // cout<<"             ZZ ";
    // for(int i = 1; i <= nBins; i++){
    //   cout<<"& $ "<<boost::format("%.2f")%histo_ZZ->GetBinContent(i)<<" $ ";
    //   total_background_binned[i-1] += histo_ZZ->GetBinContent(i);
    // }
    // cout<<"\\\\"<<endl;
    // cout<<"          WMuNu ";
    // for(int i = 1; i <= nBins; i++){
    //   cout<<"& $ "<<boost::format("%.2f")%histo_WMuNu->GetBinContent(i)<<" $ ";
    //   total_background_binned[i-1] += histo_WMuNu->GetBinContent(i);
    // }
    // cout<<"\\\\"<<endl;
    // cout<<"         WTauNu ";
    // for(int i = 1; i <= nBins; i++){
    //   cout<<"& $ "<<boost::format("%.2f")%histo_WTauNu->GetBinContent(i)<<" $ ";
    //   total_background_binned[i-1] += histo_WTauNu->GetBinContent(i);
    // }
    // cout<<"\\\\"<<endl;
    // cout<<"             WW ";
    // for(int i = 1; i <= nBins; i++){
    //   cout<<"& $ "<<boost::format("%.2f")%histo_WW->GetBinContent(i)<<" $ ";
    //   total_background_binned[i-1] += histo_WW->GetBinContent(i);
    // }
    // cout<<"\\\\"<<endl;
    
    cout<<"\\hline"<<endl;
    
    cout<<"          total ";
    for(int i = 1; i <= nBins; i++){
      cout<<"& $ "<<boost::format("%.2f")%total_background_binned[i-1]<<" $ ";
    }
    cout<<"\\\\"<<endl;
    
    cout<<"\\hline"<<endl;
    
    cout<<"           data ";
    for(int i = 1; i <= nBins; i++){
      cout<<"& $ "<<histo_data->GetBinContent(i)<<" $ ";
    }
    cout<<"\\\\"<<endl;
    
    cout<<"\\hline"<<endl;
  }
  
  if (histname == "Photon_Et_range_4" || histname == "h_photonic_recoil_4" || histname == "h_phoRecoilMt_4")
  {
    TFile* f_ZeeG_histos;
    if(histname == "Photon_Et_range_4")
      f_ZeeG_histos = new TFile("ZeeG_histos_Pt.root","RECREATE");
    else if(histname == "h_photonic_recoil_4")
      f_ZeeG_histos = new TFile("ZeeG_histos_MET.root","RECREATE");
    else if(histname == "h_phoRecoilMt_4")
      f_ZeeG_histos = new TFile("ZeeG_histos_Mt.root","RECREATE");
    f_ZeeG_histos->cd();
    histo_data->Write();
    histo_ZllG_combined->Write();
    histo_ZllG_JESUp_combined->Write();
    histo_ZllG_JESDown_combined->Write();
    histo_ZllG_PESUp_combined->Write();
    histo_ZllG_PESDown_combined->Write();
    histo_ZllG_straightUp_combined->Write();
    histo_ZllG_straightDown_combined->Write();
    histo_ZllG_twistedUp_combined->Write();
    histo_ZllG_twistedDown_combined->Write();
    histo_ZllG_gammaUp_combined->Write();
    histo_ZllG_gammaDown_combined->Write();
    histo_ZllG_qcdscale_combined->Write();
    histo_TTG->Write();
    histo_TTG_JESUp->Write();
    histo_TTG_JESDown->Write();
    histo_TTG_PESUp->Write();
    histo_TTG_PESDown->Write();
    histo_WZ->Write();
    histo_WZ_JESUp->Write();
    histo_WZ_JESDown->Write();
    histo_WZ_PESUp->Write();
    histo_WZ_PESDown->Write();
    // histo_ZGGToLLGG->Write();
    // histo_ZGGToLLGG_JESUp->Write();
    // histo_ZGGToLLGG_JESDown->Write();
    // histo_ZGGToLLGG_PESUp->Write();
    // histo_ZGGToLLGG_PESDown->Write();
    // histo_ZGGToLLGG_phoSFUp->Write();
    // histo_ZGGToLLGG_phoSFDown->Write();
    histo_jetfake->Write();
    histo_jetfake_errUp->Write();
    histo_jetfake_errDown->Write();
    f_ZeeG_histos->Close();
  }
  
  if (histname == "Photon_Et_range_4" || histname == "h_photonic_recoil_4"){
    for(int i = 1; i <= nBins; i++){
      double binWidth = histo_data->GetBinWidth(i);
      histo_data->SetBinContent(i,histo_data->GetBinContent(i)/binWidth);
      histo_data->SetBinError(i,histo_data->GetBinError(i)/binWidth);
      histo_ZllG_combined->SetBinContent(i,histo_ZllG_combined->GetBinContent(i)/binWidth);
      histo_ZllG_combined->SetBinError(i,histo_ZllG_combined->GetBinError(i)/binWidth);
      histo_TTG->SetBinContent(i,histo_TTG->GetBinContent(i)/binWidth);
      histo_TTG->SetBinError(i,histo_TTG->GetBinError(i)/binWidth);
      // histo_WWG->SetBinContent(i,histo_WWG->GetBinContent(i)/binWidth);
      // histo_WWG->SetBinError(i,histo_WWG->GetBinError(i)/binWidth);
      histo_WZ->SetBinContent(i,histo_WZ->GetBinContent(i)/binWidth);
      histo_WZ->SetBinError(i,histo_WZ->GetBinError(i)/binWidth);
      // histo_ZGGToLLGG->SetBinContent(i,histo_ZGGToLLGG->GetBinContent(i)/binWidth);
      // histo_ZGGToLLGG->SetBinError(i,histo_ZGGToLLGG->GetBinError(i)/binWidth);
      // histo_Zll_100toInf->SetBinContent(i,histo_Zll_100toInf->GetBinContent(i)/binWidth);
      // histo_Zll_100toInf->SetBinError(i,histo_Zll_100toInf->GetBinError(i)/binWidth);
      histo_jetfake->SetBinContent(i,histo_jetfake->GetBinContent(i)/binWidth);
      histo_jetfake->SetBinError(i,histo_jetfake->GetBinError(i)/binWidth);
    }
  }
  TPad *pad1 = new TPad("pad1","pad1",0.01,0.26,0.99,0.99);
  pad1->Draw(); pad1->cd();
  pad1->SetFillColor(0); pad1->SetFrameBorderMode(0); pad1->SetBorderMode(0);
  pad1->SetBottomMargin(0.);
  
  TH1F *histo_allbackgrounds = (TH1F*)histo_ZllG_combined->Clone("histo_allbackgrounds");
  histo_allbackgrounds->Add(histo_TTG);
  // histo_allbackgrounds->Add(histo_WWG);
  histo_allbackgrounds->Add(histo_WZ);
  // histo_allbackgrounds->Add(histo_Zll_100toInf);
  histo_allbackgrounds->Add(histo_jetfake);
    
  for(int i = 1; i <= nBins; i++){
    double background = histo_allbackgrounds->GetBinContent(i);
    // Add statistical errors
    double sum_binerrors_squared = 0.0;
    sum_binerrors_squared += pow(histo_ZllG_combined->GetBinError(i),2);
    sum_binerrors_squared += pow(histo_TTG->GetBinError(i),2);
    // sum_binerrors_squared += pow(histo_WWG->GetBinError(i),2);
    sum_binerrors_squared += pow(histo_WZ->GetBinError(i),2);
    // sum_binerrors_squared += pow(histo_ZGGToLLGG->GetBinError(i),2);
    // sum_binerrors_squared += pow(histo_Zll_100toInf->GetBinError(i),2);
    sum_binerrors_squared += pow(histo_jetfake->GetBinError(i),2);
    double binerror = sqrt(sum_binerrors_squared); // Include just the statistical error
    double jeserr = (fabs(jesup_shift[i-1])+fabs(jesdown_shift[i-1]))/2.0;
    double peserr = (fabs(pesup_shift[i-1])+fabs(pesdown_shift[i-1]))/2.0;
    double straighterr_ZllG = (fabs(straightup_shift_ZllG[i-1])+fabs(straightdown_shift_ZllG[i-1]))/2.0;
    double twistederr_ZllG = (fabs(twistedup_shift_ZllG[i-1])+fabs(twisteddown_shift_ZllG[i-1]))/2.0;
    double gammaerr_ZllG = (fabs(gammaup_shift_ZllG[i-1])+fabs(gammadown_shift_ZllG[i-1]))/2.0;
    double renerr = (fabs(renup_shift[i-1])+fabs(rendown_shift[i-1]))/2.0;
    double facerr = (fabs(facup_shift[i-1])+fabs(facdown_shift[i-1]))/2.0;
    double pdferr = (fabs(pdfup_shift[i-1])+fabs(pdfdown_shift[i-1]))/2.0;
    double qcdscaleerr_ZllG = fabs(qcdscale_shift_ZllG[i-1]);
    double jetfakeerr = (fabs(syst_shiftUp_jetfake[i-1])+fabs(syst_shiftDown_jetfake[i-1]))/2.0;
    if (histname == "Photon_Et_range_4" || histname == "h_photonic_recoil_4"){
      double binWidth = histo_data->GetBinWidth(i);
      jeserr /= binWidth;
      peserr /= binWidth;
      straighterr_ZllG /= binWidth;
      twistederr_ZllG /= binWidth;
      gammaerr_ZllG /= binWidth;
      renerr /= binWidth;
      facerr /= binWidth;
      pdferr /= binWidth;
      qcdscaleerr_ZllG /= binWidth;
      jetfakeerr /= binWidth;
    }
    binerror = sqrt(sum_binerrors_squared+pow(background*photon_scale_factor_unc,2)+pow(background*electron_scale_factor_unc,2)+pow(jeserr,2)+pow(peserr,2)+pow(straighterr_ZllG,2)+pow(twistederr_ZllG,2)+pow(gammaerr_ZllG,2)+pow(renerr,2)+pow(facerr,2)+pow(pdferr,2)+pow(qcdscaleerr_ZllG,2)+pow(jetfakeerr,2));
    histo_allbackgrounds->SetBinError(i,binerror);
  }
  histo_allbackgrounds->SetFillColorAlpha(kGray+1,0.6);
  histo_vector.push_back(histo_allbackgrounds);
  
  if (histname == "Photon_Et_range_4"){
    for(int i = 1; i <= nBins; i++){
      float data = histo_data->GetBinContent(i);
      float binWidth = histo_data->GetBinWidth(i);
      cout<<"data, bin "<<i<<": "<<(data*binWidth)<<endl;
    }
    for(int i = 1; i <= nBins; i++){
      float background = histo_allbackgrounds->GetBinContent(i);
      float background_err = histo_allbackgrounds->GetBinError(i);
      float binWidth = histo_allbackgrounds->GetBinWidth(i);
      cout<<"background, bin "<<i<<": "<<(background*binWidth)<<" +/- "<<(background_err*binWidth)<<endl;
    }
  }

  TH1F *histo_allbackgrounds_outline = (TH1F*)histo_allbackgrounds->Clone("histo_allbackgrounds_outline");
  histo_allbackgrounds_outline->SetFillColorAlpha(kWhite,0.0);
  histo_allbackgrounds_outline->SetLineWidth(1);
  histo_vector.push_back(histo_allbackgrounds_outline);
    
  THStack *stackHisto = new THStack("stackHisto","Title");
  stackHisto->Add(histo_WZ);
  stackHisto->Add(histo_TTG);
  stackHisto->Add(histo_jetfake);
  stackHisto->Add(histo_ZllG_combined);
  stackHisto->SetTitle("");
  
  for(int i = 0; i < int(histo_vector.size()); i++){
    histo_vector[i]->SetStats(0);
    histo_vector[i]->SetTitle("");
    histo_vector[i]->SetLineColor(kBlack);
    histo_vector[i]->GetXaxis()->SetTitle(xaxis_title);
    histo_vector[i]->GetXaxis()->SetLabelFont(42);
    histo_vector[i]->GetXaxis()->SetLabelSize(0.06);
    histo_vector[i]->GetXaxis()->SetTitleFont(42);
    histo_vector[i]->GetXaxis()->SetTitleSize(0.06);
    histo_vector[i]->GetYaxis()->SetTitle("Events / bin");
    if (histname == "Photon_Et_range_4" || histname == "h_photonic_recoil_4")
      histo_vector[i]->GetYaxis()->SetTitle("Events / GeV");
    histo_vector[i]->GetYaxis()->SetLabelFont(42);
    histo_vector[i]->GetYaxis()->SetLabelSize(0.06);
    histo_vector[i]->GetYaxis()->SetTitleFont(42);
    histo_vector[i]->GetYaxis()->SetTitleSize(0.06);
    histo_vector[i]->GetYaxis()->SetTitleOffset(0.9);
  }
  
  //Accommodate both the data and background plots
  double ymax_data = 0.0;
  double ymax_background = 0.0;
  for(int i = 1; i <= nBins; i++){
    double y_data = histo_data->GetBinContent(i);
    double y_error_data = histo_data->GetBinError(i);
    double y_high_data = y_data+y_error_data;
    if(y_high_data > ymax_data)
      ymax_data = y_high_data;
    double y_background = histo_allbackgrounds->GetBinContent(i);
    double y_error_background = histo_allbackgrounds->GetBinError(i);
    double y_high_background = y_background+y_error_background;
    if(y_high_background > ymax_background)
      ymax_background = y_high_background;
  }
  
  double ymin = 0.0003;
  double ymax = 1.3*ymax_data;
  if(ymax_background > ymax_data)
    ymax = 1.3*ymax_background;
  if (histname == "Photon_Et_range_4" || histname == "h_photonic_recoil_4"){
    pad1->SetLogy();
    // ymax *= 5;
    ymax = 5;
  }
  // else if (histname == "nJet_4"){
  //   ymax = 2.0;
  // }
  histo_data->GetYaxis()->SetRangeUser(ymin,ymax);
  if(histname == "nJet_4")
    histo_data->GetXaxis()->SetRangeUser(0,10);
  histo_data->Draw();
  stackHisto->Draw("HIST SAME");
  histo_allbackgrounds->Draw("E2 SAME");
  histo_allbackgrounds_outline->Draw("HIST SAME");
  histo_data->SetLineColor(kBlack);
  histo_data->SetMarkerColor(kBlack);
  histo_data->Draw("E0 P0 SAME");
  gPad->RedrawAxis();
  
  //Central location of leg defined to be location of leg in phoPt plot
  TLegend* leg = new TLegend(0.5+leg_xoffset,0.58075+leg_yoffset,0.885387+leg_xoffset,0.862969+leg_yoffset,"");
  leg->AddEntry(histo_data,"Data");
  leg->AddEntry(histo_ZllG_combined,"Z(ll)#gamma","F");
  leg->AddEntry(histo_jetfake,"jet#rightarrow#gamma MisID","F");
  leg->AddEntry(histo_TTG,"tt#gamma","F");
  leg->AddEntry(histo_WZ,"VV#gamma","F");
  leg->SetFillColor(kWhite);
  leg->SetShadowColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.045);
  leg->Draw();
  
  float lumiTextSize = 0.6;
  float lumiTextOffset = 0.2;
  float cmsTextSize = 0.75;
  TLatex *texS = new TLatex(0.60023,0.917173,"35.9 fb^{-1} (13 TeV)");
  texS->SetNDC();
  texS->SetTextFont(42);
  texS->SetTextSize(lumiTextSize*t_m);
  texS->Draw();
  TLatex *texS1 = new TLatex(0.13592,0.817173,"#bf{CMS} #it{Preliminary}");
  texS1->SetNDC();
  texS1->SetTextFont(42);
  texS1->SetTextSize(cmsTextSize*t_m);
  texS1->Draw();
  
  c->cd();
  TPad *pad2 = new TPad("pad2","pad2",0.01,0.01,0.99,0.26);
  pad2->Draw(); pad2->cd();
  pad2->SetFillColor(0); pad2->SetFrameBorderMode(0); pad2->SetBorderMode(0);
  pad2->SetTopMargin(0.0);
  pad2->SetBottomMargin(0.35);
  
  double max_ratio = 3.5;
  
  TH1F* Ratio = (TH1F*)histo_data->Clone("Ratio");
  TH1F* Ratio_background = (TH1F*)histo_allbackgrounds->Clone("Ratio_background");
  for(int i = 1; i <= nBins; i++){
    double y_data = histo_data->GetBinContent(i);
    double y_error_data = histo_data->GetBinError(i);
    double y_background = histo_allbackgrounds->GetBinContent(i);
    double y_error_background = histo_allbackgrounds->GetBinError(i);
    double Ratiocontent = 0.0;
    double Ratioerror = max_ratio;
    double Ratioerror_background = max_ratio;
    if(y_background > 0.){
      Ratiocontent = y_data/y_background;
      Ratioerror_background = y_error_background/y_background;
      if(y_error_data > 0.)
        Ratioerror = y_error_data/y_background;
    }
    else if(y_data > 0.){
      Ratiocontent = 3.*max_ratio;
    }
    Ratio->SetBinContent(i,Ratiocontent);
    Ratio->SetBinError(i,Ratioerror);
    Ratio_background->SetBinContent(i,1);
    Ratio_background->SetBinError(i,Ratioerror_background);
  }
  
  Ratio_background->GetYaxis()->SetRangeUser(0.0,max_ratio-0.01);
  Ratio_background->GetYaxis()->SetTitle("Data/SM");
  Ratio_background->GetYaxis()->CenterTitle();
  Ratio_background->GetYaxis()->SetLabelSize(0.14);
  Ratio_background->GetYaxis()->SetTitleSize(0.15);
  Ratio_background->GetYaxis()->SetLabelFont(42);
  Ratio_background->GetYaxis()->SetTitleFont(42);
  Ratio_background->GetYaxis()->SetTitleOffset(0.30);
  Ratio_background->GetYaxis()->SetNdivisions(305);
  Ratio_background->GetXaxis()->SetTitle(xaxis_title);
  Ratio_background->GetXaxis()->SetLabelSize(0.16);
  Ratio_background->GetXaxis()->SetTitleSize(0.18);
  Ratio_background->GetXaxis()->SetLabelFont(42);
  Ratio_background->GetXaxis()->SetTitleFont(42);
  Ratio_background->GetXaxis()->SetTitleOffset(0.9);
  Ratio_background->GetXaxis()->SetTickLength(0.05);
  Ratio_background->SetStats(0);
  Ratio->SetMarkerStyle(0);
  double xmin = histo_data->GetXaxis()->GetBinLowEdge(1);
  double xmax = histo_data->GetXaxis()->GetBinUpEdge(nBins);
  if(histname == "nJet_4"){
    Ratio_background->GetXaxis()->SetRangeUser(0,10);
    xmax = 10.0;
  }
  TLine* line = new TLine(xmin,1.,xmax,1.);
  line->SetLineStyle(2);
  line->SetLineColor(kBlack);
  gStyle->SetLineStyleString(11,"3 12");
  TLine* line0 = new TLine(xmin,0.5,xmax,0.5);
  line0->SetLineStyle(11);
  line0->SetLineColor(kBlack);
  Ratio_background->Draw("E2");
  line->Draw("SAME");
  line0->Draw("SAME");
  for(int i = 1; i <= (2*max_ratio-3); i++){
    double y_coord = 1.0 + 0.5*i;
    TLine* line_i = new TLine(xmin,y_coord,xmax,y_coord);
    line_i->SetLineStyle(11);
    line_i->SetLineColor(kBlack);
    line_i->Draw("SAME");
  }
  Ratio->Draw("E0 P0 SAME");

  double background_unc = sqrt(background_unc_sumsquares);
  
  if(histname == "Photon_Et_range_4"){
    cout<<"ZeeG region"<<endl;
    cout<<"------------------------------------"<<endl;
    cout<<"Z(ll)+Gamma: "<<int_ZllG<<" +- "<<err_ZllG<<endl;
    cout<<"tt+Gamma: "<<int_TTG<<" +- "<<err_TTG<<endl;
    // cout<<"WWG: "<<int_WWG<<" +- "<<err_WWG<<endl;
    cout<<"WZ: "<<int_WZ<<" +- "<<err_WZ<<endl;
    // cout<<"ZGGToLLGG: "<<int_ZGGToLLGG<<" +- "<<err_ZGGToLLGG<<endl;
    // cout<<"Zll: "<<int_sum_Zll<<" +- "<<err_sum_Zll<<endl;
    cout<<"Jet faking photon: "<<int_jetfake<<" +- "<<err_jetfake<<endl;
    cout<<"Total background: "<<total_background<<" +- "<<background_unc<<endl;
    cout<<"------------------------------------"<<endl;
    cout<<"Data: "<<int_data<<endl;
    cout<<"------------------------------------"<<endl;
  }

  c->SaveAs(TString("zeeg_"+plotname+".png"));
  c->SaveAs(TString("zeeg_"+plotname+".pdf"));
  delete(c);
}

void zeeg_plotter()
{
  std::vector<string> histnames;
  histnames.clear();
  std::vector<Double_t> leg_xoffsets;
  leg_xoffsets.clear();
  std::vector<Double_t> leg_yoffsets;
  leg_yoffsets.clear();
  std::vector<TString> xaxis_titles;
  xaxis_titles.clear();
  std::vector<TString> plotnames;
  plotnames.clear();

//  histnames.push_back(TString("_4"));
//  leg_xoffsets.push_back(0.);
//  leg_yoffsets.push_back(0.);
//  xaxis_titles.push_back(TString(""));
//  plotnames.push_back(TString(""));

  histnames.push_back("Photon_Et_range");
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Photon #it{p}_{T} [GeV]"));
  plotnames.push_back(TString("phoPt"));

  histnames.push_back("Photon_SCeta");
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Photon #eta"));
  plotnames.push_back(TString("phoEta"));

  histnames.push_back("Photon_SCphi");
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Photon #phi"));
  plotnames.push_back(TString("phoPhi"));

  histnames.push_back("pfMET");
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("pfMET [GeV]"));
  plotnames.push_back(TString("pfMET"));

  histnames.push_back("nJet");
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Number of Jets"));
  plotnames.push_back(TString("nJet"));

  histnames.push_back("h_photonic_recoil");
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Photonic Recoil [GeV]"));
  plotnames.push_back(TString("recoil"));

  histnames.push_back("h_dPhi_phoRecoil");
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("#Delta#phi(Photon,Recoil)"));
  plotnames.push_back(TString("dPhiPhoRecoil"));

  histnames.push_back("h_leadingLeptonPt");
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Leading Electron #it{p}_{T} [GeV]"));
  plotnames.push_back(TString("elePtLeading"));

  histnames.push_back("h_leadingLeptonEta");
  leg_xoffsets.push_back(0.15);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Leading Electron #eta"));
  plotnames.push_back(TString("eleEtaLeading"));

  histnames.push_back("h_leadingLeptonPhi");
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Leading Electron #phi"));
  plotnames.push_back(TString("elePhiLeading"));

  histnames.push_back("h_subleadingLeptonPt");
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Subleading Electron #it{p}_{T} [GeV]"));
  plotnames.push_back(TString("elePtSubleading"));

  histnames.push_back("h_subleadingLeptonEta");
  leg_xoffsets.push_back(0.15);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Subleading Electron #eta"));
  plotnames.push_back(TString("eleEtaSubleading"));

  histnames.push_back("h_subleadingLeptonPhi");
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Subleading Electron #phi"));
  plotnames.push_back(TString("elePhiSubleading"));
  
  histnames.push_back("h_dileptonPt");
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Dielectron #it{p}_{T} [GeV]"));
  plotnames.push_back(TString("dielePt"));
  
  histnames.push_back("h_dileptonM");
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Dielectron Invariant Mass [GeV]"));
  plotnames.push_back(TString("dieleM"));

  histnames.push_back("h_phoPT_over_photonicRecoil");
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Photon #it{p}_{T} / Recoil"));
  plotnames.push_back(TString("phoPtOverRecoil"));

  histnames.push_back("h_dileptonPt_over_pfMET");
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Dielectron #it{p}_{T} / pfMET"));
  plotnames.push_back(TString("dielePtOverpfMET"));

  histnames.push_back("h_min_dphijetmet");
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Min. #Delta#phi(jets,MET)"));
  plotnames.push_back(TString("dPhiJetsMET"));

  histnames.push_back("h_min_dphijetrecoil");
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Min. #Delta#phi(jets,Recoil)"));
  plotnames.push_back(TString("dphiJetsRecoil"));

  histnames.push_back("h_dPhi_leptons");
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("#Delta#phi(Leading Ele,Subleading Ele)"));
  plotnames.push_back(TString("dPhiLeadingEleSubleadingEle"));

  histnames.push_back("h_phoRecoilMt");
  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Photon-Recoil #it{M}_{T} [GeV]"));
  plotnames.push_back(TString("phoRecoilmT"));
  
  for(int i = 0; i < histnames.size(); i++){
    plot(histnames[i],leg_xoffsets[i],leg_yoffsets[i],xaxis_titles[i],plotnames[i]);
  }
}