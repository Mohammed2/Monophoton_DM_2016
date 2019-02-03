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

void plot(Double_t leg_xoffset, Double_t leg_yoffset, TString xaxis_title, TString plotname)
{
  std::vector<TH1F*> histo_vector;
  histo_vector.clear();
  
  double total_background = 0.0;
  
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
  
  TFile *f_postfit_shapes = new TFile("fitDiagnostics.root");
  TGraphAsymmErrors* tgae_data = (TGraphAsymmErrors*)((TGraphAsymmErrors*)f_postfit_shapes->Get("shapes_prefit/CR_DE/data"))->Clone("tgae_data");
  Float_t PtBins[7]={175.,200.,250., 300., 400., 600., 1000.0};
  TH1F* histo_data = new TH1F("data","data",6,PtBins);
  const int nBins = histo_data->GetXaxis()->GetNbins();
  for(int i = 1; i <= nBins; i++){
    Double_t x = 0.0;
    Double_t y = 0.0;
    tgae_data->GetPoint(i-1, x, y);
    histo_data->SetBinContent(i, y*histo_data->GetBinWidth(i));
    histo_data->SetBinError(i, tgae_data->GetErrorY(i-1)*histo_data->GetBinWidth(i));
    cout<<"data_bin "<<i<<": y="<<y<<endl;
  }
  cout<<"histo_data_integral: "<<(histo_data->Integral())<<endl;
  Float_t int_data = histo_data->Integral();
  histo_data->SetLineWidth(3);
  histo_data->SetLineColor(kWhite);
  histo_data->SetMarkerStyle(kFullSquare);
  histo_data->SetMarkerColor(kWhite);
  histo_vector.push_back(histo_data);
    
  TH1F* histo_jetfake = (TH1F*)((TH1F*)f_postfit_shapes->Get("shapes_prefit/CR_DE/QCD"))->Clone("histo_jetfake");
  // histo bin contents have been divided by bin width, so restore true normalization by multiplying by bin width
  for(int i = 1; i <= nBins; i++){
    histo_jetfake->SetBinContent(i, histo_jetfake->GetBinContent(i)*histo_jetfake->GetBinWidth(i));
    histo_jetfake->SetBinError(i, histo_jetfake->GetBinError(i)*histo_jetfake->GetBinWidth(i));
    cout<<"jetfake_bin_"<<i<<": error="<<histo_jetfake->GetBinError(i)<<endl;
  }
  Float_t int_jetfake = histo_jetfake->Integral();
  total_background += int_jetfake;
  histo_jetfake->SetFillColor(kBlue-4);
  histo_vector.push_back(histo_jetfake);
  
  TH1F* histo_ZllG_combined = (TH1F*)((TH1F*)f_postfit_shapes->Get("shapes_prefit/CR_DE/ZllG"))->Clone("histo_ZllG_combined");
  for(int i = 1; i <= nBins; i++){
    histo_ZllG_combined->SetBinContent(i, histo_ZllG_combined->GetBinContent(i)*histo_ZllG_combined->GetBinWidth(i));
    histo_ZllG_combined->SetBinError(i, histo_ZllG_combined->GetBinError(i)*histo_ZllG_combined->GetBinWidth(i));
    cout<<"ZllG_combined_bin_"<<i<<": error="<<histo_ZllG_combined->GetBinError(i)<<endl;
  }
  Float_t int_ZllG_combined = histo_ZllG_combined->Integral();
  total_background += int_ZllG_combined;
  histo_ZllG_combined->SetFillColor(kSpring-9);
  histo_vector.push_back(histo_ZllG_combined);
  
  TH1F* histo_TTG = (TH1F*)((TH1F*)f_postfit_shapes->Get("shapes_prefit/CR_DE/TTG"))->Clone("histo_TTG");
  for(int i = 1; i <= nBins; i++){
    histo_TTG->SetBinContent(i, histo_TTG->GetBinContent(i)*histo_TTG->GetBinWidth(i));
    histo_TTG->SetBinError(i, histo_TTG->GetBinError(i)*histo_TTG->GetBinWidth(i));
    cout<<"TTG_bin_"<<i<<": error="<<histo_TTG->GetBinError(i)<<endl;
  }
  Float_t int_TTG = histo_TTG->Integral();
  total_background += int_TTG;
  histo_TTG->SetFillColor(kOrange-3);
  histo_vector.push_back(histo_TTG);
  
  TH1F* histo_WZ = (TH1F*)((TH1F*)f_postfit_shapes->Get("shapes_prefit/CR_DE/WZ"))->Clone("histo_WZ");
  for(int i = 1; i <= nBins; i++){
    histo_WZ->SetBinContent(i, histo_WZ->GetBinContent(i)*histo_WZ->GetBinWidth(i));
    histo_WZ->SetBinError(i, histo_WZ->GetBinError(i)*histo_WZ->GetBinWidth(i));
    cout<<"WZ_bin_"<<i<<": error="<<histo_WZ->GetBinError(i)<<endl;
  }
  Float_t int_WZ = histo_WZ->Integral();
  total_background += int_WZ;
  histo_WZ->SetFillColor(kRed-5);
  histo_vector.push_back(histo_WZ);
  
  for(int i = 1; i <= nBins; i++){
    double binWidth = histo_data->GetBinWidth(i);
    histo_data->SetBinContent(i,histo_data->GetBinContent(i)/binWidth);
    histo_data->SetBinError(i,histo_data->GetBinError(i)/binWidth);
    histo_ZllG_combined->SetBinContent(i,histo_ZllG_combined->GetBinContent(i)/binWidth);
    histo_ZllG_combined->SetBinError(i,histo_ZllG_combined->GetBinError(i)/binWidth);
    histo_TTG->SetBinContent(i,histo_TTG->GetBinContent(i)/binWidth);
    histo_TTG->SetBinError(i,histo_TTG->GetBinError(i)/binWidth);
    histo_WZ->SetBinContent(i,histo_WZ->GetBinContent(i)/binWidth);
    histo_WZ->SetBinError(i,histo_WZ->GetBinError(i)/binWidth);
    histo_jetfake->SetBinContent(i,histo_jetfake->GetBinContent(i)/binWidth);
    histo_jetfake->SetBinError(i,histo_jetfake->GetBinError(i)/binWidth);
  }

  TPad *pad1 = new TPad("pad1","pad1",0.01,0.26,0.99,0.99);
  pad1->Draw(); pad1->cd();
  pad1->SetFillColor(0); pad1->SetFrameBorderMode(0); pad1->SetBorderMode(0);
  pad1->SetBottomMargin(0.);
  
  TH1F *histo_allbackgrounds = (TH1F*)histo_ZllG_combined->Clone("histo_allbackgrounds");
  histo_allbackgrounds->Add(histo_TTG);
  histo_allbackgrounds->Add(histo_WZ);
  histo_allbackgrounds->Add(histo_jetfake);
  for(int i = 1; i <= nBins; i++){
    double background = histo_allbackgrounds->GetBinContent(i);
    // Add bin errors
    double sum_binerrors_squared = 0.0;
    sum_binerrors_squared += pow(histo_ZllG_combined->GetBinError(i),2);
    sum_binerrors_squared += pow(histo_TTG->GetBinError(i),2);
    sum_binerrors_squared += pow(histo_WZ->GetBinError(i),2);
    sum_binerrors_squared += pow(histo_jetfake->GetBinError(i),2);
    double binerror = sqrt(sum_binerrors_squared); // Include just the statistical error
    histo_allbackgrounds->SetBinError(i,binerror);
  }
  histo_allbackgrounds->SetFillColorAlpha(kGray+1,0.6);
  histo_vector.push_back(histo_allbackgrounds);
  
  TH1F *histo_allbackgrounds_outline = (TH1F*)histo_allbackgrounds->Clone("histo_allbackgrounds_outline");
  histo_allbackgrounds_outline->SetFillColorAlpha(kWhite,0.0);
  histo_allbackgrounds_outline->SetLineWidth(1);
  histo_vector.push_back(histo_allbackgrounds_outline);
 
  // DEBUG
  // if(histname == "Photon_Et_range_24"){
  //   cout<<"Jet faking photon: "<<histo_jetfake->GetBinContent(nBins)*histo_jetfake->GetBinWidth(nBins)<<" +- "<<histo_jetfake->GetBinError(nBins)*histo_jetfake->GetBinWidth(nBins)<<endl;
  //   cout<<"Spikes: "<<histo_spikes->GetBinContent(nBins)*histo_spikes->GetBinWidth(nBins)<<" +- "<<histo_spikes->GetBinError(nBins)*histo_spikes->GetBinWidth(nBins)<<endl;
  //   cout<<"Electron faking photon: "<<histo_elefake->GetBinContent(nBins)*histo_elefake->GetBinWidth(nBins)<<" +- "<<histo_elefake->GetBinError(nBins)*histo_elefake->GetBinWidth(nBins)<<endl;
  //   cout<<"Beam halo: "<<histo_bhalo->GetBinContent(nBins)*histo_bhalo->GetBinWidth(nBins)<<" +- "<<histo_bhalo->GetBinError(nBins)*histo_bhalo->GetBinWidth(nBins)<<endl;
  //   cout<<"ZNuNu+gamma: "<<histo_ZNuNuG->GetBinContent(nBins)*histo_ZNuNuG->GetBinWidth(nBins)<<" +- "<<histo_ZNuNuG->GetBinError(nBins)*histo_ZNuNuG->GetBinWidth(nBins)<<endl;
  //   cout<<"W+gamma: "<<histo_WG->GetBinContent(nBins)*histo_WG->GetBinWidth(nBins)<<" +- "<<histo_WG->GetBinError(nBins)*histo_WG->GetBinWidth(nBins)<<endl;
  //   cout<<"GJets: "<<histo_GJets_40toInf->GetBinContent(nBins)*histo_GJets_40toInf->GetBinWidth(nBins)<<" +- "<<histo_GJets_40toInf->GetBinError(nBins)*histo_GJets_40toInf->GetBinWidth(nBins)<<endl;
  //   cout<<"Z(ll)+Gamma: "<<histo_ZllG_combined->GetBinContent(nBins)*histo_ZllG_combined->GetBinWidth(nBins)<<" +- "<<histo_ZllG_combined->GetBinError(nBins)*histo_ZllG_combined->GetBinWidth(nBins)<<endl;
  //   cout<<"tt+Gamma: "<<histo_TTG->GetBinContent(nBins)*histo_TTG->GetBinWidth(nBins)<<" +- "<<histo_TTG->GetBinError(nBins)*histo_TTG->GetBinWidth(nBins)<<endl;
  //   cout<<"t+Gamma: "<<histo_TG->GetBinContent(nBins)*histo_TG->GetBinWidth(nBins)<<" +- "<<histo_TG->GetBinError(nBins)*histo_TG->GetBinWidth(nBins)<<endl;
  //   // cout<<"WWG: "<<histo_WWG->GetBinContent(nBins)*histo_WWG->GetBinWidth(nBins)<<" +- "<<histo_WWG->GetBinError(nBins)*histo_WWG->GetBinWidth(nBins)<<endl;
  //   cout<<"Diphoton: "<<histo_diphoton->GetBinContent(nBins)*histo_diphoton->GetBinWidth(nBins)<<" +- "<<histo_diphoton->GetBinError(nBins)*histo_diphoton->GetBinWidth(nBins)<<endl;
  //   cout<<"WZ: "<<histo_WZ->GetBinContent(nBins)*histo_WZ->GetBinWidth(nBins)<<" +- "<<histo_WZ->GetBinError(nBins)*histo_WZ->GetBinWidth(nBins)<<endl;
  //   cout<<"ZZ: "<<histo_ZZ->GetBinContent(nBins)*histo_ZZ->GetBinWidth(nBins)<<" +- "<<histo_ZZ->GetBinError(nBins)*histo_ZZ->GetBinWidth(nBins)<<endl;
  //   cout<<"WMuNu: "<<histo_WMuNu->GetBinContent(nBins)*histo_WMuNu->GetBinWidth(nBins)<<" +- "<<histo_WMuNu->GetBinError(nBins)*histo_WMuNu->GetBinWidth(nBins)<<endl;
  //   cout<<"WTauNu: "<<histo_WTauNu->GetBinContent(nBins)*histo_WTauNu->GetBinWidth(nBins)<<" +- "<<histo_WTauNu->GetBinError(nBins)*histo_WTauNu->GetBinWidth(nBins)<<endl;
  //   cout<<"WW: "<<histo_WW->GetBinContent(nBins)*histo_WW->GetBinWidth(nBins)<<" +- "<<histo_WW->GetBinError(nBins)*histo_WW->GetBinWidth(nBins)<<endl;
  //   // cout<<"WZG: "<<histo_WZG->GetBinContent(nBins)*histo_WZG->GetBinWidth(nBins)<<" +- "<<histo_WZG->GetBinError(nBins)*histo_WZG->GetBinWidth(nBins)<<endl;
  //   // cout<<"WGG: "<<histo_WGG->GetBinContent(nBins)*histo_WGG->GetBinWidth(nBins)<<" +- "<<histo_WGG->GetBinError(nBins)*histo_WGG->GetBinWidth(nBins)<<endl;
  //   // cout<<"ZGGToNuNuGG: "<<histo_ZGGToNuNuGG->GetBinContent(nBins)*histo_ZGGToNuNuGG->GetBinWidth(nBins)<<" +- "<<histo_ZGGToNuNuGG->GetBinError(nBins)*histo_ZGGToNuNuGG->GetBinWidth(nBins)<<endl;
  //   cout<<"Total background: "<<histo_allbackgrounds->GetBinContent(nBins)*histo_allbackgrounds->GetBinWidth(nBins)<<" +- "<<histo_allbackgrounds->GetBinError(nBins)*histo_allbackgrounds->GetBinWidth(nBins)<<endl;
  // }
      
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
  double ymax = 1.7*ymax_data;
  if(ymax_background > ymax_data)
    ymax = 1.7*ymax_background;
  pad1->SetLogy();
  ymax *= 100;
  histo_data->GetYaxis()->SetRangeUser(ymin,ymax);
  histo_data->SetLineColor(kWhite);
  histo_data->SetMarkerColor(kWhite);
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
  leg->AddEntry(histo_WZ,"WZ","F");
  leg->SetNColumns(2);
  leg->SetFillColor(kWhite);
  leg->SetShadowColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.040);
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

  cout<<"ZeeG region"<<endl;
  cout<<"------------------------------------"<<endl;
  cout<<"Z(ll)+Gamma: "<<int_ZllG_combined<<endl;
  cout<<"tt+Gamma: "<<int_TTG<<endl;
  cout<<"WZ: "<<int_WZ<<endl;
  cout<<"Jet faking photon: "<<int_jetfake<<endl;
  cout<<"Total background: "<<total_background<<endl;
  cout<<"------------------------------------"<<endl;
  cout<<"Data: "<<int_data<<endl;
  cout<<"------------------------------------"<<endl;
  cout<<endl;

  string plotTitle = "prefit_zeeg_";
  plotTitle += plotname;
  c->SaveAs(TString(plotTitle+".png"));
  c->SaveAs(TString(plotTitle+".pdf"));
  delete(c);
}

void prefit_zeeg_plotter()
{
  std::vector<Double_t> leg_xoffsets;
  leg_xoffsets.clear();
  std::vector<Double_t> leg_yoffsets;
  leg_yoffsets.clear();
  std::vector<TString> xaxis_titles;
  xaxis_titles.clear();
  std::vector<TString> plotnames;
  plotnames.clear();

//  leg_xoffsets.push_back(0.);
//  leg_yoffsets.push_back(0.);
//  xaxis_titles.push_back(TString(""));
//  plotnames.push_back(TString(""));

  leg_xoffsets.push_back(0.);
  leg_yoffsets.push_back(0.);
  xaxis_titles.push_back(TString("Photon #it{p}_{T} [GeV]"));
  plotnames.push_back(TString("phoPt"));

  for(int i = 0; i < plotnames.size(); i++){
    plot(leg_xoffsets[i],leg_yoffsets[i],xaxis_titles[i],plotnames[i]);
  }
}