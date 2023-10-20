 #include "TFile.h"
#include "TCanvas.h"
#include "TH1.h"
#include "input/Format.C"
#include <iostream>
#include "TStyle.h"

Double_t bwconv0(Double_t *x, Double_t *par);
bool EqArr(int xx0, int xx1, int xx2, int xx3, int yy[4]);

char gg_names[10][200] ={"eta_spectrum_pbsc_hist_FG11","eta_spectrum_pbsc_hist_BG11","sigma_cut_eta_spectrum_pbsc_hist_FG11","sigma_cut_eta_spectrum_pbsc_hist_BG11",
"r_3_eta_spectrum_pbsc_hist_FG11","r_3_eta_spectrum_pbsc_hist_BG11","r_3_cnt_eta_spectrum_pbsc_hist_FG11","r_3_cnt_eta_spectrum_pbsc_hist_BG11",
"emcid_eta_spectrum_pbsc_hist_FG11","emcid_eta_spectrum_pbsc_hist_BG11"}; 
double mass = 0.139;
double Gamma = 0.00991;
void PI_FIT (int hh, int lb, int rb, double pt_min, double pt_max, int hist_N = 7, int meth = 0, bool smooth = true){
    char filename[200] = "input/HeAu_18968_Artem_newcuts.root";    ///HeAu_19829_Artem_onlyemcid Lambda_AuAu_V1_pm_g.root
    char histname_bg[200], histname_fg[200];
    sprintf(histname_bg,gg_names[2*hh+1]);///"LAM_V2_BG12"eta_prime_spectrum_nopid_cut_BG12;
    sprintf(histname_fg,gg_names[2*hh]);///
    int NVertex = 6;
    int Nrp = 1;
    double low_edge_int = 0.7;
    double high_edge_int = 0.99;
    double fit1_l = 0.04;
    double fit1_h =  0.29;
    double dif_coeff = 0.995;
    int reb = 10;
    int proj_axis = 0; /// 0 - X ; 1 - Y
    char outputname[200];
    sprintf(outputname,"output/fits/Kstar_%d_%.1f_%.1f.png",hh,pt_min,pt_max);
    //sprintf(outputname,"output/test.png");


 gStyle->SetOptStat(000);
 gStyle->SetOptFit(1111);
 //1****************start initialization************************************

 const int centr_left=lb;
 const int centr_right=rb;
 double sim;
 char slicename[200], slicename1[200],slicename2[200], slicenamebg[200], dirname[200];
 double gamma_keff[4]={1,2,1,2};
 double pt_bin[2];
 double central_bin[] = {0.,20, 20,40,40, 60,60,80,90, 93};



    TFile *inFile;
    inFile = TFile::Open(filename);
    inFile->cd("c00_z02_r00");

    TH2D *test = (TH2D*) gROOT->FindObject(histname_fg);
    TH1D *test_X = test->ProjectionX("test_X");
    TH1D *test_Y = test->ProjectionY("test_Y");

    int NbinsX = test->GetNbinsX();
    int NbinsY = test->GetNbinsY();
    double low_edgeX = test_X->GetBinLowEdge(1);
    double low_edgeY = test_Y->GetBinLowEdge(1);
    double up_edgeX = test_X->GetBinLowEdge(NbinsX+1);
    double up_edgeY = test_Y->GetBinLowEdge(NbinsY+1);

    int Nbins = NbinsX;
    double low_edge = low_edgeX;
    double up_edge = up_edgeX;
    if(proj_axis == 1)
    {
        Nbins = NbinsY;
        low_edge = low_edgeY;
        up_edge = up_edgeY;
    }
    TH2D *sum, *sumbg;
  sprintf(slicename,"#pi^{+}#pi^{-}, Pt: %.1f - %.1f GeV/c, %.0f - %.0f%%",pt_min,pt_max,central_bin[centr_left],central_bin[centr_right-1]);
  sum = new TH2D(slicename, slicename, NbinsX,low_edgeX,up_edgeX,NbinsY,low_edgeY,up_edgeY);
  sumbg = new TH2D("sumbg",  "sumbg", NbinsX,low_edgeX,up_edgeX,NbinsY,low_edgeY,up_edgeY);


  TH1F *stat = new TH1F("stat",  "stat",  256, -0.5, 255.5);
  sum->Reset(); sumbg->Reset();

  sprintf(slicename,"pK^{#minus}, Pt: %.1f #minus %.1f GeV/c, %.0f #minus %.0f%%",pt_min,pt_max,central_bin[centr_left],central_bin[centr_right-1]);
  cout<<slicename<<endl;
  sprintf(slicename1,"Slice_pt_1: %f - %f",pt_min,pt_max);
  sprintf(slicename2,"Slice_pt_2: %f - %f",pt_min,pt_max);
  sprintf(slicenamebg,"Slicebg_pt: %f - %f",pt_min,pt_max);
  TH1D* slice0, *slicebg0,* slice,* slice1,* slice2,* slicebg;
  slice = new TH1D(slicename, slicename, Nbins,low_edge,up_edge);
  slice1 = new TH1D(slicename1, slicename1, Nbins,low_edge,up_edge);
  slice2 = new TH1D(slicename2, slicename2, Nbins,low_edge,up_edge);
  slicebg = new TH1D(slicenamebg, slicenamebg, Nbins,low_edge,up_edge);

  //1****************end of initialization************************************
  //1****************start calculating histos*********************************
    double int_fg, int_bg;
    //Granicy integrirovaniia FG i BG dlia vychitaniya
    int left = slice->FindBin(low_edge_int);
    int right = slice->FindBin(high_edge_int);

    pt_bin[0] = (sum->ProjectionY(""))->FindBin(pt_min);
    pt_bin[1] = (sum->ProjectionY(""))->FindBin(pt_max);
    if(proj_axis == 1)
    {
        pt_bin[0] = (sum->ProjectionX(""))->FindBin(pt_min);
        pt_bin[1] = (sum->ProjectionX(""))->FindBin(pt_max);
    }
    cout<<pt_bin[0]<<" "<<pt_bin[1]<<endl;


  //Sobiraem histogrammu po klassam (10x1x12 = 120)
   for (int i=centr_left; i<centr_right; i++) //Centrality
    {
      for (int j=0; j<NVertex; j++) //Vertex
      {
          for (int k=0; k<Nrp; k++) //RB
          {
              sprintf(dirname, "c0%d_z0%d_r0%d", i, j, k);
              inFile->cd(dirname);

              stat->Add((TH1D*) gROOT->FindObject("PoolStatistics"));
              sum->Add((TH2D*) gROOT->FindObject(histname_fg));
              sumbg->Add((TH2D*) gROOT->FindObject(histname_bg));

	          if(proj_axis == 0)
              {
                  slice->Add(sum ->ProjectionX("",pt_bin[0],pt_bin[1]));
                  slice0 = sum ->ProjectionX("",pt_bin[0],pt_bin[1]);
                  int_fg = slice0->Integral(left,right);
                  slice1->Add(sum ->ProjectionX("",pt_bin[0],pt_bin[+1]));
                  slicebg0 = sumbg ->ProjectionX("",pt_bin[0],pt_bin[1]);
              } else
              {
                  slice->Add(sum ->ProjectionY("",pt_bin[0],pt_bin[1]));
                  slice0 = sum ->ProjectionY("",pt_bin[0],pt_bin[1]);
                  int_fg = slice0->Integral(left,right);
                  slice1->Add(sum ->ProjectionY("",pt_bin[0],pt_bin[+1]));
                  slicebg0 = sumbg ->ProjectionY("",pt_bin[0],pt_bin[1]);
              }

              slicebg0->Sumw2();
              int_bg = slicebg0->Integral(left,right);
              if (int_bg==0||int_fg==0) continue;
              slicebg0->Scale(int_fg/int_bg*dif_coeff);
              slicebg->Add(slicebg0);

              if(smooth) slice->Add(slicebg0,-1);

              sum->Add(sum,-1.0);
              sumbg->Add(sumbg,-1.0);
          }
      }
    }

   if(!smooth)
   {
       slicebg->Scale(slice->Integral(left,right)/slicebg->Integral(left,right));
       slice->Add(slicebg,-dif_coeff);
   }

   slicebg->Sumw2();
   slice->Rebin(reb);
   slicebg->Rebin(reb);
   slice1->Rebin(reb);

  //1****************end of calculating histos*********************************
  //1****************start drawing histos**************************************
  TCanvas *c1 = new TCanvas("c1","c1",1200,1200);

  cout<<"Integration range: "<<slice->GetBinCenter(left)<<" - "<<slice->GetBinCenter(right)<<endl;
  slice->SetAxisRange(fit1_l-0.005, fit1_h+0.005);
  slice1->SetAxisRange(fit1_l-0.005, fit1_h+0.005);
  slicebg->SetAxisRange(fit1_l-0.005, fit1_h+0.005);
  slicebg->SetLineColor(4);
  slicebg->SetLineColor(2);
  slice->Draw();
  slice->SetXTitle("M_{inv}");
  c1->Update();


    //1****************start fitting histos**************************************

    TF1 *parab = new TF1("parab","([0] - [1]/2.*([4]*[4]-[5]*[5]) - [2]/3.*([4]*[4]*[4]-[5]*[5]*[5]) - [3]/4.*([4]*[4]*[4]*[4]-[5]*[5]*[5]*[5])) /[6] + [1]*x + [2]*x*x + [3]*x*x*x",fit1_l,fit1_h);

    TF1 *bwconv = new TF1("bwconv",bwconv0,fit1_l,fit1_h,8);

    TF1 *mass_res = new TF1("mass_res","[0]+[1]*x+[2]*x*x+[3]*x*x*x",0.5,8.0);
    mass_res->SetParameter(0,1.320407e-01);
    mass_res->SetParameter(1,4.683030e-01);
    mass_res->SetParameter(2,-6.037768e-02);
    mass_res->SetParameter(3, 2.746577e-03);

    TF1 *gamma_bw = new TF1("gamma_bw","[0]+[1]*x+[2]*x*x",0.5,8.0);
    gamma_bw->SetParameter(0, 3.549325e+000);
    gamma_bw->SetParameter(1, 7.769890e-001);
    gamma_bw->SetParameter(2, -4.021910e-002);

    double bin_center;


    double gauss;
    gauss = (mass_res->Eval(bin_center))/1000;

    double gamma;
    gamma = (gamma_bw->Eval(bin_center))/1000;

    cout<<endl;
    cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;

    cout<<"Current bin center is "<<bin_center<<" GeV/c"<<endl;
    cout<<"Expecting gaus: "<<gauss<<endl;
    cout<<"Expecting: Gamma: "<<gamma<<endl;
    cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
    cout<<endl;

    //Gaus
    bwconv->SetParameter(3, gauss );
    bwconv->SetParLimits(3, 0.9*gauss, 1.1*gauss);
    //Gamma
    bwconv->SetParLimits(0, Gamma*0.8, Gamma*1.5);
    bwconv->SetParameter(0, Gamma );
    //Integral
    bwconv->SetParameter(2, 10000 );
    //Massa
    bwconv->SetParameter(1, mass);
    bwconv->SetParLimits(1, mass-0.5*Gamma,mass+0.5*Gamma);


    bwconv->SetParNames("#Gamma","Mass","Peak","#sigma","p0","p1","p2","p3");
    //bwconv->SetParLimits(6, -3.594e9,3.9541e03);
    //bwconv->SetParLimits(4,0, 610);


    slice->Fit("bwconv","m","",fit1_l,fit1_h);

    parab->SetParameters(bwconv->GetParameters()+4);
    parab->SetParErrors(bwconv->GetParErrors()+4);

    parab->FixParameter(4, mass+Gamma*1.5);
    parab->FixParameter(5, mass-Gamma*1.5);
    parab->FixParameter(6, Gamma*1.5*2);
    parab->SetLineColor(2);
    parab->SetLineColor(2);
    parab->Draw("same");

    //1****************end of fitting and drawing**************************************
 //1****************start calculating yeilds and etrrors*****************************
  for (int i=0; i<2000; i++)
   {
     slice2 -> SetBinContent(i, (slice->GetBinContent(i)) - (parab->Eval(slice->GetBinCenter(i))) );
   }

  slice1->SetAxisRange(fit1_l, fit1_h);
// slice2 -> Draw();

  double gran2 = slice->FindBin(bwconv->GetParameter(1)-bwconv->GetParameter(0)*2);//1.01;
  double gran3 = slice->FindBin(bwconv->GetParameter(1)+bwconv->GetParameter(0)*2);//1.028;
  double intgr = 0;
  double intgr1 = 0;

  for (int i=0; i<4000; i++)
   {
     if ( i > gran2 && i <= gran3 )
      {
   //cout<<i<<"  "<<slice1->GetBinContent(i)<<endl;
       intgr = intgr + slice2->GetBinContent(i);
       intgr1 = intgr1 + slice1->GetBinContent(i);
      }
   }

  double bins = slice1->GetBinCenter(2) - slice1->GetBinCenter(1);
  double err1 = sqrt(intgr1)/intgr;
  double err2 = (bwconv->GetParError(4))/(bwconv->GetParameter(4));
  cout<<"Err1: "<<err1<<"  Err2: "<<err2<<"  "<<bwconv->GetParError(4)<<" / "<<bwconv->GetParameter(4)<<endl;
  cout<<"Err1: "<<err1<<"  Err2: "<<err2<<"  "<<parab->GetParError(0)<<" / "<<parab->GetParameter(0)<<endl;

  cout<<"Total integral: "<<(int) intgr<<" +/- "<<	(int) sqrt( intgr1 ) <<endl;
//   cout<<"Total integral (bad): "<<(bwconv->GetParameter(2))/bins<<" +/- "<<(bwconv->GetParError(2))/bins<<endl;
   cout<<"Total Stats: "<<stat->GetBinContent(2)/1e6<<" *10^6"<<endl;

    char integral[200];
    sprintf(integral,"%.f/1e6",intgr/stat->GetBinContent(2)*1e6);
    DrawLegendTitle(0.23,0.85,0.24,0.86,integral,0.04,22);
    sprintf(integral,"#pm%.3f%%",100*(1*sqrt( intgr1 )/intgr));
    DrawLegendTitle(0.23,0.8,0.24,0.81,integral,0.04,22);

   c1->Update();
    c1->SaveAs(outputname);


 //1****************end of calculating yeilds and etrrors*****************************

}



 Double_t bwconv0(Double_t *x, Double_t *par)
 {

     // Numeric constants
     Double_t invsq2pi = 0.3989422804014; // (2 pi)^(-1/2)
     Double_t twoPi = 6.2831853071795;//2Pi

     // Control constants
     Double_t np = 100.0; // number of convolution steps
     Double_t sc = 4; // convolution extends to +-sc Gaussian sigmas

     // Variables
     Double_t xx;
     Double_t fland;
     Double_t sum = 0.0;
     Double_t xlow,xupp;
     Double_t step;
     Double_t i;

     double gran1 = (mass-Gamma*3);//1.0;
     double gran2 = (mass-Gamma*1.5);//1.01;
     double gran3 = (mass+Gamma*1.5);//1.028;
     double gran4 = (mass+Gamma*3);//1.04;


     // Range of convolution integral
     xlow = x[0] - sc * par[3];
     xupp = x[0] + sc * par[3];

     step = (xupp-xlow) / np;

     // Convolution integral of Breit and Gaussian by sum
     for(i=1.0; i<=np/2; i++) {
         xx = xlow + (i-.5) * step;
         fland = TMath::BreitWigner(xx,par[1],par[0]);
         sum += fland * TMath::Gaus(x[0],xx,par[3]);

         xx = xupp - (i-.5) * step;
         fland = TMath::BreitWigner(xx,par[1],par[0]);
         sum += fland * TMath::Gaus(x[0],xx,par[3]);
     }

     double deltax = gran3 - gran2;
     double x2 = gran3;
     double x1 = gran2;

     return (par[2] * step * sum * invsq2pi / par[3])+
            ( (par[4] - par[5]/2*(x2*x2-x1*x1) - par[6]/3*(x2*x2*x2-x1*x1*x1) - par[7]/4*(x2*x2*x2*x2-x1*x1*x1*x1)) / deltax ) + par[5]*x[0]+par[6]*x[0]*x[0]+par[7]*x[0]*x[0]*x[0];
 }
 bool EqArr(int xx0, int xx1, int xx2, int xx3, int yy[4])
 {
    bool ret = true;
    int xx[4]={xx0, xx1, xx2, xx3};
     for (int ii = 0; ii < 4; ++ii) {
         if(xx[ii]!=yy[ii]) ret=false;
     }
     return ret;
 }
