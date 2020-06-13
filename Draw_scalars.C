#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <TGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include "TFile.h"
#include <string>
#include <TMultiGraph.h>
#include <TLegend.h>
#include <TColor.h>
#include <math.h>

using namespace std;


//Int_t Draw_scalars()

Int_t Draw_scalars(const char *  	inputFile, const char* outputFile)
{

    TFile* scalar_results = new TFile(outputFile, "RECREATE");


  //  ifstream iFile(inputFile);  
    
    
    
    TString talyName = inputFile;
    ifstream talyInput;
    talyInput.open(talyName);

    if(!talyInput.is_open()){
        cout << "Could not open"<< talyName << endl;
        return kFALSE;
    }


    vector<double> vec_data; 
    float arr_scalars[48][32];
    float arr_reg[48];
    for (int stp=0; stp<32; stp++)
    {
        arr_reg[stp] = stp;
    }
    for(Int_t reg=0; reg<32; reg++){
        for(Int_t ch=0; ch<48; ch++){
            talyInput >> arr_scalars[ch][reg];
        }        
    }
    
    const EColor colours[]={kOrange,kMagenta,kBlue,kCyan,kGray,kBlack,kGreen,kYellow};

    std::vector<double> vec_counttopeak;
    std::vector<double> vec_postopeak;
    std::vector<double> vec_counttopeak1;
    std::vector<double> vec_postopeak1;
    
    TMultiGraph *mg_Count_vs_blposL0 = new TMultiGraph("mg_Count_vs_blposL0","mg_Count_vs_blposL0");
    TMultiGraph *mg_Count_vs_blposR0 = new TMultiGraph("mg_Count_vs_blposR0","mg_Count_vs_blposR0");
    TMultiGraph *mg_Count_vs_blposL1 = new TMultiGraph("mg_Count_vs_blposL1","mg_Count_vs_blposL1");
    TMultiGraph *mg_Count_vs_blposR1 = new TMultiGraph("mg_Count_vs_blposR1","mg_Count_vs_blposR1");
    
    TLegend* cntbl_legL0 = new TLegend(0.7,0.7,0.9,0.9);
    TLegend* cntbl_legL1 = new TLegend(0.7,0.7,0.9,0.9);
    TLegend* cntbl_legR0 = new TLegend(0.7,0.7,0.9,0.9);
    TLegend* cntbl_legR1 = new TLegend(0.7,0.7,0.9,0.9);
    
    TCanvas * c_count_vs_bl_L0= new TCanvas("c_count_vs_bl_L0","c_count_vs_bl_L0");
    TCanvas * c_count_vs_bl_L1 = new TCanvas("c_count_vs_bl_L1","c_count_vs_bl_L1");
    TCanvas * c_count_vs_bl_R0= new TCanvas("c_count_vs_bl_R0","c_count_vs_bl_R0");
    TCanvas * c_count_vs_bl_R1 = new TCanvas("c_count_vs_bl_R1","c_count_vs_bl_R1");
    
    TF1* fl = new TF1 ( "fl", "pol1" );
    TF1* fr = new TF1 ( "fr", "pol1" );

	fl->SetLineColor(kBlue);
    	fr->SetLineColor(kRed);

    TCanvas * c_countvsbl= new TCanvas("c_countvsbl","c_countvsbl");
  //  TCanvas * c_countvsbl_R= new TCanvas("c_countvsbl_R","c_countvsbl_R");
    
    c_countvsbl->Divide(4,4);
   // c_countvsbl_R->Divide(4,4);
    
//     TF1 * fl = new TF1("fl","[0]*TMath::Power(x,[1])",0,30);
//     TF1 * fr = new TF1("fr","[0]*TMath::Power(x,[1])",0,30);

    int cclr=0;
    
    for (int ch =16; ch< 32; ch++)
    {
        double norm_peak = *std::max_element(arr_scalars[ch],arr_scalars[ch]+32);
        double norm_peak_pos = std::max_element(arr_scalars[ch], arr_scalars[ch]+32) - arr_scalars[ch];
        
        for (int reg=0; reg< 32; reg++)
        {
            if (arr_scalars[ch][reg] > 0)
            {
                vec_counttopeak.push_back(arr_scalars[ch][reg]);
                double b0_bl = (norm_peak_pos - reg )*2;
                vec_postopeak.push_back(pow(b0_bl,2));
                //cout<<"{"<<reg<<","<<pow(b0_bl,2)<<","<<arr_scalars[ch][reg]<<"}";
            }
            if (arr_scalars[ch][reg] == *std::max_element(arr_scalars[ch],arr_scalars[ch]+32)) break;
        }
        //cout<<vec_counttopeak.size()<<" :\t";
        
        for (int reg=norm_peak_pos; reg< 32; reg++)
        {
            if (arr_scalars[ch][reg] > 0)
            {
                vec_counttopeak1.push_back(arr_scalars[ch][reg]);
                double b0_bl = (norm_peak_pos - reg )*2;
                vec_postopeak1.push_back(pow(b0_bl,2));
                //cout<<"R: {"<<reg<<","<<pow(b0_bl,2)<<","<<arr_scalars[ch][reg]<<"}";
            }
            if (arr_scalars[ch][reg] == 0) break;
        }
     //   cout<<vec_counttopeak1.size()<<" :\t";
        
        double arr_counttopeak[vec_counttopeak.size()];
        double arr_postopeak[vec_counttopeak.size()];
        double arr_counttopeak1[vec_counttopeak1.size()];
        double arr_postopeak1[vec_counttopeak1.size()];
        
        for (int i =0; i < vec_counttopeak.size(); i++)
        {
            arr_counttopeak[i] = log(vec_counttopeak[i]/norm_peak);
            arr_postopeak[i] = vec_postopeak[i];
          //  cout<<vec_counttopeak[i]<<"-"<<i<<"-"<<vec_postopeak[i]<<"\t";
            
        }
        for (int i =0; i < vec_counttopeak1.size(); i++)
        {
            arr_counttopeak1[i] = log(vec_counttopeak1[i]/norm_peak);
            arr_postopeak1[i] = vec_postopeak1[i];
          //  cout<<vec_counttopeak[i]<<"-"<<i<<"-"<<vec_postopeak[i]<<"\t";
            
        }
      //  cout<<endl<<endl;
        if (vec_counttopeak.size() > 0)
        {
            TGraph *gcount_vs_blposL = new TGraph(vec_counttopeak.size(),arr_postopeak, arr_counttopeak);
            gcount_vs_blposL->SetName(Form("Channel %d",ch));
            gcount_vs_blposL->SetTitle(Form("Channel %d",ch));
            gcount_vs_blposL->SetLineColor(kBlue);
            gcount_vs_blposL->SetLineWidth(3);
            gcount_vs_blposL->SetLineStyle(9);
            //gcount_vs_blposL->SetMarkerStyle(4);
            //gcount_vs_blposL->SetMarkerSize(1);
            gcount_vs_blposL->Fit(fl,"q");
    		gcount_vs_blposL->GetXaxis()->SetTitle("(bl_{peak} - bl) ^{2} mV");
    		gcount_vs_blposL->GetYaxis()->SetTitle("ln(Counts_{ch}/maxCount_{ch})");
            
            TGraph *gcount_vs_blposR = new TGraph(vec_counttopeak1.size(),arr_postopeak1, arr_counttopeak1);
            gcount_vs_blposR->SetName(Form("Channel %d",ch));
            gcount_vs_blposR->SetTitle(Form("Channel %d",ch));
            gcount_vs_blposR->SetLineColor(kRed);
            gcount_vs_blposR->SetLineWidth(3);
            gcount_vs_blposR->SetLineStyle(9);
            //gcount_vs_blposR->SetMarkerStyle(4);
            //gcount_vs_blposR->SetMarkerSize(1);
            gcount_vs_blposR->Fit(fr,"q");
    		gcount_vs_blposR->GetXaxis()->SetTitle("(bl_{peak} - bl) ^{2} mV");
    		gcount_vs_blposR->GetYaxis()->SetTitle("ln(Counts_{ch}/maxCount_{ch})");
            
            double slope_l =  fl->GetParameter(1);
            double slope_r =  fr->GetParameter(1);
            double var_l = -1/(2*slope_l);
            double var_r = -1/(2*slope_r);
            cout<<ch-15<<"\t"<<(sqrt(var_l)+sqrt(var_r))/2<<"\t"<<endl;
            
            if (ch-16 < 8)
            {
                mg_Count_vs_blposL0->Add(gcount_vs_blposL,"LP");
                mg_Count_vs_blposR0->Add(gcount_vs_blposR,"LP");
                cntbl_legL0->AddEntry(gcount_vs_blposL,Form("Channel %d",ch-16));
                cntbl_legR0->AddEntry(gcount_vs_blposR,Form("Channel %d",ch-16));
            }
            else 
            {
                mg_Count_vs_blposL1->Add(gcount_vs_blposL,"LP");
                mg_Count_vs_blposR1->Add(gcount_vs_blposR,"LP");
                cntbl_legL1->AddEntry(gcount_vs_blposL,Form("Channel %d",ch-16));
                cntbl_legR1->AddEntry(gcount_vs_blposR,Form("Channel %d",ch-16));
            }
            c_countvsbl->cd(ch-15);
            gcount_vs_blposL->Draw();
            gcount_vs_blposR->Draw("same");
        }
        vec_counttopeak.clear();
        vec_postopeak.clear();
        vec_counttopeak1.clear();
        vec_postopeak1.clear();
        cclr++;
        if (cclr ==8 ) cclr =0;
    }
    c_countvsbl->Write();
    c_count_vs_bl_L0->cd();
    mg_Count_vs_blposL0->GetXaxis()->SetTitle("(bl_{peak} - bl) ^{2} mV");
    mg_Count_vs_blposL0->GetYaxis()->SetTitle("ln(Counts_{ch}/maxCount_{ch})");
    //mg_Count_vs_blposL0->GetYaxis()->SetRangeUser(0.001,1);
    mg_Count_vs_blposL0->Draw("same");
    cntbl_legL0->Draw("same");
    //c_count_vs_bl_L0->SetLogy();
    c_count_vs_bl_L0->Write();
    
    c_count_vs_bl_L1->cd();
    mg_Count_vs_blposL1->GetXaxis()->SetTitle("(bl_{peak} - bl) ^{2} mV");
    mg_Count_vs_blposL1->GetYaxis()->SetTitle("ln(Counts_{ch}/maxCount_{ch})");
    //mg_Count_vs_blposL1->GetYaxis()->SetRangeUser(0.001,1);
    mg_Count_vs_blposL1->Draw("same");
    cntbl_legL1->Draw("same");
    //c_count_vs_bl_L1->SetLogy();
    c_count_vs_bl_L1->Write();
/*    mg_Count_vs_blpos0->SetLogy();
    mg_Count_vs_blpos1->SetLogy(); */   
    mg_Count_vs_blposL0->Write();
    mg_Count_vs_blposL1->Write();
    
    
    c_count_vs_bl_R0->cd();
    mg_Count_vs_blposR0->GetXaxis()->SetTitle("(bl_{peak} - bl) ^{2} mV");
    mg_Count_vs_blposR0->GetYaxis()->SetTitle("Counts_{ch}/maxCount_{ch}");
    //mg_Count_vs_blposR0->GetYaxis()->SetRangeUser(0.001,1);
    mg_Count_vs_blposR0->Draw("same");
    cntbl_legR0->Draw("same");
    //c_count_vs_bl_R0->SetLogy();
    c_count_vs_bl_R0->Write();
    
    c_count_vs_bl_R1->cd();
    mg_Count_vs_blposR1->GetXaxis()->SetTitle("(bl_{peak} - bl) ^{2} mV");
    mg_Count_vs_blposR1->GetYaxis()->SetTitle("Counts_{ch}/maxCount_{ch}");
    mg_Count_vs_blposR1->GetYaxis()->SetRangeUser(0.001,1);
    mg_Count_vs_blposR1->Draw("same");
    cntbl_legR1->Draw("same");
    //c_count_vs_bl_R1->SetLogy();
    c_count_vs_bl_R1->Write();
/*    mg_Count_vs_blpos0->SetLogy();
    mg_Count_vs_blpos1->SetLogy(); */   
    mg_Count_vs_blposR0->Write();
    mg_Count_vs_blposR1->Write();
    
    
    
    
        
    TGraph *gchannel_vs_scalars[16];
    TMultiGraph *mg_Scalars0 = new TMultiGraph("mg_Scalars0","mg_Scalars0");
    TMultiGraph *mg_Scalars1 = new TMultiGraph("mg_Scalars1","mg_Scalars1");
    TCanvas * c_counts = new TCanvas("c_counts","c_counts");
    TCanvas * c_scalers = new TCanvas("c_scalers","c_scalers");

    TLegend* sd_leg0 = new TLegend(0.7,0.7,0.9,0.9);
    TLegend* sd_leg1 = new TLegend(0.7,0.7,0.9,0.9);
    
    int clr=0;
    c_counts->Divide(4,4);
    double gaus_sigma[16];
    double gaus_index[16];
    
    double f0 = 1/ ((2*sqrt(3))*(20/1000000000));
    double vth2 = 0;
    
    
    
    for (int ch=16; ch<32; ch++)
    {
        gchannel_vs_scalars[ch-16] = new TGraph(32,arr_reg, arr_scalars[ch]);
        gchannel_vs_scalars[ch-16]->SetName(Form("Channel %d",ch));
        gchannel_vs_scalars[ch-16]->SetTitle(Form("Channel %d",ch));
        gchannel_vs_scalars[ch-16]->SetLineColor(colours[clr]);
        gchannel_vs_scalars[ch-16]->SetLineWidth(3);
        gchannel_vs_scalars[ch-16]->SetMarkerStyle(4);
        gchannel_vs_scalars[ch-16]->SetMarkerSize(1);


        double peak_x = std::max_element(arr_scalars[ch], arr_scalars[ch]+32) - arr_scalars[ch];
        TF1* f1 = new TF1 ( "f1", "gaus",peak_x,1 );
        gchannel_vs_scalars[ch-16]->Fit ( f1, "q");
        cout<<"const : "<<f1->GetParameter(0)<<"  mean : "<<f1->GetParameter(1)<<" sigma : "<<f1->GetParameter(2)*2<<endl;
        gaus_sigma[ch-16] = f1->GetParameter(2)*2;
        gaus_index[ch-16] = ch;
        
        delete f1;
        
      //  f1 = gchannel_vs_scalars[ch-16]->GetFunction ( "f1" );
     //   Double_t trck_constant = f1->GetParameter ( 0 ); //constant - c
    //    Double_t trck_slope = f1->GetParameter ( 1 ); //slope    - m

        if (ch-16 < 8)
        {
            mg_Scalars0->Add(gchannel_vs_scalars[ch-16],"LP");
            sd_leg0->AddEntry(gchannel_vs_scalars[ch-16],Form("Channel %d",ch-16));
            
        }
        else 
        {
            mg_Scalars1->Add(gchannel_vs_scalars[ch-16],"LP");
            sd_leg1->AddEntry(gchannel_vs_scalars[ch-16],Form("Channel %d",ch-16));
            
        }
        c_counts->cd(ch-15);
        gchannel_vs_scalars[ch-16]->Draw();
        clr++;
        if (clr ==8 ) clr =0;
    }
    c_counts->Write();
    sd_leg0->SetHeader("ASIC 0","C");
    sd_leg1->SetHeader("ASIC 1","C");
    
/*    
    mg_Scalars0->GetXaxis()->SetTitle("Baseline Register");
    mg_Scalars0->GetYaxis()->SetRangeUser(1,30000000);
    mg_Scalars0->GetYaxis()->SetTitle("Counts");
    mg_Scalars0->GetYaxis()->SetTitleSize(0.04);
    mg_Scalars0->GetXaxis()->SetTitleSize(0.04);
    
    mg_Scalars1->GetXaxis()->SetTitle("Baseline Register");
    mg_Scalars1->GetYaxis()->SetRangeUser(1,30000000);
    mg_Scalars1->GetYaxis()->SetTitle("Counts");
    mg_Scalars1->GetYaxis()->SetTitleSize(0.04);
    mg_Scalars1->GetXaxis()->SetTitleSize(0.04);
  */ 
    c_scalers->Divide(1,2);
    c_scalers->cd(1);
    mg_Scalars0->Draw("same,LP"); 
    sd_leg0->Draw("same");
    mg_Scalars0->Write();
    c_scalers->cd(2);
    mg_Scalars1->Draw("same,LP");
    sd_leg1->Draw("same");
    mg_Scalars1->Write();
    c_scalers->Write();
    
    TGraph * gsigma = new TGraph (16,gaus_index,gaus_sigma);
    gsigma->SetName("gsigma");
    gsigma->SetMarkerStyle(8);
    gsigma->SetMarkerSize(2);
    gsigma->GetXaxis()->SetTitle("Channel No");
    gsigma->GetYaxis()->SetTitle("#sigma mV");
    gsigma->Write();
    
   // TCanvas *can_counts;

    return 0;
}
