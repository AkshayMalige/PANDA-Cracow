#include <iostream>
#include <fstream>
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

using namespace std;


class scan {
   public:
	    	string tdc;
		int card;
		int asic;
		int channel;
		int width;
		int peak;
		float mean_bl;
		float sd;
		int status;  
	    	string file_name;

};
//Int_t Fee_scan_analysis(const char *  	inputFile, const char* outputFile)
Int_t Fee_scan_analysis2(int fn)
{

  //  const char* filename[x];
    TString talysName = "a.txt";  //a has all the file names you must read
    fstream talysInput(talysName,ios::in);
    if(!talysInput.is_open())
{
cout << "Could not open"<< talysName << endl;
return kFALSE;
}

    const Int_t max =fn;  // max is the total number of files in the list a.txt

    TString file_name[max];

    for(Int_t i=0; i<max; i++)
    {
        talysInput >> file_name[i];
    }

    talysInput.close();
    
    for (int k=0; k<max; k++)
    {
    
    	cout<<file_name[k]<<endl;

    	TString f1 = file_name[k];
	TString ot = f1 + ".root"; 

//////////////////////////////////////////////////////////////////////////////////////////

    	TFile* scan_results = new TFile(ot, "RECREATE");


    	ifstream iFile(f1);  

    	scan sc_obj;

    	vector<scan> vec_data; 
 

	    while (!iFile.eof())
	    {
	    	string tdc;
		int card;
		int asic;
		double channel;
		int width;
		int peak;
		float mean_bl;
		float sd;
		int status;  
	    	string file_name;

		iFile >> tdc >>card >> asic >> channel >> width >> peak >> mean_bl >> sd >> status >> file_name;

		sc_obj.tdc = tdc;
		sc_obj.card = card;
		sc_obj.asic = asic;
		sc_obj.channel = channel;
		sc_obj.width = width;
		sc_obj.peak = peak;
		sc_obj.mean_bl = mean_bl;
		sc_obj.sd = sd;
		sc_obj.status = status;
		sc_obj.file_name = file_name;

		vec_data.push_back(sc_obj);

	    }
	    
		TH2F* h0 = new TH2F("h0","Case 0;Baseline position;Noise Width [mV]",34,-1,33,34,-1,33);
		TH2F* h1 = new TH2F("h1","Case 1;Baseline position;Noise Width [mV]",34,-1,33,34,-1,33);
		TH2F* h2 = new TH2F("h2","Case 2;Baseline position;Noise Width [mV]",34,-1,33,34,-1,33);
		TH2F* h3 = new TH2F("h3","Case 3;Baseline position;Noise Width [mV]",34,-1,33,34,-1,33);
		TH2F* h4 = new TH2F("h4","Case 4;Baseline position;Noise Width [mV]",35,-1,34,35,-1,34);
		TH1F* h_count = new TH1F("h_count","Baseline_scan",5,0,5);
		TH1F* h1_width = new TH1F("h1_width","Baseline_scan;Width [mV]",32,0,32);
		TH1F* h1_peak = new TH1F("h1_peak","Baseline_scan;Baseline position [mV]",32,0,32);

		TH2D* channel_vs_width = new TH2D("channel_vs_width","channel_vs_width; Channel No;Noise Width [mV]",50,0,50,64,0,64);
	//	TH2D* channel_vs_position = new TH2D("channel_vs_position","channel_vs_position; Channel No;Baseline position [mV]",50,0,50,64,32,32);
	    TH2D* channel_vs_meanBL = new TH2D("channel_vs_meanBL","channel_vs_meanBL; Channel No",50,0,50,64,-32,32);
	    TH2D* channel_vs_SD = new TH2D("channel_vs_SD","channel_vs_SD; Channel No",50,0,50,10,0,5);

		double bl_clannel[vec_data.size()-1];
		double bl_mean[vec_data.size()-1];
		double bl_sd[vec_data.size()-1];
		double bl_width[vec_data.size()-1];

		std::vector<std::vector<scan>> vec_cases;
		for (int a=0; a< 5; a++)
		{
			std::vector<scan> vec_case_i;
			for ( int aa =0; aa <vec_data.size()-1; aa++)
			{
				if (a == vec_data[aa].status)
				{
					vec_case_i.push_back(vec_data[aa]);	
				}
			}
			vec_cases.push_back(vec_case_i);
		}

		double case_ch[5][vec_data.size()-1];
		double case_bl[5][vec_data.size()-1];
		double case_sd[5][vec_data.size()-1];
		double case_rms2[5][vec_data.size()-1];

	for (int st=0; st<5; st++)
	{
		std::vector<scan> vec_A;
		for(int vst=0; vst<vec_cases[st].size(); vst++)
		{
			vec_A = vec_cases[st];
			case_ch[st][vst] = (16.0*vec_A[vst].card)+vec_A[vst].channel;
            case_bl[st][vst] = vec_A[vst].mean_bl;
            case_sd[st][vst] = vec_A[vst].sd;
            case_rms2[st][vst] = vec_A[vst].width ;
            //* 0.2886 ; // [ 1/sqrt(12) ]* Vp for scqaure wave 
            //cout<<st<<"\t"<<vec_A[vst].channel<<"\t"<<vec_A[vst].sd<<endl;
		}

	}


	for ( int i =0; i <vec_data.size()-1; i++)
	{
		if (vec_data[i].status==0)
		{
				h0->Fill(vec_data[i].peak,vec_data[i].width);
				h_count->Fill(0);

		}
		else if (vec_data[i].status==1)
		{
				h1->Fill(vec_data[i].peak,vec_data[i].width);
				h_count->Fill(1);
				h1_width->Fill(vec_data[i].width);
				h1_peak->Fill(vec_data[i].peak);

		}
		else if (vec_data[i].status==2)
		{
				h2->Fill(vec_data[i].peak,vec_data[i].width);
				h_count->Fill(2);


		}
		else if (vec_data[i].status==3)
		{
				h3->Fill(vec_data[i].peak,vec_data[i].width);
				h_count->Fill(3);


		}
		else if (vec_data[i].status==4)
		{
				h4->Fill(vec_data[i].peak,vec_data[i].width);
				h_count->Fill(4);


		}

		channel_vs_width->Fill(((16*vec_data[i].card)+vec_data[i].channel),vec_data[i].width);
	//	channel_vs_position->Fill(((16.0*vec_data[i].card)+vec_data[i].channel),vec_data[i].peak);
		channel_vs_SD->Fill(((16*vec_data[i].card)+vec_data[i].channel),vec_data[i].sd);
		channel_vs_meanBL->Fill(((16*vec_data[i].card)+vec_data[i].channel),vec_data[i].mean_bl);
		bl_clannel[i]  = ((16*vec_data[i].card)+vec_data[i].channel);
		bl_mean[i]     = vec_data[i].mean_bl;
		bl_sd[i]       = vec_data[i].sd;
		bl_width[i]      = vec_data[i].width;
	}

		TGraph *gbl = new TGraph(vec_data.size()-1,bl_clannel, bl_mean);
		gbl->Draw("AP");
		gbl->SetMarkerSize(2);
		gbl->SetMarkerStyle(8);
        gbl->SetMarkerColor(kRed);
		gbl->GetXaxis()->SetRangeUser(16,31.5);
		gbl->GetYaxis()->SetRangeUser(-32,32);
		gbl->GetXaxis()->SetTitle("Channel No");
		gbl->GetYaxis()->SetTitle("Bl_{mean} mV");
		gbl->GetYaxis()->SetTitleSize(0.04);
		gbl->GetXaxis()->SetTitleSize(0.04);
	    gbl->GetXaxis()->SetNdivisions(520);
		gbl->SetName("gmean_BL");	
		gbl->Write();
        
        TGraph *gsd = new TGraph(vec_data.size()-1,bl_clannel, bl_sd);
		gsd->Draw("AP");
		gsd->SetMarkerSize(2);
		gsd->SetMarkerStyle(8);
        gsd->SetMarkerColor(kRed);
		gsd->GetXaxis()->SetRangeUser(16,31.5);
		gsd->GetYaxis()->SetRangeUser(0,5);
		gsd->GetXaxis()->SetTitle("Channel No");
		gsd->GetYaxis()->SetTitle("#sigma mV");
		gsd->GetYaxis()->SetTitleSize(0.04);
		gsd->GetXaxis()->SetTitleSize(0.04);
	    gsd->GetXaxis()->SetNdivisions(520);
		gsd->SetName("gStd_dev");	
		gsd->Write();
        
        TGraph *gwidth = new TGraph(vec_data.size()-1,bl_clannel, bl_width);
		gwidth->Draw("AP");
		gwidth->SetMarkerSize(2);
		gwidth->SetMarkerStyle(8);
        gwidth->SetMarkerColor(kRed);
		gwidth->GetXaxis()->SetRangeUser(16,31.5);
		gwidth->GetYaxis()->SetRangeUser(0,64);
		gwidth->GetXaxis()->SetTitle("Channel No");
		gwidth->GetYaxis()->SetTitle("Full width mV");
		gwidth->GetYaxis()->SetTitleSize(0.04);
		gwidth->GetXaxis()->SetTitleSize(0.04);
	    gwidth->GetXaxis()->SetNdivisions(520);
		gwidth->SetName("gFull_width");	
		gwidth->Write();
        

	    int arr_marker_style[5] = {4,8,27,28,30};
	    
	    TGraph *gchannel_vs_meanBL[5];
	    TMultiGraph *mg_BL = new TMultiGraph("mg_BL","mg_BL");
	    TCanvas * c_bl = new TCanvas("c_bl","c_bl");
	    c_bl->cd();
	    for (int cs=0; cs<5; cs++)
	    {
            gchannel_vs_meanBL[cs] = new TGraph(vec_cases[cs].size(),case_ch[cs], case_bl[cs]);
            gchannel_vs_meanBL[cs]->SetMarkerSize(2);
            gchannel_vs_meanBL[cs]->SetMarkerColor(kRed);
            gchannel_vs_meanBL[cs]->SetMarkerStyle(arr_marker_style[cs]);
            gchannel_vs_meanBL[cs]->GetXaxis()->SetRangeUser(15.5,31.5);
            gchannel_vs_meanBL[cs]->GetYaxis()->SetRangeUser(-32,32);
            gchannel_vs_meanBL[cs]->SetName(Form("SD Case %d",cs));
            gchannel_vs_meanBL[cs]->GetXaxis()->SetTitle("Channel No");
            gchannel_vs_meanBL[cs]->GetYaxis()->SetTitle("Bl_{mean} mV");
            gchannel_vs_meanBL[cs]->GetYaxis()->SetTitleSize(0.04);
            gchannel_vs_meanBL[cs]->GetXaxis()->SetTitleSize(0.04);
            //gchannel_vs_meanBL[cs]->SetTitle("SD");
            mg_BL->Add(gchannel_vs_meanBL[cs],"AP");
	    }    
	    mg_BL->GetXaxis()->SetTitle("Channel No");
	    mg_BL->GetYaxis()->SetTitle("Bl_{mean} mV");
	    mg_BL->GetYaxis()->SetTitleSize(0.04);
	    mg_BL->GetXaxis()->SetTitleSize(0.04);
	    mg_BL->GetXaxis()->SetRangeUser(15.5,31.5); 
	    mg_BL->GetYaxis()->SetRangeUser(-32,32);
	    mg_BL->GetXaxis()->SetNdivisions(520);
	    mg_BL->Draw("AP");    
	    TLegend* bl_leg = new TLegend(0.7,0.7,0.9,0.9);
	    for(int c=0; c<5; c++)
	    {
		bl_leg->AddEntry(gchannel_vs_meanBL[c],Form("Case %d",c));

	    }
	    bl_leg->Draw();
	    mg_BL->Write();
	    c_bl->SetGridx();
	    c_bl->Write();
	    

		TString s_bl = f1 + "_BL.png"; 
		c_bl->SaveAs(s_bl,"png");


	    TGraph *gchannel_vs_SD[5];
	    TMultiGraph *mg_SD = new TMultiGraph("mg_SD","mg_SD");
	    TCanvas * c_sd = new TCanvas("c_sd","c_sd");
	    c_sd->cd();
	    for (int cs=0; cs<5; cs++)
	    {
            gchannel_vs_SD[cs] = new TGraph(vec_cases[cs].size(),case_ch[cs], case_sd[cs]);
            gchannel_vs_SD[cs]->SetMarkerSize(2);
            gchannel_vs_SD[cs]->SetMarkerColor(kRed);
            gchannel_vs_SD[cs]->SetMarkerStyle(arr_marker_style[cs]);
            // gchannel_vs_SD[cs]->SetMarkerStyle(8);
            gchannel_vs_SD[cs]->GetXaxis()->SetRangeUser(15.5,31.5);
            gchannel_vs_SD[cs]->SetName(Form("SD Case %d",cs));
            gchannel_vs_SD[cs]->GetXaxis()->SetTitle("Channel No");
            gchannel_vs_SD[cs]->GetYaxis()->SetTitle("#sigma mV");
            gchannel_vs_SD[cs]->GetYaxis()->SetTitleSize(0.04);
            gchannel_vs_SD[cs]->GetXaxis()->SetTitleSize(0.04);
            //gchannel_vs_SD[cs]->SetTitle("SD");
            mg_SD->Add(gchannel_vs_SD[cs],"AP");
	    }    
	    mg_SD->GetXaxis()->SetTitle("Channel No");
	    mg_SD->GetYaxis()->SetTitle("#sigma mV");
	    mg_SD->GetYaxis()->SetTitleSize(0.04);
	    mg_SD->GetXaxis()->SetTitleSize(0.04);
	    mg_SD->GetXaxis()->SetRangeUser(15.5,31.5);
	    mg_SD->GetYaxis()->SetRangeUser(0,5);
	    mg_SD->GetXaxis()->SetNdivisions(520);
	    mg_SD->Draw("AP");    
	    TLegend* sd_leg = new TLegend(0.7,0.7,0.9,0.9);
	    for(int c=0; c<5; c++)
	    {
		sd_leg->AddEntry(gchannel_vs_SD[c],Form("Case %d",c));

	    }
	   sd_leg->Draw();
	    mg_SD->Write();
	    c_sd->SetGridx();
	    c_sd->Write();
		TString s_sd = f1 + "_SD.png"; 
		c_sd->SaveAs(s_sd,"png");

	    

	    TGraph *gchannel_vs_RMS[5];
	    TMultiGraph *mg_RMS = new TMultiGraph("mg_Width","mg_Width");
	    TCanvas * c_rms = new TCanvas("c_rms","c_rms");
	    c_rms->cd();
	    for (int cs=0; cs<5; cs++)
	    {
            gchannel_vs_RMS[cs] = new TGraph(vec_cases[cs].size(),case_ch[cs], case_rms2[cs]);
            gchannel_vs_RMS[cs]->SetMarkerSize(2);
            gchannel_vs_RMS[cs]->SetMarkerColor(kRed);
            gchannel_vs_RMS[cs]->SetMarkerStyle(arr_marker_style[cs]);
            gchannel_vs_RMS[cs]->GetXaxis()->SetRangeUser(15.5,31.5);
            gchannel_vs_RMS[cs]->SetName(Form("SD Case %d",cs));
            gchannel_vs_RMS[cs]->GetXaxis()->SetTitle("Channel No");
            gchannel_vs_RMS[cs]->GetYaxis()->SetTitle("Full width [mV]");
            gchannel_vs_RMS[cs]->GetYaxis()->SetTitleSize(0.04);
            gchannel_vs_RMS[cs]->GetXaxis()->SetTitleSize(0.04);
            //gchannel_vs_RMS[cs]->SetTitle("RMS");
            mg_RMS->Add(gchannel_vs_RMS[cs],"AP");
	    }    
	    mg_RMS->GetXaxis()->SetTitle("Channel No");
	    mg_RMS->GetYaxis()->SetTitle("Full width [mV]");
	    mg_RMS->GetYaxis()->SetTitleSize(0.04);
	    mg_RMS->GetXaxis()->SetTitleSize(0.04);
	    mg_RMS->GetXaxis()->SetRangeUser(15.5,31.5);
	    mg_RMS->GetYaxis()->SetRangeUser(0,64);
	    mg_RMS->GetXaxis()->SetNdivisions(520);
	    mg_RMS->Draw("AP");    
	    TLegend* rms_leg = new TLegend(0.7,0.7,0.9,0.9);
	    for(int c=0; c<5; c++)
	    {
            rms_leg->AddEntry(gchannel_vs_RMS[c],Form("Case %d",c));
	    }
	    
	    rms_leg->Draw();
	    mg_RMS->Write();
	    c_rms->SetGridx();
	    c_rms->Write();  
		TString s_width = f1 + "_Width.png"; 
		c_rms->SaveAs(s_width,"png");  


		h_count->GetXaxis()->SetBit(TAxis::kLabelsHori);
	   	h_count->GetXaxis()->SetBinLabel(1,"Case 0");
	   	h_count->GetXaxis()->SetBinLabel(2,"Case 1");
		h_count->GetXaxis()->SetBinLabel(3,"Case 2");
	   	h_count->GetXaxis()->SetBinLabel(4,"Case 3");
	   	h_count->GetXaxis()->SetBinLabel(5,"Case 4");
		h_count->SetLineWidth ( 3 );
	    h_count->SetLineColor ( kTeal-8 );
	    h_count->SetFillColor ( kTeal-8 );
		h_count->GetXaxis()->SetLabelSize ( 0.05 );
	    h_count->GetYaxis()->SetLabelSize ( 0.05 );
	    h_count->GetXaxis()->SetTitleSize ( 0.05 );

		h0->GetXaxis()->SetLabelSize ( 0.05 );
	    h0->GetYaxis()->SetLabelSize ( 0.05 );
	    h0->GetXaxis()->SetTitleSize ( 0.05 );

		h1->GetXaxis()->SetLabelSize ( 0.05 );
	    h1->GetYaxis()->SetLabelSize ( 0.05 );
	    h1->GetXaxis()->SetTitleSize ( 0.05 );

		h1->GetXaxis()->SetLabelSize ( 0.05 );
	    h1->GetYaxis()->SetLabelSize ( 0.05 );
	    h1->GetXaxis()->SetTitleSize ( 0.05 );

		h2->GetXaxis()->SetLabelSize ( 0.05 );
	    h2->GetYaxis()->SetLabelSize ( 0.05 );
	    h2->GetXaxis()->SetTitleSize ( 0.05 );

		h3->GetXaxis()->SetLabelSize ( 0.05 );
	    h3->GetYaxis()->SetLabelSize ( 0.05 );
	    h3->GetXaxis()->SetTitleSize ( 0.05 );

		h4->GetXaxis()->SetLabelSize ( 0.05 );
	    h4->GetYaxis()->SetLabelSize ( 0.05 );
	    h4->GetXaxis()->SetTitleSize ( 0.05 );

		h1_width->GetXaxis()->SetLabelSize ( 0.05 );
	    h1_width->GetYaxis()->SetLabelSize ( 0.05 );
	    h1_width->GetXaxis()->SetTitleSize ( 0.05 );
		h1_width->SetLineWidth ( 2 );

		h1_peak->GetXaxis()->SetLabelSize ( 0.05 );
	    h1_peak->GetYaxis()->SetLabelSize ( 0.05 );
	    h1_peak->GetXaxis()->SetTitleSize ( 0.05 );
		h1_peak->SetLineWidth ( 2 );

		//channel_vs_meanBL->SetMarkerSize(3);
		channel_vs_meanBL->SetMarkerStyle(8);
		//channel_vs_meanBL->GetXaxis()->SetRangeUser(16,32);
		channel_vs_meanBL->GetYaxis()->SetTitle("Bl_{mean} = #Sigma _{i} B_{i} * C_{i} / #Sigma_{i} C_{i} ");
		channel_vs_meanBL->GetYaxis()->SetTitleSize(0.04);
		channel_vs_meanBL->GetXaxis()->SetTitleSize(0.04);
		//channel_vs_SD->SetMarkerSize(3);
		channel_vs_SD->SetMarkerStyle(8);
	//	channel_vs_SD->GetXaxis()->SetRangeUser(16,32);
		channel_vs_SD->GetYaxis()->SetTitle("#sigma B^{2} = #Sigma_{i} (B_{i} - B_{m})^{2} * C_{i} / #Sigma_{i} C_{i} ");
		channel_vs_SD->GetYaxis()->SetTitleSize(0.04);
		channel_vs_SD->GetXaxis()->SetTitleSize(0.04);
        
        
        TH1D* proj_sd = new TH1D("proj_sd","proj_sd",5,0,5);
        proj_sd = channel_vs_SD->ProjectionY("",16,32);

		h0->Write();
		h1->Write();
		h2->Write();
		h3->Write();
		h4->Write();
		h_count->Write();
		h1_width->Write();
		h1_peak->Write();
	//	channel_vs_position->Write();
		channel_vs_width->Write();
		channel_vs_meanBL->Write();
		channel_vs_SD->Write();
        proj_sd->Write();
		
		scan_results->Close();
	}

	return 0;
}
