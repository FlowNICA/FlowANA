struct coord
{
    float xlow;
    float ylow;
    float xup;
    float yup;
};

void ReadFlow(const char *inFileName)
{
    // Setting up global variables for the plot
    gROOT->SetStyle("Pub");
    gROOT->ForceStyle();
    gStyle->SetPalette(kDarkRainBow);
    gStyle->SetErrorX(0);
    gStyle->SetTitleSize(0.05);

    const double binpt [11] = {0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.8,2.3,2.8,4.0};
    const TString hName [4] = {TString("h^{#pm}"),TString("#pi^{#pm}"),TString("K^{#pm}"),TString("p+#bar{p}")};
    const TString hCent [6] = {TString("0-10%"),TString("10-20%"),TString("20-30%"),TString("30-40%"),TString("40-50%"),TString("50-60%")};
    const std::pair<Double_t,Double_t> v2range[6][4] = {{{-0.023,0.06},{-0.023,0.06},{-0.023,0.06},{-0.023,0.06}},
                                                        {{-0.005,0.09},{-0.005,0.09},{-0.005,0.09},{-0.005,0.09}},
                                                        {{-0.023,0.09},{-0.023,0.09},{-0.023,0.09},{-0.023,0.09}},
                                                        {{-0.051,0.16},{-0.051,0.16},{-0.051,0.16},{-0.051,0.16}},
                                                        {{-0.023,0.16},{-0.023,0.16},{-0.023,0.16},{-0.023,0.16}},
                                                        {{-0.023,0.16},{-0.023,0.16},{-0.023,0.16},{-0.023,0.16}}};

    TFile *fi = new TFile(inFileName,"read");
    TH1F *hv2[4][6][10][4]; //idet-icent-ipt-ipid
    TH1F *hv22[4][10][4];
    TH1F *hpt[6][10][4];

    for (int idet=0;idet<4;idet++)
    {
        for (int icent=0;icent<6;icent++)
        {
            for (int ipt=0; ipt<10; ipt++)
            {
                for (int ipid=0;ipid<4;ipid++)
                {
                    if (idet == 0)
                    {
                        hpt[icent][ipt][ipid] = (TH1F*) fi->Get(Form("hpt_%i_%i_%i",icent,ipt,ipid));
                    }
                    if (icent==0)
                    {
                        hv22[idet][ipt][ipid] = (TH1F*) fi->Get(Form("hv22_%i_%i_%i",idet,ipt,ipid));
                    }
                    hv2[idet][icent][ipt][ipid] = (TH1F*) fi->Get(Form("hv2_%i_%i_%i_%i",idet,icent,ipt,ipid));
                }
            }
        }
    }
    std::vector<Double_t> vPt, vV2tpc, vV2fhcal;
    std::vector<Double_t> ePt, eV2tpc, eV2fhcal;
    int det1 = 0, det2 = 3;
    TGraphErrors *grv2tpc[6][4], *grv2fhcal[6][4]; //icent-ipid

    for (int icent=0;icent<6;icent++)
    {
        for (int ipid=0;ipid<4;ipid++)
        {
            for (int ipt=0;ipt<10;ipt++)
            {
                vPt.push_back(hpt[icent][ipt][ipid]->GetMean());
                ePt.push_back(0.);
                vV2tpc.push_back(hv2[det1][icent][ipt][ipid]->GetMean());
                eV2tpc.push_back(hv2[det1][icent][ipt][ipid]->GetMeanError());
                vV2fhcal.push_back(hv2[det2][icent][ipt][ipid]->GetMean());
                eV2fhcal.push_back(hv2[det2][icent][ipt][ipid]->GetMeanError());
            }
            grv2tpc[icent][ipid] = new TGraphErrors(vPt.size(),&vPt[0],&vV2tpc[0],&ePt[0],&eV2tpc[0]);
            grv2fhcal[icent][ipid] = new TGraphErrors(vPt.size(),&vPt[0],&vV2fhcal[0],&ePt[0],&eV2fhcal[0]);

            grv2tpc[icent][ipid]->GetXaxis()->SetTitle("p_{T}, [GeV/c]");
            grv2tpc[icent][ipid]->GetYaxis()->SetTitle("v_{2}");
            grv2fhcal[icent][ipid]->GetXaxis()->SetTitle("p_{T}, [GeV/c]");
            grv2fhcal[icent][ipid]->GetYaxis()->SetTitle("v_{2}");

            grv2tpc[icent][ipid]  ->SetMarkerStyle(20);
            grv2fhcal[icent][ipid]->SetMarkerStyle(26);
            
            grv2tpc[icent][ipid]  ->GetYaxis()->SetRangeUser(v2range[icent][ipid].first,v2range[icent][ipid].second);
            grv2fhcal[icent][ipid]->GetYaxis()->SetRangeUser(v2range[icent][ipid].first,v2range[icent][ipid].second);

            grv2tpc[icent][ipid]  ->GetXaxis()->SetLimits(0.,3.49);
            grv2fhcal[icent][ipid]->GetXaxis()->SetLimits(0.,3.49);

            grv2tpc[icent][ipid]  ->SetMarkerSize(1.6);
            grv2fhcal[icent][ipid]->SetMarkerSize(1.6);

            vPt.clear(); ePt.clear();
            vV2tpc.clear(); eV2tpc.clear();
            vV2fhcal.clear(); eV2fhcal.clear();
        }
    }

    TH1F *hptTmp;
    TGraphErrors *grv22tpc[4], *grv22fhcal[4];
    for (int ipid=0;ipid<4;ipid++) {
        for (int ipt = 0; ipt < 10; ipt++) {
            hptTmp = (TH1F*) hpt[1][ipt][ipid]->Clone();
            hptTmp->Add(hpt[2][ipt][ipid]);
            hptTmp->Add(hpt[3][ipt][ipid]);
            vPt.push_back(hptTmp->GetMean());
            ePt.push_back(0.);
            vV2tpc.push_back(hv22[det1][ipt][ipid]->GetMean());
            eV2tpc.push_back(hv22[det1][ipt][ipid]->GetMeanError());
            vV2fhcal.push_back(hv22[det2][ipt][ipid]->GetMean());
            eV2fhcal.push_back(hv22[det2][ipt][ipid]->GetMeanError());
            delete hptTmp;
        }
        grv22tpc[ipid] = new TGraphErrors(vPt.size(),&vPt[0],&vV2tpc[0],&ePt[0],&eV2tpc[0]);
        grv22fhcal[ipid] = new TGraphErrors(vPt.size(),&vPt[0],&vV2fhcal[0],&ePt[0],&eV2fhcal[0]);

        grv22tpc[ipid]->GetXaxis()->SetTitle("p_{T}, [GeV/c]");
        grv22tpc[ipid]->GetYaxis()->SetTitle("v_{2}");
        grv22fhcal[ipid]->GetXaxis()->SetTitle("p_{T}, [GeV/c]");
        grv22fhcal[ipid]->GetYaxis()->SetTitle("v_{2}");

        grv22tpc[ipid]  ->SetMarkerStyle(20);
        grv22fhcal[ipid]->SetMarkerStyle(26);

        grv22tpc[ipid]  ->GetYaxis()->SetRangeUser(v2range[1][ipid].first,v2range[1][ipid].second);
        grv22fhcal[ipid]->GetYaxis()->SetRangeUser(v2range[1][ipid].first,v2range[1][ipid].second);

        grv22tpc[ipid]  ->GetXaxis()->SetLimits(0.,3.49);
        grv22fhcal[ipid]->GetXaxis()->SetLimits(0.,3.49);

        grv22tpc[ipid]  ->SetMarkerSize(1.6);
        grv22fhcal[ipid]->SetMarkerSize(1.6);

        vPt.clear(); ePt.clear();
        vV2tpc.clear(); eV2tpc.clear();
        vV2fhcal.clear(); eV2fhcal.clear();
    }

    TCanvas *canv[6];
    TLegend *leg[6][4];
    TLine   *line0 = new TLine();
    coord lc[4] = {{0.2,0.7,0.55,0.89},{0.2,0.82,0.55,0.89},{0.2,0.82,0.55,0.89},{0.2,0.82,0.55,0.89}};
    line0->SetLineStyle(2);
    line0->SetLineColor(4);
    line0->SetLineWidth(1);
    for (int icent=0; icent<6; icent++)
    {
        canv[icent] = new TCanvas(Form("canv_%i",icent),Form("canv cent %i",icent),900,900);
        canv[icent] ->Divide(2,2);
    }
    for (int icent=0; icent<6; icent++)
    {
        for (int ipid=0;ipid<4;ipid++)
        {
            canv[icent]->cd(ipid+1);
            grv2tpc[icent][ipid]->Draw("AP PLC PMC");
            grv2fhcal[icent][ipid]->Draw("P PLC PMC");
            leg[icent][ipid] = new TLegend(lc[ipid].xlow,lc[ipid].ylow,lc[ipid].xup,lc[ipid].yup);
            leg[icent][ipid]->SetBorderSize(0.);
            leg[icent][ipid]->SetHeader(Form("%s, %s",hCent[icent].Data(),hName[ipid].Data()),"C");
            if (ipid == 0)
            {
                leg[icent][ipid]->AddEntry(grv2tpc[icent][ipid],"v2{TPC}","p");
                leg[icent][ipid]->AddEntry(grv2fhcal[icent][ipid],"v2{FHCal}","p");
            }
            leg[icent][ipid]->Draw();

            line0->DrawLine(0.,0.,3.49,0.);
        }
        // canv[icent]->SaveAs(Form("/home/peter/Documents/WorkLocal/STAR/Pics/Models/EPflow/v2_cent%i.png",icent));
            if (icent == 0)
                canv[icent]->Print(Form("/home/peter/Documents/WorkLocal/STAR/Pics/Models/EPflow/v2.pdf("),"pdf");
            if (icent == 5)
                canv[icent]->Print(Form("/home/peter/Documents/WorkLocal/STAR/Pics/Models/EPflow/v2.pdf)"),"pdf");
            if (icent != 0 && icent != 5)
                canv[icent]->Print(Form("/home/peter/Documents/WorkLocal/STAR/Pics/Models/EPflow/v2.pdf"),"pdf");
    }

    TCanvas *canv2 = new TCanvas("canv2","10-40%", 900, 900);
    TLegend *leg2[4];
    canv2->Divide(2,2);
    for (int ipid=0;ipid<4;ipid++) {
        canv2->cd(ipid + 1);
        grv22tpc[ipid]->Draw("AP PLC PMC");
        grv22fhcal[ipid]->Draw("P PLC PMC");

        leg2[ipid] = new TLegend(lc[ipid].xlow,lc[ipid].ylow,lc[ipid].xup,lc[ipid].yup);
        leg2[ipid]->SetBorderSize(0.);
        leg2[ipid]->SetHeader(Form("10-40%, %s",hName[ipid].Data()),"C");
        if (ipid == 0)
        {
            leg2[ipid]->AddEntry(grv22tpc[ipid],"v2{TPC}","p");
            leg2[ipid]->AddEntry(grv22fhcal[ipid],"v2{FHCal}","p");
        }
        leg2[ipid]->Draw();

        line0->DrawLine(0.,0.,3.49,0.);
    }
    canv2->SaveAs("/home/peter/Documents/WorkLocal/STAR/Pics/Models/EPflow/v2_cent1040.pdf");
    
}