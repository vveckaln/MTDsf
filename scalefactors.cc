#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TROOT.h"
#include "TStyle.h"
#include <string>
#include "tdrstyle.C"


void save(TH2 * hist, const char * tag) 
{
    TLatex latex;
    latex.SetNDC();
    latex.SetTextAngle(0);
    latex.SetTextColor(kBlack);
    latex.SetTextFont(42);
    latex.SetTextAlign(31);
    latex.SetTextSize(0.03);

    TCanvas * can = new TCanvas;
    hist -> Draw("COLZ");
    latex.DrawLatex(0.6, 0.93, tag);
    const string dir = "sfpics";
    can -> SaveAs((dir + "/" + hist -> GetName() + ".png").c_str());
    can -> SaveAs((dir + "/" + hist -> GetName() + ".pdf").c_str());
    delete can;
}




// TH2F * prepareHisto(std::string nameHisto, int nbin_Pt, float * bins_Pt, float * sf_Pt, int nbin_Eta, float * bins_Eta, float *sf_Eta) 
// {
//     TH2F * eff = new TH2F(nameHisto.c_str(), "", nbin_Pt, bins_Pt, nbin_Eta, bins_Eta);
//     for(unsigned char ptind = 0; ptin <= nbin_Pt + 1; ptind ++) 
//       { 
//         for(unsinged char etaind = 0; etaind <= nbin_Eta + 1; etaind ++) 
// 	  { 
//             eff -> SetBinContent(ptind, etaind, sf_Pt[ptind]*sf_Eta[etaind]); 
//         }
//     }
//     return eff;
// }





//int run(std::string name) {
int scalefactors(const char * name)
{
  enum  MTDsections     {BARREL = 0, ENDCAP = 1}; 

  gROOT -> SetBatch(kTRUE);
  gStyle -> SetOptStat(0);
  TFile * file = TFile::Open(name, "update");
    
  {
    const unsigned char   nbin_Pt                                                 = 1;
    const float           bins_Pt[nbin_Pt + 1]                                    = {20.0, 200.0};
    const unsigned char   nbin_Eta                                                = 3;
    const float           bins_Eta[nbin_Eta + 1]                                  = {-3.5, -1.5, 1.5, 3.5};

    const float           SF_eff_Pt_BTagSFLoose_barrel[nbin_Pt ]                  = {1.005};
    const float           SF_eff_Pt_BTagSFLoose_endcap[nbin_Pt]                   = {1.011};
    
    const float           SF_frate_Pt_BTagSFLoose_barrel[nbin_Pt]                 = {1.0};
    const float           SF_frate_Pt_BTagSFLoose_endcap[nbin_Pt]                 = {1.0};


    const float           SF_eff_Pt_BTagSFMedium_barrel[nbin_Pt]                  = {1.023};
    const float           SF_eff_Pt_BTagSFMedium_endcap[nbin_Pt]                  = {1.048};


    const float           SF_frate_Pt_BTagSFMedium_barrel[nbin_Pt]                = {1.0};
    const float           SF_frate_Pt_BTagSFMedium_endcap[nbin_Pt]                = {1.0};

    const float           SF_eff_Pt_BTagSFTight_barrel[nbin_Pt]                   = {1.019};
    const float           SF_eff_Pt_BTagSFTight_endcap[nbin_Pt]                   = {1.090};
    
    const float           SF_frate_Pt_BTagSFTight_barrel[nbin_Pt]                 = {1.0};
    const float           SF_frate_Pt_BTagSFTight_endcap[nbin_Pt]                 = {1.0};

    const float           SF_eff_Pt_TauSF_barrel[nbin_Pt]                         = {1.05};
    const float           SF_eff_Pt_TauSF_endcap[nbin_Pt]                         = {1.093};

    const float           SF_frate_Pt_TauSF_barrel[nbin_Pt]                       = {1.0};
    const float           SF_frate_Pt_TauSF_endcap[nbin_Pt]                       = {1.0};

    const float           SF_eff_Pt_EleIDSF_barrel[nbin_Pt]                       = {1.033};
    const float           SF_eff_Pt_EleIDSF_endcap[nbin_Pt]                       = {1.000};

    const float           SF_frate_Pt_EleIDSF_barrel[nbin_Pt]                     = {1.0};
    const float           SF_frate_Pt_EleIDSF_endcap[nbin_Pt]                     = {1.0};
      
    const unsigned char   nhistos                                                 = 10;
    const float         * values[nhistos][2] =
      {
	{
	  SF_eff_Pt_BTagSFLoose_barrel,
	  SF_eff_Pt_BTagSFLoose_endcap
	},
	{
    
	  SF_frate_Pt_BTagSFLoose_barrel,
	  SF_frate_Pt_BTagSFLoose_endcap
	},
	{
	  SF_eff_Pt_BTagSFMedium_barrel,
	  SF_eff_Pt_BTagSFMedium_endcap
	},
	{
	  SF_frate_Pt_BTagSFMedium_barrel,
	  SF_frate_Pt_BTagSFMedium_endcap
	},
	{
	  SF_eff_Pt_BTagSFTight_barrel,
	  SF_eff_Pt_BTagSFTight_endcap
	},
	{
	  SF_frate_Pt_BTagSFTight_barrel,
	  SF_frate_Pt_BTagSFTight_endcap
	},
	{

	  SF_eff_Pt_TauSF_barrel,
	  SF_eff_Pt_TauSF_endcap           
	},
	{
	  SF_frate_Pt_TauSF_barrel,
	  SF_frate_Pt_TauSF_endcap
	},
	{
	  SF_eff_Pt_EleIDSF_barrel,
	  SF_eff_Pt_EleIDSF_endcap
	},
	{
	  SF_frate_Pt_EleIDSF_barrel,
	  SF_frate_Pt_EleIDSF_endcap
	}
      };
  
    const char * hnames[nhistos] = 
      {
	"eff_BTagSFLoose",
	"frate_BTagSFLoose",
	"eff_BTagSFMedium",
	"frate_BTagSFMedium",
	"eff_BTagSFTight",
	"frate_BTagSFTight",
	"eff_TauSF",
	"frate_TauSF",
	"eff_EleIDSF",
	"frate_EleIDSF"
      };

    const char * htitles[nhistos] =
      {
	"b-tagging efficiency loose Scale Factor",
	"b-tagging fakerate loose Scale Factor",
	"b-tagging efficiency medium Scale Factor"
	"b-tagging fakerate medium Scale Factor",
	"b-tagging efficiency tight Scale Factor",
	"b-tagging fakerate tight Scale Factor",
	"Tau efficiency Scale Factor",
	"Tau fakerate Scale Factor",
	"Electron identification efficiency Scale Factor",
	"Electron identification fakerate Scale Factor"
      };


    for (unsigned char hind = 0; hind < nhistos; hind ++)
      {
	TH2 * hist = new TH2F(hnames[hind], ";p_{T} [GeV]; #eta", nbin_Pt, bins_Pt, nbin_Eta, bins_Eta);
	for (unsigned char ptind = 0; ptind < nbin_Pt; ptind ++)
	  {
	    for (unsigned char etaind = 0; etaind < nbin_Eta; etaind ++)
	      {
		const unsigned char MTSECT = fabs(bins_Eta[etaind]) <= 1.5 and fabs(bins_Eta[etaind + 1]) <= 1.5 ? BARREL : ENDCAP;
		hist -> Fill(bins_Pt[ptind],  bins_Eta[etaind], values[hind][MTSECT][0]);
		//	printf("pt %u %f  et %u %f hind %u MTSECT %u %f\n", ptind, bins_Pt[ptind], etaind, bins_Eta[etaind], hind, MTSECT, values[hind][MTSECT][0]);
	      }
	  }
	//	getchar();
	hist -> GetZaxis() -> SetRangeUser(1.0, 1.07);
	save(hist, htitles[hind]);
	hist -> Write();
      }

  }
 
  //=================================================================================================//
  {
    const unsigned char   nbin_Pt                                                   = 3;
    const float           bins_Pt[nbin_Pt + 1]                                          = {10.0, 20.0, 40.0, 200.0};
    const unsigned char   nbin_Eta                                                  = 3;
    const float           bins_Eta[nbin_Eta +1]                                     = {-2.8, -1.5, 1.5, 2.8};

    const unsigned char nbin_Pt_pr = 5;
    const float           SF_eff_Pt_EleIsoSF_barrel[nbin_Pt_pr]                        = {1.019, 1.019, 1.019, 1.026, 1.026};
    const float           SF_eff_Pt_EleIsoSF_endcap[nbin_Pt_pr]                        = {1.033, 1.033, 1.035, 1.045, 1.045};
    const float           SF_frate_Pt_EleIsoSF_barrel[nbin_Pt_pr]                      = {1.0,   1.0,   1.0,   1.0,   1.0};
    const float           SF_frate_Pt_EleIsoSF_endcap[nbin_Pt_pr]                      = {1.0,   1.0,   1.0,   1.0,   1.0};

    const float           SF_eff_Pt_MuonIsoSF_barrel[nbin_Pt_pr]                       = {1.019, 1.019, 1.019, 1.026, 1.026};
    const float           SF_eff_Pt_MuonIsoSF_endcap[nbin_Pt_pr]                       = {1.033, 1.033, 1.035, 1.045, 1.045};
    const float           SF_frate_Pt_MuonIsoSF_barrel[nbin_Pt_pr]                     = {1.0,   1.0,   1.0,   1.0,   1.0};
    const float           SF_frate_Pt_MuonIsoSF_endcap[nbin_Pt_pr]                     = {1.0,   1.0,   1.0,   1.0,   1.0};

    const float           SF_eff_Pt_PhotonIsoSF_barrel[nbin_Pt_pr]                     = {1.019, 1.019, 1.019, 1.026, 1.026};
    const float           SF_eff_Pt_PhotonIsoSF_endcap[nbin_Pt_pr]                     = {1.033, 1.033, 1.035, 1.045, 1.045};
    const float           SF_frate_Pt_PhotonIsoSF_barrel[nbin_Pt_pr]                   = {1.0,   1.0,   1.0,   1.0,   1.0};
    const float           SF_frate_Pt_PhotonIsoSF_endcap[nbin_Pt_pr]                   = {1.0,   1.0,   1.0,   1.0,   1.0};

    const unsigned char   nhistos                                                   = 6;

    const float         * values[nhistos][2] =
      {
	{
	  SF_eff_Pt_EleIsoSF_barrel,
	  SF_eff_Pt_EleIsoSF_endcap
	},
	{
	  SF_frate_Pt_EleIsoSF_barrel,
	  SF_frate_Pt_EleIsoSF_endcap
	},
	{
	  SF_eff_Pt_MuonIsoSF_barrel,
	  SF_eff_Pt_MuonIsoSF_endcap
	},
	{
	  SF_frate_Pt_MuonIsoSF_barrel,
	  SF_frate_Pt_MuonIsoSF_endcap
	},
	{
	  SF_eff_Pt_PhotonIsoSF_barrel,
	  SF_eff_Pt_PhotonIsoSF_endcap
	},
	{
	  SF_frate_Pt_PhotonIsoSF_barrel,
	  SF_frate_Pt_PhotonIsoSF_endcap
	}
      };

    const char * hnames[nhistos] = 
      {
	"eff_EleIsoSF",
	"frate_EleIsoSF",
	"eff_MuonIsoSF",
	"frate_MuonIsoSF",
	"eff_PhotonIsoSF",
	"frate_PhotonIsoSF"
      };

    const char * htitles[nhistos] = 
      {
	"Electron isolation efficiency Scale Factor",
	"Electron isolation fakerate Scale Factor",
	"Muon isolation efficiency Scale Factor",
	"Muon isolation fakerate Scale Factor",
	"Photon isolation efficiency Scale Factor",
	"Photon isolation fakerate Scale Factor"
      };

    for (unsigned char hind = 0; hind < nhistos; hind ++)
      {
	TH2 * hist = new TH2F(hnames[hind], ";p_{T} [GeV]; #eta", nbin_Pt, bins_Pt, nbin_Eta, bins_Eta);
	for (unsigned char ptind = 0; ptind < nbin_Pt; ptind ++)
	  {
	    for (unsigned char etaind = 0; etaind < nbin_Eta; etaind ++)
	      {
		const unsigned char MTSECT = fabs(bins_Eta[etaind]) <= 1.5 and fabs(bins_Eta[etaind + 1]) <= 1.5 ? BARREL : ENDCAP;
		hist -> Fill(bins_Pt[ptind],  bins_Eta[etaind], values[hind][MTSECT][ptind]);
	      }
	  }
	hist -> GetZaxis() -> SetRangeUser(1.0, 1.07);
	save(hist, htitles[hind]);
	hist -> Write();
      }

  }
  
  file -> Close();

  return 0;
}




