#ifndef CONDCORE_SISTRIPPLUGINS_SISTRIPTKMAPS_H
#define CONDCORE_SISTRIPPLUGINS_SISTRIPTKMAPS_H

// CMSSW includes
#include "CalibTracker/StandaloneTrackerTopology/interface/StandaloneTrackerTopology.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "DataFormats/TrackerCommon/interface/PixelBarrelName.h"
#include "DataFormats/TrackerCommon/interface/PixelEndcapName.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// ROOT includes
#include "TArrow.h"
#include "TPaletteAxis.h"
#include "TGaxis.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TLatex.h"
#include "TH2Poly.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TProfile2D.h"
#include "TStyle.h"

// STL includes
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <boost/tokenizer.hpp>
#include <boost/range/adaptor/indexed.hpp>

#define MYOUT LogDebug("SiStripTkMaps")

/*--------------------------------------------------------------------
/ Ancillary class to build SiStrip Tracker maps
/--------------------------------------------------------------------*/
class SiStripTkMaps {
public:
  SiStripTkMaps(const char* option)
      : m_option{option},
        m_trackerTopo{StandaloneTrackerTopology::fromTrackerParametersXMLFile(
            edm::FileInPath("Geometry/TrackerCommonData/data/PhaseI/trackerParameters.xml").fullPath())} {}

  ~SiStripTkMaps() {}

  //============================================================================
  void bookMap(const std::string mapTitle, const std::string zAxisTitle) {
    double minx = 0xFFFFFF, maxx = -0xFFFFFF, miny = 0xFFFFFF, maxy = -0xFFFFFFF;
    readVertices(minx, maxx, miny, maxy);

    // set the titles
    m_zAxisTitle = zAxisTitle;
    m_mapTitle = mapTitle;

    TGaxis::SetMaxDigits(2);

    // margin of the box
    static constexpr int margin = 5;
    m_trackerMap =
        new TH2Poly("Tracker Map", m_mapTitle.c_str(), minx - margin, maxx + margin, miny - margin, maxy + margin);
    m_trackerMap->SetFloat();
    m_trackerMap->SetOption(m_option);
    m_trackerMap->SetStats(false);
    m_trackerMap->GetZaxis()->SetLabelSize(0.03);
    m_trackerMap->GetZaxis()->SetTitleOffset(0.5);
    m_trackerMap->GetZaxis()->SetTitleSize(0.05);
    m_trackerMap->GetZaxis()->SetTitle(m_zAxisTitle.c_str());
    m_trackerMap->GetZaxis()->CenterTitle();

    for (const auto& pair : m_bins) {
      m_trackerMap->AddBin(pair.second->Clone());
    }
  }

  //============================================================================
  void fill(long rawid, double val) {
    m_trackerMap->Fill(TString::Format("%ld", rawid), val);
    m_values.push_back(val);
  }

  //============================================================================
  void drawMap(TCanvas& canvas, std::string option = "") {
    static constexpr float tmargin_ = 0.08;
    static constexpr float bmargin_ = 0.02;
    static constexpr float lmargin_ = 0.02;
    static constexpr float rmargin_ = 0.08;

    canvas.cd();
    adjustCanvasMargins(canvas.cd(), tmargin_, bmargin_, lmargin_, rmargin_);
    canvas.Update();

    m_trackerMap->SetTitle("");
    if (!option.empty()) {
      m_trackerMap->Draw(option.c_str());
    } else {
      m_trackerMap->Draw();
    }

    canvas.SetFrameLineColor(0);
    gPad->Update();
    TPaletteAxis* palette = (TPaletteAxis*)m_trackerMap->GetListOfFunctions()->FindObject("palette");
    if (palette != nullptr) {
      palette->SetLabelSize(0.02);
      palette->SetX1NDC(1 - rmargin_);
      palette->SetX2NDC(1 - rmargin_ + lmargin_);
    }
  }

  //============================================================================
  void dressMap(TCanvas& canv) {
    static constexpr int wH_ = 3000;
    static constexpr int hH_ = 850;
    if (canv.GetWindowHeight() != hH_ && canv.GetWindowWidth() != wH_) {
      canv.SetWindowSize(wH_, hH_);
    }

    std::array<std::string, 12> barrelNames = {
        {"TIB L2", "TIB L1", "TIB L4", "TIB L3", "TOB L2", "TOB L1", "TOB L4", " TOB L3", "TOB L6", "TOB L5"}};
    std::array<std::string, 4> endcapNames = {{"TID", "TEC", "TID", "TEC"}};
    std::array<std::string, 24> disknumbering = {{"+1", "+2", "+3", "+1", "+2", "+3", "+4", "+5",
                                                  "+6", "+7", "+8", "+9", "-1", "-2", "-3", "-1",
                                                  "-2", "-3", "-4", "-5", "-6", "-7", "-8", "-9"}};

    static constexpr std::array<float, 12> b_coordx = {
        {0.1, 0.1, 0.26, 0.26, 0.41, 0.41, 0.56, 0.56, 0.725, 0.725, 0.05, 0.17}};
    static constexpr std::array<float, 12> b_coordy = {
        {0.70, 0.45, 0.70, 0.45, 0.70, 0.46, 0.70, 0.46, 0.70, 0.46, 0.85, 0.85}};

    static constexpr std::array<float, 4> e_coordx = {{0.01, 0.21, 0.01, 0.21}};
    static constexpr std::array<float, 4> e_coordy = {{0.89, 0.89, 0.17, 0.17}};

    static constexpr std::array<float, 24> n_coordx = {{0.01,  0.087, 0.165, 0.227, 0.305, 0.383, 0.461, 0.539,
                                                        0.616, 0.694, 0.772, 0.850, 0.01,  0.087, 0.165, 0.227,
                                                        0.305, 0.383, 0.461, 0.539, 0.617, 0.695, 0.773, 0.851}};

    static constexpr std::array<float, 24> n_coordy = {{0.85, 0.85, 0.85, 0.85, 0.85, 0.85, 0.85, 0.85,
                                                        0.85, 0.85, 0.85, 0.85, 0.13, 0.13, 0.13, 0.13,
                                                        0.13, 0.13, 0.13, 0.13, 0.13, 0.13, 0.13, 0.13}};

    canv.cd();
    for (const auto& name : barrelNames | boost::adaptors::indexed(0)) {
      auto ltx = TLatex();
      ltx.SetTextFont(62);
      ltx.SetTextSize(0.035);
      ltx.SetTextAlign(11);
      ltx.DrawLatexNDC(b_coordx[name.index()], b_coordy[name.index()], name.value().c_str());
    }

    for (const auto& name : endcapNames | boost::adaptors::indexed(0)) {
      auto ltx = TLatex();
      ltx.SetTextFont(62);
      ltx.SetTextSize(0.05);
      ltx.SetTextAlign(11);
      ltx.DrawLatexNDC(e_coordx[name.index()], e_coordy[name.index()], name.value().c_str());
    }

    for (const auto& name : disknumbering | boost::adaptors::indexed(0)) {
      auto ltx = TLatex();
      ltx.SetTextFont(62);
      ltx.SetTextSize(0.035);
      ltx.SetTextAlign(11);
      ltx.DrawLatexNDC(n_coordx[name.index()], n_coordy[name.index()], name.value().c_str());
    }

    auto ltx = TLatex();
    ltx.SetTextFont(62);
    ltx.SetTextSize(0.045);
    ltx.SetTextAlign(11);
    ltx.DrawLatexNDC(gPad->GetLeftMargin(), 1 - gPad->GetTopMargin() + 0.03, m_mapTitle.c_str());

    // barrel axes
    float phiX1 = 0.09;
    float phiX2 = 0.23;
    float phiY1 = 0.24;
    float phiY2 = phiY1;
    auto arrowPhi = TArrow();
    arrowPhi.SetLineColor(kBlue);
    arrowPhi.SetLineWidth(2);
    arrowPhi.SetOption("|>");
    arrowPhi.SetArrowSize(10);
    arrowPhi.DrawLineNDC(phiX1, phiY1, phiX2, phiY2);
    //arrowPhi.SetNDC(kTRUE);
    //arrowPhi.DrawArrow(0.09,0.25,0.23,0.25,0.1,"|>");

    float zX1 = phiX2;
    float zX2 = zX1;
    float zY1 = phiY1;
    float zY2 = 0.45;
    auto arrowZ = TArrow();
    arrowZ.SetLineColor(kBlue);
    arrowZ.SetLineWidth(2);
    arrowZ.SetOption("|>");
    arrowZ.SetArrowSize(10);
    arrowZ.DrawLineNDC(zX1, zY1, zX2, zY2);
    //arrowZ.SetNDC(kTRUE);
    //arrowZ.DrawArrow(0.09,0.25,0.09,0.45,0.1,"|>");

    auto textPhi = TLatex();
    textPhi.SetTextSize(0.04);
    textPhi.SetTextAlign(11);
    textPhi.SetTextColor(kBlue);
    textPhi.DrawLatexNDC(phiX1 - 0.01, phiY1 + 0.01, "#phi");

    auto textZ = TLatex();
    textZ.SetTextSize(0.04);
    textZ.SetTextAlign(11);
    textZ.SetTextColor(kBlue);
    textZ.DrawLatexNDC(zX1 + 0.005, zY2, "z");

    // endcap axes
    float xX1 = 0.85;
    float xX2 = 0.89;
    float xY1 = 0.83;
    float xY2 = xY1;
    auto arrowX = TArrow();
    arrowX.SetLineColor(kBlue);
    arrowX.SetLineWidth(2);
    arrowX.SetOption("|>");
    arrowX.SetArrowSize(10);
    arrowX.DrawLineNDC(xX1, xY1, xX2, xY2);
    //arrowX.SetNDC(kTRUE);
    //arrowX.DrawArrow(0.09,0.25,0.23,0.25,0.1,"|>");

    float yX1 = xX2;
    float yX2 = yX1;
    float yY1 = xY1;
    float yY2 = 0.95;
    auto arrowY = TArrow();
    arrowY.SetLineColor(kBlue);
    arrowY.SetLineWidth(2);
    arrowY.SetOption("|>");
    arrowY.SetArrowSize(10);
    arrowY.DrawLineNDC(yX1, yY1, yX2, yY2);
    //arrowY.SetNDC(kTRUE);
    //arrowY.DrawArrow(0.09,0.25,0.09,0.45,0.1,"|>");

    auto textX = TLatex();
    textX.SetTextSize(0.04);
    textX.SetTextAlign(11);
    textX.SetTextColor(kBlue);
    textX.DrawLatexNDC(xX1, xY1 - 0.03, "x");

    auto textY = TLatex();
    textY.SetTextSize(0.04);
    textY.SetTextAlign(11);
    textY.SetTextColor(kBlue);
    textY.DrawLatexNDC(yX1 - 0.01, yY2, "y");

    canv.Modified();
    canv.Update();  // make sure it's really (re)drawn
  }

  //============================================================================
  const TH2Poly* getTheMap() { return m_trackerMap; }

  //============================================================================
  inline const std::string& getTheMapTitle() { return m_mapTitle; }

  //============================================================================
  inline const std::string& getTheZAxisTitle() { return m_zAxisTitle; }

  //============================================================================
  inline const std::vector<unsigned int>& getTheFilledIds() { return m_detIdVector; }

  //============================================================================
  inline const std::vector<double>& getTheFilledValues() { return m_values; }

  //============================================================================
  void setZAxisRange(double xmin, double xmax) { m_trackerMap->GetZaxis()->SetRangeUser(xmin, xmax); }

  //============================================================================
  void adjustCanvasMargins(TVirtualPad* pad, float top, float bottom, float left, float right) {
    if (top > 0) {
      pad->SetTopMargin(top);
    }
    if (bottom > 0) {
      pad->SetBottomMargin(bottom);
    }
    if (left > 0) {
      pad->SetLeftMargin(left);
    }
    if (right > 0) {
      pad->SetRightMargin(right);
    }
  }

  //============================================================================
  void readVertices(double& minx, double& maxx, double& miny, double& maxy) {
    std::ifstream in;

    in.open(edm::FileInPath("DQM/SiStripMonitorClient/data/Geometry/tracker_map_bare").fullPath().c_str());

    if (!in.good()) {
      throw cms::Exception("FileError") << "Problem opening corner file!!" << std::endl;
      return;
    }

    while (in.good()) {
      long detid = 0;
      double x[5], y[5];

      std::string line;
      std::getline(in, line);

      TString string(line);
      TObjArray* array = string.Tokenize(" ");
      int ix{0}, iy{0};
      bool isPixel{false};
      for (int i = 0; i < array->GetEntries(); ++i) {
        if (i == 0) {
          detid = static_cast<TObjString*>(array->At(i))->String().Atoll();

          // Drop Pixel Data
          DetId detId(detid);
          if (detId.subdetId() == PixelSubdetector::PixelBarrel || detId.subdetId() == PixelSubdetector::PixelEndcap) {
            isPixel = true;
            break;
          }
        } else {
          if (i % 2 == 0) {
            x[ix] = static_cast<TObjString*>(array->At(i))->String().Atof();
            if (x[ix] < minx) {
              minx = x[ix];
            }
            if (x[ix] > maxx) {
              maxx = x[ix];
            }
            ++ix;
          } else {
            y[iy] = static_cast<TObjString*>(array->At(i))->String().Atof();
            if (y[iy] < miny) {
              miny = y[iy];
            }
            if (y[iy] > maxy) {
              maxy = y[iy];
            }
            ++iy;
          }  // else
        }    // else
      }      // loop on entries

      if (isPixel) {
        continue;
      }

      m_bins[detid] = std::make_shared<TGraph>(ix, x, y);
      m_bins[detid]->SetName(TString::Format("%ld", detid));
      m_bins[detid]->SetTitle(TString::Format("Module ID=%ld", detid));
      m_detIdVector.push_back(detid);
    }
  }

private:
  Option_t* m_option;
  std::string m_mapTitle = "";
  std::string m_zAxisTitle = "";
  double m_axmin, m_axmax;
  std::map<long, std::shared_ptr<TGraph> > m_bins;
  std::vector<unsigned int> m_detIdVector;
  std::vector<double> m_values;
  TrackerTopology m_trackerTopo;
  TH2Poly* m_trackerMap{nullptr};
};

#endif
