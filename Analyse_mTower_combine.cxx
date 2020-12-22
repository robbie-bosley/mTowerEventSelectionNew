#include "TString.h"
#include "TList.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TROOT.h"
#include "TChain.h"
#include "TMath.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TGaxis.h"
#include "TF1.h"
#include "TF2.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TVector3.h"
#include "TRandom3.h"
#include "TParameter.h"


#include "classes/mTowerHit.h"
#include "classes/mTowerClusterRobbie.h"
#include "classes/mTowerEvent.h"
#include "classes/mTowerChipRobbie.h"

//#include "/eos/project/m/mtower/public/hiroki/11_id_pileup/fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequence.hh"

#include "/eos/project/m/mtower/public/hiroki/00_util/util.hh"
#include "/eos/project/m/mtower/public/hiroki/00_util/mTowerUtil.h"
#include "/eos/project/m/mtower/public/hiroki/00_util/load_hot_pixels.h"
#include "/eos/project/m/mtower/public/hiroki/00_util/load_module_thickness.h"
#include "/eos/project/m/mtower/public/hiroki/00_util/load_inclination_parameters.h"
#include "/eos/project/m/mtower/public/hiroki/00_util/load_alignment_parameters.h"

#include "/eos/project/m/mtower/public/hiroki/11_id_pileup/cellinfo.h"

#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <tuple>
#include <fstream>
#include <string>

using namespace std;
using namespace fastjet;
using namespace ROOT::Math;

//conversion tables
const std::map< Int_t, Int_t > chipid2lanerobbie_lut = {
  { 0,40},{ 1,39},{ 2,42},{ 3,41},{ 4,44},{ 5,43},{ 6,46},{ 7,45},
  { 8,48},{ 9,47},{10,50},{11,49},{12,52},{13,51},{14,54},{15,53},
  {16,38},{17,55},{18,36},{19,37},{20,32},{21,35},{22,34},{23,33},
  {24,64},{25,63},{26,66},{27,65},{28,68},{29,67},{30,70},{31,69},
  {32,72},{33,71},{34,74},{35,73},{36,76},{37,75},{38,78},{39,77},
  {40,62},{41,79},{42,60},{43,61},{44,56},{45,59},{46,58},{47,57}
};

const std::map< Int_t, Int_t > lane2chipidrobbie_lut = {
  {40, 0},{39, 1},{42, 2},{41, 3},{44, 4},{43, 5},{46, 6},{45, 7},
  {48, 8},{47, 9},{50,10},{49,11},{52,12},{51,13},{54,14},{53,15},
  {38,16},{55,17},{36,18},{37,19},{32,20},{35,21},{34,22},{33,23},
  {64,24},{63,25},{66,26},{65,27},{68,28},{67,29},{70,30},{69,31},
  {72,32},{71,33},{74,34},{73,35},{76,36},{75,37},{78,38},{77,39},
  {62,40},{79,41},{60,42},{61,43},{56,44},{59,45},{58,46},{57,47}
};

const std::map< Int_t, Int_t > chipid2layerrobbie_lut = {
  { 0,22},{ 1,22},{ 2,20},{ 3,20},{ 4,18},{ 5,18},{ 6,16},{ 7,16},
  { 8,14},{ 9,14},{10,12},{11,12},{12,10},{13,10},{14, 8},{15, 8},
  {16, 6},{17, 6},{18, 4},{19, 4},{20, 0},{21, 0},{22, 2},{23, 2},
  {24,23},{25,23},{26,21},{27,21},{28,19},{29,19},{30,17},{31,17},
  {32,15},{33,15},{34,13},{35,13},{36,11},{37,11},{38, 9},{39, 9},
  {40, 7},{41, 7},{42, 5},{43, 5},{44, 1},{45, 1},{46, 3},{47, 3}
};

const std::map< Int_t, Int_t > lane2layerrobbie_lut = {
  {40,22},{39,22},{42,20},{41,20},{44,18},{43,18},{46,16},{45,16},
  {48,14},{47,14},{50,12},{49,12},{52,10},{51,10},{54, 8},{53, 8},
  {38, 6},{55, 6},{36, 4},{37, 4},{32, 0},{35, 0},{34, 2},{33, 2},
  {64,23},{63,23},{66,21},{65,21},{68,19},{67,19},{70,17},{69,17},
  {72,15},{71,15},{74,13},{73,13},{76,11},{75,11},{78, 9},{77, 9},
  {62, 7},{79, 7},{60, 5},{61, 5},{56, 1},{59, 1},{58, 3},{57, 3}
};

const std::map<int,bool> layer2isInvrobbie_lut = {
  { 0, kFALSE}, { 1, kTRUE}, { 2, kFALSE}, { 3, kTRUE}, 
  { 4, kFALSE}, { 5, kTRUE}, { 6, kFALSE}, { 7, kTRUE}, 
  { 8, kFALSE}, { 9, kTRUE}, {10, kFALSE}, {11, kTRUE}, 
  {12, kFALSE}, {13, kTRUE}, {14, kFALSE}, {15, kTRUE}, 
  {16, kFALSE}, {17, kTRUE}, {18, kFALSE}, {19, kTRUE}, 
  {20, kFALSE}, {21, kTRUE}, {22, kFALSE}, {23, kTRUE} 
};

bool IsLeftChip(int lane){
  int layerNr = lane2layerrobbie_lut.at(lane);
  bool isInv  = layer2isInvrobbie_lut.at(layerNr);
  int chipid = lane2chipidrobbie_lut.at(lane); 
  bool isOdd;
  if (chipid%2 == 1){ isOdd = kTRUE;}
  else {isOdd = kFALSE;}
  bool isLeft = (bool)(isOdd != isInv);
  return isLeft;
}

//for hit maps
int lane2padhitmap(int lane){
  int layerNr = lane2layerrobbie_lut.at(lane);
  bool isLeft = IsLeftChip(lane);
  
  int padid;
  if (layerNr < 6) {padid = layerNr+1;}
  else if (layerNr < 12) {padid = layerNr+7;}
  else if (layerNr < 18) {padid = layerNr+13;}
  else {padid = layerNr+19;}
  if (isLeft) {padid += 6;}
 
  return padid; 
}

//for hit maps
int lane2padoccupancy(int lane){
  int layerNr = lane2layerrobbie_lut.at(lane);
  int padid = layerNr+1;
  return padid; 
}

void Analyse_mTower(int run)
{
  //USER:HIROKI
  Double_t    energy_in   = 5.0; 
  Double_t  wBinBulk      = 0.50; //[mm]                                                                                                                                                                
  Int_t     nLayerBulk    = 24;
  Int_t     minNlayerBulk = 3;
  Double_t  W0Bulk        = 2.;
  Double_t  thDistBulk    = 5.;
  
  Int_t     ith           = 0;
  Int_t     ndiv          = 1;

  std::ios oldState(nullptr);
  oldState.copyfmt( std::cout );
  TH1::StatOverflows(kTRUE);
  //gStyle->SetOptTitle(0);                                                                                                                                                                         
  //gStyle->SetOptStat(0);                                                                                                                                                                              
  //gStyle->SetOptFit(0);                                                                                                                                                                              
  gStyle->SetPalette(kRainBow);
  gStyle->SetLineScalePS(1.0);

  //USER:ROBBIE
  //-------------------------------------------------------
  //macro to read the mTower TB data
  //using event classes: mTowerEvent, mTowerHit, mTowerChipRobbie, mTowerClusterRobbie
  //detailed code description can be found in the bachelor thesis of Aart van Bochove
  //authors: N. van der Kolk, A. van Bochove, R Bosley
  //-------------------------------------------------------
  
  //-------------------------------------------------------
  //Change this few variables to what applies for your analysis
  //-------------------------------------------------------
  
  TString fileLocationOutputTree = "./";
  TString fileLocation = "/eos/project/m/mtower/Data/Data_TB_February_2020/mTower_Data_DESY_Feb_2020_raw1/"; //The location of the selected data
  TString maskingFileLocation = "/eos/project/m/mtower/public/analysis_fw_sample_HY/data/hotpixel_TB/"; //The location of the masking .txt files

  //USER:HIROKI
  TString s_thickness 	= "/eos/project/m/mtower/public/hiroki/00_util/setting/module_thickness_Feb_2020.txt"     ;
  TString s_alignment 	= "/eos/project/m/mtower/public/hiroki/00_util/setting/chip_correction_Oct_2020.txt"      ;
  TString s_inclination	= "/eos/project/m/mtower/public/hiroki/00_util/setting/detector_inclination_Oct_2020.txt" ;

  //USER:ROBBIE
 
  bool UND = false; //is true if an earlier made TTree is used, so if not using the original data.
  bool testing = false; //Testing outputfile, run over small amount of events.

  bool DB = false; //debug print statements
  bool HP = false; //find hot pixels
  bool CT = true; //Create TTree
  bool C2 = true; //Criterion 2. Every cluster in the first layer with hits behind it in layer 2 is accepted. Only events with 1 accepted cluster are accepted for analysis.
  bool C4 = true; //Criterion 4. Only events where there are no hits outside of a certain area in the second layer centrated around the accepted cluster (C2) in the first layer are accepted for analysis. NEEDS C2
  bool C6 = true; // Criterion 6. Only events with clusters which are not within a certain number of pixels of the layer border are accepted for analysis. NEEDS C2

  // Activate some debug tools?
  bool CheckRejects = false;
  bool CheckThirdLayer = false;
  bool CheckAccepts = false;

  int nPixelRadiusC2 = 10; //The search area behind a cluster is a circle centered around the average coordinate of the cluster. nPixelRadiusC2 is the radius of that circle.
  int nPixelBorderC6 = 50; //The width of the border in which to reject clusters. If a 'particle' cluster is found in layer0 within this many pixels of the layer's border, that event will be rejected. 
  int nPixelRadiusC4 = 120; //Radius of circle for criterion C4

  int nParts = 1; //In how many parts is this run done?
  int part = 1; //Which part is this?

  const int maxNChips = 48; //48 chips is the maximum size of mTower
  const int laneOffset = 32; 
  std::cout << "*** lane offset used is "<<laneOffset<<endl;
  //for the mTower cosmics test and 2019 and 2020 test beam the offset is 32, for the testbeam data of 2018 there is no offset
  const int rowsPerChip = 512; //512 for ALPIDE chip
  const int columnsPerChip = 1024; //1024 for ALPIDE chip
  double nPixelsGap = 5; //Width of gap between chips in units of pixels
  const int laneNumber[6] = {0,3,27,24,2,1}; //corresponds 'lanecode' with laneNumber-laneOffset. Element i is a chip in the (i/2)th layer. Lane 32, 35 in layer 0, 59, 56 in l1, 34, 33 in l2

  //USER:HIROKI
  Double_t dcenter =  5.;
  Double_t dlim    = 12.;
  Double_t wBinLim = 2.*wBinBulk;

  //############################################################                                                                                                                                       
  // basic conditions                                                                                                                                                                                          
  //############################################################                                                                                                                                         
  Double_t energy  = energy_in ;
  TString  suffix0 = Form("%02d_GeV", (Int_t)( 10.*energy ));
  TString  suffix  = Form("%02d_GeV_b%03d_n%02d_m%02d_w%03d_d%03d",
      (Int_t)( 10.*energy       ),
      (Int_t)(100.*wBinBulk     ),
      (Int_t)(     nLayerBulk   ),
      (Int_t)(     minNlayerBulk),
      (Int_t)( 10.*W0Bulk       ),
      (Int_t)( 10.*thDistBulk   )
      );

  //USER:ROBBIE
  //-------------------------------------------------------
  //Setting up some variables automatically
  //-------------------------------------------------------

  if (!(C2)) {C4 = false; C6 = false;} //If C2 is not applied, all of these criteria cannot be applied
  TString beamEnergy; //for masking inputfile
  if (run == 1309 || run == 1310 || run == 1346 || run == 1375 || run == 1376) {
    beamEnergy = "5R8";
    energy_in = 5.8;
  }
  else if (run == 1245 || run == 1250 || run == 1261 || run == 1308 || run == 1333 || run == 1339 || run == 1413) {
    beamEnergy = "5R0";
    energy_in = 5.0;
  }
  else if (run == 1257 || run == 1272 || run == 1274 || run == 1275 || run == 1338 || run == 1345) {
    beamEnergy = "4R0";
    energy_in = 4.0;
  }
  else if (run == 1335 || run == 1341 || run == 1262) {
    beamEnergy = "3R0";
    energy_in = 3.0;
  }
  else if (run == 1276 || run == 1337 || run == 1344) {
    beamEnergy = "2R0";
    energy_in = 2.0;
  }
  else if (run == 1336 || run == 1343 || run == 1263) {
    beamEnergy = "1R0";
    energy_in = 1.0;
  }
  else {
    beamEnergy = "0";
    energy_in = 0.0;
  }

  //later the rows and columns will be translated to an absolute coordinate system with x and y. This are the maximum values of those coordinates
  int maxX = 2*rowsPerChip + nPixelsGap - 1;
  int maxY = columnsPerChip - 1;

  //set plain style for histograms
  gROOT->SetStyle("Plain");

  //-------------------------------------------------------
  //extra masking
  //-------------------------------------------------------

  TH3F* hMaskPtn = new TH3F("hMaskPtn","Hit mask",maxNChips, 0, maxNChips, columnsPerChip, 0, columnsPerChip, rowsPerChip, 0, rowsPerChip); //Histogram with pixels to mask: lane, column, row

  //USER:HIROKI
  //############################################################                                                                                                                        
  // load layer z position                                                                                                                                              
  // load calibration factor                                                                                                                                                                     
  // load inclination function                                                                                                                                                                  
  //############################################################                                                                                                                                
  auto aposz    = mTowerThickness::get_sensor_position( s_thickness.Data() );
  auto paralign = mTowerAlignment::get_sensor_position( s_alignment.Data() );
  auto parinc   = mTowerInclination::get_inclination_function( s_inclination.Data(), energy );

  
  
  //USER:ROBBIE
  if (beamEnergy != "0")
    {
      ifstream in;
      TString maskingFileName = "outMaskPtn_TB_";
      maskingFileName += beamEnergy;
      maskingFileName += "_GeV.txt";
      in.open(Form(maskingFileName,maskingFileLocation.Data()));
      
      Float_t chip_id, nr_lane, hot_pixel_column, hot_pixel_row, pixel_entry;
      Double_t difference, average, std_dev;
      Int_t nlines = 0;
      
      while (1)
	{
	  //in >> chip_id >> nr_lane >> hot_pixel_column >> hot_pixel_row >> pixel_entry >> difference >> average >> std_dev;
	  in >> chip_id >> nr_lane >> hot_pixel_column >> hot_pixel_row >> pixel_entry;
	  if (!in.good()) break;
	  
	  hMaskPtn->Fill(nr_lane-laneOffset,hot_pixel_column,hot_pixel_row);
	  nlines++;
	}
      
      in.close();
    }

  //-------------------------------------------------------
  //create outputfile
  //-------------------------------------------------------

  TString baseName = "Run_";
  TString outputFileName = "";

  if (testing)
    {
      outputFileName = "results_";
      outputFileName += baseName;
      outputFileName += run;
      outputFileName += "_test.root";
    }
  else
    {
      outputFileName = "results_";
      outputFileName += baseName;
      outputFileName += run;
      if (C2)
	{
	  outputFileName += "_C2_";
	  outputFileName += nPixelRadiusC2;
	  outputFileName += "p";
	}
      if (C4)
	{
	  outputFileName += "_C4_";
	  outputFileName += nPixelRadiusC4;
	  outputFileName += "p";
	}
      if (C6)
	{
	  outputFileName += "_C6_";
	  outputFileName += nPixelBorderC6;
	  outputFileName += "pixels";
	}
      if (beamEnergy == "0") outputFileName += "_noExtraMask";
      if (nParts != 1)
	{
	  outputFileName += "_part";
	  outputFileName += part;
	}
      outputFileName += "C246_3layerclusteringfixed.root";
    }
  
  std::cout << "*** Name of outputfile = " << outputFileName << endl;
  
  TFile* outputFile = new TFile(outputFileName,"recreate");
  outputFile->cd();

  //USER:HIROKI
  //############################################################
  // output object definition                            
  //############################################################                                                                           
  TH2D* hBulkSize  = new TH2D("hBulkSize" ,"hBulkSize" , 600 ,-0.5,599.5 ,100, -0.5, 99.5);//all                                                                                                     
  TH2D* hBulkSize0 = new TH2D("hBulkSize0","hBulkSize0", 600 ,-0.5,599.5 ,100, -0.5, 99.5);//<lim                                                                                           
  TH2D* hBulkSize1 = new TH2D("hBulkSize1","hBulkSize1", 600 ,-0.5,599.5 ,100, -0.5, 99.5);//<center                                                                             

  const Int_t nla = 4;
  Int_t max_layer[nla] = {6,12,18,24};
  TH3D** hNhit0    ; InitObj1D<TH3D*>(hNhit0    , nla);// <lim                                                                                                                    
  for(Int_t i=0; i<nla; i++)
    hNhit0[i] = new TH3D(
        Form("hNhit0_%02d",max_layer[i]),
        Form("hNhit0_%02d;Ncore;Ncore^{center};NHitsTot(ila<%d)",max_layer[i],max_layer[i]),
        4, -0.5, 3.5,/*Nshower*/
        4, -0.5, 3.5,/*Nshower_central_region*/
        500, 0, 5000 );

  TH3D** hNhit1    ; InitObj1D<TH3D*>(hNhit1    , nla);// <center                                                                                                                    
  for(Int_t i=0; i<nla; i++)
    hNhit1[i] = new TH3D(
        Form("hNhit1_%02d",max_layer[i]),
        Form("hNhit1_%02d;Ncore;Ncore^{center};NHitsTot(ila<%d)",max_layer[i],max_layer[i]),
        4, -0.5, 3.5,/*Nshower*/
        4, -0.5, 3.5,/*Nshower_central_region*/
        500, 0, 5000 );

  TH1D** hNhitFinal    ; InitObj1D<TH1D*>(hNhitFinal    , nla);
  for(Int_t i=0; i<nla; i++)
    hNhitFinal[i] = new TH1D(
        Form("hNhitFinal_%02d",max_layer[i]),
        Form("hNhitFinal_%02d;NHitsTot(ila<%d)",max_layer[i],max_layer[i]),
	500, 0, 5000 );

  TH3D** hNclus0    ; InitObj1D<TH3D*>(hNclus0    , nla);// <lim                                                                                                                    
  for(Int_t i=0; i<nla; i++)
    hNclus0[i] = new TH3D(
        Form("hNclus0_%02d",max_layer[i]),
        Form("hNclus0_%02d;Ncore;Ncore^{center};NClustersTot(ila<%d)",max_layer[i],max_layer[i]),
        4, -0.5, 3.5,/*Nshower*/
        4, -0.5, 3.5,/*Nshower_central_region*/
        500, 0, 1000 );

  TH3D** hNclus1    ; InitObj1D<TH3D*>(hNclus1    , nla);// <center                                                                                                                          
  for(Int_t i=0; i<nla; i++)
    hNclus1[i] = new TH3D(
        Form("hNclus1_%02d",max_layer[i]),
        Form("hNclus1_%02d;Ncore;Ncore^{center};NClustersTot(ila<%d)",max_layer[i],max_layer[i]),
        4, -0.5, 3.5,/*Nshower*/
        4, -0.5, 3.5,/*Nshower_central_region*/
	500, 0, 1000 );

  TH1D** hNclusFinal    ; InitObj1D<TH1D*>(hNclusFinal    , nla);
  for(Int_t i=0; i<nla; i++)
    hNclusFinal[i] = new TH1D(
        Form("hNclusFinal_%02d",max_layer[i]),
        Form("hNclusFinal_%02d;NClustersTot(ila<%d)",max_layer[i],max_layer[i]),
        500, 0, 1000 );
  
  //USER:ROBBIE
  //-------------------------------------------------------
  //define histograms
  //-------------------------------------------------------
  //output histogram list
  TList* histogramList = new TList(); 

  //number of hits distributions
  TH1D* hHitsDistributionSelection = new TH1D("hitDistrSel","", 5000,0,5000);
  hHitsDistributionSelection -> Sumw2();
  hHitsDistributionSelection -> SetTitle("Distribution of the number of hits per event;number of hits;#");
  histogramList -> Add(hHitsDistributionSelection);

  TH1D* hHitsDistribution = new TH1D("hitsDistribution","Hits Distribution",5000,0,5000); //before selection
  hHitsDistribution -> Sumw2();
  hHitsDistribution -> SetTitle("Distribution of the number of hits per event;number of hits;#");
  histogramList -> Add(hHitsDistribution);
  
  TH2I* hHitMap[maxNChips];
  TH1D* hClusterSize[maxNChips]; 
  TH1D* hClusterSizeLayer[3];
  TH1D* hNClusters[maxNChips];
  TH1D* hNClustersLayer[3];
  TH1D* hNHitsLayer[maxNChips/2];
  TH2D* hNClustersLayer0v1;
  TH2D* hNHitsLayer0v1;
  TH2D* hSidevsNHits;
  
  for (int lane = 0; lane < maxNChips; lane++)
    {
      TString name = "hitMap_lane";
      name += lane;
      hHitMap[lane] = new TH2I(name,"",columnsPerChip,0,columnsPerChip,rowsPerChip,0,rowsPerChip);
      hHitMap[lane] -> SetTitle("Hit Map;column;row");
      histogramList -> Add(hHitMap[lane]);
      
      TString nameCS = "clusterSize_lane";
      nameCS += lane;
      hClusterSize[lane] = new TH1D(nameCS,"",100,0,100);
      hClusterSize[lane] -> Sumw2();
      hClusterSize[lane] -> SetTitle(";cluster size;#");
      histogramList -> Add(hClusterSize[lane]);
      
      //hNClusters
      TString nameNC = "nClusters_lane";
      nameNC += lane;
      hNClusters[lane] = new TH1D(nameNC,"",100,0,100);
      hNClusters[lane] -> Sumw2();
      hNClusters[lane] -> SetTitle(";number of clusters;#");
      histogramList -> Add(hNClusters[lane]);
    }
  
  for (int layer = 0; layer < 3; layer++)
    {
      TString nameCSL = "clusterSize_layer";
      nameCSL += layer;
      hClusterSizeLayer[layer] = new TH1D(nameCSL,"",100,0,100);
      hClusterSizeLayer[layer] -> Sumw2();
      hClusterSizeLayer[layer] -> SetTitle(";cluster size;#");
      histogramList -> Add(hClusterSizeLayer[layer]);
      
      //hNClusters
      TString nameNCL = "nClustersLayer_layer";
      nameNCL += layer;
      hNClustersLayer[layer] = new TH1D(nameNCL,"",100,0,100);
      hNClustersLayer[layer] -> Sumw2();
      hNClustersLayer[layer] -> SetTitle(";number of clusters;#");
      histogramList -> Add(hNClustersLayer[layer]);
      
      TString nameNHL = "nHitsLayer_layer";
      nameNHL += layer;
      hNHitsLayer[layer] = new TH1D(nameNHL,"",100,0,100);
      hNHitsLayer[layer] -> Sumw2();
      hNHitsLayer[layer] -> SetTitle(";number of hits;#");
      histogramList -> Add(hNHitsLayer[layer]);
    }
  
  //hNClusters
  TString nameNCLV = "nClustersLayer_layer0v1";
  hNClustersLayer0v1 = new TH2D(nameNCLV,"",100,0,100,100,0,100);
  hNClustersLayer0v1 -> Sumw2();
  hNClustersLayer0v1 -> SetTitle(";number of clusters;#");
  histogramList -> Add(hNClustersLayer0v1);

  //nHits per layer
  TString nameNHLV = "nHitsLayer_layer0v1";
  hNHitsLayer0v1 = new TH2D(nameNHLV,"",100,0,100,100,0,100);
  hNHitsLayer0v1 -> Sumw2();
  hNHitsLayer0v1 -> SetTitle(";number of hits;#");
  histogramList -> Add(hNHitsLayer0v1);
  
  //Distance from nearest side vs nhits
  TString nameNSVNH = "r_side_vs_nhits";
  hSidevsNHits = new TH2D(nameNSVNH, "", 512, 0, 512, 500, 0, 5000);
  hSidevsNHits -> Sumw2();
  hSidevsNHits -> SetTitle(";Number of hits;Distance from side");
  histogramList -> Add(hSidevsNHits);
  
  //-------------------------------------------------------
  //read the data
  //-------------------------------------------------------
  
  //open the input root file
  fileLocation += baseName;
  fileLocation += run;
  if (UND) fileLocation += "_eventsLeft_normalwindow_C246";
  else
    {
      fileLocation += "/rootout_raw/conv_Run_";
      fileLocation += run;
    }
  fileLocation +=".root";
  std::cout<<endl<<"*** Reading file: "<<fileLocation<<endl<<endl<<"Ignore the following warning if there is one."<<endl;
  TFile* inputFile = TFile::Open(fileLocation);
  std::cout<<endl;

  //USER:HIROKI
  //############################################################                                                                                                                                 
  // get pixels to be masked                                                                                                                                                                    
  //############################################################                                                                                                                                   
  TString s_file, s_ext, s_path ; s_file = GetFileName (fileLocation, s_path, s_ext );
  TString stmp      = Form("outMaskPtn_TB_%dR%d_GeV.txt", (Int_t)energy, (Int_t)(energy*10)%10 );
  TString s_maskptn = Form("%s/hotpixel_TB/%s",s_path.Data(),stmp.Data());
  vector<pair<Int_t,Int_t>> vmask[nchip];
  cout<< "Load hot pixels from "<< s_maskptn.Data() <<"." <<endl;
  for(Int_t ich=0; ich<nchip; ich++) vmask[ich] = mTowerHotPixels::get_hot_pixels( s_maskptn.Data(), chipid2lane_lut.at( ich ) );

  //USER:ROBBIE
  inputFile->cd(); //set to current directory
  
  //get the tree
  TTree* frames = (TTree*)inputFile->Get("Frames");
  //get the branches of the TTree
  int runNumber, fileNumber, eventIndex, nHits, dataSize, eventNumberOriginal;
  frames->SetBranchAddress("runNumber",&runNumber);
  frames->SetBranchAddress("fileNumber",&fileNumber);
  frames->SetBranchAddress("eventNumber",&eventIndex); //Note that this does not run from 0 to nEvents but resets to 0 sometimes
  frames->SetBranchAddress("nHits",&nHits);
  frames->SetBranchAddress("dataSize",&dataSize);
  vector<Int_t>* vlane = new vector<Int_t>();
  vector<Int_t>* vcolumn = new vector<Int_t>();
  vector<Int_t>* vrow = new vector<Int_t>();
  //USER:HIROKI
  Long64_t packetState  ;
  vector<Int_t>*vst_lane   = new vector<Int_t>();
  vector<Int_t>*vst_error  = new vector<Int_t>();
  vector<Int_t>*vst_roflag = new vector<Int_t>();
  //USER:ROBBIE
  frames->SetBranchAddress("lane",&vlane);
  frames->SetBranchAddress("column",&vcolumn);
  frames->SetBranchAddress("row",&vrow);
  frames->SetBranchAddress("packetState" , &packetState );
  frames->SetBranchAddress("st_lane"     , &vst_lane    );
  frames->SetBranchAddress("st_error"    , &vst_error   );
  frames->SetBranchAddress("st_roflag"   , &vst_roflag  );

  int nEvents = frames->GetEntries();

  int nEventsOriginal = nEvents;
  if (UND)
    {
      frames->SetBranchAddress("eventNumberOriginal",&eventNumberOriginal);
      frames->GetEntry(nEvents-1);
      nEventsOriginal = eventNumberOriginal;
    }

  //number of hits vs event number histogram
  outputFile->cd(); //set outputFile as current directory

  TH1I* hHitsvsEvent = new TH1I("hitsvsEvent","",nEventsOriginal,0,nEventsOriginal);
  hHitsvsEvent -> SetTitle("Number of hits per event;event number;number of hits");
  histogramList -> Add(hHitsvsEvent);

  inputFile->cd(); //set inputFile as current directory

  //-------------------------------------------------------
  // Create TTree with selected events
  //-------------------------------------------------------

  //creat the file and the tree
  TFile* outputTree;
  int eventNumberNew;
  int eventNumberOld;
  fileLocationOutputTree += baseName;
  fileLocationOutputTree += run;
  fileLocationOutputTree += "_eventsLeft_normalwindow_C246_layer2fixedclusters.root";
  if (!(CT)) fileLocationOutputTree = "tobedeleted"; //Otherwise a previously made file might be deleted later
  outputTree = new TFile(fileLocationOutputTree,"recreate");
  outputTree->cd();
  TTree Frames("Frames", "Focal frames");
  
  //USER:HIROKI
  TTree* outtree = new TTree("tree","Event Selection");
  Int_t     out_event         ;
  Int_t     out_runNumber     ;
  Int_t     out_nHitsTot      ;
  Int_t     out_nClustersTot  ;
  Int_t     out_status        ;
  Double_t  out_showerX       ;
  Double_t  out_showerY       ;
  
  //USER:ROBBIE

  TTree* prepreparedtree = new TTree("prepreparedselections","Pre-prepared Event Selections");
  Int_t  prepreparedEvent;
  Int_t  prepreparedHiroki;
  Int_t  prepreparedRobbie;

  if (CT) //if CT, create the branches of the tree
    {
      eventNumberNew = -1;
      //Frames.Branch("runNumber",&runNumber,"r/I");
      //Frames.Branch("fileNumber",&fileNumber,"f/I");
      //Frames.Branch("eventNumber",&eventNumberNew,"e/I");
      //Frames.Branch("eventNumberOriginal",&eventNumberOld,"eo/I"); //This number does not reset to 0 sometimes, it goes from 0 to nEvents
      Frames.Branch("nHits",&nHits,"n/I");
      //Frames.Branch("dataSize",&dataSize,"d/I");
      //Frames.Branch("lane",&vlane);
      //Frames.Branch("column",&vcolumn);
      //Frames.Branch("row",&vrow);

      //USER:HIROKI

      outtree->Branch("event"          , &out_event        ,"out_event/I"         );
      outtree->Branch("runNumber"      , &out_runNumber    ,"out_runNumber/I"     );
      outtree->Branch("nHitsTot"       , &out_nHitsTot     ,"out_nHitsTot/I"      );
      outtree->Branch("nClustersTot"   , &out_nClustersTot ,"out_nClustersTot/I"  );
      outtree->Branch("status"         , &out_status       ,"out_status/I"        );
      outtree->Branch("showerX"        , &out_showerX      ,"out_showerX/D"       );
      outtree->Branch("showerY"        , &out_showerY      ,"out_showerY/D"       );
      //USER:ROBBIE
      prepreparedtree->Branch("event", &prepreparedEvent, "event/I");
      prepreparedtree->Branch("HirokiStatus", &prepreparedHiroki, "HirokiStatus/I");
      prepreparedtree->Branch("RobbieStatus", &prepreparedRobbie, "RobbieStatus/I");
    }
  else //else delete the file again. This making and deleting is needed because Frames can not be created in an if statement.
    {
      remove(fileLocationOutputTree);
    }

  //USER:HIROKI
    //############################################################                                                                                                                                                                    
  // initialize layer/chip object                                                                                                                                                                                                  
  //############################################################                                                                                                                                                             
  mTowerLayer* clayers [nlayer];
  mTowerChip*  cchips  [nchip ] ;
  for(Int_t ila=0; ila<nlayer; ila++) clayers[ila] = new mTowerLayer( ila );
  for(Int_t ich=0; ich<nchip ; ich++) cchips [ich] = new mTowerChip ( ich );
  for(Int_t ich=0; ich<nchip ; ich++) cchips [ich]->set_method( mTowerChip::k8N ) ;
  for(Int_t ila=0; ila<nlayer; ila++) clayers[ila]->set_zpos( aposz[ila] );

  //############################################################                                                                                                                                                                       
  // variables for the temporary-use                                                                                                                                                                                                  
  //############################################################                                                                                                                                                                       
  Int_t nbin = 1;
  nbin = (Int_t)(16./wBinBulk);
  TH2D* htmpBulk = new TH2D("htmpBulk", "htmpBulk", 2*nbin, -1.*wBinBulk*nbin, wBinBulk*nbin, 2*nbin, -1.*wBinBulk*nbin, wBinBulk*nbin);

  nbin = (Int_t)(16./wBinLim );
  TH2D* htmpLim  = new TH2D("htmpLim" , "htmpLim ", 2*nbin, -1.*wBinLim *nbin, wBinLim *nbin, 2*nbin, -1.*wBinLim *nbin, wBinLim *nbin);

  //############################################################                                                                                                                                                                    
  // [cell-id, cell-info] map                                                                                                                                                                                                            
  //############################################################                                                                                                                                                          
  std::map<Int_t, cellinfo> mcellsBulk{};
  std::map<Int_t, cellinfo> mcellsLim {};

  //USER:ROBBIE
  inputFile->cd();

  //------------------------------------------------------- 
  // Loop over events 
  //-------------------------------------------------------

  //USER:HIROKI
  Int_t eff_den1 = 0; Int_t eff_neu1 = 0;
  
  //USER:ROBBIE
  std::cout<<"*** Loop over events in input file"<<endl;
  int nReadEvents = 0; //count how many events the loop has read
  int nEmptyEvents = 0; //count the number of "events" that have less than 1 hits
    
  //The events which are looped over
  int minEvent = nEvents/nParts*(part-1);
  int maxEvent = nEvents/nParts*part;  
  if (testing)
    {
      minEvent = 0;
      maxEvent = 100;
    }

  Bool_t isFlushTree = kFALSE ;

  for (int event = minEvent; event < maxEvent; event++)
    //for (int event = minEvent; event < 1000; event++)
    {
      prepreparedEvent = event;
      frames->GetEntry(event);
      if(( event %  10000 == 0) && DB) cout<<"event = "<< event <<endl;

      if (CT)
	{
	  if (UND) eventNumberOld = eventNumberOriginal;
	  else eventNumberOld = event;
	}

      //USER:HIROKI
      Double_t  nHitsTotNew    [nla] = {0};
      Double_t  nClustersTotNew[nla] = {0};
      mcellsBulk.clear(); htmpBulk->Reset("ICESM");
      mcellsLim.clear() ; htmpLim ->Reset("ICESM");
      
      //* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
      // check inclination parameter exists
      //* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
      UInt_t runPeriod = ( runNumber < 1413 )? 0 : 1 ;
      if(runPeriod >= parinc.size() ){ 
	if (DB) cerr<<" === wrong run number ===  or  === failure on loading inclination parameter ===" <<endl; 
	return; 
      }
      

      //* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *                                                                                                                                                              
      // check event status                                                                                                                                                                                                          
      //* * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      Bool_t   isgood = check_event_quality( vlane, vst_lane, vst_error, vst_roflag );
      if( !isgood ) continue;

      //* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *                                                                                                                                                                 
      // set trans parameters                                                                                                                                                                                                               
      //* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *                                                                                                                                                                        
      Double_t ptmp[4] = {0};
      for(Int_t ich=0; ich<nchip ; ich++){
	auto  nrlayer = chipid2layer_lut.at( ich );
	Int_t layer   = get<0>(nrlayer);
	ptmp[0] = get<0>(paralign.at(ich)) - get<0>(parinc[runPeriod])->Eval( aposz[layer] );
	ptmp[1] = get<1>(paralign.at(ich)) - get<1>(parinc[runPeriod])->Eval( aposz[layer] );
	ptmp[2] = get<2>(paralign.at(ich)) ;
	ptmp[3] = get<3>(paralign.at(ich)) ;
	cchips [ich] ->set_trans( ptmp[0], ptmp[1], ptmp[2], ptmp[3] );
      }
      
      //* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *                                                                                                                                                                       
      // reset (except transform parameters)                                                                                                                                                                                             
      //* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *                                                                                                                                                                  
      for(Int_t ich=0; ich<nchip ; ich++) cchips [ich]->clear_chip() ;
      for(Int_t ila=0; ila<nlayer; ila++) clayers[ila]->clear_layer();

      //USER:ROBBIE
      //create mTowerEvent
      mTowerEvent* currentEvent = new mTowerEvent(runNumber,eventIndex);
      currentEvent->setNHits(nHits);
      currentEvent->setNChips(maxNChips);
      TObjArray* hitList = currentEvent->getHits();
      if (DB)
	{
	  if (DB) std::cout<<endl<<"(DB) Run: "<<runNumber<<", event: "<<event<<"/"<<nEvents-1<<", number of hits: "<<nHits;
	  if (UND) std::cout<<", original event number: "<<eventNumberOriginal;
	  if (DB) std::cout<<endl<<"(DB) Loop over hits in event "<<event<<" to add hits to hitlist and apply extra mask"<<endl;
	}

      int nHitsMasked = 0; //number of extra hits masked
      for (int hit=0; hit<nHits; hit++)
	{
	  mTowerHit* currentHit = new mTowerHit();
	  Int_t lane = vlane->at(hit);
	  Int_t col = vcolumn->at(hit);
	  Int_t row = vrow->at(hit);
	  Int_t chipid  = lane2chipidrobbie_lut.at( lane );
	  if (hMaskPtn->GetBinContent(lane-laneOffset+1,col+1,row+1) == 0) //extra masking
	    {
 	      currentHit->setCoordinates(lane, col, row);
	      hitList->Add(currentHit);
	      //USER:HIROKI
	      cchips[chipid]->add_hit( col, row );
	    }
	  //USER:ROBBIE
	  else nHitsMasked++;
	}
      if (DB) std::cout<<"(DB) "<<nHitsMasked<<" extra hits masked in this event"<<endl;
     
      //------------------------------------------------------- 
      //Read the event and do analysis
      //------------------------------------------------------- 

      nReadEvents++;
      
      int entries = hitList->GetEntries();
      
      if (entries < 1) nEmptyEvents++;
      else
      	{
	  //Initialising variables
	  double meanCol[maxNChips] = {0.0};
	  double mean2Col[maxNChips] = {0.0};
	  double meanRow[maxNChips] = {0.0};
	  double mean2Row[maxNChips] = {0.0};
	  int nHitsPerLane[maxNChips]={0};
	  mTowerChipRobbie* hitsInChip[maxNChips];
	  vector<vector<int>> hitsInLayer(maxNChips/2, vector<int>{});
	  for (int l = 0; l<maxNChips ; l++)
	    {
	      hitsInChip[l] = new mTowerChipRobbie(l);
	      hitsInChip[l]->setLane(l+laneOffset);
	    }

	  //loop over all entries
	  if (DB) std::cout<<"(DB) Loop over hits in event "<<event<<" to get properties"<<endl;
	  for (int hit = 0; hit < entries; hit++)
	    {
	      mTowerHit* currentHit = (mTowerHit*)hitList->At(hit);
	      int lane = currentHit->getLane();
	      int row = currentHit->getRow();
	      int column = currentHit->getColumn();

	      if ((lane-laneOffset) <maxNChips &&(lane-laneOffset)>-1 )
		{	  
		  hHitMap[lane-laneOffset]->Fill(column,row,1);
		  nHitsPerLane[lane-laneOffset]++;
		  //hHitsPerLaneNoSel->Fill(lane-laneOffset);
		  meanCol[lane-laneOffset] += column;
		  mean2Col[lane-laneOffset] += column*column;
		  meanRow[lane-laneOffset] += row;
		  mean2Row[lane-laneOffset] += row*row;
		  hitsInLayer[lane2layerrobbie_lut.at(lane)].push_back(hit);
		  hitsInChip[lane-laneOffset]->AddHit(currentHit);
		}
	      else
		{
		  std::cout<<"lane number of hit "<<hit<<" out of range: "<<lane<<endl;
		}
	    } //loop over entries
	  
	  //calculate occupancy per chip
	  /*for (int l = 0; l<maxNChips ; l++)
	    {
	      if (hitsInChip[l]->getNHits() > 0)
		{
		  //hEventsPerLane->Fill(l);
		  hOccupancy[l] -> Fill(double(hitsInChip[l]->getNHits())/double(rowsPerChip*columnsPerChip));
		}
		}*/
	  
	  
	  //Fill raw data histogram
	  hHitsDistribution->Fill(entries);

	  //------------------------------------------------------- 
	  //Clustering
	  //-------------------------------------------------------

	  if (DB) std::cout<<"(DB) CLUSTERING"<<endl;
	  vector<vector<int>> vClusters_layer0(4, vector<int> {}); //Vector with for every cluster: lanecode, meanX, meanY, clustersize
	  vector<vector<int>> vClusters_layer1(4, vector<int> {}); //Vector with for every cluster: lanecode, meanX, meanY, clustersize
	  vector<vector<int>> vClusters_layer2(4, vector<int> {}); //Vector with for every cluster: lanecode, meanX, meanY, clustersize
	  int nClusters = 0;
	  int nClusters_layer1 = 0;
	  int nClusters_layer2 = 0;
	  double minsep = 0;
	  	 
	  hNHitsLayer[0]->Fill((hitsInChip[laneNumber[0]]->getNHits())+(hitsInChip[laneNumber[1]]->getNHits()));
	  hNHitsLayer[1]->Fill((hitsInChip[laneNumber[2]]->getNHits())+(hitsInChip[laneNumber[3]]->getNHits()));
	  hNHitsLayer0v1->Fill(((hitsInChip[laneNumber[0]]->getNHits())+(hitsInChip[laneNumber[1]]->getNHits())), ((hitsInChip[laneNumber[2]]->getNHits())+(hitsInChip[laneNumber[3]]->getNHits())));
     
	  for (int lanecode = 0; lanecode < 2; lanecode++) //loop over chips in first layer to find all clusters
	    {
	      int l = laneNumber[lanecode];
	      if (hitsInChip[l]->getNHits()>0)
		{
		  //get the array of clusters
		  hitsInChip[l]->Clusterize(); // This is the main clustering workhorse. You can find it in mTowerChipRobbie.cxx in the ./classes/ folder
		  TObjArray* clusterlist = hitsInChip[l]->getClusters();
		  for (int c = 0;c<clusterlist->GetEntries();c++) //loop over clusters
		    {
		      mTowerClusterRobbie* cluster = (mTowerClusterRobbie*) clusterlist->At(c);
		      if (cluster)
			{
			  nClusters++;

			  //get properties of cluster
			  if (DB) std::cout<<"(DB) layer0 cluster id = "<<cluster->getId()<<", and lane = "<<cluster->getLane()<<", has "<<cluster->getNHits()<<" hits. This is cluster nr "<<nClusters-1<<endl;
			  vClusters_layer0[0].push_back(lanecode); //lanecode
			  
			  //get mean coordinate of cluster.
			  int meanRow = round(cluster -> getMeanRow());
			  int meanColumn = round(cluster -> getMeanColumn());
			  
			  //change chip dependent coordinates to absolute coordinates   
			  int meanX = meanRow;
			  int meanY = meanColumn;
			  if (IsLeftChip(l+laneOffset)) //cluster is in the left chip
			    {
			      meanX = -meanX - 1 - nPixelsGap;
			      meanY = columnsPerChip - 1 - meanY;
			    }
			  meanX += rowsPerChip + nPixelsGap;

			  vClusters_layer0[1].push_back(meanX); //meanX
			  vClusters_layer0[2].push_back(meanY); //meanY

			  vClusters_layer0[3].push_back(cluster->getNHits()); //clustersize
			  hClusterSize[l]->Fill(cluster->getNHits());
			  hClusterSizeLayer[lane2layerrobbie_lut.at(cluster->getLane())]->Fill(cluster->getNHits()); // Put both in one layer
			  if (DB) std::cout<<"(DB) made it to end of if(cluster)"<<endl;
			}
		    }
		}
	      hNClusters[l]->Fill(nClusters);
	      hNClustersLayer[lane2layerrobbie_lut.at(l+laneOffset)]->Fill(nClusters);
	      if (DB) std::cout<<"(DB) made it to end of layer 0"<<endl;
	    }
	  
	  //Clustering in second layer
	  for (int lanecode = 2; lanecode < 4; lanecode++) //loop over chips in second layer to find all clusters
	    {
	      int l = laneNumber[lanecode];
	      if (hitsInChip[l]->getNHits()>0)
		{
		  //get the array of clusters
		  hitsInChip[l]->Clusterize();
		  TObjArray* clusterlist = hitsInChip[l]->getClusters();
		  for (int c = 0;c<clusterlist->GetEntries();c++) //loop over clusters
		    {
		      mTowerClusterRobbie* cluster = (mTowerClusterRobbie*) clusterlist->At(c);
		      if (cluster)
			{
			  nClusters_layer1++;
			  
			  //get properties of cluster
			  if (DB) std::cout<<"(DB) layer1 cluster id = "<<cluster->getId()<<", and lane = "<<cluster->getLane()<<", has "<<cluster->getNHits()<<" hits. This is cluster nr "<<nClusters_layer1-1<<endl;
			  vClusters_layer1[0].push_back(lanecode); //lanecode
			  
			  //get mean coordinate of cluster.
			  int meanRow = round(cluster -> getMeanRow());
			  int meanColumn = round(cluster -> getMeanColumn());
			  
			  //change chip dependent coordinates to absolute coordinates   
			  int meanX = meanRow;
			  int meanY = meanColumn;
			  if (IsLeftChip(l+laneOffset)) //cluster is in the left chip
			    {
			      meanX = -meanX - 1 - nPixelsGap;
			      meanY = columnsPerChip - 1 - meanY;
			    }
			  meanX += rowsPerChip + nPixelsGap;

			  vClusters_layer1[1].push_back(meanX); //meanX
			  vClusters_layer1[2].push_back(meanY); //meanY
			  vClusters_layer1[3].push_back(cluster->getNHits()); //clustersize
			  hClusterSize[l]->Fill(cluster->getNHits());
			  hClusterSizeLayer[lane2layerrobbie_lut.at(l+laneOffset)]->Fill(cluster->getNHits()); // Put both in one layer
			  if (DB) std::cout<<"(DB) end of clustering"<<endl;
			}
		    }
		}
	      hNClusters[l]->Fill(nClusters_layer1);
	      hNClustersLayer[lane2layerrobbie_lut.at(l+laneOffset)]->Fill(nClusters_layer1);
	      hNClustersLayer0v1->Fill(nClusters, nClusters_layer1);
	      if (DB) std::cout<<"(DB) end of fillin cluster info"<<endl;
	    }

	  for (int lanecode = 4; lanecode < 6; lanecode++) //loop over chips in second layer to find all clusters
	    {
	      int l = laneNumber[lanecode];
	      if (hitsInChip[l]->getNHits()>0)
		{
		  //get the array of clusters
		  hitsInChip[l]->Clusterize();
		  TObjArray* clusterlist = hitsInChip[l]->getClusters();
		  for (int c = 0;c<clusterlist->GetEntries();c++) //loop over clusters
		    {
		      mTowerClusterRobbie* cluster = (mTowerClusterRobbie*) clusterlist->At(c);
		      if (cluster)
			{
			  nClusters_layer2++;
			  
			  //get properties of cluster
			  if (DB) std::cout<<"(DB) layer2cluster id = "<<cluster->getId()<<", and lane = "<<cluster->getLane()<<", has "<<cluster->getNHits()<<" hits. This is cluster nr "<<nClusters_layer2-1<<endl;
			  vClusters_layer2[0].push_back(lanecode); //lanecode
			  
			  //get mean coordinate of cluster.
			  int meanRow = round(cluster -> getMeanRow());
			  int meanColumn = round(cluster -> getMeanColumn());
			  
			  //change chip dependent coordinates to absolute coordinates   
			  int meanX = meanRow;
			  int meanY = meanColumn;
			  if (IsLeftChip(l+laneOffset)) //cluster is in the left chip
			    {
			      meanX = -meanX - 1 - nPixelsGap;
			      meanY = columnsPerChip - 1 - meanY;
			    }
			  meanX += rowsPerChip + nPixelsGap;

			  vClusters_layer2[1].push_back(meanX); //meanX
			  vClusters_layer2[2].push_back(meanY); //meanY
			  vClusters_layer2[3].push_back(cluster->getNHits()); //clustersize

			  hClusterSize[l]->Fill(cluster->getNHits());
			  if (DB) std::cout<<"(DB) got this far."<<endl;
			  hClusterSizeLayer[lane2layerrobbie_lut.at(l+laneOffset)]->Fill(cluster->getNHits()); // Put both in one layer
			  if (DB) std::cout<<"(DB) got further."<<endl;
			  if (DB) std::cout<<"(DB) end of clustering"<<endl;
			}
		    }
		}
	      hNClusters[l]->Fill(nClusters_layer2);
	      hNClustersLayer[lane2layerrobbie_lut.at(l+laneOffset)]->Fill(nClusters_layer2);
	      hNClustersLayer0v1->Fill(nClusters, nClusters_layer2);
	      if (DB) std::cout<<"(DB) end of fillin cluster info"<<endl;
	    }

	  //------------------------------------------------------- 
	  //Event selection
	  //-------------------------------------------------------	      	      

	  bool eventRejected = false; //variable for event criteria
	  int acceptedCluster = -1; //the index of the accepted cluster - for now no cluster accepted
	  int meanXAC; //the mean x of the accepted cluster
	  int meanYAC; //the mean y of the accepted cluster
	  int meanXASC = 0; // the mean x of the layer1 cluster used to verify the particle in layer0.
	  int meanYASC = 0; // the mean y of the layer1 cluster used to verify the particle in layer0.
	  int meanXATC = 0; // the mean x of the layer2 cluster used to verify the particle in layer0.
	  int meanYATC = 0; // the mean y of the layer2 cluster used to verify the particle in layer0.

	  //------------------------------------------------------- 
	  //Criteria C2
	  //-------------------------------------------------------	      	      
	  
	  if (C2 && !(eventRejected))
	    {
	      int nparticles = 0;
	      for (int c = 0; c < nClusters; c++) //loop over clusters for criteria C2
		{
		  if (eventRejected) break;
		  if (vClusters_layer0[0][c] < 2) //if cluster is in first layer
		    {
		      int meanXC = vClusters_layer0[1][c]; //mean x of this cluster
		      int meanYC = vClusters_layer0[2][c]; //mean y of this cluster
		      bool particle = false; // Flag for if this event contains a particle or not.
		      bool doubled = false; // Flag for if this event contains two 'particle' clusters that were identified by the same layer1 cluster.
		      if (DB) std::cout<<"(DB) nClusters_layer2 ="<<nClusters_layer2<<endl;

		      //Instead of looping over hits in second layer, instead loop over clusters in second layer.
		      for (int c2 = 0; c2 < nClusters_layer1; c2++)
			{
			  int meanXC2 = vClusters_layer1[1][c2];
			  int meanYC2 = vClusters_layer1[2][c2];
			  if (DB) std::cout<<"(DB)HERE?"<<endl;
			  if (vClusters_layer1[3][c2] > 1) { //Require that the layer 1 cluster is more than 1 pixel large
			    if (pow(pow(meanXC-meanXC2,2)+pow(meanYC-meanYC2,2),0.5)<nPixelRadiusC2) //Is mean of the layer 1 cluster within shadow of layer 0 cluster?
			      {
				if ((meanXC2 == meanXASC) && (meanYC2 == meanYASC)) { // If true, then we assume this cluster is part of the same particle track as the already accepted particle, so we skip right over it as we already know about this particle!
				  doubled = true;
				  continue;
				}
				particle = true;
				if (CheckRejects || CheckAccepts) {
				  std::cout << "EVENT " << eventIndex << "  Found a particle in layer0 at (" << meanXC << "," << meanYC << "), from a layer1 cluster at (" << meanXC2 << "," << meanYC2 << ")" << endl;
				}
				if (DB) std::cout<<"(DB) cluster nr "<<c<<" accepted by criterion C2"<<endl;
				if (acceptedCluster != -1) //If there is already an accepted cluster
				  {
				    eventRejected = true;
				    if (DB) std::cout<<"(DB) more then 1 accepted cluster in the first layer, so event rejected by criterion C2"<<endl;
				    if (eventRejected && CheckRejects) {
				      std::cout << "+++++++++EVENTREJECTION++++++++++" << endl;
				      std::cout << "EVENT NUMBER " << eventIndex << " WAS REJECTED BY C2, as there was more than one particle. There were " << nHits << " hits in the event." << endl;
				    }
				  }
				acceptedCluster = c;
				nparticles++;
				meanXAC = meanXC;
				meanYAC = meanYC;
				meanXASC = meanXC2;
				meanYASC = meanYC2;
				break; //stop loop over clusters
			      }
			  }
			}
		      if ((!particle) && (!eventRejected) && (!doubled)) { // Loop over clusters in third layer, do we find any evidence here of a throughgoing particle?
			for (int c3 = 0; c3 < nClusters_layer2; c3++)
			  {
			    int meanXC3 = vClusters_layer2[1][c3];
			    int meanYC3 = vClusters_layer2[2][c3];
			    if ((meanXC3 != meanXATC) && (meanYC3 != meanXATC)) { // This statement avoids doubling up particles being identified by the same layer2 cluster.
			      if (vClusters_layer2[3][c3] > 1) { //Require that the layer 2 cluster is more than 1 pixel large
				if (pow(pow(meanXC-meanXC3,2)+pow(meanYC-meanYC3,2),0.5)<(nPixelRadiusC2)) //Is mean of the layer 2 cluster within shadow of layer 0 cluster?
				  {
				    if (CheckRejects || CheckThirdLayer || CheckAccepts) {
				      std::cout << "EVENT " << eventIndex << "  Found a particle in layer0 at (" << meanXC << "," << meanYC << "), from a layer2 cluster at (" << meanXC3 << "," << meanYC3 << ")" << endl;
				    }
				    if (DB) std::cout<<"(DB) cluster nr "<<c<<" accepted by criterion C2"<<endl;
				    if (acceptedCluster != -1) //If there is already an accepted cluster
				      {
					eventRejected = true;
					if (DB) std::cout<<"(DB) more then 1 accepted cluster in the first layer, so event rejected by criterion C2"<<endl;
					if (eventRejected && CheckRejects) {
					  std::cout << "+++++++++EVENTREJECTION++++++++++" << endl;
					  std::cout << "EVENT NUMBER " << eventIndex << " WAS REJECTED BY C2, as there was more than one particle. There were " << nHits << " hits in the event." << endl;
					}
				      }
				    acceptedCluster = c;
				    nparticles++;
				    meanXAC = meanXC;
				    meanYAC = meanYC;
				    meanXATC = meanXC3;
				    meanYATC = meanYC3;
				    break; //stop loop over clusters
				  }
			      }
			    }
			  }
		      } else if (!doubled) { // This loop will just find a corresponding cluster in layer2, it won't accept or reject as we already have a particle we are happy with!
			for (int c3 = 0; c3 < nClusters_layer2; c3++)
			  {
			    int meanXC3 = vClusters_layer2[1][c3];
			    int meanYC3 = vClusters_layer2[2][c3];
			    if ((meanXC3 != meanXATC) && (meanYC3 != meanXATC)) {
			      if (vClusters_layer2[3][c3] > 1) { //Require that the layer 2 cluster is more than 1 pixel large
				if (pow(pow(meanXC-meanXC3,2)+pow(meanYC-meanYC3,2),0.5)<(nPixelRadiusC2)) //Is mean of the layer 2 cluster within shadow of layer 0 cluster?
				  {
				    if (CheckRejects || CheckThirdLayer) {
				      std::cout << "EVENT " << eventIndex << "  Particle in layer0 at (" << meanXC << "," << meanYC << "), verified by a layer2 cluster at (" << meanXC3 << "," << meanYC3 << ")" << endl;
				    }
				    meanXATC = meanXC3;
				    meanYATC = meanYC3;
				    break; //stop loop over clusters
				  }
			      }
			    }
			  }
		      } 
		      if (DB && c != acceptedCluster) std::cout<<"(DB) cluster nr "<<c<<" ignored by criterion C2"<<endl;
		    }//if cluster in first layer
		  if (DB) std::cout << "Event " << eventIndex << ", PARTICLES " << nparticles << endl;
		} //loop over clusters for criterion C2

	      if (acceptedCluster == -1) //If there is no cluster accepted
		{
		  eventRejected = true;
		  if (DB) std::cout<<"(DB) 0 accepted clusters in the first layer, so event rejected"<<endl;
		  if (eventRejected && CheckRejects) {
		    std::cout << "+++++++++EVENTREJECTION++++++++++" << endl;
		    std::cout << "EVENT NUMBER " << eventIndex << " WAS REJECTED BY C2, as there were no particles. There were " << nHits << " hits in the event." << endl;
		  }
		}
	      if (nparticles == 0) {
		eventRejected = true;
		if (eventRejected && CheckRejects) {
		  std::cout << "+++++++++EVENTREJECTION++++++++++" << endl;
		  std::cout << "EVENT NUMBER " << eventIndex << " WAS REJECTED BY C2, as there were no particles. There were " << nHits << " hits in the event." << endl;
		}
	      }
	    } //if C2

	  // Find the minimum separation from the sides.
	  if (meanXAC > 512) {
	    float xsep = 1024.0 - meanXAC;
	    if (meanYAC > 512) {
	      float ysep = 1024.0 - meanYAC;
	      minsep = std::min(xsep, ysep);
	    } else if (meanYAC <= 512) {
	      float ysep = 0.0 + meanYAC;
	      minsep = std::min(xsep, ysep);
	    } else {
	      if (DB) std::cout<<"(DB) SOMETHING WENT WRONG WITH THE SEPARATION: YAC IS NEITHER BIGGER THAN, OR SMALLER/EQUAL TO, 512.0"<<endl;
	      minsep = 0.0;
	    }
	  } else if (meanXAC <= 512) {
	    float xsep = 0.0+meanXAC;
	    if (meanYAC > 512) {
	      float ysep = 1024.0 - meanYAC;
	      minsep = std::min(xsep, ysep);
	    } else if (meanYAC <= 512) {
	      float ysep = 0.0 + meanYAC;
	      minsep = std::min(xsep, ysep);
	    } else {
	      if (DB) std::cout<<"(DB) SOMETHING WENT WRONG WITH THE SEPARATION: YAC IS NEITHER BIGGER THAN, OR SMALLER/EQUAL TO, 512.0"<<endl;
	      minsep = 0.0;
	    }
	  } else {
	    if (DB) std::cout<<"(DB) SOMETHING WENT WRONG WITH THE SEPARATION: XAC IS NEITHER BIGGER THAN, OR SMALLER/EQUAL TO, 512.0"<<endl;
	    minsep = 0.0;
	  }
	  
	  //------------------------------------------------------- 
	  //Criterion C6
	  //-------------------------------------------------------	      	      
	  if ((C6) && !(eventRejected)) 
	    {
	      if ((meanXAC < nPixelBorderC6) || (meanXAC > (columnsPerChip - nPixelBorderC6)) || (meanYAC < nPixelBorderC6) || (meanYAC > (columnsPerChip - nPixelBorderC6)))
		{
		  eventRejected = true;
		  if (DB) std::cout<<"(DB) event rejected by C6 criterion"<<endl;
		}
	      if (eventRejected && CheckRejects) {
		std::cout << "+++++++++EVENTREJECTION++++++++++" << endl;
		std::cout << "EVENT NUMBER " << eventIndex << " WAS REJECTED BY C6, AS IT FELL TOO CLOSE TO THE BORDER OF THE CHIP. There were " << nHits << " hits in the event." << endl;
	      }
	    }
	  
	  //------------------------------------------------------- 
	  //Criterion C4
	  //-------------------------------------------------------	      	      
	  if ((C4) && !(eventRejected))
	    {
	      for (int c2 = 0; c2 < nClusters_layer1; c2++)
		{
		  int meanXC2 = vClusters_layer1[1][c2]; //mean x of this cluster
		  int meanYC2 = vClusters_layer1[2][c2]; //mean y of this cluster
		  if ((pow(pow(meanXAC-meanXC2,2)+pow(meanYAC-meanYC2,2),0.5) > nPixelRadiusC4) && (vClusters_layer1[3][c2] > 1)) { //Is mean of the layer 1 cluster within shadow of layer 0 cluster? Is the cluster size greater than 1?
		    eventRejected = true;
		    if (DB) std::cout<<"(DB) event rejected by C4 criterion"<<endl;
		    if (eventRejected && CheckRejects && (nHits > 2000)) {
		      std::cout << "+++++++++EVENTREJECTION++++++++++" << endl;
		      std::cout << "EVENT NUMBER " << eventIndex << " WAS REJECTED BY C4. There were " << nHits << " hits in the event." << endl;
		    }
		  }
		} //Loop over layer 1 clusters
	    } //C4 loop
	  
	  
	  if (!eventRejected && CheckAccepts && ((nHits >= 2000) || (nHits <= 800))) {
	    std::cout << "++++++++EVENTREJECTION+++++++++++" << endl;
	    std::cout << "WARNING WARNING WARNING: Event " << eventIndex << " was ACCEPTED despite being outside the 'signal' range, with " << nHits << " hits." << endl;
	  }

	  if (eventRejected) {
	    prepreparedRobbie = 0;
	  } else {
	    prepreparedRobbie = 1;
	  }

	  if (!(eventRejected))
	    {
	      if (DB) std::cout << "Filling histograms" << endl;

	      if (UND) hHitsvsEvent->Fill(eventNumberOriginal,entries);
	      else hHitsvsEvent->Fill(event,entries);
	      hHitsDistributionSelection -> Fill(entries);

	      hSidevsNHits -> Fill(minsep, entries);
     
	      for (int l = 0; l<maxNChips ; l++)
		{
		  if (nHitsPerLane[l] > 0)
		    {
		      //calculate mean
		      meanCol[l] = meanCol[l]/nHitsPerLane[l];
		      //hMeanCol[l]->Fill(meanCol[l]);
		      meanRow[l] = meanRow[l]/nHitsPerLane[l];
		      //hMeanRow[l]->Fill(meanRow[l]);
		      //calculate spread
		      mean2Col[l] = mean2Col[l]/nHitsPerLane[l];
		      mean2Row[l] = mean2Row[l]/nHitsPerLane[l];
		      double spreadCol =  TMath::Sqrt( mean2Col[l] - ( meanCol[l]* meanCol[l]) );
		      double spreadRow =  TMath::Sqrt( mean2Row[l] - ( meanRow[l]* meanRow[l]) );
		      //hSpreadCol[l]->Fill(spreadCol);
		      //hSpreadRow[l]->Fill(spreadRow);
		    }
		  delete hitsInChip[l];
		}		      

	      if (CT) //Fill TTree
		{
		  outputTree->cd();
		  eventNumberNew++;
		  Frames.Fill();
		  inputFile->cd();
		}

	    } //if the event is not rejected
	} //end of (Robbie's) event selection
      
      //USER:HIROKI
      //* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *                                                                                                                                                                        
      // fill chips into layer                                                                                                                                                                                                              
      //* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *                                                                                                                                                                        
      for(Int_t ich=0; ich<nchip ; ich++){
	auto      ilayer = get<0>( chipid2layer_lut.at(ich) );
	clayers[ilayer]->add_chip( cchips[ich] );
      }
      
      
      
      //* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *                                                                                                                                                                          
      //total number of hits/clusters                                                                                                                                                                                                     
      //* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *                                                                                                                                                                       
      for(Int_t i=0; i<nla; i++){
	nClustersTotNew[i] = 0;
	nHitsTotNew    [i] = 0;
      }
      for(Int_t ich=0; ich<nchip ; ich++){
	auto lane  = cchips[ich]->get_nrlane();
	auto nclus = cchips[ich]->get_nClusters();
	auto nhits = cchips[ich]->get_nPixels();
	auto scale = cchips[ich]->get_scale();
	for(Int_t i=0; i<nla; i++)
	  if( get<0>( lane2layer_lut.at(lane)) < max_layer[i] ){
	    nClustersTotNew[i] += (Double_t)nclus         ;
	    nHitsTotNew    [i] += (Double_t)nhits * scale ;
	  }
      }
      
      //* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *                                                                                                                                                                        
      // fill hit positon                                                                                                                                                                                                                
      //* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *                                                                                                                                                                     
      for(Int_t ila=0; ila<nlayer; ila++){
	for(UInt_t ic=0; ic<clayers[ila]->get_nChips(); ic++){
	  auto chip   = clayers[ila]->get_chip(ic);
	  
	  for (UInt_t iclus=0; iclus<chip->get_nClusters(); iclus++){
	    auto xy   = chip->get_cluster(iclus)->get_pos_glob_in_mm();
	    Double_t x = get<0>(xy) ;
	    Double_t y = get<1>(xy) ;
	    
	    //Bulk                                                                                                                                                                                                                        
	    if( ila < nLayerBulk ){
	      if( TMath::Abs(x) < dlim && TMath::Abs(y) < dlim ){
		Int_t    icell = htmpBulk->FindBin(x,y);
		Double_t xcell = htmpBulk->GetXaxis()->GetBinCenter( htmpBulk->GetXaxis()->FindBin(x) );
		Double_t ycell = htmpBulk->GetYaxis()->GetBinCenter( htmpBulk->GetYaxis()->FindBin(y) );
		if( mcellsBulk.find(icell) == mcellsBulk.end() ){
		  mcellsBulk.emplace( icell, cellinfo(xcell, ycell) );
		  mcellsBulk.at(icell).Add(x, y, ila);
		}else{
		  mcellsBulk.at(icell).Add(x, y, ila);
		}
	      }
	      else {
		Int_t    icell = htmpLim->FindBin(x,y);
		Double_t xcell = htmpLim->GetXaxis()->GetBinCenter( htmpLim->GetXaxis()->FindBin(x) );
		Double_t ycell = htmpLim->GetYaxis()->GetBinCenter( htmpLim->GetYaxis()->FindBin(y) );
		if( mcellsLim.find(icell) == mcellsLim.end() ){
		  mcellsLim.emplace( icell, cellinfo(xcell, ycell) );
		  mcellsLim.at(icell).Add(x, y, ila);
		}else{
		  mcellsLim.at(icell).Add(x, y, ila);
		}
	      }
	    }
	    
	  }//cluster                                                                                                                                                                                                                      
	}//chip                                                                                                                                                                                                                      
      }//layer          
      
      //* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *                                                                                                                                                                     
      // cell removal                                                                                                                                                                                                                 
      //* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *                                                                                                                                                                
      
      //Bulk                                                                                                                                                                                                                           
      vector<Int_t> vCellRemoveBulk;
      vCellRemoveBulk.clear();
      for ( auto it = mcellsBulk.begin(); it != mcellsBulk.end(); it++ ) {
	auto icell = it->first  ;
	auto info  = it->second ;
	if( info.GetNlayer() < minNlayerBulk || info.GetN() <= W0Bulk )
	  vCellRemoveBulk.push_back( icell );
      }
      for( UInt_t i=0; i<vCellRemoveBulk.size(); i++ )
	mcellsBulk.erase(vCellRemoveBulk.at(i));
      
      //Lim                                                                                                                                                                                                                       
      vector<Int_t> vCellRemoveLim;
      vCellRemoveLim.clear();
      for ( auto it = mcellsLim.begin(); it != mcellsLim.end(); it++ ) {
	auto icell = it->first  ;
	auto info  = it->second ;
	if( info.GetNlayer() < minNlayerBulk || info.GetN() <= W0Bulk )
	  vCellRemoveLim.push_back( icell );
      }
      for( UInt_t i=0; i<vCellRemoveLim.size(); i++ )
	mcellsLim.erase(vCellRemoveLim.at(i));

      //============================================================                                                                                                                                                            
      //* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *                                                                                                                                                                
      // pileup suspicious event (stage1: Bulk)                                                                                                                                                                                 
      //* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *                                                                                                                                                         
      //============================================================                                                                                                                                                          
      if (DB) cerr<<"=================================================="<<endl;
      vector<tuple<Int_t, Double_t, Double_t, Int_t>> vShowerBulk ;
      vShowerBulk.clear();
      {//FJ start                                                                                                                                                                                                          
	vector<PseudoJet> particles;
	particles.clear();
	//Bulk                                                                                                                                                                                                                     
	for ( auto it1 = mcellsBulk.begin(); it1 != mcellsBulk.end(); it1++ ) {
	  //cell center                                                                                                                                                                                                               
	  TVector3 vtmp;
	  //vtmp.SetPtEtaPhi( (it1->second).GetN(), (it1->second).GetXCenter()/10., (it1->second).GetYCenter()/10. );                                                                                            
	  vtmp.SetPtEtaPhi( (it1->second).GetN(), (it1->second).GetX()/10., (it1->second).GetY()/10. );
	  particles.push_back( PseudoJet( vtmp.Px(), vtmp.Py(), vtmp.Pz(), vtmp.Mag()) );
	}
	//Lim                                                                                                                                                                                                    
	for ( auto it1 = mcellsLim.begin(); it1 != mcellsLim.end(); it1++ ) {
	  //cell center                                                                                                                                                                                     
	  TVector3 vtmp;
	  //vtmp.SetPtEtaPhi( (it1->second).GetN(), (it1->second).GetXCenter()/10., (it1->second).GetYCenter()/10. );                                                                                       
	  vtmp.SetPtEtaPhi( (it1->second).GetN(), (it1->second).GetX()/10., (it1->second).GetY()/10. );
	  particles.push_back( PseudoJet( vtmp.Px(), vtmp.Py(), vtmp.Pz(), vtmp.Mag()) );
	}
	
	Double_t R = thDistBulk / 10.;
	JetDefinition jet_def(antikt_algorithm, R);
	ClusterSequence cs(particles, jet_def);
	vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());
	
	// print the jets                                                                                                                                                                                            
	if (DB) cerr <<"event "<< event <<"  (Bulk), "<< nHitsTotNew[nla-1] << endl;
	for (unsigned i = 0; i < jets.size(); i++) {
	  if (DB) cerr << "jet " <<setw(3)<< i << ": "
	       << "  pt: "<< setw(8) << jets[i].pt()
	       << "  x : "<< setw(8) << 10.*jets[i].rap()
	       << "  y : "<< setw(8) << 10.*jets[i].phi_std()
	       << endl;
	  vector<PseudoJet> constituents = jets[i].constituents();
	  for (unsigned j = 0; j < constituents.size(); j++)
	    if (DB) cerr << "    constituent " <<setw(3)<< j << "'s "
		 << "  pt: "<< setw(8) << constituents[j].pt()
		 << "  x : "<< setw(8) << 10.*constituents[j].rap()
		 << "  y : "<< setw(8) << 10.*constituents[j].phi_std()
		 << endl;
	}
	for (unsigned i = 0; i < jets.size(); i++) {
	  if( jets[i].constituents().size() < 2 ) continue;//TODO                                                                                                                                    
	  
	  vShowerBulk.push_back(make_tuple(
					   (Int_t)(jets[i].pt()),
					   10.*jets[i].rap(),
					   10.*jets[i].phi_std(),
					   (Int_t)(jets[i].constituents().size())
					   ));
	}
      }//FJ end
      
      // print shower list                                                                                                                                                                                         
      for( UInt_t i=0; i<vShowerBulk.size(); i++)
	if (DB) cerr
	  <<setw(3) << get<0>(vShowerBulk.at(i))
	  <<setw(3) << get<3>(vShowerBulk.at(i))
	  <<setw(20)<< get<1>(vShowerBulk.at(i))
	  <<setw(20)<< get<2>(vShowerBulk.at(i))
	  << endl;

      //* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *                                                                                                                                               
      // extract stage1 status                                                                                                                                                                  
      //* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *                                                                                                                                
      Bool_t status0 = kFALSE;
      Bool_t status1 = kFALSE;
      {
	Int_t nCenter = 0;
	Int_t nLim    = 0;
	Int_t nAll    = 0;
	for( UInt_t i=0; i<vShowerBulk.size(); i++){
	  Int_t    a = get<0>(vShowerBulk.at(i)) ;
	  Int_t    n = get<3>(vShowerBulk.at(i)) ;
	  Double_t x = get<1>(vShowerBulk.at(i)) ;
	  Double_t y = get<2>(vShowerBulk.at(i)) ;
	  
	  nAll++;
	  if( TMath::Abs(x) < dcenter && TMath::Abs(y) < dcenter ) nCenter++;
	  if( TMath::Abs(x) < dlim    && TMath::Abs(y) < dlim    ) nLim++   ;
	}
	if( nLim   ==nAll && nLim    == 1 ) status0 = kTRUE;
	if( nCenter==nAll && nCenter == 1 ) status1 = kTRUE;
	
	//* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *                                                                                                                                                     
	// fill nHit distributions                                                                                                                                                                  
	//* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *                                                                                                                             
	for(Int_t i=0; i<nla; i++){
	  hNhit0 [i] ->Fill( nAll, nLim   , nHitsTotNew[i]     );
	  hNhit1 [i] ->Fill( nAll, nCenter, nHitsTotNew[i]     );
	  hNclus0[i] ->Fill( nAll, nLim   , nClustersTotNew[i] );
	  hNclus1[i] ->Fill( nAll, nCenter, nClustersTotNew[i] );
	}

	if(nAll == 0) eff_den1++;
	eff_neu1++;
	
	// event of Interest                                                                                                                                                         
	for( UInt_t i=0; i<vShowerBulk.size(); i++){
	  Int_t    a = get<0>(vShowerBulk.at(i)) ;
	  Double_t x = get<1>(vShowerBulk.at(i)) ;
	  Double_t y = get<2>(vShowerBulk.at(i)) ;
	  Int_t    n = get<3>(vShowerBulk.at(i)) ;
	  
	  if( kTRUE   ) hBulkSize ->Fill( a, n );
	  if( status0 ) hBulkSize0->Fill( a, n );
	  if( status1 ) hBulkSize1->Fill( a, n );
	}
      }
      if (DB) cerr<<"=================================================="<<endl;
      
      //* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *                                                                                                                    
      // set tree variables                                                                                                                                                                         
      //* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *                                                                                                         
      out_event         = event                     ;
      out_runNumber     = runNumber               ;
      out_nHitsTot      = nHitsTotNew    [nla-1]  ;
      out_nClustersTot  = nClustersTotNew[nla-1]  ;
      out_showerX = (vShowerBulk.size()==1)? get<1>(vShowerBulk.at(0)) : -999 ;
      out_showerY = (vShowerBulk.size()==1)? get<2>(vShowerBulk.at(0)) : -999 ;
      out_status  = 0;
      if(status0) out_status++;
      if(status1) out_status++;

      //============================================================                                                                                                                                                   
      //* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *                                                                                                                                                                
      // pileup suspicious event (stage2: Core)                                                                                                                                                                                   
      //* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *                                                                                                                                                             
      //============================================================                                                                                                                                                               
      if( !status1 ){ outtree->Fill(); continue; }
      
      Bool_t isSelect = (vShowerBulk.size() < 1 )? kFALSE : kTRUE ;
      
      vector<pair<Double_t, Double_t>> vclus; vclus.clear();
      
      for( UInt_t i=0; i<vShowerBulk.size(); i++){
	Double_t xs = get<1>(vShowerBulk.at(i)) ;
	Double_t ys = get<2>(vShowerBulk.at(i)) ;
	
	Int_t nclus = 0;
	for(Int_t ila=0; ila<2; ila++){
	  for(UInt_t ic=0; ic<clayers[ila]->get_nChips(); ic++){
	    auto chip   = clayers[ila]->get_chip(ic);
	    for (UInt_t iclus=0; iclus<chip->get_nClusters(); iclus++){
	      auto xy   = chip->get_cluster(iclus)->get_pos_glob_in_mm();
	      Double_t xc = get<0>(xy) ;
	      Double_t yc = get<1>(xy) ;
	      Double_t rc = 0.;
	      rc += TMath::Power(xc-xs,2);
	      rc += TMath::Power(yc-ys,2);
	      rc =  TMath::Sqrt(rc);
	      if( rc > thDistBulk ) continue;
	      nclus++;
	      
	      //if (DB) cerr << iev <<" "<< rc <<endl;                                                                                                                                                              
	      
	      isSelect &= ( rc < 2.*wBinBulk );
	      isSelect &= (( ila==0 && nclus>1 )? kFALSE : kTRUE );
	      vclus.push_back( xy );
	    }
	  }
	}
	isSelect &= ( nclus >= 2 );
	if(!isSelect)break;

	for (UInt_t iclus=0    ; iclus<vclus.size(); iclus++){
	  for (UInt_t jclus=iclus; jclus<vclus.size(); jclus++){
	    Double_t x0 = get<0>(vclus.at(iclus)) ;
	    Double_t y0 = get<1>(vclus.at(iclus)) ;
	    Double_t x1 = get<0>(vclus.at(jclus)) ;
	    Double_t y1 = get<1>(vclus.at(jclus)) ;
	    Double_t r = 0.;
	    r += TMath::Power(x1-x0,2);
	    r += TMath::Power(y1-y0,2);
	    r =  TMath::Sqrt(r);
	    isSelect &= ( r < 0.5 );
	    if(!isSelect)break;
	  }
	  if(!isSelect)break;
	}
      }
      
      if( isSelect ){
	for(Int_t i=0; i<nla; i++){
	  hNhitFinal [i]->Fill( nHitsTotNew    [i] );
	  hNclusFinal[i]->Fill( nClustersTotNew[i] );
	}
	//if(nHitsTotNew[3] > 1900) cout<<iev<<endl;                                                                                                                            
	prepreparedHiroki = 1;
      } else {
	prepreparedHiroki = 0;
      }
      
      if (DB) cerr<<endl;


    //* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *                                                                                                                                                   
    // fill tree branches                                                                                                                                                                                        
    //* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *                                                                                                                                                
    if(isSelect) out_status++;
    outtree->Fill();
    prepreparedtree->Fill();
    if( isFlushTree ){ outtree->AutoSave("FlushBaskets") ; outtree->AutoSave("FlushBaskets") ; isFlushTree = kFALSE; }


    //USER:ROBBIE  
    delete currentEvent;
    } //end of loop over events

  //USER:HIROKI
  Double_t eff = 1. - (Double_t)eff_den1/(Double_t)eff_neu1 ;
  if (DB) cout<<"efficiency(Bulk) = "<< /*setprecision(2) <<*/ eff << endl;
  if (DB) cout.copyfmt( oldState );
  
  //############################################################                                                                                                                                                       
  // write output objects                                                                                                                                                                                                       
  //############################################################                                                                                                                                                                  
  TParameter<Double_t>* peff = new TParameter<Double_t>( "efficiency", eff );
  //TNamed                nfile(  "file"      , inputfilename );
  //TNamed                ncond(  "condition" , suffix.Data() );
  
  //USER:ROBBIE
  //(Optional) Find hot pixels and print them
  if (HP)
    {
      std::cout<<"*** Find hot pixels"<<endl;
      for (int l = 0;l<maxNChips;l++)
	{
	  int nHitsLane = hHitMap[l]->GetEntries();
	  if (nHitsLane > 0)
	    {
	      for (int c = 1; c < columnsPerChip+1; c++) //shift by 1 because bins start counting at 1
		{
		  for (int r = 1;r < rowsPerChip+1;r++)
		    {
		      int nHitsBin = hHitMap[l]->GetBinContent(c,r);
		      double fractionBin = (double)nHitsBin/(double)nHitsLane;
		      if (fractionBin > 0.001) std::cout<<"Lane "<<l<<": Bin ("<<c<<", "<<r<<") has "<<nHitsBin<<" entries, "<<fractionBin<<endl;
		    }
		}
	    }
	}
    } //if HP

  if (DB) std::cout<<"*** There were "<<nReadEvents<<" events read by the event loop"<<endl;
  if (DB) std::cout<<"*** There were "<<nEmptyEvents<<" events with less than 1 hits"<<endl;
    
  //----------------------------------------------------------------------------------
 
  //Write the output to disk
  outputFile->cd();
  outputFile->Write();
  outputFile->Close();
  
  if (DB) std::cout<<" I managed to get all the way to saving the tree!";

  if (CT)
    {
      outputTree->cd();
      outputTree->Write();

      //USER:HIROKI
      peff->Write();
      hBulkSize ->Write();
      hBulkSize0->Write();
      hBulkSize1->Write();
      for(Int_t i=0; i<nla; i++) hNhit0      [i] ->Write();
      for(Int_t i=0; i<nla; i++) hNhit1      [i] ->Write();
      for(Int_t i=0; i<nla; i++) hNhitFinal  [i] ->Write();
      for(Int_t i=0; i<nla; i++) hNclus0     [i] ->Write();
      for(Int_t i=0; i<nla; i++) hNclus1     [i] ->Write();
      for(Int_t i=0; i<nla; i++) hNclusFinal [i] ->Write();
      
      //USER:ROBBIE
      outputTree->Close();
    }
}
