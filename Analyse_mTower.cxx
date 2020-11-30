#include "TString.h"
#include "TList.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TROOT.h"
#include "TChain.h"
#include "TMath.h"

#include "classes/mTowerHit.h"
#include "classes/mTowerCluster.h"
#include "classes/mTowerEvent.h"
#include "classes/mTowerChip.h"

#include <iostream>
#include <fstream>
#include <string>


//conversion tables
const std::map< Int_t, Int_t > chipid2lane_lut = {
  { 0,40},{ 1,39},{ 2,42},{ 3,41},{ 4,44},{ 5,43},{ 6,46},{ 7,45},
  { 8,48},{ 9,47},{10,50},{11,49},{12,52},{13,51},{14,54},{15,53},
  {16,38},{17,55},{18,36},{19,37},{20,32},{21,35},{22,34},{23,33},
  {24,64},{25,63},{26,66},{27,65},{28,68},{29,67},{30,70},{31,69},
  {32,72},{33,71},{34,74},{35,73},{36,76},{37,75},{38,78},{39,77},
  {40,62},{41,79},{42,60},{43,61},{44,56},{45,59},{46,58},{47,57}
};

const std::map< Int_t, Int_t > lane2chipid_lut = {
  {40, 0},{39, 1},{42, 2},{41, 3},{44, 4},{43, 5},{46, 6},{45, 7},
  {48, 8},{47, 9},{50,10},{49,11},{52,12},{51,13},{54,14},{53,15},
  {38,16},{55,17},{36,18},{37,19},{32,20},{35,21},{34,22},{33,23},
  {64,24},{63,25},{66,26},{65,27},{68,28},{67,29},{70,30},{69,31},
  {72,32},{71,33},{74,34},{73,35},{76,36},{75,37},{78,38},{77,39},
  {62,40},{79,41},{60,42},{61,43},{56,44},{59,45},{58,46},{57,47}
};

const std::map< Int_t, Int_t > chipid2layer_lut = {
  { 0,22},{ 1,22},{ 2,20},{ 3,20},{ 4,18},{ 5,18},{ 6,16},{ 7,16},
  { 8,14},{ 9,14},{10,12},{11,12},{12,10},{13,10},{14, 8},{15, 8},
  {16, 6},{17, 6},{18, 4},{19, 4},{20, 0},{21, 0},{22, 2},{23, 2},
  {24,23},{25,23},{26,21},{27,21},{28,19},{29,19},{30,17},{31,17},
  {32,15},{33,15},{34,13},{35,13},{36,11},{37,11},{38, 9},{39, 9},
  {40, 7},{41, 7},{42, 5},{43, 5},{44, 1},{45, 1},{46, 3},{47, 3}
};

const std::map< Int_t, Int_t > lane2layer_lut = {
  {40,22},{39,22},{42,20},{41,20},{44,18},{43,18},{46,16},{45,16},
  {48,14},{47,14},{50,12},{49,12},{52,10},{51,10},{54, 8},{53, 8},
  {38, 6},{55, 6},{36, 4},{37, 4},{32, 0},{35, 0},{34, 2},{33, 2},
  {64,23},{63,23},{66,21},{65,21},{68,19},{67,19},{70,17},{69,17},
  {72,15},{71,15},{74,13},{73,13},{76,11},{75,11},{78, 9},{77, 9},
  {62, 7},{79, 7},{60, 5},{61, 5},{56, 1},{59, 1},{58, 3},{57, 3}
};

const std::map<int,bool> layer2isInv_lut = {
  { 0, kFALSE}, { 1, kTRUE}, { 2, kFALSE}, { 3, kTRUE}, 
  { 4, kFALSE}, { 5, kTRUE}, { 6, kFALSE}, { 7, kTRUE}, 
  { 8, kFALSE}, { 9, kTRUE}, {10, kFALSE}, {11, kTRUE}, 
  {12, kFALSE}, {13, kTRUE}, {14, kFALSE}, {15, kTRUE}, 
  {16, kFALSE}, {17, kTRUE}, {18, kFALSE}, {19, kTRUE}, 
  {20, kFALSE}, {21, kTRUE}, {22, kFALSE}, {23, kTRUE} 
};

bool IsLeftChip(int lane){
  int layerNr = lane2layer_lut.at(lane);
  bool isInv  = layer2isInv_lut.at(layerNr);
  int chipid = lane2chipid_lut.at(lane); 
  bool isOdd;
  if (chipid%2 == 1){ isOdd = kTRUE;}
  else {isOdd = kFALSE;}
  bool isLeft = (bool)(isOdd != isInv);
  return isLeft;
}

//for hit maps
int lane2padhitmap(int lane){
  int layerNr = lane2layer_lut.at(lane);
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
  int layerNr = lane2layer_lut.at(lane);
  int padid = layerNr+1;
  return padid; 
}

void Analyse_mTower(int run)
{
  //-------------------------------------------------------
  //macro to read the mTower TB data
  //using event classes: mTowerEvent, mTowerHit, mTowerChip, mTowerCluster
  //detailed code description can be found in the bachelor thesis of Aart van Bochove
  //authors: N. van der Kolk, A. van Bochove
  //-------------------------------------------------------
  
  //-------------------------------------------------------
  //Change this few variables to what applies for your analysis
  //-------------------------------------------------------
  
  TString fileLocationOutputTree = "./";
  TString fileLocation = "/eos/project/m/mtower/Data/Data_TB_February_2020/mTower_Data_DESY_Feb_2020_raw1/"; //The location of the selected data
  TString maskingFileLocation = "/eos/project/m/mtower/public/analysis_fw_sample_HY/data/hotpixel_TB/"; //The location of the masking .txt files
 
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

  //-------------------------------------------------------
  //Setting up some variables automatically
  //-------------------------------------------------------

  if (!(C2)) {C4 = false; C6 = false;} //If C2 is not applied, all of these criteria cannot be applied
  TString beamEnergy; //for masking inputfile
  if (run == 1309 || run == 1310 || run == 1346 || run == 1375 || run == 1376) beamEnergy = "5R8";
  else if (run == 1245 || run == 1250 || run == 1261 || run == 1308 || run == 1333 || run == 1339 || run == 1413) beamEnergy = "5R0";
  else if (run == 1257 || run == 1272 || run == 1274 || run == 1275 || run == 1338 || run == 1345) beamEnergy = "4R0";
  else if (run == 1335 || run == 1341 || run == 1262) beamEnergy = "3R0";
  else if (run == 1276 || run == 1337 || run == 1344) beamEnergy = "2R0";
  else if (run == 1336 || run == 1343 || run == 1263) beamEnergy = "1R0";
  else beamEnergy = "0";

  //later the rows and columns will be translated to an absolute coordinate system with x and y. This are the maximum values of those coordinates
  int maxX = 2*rowsPerChip + nPixelsGap - 1;
  int maxY = columnsPerChip - 1;

  //set plain style for histograms
  gROOT->SetStyle("Plain");

  //-------------------------------------------------------
  //extra masking
  //-------------------------------------------------------

  TH3F* hMaskPtn = new TH3F("hMaskPtn","Hit mask",maxNChips, 0, maxNChips, columnsPerChip, 0, columnsPerChip, rowsPerChip, 0, rowsPerChip); //Histogram with pixels to mask: lane, column, row

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
  frames->SetBranchAddress("lane",&vlane);
  frames->SetBranchAddress("column",&vcolumn);
  frames->SetBranchAddress("row",&vrow);
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
    }
  else //else delete the file again. This making and deleting is needed because Frames can not be created in an if statement.
    {
      remove(fileLocationOutputTree);
    }

  inputFile->cd();

  //------------------------------------------------------- 
  // Loop over events 
  //-------------------------------------------------------
  
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

  for (int event = minEvent; event < maxEvent; event++)
    {
      frames->GetEntry(event);

      if (CT)
	{
	  if (UND) eventNumberOld = eventNumberOriginal;
	  else eventNumberOld = event;
	}

      //create mTowerEvent
      mTowerEvent* currentEvent = new mTowerEvent(runNumber,eventIndex);
      currentEvent->setNHits(nHits);
      currentEvent->setNChips(maxNChips);
      TObjArray* hitList = currentEvent->getHits();
      if (DB)
	{
	  std::cout<<endl<<"(DB) Run: "<<runNumber<<", event: "<<event<<"/"<<nEvents-1<<", number of hits: "<<nHits;
	  if (UND) std::cout<<", original event number: "<<eventNumberOriginal;
	  std::cout<<endl<<"(DB) Loop over hits in event "<<event<<" to add hits to hitlist and apply extra mask"<<endl;
	}

      int nHitsMasked = 0; //number of extra hits masked
      for (int hit=0; hit<nHits; hit++)
	{
	  mTowerHit* currentHit = new mTowerHit();
	  Int_t lane = vlane->at(hit);
	  Int_t col = vcolumn->at(hit);
	  Int_t row = vrow->at(hit);
	  if (hMaskPtn->GetBinContent(lane-laneOffset+1,col+1,row+1) == 0) //extra masking
	    {
 	      currentHit->setCoordinates(lane, col, row);
	      hitList->Add(currentHit);
	    }
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
	  mTowerChip* hitsInChip[maxNChips];
	  vector<vector<int>> hitsInLayer(maxNChips/2, vector<int>{});
	  for (int l = 0; l<maxNChips ; l++)
	    {
	      hitsInChip[l] = new mTowerChip(l);
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
		  hitsInLayer[lane2layer_lut.at(lane)].push_back(hit);
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
		  hitsInChip[l]->Clusterize(); // This is the main clustering workhorse. You can find it in mTowerChip.cxx in the ./classes/ folder
		  TObjArray* clusterlist = hitsInChip[l]->getClusters();
		  for (int c = 0;c<clusterlist->GetEntries();c++) //loop over clusters
		    {
		      mTowerCluster* cluster = (mTowerCluster*) clusterlist->At(c);
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
			  hClusterSizeLayer[lane2layer_lut.at(cluster->getLane())]->Fill(cluster->getNHits()); // Put both in one layer
			  if (DB) std::cout<<"(DB) made it to end of if(cluster)"<<endl;
			}
		    }
		}
	      hNClusters[l]->Fill(nClusters);
	      hNClustersLayer[lane2layer_lut.at(l+laneOffset)]->Fill(nClusters);
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
		      mTowerCluster* cluster = (mTowerCluster*) clusterlist->At(c);
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
			  hClusterSizeLayer[lane2layer_lut.at(l+laneOffset)]->Fill(cluster->getNHits()); // Put both in one layer
			  if (DB) std::cout<<"(DB) end of clustering"<<endl;
			}
		    }
		}
	      hNClusters[l]->Fill(nClusters_layer1);
	      hNClustersLayer[lane2layer_lut.at(l+laneOffset)]->Fill(nClusters_layer1);
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
		      mTowerCluster* cluster = (mTowerCluster*) clusterlist->At(c);
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
			  hClusterSizeLayer[lane2layer_lut.at(l+laneOffset)]->Fill(cluster->getNHits()); // Put both in one layer
			  if (DB) std::cout<<"(DB) got further."<<endl;
			  if (DB) std::cout<<"(DB) end of clustering"<<endl;
			}
		    }
		}
	      hNClusters[l]->Fill(nClusters_layer2);
	      hNClustersLayer[lane2layer_lut.at(l+laneOffset)]->Fill(nClusters_layer2);
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
	} //end of event selection
      delete currentEvent;
    } //end of loop over events

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

  std::cout<<"*** There were "<<nReadEvents<<" events read by the event loop"<<endl;
  std::cout<<"*** There were "<<nEmptyEvents<<" events with less than 1 hits"<<endl;
    
  //----------------------------------------------------------------------------------
 
  //Write the output to disk
  outputFile->cd();
  outputFile->Write();
  outputFile->Close();
  
  std::cout<<" I managed to get all the way to saving the tree!";

  if (CT)
    {
      outputTree->cd();
      outputTree->Write();
      outputTree->Close();
    }
}
