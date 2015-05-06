//------------------------------------------------------------//
//                                                            //
//  PROGRAM FOR "PET" USING DT5740 AND MATRIX 1 (64 crystals) //
//                                                            //
//------------------------------------------------------------//

//NB new version with Short_t and ULong64_t pulizia

// compile with 
// g++ -o PETangerlogic PETangerlogic.cpp `root-config --cflags --glibs` -lSpectrum -lMLP -lTreePlayer -lMinuit 

// The macro wants to take the data acquired from the matrix and be able to assign photoelectric events to each of the 64 crystals.
// Therefore, it is necessary to recognize the interaction crystal among the 4 connected to the MPPC taken into consideration, and 
// discriminate the photoelectric events from comptons.
//
// Just to clarify, the input is a series of events, each one with two time tags (for the moment they are not used, but they will be 
// important in real PET reconstruction) and 16 numbers, i.e. the integrated charge on the 16 MPPC. Since each MPPC have 4 crystals coupled
// it is not know a priori which crystal scintillated
// 
// The "standard" way would be to select from the 16 MPPC spectra only the photoelectric events, then compute the anger logic and have 
// a flood map with 64 nice peaks. This cannot be done if the output of the 4 crystals connected to each mppc is not homogeneous enough, or 
// if the alignment is not good. As it turns out, the current status of the matrix or of the alignment is quite bad so in many channels it is not
// possible to find a single peak. Besides, a more general approach is always preferable.
//
// So the procedure is done as follow:
//
// 1. First all data are used to calculate a flood histogram, without any data cut. Why?
// What can happen when gammas arrive in the matrix is that either
// (a) they deposit all the energy in a single crystal with a nice photoelectric event
// (b) they deposit some of their energy in a single crystal but with some other process, not photoelectric (compton scattering), and then 
// escape from the matrix
// (c) they deposit some energy in one crystal (compton scattering), then interact with another crystal (either in the same MPPC or not)
//     either with compton scattering of with photoelectric. In this case they can eventually escape from the matrix or not.
// And of course this is generally true not only for the 511KeV (of 662 KeV in our tests) we are intersted in, but also for any other gamma
// coming from scattering outside of the matrix of from background radiation
// It is clear that in PET we only want to select events corresponding to 511KeV gammas that are in case (a). 
// If we try to discriminate photoelectric events in the MPPC spectra, we fail, simply because the 4 crystals connected to the same MPPC could have 
// a very different effective light output, either because their conditions are not the same (polishing, damages, damaged wrapping and so on) or because they
// are not perfectly aligned to the MPPC.
// Instead, by calculating an anger logic position for each event, what we get is the interaction position of the gammas, calculated from the light shared 
// to all the 16 MPPCs, and this automatically removes the gammas in case (c), because a weighted position of these gammas will place them in a position
// that does not correspond to one of the 64 peaks (of course, this is true with some efficiency, I would like to see simulations on this issue). Therefore, 
// in a flood histogram without any selection, we are able to see the 64 peaks corresponding to the events of type (a) and (b).
//
// 2. Find the 64 peaks in this flood histogram. They correspond to events that are "contained in one crystal". 
// 3. Fit the peaks, find a position and a sigma (2d, of course)
// 4. For each MPPC, select only the events in this 4 peaks (+/- 2 sigma, for example, or in some other way) and flush them in 4 histogram, one for each crystal
// 5. Plot the 64 crystal histrograms. Now a photopeak should be clearly visible
// 6. Fit the photopeak, select the events in the peak. The final result will be gammas of type (a)


//TODO make the algorithm that looks for the photopeak in the single crystal spectra more solid
//TODO prevent elliptic cuts to cross each other



#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TString.h"
#include "TApplication.h"
#include "TLegend.h"
#include "TTree.h"
#include "TFile.h"
#include "TF2.h"
#include "TSpectrum.h"
#include "TSpectrum2.h"
#include "TTreeFormula.h"
#include "TMath.h"
#include "TChain.h"
#include "TCut.h"
#include "TLine.h"
#include "TError.h"
#include "TEllipse.h"
#include "TMinuit.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>

#include <stdlib.h> 
#include <stdio.h> 
#include <unistd.h>
#include <cmath> 

#include "TRandom1.h"

//#define Params.histomax 16384
#define ENERGY_RESOLUTION 0.12
#define ENERGY_RESOLUTION_SATURATION_CORRECTION 0.35
//----------------------------------------------------------//
//  Creating a struct for input                             //
//----------------------------------------------------------//

struct Input_t
{
  int DigitizerChannel;					//num of the digitizer channel
  TString LabelName;					//name of the MPPC channel
  int CanvasCd;						//position in the 4x4 Canvases, so that the visualization is alwais in the same order
  bool Data;						//this channel have data (true) or no (false)
  int BigCanvasPosition;				//position in the big canvas fo 64 spectra
  bool Module[2];					//module 0 and 1 true or false
  double param0,param1;					//saturation parameters
  
  Float_t *xpeaks;					//array with x positions of peaks in the 2d flood histogram generated accepting all the events
  Float_t *ypeaks;					//array with y positions of peaks in the 2d flood histogram generated accepting all the events
  Float_t *zpeaks;					//array with z positions of peaks in the 2d flood histogram generated accepting all the events
  Double_t *xsigma;					//array with x sigmas of peaks in the 2d flood histogram generated accepting all the events
  Double_t *ysigma;					//array with y sigmas of peaks in the 2d flood histogram generated accepting all the events
  
  float fit2DmeanX[4],fit2DmeanY[4];			//results of the 2d fits for the 4 crystals, means on X and Y
  float fit2DsigmaX[4],fit2DsigmaY[4];  		//results of the 2d fits for the 4 crystals, sigmas on X and Y
  float CrystalMean[4],CrystalSigma[4];			//for each MPPC, 4 crystals. Once the single spectra are extracted, a fit finds out the mean and sigma of the photoelectric peaks
  float LightSharingCrystalMean,LightSharingCrystalSigma;
  float SumCrystalMean[4],SumCrystalSigma[4];
  float ellipseCenterX[4],ellipseCenterY[4];		//values that will be actually used in the elliptic cuts. could be different from the ones above if the ellipes intercept
  float ellipseWidthX[4],ellipseWidthY[4];		//values that will be actually used in the elliptic cuts. could be different from the ones above if the ellipes intercept
  
  //histograms and canvases
  TH1F *channel;					//the histogram of this channel MPPC spectrum -
  
  TH1F *channelClone;					//the histogram of this channel MPPC spectrum - clone
  TH1F *channelOriginal;				//the histogram of this channel MPPC spectrum - no cuts
  //plots in photons
  TH1F *LYchannel;					//same but scaled in Ph/MeV
  TH1F *LYchannelTemp;					//same but scaled in Ph/MeV
  TH1F *LYchannelOriginal;				//same but scaled in Ph/MeV
  //plots in pixels fired
  TH1F *PFchannel;					//same but scaled in Pixels Fired
  TH1F *PFchannelTemp;					//same but scaled in Pixels Fired
  TH1F *PFchannelOriginal;				//same but scaled in Pixels Fired
  
  
  TCanvas *SpectraCanvas;				//canvas for the 2x2 spectra of the 4 crystals in this mppc
  TCanvas *SumSpectraCanvas;
  TH1F *crystal[4];					//4 histograms for the 4 spectra of the crystals in this mppc
  TH1F *crystalClone[4];
  TH1F *crystalHighlight[4];
  TH1F *SumcrystalHighlight[4];
  TH1F *crystalSum[4];					//spectra of the 4 crystals, but using the 16ch sum instead of just trigger mppc
  TH2F *recFlood;					//histogram for the 4x4 flood reconstruction
  TCanvas *recChannelCanvas;
  TCanvas *MPPCCanvas;
  TCanvas *LightSharingCanvas;				//a canvas for the lightsharing when you cut only on the TriggerChannel info;
  TCanvas *LightSharingCanvasPE;			//a canvas for the lightsharing when you cut on the TriggerChannel info and on the PE events;
  TCanvas *LYLightSharingCanvasPE;			//a canvas for the lightsharing when you cut on the TriggerChannel info and on the PE events;

  double lightSharingData[32];				//array to hold the peak positions, in LY, of the sharing data. 32 is just for simmetry..
  
  TH1F *spectraTag[32];					//the 32 spectra when tagging on this channel
  TH1F *spectraTagPE[32];				//the 32 spectra when tagging on this channel and cutting on the PE
  TH1F *LYspectraTagPE[32];
  TF2 *f2;						//2d function to fit the peaks in the flood histos
  TCanvas *CrystalLightSharingCanvas[4];		//4 canvas, one per crystal, to show the light sharing
  TH1F *CrystalLightSharing[4][32];			//4 crystal by 32 by 16 channels...
  TCanvas *CrystalDOICanvas[4];				//4 canvas, one per crystal, to play with DOI
  TCanvas *SumCrystalDOICanvas[4];
  TH1F *DOIRatio[4];					//4 spectra for the DOI ratio (trigger channel / all channels), one per crystal
  TH1F *SumDOIRatio[4];
  TH1F *SumDOIRatioPeak[4];
  double SumLYMean,SumLYSigma;
  
} Input[32];



//----------------------------------------------------------//
//  Creating a struct for parameters                        //
//----------------------------------------------------------//


struct Params_t 
{
  std::string saturationFileName;
  bool coincidence;
  bool deepAnalysis;
  bool doiAnalysis;
  bool saveAnalysisTree;
  bool LightSharingAnalysis;
  bool ChannelOn[32];
  int ChannelModule[32];
  std::string ChannelLabel[32]; 
  int ChannelPosition[32];
  int histomax;
  double MppcPDE;
  double MppcGain;
  double SourceEnergy;
  double AdcChargeBinning;
  double LYhistomax;
  double PFhistomax;
  int binning;
  int LYbinning;
  int PFbinning;
  int histo2DglobalBin;
  int histo2DchannelBin;
  double PEmin;
  double PEmax;
  double lipCutLevel;
  bool CutGrid;
  bool OnlyLyAnalysys;
} Params;

//ttree here, or the functions will never see it
TTree *t1;

//----------------------------------------------------------//
//  Forward declarations                                    //
//----------------------------------------------------------//
int parseConfigFile(std::istream& FileHeader, Params_t &Params);
int printConfigFile(Params_t &Params);
int makeCanvas(TCanvas * Canvas,std::string titleOrName,int ModuleNum,int width,int heigth);

// int makeHisto2D(TTree *t1,TH2F & Histo2D,TCanvas * Canvas,std::string Name,std::string Title,std::string xTitle,std::string yTitle,std::string zTitle,int ModuleNum)
// {
//   //then the first histos...
//   // the flood histo with no cuts...
//   std::stringstream histoStream;
//   std::string histoString;
//   
//   histoStream << titleOrName << ModuleNum;
//   histoString = histoStream.str();
//   //recClean[k] = new TH2F(histoString.c_str(),histoString.c_str(),1000,-7,7,1000,-7,7);
//   Name
//   Histo2D->GetXaxis()->SetTitle(xTitle);
//   Histo2D->GetYaxis()->SetTitle(yTitle);
//   Histo2D->GetZaxis()->SetTitle(zTitle);
//   Histo2D->SetTitle(Title);
//   varStream << "FloodY_" << ModuleNum << ":FloodX_" << ModuleNum << ">>" << histoString ;
//   varString = varStream.str();
//   Canvas[ModuleNum]->cd();
//   t1->Draw(varString.c_str(),"","COLZ");
// } 



//----------------------------------------------------------//
//  Main program                                            //
//----------------------------------------------------------//
int main (int argc, char** argv)
{
  
  //----------------------------------------------------------//
  //  Random stuff...                                         //
  //----------------------------------------------------------//
  //binning choice
  //int Params.histomax = 8196;
  //int BINNING_DIV = 32;
  //define colors of the histograms
  int colors[4] = {2,4,5,6};
  bool ModuleOn[2] = {false,false};
  gStyle->SetOptFit(0001);
  bool correctingForSaturation = false;
  gErrorIgnoreLevel = kError;
  
  //----------------------------------------------------------//
  //  Fix for LIP chip                                        //
  //----------------------------------------------------------//
  //the output from LIP chip suffers from a minor problem when the number of channels triggering is low. It can reco the event in a "grid"
  //all the events in this grid have to be discarded, so if this is the case the user can switch on this cut from the config file and decide the 
  //width of the grid to cut out. Here we define some useful string and streams
  std::stringstream LipCutStream[2];
  std::string lipCutString[2];
  
  
  //----------------------------------------------------------//
  //  Parse the config file                                   //
  //----------------------------------------------------------//
  std::fstream configFile("analysis_config.txt");
  parseConfigFile(configFile,Params);
  printConfigFile(Params);
  //toggle deep analysis
  if(Params.deepAnalysis)
  {
    std::cout << "Software will perform deep analysis" << std::endl;	 
  }
  if(Params.doiAnalysis)
  {
    std::cout << "Software will perform DOI analysis" << std::endl;	 
  }
  
  //calculate scaling factors
  //scaling factor LY: the channels already contain the real ADC counts, now we translate them in LY using the conversion factor derived from the params. It will be the same of course for each chN
  double ChtoLYscalingFactor = (Params.AdcChargeBinning)/(Params.SourceEnergy*Params.MppcGain*Params.MppcPDE*1.6*TMath::Power(10,-19));
  //scaling factor PF: the channels already contain the real ADC counts, now we translate them in Pixels Fired using the conversion factor derived from the params. It will be the same of course for each chN
  double ChtoPFscalingFactor = (Params.AdcChargeBinning)/(Params.MppcGain*1.6*TMath::Power(10,-19));
  
  //----------------------------------------------------------//
  //  Fill input struct with params                           //
  //----------------------------------------------------------//
  for(int i=0; i<32; i++)
  {
    for (int j = 0 ; j < 4 ; j++) //Initialize these data 
    {
      Input[i].fit2DmeanX[j]	= 0;
      Input[i].fit2DmeanY[j]	= 0;
      Input[i].fit2DsigmaX[j]	= 0;
      Input[i].fit2DsigmaY[j]	= 0;
      Input[i].CrystalMean[j] = 0;
      Input[i].CrystalSigma[j] = 0;
      Input[i].ellipseCenterX[j]	= 0;
      Input[i].ellipseCenterY[j]	= 0;
      Input[i].ellipseWidthX[j]	= 0;
      Input[i].ellipseWidthY[j]	= 0;
      
      
    }
    Input[i].DigitizerChannel = i;
    Input[i].LabelName = Params.ChannelLabel[i];
    //int t = (i / 2);
    //Input[i].CanvasCd =  13-t + 2 * (t % 4);			//this way 4x4 canvases start with A1 in bottom left position
    Input[i].CanvasCd = Params.ChannelPosition[i];		//pass the position param to the input struct
    //Input[i].Data = !(i % 2); 				//this way only even channels have data
    Input[i].Data = Params.ChannelOn[i];			//pass the on off param to the input struct
    Input[i].param0=1;						//Initialize saturation params
    Input[i].param1=1;
    if(Params.ChannelModule[i] == 0)
    {
      Input[i].Module[0] = true;
      Input[i].Module[1] = false;
    }
    else
    {
      Input[i].Module[0] = false;
      Input[i].Module[1] = true;
    }
  }
  //so, which modules are on? the ones that have at least one channel on...
  //find out which one is on, and write it in ModuleOn[2]
  for (int i = 0 ; i < 32 ; i++)
  {
    ModuleOn[0] |= Input[i].Data && Input[i].Module[0];
    ModuleOn[1] |= Input[i].Data && Input[i].Module[1];
  }
  
  
  
  
  
  
  
  //----------------------------------------------------------//
  //  Saturation correction - import data if any              //
  //----------------------------------------------------------//
  int iSaturation=0;
  std::ifstream inFile;
  inFile.open(Params.saturationFileName.c_str(),std::ios::in);
  if(inFile.is_open()) 
  {
    correctingForSaturation = true;
    //HISTOMAX = 8196;
    //Params.binning_div = 32;
    std::cout << "Taking Saturation Data from " << Params.saturationFileName << std::endl;	
    while(!inFile.eof())
    {
	//dirty hack to accept both "max adc channels" and "max charge collected" data. It's a terrible implementation and we should really just use charges
        double p0,p1;
        inFile>> p0 >> p1;
	if (p0 > 10.0) //then it means the saturation file reports the p0 and p1 in terms of adc counts in the charge integration plot (i bet you don't collect 10 Coulombs in our setup...)
	{
	  Input[iSaturation].param0 = p0;
	  Input[iSaturation].param1 = p1;
	}
	else //then it means the saturation file reports the p0 and p1 in terms of total charge collected, and we need to convert them back to the adc scale
	{
	  Input[iSaturation].param0 = p0 / Params.AdcChargeBinning;  
	  Input[iSaturation].param1 = p1 / Params.AdcChargeBinning;
	}
        iSaturation++;
    }
  }
  
  //----------------------------------------------------------//
  //  Check input args                                        //
  //----------------------------------------------------------//
  if(argc == 0)
  {
    std::cout << "ERROR: YOU NEED TO PROVIDE AT LEAST AN INPUT FILE!" << std::endl;
  }
  else//if there is at least one file, go!
  {
    //play with strings for the output files...
    std::string firstFileName = argv[1];
    std::size_t found = firstFileName.find_first_of("_",firstFileName.find_first_of("_")+1);
    std::string fileRoot;
    fileRoot = firstFileName.substr(found+1,(firstFileName.length()-(found+1)) -5 );
    //std::cout << "fileRoot " << fileRoot << std::endl;
    if (Params.deepAnalysis)
      fileRoot += "_deep_";
    //std::string fileRoot;
    if (correctingForSaturation)
      fileRoot += "_saturation_corrected.root";
    else
      fileRoot += ".root";
    
    // functions for fitting...
    TF1 *gauss = new TF1("gauss",  "[0]*exp(-0.5*((x-[1])/[2])**2)",0,20000); //quite incredibly, this means that [2] can be negative... so the fix will be to take always the module of [2]...
    TF1 *LYgauss = new TF1("LYgauss",  "[0]*exp(-0.5*((x-[1])/[2])**2)",0,60000); //quite incredibly, this means that [2] can be negative... so the fix will be to take always the module of [2]...
    TF1 *LYSumgauss = new TF1("LYSumgauss",  "[0]*exp(-0.5*((x-[1])/[2])**2)",0,60000); //quite incredibly, this means that [2] can be negative... so the fix will be to take always the module of [2]...
    TF1 *LSgauss = new TF1("LSgauss",  "[0]*exp(-0.5*((x-[1])/[2])**2)",0,60000); //quite incredibly, this means that [2] can be negative... so the fix will be to take always the module of [2]...
    
    //----------------------//
    //  TChain              //
    //----------------------//
    
    //tchain to merge the input root ttrees
    TChain *chain =  new TChain("adc");
    for (int i = 1 ; i < argc ; i++)
    {
      std::cout << "Adding file " << argv[i] << std::endl;
      chain->Add(argv[i]);
    }
    
    //----------------------//
    //  TTree               //
    //----------------------//
    //and now we make the real ttree used in this program
    
    //original ones
    //Double_t xmppc[32]={-4.8,0,-1.6,0,1.6,0,4.8,0,-4.8,0,-1.6,0,1.6,0,4.8,0,-4.8,0,-1.6,0,1.6,0,4.8,0,-4.8,0,-1.6,0,1.6,0,4.8,0};
    //Double_t ymppc[32]={-4.8,0,-4.8,0,-4.8,0,-4.8,0,-1.6,0,-1.6,0,-1.6,0,-1.6,0,1.6,0,1.6,0,1.6,0,1.6,0,4.8,0,4.8,0,4.8,0,4.8,0};
    
    //temporary for zheng - nino (but they are actually the same as original ones...)
    //Double_t xmppc[32]={-4.8,0,-1.6,0,1.6,0,4.8,0,-4.8,0,-1.6,0,1.6,0,4.8,0,-4.8,0,-1.6,0,1.6,0,4.8,0,-4.8,0,-1.6,0,1.6,0,4.8,0};
    //Double_t ymppc[32]={-4.8,0,-4.8,0,-4.8,0,-4.8,0,-1.6,0,-1.6,0,-1.6,0,-1.6,0,1.6,0,1.6,0,1.6,0,1.6,0,4.8,0,4.8,0,4.8,0,4.8,0};
    
    //mapping for the new boards
    Double_t xmppc[32]={-4.8,-4.8,-4.8,-4.8,-1.6,-1.6,-1.6,-1.6,1.6,1.6,1.6,1.6,4.8,4.8,4.8,4.8,-4.8,-4.8,-4.8,-4.8,-1.6,-1.6,-1.6,-1.6,1.6,1.6,1.6,1.6,4.8,4.8,4.8,4.8};
    Double_t ymppc[32]={4.8,1.6,-1.6,-4.8,4.8,1.6,-1.6,-4.8,4.8,1.6,-1.6,-4.8,4.8,1.6,-1.6,-4.8,4.8,1.6,-1.6,-4.8,4.8,1.6,-1.6,-4.8,4.8,1.6,-1.6,-4.8,4.8,1.6,-1.6,-4.8};
    
    //variables for the flood histogram computation
    Double_t posx;
    Double_t posy;
    //create the new ttree
    std::stringstream snames[32],stypes[32];
    std::string names[32],types[32];
    
    t1 = new TTree("adc","adc"); 
    
    //variables
    ULong64_t ExtendedTimeTag,t1_ExtendedTimeTag;
    ULong64_t DeltaTimeTag,t1_DeltaTimeTag;
    int TriggerChannel[2];
    Short_t charge[32];
    Short_t t1_charge[32];
    ULong64_t counter = 0;
    ULong64_t badEvents = 0;
    float floodx[2],floody[2],firstonsecond[2],ratiotriggeronall[2];
    Bool_t badevent;
    
    //branches
    TBranch *b_ExtendedTimeTag;
    TBranch *b_DeltaTimeTag;
    TBranch *b_charge[32];
    TBranch *b_TriggerChannel;
    TBranch *b_floodx;
    TBranch *b_floody;
    TBranch *b_firstonsecond;
    chain->SetBranchAddress("ExtendedTimeTag", &ExtendedTimeTag, &b_ExtendedTimeTag);
    chain->SetBranchAddress("DeltaTimeTag", &DeltaTimeTag, &b_DeltaTimeTag);
    for(int i=0; i<32; i++)
    {
      snames[i].str(std::string());
      stypes[i].str(std::string());
      snames[i] << "ch" << i;
      names[i] = snames[i].str();
      chain->SetBranchAddress(names[i].c_str(), &charge[i], &b_charge[i]);
    }
    
    //create the second ttree, t1
    //first 2 branches of the ttree
    t1->Branch("ExtendedTimeTag",&t1_ExtendedTimeTag,"ExtendedTimeTag/l"); 	//absolute time tag of the event
    t1->Branch("DeltaTimeTag",&t1_DeltaTimeTag,"DeltaTimeTag/l"); 			//delta time from previous event
    //branches of the 32 channels data
    for (int i = 0 ; i < 32 ; i++){
      //empty the stringstreams
      snames[i].str(std::string());
      stypes[i].str(std::string());
      t1_charge[i] = 0;
      snames[i] << "ch" << i;
      stypes[i] << "ch" << i << "/S";  
      names[i] = snames[i].str();
      types[i] = stypes[i].str();
      t1->Branch(names[i].c_str(),&t1_charge[i],types[i].c_str());
    }
    
    
    for (int k = 0 ; k < 2 ; k++)
    {
      if(ModuleOn[k])
      {
	std::stringstream nameStream,typeStream;
	std::string nameString,typeString;
	
	nameStream << "TriggerChannel_" << k;
	nameString = nameStream.str();
	typeStream << "TriggerChannel_" << k << "/I";
	typeString = typeStream.str();
	t1->Branch(nameString.c_str(),&TriggerChannel[k],typeString.c_str());		//TriggerChannel[2] = channel where the signal is the highest for module [0] and [1]
	
	nameStream.str(std::string());
	typeStream.str(std::string());
	nameStream << "FloodX_" << k;
	nameString = nameStream.str();
	typeStream << "FloodX_" << k << "/F";
	typeString = typeStream.str();
	t1->Branch(nameString.c_str(),&floodx[k],typeString.c_str());					//x position in the complete flood histogram
	
	nameStream.str(std::string());
	typeStream.str(std::string());
	nameStream << "FloodY_" << k;
	nameString = nameStream.str();
	typeStream << "FloodY_" << k << "/F";
	typeString = typeStream.str();
	t1->Branch(nameString.c_str(),&floody[k],typeString.c_str());					//y position in the complete flood histogram
	
	nameStream.str(std::string());
	typeStream.str(std::string());
	nameStream << "FirstOnSecond_" << k;
	nameString = nameStream.str();
	typeStream << "FirstOnSecond_" << k << "/F";
	typeString = typeStream.str();
	t1->Branch(nameString.c_str(),&firstonsecond[k],typeString.c_str());					//ratio highest signal to second highest 
	
	nameStream.str(std::string());
	typeStream.str(std::string());
	nameStream << "RatioTriggerOnAll_" << k;
	nameString = nameStream.str();
	typeStream << "RatioTriggerOnAll_" << k << "/F";
	typeString = typeStream.str();
	t1->Branch(nameString.c_str(),&ratiotriggeronall[k],typeString.c_str());					//ratio highest signal to second highest 
	
	
	
      }
    }
    
    t1->Branch("BadEvent",&badevent,"BadEvent/O");
    
    std::cout << "Filling the TTree for the analysis... " << std::endl;
    Int_t nevent = chain->GetEntries();
    long long int GoodCounter = 0;
    
    for (Int_t i=0;i<nevent;i++) { //loop on all the entries of tchain
      chain->GetEvent(i);              //read complete accepted event in memory
      double maxCharge[2] = {0,0};
      double secondCharge[2] = {0,0};
      float columnsum[2]={0,0};
      float rowsum[2]={0,0};
      float total[2]={0,0};
      int badCharge = 0;
      badevent = false;
      
      t1_ExtendedTimeTag = ExtendedTimeTag;
      t1_DeltaTimeTag = DeltaTimeTag;
      for (int i = 0 ; i < 32 ; i++) //loop on channels
      {
	if(Input[i].Data) //if channel is on
	{
	  if(correctingForSaturation /*&& i != 6*/) //temporary 
	  {
	    if(charge[i] > Input[i].param0)
	    {
	      badCharge++;
	      //std::cout << "BadCharge " << (int)round(-Input[i].param0 * TMath::Log(1.0 - ( charge[i]/Input[i].param0 ))) << std::endl;
	    }
	    t1_charge[i] = (int)round(-Input[i].param0 * TMath::Log(1.0 - ( charge[i]/Input[i].param0 )));	//this way, the param0 input has to be calculated in terms of number of bins in the binning of this very measurement. Usually we acquire at 156fC binning
	    if(charge[i] > Input[i].param0)
	    {
	      //std::cout << "BadCharge " << t1_charge[i] << std::endl;
	    }
	    
	  }
	  else
	    t1_charge[i] = charge[i];
	}
	else 
	  t1_charge[i] = 0;
	for( int k = 0 ; k < 2 ; k++) //cycle on modules
	{
	  if(Input[i].Module[k]) //if channel is in this module
	  {
	    if(Input[i].Data) //if channel is on
	    {
	      if (t1_charge[i] > maxCharge[k])
	      {
		maxCharge[k] = t1_charge[i];
		TriggerChannel[k] = i;
	      }
	      total[k] += t1_charge[i];
	      rowsum[k] += t1_charge[i]*xmppc[i];
	      columnsum[k] += t1_charge[i]*ymppc[i];
	    }
	  }
	}
      }
      
      if(badCharge)
      {
	badevent = true;
	badEvents++;
      }
      
      //find second highest charge
      for (int k = 0 ; k < 2 ; k++)
      {
	for (int i = 0 ; i < 32 ; i++)
	{
	  if(Input[i].Module[k]) //if channel is in this module
	  {
	    if(Input[i].Data) //if channel is on
	    {
	      if (t1_charge[i] != maxCharge[k])//TODO here actually I'm not considering the case (unlikely) there are two channels with same value = maxcharge
	      {
		if (t1_charge[i] > secondCharge[k])
		{
		  secondCharge[k] = t1_charge[i];
		}
	      }
	    }
	  }
	}
      }
      for (int k = 0 ; k < 2 ; k++)
      {
	if(ModuleOn[k])
	{
	  //compute flood x and y
	  floodx[k]=rowsum[k]/total[k];
	  floody[k]=columnsum[k]/total[k];
	  //compute first on second ratio
	  firstonsecond[k] = maxCharge[k] / secondCharge[k];
	}
      }
      
      //compute the ratio trigger/sum 
      for (int k = 0 ; k < 2 ; k++)
      {
	if(ModuleOn[k])
	{
	  ratiotriggeronall[k] = maxCharge[k]/total[k];
	}
      }
      
      
      if(!badevent)
      {
	t1->Fill();//fills the tree with the data only if it's a good event
	GoodCounter++;
      }
      //counter to give a feedback to the user
      counter++;
      
      int perc = ((100*counter)/nevent); //should strictly have not decimal part, written like this...
      if( (perc % 10) == 0 )
      {
	std::cout << "\r";
	std::cout << perc << "% done... ";
	//std::cout << counter << std::endl;
      }
    }
    std::cout << std::endl;
    
    std::cout << "Tot events = \t" << counter << std::endl;
    std::cout << "Accepted events = \t" << GoodCounter << std::endl;
    std::cout << "Bad events = \t" << badEvents << std::endl;
    //-----------------------------
//     TFile* fTree = new TFile("tempSaturation.root","recreate");
//     fTree->cd();
//     t1->Write();
//     fTree->Close();
//     std::cout << "DONE" << std::endl;
    //-----------------------------
    

//     
    //----------------------------------------------------------//
    //   Useful arrays                                          //
    //----------------------------------------------------------//
    double Pedestal[32] = //in case the pedestal is not in zero... but we tested that it is in zero, so...
    { 
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0
    };
    //position of each crystal in the big canvas
    int BigCanvasMapping[64] = 
    {
      49,
      50,
      57,
      58,
      51,
      52,
      59,
      60,
      53,
      54,
      61,
      62,
      55,
      56,
      63,
      64,
      33,
      34,
      41,
      42,
      35,
      36,
      43,
      44,
      37,
      38,
      45,
      46,
      39,
      40,
      47,
      48,
      17,
      18,
      25,
      26,
      19,
      20,
      27,
      28,
      21,
      22,
      29,
      30,
      23,
      24,
      31,
      32,
      1,
      2,
      9,
      10,
      3,
      4,
      11,
      12,
      5,
      6,
      13,
      14,
      7,
      8,
      15,
      16
    };

    
    //----------------------------------------------------------//
    //  RUN ON THE CHANNELS                                     //
    //----------------------------------------------------------//
    
    //Create the canvases etc 
    
    TCanvas *FloodCanvasClean[2], *CanvasFirstOnSecond[2], *FloodCanvasCutOnRatio[2], *RecCanvasFull[2],*CanvasOriginal[2], *Canvas1[2];
    TCanvas *LYCanvasOriginal[2], *LYCanvas1[2];
    TCanvas *PFCanvasOriginal[2], *PFCanvas1[2];
    TH2F *recClean[2], *recCutOnRatio[2], *recFinal[2];
    TH1F *HistoFirstOnSecond[2];
    TCanvas *FloodCanvas[2]; 
    TCanvas *BigSpectraCanvas[2];
    TCanvas *SumBigSpectraCanvas[2];
    
    int crystalCounter[2] = {0,0};
    int crystalCounter2[2] = {0,0};
    
    bool peak2Dfound[2][64];
    for(int k = 0; k < 2 ; k++) //run on modules (0 and 1)
    {
      for(int i = 0; i < 64 ; i++)//run on all possible crystals (64)
      {
	peak2Dfound[k][i] = false; 
      }
    }
    
    std::ofstream outputCutFile[2];
    
    
    for (int k = 0; k < 2 ; k++) //run on modules (0 and 1)
    {
      if(ModuleOn[k]) //if the module is set as ON
      {
	
	std::cout << "Analyzing Module " << k << "..." << std::endl;
	//a string and a stream to rule them all...
	std::stringstream histoStream,varStream,cutStream;
	std::string histoString,varString,cutString;
	
	//we want to save in a text file the cuts for each crystal area in the flood plot
	//
	std::stringstream outputCutFileStream;
        outputCutFileStream << "Module_" << k << "_Cuts.dat";
	outputCutFile[k].open(outputCutFileStream.str().c_str(), std::ofstream::out); 
	
	//Canvas
	//not really needed to use external function, but it was done to test something...
	//create the canvas
	FloodCanvasClean[k]		= new TCanvas();
	CanvasFirstOnSecond[k]		= new TCanvas();
	FloodCanvasCutOnRatio[k]	= new TCanvas();
	RecCanvasFull[k]		= new TCanvas();
	CanvasOriginal[k]		= new TCanvas();
	Canvas1[k]			= new TCanvas();
	LYCanvasOriginal[k]		= new TCanvas();
	LYCanvas1[k]			= new TCanvas();
	PFCanvasOriginal[k]		= new TCanvas();
	PFCanvas1[k]			= new TCanvas();
	FloodCanvas[k]			= new TCanvas();
	BigSpectraCanvas[k]		= new TCanvas();
	SumBigSpectraCanvas[k]		= new TCanvas();
	
	//give them names and dimensions
	makeCanvas(FloodCanvasClean[k],"FloodCanvasClean-Module",k,1000,1000);
	makeCanvas(CanvasFirstOnSecond[k],"CanvasFirstOnSecond-Module",k,800,600);
	makeCanvas(FloodCanvasCutOnRatio[k],"FloodCanvasCutOnRatio-Module",k,1200,800);
	makeCanvas(RecCanvasFull[k],"RecCanvasFull-Module",k,1000,1000);
	makeCanvas(CanvasOriginal[k],"AdcCountsOriginal-Module",k,1200,800);
	makeCanvas(Canvas1[k],"AdcCounts_Trigger-Module",k,1200,800);
	makeCanvas(LYCanvasOriginal[k],"LightYield_Trigger-Module",k,1200,800);
	makeCanvas(LYCanvas1[k],"SumLightYield_Trigger-Module",k,1200,800);
	makeCanvas(PFCanvasOriginal[k],"PixelsFired_Trigger-Module",k,1200,800);
	makeCanvas(PFCanvas1[k],"SumPixelsFired_Trigger-Module",k,1200,800);
	makeCanvas(FloodCanvas[k],"FloodCanvas-Module",k,1200,800);
	makeCanvas(BigSpectraCanvas[k],"BigSpectraCanvas-Module",k,1200,800);
	makeCanvas(SumBigSpectraCanvas[k],"SumBigSpectraCanvas-Module",k,1200,800);
	
	//divide the ones that need it
	CanvasOriginal[k]->Divide(4,4);
	Canvas1[k]->Divide(4,4);
	LYCanvasOriginal[k]->Divide(4,4);
	LYCanvas1[k]->Divide(4,4);
	PFCanvasOriginal[k]->Divide(4,4);
	PFCanvas1[k]->Divide(4,4);
	FloodCanvas[k]->Divide(4,4);
	BigSpectraCanvas[k]->Divide(8,8);
	SumBigSpectraCanvas[k]->Divide(8,8);
	
	// then the first histos...
	// the flood histo with "no cuts"...
	histoStream << "recClean-Module" << k;
	histoString = histoStream.str();
	recClean[k] = new TH2F(histoString.c_str(),histoString.c_str(),Params.histo2DglobalBin,-7,7,Params.histo2DglobalBin,-7,7);
	recClean[k]->GetXaxis()->SetTitle("Flood X");
	recClean[k]->GetYaxis()->SetTitle("Flood Y");
	recClean[k]->GetZaxis()->SetTitle("N");
	recClean[k]->SetTitle("Reconstruction of entire dataset");
	varStream << "FloodY_" << k << ":FloodX_" << k << ">>" << histoString ;
	varString = varStream.str();
	FloodCanvasClean[k]->cd();
	
	// Cut of grid, useful for Lip chip input
	if(Params.CutGrid){
	LipCutStream[k] <<  "((FloodX_" 
	          << k
	          << " < (-4.8 -"
		  << Params.lipCutLevel
		  << "))"
	          << " | "
	          <<  "((FloodX_" 
		  << k
		  << " > (-4.8 +"
		  << Params.lipCutLevel
		  << "))"
		  << " && "
		  <<  "(FloodX_" 
		  << k
		  << " < (-1.6 -"
		  << Params.lipCutLevel
		  << ")))"
		  << " | "
		  <<  "((FloodX_" 
		  << k
		  << " > (-1.6 +"
		  << Params.lipCutLevel
		  << "))"
		  << " && "
		  <<  "(FloodX_" 
		  << k
		  << " < (1.6 -"
		  << Params.lipCutLevel
		  << ")))"
		  << " | "
		  <<  "((FloodX_" 
		  << k
		  << " > (1.6 +"
		  << Params.lipCutLevel
		  << "))"
		  << " && "
		  <<  "(FloodX_" 
		  << k
		  << " < (4.8 -"
		  << Params.lipCutLevel
		  << ")))"
		  << " | "
		  <<  "(FloodX_" 
		  << k
		  << " > (4.8 +"
		  << Params.lipCutLevel
		  << ")))"
		  << " && "
		  <<  "((FloodY_" 
		  << k
		  << " < (-4.8 -"
		  << Params.lipCutLevel
		  << "))"
		  << " | "
		  <<  "((FloodY_" 
		  << k
		  << " > (-4.8 +"
		  << Params.lipCutLevel
		  << "))"
		  << " && "
		  <<  "(FloodY_" 
		  << k
		  << " < (-1.6 -"
		  << Params.lipCutLevel
		  << ")))"
		  << " | "
		  <<  "((FloodY_" 
		  << k
		  << " > (-1.6 +"
		  << Params.lipCutLevel
		  << "))"
		  << " && "
		  <<  "(FloodY_" 
		  << k
		  << " < (1.6 -"
		  << Params.lipCutLevel
		  << ")))"
		  << " | "
		  <<  "((FloodY_" 
		  << k
		  << " > (1.6 +"
		  << Params.lipCutLevel
		  << "))"
		  << " && "
		  <<  "(FloodY_" 
		  << k
		  << " < (4.8 -"
		  << Params.lipCutLevel
		  << ")))"
		  << " | "
		  <<  "(FloodY_" 
		  << k
		  << " > (4.8 +"
		  << Params.lipCutLevel
		  << ")))";
	}
	else
	  LipCutStream[k] << "";
	
	lipCutString[k] = LipCutStream[k].str();
	
	TCut lipCut = lipCutString[k].c_str();
	TCut badEventCut = "BadEvent == 0";
	t1->Draw(varString.c_str(),lipCut+badEventCut,"COLZ");
	histoStream.str(std::string());
	varStream.str(std::string());
	
	//draw the mppc real dimensions on the 2d plot
	TLine *line[10];
	line[0] = new TLine(-6.4,-6.4,-6.4,6.4);
	line[1] = new TLine(-3.2,-6.4,-3.2,6.4);
	line[2] = new TLine(0,-6.4,0,6.4);
	line[3] = new TLine(3.2,-6.4,3.2,6.4);
	line[4] = new TLine(6.4,-6.4,6.4,6.4);
	line[5] = new TLine(-6.4,-6.4,6.4,-6.4);
	line[6] = new TLine(-6.4,-3.2,6.4,-3.2);
	line[7] = new TLine(-6.4,0,6.4,0);
	line[8] = new TLine(-6.4,3.2,6.4,3.2);
	line[9] = new TLine(-6.4,6.4,6.4,6.4);
	for ( int iLine = 0 ; iLine < 10 ; iLine++)
	{
	  line[iLine]->SetLineColor(kRed);
	  line[iLine]->SetLineWidth(2);
	  line[iLine]->Draw();
	}

	//then, plot the firstonsecond, to see what it looks like... //not useful
	histoStream << "HistoFirstOnSecond-Module" << k;
	histoString = histoStream.str();
	HistoFirstOnSecond[k] = new TH1F(histoString.c_str(),histoString.c_str(),400,0,10);
	HistoFirstOnSecond[k]->GetXaxis()->SetTitle("first charge / second charge");
	HistoFirstOnSecond[k]->GetYaxis()->SetTitle("N");
	HistoFirstOnSecond[k]->SetTitle("Ratio highest MPPC charge / second highest");
	varStream << "FirstOnSecond_" << k << ">>" << histoString ;
	varString = varStream.str();
	CanvasFirstOnSecond[k]->cd();
	t1->Draw(varString.c_str());
	histoStream.str(std::string());
	varStream.str(std::string());
	
	//not useful either
	histoStream << "recCutOnRatio-Module" << k;
	histoString = histoStream.str();
	recCutOnRatio[k] = new TH2F(histoString.c_str(),histoString.c_str(),Params.histo2DglobalBin,-7,7,Params.histo2DglobalBin,-7,7);
	recCutOnRatio[k]->GetXaxis()->SetTitle("Flood X");
	recCutOnRatio[k]->GetYaxis()->SetTitle("Flood Y");
	recCutOnRatio[k]->GetZaxis()->SetTitle("N");
	recCutOnRatio[k]->SetTitle("Reconstruction of dataset, cut on first/second > 1.7");
	varStream << "FloodY_" << k << ":FloodX_" << k << ">>" << histoString ;
	varString = varStream.str();
	cutStream << "FirstOnSecond_" << k << "> 1.7";
	cutString = cutStream.str();
	FloodCanvasCutOnRatio[k]->cd();
	t1->Draw(varString.c_str(),cutString.c_str(),"COLZ");
	histoStream.str(std::string());
	varStream.str(std::string());
	cutStream.str(std::string());
	
	histoStream << "recFinal-Module" << k;
	histoString = histoStream.str();
	recFinal[k] = new TH2F(histoString.c_str(),histoString.c_str(),Params.histo2DglobalBin,-7,7,Params.histo2DglobalBin,-7,7);
	
	//REAL DEAL
	//So now first thing, 32 flood histograms are generated, without any selection except for the requirement that for each MPPC only the 
	//events tagged as TriggerChannel in that MPPC are used. Everything is taken from the "zero cut" floodx and floody branches
	//in the ttree
	//Flood canvas, the initial flood with no cut. We divide it in 32 canvases, cause finding 4 peaks is already quite challenging...
	//f2->SetParameters(14,-2,0.3,0.76,0.2);
	
	//for each channel, a reconstruction flood is done (then in reality only 16 are filled and plotted, because only 
	//16 channels have data, but i leave them more generally like this
	//TH2F *recFlood[32];
	//strings and stuff
	std::stringstream recFloodStream[32], drawRecStream[32],drawCutStream[32],SpectraCanvasStream[32];
	std::string recFloodString[32], drawRecString[32],drawCutString[32],SpectraCanvasString[32];
	//TSpectrum2 is used to find the peak positions
	TSpectrum2 *s_peak[32];
	
// // 	//La mega stringa 
// // 	std::stringstream MegaStream;
// // 	std::string MegaString;

	//Peak position and sigma will be stored in these variables
	//Double_t mean[32],sigma[32];
	//loops on the TH1F array, fills them with data from the tree
	std::cout << "Filling histograms..." << std::endl;

	for(int i=0; i<32; i++)// run on all the channels 
	{
	  if(Input[i].Module[k])//if the channel belongs to this module
	  {
	    if(Input[i].Data) //do it only if the channel have data
	    {
	      std::cout << "Analyzing channel " << i << "..." << std::endl;
	      
 
	      //-------------------------------------------//
	      //                 ADC counts                //
	      //-------------------------------------------//
	      
	      //------------------(1)----------------------//
	      //plot the original spectra of the 16 mppc, without TriggerChannel selection
	      CanvasOriginal[k]->cd(Input[i].CanvasCd);//select which subcanvas is used for drawing this histo
	      std::stringstream channelOriginalStream,channelOriginalVarStream;
	      std::string channelOriginalString,channelOriginalVarString;
	      channelOriginalStream << "Original channel " << i << " - " << Input[i].LabelName << " - Module " << k;
	      channelOriginalString = channelOriginalStream.str();
	      channelOriginalVarStream << names[i] << " >> " << channelOriginalString;
	      channelOriginalVarString = channelOriginalVarStream.str();
	      Input[i].channelOriginal = new TH1F(channelOriginalString.c_str(),"",Params.binning,1,Params.histomax+1); //TH1F for this channel created here
	      Input[i].channelOriginal->SetTitle(channelOriginalString.c_str());
	      gPad-> SetLogy();
	      t1->Draw(channelOriginalVarString.c_str());
	      
	      //------------------(2)---------------------//
	      //now plot the same spectra, but using the "max charge" info -> TriggerChannel
	      char channelName[200]; 
	      char FillTheHisto[220]; 
	      char FillTheHistoTag[220]; 
	      std::stringstream MixedStream,MixedVarStream,MixedCutStream;
	      std::string MixedString,MixedVarString,MixedCutString;
	      sprintf(channelName,"Channel %d - %s Module %d - Spectrum",i,Input[i].LabelName.Data(),k);
	      sprintf(FillTheHisto,"%s>>%s",names[i].c_str(),channelName); 
	      sprintf(FillTheHistoTag,"TriggerChannel_%d == %d",k,i); 
	      Input[i].channel = new TH1F(channelName,"",Params.binning,1,Params.histomax+1);
	      Input[i].channel->SetTitle(channelName);
	      Input[i].channel->GetXaxis()->SetTitle("ADC Channels");
	      Input[i].channel->GetYaxis()->SetTitle("Counts");
	      Canvas1[k]->cd(Input[i].CanvasCd);//select which subcanvas is used for drawing this histo  
	      t1->Draw(FillTheHisto,FillTheHistoTag);//Draw the histo taking data from the tree. For example, the real command would look like 
	      //clone the channel histo (to be used in the mixed plot)
	      Input[i].channelClone = (TH1F*)Input[i].channel->Clone("clone"); //clones the histo
	      //     if (i!=ch) channel[i]->SetStats(0);
	      //     else channel[i]->SetStats(1);
	      //if the channel that is being filled in this step of the loop is the trigger 
	      //channel, fit its histogram and get the mean and sigma of the photopeak 
	      // 	if (Input[i].LabelName == TriggerChannel)
	      // 	{
	      // 	  TriggerChannelID = i;
	      Input[i].channel->SetStats(1);
	      //-------------------------------------------//
	      //             END of ADC counts             //
	      //-------------------------------------------//
	      
	      
	      //-------------------------------------------//
	      //               Pixels Fired                //
	      //-------------------------------------------//
	      //1)
	      //plot the trigger spectra of the 16 mppc, with TriggerChannel selection but rescaled to get pixels fired
	      PFCanvasOriginal[k]->cd(Input[i].CanvasCd);//select which subcanvas is used for drawing this histo
	      std::stringstream PFchannelOriginalStream,PFchannelOriginalStreamVar,PFchannelOriginalStreamCut;
	      std::string PFchannelOriginalString,PFchannelOriginalStringVar,PFchannelOriginalStringCut;
	      PFchannelOriginalStream << "MPPC " << Input[i].LabelName << " - Module " << k;
	      PFchannelOriginalString = PFchannelOriginalStream.str();
	      PFchannelOriginalStreamVar << ChtoPFscalingFactor << " * " << names[i] << " >> " << PFchannelOriginalString;
	      PFchannelOriginalStringVar = PFchannelOriginalStreamVar.str();
	      PFchannelOriginalStreamCut << "TriggerChannel_" << k << " == " << i;
	      PFchannelOriginalStringCut = PFchannelOriginalStreamCut.str();
	      Input[i].PFchannelOriginal = new TH1F(PFchannelOriginalString.c_str(),"",(int)((1.0/1.5)*Params.PFbinning),1,Params.PFhistomax+1); //TH1F for this channel created here
	      //Input[i].PFchannelOriginal->SetTitle(PFchannelOriginalString.c_str());
// 	      gPad-> SetLogy();
	      //Input[i].channelOriginal->SetLogy();
	      //t1->Draw(PFchannelOriginalVarString.c_str());
// 	      for(int iBin=0;iBin<Input[i].channelOriginal->GetNbinsX();iBin++)
// 	      {
// 		Input[i].PFchannelOriginal->SetBinContent(iBin,Input[i].channelOriginal->GetBinContent(iBin));
// 	      }
	      Input[i].PFchannelOriginal->SetTitle(PFchannelOriginalString.c_str());
	      Input[i].PFchannelOriginal->GetXaxis()->SetTitle("N of Pixels Fired");
	      Input[i].PFchannelOriginal->GetYaxis()->SetTitle("Counts");
	      Input[i].PFchannelOriginal->GetXaxis()->SetLabelSize(0.06);
	      Input[i].PFchannelOriginal->GetYaxis()->SetLabelSize(0.06);
	      Input[i].PFchannelOriginal->GetXaxis()->SetNdivisions(505);
	      Input[i].PFchannelOriginal->GetYaxis()->SetNdivisions(505);
	      Input[i].PFchannelOriginal->GetYaxis()->SetTitleSize(0.06);
	      Input[i].PFchannelOriginal->GetYaxis()->SetTitleOffset(0.7);
	      Input[i].PFchannelOriginal->GetXaxis()->SetTitleSize(0.06);
	      Input[i].PFchannelOriginal->GetXaxis()->SetTitleOffset(0.7);
	      //std::cout << PFchannelOriginalStringVar << std::endl;
	      //std::cout << PFchannelOriginalStringCut << std::endl;
	      t1->Draw(PFchannelOriginalStringVar.c_str(),PFchannelOriginalStringCut.c_str());
	      //2)
	      //plot the SUM spectra of the 16 mppc, with TriggerChannel selection and rescaled to pixel fired
	      PFCanvas1[k]->cd(Input[i].CanvasCd);//select which subcanvas is used for drawing this histo
	      std::stringstream PFchannelStream,PFchannelStreamVar,PFchannelStreamCut;
	      std::string PFchannelString,PFchannelStringVar,PFchannelStringCut,tempHistoPF;
	      //tempHistoPF = "temporary histogram";
	      PFchannelStream << "Sum Pixels Fired for Trigger in " << Input[i].LabelName << " - Module " << k;
	      PFchannelString = PFchannelStream.str();	      
	      //creates a string for the variables to draw, adding only the active channels
	      int PFIsNotFirst = 0;
	      //we add at the beginning the ChtoPFscalingFactor
	      PFchannelStreamVar << ChtoPFscalingFactor <<  " * ("; //
	      for(int PFIndex = 0 ; PFIndex < 32 ; PFIndex++)
	      {
		if(Input[PFIndex].Module[k])//if the channel belongs to this module
		{
		  if(Input[PFIndex].Data)//if channel has data
		  {
		    if(PFIsNotFirst) PFchannelStreamVar << "+";
		    PFIsNotFirst++;
		    PFchannelStreamVar << "ch" << PFIndex;
		  }
		}
	      }
	      PFchannelStreamVar <<  ") >> " << PFchannelString;
	      PFchannelStringVar = PFchannelStreamVar.str(); //useful later 
	      //fill this histogram only if the TriggerChannel is indentified as this channel
	      PFchannelStreamCut << "TriggerChannel_" << k << " == " << i;
	      PFchannelStringCut = PFchannelStreamCut.str();
	      //PFchannelOriginalVarString = PFchannelOriginalVarStream.str();
	      //Input[i].PFchannelTemp = new TH1F(tempHistoPF.c_str(),"",Params.binning,0,Params.histomax); //TH1F for this channel created here
	      Input[i].PFchannel = new TH1F(PFchannelString.c_str(),"",Params.PFbinning,1,1.5*Params.PFhistomax+1); //TH1F for this channel created here
	      //Input[i].PFchannelOriginal->SetTitle(PFchannelOriginalString.c_str());
	      //gPad-> SetLogy();
	      //Input[i].channelOriginal->SetLogy();
	      //t1->Draw(PFchannelOriginalVarString.c_str());
// 	      for(int iBin=0;iBin<Input[i].channel->GetNbinsX();iBin++)
// 	      {
// 		Input[i].PFchannel->SetBinContent(iBin,Input[i].PFchannelTemp->GetBinContent(iBin));
// 	      }
	      Input[i].PFchannel->SetTitle(PFchannelString.c_str());
	      Input[i].PFchannel->GetXaxis()->SetTitle("N of Pixels Fired");
	      Input[i].PFchannel->GetYaxis()->SetTitle("Counts");
	      Input[i].PFchannel->GetXaxis()->SetLabelSize(0.06);
	      Input[i].PFchannel->GetYaxis()->SetLabelSize(0.06);
	      Input[i].PFchannel->GetXaxis()->SetNdivisions(505);
	      Input[i].PFchannel->GetYaxis()->SetNdivisions(505);
	      Input[i].PFchannel->GetYaxis()->SetTitleSize(0.06);
	      Input[i].PFchannel->GetYaxis()->SetTitleOffset(0.7);
	      Input[i].PFchannel->GetXaxis()->SetTitleSize(0.06);
	      Input[i].PFchannel->GetXaxis()->SetTitleOffset(0.7);
	      //Input[i].PFchannel->Draw();
	      //std::cout << PFchannelStringVar << std::endl;
	      //std::cout << PFchannelStringCut << std::endl;
	      t1->Draw(PFchannelStringVar.c_str(),PFchannelStringCut.c_str());
	      //-------------------------------------------//
	      //            END of Pixels Fired            //
	      //-------------------------------------------//
	      
	      
	      
	      //-------------------------------------------//
	      //                   Photons                 //
	      //-------------------------------------------//
	      //1)
	      //plot the trigger spectra of the 16 mppc, with TriggerChannel selection but rescaled to get LY
	      LYCanvasOriginal[k]->cd(Input[i].CanvasCd);//select which subcanvas is used for drawing this histo
	      std::stringstream LYchannelOriginalStream,LYchannelOriginalStreamVar,LYchannelOriginalStreamCut;
	      std::string LYchannelOriginalString,LYchannelOriginalStringVar,LYchannelOriginalStringCut;
	      LYchannelOriginalStream << "Spectrum of MPPC " << Input[i].LabelName << " - Module " << k;
	      LYchannelOriginalString = LYchannelOriginalStream.str();
	      LYchannelOriginalStreamVar << ChtoLYscalingFactor << " * " << names[i] << " >> " << LYchannelOriginalString;
	      LYchannelOriginalStringVar = LYchannelOriginalStreamVar.str();
	      LYchannelOriginalStreamCut << "TriggerChannel_" << k << " == " << i;
	      LYchannelOriginalStringCut = LYchannelOriginalStreamCut.str();
	      Input[i].LYchannelOriginal = new TH1F(LYchannelOriginalString.c_str(),"",(int)((1.0/1.5)*Params.LYbinning),1,Params.LYhistomax+1); //TH1F for this channel created here
	      //Input[i].LYchannelOriginal->SetTitle(LYchannelOriginalString.c_str());
// 	      gPad-> SetLogy();
	      //Input[i].channelOriginal->SetLogy();
	      //t1->Draw(LYchannelOriginalVarString.c_str());
// 	      for(int iBin=0;iBin<Input[i].channelOriginal->GetNbinsX();iBin++)
// 	      {
// 		Input[i].LYchannelOriginal->SetBinContent(iBin,Input[i].channelOriginal->GetBinContent(iBin));
// 	      }
	      Input[i].LYchannelOriginal->SetTitle(LYchannelOriginalString.c_str());
	      Input[i].LYchannelOriginal->GetXaxis()->SetTitle("Ph/MeV");
	      Input[i].LYchannelOriginal->GetYaxis()->SetTitle("Counts");
	      Input[i].LYchannelOriginal->GetXaxis()->SetLabelSize(0.06);
	      Input[i].LYchannelOriginal->GetYaxis()->SetLabelSize(0.06);
	      Input[i].LYchannelOriginal->GetXaxis()->SetNdivisions(505);
	      Input[i].LYchannelOriginal->GetYaxis()->SetNdivisions(505);
	      Input[i].LYchannelOriginal->GetYaxis()->SetTitleSize(0.06);
	      Input[i].LYchannelOriginal->GetYaxis()->SetTitleOffset(0.7);
	      Input[i].LYchannelOriginal->GetXaxis()->SetTitleSize(0.06);
	      Input[i].LYchannelOriginal->GetXaxis()->SetTitleOffset(0.7);
	      //std::cout << LYchannelOriginalStringVar << std::endl;
	      //std::cout << LYchannelOriginalStringCut << std::endl;
	      t1->Draw(LYchannelOriginalStringVar.c_str(),LYchannelOriginalStringCut.c_str());
	      //and now let's look for a photopeak. It makes sense only for 4x4 matrices, but in the 8x8 case it will just not work, no big deal..
	      TSpectrum *sLY;
	      sLY = new TSpectrum(5);
	      //Input[i].crystal[j]->Smooth(10);
	      //Input[i].crystal[j]->Rebin(2);
	      //Input[i].SpectraCanvas->cd(j+1);
	      Int_t LYCrystalPeaksN = sLY->Search(Input[i].LYchannelOriginal,1,"goff",0.5); //TODO pass to "goff"
	      Float_t *LYCrystalPeaks = sLY->GetPositionX();
	      Float_t *LYCrystalPeaksY = sLY->GetPositionY();
	      float LYmaxPeak = 0.0;
	      int LYpeakID = 0;
	      for (int LYpeakCounter = 0 ; LYpeakCounter < LYCrystalPeaksN ; LYpeakCounter++ )
	      {
		if(LYCrystalPeaks[LYpeakCounter] > LYmaxPeak)
		{
		  LYmaxPeak = LYCrystalPeaks[LYpeakCounter];
		  LYpeakID = LYpeakCounter;
		}
	      }
	      float LYenergyResolution;
	      if (correctingForSaturation)
		LYenergyResolution = ENERGY_RESOLUTION_SATURATION_CORRECTION;
	      else
		LYenergyResolution = ENERGY_RESOLUTION;
	      float LYpar0 = LYCrystalPeaksY[LYpeakID];
	      float LYpar1 = LYCrystalPeaks[LYpeakID];
	      float LYpar2 = (LYCrystalPeaks[LYpeakID]*LYenergyResolution)/2.35;
	      //std::cout << "Crystal " << BigCanvasMapping[crystalCounter[k]] << " = " << par0 << " " << par1 << " " << par2 << std::endl;
	      LYgauss->SetParameter(0,LYpar0);
	      LYgauss->SetParameter(1,LYpar1);
	      LYgauss->SetParameter(2,LYpar2); //expected FWHM en res = 12%
	      Input[i].LYchannelOriginal->Fit("LYgauss","QN","",LYpar1-1.2*LYpar2,LYpar1+2.0*LYpar2);
	      //store the mean and sigma in the data struct
	      Input[i].LightSharingCrystalMean = LYgauss->GetParameter(1);
	      Input[i].LightSharingCrystalSigma = std::abs(LYgauss->GetParameter(2));
	      //2)
	      //plot the SUM spectra of the 16 mppc, with TriggerChannel selection and rescaled to LY
	      LYCanvas1[k]->cd(Input[i].CanvasCd);//select which subcanvas is used for drawing this histo
	      std::stringstream LYchannelStream,LYchannelStreamVar,LYchannelStreamCut;
	      std::string LYchannelString,LYchannelStringVar,LYchannelStringCut,tempHistoLY;
	      //tempHistoLY = "temporary histogram";
	      LYchannelStream << "Sum Spectrum for Trigger in " << Input[i].LabelName << " - Module " << k;
	      LYchannelString = LYchannelStream.str();	      
	      //creates a string for the variables to draw, adding only the active channels
	      int LYIsNotFirst = 0;
	      //we add at the beginning the ChtoLYscalingFactor
	      LYchannelStreamVar << ChtoLYscalingFactor <<  " * ("; //
	      for(int LYIndex = 0 ; LYIndex < 32 ; LYIndex++)
	      {
		if(Input[LYIndex].Module[k])//if the channel belongs to this module
		{
		  if(Input[LYIndex].Data)//if channel has data
		  {
		    if(LYIsNotFirst) LYchannelStreamVar << "+";
		    LYIsNotFirst++;
		    LYchannelStreamVar << "ch" << LYIndex;
		  }
		}
	      }
	      LYchannelStreamVar <<  ") >> " << LYchannelString;
	      LYchannelStringVar = LYchannelStreamVar.str(); //useful later 
	      //fill this histogram only if the TriggerChannel is indentified as this channel
	      LYchannelStreamCut << "TriggerChannel_" << k << " == " << i;
	      LYchannelStringCut = LYchannelStreamCut.str();
	      //LYchannelOriginalVarString = LYchannelOriginalVarStream.str();
	      //Input[i].LYchannelTemp = new TH1F(tempHistoLY.c_str(),"",Params.binning,0,Params.histomax); //TH1F for this channel created here
	      Input[i].LYchannel = new TH1F(LYchannelString.c_str(),"",2*Params.LYbinning,1,1.5*Params.LYhistomax+1); //TH1F for this channel created here
	      //Input[i].LYchannelOriginal->SetTitle(LYchannelOriginalString.c_str());
	      //gPad-> SetLogy();
	      //Input[i].channelOriginal->SetLogy();
	      //t1->Draw(LYchannelOriginalVarString.c_str());
// 	      for(int iBin=0;iBin<Input[i].channel->GetNbinsX();iBin++)
// 	      {
// 		Input[i].LYchannel->SetBinContent(iBin,Input[i].LYchannelTemp->GetBinContent(iBin));
// 	      }
	      Input[i].LYchannel->SetTitle(LYchannelString.c_str());
	      Input[i].LYchannel->GetXaxis()->SetTitle("Ph/MeV");
	      Input[i].LYchannel->GetYaxis()->SetTitle("Counts");
	      Input[i].LYchannel->GetXaxis()->SetLabelSize(0.06);
	      Input[i].LYchannel->GetYaxis()->SetLabelSize(0.06);
	      Input[i].LYchannel->GetXaxis()->SetNdivisions(505);
	      Input[i].LYchannel->GetYaxis()->SetNdivisions(505);
	      Input[i].LYchannel->GetYaxis()->SetTitleSize(0.06);
	      Input[i].LYchannel->GetYaxis()->SetTitleOffset(0.7);
	      Input[i].LYchannel->GetXaxis()->SetTitleSize(0.06);
	      Input[i].LYchannel->GetXaxis()->SetTitleOffset(0.7);
	      //Input[i].LYchannel>-Draw();
	      //std::cout << LYchannelStringVar << std::endl;
	      //std::cout << LYchannelStringCut << std::endl;
	      t1->Draw(LYchannelStringVar.c_str(),LYchannelStringCut.c_str());
	      //find the peak for the sum plots
	      TSpectrum *sSumLY;
	      sSumLY = new TSpectrum(5);
	      
	      Int_t SumLYPeaks = sSumLY->Search(Input[i].LYchannel,1,"goff",0.5); //TODO pass to "goff"
	      Float_t *LYSumPeaks = sSumLY->GetPositionX();
	      Float_t *LYSumPeaksY = sSumLY->GetPositionY();
	      
	      LYmaxPeak = 0.0;
	      LYpeakID = 0;
	      for (int LYpeakCounter = 0 ; LYpeakCounter < SumLYPeaks ; LYpeakCounter++ )
	      {
		if(LYSumPeaks[LYpeakCounter] > LYmaxPeak)
		{
		  LYmaxPeak = LYSumPeaks[LYpeakCounter];
		  LYpeakID = LYpeakCounter;
		}
	      }
	      float LYSumpar0 = LYSumPeaksY[LYpeakID];
	      float LYSumpar1 = LYSumPeaks[LYpeakID];
	      float LYSumpar2 = (LYSumPeaks[LYpeakID]*LYenergyResolution)/2.35;
	      //std::cout << "Crystal " << BigCanvasMapping[crystalCounter[k]] << " = " << par0 << " " << par1 << " " << par2 << std::endl;
	      LYSumgauss->SetParameter(0,LYSumpar0);
	      LYSumgauss->SetParameter(1,LYSumpar1);
	      LYSumgauss->SetParameter(2,LYSumpar2); //expected FWHM en res = 12%
	      Input[i].LYchannel->Fit("LYSumgauss","Q","",LYSumpar1-3*LYSumpar2,LYSumpar1+2.0*LYSumpar2);
	      //store the mean and sigma in the data struct
	      Input[i].SumLYMean = LYSumgauss->GetParameter(1);
	      Input[i].SumLYSigma = std::abs(LYSumgauss->GetParameter(2));
	      
	      
	      
	      //-------------------------------------------//
	      //             END of Photons                //
	      //-------------------------------------------//
	      
	      
	      
	      //now, knowing which channel is the trigger, fill the other spectra only with events that are trigger in TriggerChannel
	      //it means you need to run on the 16-32 channels again..
	      std::stringstream lightShStream;
	      std::string lightShString;
	      lightShStream << "Channel " << i << " - " << Input[i].LabelName << " - Light Sharing " << " - Module " << k;
	      lightShString = lightShStream.str();
	      Input[i].LightSharingCanvas = new TCanvas(lightShString.c_str(),"",1200,800);
	      Input[i].LightSharingCanvas->Divide(4,4);
	      for (int kTag = 0; kTag < 32 ; kTag++)
	      {
		if(Input[kTag].Module[k])//if the channel belongs to this module
		{
		  if(Input[kTag].Data)//and if the channel is on
		  {
		    std::stringstream spectraTagStream,spectraTagVarStream,spectraTagCutStream;
		    std::string spectraTagString,spectraTagVarString,spectraTagCutString;
		    spectraTagStream << "Channel " << kTag << " - " << Input[kTag].LabelName << " - Trigger on " << i << " - Module " << k;
		    spectraTagString = spectraTagStream.str();
		    
		    spectraTagVarStream << names[kTag] << " >> " << spectraTagString;
		    spectraTagVarString = spectraTagVarStream.str();
		    
		    spectraTagCutStream << "TriggerChannel_" << k << " == " << i;
		    spectraTagCutString = spectraTagCutStream.str();
		    
		    Input[i].spectraTag[kTag] = new TH1F(spectraTagString.c_str(),spectraTagString.c_str(),Params.binning,1,Params.histomax+1);
		    Input[i].LightSharingCanvas->cd(Input[kTag].CanvasCd);
		    t1->Draw(spectraTagVarString.c_str(),spectraTagCutString.c_str());
		  }
		}
	      }
	      
	      
	      //the same plot, but cutting on the "photopeak", i.e. the standard procedure that makes no sense for most of the channels when the crystal are 64
	      std::stringstream lightShStreamPE;
	      std::string lightShStringPE;
	      lightShStreamPE << "Channel " << i << " - " << Input[i].LabelName << " - Light Sharing PE" << " - Module " << k;
	      lightShStringPE = lightShStreamPE.str();
	      Input[i].LightSharingCanvasPE = new TCanvas(lightShStringPE.c_str(),"",1200,800);
	      Input[i].LightSharingCanvasPE->Divide(4,4);
	      for (int kTag = 0; kTag < 32 ; kTag++)
	      {
		if(Input[kTag].Module[k])//if the channel belongs to this module
		{
		  if(Input[kTag].Data)//and if the channel is on
		  {
		    std::stringstream spectraTagStreamPE,spectraTagVarStreamPE,spectraTagCutStreamPE;
		    std::string spectraTagStringPE,spectraTagVarStringPE,spectraTagCutStringPE;
		    spectraTagStreamPE << "Channel " << kTag << " - " << Input[kTag].LabelName << " - Trigger on " << i << " - Module " << k << "PE CUT";
		    spectraTagStringPE = spectraTagStreamPE.str();
		    spectraTagVarStreamPE << names[kTag] << " >> " << spectraTagStringPE;
		    spectraTagVarStringPE = spectraTagVarStreamPE.str();
		    spectraTagCutStreamPE << "TriggerChannel_" << k << " == " << i << "&& ch" << i << ">" << Params.PEmin << "&& ch" << i << "<" << Params.PEmax;
		    spectraTagCutStringPE = spectraTagCutStreamPE.str();
		    Input[i].spectraTagPE[kTag] = new TH1F(spectraTagStringPE.c_str(),spectraTagStringPE.c_str(),Params.binning,1,Params.histomax+1);
		    Input[i].LightSharingCanvasPE->cd(Input[kTag].CanvasCd);
		    t1->Draw(spectraTagVarStringPE.c_str(),spectraTagCutStringPE.c_str());
		  }
		}
	      }
	      
	      
	      //the same plot again, but in LY scale and cutting on the real photopeak (from fit)
	      std::stringstream LYlightShStreamPE;
	      std::string LYlightShStringPE;
	      LYlightShStreamPE << "LY Sharing - Channel " << i << " - " << Input[i].LabelName << " - Light Sharing LY" << " - Module " << k;
	      LYlightShStringPE = LYlightShStreamPE.str();
	      Input[i].LYLightSharingCanvasPE = new TCanvas(LYlightShStringPE.c_str(),"",1200,800);
	      Input[i].LYLightSharingCanvasPE->Divide(4,4);
	      for (int kTag = 0; kTag < 32 ; kTag++)
	      {
		if(Input[kTag].Module[k])//if the channel belongs to this module
		{
		  if(Input[kTag].Data)//and if the channel is on
		  {
		    std::stringstream LYspectraTagStreamPE,LYspectraTagVarStreamPE,LYspectraTagCutStreamPE;
		    std::string LYspectraTagStringPE,LYspectraTagVarStringPE,LYspectraTagCutStringPE;
		    if(kTag == i)
		      LYspectraTagStreamPE << "Light in Trigger Channel " << Input[i].LabelName << " - Module " << k;
		    else
		      LYspectraTagStreamPE << "Light Shared to Channel " << Input[kTag].LabelName << " - Trigger on " << Input[i].LabelName << " - Module " << k;
		    LYspectraTagStringPE = LYspectraTagStreamPE.str();
		    LYspectraTagVarStreamPE << ChtoLYscalingFactor << " * " << names[kTag] << " >> " << LYspectraTagStringPE;
		    LYspectraTagVarStringPE = LYspectraTagVarStreamPE.str();
		    LYspectraTagCutStreamPE << "TriggerChannel_" << k << " == " << i << "&& (( ch" << i << " * " << ChtoLYscalingFactor << " ) >" << Input[i].LightSharingCrystalMean - 1.5*Input[i].LightSharingCrystalSigma << " ) && ((ch" << i << "*"<< ChtoLYscalingFactor << ") <" << Input[i].LightSharingCrystalMean + 2.5*Input[i].LightSharingCrystalSigma << ")";
		    LYspectraTagCutStringPE = LYspectraTagCutStreamPE.str() ;
		    //std::cout << LYspectraTagVarStringPE << std::endl;
		    //std::cout << LYspectraTagCutStringPE << std::endl;
		    Input[i].LYspectraTagPE[kTag] = new TH1F(LYspectraTagStringPE.c_str(),LYspectraTagStringPE.c_str(),Params.binning*4,1,Params.LYhistomax+1);
		    Input[i].LYLightSharingCanvasPE->cd(Input[kTag].CanvasCd);
		    //Input[i].LYspectraTagPE[kTag]->SetTitle("Ph/MeV");
		    Input[i].LYspectraTagPE[kTag]->GetXaxis()->SetTitle("Ph/MeV");
		    Input[i].LYspectraTagPE[kTag]->GetYaxis()->SetTitle("Counts");
		    Input[i].LYspectraTagPE[kTag]->GetXaxis()->SetLabelSize(0.06);
		    Input[i].LYspectraTagPE[kTag]->GetYaxis()->SetLabelSize(0.06);
		    Input[i].LYspectraTagPE[kTag]->GetXaxis()->SetNdivisions(505);
		    Input[i].LYspectraTagPE[kTag]->GetYaxis()->SetNdivisions(505);
		    Input[i].LYspectraTagPE[kTag]->GetYaxis()->SetTitleSize(0.06);
		    Input[i].LYspectraTagPE[kTag]->GetYaxis()->SetTitleOffset(0.7);
		    Input[i].LYspectraTagPE[kTag]->GetXaxis()->SetTitleSize(0.06);
		    Input[i].LYspectraTagPE[kTag]->GetXaxis()->SetTitleOffset(0.7);
		    t1->Draw(LYspectraTagVarStringPE.c_str(),LYspectraTagCutStringPE.c_str());
		    if(Params.LightSharingAnalysis)
		    {
		      //now we fit each one of these spectra, with a simple gauss function, to find the peak positions
		      //and we store them in the array
// 		      TSpectrum *sLS;
// 		      sLS = new TSpectrum(1);
// 		      //Input[i].crystal[j]->Smooth(10);
// 		      //Input[i].crystal[j]->Rebin(2);
// 		      //Input[i].SpectraCanvas->cd(j+1);
// 		      Int_t LSCrystalPeaksN = sLS->Search(Input[i].LYspectraTagPE[kTag],1,"",0.2); //TODO pass to "goff"
// 		      Float_t *LSCrystalPeaks = sLS->GetPositionX();
// 		      Float_t *LSCrystalPeaksY = sLS->GetPositionY();
// 		      float LSmaxPeak = 0.0;
// 		      int LSpeakID = 0;
// 		      for (int LSpeakCounter = 0 ; LSpeakCounter < LSCrystalPeaksN ; LSpeakCounter++ )
// 		      {
// 			if(LSCrystalPeaks[LSpeakCounter] > LSmaxPeak)
// 			{
// 			  LSmaxPeak = LSCrystalPeaks[LSpeakCounter];
// 			  LSpeakID = LSpeakCounter;
// 			}
// 		      }
// 		      float LSenergyResolution;
// 		      if (correctingForSaturation)
// 			LSenergyResolution = ENERGY_RESOLUTION_SATURATION_CORRECTION;
// 		      else
// 			LSenergyResolution = ENERGY_RESOLUTION;
// 		      float LSpar0 = LSCrystalPeaksY[LSpeakID];
// 		      float LSpar1 = LSCrystalPeaks[LSpeakID];
// 		      float LSpar2 = (LSCrystalPeaks[LSpeakID]*LSenergyResolution)/2.35;
// 		      //std::cout << "Crystal " << BigCanvasMapping[crystalCounter[k]] << " = " << par0 << " " << par1 << " " << par2 << std::endl;
// 		      LSgauss->SetParameter(0,LSpar0);
// 		      LSgauss->SetParameter(1,LSpar1);
// 		      LSgauss->SetParameter(2,LSpar2); //expected FWHM en res = 12%
// 		      //Input[i].LYchannelOriginal->Fit("LYgauss","Q","",LYpar1-1.2*LYpar2,LYpar1+2.0*LYpar2);
// 		      Input[i].LYspectraTagPE[kTag]->Fit("LSgauss","","",LSpar1-3*LSpar2,LSpar1+2.0*LSpar2);
// 		      Input[i].lightSharingData[kTag] = LSgauss->GetParameter(1);
		      //let's trust the mean of the histo...
		      //Input[i].lightSharingData[kTag] = Input[i].LYspectraTagPE[kTag]->GetMean();
		      //std::cout << kTag << "\t" << Input[i].lightSharingData[kTag] << std::endl;
		      
		      //make it simple: they should be single peaks, so let's just take the position of the max as fit starting point
		      float LSenergyResolution;
		      if (correctingForSaturation)
			LSenergyResolution = ENERGY_RESOLUTION_SATURATION_CORRECTION;
		      else
			LSenergyResolution = ENERGY_RESOLUTION;
		      int binmax = Input[i].LYspectraTagPE[kTag]->GetMaximumBin();
                      double maxPos = Input[i].LYspectraTagPE[kTag]->GetXaxis()->GetBinCenter(binmax);
		      
		      float LSpar1 = maxPos;
// 		      float LSpar1 = LSCrystalPeaks[LSpeakID];
		      float LSpar2 = (maxPos*LSenergyResolution)/2.35;
// 		      LSgauss->SetParameter(0,LSpar0);
		      LSgauss->SetParameter(1,LSpar1);
		      LSgauss->SetParameter(2,LSpar2); //expected FWHM en res = 12%
		      //double LSfitMin = 0;
		      //if (LSpar1-3*LSpar2 > 0 )
			//LSfitMin = LSpar1-3*LSpar2;
		      //double LSfitMax = LSpar1+2.0*LSpar2;
		      Input[i].LYspectraTagPE[kTag]->Fit("LSgauss","QR");
		      std::string converged = "CONVERGED";
		      std::string status = (std::string) gMinuit->fCstatu;
		      //std::cout << "STATUS = " << status << std::endl;
		      //std::cout << "COMPARE = " << status.compare(converged) << std::endl;
		      if(status.compare(converged) == 1) // why 1????? bloody ROOT..
		      {
			//std::cout << "Entering cycle" << std::endl;
			if(LSgauss->GetParameter(1) >= 0) //no negative light sharing please 
			  Input[i].lightSharingData[kTag] = LSgauss->GetParameter(1);
			else
			  Input[i].lightSharingData[kTag] = 0;
		      }
		      else
		      {
			Input[i].lightSharingData[kTag] = maxPos;
		      }
		      
		    }
		  }
		}
	      }
	      
	      
	      //generate a flood histogram for this channel
	      //a modification was needed to adapt to LIP output: they get data only for the channels above a fixed threshold, this results in 
	      //a number of channels with data != 0 that is generally < 16. When the channels with data are really a few, the anger logic reco 
	      //will "accumulate" in lines corresponding to the x and y positions of the mppcs. This clearly happens when the 16 channels for a given
	      //event have data only in a column or in a row (or worse, only in the trigger channel, which would result in reconstructing the event in
	      //the real center of the mppc. This will prevent the algorithm to find the 4 peaks (or maybe it can be made in such a way that the 4 peaks
	      //are found anyway, but honestly i don't have any idea how, for the moment). So we have to set a rule to avoid this events to even enter 
	      //our scatter plot
	      recFloodStream[i] << "recFlood" << i << "_Module_" << k;
	      recFloodString[i] = recFloodStream[i].str();								//histogram name
	      Input[i].recFlood = new TH2F(recFloodString[i].c_str(),recFloodString[i].c_str(),Params.histo2DchannelBin,-7,7,Params.histo2DchannelBin,-7,7);		//create it, with 250x250 bins
	      Input[i].recFlood->GetXaxis()->SetTitle("FloodX");								//axis
	      Input[i].recFlood->GetYaxis()->SetTitle("FloodY");								//axis
	      Input[i].recFlood->GetZaxis()->SetTitle("N");									//axis
	      //draw vars
	      drawRecStream[i] << "FloodY_" << k << ":FloodX_" << k << " >> " << recFloodString[i];
	      drawRecString[i] = drawRecStream[i].str();								//variable to pass to the Draw
	      drawCutStream[i] << "TriggerChannel_" << k << " == " << i;						//origina cut to pass to the Draw
	      //now the cuts needed to adapt to LIP output FIXME hardcoded numbers, can be made better with a cycle
	      if(Params.CutGrid)
		drawCutStream[i] << " && " << lipCutString[k];
	      drawCutString[i] = drawCutStream[i].str();								
	      //canvas choice      
	      FloodCanvas[k]->cd(Input[i].CanvasCd);									//select which subcanvas is used for drawing this histo  
	      //Draw command
	      t1->Draw(drawRecString[i].c_str(),drawCutString[i].c_str(),"COLZ");
	      
	      //then write then same in a canvas, so it's nicer when saved separately in the root file..
	      std::string C_name = "Canvas" + recFloodString[i]; 
	      Input[i].recChannelCanvas = new TCanvas(C_name.c_str(),"",800,800);
	      std::stringstream C_varStream; 
	      C_varStream  << "FloodY_" << k << ":FloodX_" << k << " >> " << C_name; 
	      std::string C_var;
	      C_var = C_varStream.str();
	      t1->Draw(C_var.c_str(),drawCutString[i].c_str(),"COLZ");
	      
	      //peak search -FIXME make it work properly...
	      s_peak[i] = new TSpectrum2(4,1);
	      Int_t nfound = s_peak[i]->Search(Input[i].recFlood,1,"col",0.3);	
	      
	      
	      Double_t tempXsigma[4] = {0.2,0.2,0.2,0.2};
	      Double_t tempYsigma[4] = {0.2,0.2,0.2,0.2};
	      //store what has been found in the data struct
	      Input[i].xpeaks = s_peak[i]->GetPositionX();
	      Input[i].ypeaks = s_peak[i]->GetPositionY();
	      
	      //temporary sigmas
	      Input[i].xsigma = &tempXsigma[0];
	      Input[i].ysigma = &tempYsigma[0];
	      // 	{
	      //real fit of the 2d peaks
	      //std::cout << "Ch " << i << std::endl;
	      for ( int j = 0 ; j < nfound ; j++)
	      {
		
		Input[i].f2 = new TF2("f2","[0]*TMath::Gaus(x,[1],[2])*TMath::Gaus(y,[3],[4])",Input[i].xpeaks[j]-2.0*Input[i].xsigma[j],Input[i].xpeaks[j]+2.0*Input[i].xsigma[j],Input[i].ypeaks[j]-2.0*Input[i].ysigma[j],Input[i].ypeaks[j]+2.0*Input[i].ysigma[j]);
		Input[i].f2->SetParameters(14,Input[i].xpeaks[j],Input[i].xsigma[j],Input[i].ypeaks[j],Input[i].ysigma[j]);
		Input[i].recFlood->Fit("f2","QNR");
		//Input[i].f2->Draw("cont1 same");
// 		Input[i].fit2DmeanX[j] = Input[i].f2->GetParameter(1);
// 		Input[i].fit2DmeanY[j] = Input[i].f2->GetParameter(3);
		//make it work by kicking its ass... - FIXME
		Input[i].fit2DmeanX[j] = Input[i].xpeaks[j];
		Input[i].fit2DmeanY[j] = Input[i].ypeaks[j];
		//we insert the 1.5 sigma limit here, it makes things easier
		Input[i].fit2DsigmaX[j] = 2.5*std::abs(Input[i].f2->GetParameter(2));
		Input[i].fit2DsigmaY[j] = 2.5*std::abs(Input[i].f2->GetParameter(4));
		
		//std::cout  << Input[i].xpeaks[j] << " " << Input[i].ypeaks[j] << " " <<  Input[i].fit2DmeanX[j] << " " << Input[i].fit2DmeanY[j ]<< " " << Input[i].fit2DsigmaX[j] << " " << Input[i].fit2DsigmaY[j] << std::endl;
	      }
	      // 	}
	      //now, sort them to give a spatial meaning to this fitting
	      //the idea is
	      //we sort them from the highest y to the lowest, then 
	      //we check first and second in term of x value, and sort them from lowest to highest
	      //and the same for third and fourth
	      //in this way, they should be ordered in such a way that, when running on j, they are 
	      //mapped in the correct subcanvas in the BigSpectraCanvas 
	      float swapMeanX,swapMeanY;
	      float swapsigmaX,swapSigmaY;
	      //first, order from highest y to lowest
	      for ( int jSwap = 0 ; jSwap < 4 - 1 ; jSwap++)
	      {
		for(int kSwap = (jSwap+1); kSwap < 4; kSwap++)   // rest of the elements
		{
		  if (Input[i].fit2DmeanY[jSwap] < Input[i].fit2DmeanY[kSwap])          // descending order
		  {
		    swapMeanX 		= Input[i].fit2DmeanX[jSwap];
		    swapMeanY		= Input[i].fit2DmeanY[jSwap];
		    swapsigmaX		= Input[i].fit2DsigmaX[jSwap];
		    swapSigmaY		= Input[i].fit2DsigmaY[jSwap];
		    
		    Input[i].fit2DmeanX[jSwap]	 =   Input[i].fit2DmeanX[kSwap];	
		    Input[i].fit2DmeanY[jSwap]	 =   Input[i].fit2DmeanY[kSwap];	
		    Input[i].fit2DsigmaX[jSwap]	 =   Input[i].fit2DsigmaX[kSwap];
		    Input[i].fit2DsigmaY[jSwap]	 =   Input[i].fit2DsigmaY[kSwap];
		    
		    Input[i].fit2DmeanX[kSwap]	 =   swapMeanX ;
		    Input[i].fit2DmeanY[kSwap]	 =   swapMeanY	;
		    Input[i].fit2DsigmaX[kSwap]	 =   swapsigmaX ;
		    Input[i].fit2DsigmaY[kSwap]	 =   swapSigmaY ;
		  }
		}
	      }
	      //now, first and second, on x
	      if (Input[i].fit2DmeanX[0] > Input[i].fit2DmeanX[1])          // ascending order
	      {
		swapMeanX 		= Input[i].fit2DmeanX[0];
		swapMeanY		= Input[i].fit2DmeanY[0];
		swapsigmaX		= Input[i].fit2DsigmaX[0];
		swapSigmaY		= Input[i].fit2DsigmaY[0];
		
		Input[i].fit2DmeanX[0]	 =   Input[i].fit2DmeanX[1];	
		Input[i].fit2DmeanY[0]	 =   Input[i].fit2DmeanY[1];	
		Input[i].fit2DsigmaX[0]	 =   Input[i].fit2DsigmaX[1];
		Input[i].fit2DsigmaY[0]	 =   Input[i].fit2DsigmaY[1];
		
		Input[i].fit2DmeanX[1]	 =   swapMeanX ;
		Input[i].fit2DmeanY[1]	 =   swapMeanY	;
		Input[i].fit2DsigmaX[1]	 =   swapsigmaX ;
		Input[i].fit2DsigmaY[1]	 =   swapSigmaY ;
	      }
	      //finally, third and fourth, on ascending x
	      if (Input[i].fit2DmeanX[2] > Input[i].fit2DmeanX[3])          // ascending order
	      {
		swapMeanX 		= Input[i].fit2DmeanX[2];
		swapMeanY		= Input[i].fit2DmeanY[2];
		swapsigmaX		= Input[i].fit2DsigmaX[2];
		swapSigmaY		= Input[i].fit2DsigmaY[2];
		
		Input[i].fit2DmeanX[2]	 =   Input[i].fit2DmeanX[3];	
		Input[i].fit2DmeanY[2]	 =   Input[i].fit2DmeanY[3];	
		Input[i].fit2DsigmaX[2]	 =   Input[i].fit2DsigmaX[3];
		Input[i].fit2DsigmaY[2]	 =   Input[i].fit2DsigmaY[3];
		
		Input[i].fit2DmeanX[3]	 =   swapMeanX ;
		Input[i].fit2DmeanY[3]	 =   swapMeanY	;
		Input[i].fit2DsigmaX[3]	 =   swapsigmaX ;
		Input[i].fit2DsigmaY[3]	 =   swapSigmaY ;
	      }
	      
	      
	      //now that they are sorted, find a way to avoid cuts to intersect
	      //one possible way: - TODO this is very dirty and can be done in a much more general way
	      //first, forget about the ellipes and use circles...
	      for ( int j = 0 ; j < 4 ; j++)
	      {
		float sx,sy;
		sx = Input[i].fit2DsigmaX[j];
		sy = Input[i].fit2DsigmaY[j];
		Input[i].fit2DsigmaX[j] = (sx+sy)/2.0;
		Input[i].fit2DsigmaY[j] = (sx+sy)/2.0;
	      }
	      //now for each peak, check if any other peak is closer than the sum of the relative circles
	      for ( int j = 0 ; j < 4 ; j++)// run on all peaks
	      {
		for ( int jOther = 0 ; jOther < 4 ; jOther++)//run on all peaks again, but...
	        {
		  if (jOther != j) //...do it only if it's not the same peak
		  {
		    float distance = TMath::Sqrt( TMath::Power(Input[i].fit2DmeanX[j]-Input[i].fit2DmeanX[jOther],2) + TMath::Power(Input[i].fit2DmeanY[j]-Input[i].fit2DmeanY[jOther],2) );
		    float sumOfRadii = (Input[i].fit2DsigmaX[j]+Input[i].fit2DsigmaX[jOther]);
		    if ( distance < sumOfRadii )
		    {
		      //std::cout << "WARNING: Peaks of Module " << k << " Channel " << i << " are overlapping!" << std::endl;
		      Input[i].fit2DsigmaX[j] = distance * ( Input[i].fit2DsigmaX[j] / sumOfRadii );
		      Input[i].fit2DsigmaX[jOther] = distance * ( Input[i].fit2DsigmaX[jOther] / sumOfRadii );
		      Input[i].fit2DsigmaY[j] = Input[i].fit2DsigmaX[j];
		      Input[i].fit2DsigmaY[jOther] = Input[i].fit2DsigmaX[jOther];
		    }
		  }
		}
	      }
	      //for the moment, trivial
	      //std::cout << "Fit values" << std::endl;
	      for ( int j = 0 ; j < 4 ; j++)
	      {
		Input[i].ellipseCenterX[j] = Input[i].fit2DmeanX[j];
		Input[i].ellipseCenterY[j] = Input[i].fit2DmeanY[j];
		Input[i].ellipseWidthX[j] = Input[i].fit2DsigmaX[j];
		Input[i].ellipseWidthY[j] = Input[i].fit2DsigmaY[j];
		TEllipse *ellipses;
		ellipses = new TEllipse(Input[i].ellipseCenterX[j],Input[i].ellipseCenterY[j],Input[i].ellipseWidthX[j],Input[i].ellipseWidthY[j]);
	        //ellipses->SetFillColor(42);
	        ellipses->SetFillStyle(4001);
	        ellipses->SetLineColor(kRed);
	        ellipses->SetLineWidth(2);
	        ellipses->Draw();
		//std::cout  <<  Input[i].fit2DmeanX[j] << " " << Input[i].fit2DmeanY[j ]<< " " << Input[i].fit2DsigmaX[j] << " " << Input[i].fit2DsigmaY[j] << std::endl;
	      }
	      
	      
	      //Plot the 4 crystal spectra - TODO axis names and plot name! 
	      //create the 2x2 canvas for this channel
	      SpectraCanvasStream[i] << "SpectraCanvas" << i << "_Module_" << k;
	      SpectraCanvasString[i] = SpectraCanvasStream[i].str();
	      Input[i].SpectraCanvas = new TCanvas(SpectraCanvasString[i].c_str(),"",800,600);
	      Input[i].SpectraCanvas->Divide(2,2);
	      //then create another one just the same, but for spectra that will use the sum of 16ch 
	      std::stringstream SumSpectraCanvasStream;
	      std::string SumSpectraCanvasString;
	      SumSpectraCanvasStream << "SumSpectraCanvas" << i << "_Module_" << k;
	      SumSpectraCanvasString = SumSpectraCanvasStream.str();
	      Input[i].SumSpectraCanvas = new TCanvas(SumSpectraCanvasString.c_str(),"",800,600);
	      Input[i].SumSpectraCanvas->Divide(2,2);
	      
	      std::stringstream drawSpectrumStream[4],drawSpectrumCutStream[4],crystalStream[4];
	      std::string drawSpectrumString[4],drawSpectrumCutString[4],crystalString[4];
	      //crystalCounter = 0;
	      //draw single crystal spectra
	      for ( int j = 0 ; j < 4 ; j++) //should always run to 4
	      {
		
		//create the crystal histos with Draw
		crystalStream[j] << "Crystal " << BigCanvasMapping[crystalCounter[k]]  << " Module " << k << " - Spectrum";
		crystalString[j] = crystalStream[j].str();
		Input[i].crystal[j] = new TH1F(crystalString[j].c_str(),crystalString[j].c_str(),Params.binning,1,Params.histomax+1); 				//TH1F for this channel created here
		
		//first the variable for the "simple" crystal spectrum (i.e. using just the info from the trigger mppc to build the spectrum)
		drawSpectrumStream[j] << "ch" << i << ">>" << crystalString[j];
		drawSpectrumString[j] = drawSpectrumStream[j].str();
		
		//then the variable for the sum spectrum
		//create the crystal histos with Draw
		std::stringstream SumHistoNameStream,SumSpectraVarStream,SumSpectraCutStream;
		std::string SumHistoNameString,SumSpectraVarString,SumSpectraCutString,SumSpectraGeneralString;
		SumHistoNameStream << "Crystal " << BigCanvasMapping[crystalCounter[k]]  << " Module " << k << " - Sum Spectrum";
		SumHistoNameString = SumHistoNameStream.str();
		Input[i].crystalSum[j] = new TH1F(SumHistoNameString.c_str(),SumHistoNameString.c_str(),2*Params.binning,1,1.8*Params.histomax+1); //in general the x-axis needs to be around 2-3 time longer for sum plots
		Input[i].crystalSum[j]->GetXaxis()->SetTitle("Charge [ADC Channels]");
		Input[i].crystalSum[j]->GetYaxis()->SetTitle("Counts");
		
		
		int SumIsNotFirst = 0;
		SumSpectraVarStream << "("; //
		for(int SumIndex = 0 ; SumIndex < 32 ; SumIndex++)
		{
		  if(Input[SumIndex].Module[k])//if the channel belongs to this module
		  {
		    if(Input[SumIndex].Data)//if channel has data
		    {
		      if(SumIsNotFirst) SumSpectraVarStream << "+";
		      SumIsNotFirst++;
		      SumSpectraVarStream << "ch" << SumIndex;
		    }
		  }
		}
		SumSpectraVarStream <<  ")";
		SumSpectraGeneralString = SumSpectraVarStream.str(); //useful later
		SumSpectraVarStream <<  " >> " << SumHistoNameString;
		SumSpectraVarString = SumSpectraVarStream.str();
		
		//then we define the string that identifies the crystal in the flood histogram
		std::stringstream crystalFloodCut;
		crystalFloodCut << 
		"TriggerChannel_" <<
		k <<
		" == " << 
		i <<
		" && (( TMath::Power(FloodX_" <<
		k << 
		" - " <<
		Input[i].ellipseCenterX[j] <<
		",2) / TMath::Power(" <<
		Input[i].ellipseWidthX[j]  <<
		",2))  + ( TMath::Power(FloodY_" << 
		k << 
		" - " <<
		Input[i].ellipseCenterY[j] <<
		",2) / TMath::Power(" <<
		Input[i].ellipseWidthY[j]  <<
		",2))) < 1";
		
		//and we assign it to the cut strings for both plots
		//drawSpectrumCutStream[j] = crystalFloodCut; //reduntant
		drawSpectrumCutString[j] = crystalFloodCut.str();
		//SumSpectraCutStream = crystalFloodCut; //reduntant
		SumSpectraCutString = crystalFloodCut.str();
		
		//sum spectrum
		Input[i].SumSpectraCanvas->cd(j+1);
		if(!( (Input[i].fit2DmeanX[j] == 0) && (Input[i].fit2DmeanY[j] == 0) && (Input[i].fit2DsigmaX[j] == 0) && (Input[i].fit2DsigmaY[j] == 0) )) 
		{
		  //if this is true, a peak was found, then draw it 
		  // but also save the info that the peak was found
		  peak2Dfound[k][crystalCounter[k]] = true;
		  //std::cout << SumSpectraVarString << std::endl;
		  //std::cout << SumSpectraCutString << std::endl;
		  t1->Draw(SumSpectraVarString.c_str(),SumSpectraCutString.c_str());  
		}
		
		//plot the simple spectrum
		Input[i].SpectraCanvas->cd(j+1);
		if(!( (Input[i].fit2DmeanX[j] == 0) && (Input[i].fit2DmeanY[j] == 0) && (Input[i].fit2DsigmaX[j] == 0) && (Input[i].fit2DsigmaY[j] == 0) )) 
		{
		  t1->Draw(drawSpectrumString[j].c_str(),drawSpectrumCutString[j].c_str());
		}
		//clone the crystal spectrum (to be used later in the mixed plots)
		Input[i].crystalClone[j] = (TH1F*)Input[i].crystal[j]->Clone("clone"); //clones the histo
		
		
		
		//fitting on the simple spectrum
		//t1->Draw(drawSpectrumString[j].c_str(),drawSpectrumCutString[j].c_str(),"same");
		if(!( (Input[i].fit2DmeanX[j] == 0) && (Input[i].fit2DmeanY[j] == 0) && (Input[i].fit2DsigmaX[j] == 0) && (Input[i].fit2DsigmaY[j] == 0) )) 
		{
		  
		  //find peaks in each crystal spectrum, with TSpectrum
		  TSpectrum *s;
		  s = new TSpectrum(5);
		  //Input[i].crystal[j]->Smooth(10);
		  //Input[i].crystal[j]->Rebin(2);
		  Input[i].SpectraCanvas->cd(j+1);
		  Int_t CrystalPeaksN = s->Search(Input[i].crystal[j],1,"goff",0.5); //TODO pass to "goff"
		  Float_t *CrystalPeaks = s->GetPositionX();
		  Float_t *CrystalPeaksY = s->GetPositionY();
		  float maxPeak = 0.0;
		  int peakID = 0;
		  for (int peakCounter = 0 ; peakCounter < CrystalPeaksN ; peakCounter++ )
		  {
		    if(CrystalPeaks[peakCounter] > maxPeak)
		    {
		      maxPeak = CrystalPeaks[peakCounter];
		      peakID = peakCounter;
		    }
		  }
		  //std::cout << CrystalPeaks[0] << std::endl;
		  //std::cout << CrystalPeaksY[0] << std::endl;
		  //fit the spectra - TODO use the gaussian plus fermi?
		  float energyResolution;
		  if (correctingForSaturation)
		    energyResolution = ENERGY_RESOLUTION_SATURATION_CORRECTION;
		  else
		    energyResolution = ENERGY_RESOLUTION;
		  float par0 = CrystalPeaksY[peakID];
		  float par1 = CrystalPeaks[peakID];
		  float par2 = (CrystalPeaks[peakID]*energyResolution)/2.35;
		  //std::cout << "Crystal " << BigCanvasMapping[crystalCounter[k]] << " = " << par0 << " " << par1 << " " << par2 << std::endl;
		  gauss->SetParameter(0,par0);
		  gauss->SetParameter(1,par1);
		  gauss->SetParameter(2,par2); //expected FWHM en res = 12%
		  
		  Input[i].crystal[j]->Fit("gauss","Q","",par1-2.0*par2,par1+2.0*par2);
		  //store the mean and sigma in the data struct
		  Input[i].CrystalMean[j] = gauss->GetParameter(1);
		  Input[i].CrystalSigma[j] = std::abs(gauss->GetParameter(2));
		  //std::cout << Input[i].CrystalMean[j] << std::endl;
		  //std::cout << Input[i].CrystalSigma[j] << std::endl;
		}
		
		
		//fit on sum spectra
		if(!( (Input[i].fit2DmeanX[j] == 0) && (Input[i].fit2DmeanY[j] == 0) && (Input[i].fit2DsigmaX[j] == 0) && (Input[i].fit2DsigmaY[j] == 0) )) 
		{
		  
		  //find peaks in each crystal spectrum, with TSpectrum
		  TSpectrum *s;
		  s = new TSpectrum(5);
		  //Input[i].crystal[j]->Smooth(10);
		  //Input[i].crystal[j]->Rebin(2);
		  Input[i].SumSpectraCanvas->cd(j+1);
		  Int_t CrystalPeaksN = s->Search(Input[i].crystalSum[j],1,"goff",0.5); 
		  Float_t *CrystalPeaks = s->GetPositionX();
		  Float_t *CrystalPeaksY = s->GetPositionY();
		  float maxPeak = 0.0;
		  int peakID = 0;
		  for (int peakCounter = 0 ; peakCounter < CrystalPeaksN ; peakCounter++ )
		  {
		    if(CrystalPeaks[peakCounter] > maxPeak)
		    {
		      maxPeak = CrystalPeaks[peakCounter];
		      peakID = peakCounter;
		    }
		  }
		  //std::cout << CrystalPeaks[0] << std::endl;
		  //std::cout << CrystalPeaksY[0] << std::endl;
		  //fit the spectra - TODO use the gaussian plus fermi?
		  float energyResolution;
		  if (correctingForSaturation)
		    energyResolution = ENERGY_RESOLUTION_SATURATION_CORRECTION;
		  else
		    energyResolution = ENERGY_RESOLUTION;
		  float par0 = CrystalPeaksY[peakID];
		  float par1 = CrystalPeaks[peakID];
		  float par2 = (CrystalPeaks[peakID]*energyResolution)/2.35;
		  gauss->SetParameter(0,par0);
		  gauss->SetParameter(1,par1);
		  gauss->SetParameter(2,par2); //expected FWHM en res = 12%
		  Input[i].crystalSum[j]->Fit("gauss","Q","",par1-1.0*par2,par1+2.0*par2);
		  //store the mean and sigma in the data struct
		  Input[i].SumCrystalMean[j] = gauss->GetParameter(1);
		  Input[i].SumCrystalSigma[j] = std::abs(gauss->GetParameter(2));
		  //std::cout << Input[i].CrystalMean[j] << std::endl;
		  //std::cout << Input[i].CrystalSigma[j] << std::endl;
		}
		
		
		
		
		
		
		//and finally, the flood histogram with just photoelectric events
		std::stringstream InvarStream,IncutStream;
		std::string InvarString,IncutString;
		//IMPORTANT - the "+" here is fundamental!
		InvarStream << "FloodY_" << k << ":FloodX_" << k << " >> +" << "recFinal-Module" << k;
		InvarString = InvarStream.str();
		//the cut is the same used for each crystal spectrum, adding the cut on the spectrum itself
		//here the choice is one sigma on the left, 2 sigmas on the right
		//TODO does it make sense to allow external control on these cuts?
		IncutStream << drawSpectrumCutString[j] << " && ch" << i << " > " << Input[i].CrystalMean[j] - Input[i].CrystalSigma[j] << " && ch" << i << " < " <<  Input[i].CrystalMean[j] + 3.0*Input[i].CrystalSigma[j] ;
		IncutString = IncutStream.str();
		
		
		//build the same string but for the sum plots
		std::stringstream SumInvarStream,SumIncutStream;
		std::string SumInvarString,SumIncutString;
		//IMPORTANT - the "+" here is fundamental!
		SumInvarStream << "FloodY_" << k << ":FloodX_" << k << " >> +" << "recFinal-Module" << k;
		SumInvarString = SumInvarStream.str();
		//the cut is the same used for each crystal spectrum, adding the cut on the spectrum itself
		//here the choice is one sigma on the left, 2 sigmas on the right
		//TODO does it make sense to allow external control on these cuts?
		SumIncutStream << SumSpectraCutString << " && " << SumSpectraGeneralString << " > " << Input[i].SumCrystalMean[j] - 1.5*Input[i].SumCrystalSigma[j] << " &&" << SumSpectraGeneralString << " < " <<  Input[i].SumCrystalMean[j] + 3.0*Input[i].SumCrystalSigma[j] ;
		SumIncutString = SumIncutStream.str();
		
		
		//write this cut in the text file 
		if(peak2Dfound[k][crystalCounter[k]] == true)
		{
		  outputCutFile[k] << SumIncutString << std::endl;
		}
		
		//now a different cut, from the peak up
		std::stringstream PeakCutStream;
		std::string PeakCutString;
		
		PeakCutStream << SumSpectraCutString << " && " << SumSpectraGeneralString << " > " << Input[i].SumCrystalMean[j] << " &&" << SumSpectraGeneralString << " < " <<  Input[i].SumCrystalMean[j] + 3.0*Input[i].SumCrystalSigma[j] ;
		PeakCutString = PeakCutStream.str();
		
		
		//Update the mega stringstream
		// 	if(crystalCounter != 63)
		// 	  MegaStream << "(" << cutString << " ) | ";
		// 	else
		// 	  MegaStream << "(" << cutString << " )";
		//std::cout << cutString << std::endl;
		RecCanvasFull[k]->cd(); //select the proper canvas
		t1->Draw(InvarString.c_str(),IncutString.c_str(),"COLZ"); //plot with Draw
		
		
		
		//at this point let's try to put color on the photopeak selection
		//we have to draw another plot but with the last 
		std::stringstream cloneStream,cloneVarStream;
		std::string cloneString,cloneVarString;
		cloneStream << "Clone-" << crystalString[j];
		cloneString = cloneStream.str();
		Input[i].crystalHighlight[j] = new TH1F(cloneString.c_str(),"",Params.binning,1,Params.histomax+1);
		Input[i].crystalHighlight[j]->SetFillColor(kGreen);
		Input[i].crystalHighlight[j]->SetFillStyle(3001);
		cloneVarStream << "ch" << i << ">>" << cloneString;
		cloneVarString = cloneVarStream.str();
		Input[i].SpectraCanvas->cd(j+1);
		if(!( (Input[i].fit2DmeanX[j] == 0) && (Input[i].fit2DmeanY[j] == 0) && (Input[i].fit2DsigmaX[j] == 0) && (Input[i].fit2DsigmaY[j] == 0) )) 
		{
		  t1->Draw(cloneVarString.c_str(),IncutString.c_str(),"same");
		}
		
		//same thing but for the sum plots
		//we have to draw another plot but with the last 
		std::stringstream SumcloneStream,SumcloneVarStream;
		std::string SumcloneString,SumcloneVarString;
		SumcloneStream << "CloneSum-" << SumHistoNameString;
		SumcloneString = SumcloneStream.str();
		Input[i].SumcrystalHighlight[j] = new TH1F(SumcloneString.c_str(),"",2*Params.binning,1,1.8*Params.histomax+1);
		Input[i].SumcrystalHighlight[j]->SetFillColor(kGreen);
		Input[i].SumcrystalHighlight[j]->SetFillStyle(3001);
		SumcloneVarStream << SumSpectraGeneralString << ">>" << SumcloneString;
		SumcloneVarString = SumcloneVarStream.str();
		Input[i].SumSpectraCanvas->cd(j+1);
		if(!( (Input[i].fit2DmeanX[j] == 0) && (Input[i].fit2DmeanY[j] == 0) && (Input[i].fit2DsigmaX[j] == 0) && (Input[i].fit2DsigmaY[j] == 0) )) 
		{
		  t1->Draw(SumcloneVarString.c_str(),SumIncutString.c_str(),"same");
		}
		
		
		//do this only if you want the "deep" analysis
		if(Params.deepAnalysis)
		{
		  
		  //now, real fun. for each of the 4 crystals, take events that are in the photopeak and fill 16 spectra for the light sharing...
		  std::stringstream CrystalLightSharingCanvasStream;
		  std::string CrystalLightSharingCanvasString;
		  CrystalLightSharingCanvasStream << "Crystal " << BigCanvasMapping[crystalCounter[k]] << " Module "<< k <<  " - Light sharing";
		  CrystalLightSharingCanvasString = CrystalLightSharingCanvasStream.str();
		  Input[i].CrystalLightSharingCanvas[j] = new TCanvas(CrystalLightSharingCanvasString.c_str(),"",1200,800);
		  Input[i].CrystalLightSharingCanvas[j]->Divide(4,4);
		  
		  for (int kDeep = 0; kDeep < 32 ; kDeep++)
		  {
		    if(Input[kDeep].Module[k])//if the channel belongs to this module
		    {
		      if(Input[kDeep].Data)
		      {
			std::stringstream lightTagStream,lightTagVarStream,lightTagCutStream;
			std::string lightTagString,lightTagVarString,lightTagCutString;
			lightTagStream << "Channel " << kDeep << " - " << Input[kDeep].LabelName << " Module " << k << " - Trigger on Channel " << i << " crystal " << BigCanvasMapping[crystalCounter[k]];
			
			lightTagString = lightTagStream.str();
			
			lightTagVarStream << names[kDeep] << " >> " << lightTagString;
			lightTagVarString = lightTagVarStream.str();
			
			//lightTagCutStream << "TriggerChannel == " << i;
			//lightTagCutString = lightTagCutStream.str();
			
			Input[i].CrystalLightSharing[j][kDeep] = new TH1F(lightTagString.c_str(),lightTagString.c_str(),Params.binning,1,Params.histomax+1);
			Input[i].CrystalLightSharingCanvas[j]->cd(Input[kDeep].CanvasCd);
			//std::cout << lightTagVarString << std::endl;
			//std::cout << cutString << std::endl;
			//std::cout << std::endl;
			t1->Draw(lightTagVarString.c_str(),IncutString.c_str());
			
		      }
		    }
		  }
		  
		}
		
		//do this only if you want the "doi" analysis
		if(Params.doiAnalysis)
		{
		  std::stringstream CrystalDOICanvasStream;
		  std::string CrystalDOICanvasString;
		  CrystalDOICanvasStream << "Crystal " << BigCanvasMapping[crystalCounter[k]] << " Module "<< k <<  " - DOI";
		  CrystalDOICanvasString = CrystalDOICanvasStream.str();
		  Input[i].CrystalDOICanvas[j] = new TCanvas(CrystalDOICanvasString.c_str(),"",800,600);
		  
		  //name
		  std::stringstream doiStream,doiVarStream,doiCutStream;
		  std::string doiString,doiVarString,doiCutString;
		  doiStream << "Crystal " << BigCanvasMapping[crystalCounter[k]] << " Module " << k << " - DOI plot";
		  doiString = doiStream.str();
		  Input[i].DOIRatio[j] = new TH1F(doiString.c_str(),"",500,0,1);
		  Input[i].DOIRatio[j]->GetXaxis()->SetTitle("TriggerCh/AllCh");
		  Input[i].DOIRatio[j]->GetYaxis()->SetTitle("N");
		  Input[i].DOIRatio[j]->SetTitle(doiString.c_str());
		  Input[i].DOIRatio[j]->GetXaxis()->SetRangeUser(0,1);
		  //variables
		  //this plot is the simple trigger channel / all channels
		  //we need to define a string for this variable, so 
		  //first we have to find out which channels are active
		  int isNotFirst = 0;
		  doiVarStream << "(ch" << i << "/("; //
		  for(int doiIndex = 0 ; doiIndex < 32 ; doiIndex++)
		  {
		    if(Input[doiIndex].Module[k])//if the channel belongs to this module
		    {
		      if(Input[doiIndex].Data)//if channel has data
		      {
			if(isNotFirst) doiVarStream << "+";
			isNotFirst++;
			doiVarStream << "ch" << doiIndex;
		      }
		    }
		  }
		  doiVarStream <<  ")) >> " << doiString;
		  doiVarString = doiVarStream.str();
		  //std::cout << doiVarString << std::endl;  //debug
		  //cut
		  //the cut is the same of the flood histogram with just photoelectric events (i.e. only good events are analyzed)
		  //we should do doiCutString = IncutString; but it's really redundant...
		  
		  Input[i].CrystalDOICanvas[j]->cd();
		  t1->Draw(doiVarString.c_str(),IncutString.c_str());
		  
		  
		  //repeat for the sum plots
		  std::stringstream SumCrystalDOICanvasStream;
		  std::string SumCrystalDOICanvasString;
		  SumCrystalDOICanvasStream << "Sum Crystal " << BigCanvasMapping[crystalCounter[k]] << " Module "<< k <<  " - DOI";
		  SumCrystalDOICanvasString = SumCrystalDOICanvasStream.str();
		  Input[i].SumCrystalDOICanvas[j] = new TCanvas(SumCrystalDOICanvasString.c_str(),"",800,600);
		  
		  //name
		  std::stringstream SumdoiStream,SumdoiVarStream,SumdoiCutStream;
		  std::string SumdoiString,SumdoiVarString,SumdoiCutString;
		  SumdoiStream << "Sum Crystal " << BigCanvasMapping[crystalCounter[k]] << " Module " << k << " - DOI plot";
		  SumdoiString = SumdoiStream.str();
		  Input[i].SumDOIRatio[j] = new TH1F(SumdoiString.c_str(),"",500,0,1); //this one selects the "entire photopeak"
		  
		  Input[i].SumDOIRatio[j]->GetXaxis()->SetTitle("TriggerCh/AllCh");
		  Input[i].SumDOIRatio[j]->GetYaxis()->SetTitle("Counts");
		  Input[i].SumDOIRatio[j]->SetTitle(doiString.c_str());
		  Input[i].SumDOIRatio[j]->GetXaxis()->SetRangeUser(0,1);
		  
		  SumdoiStream << " - Peak";
		  Input[i].SumDOIRatioPeak[j] = new TH1F(SumdoiStream.str().c_str(),"",500,0,1); //this selects to the peak up
		  //variables
		  //this plot is the simple trigger channel / all channels
		  //we need to define a string for this variable, so 
		  //first we have to find out which channels are active
		  isNotFirst = 0;
		  SumdoiVarStream << "(ch" << i << "/("; //
		  for(int doiIndex = 0 ; doiIndex < 32 ; doiIndex++)
		  {
		    if(Input[doiIndex].Module[k])//if the channel belongs to this module
		    {
		      if(Input[doiIndex].Data)//if channel has data
		      {
			if(isNotFirst) SumdoiVarStream << "+";
			isNotFirst++;
			SumdoiVarStream << "ch" << doiIndex;
		      }
		    }
		  }
		  SumdoiVarStream <<  ")) >> " ;
		  std::string temporaryString,SumdoiVarStringPeak;
		  temporaryString = SumdoiVarStream.str();
		  
		  SumdoiVarString = temporaryString + SumdoiString;
		  
		  SumdoiVarStringPeak = temporaryString + SumdoiStream.str();
		  //std::cout << doiVarString << std::endl;  //debug
		  //cut
		  //the cut is the same of the flood histogram with just photoelectric events (i.e. only good events are analyzed)
		  //we should do doiCutString = IncutString; but it's really redundant...
		  
		  Input[i].SumCrystalDOICanvas[j]->cd();
		  Input[i].SumDOIRatio[j]->SetLineColor(1);
		  Input[i].SumDOIRatio[j]->SetFillColor(1);
		  Input[i].SumDOIRatio[j]->SetFillStyle(3001);
		  Input[i].SumDOIRatioPeak[j]->SetLineColor(2);
		  Input[i].SumDOIRatioPeak[j]->SetFillColor(2);
		  Input[i].SumDOIRatioPeak[j]->SetFillStyle(3001);
		  t1->Draw(SumdoiVarString.c_str(),SumIncutString.c_str());
		  t1->Draw(SumdoiVarStringPeak.c_str(),PeakCutString.c_str(),"same");
		  
		  
		}
		
		
		//and finally the big spectra canvas (aka the 64 crystals)
		BigSpectraCanvas[k]->cd(BigCanvasMapping[crystalCounter[k]]);
		std::stringstream CrystalTitleStream;
		std::string CrystalTitleString;
		CrystalTitleStream << "Crystal " << BigCanvasMapping[crystalCounter[k]];
		CrystalTitleString = CrystalTitleStream.str();
		Input[i].crystal[j]->SetTitle(CrystalTitleString.c_str());
		Input[i].crystal[j]->Draw();
		Input[i].crystalHighlight[j]->Draw("same");
		
		
		SumBigSpectraCanvas[k]->cd(BigCanvasMapping[crystalCounter[k]]);
// 		std::stringstream CrystalTitleStream;
// 		std::string CrystalTitleString;
// 		CrystalTitleStream << "Crystal " << BigCanvasMapping[crystalCounter[k]];
// 		CrystalTitleString = CrystalTitleStream.str();
		Input[i].crystalSum[j]->SetTitle(CrystalTitleString.c_str());
		Input[i].crystalSum[j]->Draw();
		Input[i].SumcrystalHighlight[j]->Draw("same");
		
		//update the crystal counter even if the peak was not found by TSpectrum2
		crystalCounter[k]++;
	      }
	      
	      
	      
	      
	      
	      //the "mixed plots"
	      MixedStream << "Channel " << i << " Module "<< k << " - Mixed Plot";
	      MixedString = MixedStream.str();
	      Input[i].MPPCCanvas = new TCanvas(MixedString.c_str(),"",800,600);
	      Input[i].MPPCCanvas->cd();
	      TLegend *legend = new TLegend(0.13,0.7,0.39,0.89,"");
	      legend->SetFillStyle(0);
	      Input[i].channelClone->SetLineColor(1);
	      Input[i].channelClone->SetFillColor(1);
	      Input[i].channelClone->SetFillStyle(3001);
	      legend->AddEntry(Input[i].channelClone,"MPPC spectrum","f");
	      Input[i].channelClone->Draw(); 
	      for (int j = 0 ; j < 4 ; j++)
	      {
		std::stringstream legendStream;
		std::string legendName;
		legendStream << "Crystal " <<  BigCanvasMapping[crystalCounter2[k]] << " Module " << k;
		legendName = legendStream.str();
		Input[i].crystalClone[j]->SetLineColor(j+2);
		Input[i].crystalClone[j]->SetFillColor(j+2);
		Input[i].crystalClone[j]->SetFillStyle(3001);
		legend->AddEntry(Input[i].crystalClone[j],legendName.c_str(),"f");
		Input[i].crystalClone[j]->Draw("same");
		crystalCounter2[k]++;
	      }
	      legend->Draw();
	      
	    }
	    
	  }
	}
	outputCutFile[k].close();
      }
    }
    
    //     //Finally, put only the interesting events in another ttree
    //     TFile* fTreePET = new TFile("prova.root","recreate");
    //     TTree *t2 = new TTree("events","events"); 
    //     
    //     //the mega string of all the cuts (from the final draw)
    //     MegaString = MegaStream.str();
    //     std::cout << MegaString << std::endl;
    //     
    //     //increase the limits in TFormula, because the cut string is way too long for the defaults 
    //     TFormula::SetMaxima(100000,1000,1000000);
    //     //create the formula from the cuts
    //     TTreeFormula* Formula = new TTreeFormula("Formula",MegaString.c_str(),t1);
    //     
    //     long long int TimeTag;
    //     int CrystalID;
    //     t2->Branch("TimeTag",&TimeTag,"TimeTag/L"); 			//absolute time tag of the event
    //     t2->Branch("CrystalID",&CrystalID,"CrystalID/I"); 			//crystal where the interaction occurred
    //     
    //     //loop on entries
    //     Long64_t nentries = t1->GetEntriesFast();   //nentries num eventi
    //     for (Long64_t jentry=0; jentry<nentries;jentry++)
    //     {
    //       t1->GetEntry(jentry);
    //       if( Formula->EvalInstance() )
    //       {
    // 	TimeTag = ExtendedTimeTag;
    // 	CrystalID = 0; //TODO
    // 	t2->Fill();
    //       }
    //     }
    //     
    //     fTreePET->cd();
    //     t2->Write();
    //     fTreePET->Close();
    /*
     * for(int i = 0 ; i < 32 ; i++)
     * {
     *   char channelName[50]; //a string will be needed for the TH1F creation inside the following loop. Root likes strings...
     *   char channelNameHighlight[80];
     *   char FillTheHisto[70]; //a string will be needed for filling TH1F from the tree data, inside the following loop
     *   char FillTheHistoTag[70]; //a string will be needed for the cut conditions
     *   char FillTheHistoTagHighlight[100];
     *   char FillTheHistoTagHighlightCut[100];
     *   
     *   
     *   //if (Input[i].LabelName == TriggerChannel) sprintf(channelName,"Channel %d - %s - TriggerChannel",i,Input[i].LabelName.Data()); //name of the histogram
     *   sprintf(channelName,"Channel %d - %s",i,Input[i].LabelName.Data()); //name of the histogram
     *   sprintf(channelNameHighlight,"Channel Highlight %d - %s",i,Input[i].LabelName.Data());
     *   
     *   sprintf(FillTheHistoTagHighlight,"%s>>%s",names[i].c_str(),channelNameHighlight); //string to pass to Draw(), to fill the histogram
     *   sprintf(FillTheHisto,"%s>>%s",names[i].c_str(),channelName); //string to pass to Draw(), to fill the histogram
     *   
     *   sprintf(FillTheHistoTag,"TriggerChannel == %d",i); //string to pass to Draw(), to set the "tag"
     *   
     *   channel[i] = new TH1F(channelName,"",1024,0,4096); //TH1F for this channel created here
     *   channel[i]->SetTitle(channelName);
     *   //channel[i]->SetTitleSize(0.9);
     *   channel[i]->GetXaxis()->SetTitle("ADC Channels");
     *   //     channel[i]->GetXaxis()->SetLabelSize(0.06);
     *   //     channel[i]->GetYaxis()->SetLabelSize(0.06);
     *   //     channel[i]->GetXaxis()->SetLabelOffset(-0.002);
     *   //     channel[i]->GetXaxis()->SetNdivisions(506);
     *   //     channel[i]->GetYaxis()->SetNdivisions(505);
     *   
     *   channel[i]->GetYaxis()->SetTitle("Counts");
     *   //     channel[i]->GetYaxis()->SetTitleSize(0.07);
     *   //     channel[i]->GetYaxis()->SetTitleOffset(0.5);
     *   //     channel[i]->GetXaxis()->SetTitleSize(0.06);
     *   //     channel[i]->GetXaxis()->SetTitleOffset(0.7);
     *   if (Input[i].Data) 
     *     Canvas1->cd(Input[i].CanvasCd);//select which subcanvas is used for drawing this histo  
     *     else 
     * Canvas2->cd(Input[i].CanvasCd);//select which subcanvas is used for drawing this histo  
     * 
     * t1->Draw(FillTheHisto,FillTheHistoTag);//Draw the histo taking data from the tree. For example, the real command would look like 
     * //     if (i!=ch) channel[i]->SetStats(0);
     * //     else channel[i]->SetStats(1);
     * 
     * //if the channel that is being filled in this step of the loop is the trigger 
     * //channel, fit its histogram and get the mean and sigma of the photopeak 
     * // 	if (Input[i].LabelName == TriggerChannel)
     * // 	{
     * // 	  TriggerChannelID = i;
     * channel[i]->SetStats(1);
     * // 	  //gStyle->SetStatX(1);
     * // 	  //gStyle->SetStatY(1);
     * // 	  //gStyle->SetStatW(1);
     * // 	  //gStyle->SetStatH(1);
     * // 	  TF1 *gauss = new TF1("gauss",  "[0]*exp(-0.5*((x-[1])/[2])**2)",fitMin,fitMax);
     * if(!cutFile.is_open())
     * {
     *  channel[i]->Fit("gaus","","",2400,3200);
     *  mean[i] = gaus->GetParameter(1);
     *  sigma[i] = gaus->GetParameter(2);
  }
  
  if(cutFile.is_open())
  {
  setCutLow[i] = Input[i].lowCut;
  setCutHigh[i] =Input[i].highCut;
  }
  else
  {
  setCutLow[i] = (mean[i]-HowManySigma*sigma[i]);
  setCutHigh[i] =(mean[i]+HowManySigma*sigma[i]);
  }
  std::cout << setCutLow[i] << " " << setCutHigh[i] << std::endl;
  
  sprintf(FillTheHistoTagHighlightCut,"TriggerChannel == %d && %s > %d && %s < %d",i,names[i].c_str(),setCutLow[i],names[i].c_str(),setCutHigh[i]); //string to pass to Draw(), to set the "tag"
  
  hsigMppc[i] = new TH1F(channelNameHighlight,"",1024,0,4096); //TH1F for this channel created here
  hsigMppc[i]->SetFillColor(3);//set fill color
  t1->Draw(FillTheHistoTagHighlight,FillTheHistoTagHighlightCut,"same");
  //hsigMppc[i]->GetXaxis()->SetRange(setCutLow[i],setCutHigh[i]);//set the range according to the selection
  
  //hsigMppc[i]->Draw("same");//draws on the same histo
  //gStyle->SetOptFit(0001);
  // 	  
  // 	  //highlighting the events that we use for the selection
  //  TH1 *hsigMppc = (TH1*)channel[i]->Clone("hsigMppc");//clones histo
  //  hsigMppc->GetXaxis()->SetRange((mean[i] - HowManySigma*sigma[i]),(mean[i] + HowManySigma*sigma[i]));//set the range according to the selection
  //  hsigMppc->SetFillColor(3);//set fill color
  //  hsigMppc->Draw("same");//draws on the same histo
  // 	}
  
  //delete channelName;
  //delete FillTheHisto;
  //delete channelCutName;
  //delete FillTheHistoCutVar;
  //delete FillTheHistoCutCuts;
  printf("%s %s %d %d\n",FillTheHisto,FillTheHistoTag,setCutLow[i],setCutHigh[i]);
  }*/
    
    //save lightsharing data, if needed
    if(Params.LightSharingAnalysis)
    {
      //make the folder
      std::string folderName = "lightSharing";
      std::string command = std::string("mkdir -p ") + folderName;
      system(command.c_str());
      
      //out file of the sums LY
      std::stringstream outFileLYSStream;
      std::string outFileLYName;
      outFileLYSStream << "./" << folderName << "/SumLY.dat";
      outFileLYName = outFileLYSStream.str();
      std::ofstream ofsLY(outFileLYName.c_str(), std::ofstream::out);
      
      
      for (int k = 0 ; k < 2 ; k++)
      {
	if(ModuleOn[k])
	{
	  for (int i = 0 ; i < 32 ; i++) //run on all channels
	  {
	    if(Input[i].Module[k])//if the channel belongs to this module
	    {
	      if(Input[i].Data) //if channel has data
	      {
		ofsLY << i << " " << Input[i].SumLYMean << std::endl;
		std::stringstream outFileSStream;
		std::string outFileName;
		outFileSStream << "./" << folderName << "/Ch" << i << "-" <<Input[i].LabelName << ".dat";
		outFileName = outFileSStream.str();

		std::ofstream ofs(outFileName.c_str(), std::ofstream::out);
		for (int kTag = 0; kTag < 32 ; kTag++)
		{
		  if(Input[kTag].Module[k])//if the channel belongs to this module
		  {
		    if(Input[kTag].Data)//and if the channel is on
		    {
		      ofs << kTag << "\t" << Input[i].lightSharingData[kTag] << std::endl;
		    }
		  }
		}
		ofs.close();
		
	      }
	    }
	  }
	}
      }
      ofsLY.close();
    }
    
    
    if(Params.saveAnalysisTree)
    {
      TFile* fTree = new TFile("temp.root","recreate");
      fTree->cd();
      t1->Write();
      fTree->Close();
      std::cout << "Temporary tree saved to temp.root file" << std::endl;
    }
    
    
    //now save the histograms in a root file
    TString fileRootHisto;
    //     if(treeRoot == "")
    fileRootHisto = "histo_" + fileRoot;
    TFile* fHisto = new TFile(fileRootHisto,"recreate"); //open the root file
    std::cout << "Writing results to ROOT file " << fileRootHisto << std::endl;
    TDirectory *directory[2][33][5]; //TDirectory 
    
    for (int k = 0 ; k < 2 ; k++)
    {
      if(ModuleOn[k])
      {
	int crystalCounter = 0;
	std::stringstream ModuleDirStream;
	std::string ModuleDirString;
	ModuleDirStream << "Module " << k;
	ModuleDirString = ModuleDirStream.str();
	directory[k][0][0] = fHisto->mkdir(ModuleDirString.c_str());
	directory[k][0][0]->cd();
	
	CanvasOriginal[k]->Write();
	Canvas1[k]->Write();
	PFCanvasOriginal[k]->Write();
	PFCanvas1[k]->Write();
	LYCanvasOriginal[k]->Write();
	LYCanvas1[k]->Write();
	FloodCanvasClean[k]->Write();
	//CanvasFirstOnSecond[k]->Write(); //not really needed
	//FloodCanvasCutOnRatio[k]->Write(); //not really needed
	FloodCanvas[k]->Write();
	//recFinal->Write();  
	RecCanvasFull[k]->Write();
	BigSpectraCanvas[k]->Write();
	SumBigSpectraCanvas[k]->Write();
	for (int i = 0 ; i < 32 ; i++) //run on all channels
	{
	  if(Input[i].Module[k])//if the channel belongs to this module
	  {
	    if(Input[i].Data)
	    {
	      std::stringstream DirStream;
	      std::string DirString;
	      DirStream << "Channel " << i << " - " << Input[i].LabelName;
	      DirString = DirStream.str();
	      directory[k][i+1][0] = directory[k][0][0]->mkdir(DirString.c_str());
	      directory[k][i+1][0]->cd();
	      //Input[i].recFlood->Write();
	      Input[i].channel->Write();
	      Input[i].LightSharingCanvas->Write();
	      Input[i].LightSharingCanvasPE->Write();
	      Input[i].LYLightSharingCanvasPE->Write();
	      Input[i].recChannelCanvas->Write();
	      Input[i].SpectraCanvas->Write();
	      Input[i].SumSpectraCanvas->Write();
	      Input[i].MPPCCanvas->Write();
	      
	      if(Params.deepAnalysis | Params.doiAnalysis)
	      {
		for(int j = 0 ; j < 4 ; j++)
		{
		  std::stringstream SubDirStream;
		  std::string SubDirString;
		  SubDirStream << "Crystal " <<  BigCanvasMapping[crystalCounter];
		  SubDirString = SubDirStream.str();
		  directory[k][i+1][j+1] = directory[k][i+1][0]->mkdir(SubDirString.c_str());
		  crystalCounter++;
		}
	      }
	      
	      if(Params.deepAnalysis)
	      {
		for(int j = 0 ; j < 4 ; j++)
		{
		  directory[k][i+1][j+1]->cd();
		  Input[i].crystal[j]->Write();
		  Input[i].CrystalLightSharingCanvas[j]->Write();
		}
	      }
	      if(Params.doiAnalysis)
	      {
		for(int j = 0 ; j < 4 ; j++)
		{
		  directory[k][i+1][j+1]->cd();
		  Input[i].CrystalDOICanvas[j]->Write();
		  Input[i].SumCrystalDOICanvas[j]->Write();
		}
	      }
	      
	      
	    }
	  }
	}
      }
    }
    //   for (int i = 0 ; i < 32 ; i++) //write the 64 crystal spectra
    //   {
    //     if(Input[i].Data)
    //     {
    //       
    //     }
    //   }
    //save analysis tree if the user wants it
    
    fHisto->Close();
  }
}








int parseConfigFile(std::istream& FileHeader,Params_t &Params)
{
  // defines the strings to look for
  //   std::string FileName = "name of data file := ";       
  //   std::string NumberFormat = "!number format := ";             
  //   std::string SizeX = "!matrix size [1] := ";         
  //   std::string SizeY = "!matrix size [2] := ";
  //   std::string SizeZ = "!matrix size [3] := "; 
  //   std::string SpacingX = "scaling factor (mm/pixel) [1] := ";	 
  //   std::string SpacingY = "scaling factor (mm/pixel) [2] := ";	 
  //   std::string SpacingZ = "scaling factor (mm/pixel) [3] := ";
  std::stringstream ChIdentifierStream[32];
  std::string ChIdentifierString[32];
  std::string saturation 		= "!Saturation File Name =";
  std::string analysisString 		= "!Deep Analysis =";
  std::string coincidenceString 	= "!Coincidence =";
  std::string histoMaxString 		= "!Histo 1D Max =";
  std::string histoBinString 		= "!Histo 1D Bin =";
  std::string histo2DBinString 		= "!Histo 2D-global Bin =";
  std::string histo2DChBinString 	= "!Histo 2D-channel Bin =";
  std::string PEminString 		= "!Photoelectric Peak Min =";
  std::string PEmaxString 		= "!Photoelectric Peak Max =";
  std::string GridCutString 		= "!Grid Cut Width =";
  std::string SaveTreeString 		= "!Save Analysis TTree =";
  std::string doiAnalysisString 	= "!DOI Analysis =";
  std::string LightSharingAnalysisString= "!Light Sharing Analysis =";
  std::string MppcPDEString		= "!Photodetector PDE =";
  std::string MppcGainString		= "!Photodetector Gain =";
  std::string SourgeEnergyString	= "!Gamma Source Energy =";
  std::string AdcChargeBinningString	= "!ADC charge binning =";
  std::string OnlyLyAnalysysString	= "!Only LY Analysis = ";
  
  //std::string analysisChoice;
  for (int i = 0; i < 32 ; i++)
  {
    ChIdentifierStream[i] << "!Channel [" << i << "] =";
    ChIdentifierString[i] = ChIdentifierStream[i].str(); //set the search channel strings
    Params.ChannelLabel[i] = "VOID"; //Initialize to VOID, then only the specified will be changed
    Params.ChannelOn[i] = false; //Initialize to off all channels 
    
  }
  
  Params.deepAnalysis = false;
  Params.doiAnalysis = false;
  Params.coincidence = false;
  Params.LightSharingAnalysis = false;
  Params.OnlyLyAnalysys = false;
  Params.saturationFileName = "";
  Params.histomax= 8192;
  Params.LYhistomax = (8192 * 156e-15) / (1.25e6 * 0.3 * 0.662 * 1.6e-19);	//if nothing is specified, we assume the gain 1.25e6, pde 0.3, cesium source 0.662 MeV, charge binning 156fC
  Params.PFhistomax = (8192 * 156e-15) / (1.25e6 * 1.6e-19);
  Params.MppcPDE = 0.3;
  Params.MppcGain = 1.25e6;
  Params.SourceEnergy = 0.662;
  Params.AdcChargeBinning = 156e-15;
  Params.binning = 256;
  Params.LYbinning = (int) round((256 * 156e-15) / (1.25e6 * 0.3 * 0.662 * 1.6e-19));
  Params.PFbinning = (int) round((256 * 156e-15) / (1.25e6 * 1.6e-19));
  Params.histo2DglobalBin = 1000;
  Params.histo2DchannelBin = 150;
  Params.PEmin = 2000;
  Params.PEmin = 2400;
  Params.lipCutLevel = 0.;  //default to 0, and if it remains 0 after reading the config file, it means "don't cut the grid"
  Params.CutGrid = false;
  
  //read the config file
  
  std::string s1; 
  //looks for the strings defined above, assing the rest of the line to the correct variable      
  while(getline(FileHeader,s1))	
  {
    unsigned int i;
    //look for global parameters
    i = s1.find(saturation); 
    if (i==0)  
    {
      Params.saturationFileName.assign(s1, saturation.size() , s1.size());	
      Params.saturationFileName.erase(remove_if(Params.saturationFileName.begin(), Params.saturationFileName.end(), isspace), Params.saturationFileName.end());
    }
    i = s1.find(analysisString);
    if (i==0)  
    {
      std::string analysis;
      analysis.assign(s1, analysisString.size() , s1.size());
      analysis.erase(remove_if(analysis.begin(), analysis.end(), isspace), analysis.end());
      if(analysis == "yes" |analysis == "Yes" |  analysis == "YES" | analysis == "on" | analysis == "On" | analysis == "ON" | analysis == "1" ){
	Params.deepAnalysis = true;
      } 
    }
    i = s1.find(coincidenceString);
    if (i==0)  
    {
      std::string coincidence;
      coincidence.assign(s1, coincidenceString.size() , s1.size());	
      coincidence.erase(remove_if(coincidence.begin(), coincidence.end(), isspace), coincidence.end());
      if(coincidence == "yes" |coincidence == "Yes" |  coincidence == "YES" | coincidence == "on" | coincidence == "On" | coincidence == "ON" | coincidence == "1" ){
	Params.coincidence = true;
      }
    }
    i = s1.find(SaveTreeString);
    if (i==0)  
    {
      std::string analysis;
      analysis.assign(s1, SaveTreeString.size() , s1.size());
      analysis.erase(remove_if(analysis.begin(), analysis.end(), isspace), analysis.end());
      if(analysis == "yes" |analysis == "Yes" |  analysis == "YES" | analysis == "on" | analysis == "On" | analysis == "ON" | analysis == "1" )
      {
	Params.saveAnalysisTree = true;
      } 
    }
    i = s1.find(doiAnalysisString);
    if (i==0)  
    {
      std::string analysis;
      analysis.assign(s1, doiAnalysisString.size() , s1.size());
      analysis.erase(remove_if(analysis.begin(), analysis.end(), isspace), analysis.end());
      if(analysis == "yes" |analysis == "Yes" |  analysis == "YES" | analysis == "on" | analysis == "On" | analysis == "ON" | analysis == "1" )
      {
	Params.doiAnalysis = true;
      } 
    }
    i = s1.find(LightSharingAnalysisString);
    if (i==0)  
    {
      std::string analysis;
      analysis.assign(s1, LightSharingAnalysisString.size() , s1.size());
      analysis.erase(remove_if(analysis.begin(), analysis.end(), isspace), analysis.end());
      if(analysis == "yes" |analysis == "Yes" |  analysis == "YES" | analysis == "on" | analysis == "On" | analysis == "ON" | analysis == "1" )
      {
	Params.LightSharingAnalysis = true;
      } 
    }
    
    i = s1.find(OnlyLyAnalysysString);
    if (i==0)  
    {
      std::string analysis;
      analysis.assign(s1, OnlyLyAnalysysString.size() , s1.size());
      analysis.erase(remove_if(analysis.begin(), analysis.end(), isspace), analysis.end());
      if(analysis == "yes" |analysis == "Yes" |  analysis == "YES" | analysis == "on" | analysis == "On" | analysis == "ON" | analysis == "1" )
      {
	Params.OnlyLyAnalysys = true;
      } 
    }

    i = s1.find(histoMaxString); 
    if (i==0)  
    {
      std::string maxHisto;
      maxHisto.assign(s1, histoMaxString.size() , s1.size());	
      maxHisto.erase(remove_if(maxHisto.begin(), maxHisto.end(), isspace), maxHisto.end());
      Params.histomax = atoi(maxHisto.c_str());
    }
    i = s1.find(histoBinString); 
    if (i==0)  
    {
      std::string binHisto;
      binHisto.assign(s1, histoBinString.size() , s1.size());	
      binHisto.erase(remove_if(binHisto.begin(), binHisto.end(), isspace), binHisto.end());
      Params.binning = atoi(binHisto.c_str());
    }
    
    
    i = s1.find(histo2DBinString); 
    if (i==0)  
    {
      std::string tempString;
      tempString.assign(s1, histo2DBinString.size() , s1.size());	
      tempString.erase(remove_if(tempString.begin(), tempString.end(), isspace), tempString.end());
      Params.histo2DglobalBin = atoi(tempString.c_str());
    }
    
    
    
    i = s1.find(histo2DChBinString); 
    if (i==0)  
    {
      std::string tempString;
      tempString.assign(s1, histo2DChBinString.size() , s1.size());	
      tempString.erase(remove_if(tempString.begin(), tempString.end(), isspace), tempString.end());
      Params.histo2DchannelBin = atoi(tempString.c_str());
    }
    
    
    
    
    
    
    
    i = s1.find(MppcPDEString); 
    if (i==0)  
    {
      std::string tempString;
      tempString.assign(s1, MppcPDEString.size() , s1.size());	
      tempString.erase(remove_if(tempString.begin(), tempString.end(), isspace), tempString.end());
      Params.MppcPDE = atof(tempString.c_str());
    }
    
    i = s1.find(MppcGainString); 
    if (i==0)  
    {
      std::string tempString;
      tempString.assign(s1, MppcGainString.size() , s1.size());	
      tempString.erase(remove_if(tempString.begin(), tempString.end(), isspace), tempString.end());
      Params.MppcGain = atof(tempString.c_str());
    }
    
    i = s1.find(SourgeEnergyString); 
    if (i==0)  
    {
      std::string tempString;
      tempString.assign(s1, SourgeEnergyString.size() , s1.size());	
      tempString.erase(remove_if(tempString.begin(), tempString.end(), isspace), tempString.end());
      Params.SourceEnergy = atof(tempString.c_str());
    }
    
    i = s1.find(AdcChargeBinningString); 
    if (i==0)  
    {
      std::string tempString;
      tempString.assign(s1, AdcChargeBinningString.size() , s1.size());	
      tempString.erase(remove_if(tempString.begin(), tempString.end(), isspace), tempString.end());
      Params.AdcChargeBinning = atof(tempString.c_str());
    }
    
    
    
    i = s1.find(PEminString); 
    if (i==0)  
    {
      std::string tempString;
      tempString.assign(s1, PEminString.size() , s1.size());	
      tempString.erase(remove_if(tempString.begin(), tempString.end(), isspace), tempString.end());
      Params.PEmin = atof(tempString.c_str());
    }
    i = s1.find(PEmaxString); 
    if (i==0)  
    {
      std::string tempString;
      tempString.assign(s1, PEmaxString.size() , s1.size());	
      tempString.erase(remove_if(tempString.begin(), tempString.end(), isspace), tempString.end());
      Params.PEmax = atof(tempString.c_str());
    }
    i = s1.find(GridCutString); 
    if (i==0)  
    {
      std::string tempString;
      tempString.assign(s1, GridCutString.size() , s1.size());	
      tempString.erase(remove_if(tempString.begin(), tempString.end(), isspace), tempString.end());
      Params.lipCutLevel = atof(tempString.c_str());
      if(Params.lipCutLevel != 0)
      {
	Params.CutGrid = true;
      }
    }
    
    
    
    //channel parameters
    for (int j = 0 ; j < 32 ; j++)
    {
      i = s1.find(ChIdentifierString[j]); //find the info about the channel
      if (i==0) //if the channel is in the config file  
      {
	//first, turn on the channel
	Params.ChannelOn[j] = true;
	//then, tokenize the string
	std::string parametersString,noLeadingSpace;
	parametersString.assign(s1, ChIdentifierString[j].size() , s1.size());
	noLeadingSpace = parametersString.substr(parametersString.find_first_not_of(" "),parametersString.length()-parametersString.find_first_not_of(" "));
	
	//std::cout << "noLeadingSpace " << noLeadingSpace << std::endl;
	std::string::size_type pos = noLeadingSpace.find_first_of(" ,.-\t");
	std::string token = noLeadingSpace.substr(0, pos);
	std::string::size_type pos1 = noLeadingSpace.find_first_of(" ,.-\t",pos+1);
	std::string token1 = noLeadingSpace.substr(pos+1, pos1-1);
	std::string::size_type pos2 = noLeadingSpace.find_first_of(" ,.-\t",pos1+1);
	std::string token2 = noLeadingSpace.substr(pos1+1, pos2-1);
	//std::cout << token << " " << token1 << " " << token2 << std::endl;
	token.erase(remove_if(token.begin(), token.end(), isspace), token.end());
	token1.erase(remove_if(token1.begin(), token1.end(), isspace), token1.end());
	token2.erase(remove_if(token2.begin(), token2.end(), isspace), token2.end());
	//std::cout << token << " " << token1 << " " << token2 << std::endl;
	//finally assign to the Params
	Params.ChannelModule[j] = atoi(token.c_str());
	Params.ChannelLabel[j] = token1;
	Params.ChannelPosition[j] = atoi(token2.c_str());
      }
    }
    
  }
  //now we can calculate the rescaling for the plots where light is expressed in terms of Ph/MeV (rather then in ADC counts)
  
  Params.LYhistomax 	= (Params.histomax*Params.AdcChargeBinning)/(Params.SourceEnergy*Params.MppcGain*Params.MppcPDE*1.6*TMath::Power(10,-19));
  //Params.LYbinning	= (int) round((Params.binning*Params.AdcChargeBinning)/(Params.SourceEnergy*Params.MppcGain*Params.MppcPDE*1.6*TMath::Power(10,-19))); 
  Params.LYbinning	= Params.binning; 
  
  
  Params.PFhistomax	=  (Params.histomax*Params.AdcChargeBinning)/(Params.MppcGain*1.6*TMath::Power(10,-19));
  //Params.PFbinning	= (int) round((Params.binning*Params.AdcChargeBinning)/(Params.MppcGain*1.6*TMath::Power(10,-19)));
  Params.PFbinning	= Params.binning;
  
  return 0;
}

int printConfigFile(Params_t &Params)
{
  std::cout << "/***********************************************/" << std::endl;
  std::cout << "|                                               |" << std::endl;
  std::cout << "|          Printing config parameters           |" << std::endl;
  std::cout << "|                                               |" << std::endl;
  std::cout << "/***********************************************/" << std::endl;
  std::cout << "Saturation File =\t" 			<< Params.saturationFileName 		<< std::endl;
  std::cout << "Coincidence \t=\t" 			<< Params.coincidence 			<< std::endl;
  std::cout << "Deep analysis \t=\t" 			<< Params.deepAnalysis 			<< std::endl;
  std::cout << "DOI analysis \t=\t" 			<< Params.doiAnalysis 			<< std::endl;
  std::cout << "Save analysis TTree \t=\t" 		<< Params.saveAnalysisTree 		<< std::endl;
  std::cout << "Perform Grid Cut \t=\t" 		<< Params.CutGrid 			<< std::endl;
  if(Params.CutGrid)
    std::cout << "Grid Cut Width \t=\t" 		<< Params.lipCutLevel 			<< std::endl;
  std::cout << "Histos 1D Max \t=\t" 			<< Params.histomax 			<< std::endl;
  std::cout << "Histos 1D Bins \t=\t" 			<< Params.binning 			<< std::endl;
  std::cout << "Histos 2D-global Bins \t=\t" 		<< Params.histo2DglobalBin 		<< std::endl;
  std::cout << "Histos 2D-channel Bins \t=\t" 		<< Params.histo2DchannelBin 		<< std::endl;
  std::cout << "Histo 1D LY max \t=\t" 			<< Params.LYhistomax 			<< std::endl;
  
  std::cout << "Photodetector PDE \t=\t" 		<< Params.MppcPDE 			<< std::endl;
  std::cout << "Photodetector Gain \t=\t" 		<< Params.MppcGain 			<< std::endl;
  std::cout << "Gamma Source Energy \t=\t" 		<< Params.SourceEnergy 			<< std::endl;
  std::cout << "ADC charge binning \t=\t" 		<< Params.AdcChargeBinning 		<< std::endl;
  
  std::cout << "Photoelectric hard coded cut Min (optional) \t=\t" 		<< Params.PEmin 		<< std::endl;
  std::cout << "Photoelectric hard coded cut Max (optional) \t=\t" 		<< Params.PEmax 		<< std::endl;
  for (int i = 0 ; i < 32 ; i++){
    std::cout << "Channel [" << i << "] \t=\t"	  	<< Params.ChannelOn[i] 		<< "\t"
    << Params.ChannelModule[i] 		<< "\t"
    << Params.ChannelLabel[i] 		<< "\t"
    << Params.ChannelPosition[i] 		<< std::endl; 
  }
  return 0;
}


int makeCanvas(TCanvas * Canvas,std::string titleOrName,int ModuleNum,int width,int heigth)
{
  //std::cout << Input[4].LabelName << " " << Params.ChannelLabel[4] << " " << t1->GetEntries() <<  std::endl;
  std::stringstream CanvasStream;
  std::string CanvasString;

  CanvasStream << titleOrName << ModuleNum;
  CanvasString = CanvasStream.str();
  Canvas->SetName(CanvasString.c_str());
  Canvas->SetTitle(CanvasString.c_str());
  Canvas->SetCanvasSize(width,heigth);
//   CanvasStream.str(std::string());
  
}


