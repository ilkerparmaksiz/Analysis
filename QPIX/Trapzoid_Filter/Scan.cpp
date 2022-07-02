/************************************************************************
 *
 *  Filename: Scan.cpp
 *s
 *  Description:
 *
 *	Author(s):
 *     Michael T. Febbraro
 *     David Walter
 *
 *  Creation Date: 9/25/2016
 *  Last modified: 8/24/2018
 *
 *  To compile: g++ -O3 -pedantic -o Scan.exe `root-config --cflags --libs` -lSpectrum NewScan.cpp
 *      - if errors copiling in Mac OSX
 *        - remove -03 option
 *        - clang++ instead of g++
 *        - $(root-config --cflags --libs) instead of 'root-config --cflags --libs'
 *
 *
 *  If "error while loading shared libraries: libcore.so: ..." occurs, type
 *  "source `root-config --prefix`/bin/thisroot.sh"
 *
 * -----------------------------------------------------
 * 	Nuclear Astrophysics, Physics Division
 *  Oak Ridge National Laboratory
 *
 */


#include <iostream>
#include <fstream>
#include <string>
#include <signal.h>
#include "PulseAnalysis.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TMath.h"
#include "TF1.h"

using namespace std;

typedef struct
{
  Float_t l;               // Long integral
  Float_t s;               // Short integral
  Float_t amp;             // Amplitude
  Float_t cfd;             // Trigger time
  Float_t psd;             // PSD parameter s/l
  Float_t trg;             // Det trigger
  Float_t pp;              // Position of max amplitude

} Can;

int PeakFinder(float* trace, int length, int window, int rangeL, int rangeH, float threshold, vector <float> *pa, vector <float> *pt, vector <float> *pi) {
  int i, j, k;
  bool toggle = 0;
  float peak = 0;
  float integral;
  int loc, npeaks;

  pa->clear();
  pi->clear();
  pt->clear();
  npeaks = 0;

  for (i = rangeL; i < rangeH; i++)
  {
    if (trace[i] > peak)
    {
      peak = trace[i];
      loc = i;
      j = 0;
    }
    else if (trace[i] > threshold)
    {
      i += (int)(window/2);
      pa->push_back(peak);
      pt->push_back((float)loc);
      integral = 0;
      for (k = loc - 5; k <= loc + 5; k++) { integral += trace[k];}
      pi->push_back(integral);
      peak = 0;
      npeaks++;
    }
  }
  return npeaks;
}


int Trapezoidal_Filter(float* trace, int length, float* filtered, int top, int rise) {
  int i, j, k;
  for (i = 0; i < length; i++) {
    if (i - rise - top >= 0)
    {
      filtered[i] = 0;
      for (j = i - rise; j <= i; j++) {
        filtered[i] += trace[j];
      }
      for (j = i - rise - top; j <= i - top; j++) {
        filtered[i] -= trace[j];
      }
      filtered[i] /= (float)rise;
    }
    else {
      filtered[i] = trace[i];
    }
  }

  return 0;
}

// Exponential moving average filter (low-pass filter)
int EMA_Filter(float* trace, int length, float* filtered, float alpha) {
  int i;
  filtered[0] = 0;
  for (i = 1; i < length; i++) {
    filtered[i] = alpha*trace[i] + (1-alpha)*filtered[i-1];
  }

  // Shift by 1 for better output / input pulse alignment
  for (i = 0; i < length - 1; i++) {
    filtered[i] = filtered[i+1];
  }

  return 0;
}

int SNIP(float* trace, int length, float* filtered, int iterations) {
  int i, j;
  float background[5000];
  for (i = 0; i < length; i++) {background[i] = trace[i];}

  for (i = 1; i < iterations; i++) {
    filtered[0] = background[i];
    filtered[length-1] = background[length-1];
    for (j = 1; j < length - 1; j++)
    {
        if (background[j] <= 0.5*(background[j+1] + background[j-1])) {filtered[j] = background[j];}
        else {filtered[j] = 0.5*(background[j+1] + background[j-1]);}
    }
    for (j = 0; j < length; j++) {background[j] = filtered[j];}
  }
  return 0;
}

int Scan (){


  /** ----------------------------------------------------
   *	Variable declairation
   *   ----------------------------------------------------
   */
  static Can det0, det1;

  bool	beamON,
    trg,
    data;

  float	X, offset, sigma, mu;
  int	multi;

  ifstream fp[16];

  string 	line;

  int i,j,k, pposition,
    Tracelength = 5000;

  float integral;

  int RF_offset = 0;

  float pulse[5000],
    CMAtrace[5000],
    SG_pulse[5000],
    SGderv2_pulse[5000],
    baseline[5000],
    filtered[5000];

  Float_t amplitude,
    risetime,
    falltime,
    width,
    CFD,
    tac,
    paraL,
    paraS,
    runtime,
    bg_avg,
    RF, RF_amp, ebit, meas;

  vector<float> pa, pt, pi;

  bool edge;

  Float_t xy_, x_, y_, x_2, m_, b_;
  int ii;

  Int_t gate;

  // For SG filtered pulse
  Float_t trace_min, trace_max;
  int trace_min_loc, trace_max_loc;
  bool zero_cross;

  char 	filename[250],
    prompt[10],
    openfile[250],
    prefix[250],
    runnum[250],
    interrputPrompt;

  Float_t trgtime, prevtime, difftime;
  Float_t prevtrgtime[10];
  long	TEvt = 0;

  uint64_t buffer64;
  uint32_t buffer32;
  uint16_t buffer16;


  TF1 *f1 = new TF1("f1","gaus",0,350);

  /** ----------------------------------------------------
   *	Calibrations and threshold
   *   ----------------------------------------------------
   */

  float cal[16] =
    {   0.0097,
        0.0092,
        0.0207,
        0.0148,
        0.0206,
        0.0263,
        0.0325,
        0.0376,
        0.031,
        0.0236,
        1.0, 1., 1., 1., 1., 1.
    }; // calibration (keVee / (bit * sample)) from manual check on calibration


  float threshold[16] =
    {   15940,
  	15900,
  	15840,
  	15930,
  	15920,
  	15970,
  	15930,
  	15970,
  	15950,
  	15850,
  	15750,
  	15850,
  	158250,
  	15700,
  	15700,
  	15700
    };

  // Upper electron band gate
  // Extracted from 'psdgate.cpp' using 5-sigma
  float egate[5][3] =
    {
      {2.207, 8.76E-6, 0.108},
      {1.964, 1.27E-6, 0.132},
      {2.094, 4.74E-6, 0.122},
      {2.13E-3, -3.24E-6, 0.397},
      {1.302, -2.05E-5, 0.252}
    };


  /** ----------------------------------------------------
   *	Get functions
   *   ----------------------------------------------------
   */

  PulseAnalysis *Analysis = new PulseAnalysis();


  /** ----------------------------------------------------
   *	Program start...
   *   ----------------------------------------------------
   */

  cout << " ------------------------------------------------ " << endl;
  cout << " | Scan.cpp                                      |" << endl;
  cout << " |   Experiment: ORNL-UTA flash lamp             |" << endl;
  cout << " |   Date: May 29th 2021                         |" << endl;
  cout << " ------------------------------------------------ " << endl;

  cout << "Root file name to be created: ";
  cin >> filename;

  cout << "Run file prefix ('../run#'): ";
  cin >> prefix;

  TFile *ff = new TFile(filename,"RECREATE");

  TTree *tt = new TTree("T","RunData");

  // ODeSA
  tt->Branch("d0",&det0,"l:s:amp:cfd:psd:trg:pp");
  tt->Branch("d1",&det1,"l:s:amp:cfd:psd:trg:pp");
  tt->Branch("meas",&meas,"meas/F");

  // RF Sinewave - channel 13

  // RF logic - channel 14
  tt->Branch("RF",&RF,"RF");

  // Ebit - channel 15
  tt->Branch("ebit",&ebit,"ebit");

  tt->Branch("runtime",&runtime,"runtime/F");     // Runtime in ms

  TH1F *trace0 = new TH1F("trace0","Trace for channel 0",2001,0,8000);
  //tt->Branch("trace0","TH1F", &trace0);
  TH1F *trace1 = new TH1F("trace1","Trace for channel 1",2001,0,8000);
  //tt->Branch("trace1","TH1F", &trace1);
  TH1F *traceTF = new TH1F("traceTF","Trace for channel 1 Trapezoidal Filter",2001,0,8000);
  //tt->Branch("traceTF","TH1F", &traceTF);
  TH1F *traceSNIP = new TH1F("traceSNIP","Trace for SNIP",2001,0,8000);
  //tt->Branch("traceSNIP","TH1F", &traceSNIP);
  TH1F *traceCMA = new TH1F("traceCMA","Trace for CMA",2001,0,8000);
  //tt->Branch("traceCMA","TH1F", &traceCMA);

  TH1F *trace0C = new TH1F("trace0C","Corrected trace for channel 0",600,0,2400);
  //tt->Branch("trace0C","TH1F", &trace0C);

  TH1F *PE = new TH1F("PE","PE distribution",2001,0,10000);
  TH1F *PA = new TH1F("PA","PA distribution",2001,0,2000);

  const int numFiles = 3;

  // Open files
  for (i = 0; i < numFiles; i++)
  {
    /**
    sprintf(openfile, "%s_wave%d.txt",prefix,i);
    cout << " Opening file: " << openfile;
    fp[i].open(openfile, std::ifstream::in);

    if(fp[i].is_open()) {cout << " - Open!" << endl;}
    else {{cout << " - Failed!" << endl;} }
    */

    sprintf(openfile, "%s_wave%d.dat",prefix,i);
    cout << " Opening file: " << openfile;
    fp[i].open(openfile, std::ifstream::in | std::ifstream::binary);

    if(fp[i].is_open()) {cout << " - Open!" << endl;}
    else {{cout << " - Failed!" << endl;} }
  }

  data = 1;
  runtime = 0;
  prevtime = 0;

  /** ----------------------------------------------------
   *	Process liquid scintillator det events
   *   ----------------------------------------------------
   */

  while (data)
    {
      multi = 0;
      X = -1;
      beamON = 0;
      for (j = 0; j < numFiles; j++)
	{
	 if(j > -1)
	 {
	  //if(!fp[j].is_open()){data = 0; cout << "Could not open file!!" << endl;}

    	// Stop after nth events...
    	//if (TEvt > 1000000) {data = 0;}

	  if(fp[j].is_open())
	    {

	      trace_min = 0;

/**
        // Ascii parsing
        if (!getline(fp[j], line)) {data=0; break; } // Record length
	      getline(fp[j], line); // Channel number
	      getline(fp[j], line); // Event number
	      getline(fp[j], line); // Trigger time stamp
*/

        // Binary parsing
        if (!fp[j].read((char*)&buffer32, 4)) {data = 0; break;} //cout << buffer32 << endl;
        Tracelength = (buffer32 - 16)/2; //cout << Tracelength << endl;
        fp[j].read((char*)&buffer32, 4); //cout << buffer32 << endl;
        fp[j].read((char*)&buffer32, 4); //cout << buffer32 << endl;
        fp[j].read((char*)&buffer32, 4); //cout << buffer32 << endl;

        /**
        // Trigger time in 2 ADC clock cycles ( 8 ns )
	      if (j == 0)
            trgtime = atof(line.substr(line.find_last_of(" "),line.length()).c_str());
        */

        // Trigger time in 2 ADC clock cycles ( 8 ns )
        if (j == 0)
            trgtime = buffer32;


	      // Reset variables
	      CFD = -1;
	      amplitude = -1;
	      paraL = 0;
	      paraS = 0;
	      trg = 0;
	      pposition = -1;
        meas = 0;

        if (j == 1)
        {
          pa.clear();
          pt.clear();
          pi.clear();
        }

	      // Get traces
	      for (i = 0; i < Tracelength; i++)
		    {
          /**
          if (!getline(fp[j], line)) {data = 0; break;}
          if (j != 0) {pulse[i] = 16383 - atof(line.c_str()); }
          else {pulse[i] = atof(line.c_str());}

          if (pulse[i] > (16383 - threshold[j])) {trg = 1;}
          */

          if (!fp[j].read((char*)&buffer16, 2)) {data = 0; break;}
          if (j < 12) {pulse[i] = 16383 - (float)buffer16; }
          else {pulse[i] = (float)buffer16;}

          // Added raw traces
          //if (j==0) {trace0->SetBinContent(i, pulse[i]);}
          //if (j==1) {trace1->SetBinContent(i, pulse[i]);}
		    }

	      /** Liquid can processing **/
	      if(Tracelength > 1)
		{


      // Process trace
      if (j < 3) {

        // SiPM
        if (j == 0) {
            Analysis->Baseline_restore(pulse, baseline, Tracelength, 10, 2);
            Analysis->Parameters(pulse, Tracelength, 3, &CFD, &amplitude, &risetime, &falltime, &width);
            for (i = 0; i < Tracelength; i++) { paraL += pulse[i]; }
        }

        // Photodiode
        if (j == 1) {
            Analysis->Baseline_restore(pulse, baseline, Tracelength, 10, 2);
            Analysis->Parameters(pulse, Tracelength, 3, &CFD, &amplitude, &risetime, &falltime, &width);
            for (i = 0; i < Tracelength; i++) { paraL += pulse[i]; }
        }

        // Step logic signal
        if (j == 2) {
            meas = pulse[0];
        }





        // SiPM
        /**
        if (j == 0) {
            Trapezoidal_Filter(pulse, Tracelength, baseline, 3, 3);

            for (i = 0; i < Tracelength; i++)
            {
              pulse[i] = baseline[i];
              traceTF->SetBinContent(i, pulse[i]);
            }

            SNIP(pulse, Tracelength, filtered, 50);

            for (i = 0; i < Tracelength; i++) {traceSNIP->SetBinContent(i, filtered[i]); pulse[i] -= filtered[i];}

            Analysis->CMA_Filter(pulse, Tracelength, CMAtrace, 3, 0, 60);

            for (i = 0; i < Tracelength; i++) {pulse[i] -= CMAtrace[i]; traceCMA->SetBinContent(i, pulse[i]);}

            PeakFinder(pulse, Tracelength, 10, 50, 1950, 50, &pa, &pt, &pi);

            for (i = 0; i < pa.size(); i++)
            {
              PE->Fill(pi[i]);
              PA->Fill(pa[i]);
              paraL += pa[i];
            }
            paraS = (float)pa.size();

            for (i = 0; i < Tracelength; i++) { paraL += pulse[i]; }

        }
        */

        //Analysis->CMA_Filter(pulse, Tracelength, CMAtrace, 10, pulse[0], 3.5 );
        //Analysis->PSD_Integration(pulse, Tracelength, 50, 500, 12, 1, &paraL, &paraS);
        //Trapezoidal_Filter(pulse, Tracelength, baseline, 10, 10);
        //PeakFinder(pulse, Tracelength, 5, 100, 900, 5, &pa, &pt);

      }

    }

  }
	    }


	  switch(j) {
	  case 0 :
	    det0.amp = amplitude;
	    det0.l = paraL;
	    det0.s = paraS;
	    det0.cfd = CFD;
	    det0.trg = trg;
	    det0.psd = det0.s / det0.l;
      det0.pp = pposition;
	    if (det0.trg){ det0.trg = 1;}
	    else {det0.trg = 0;}

	    // Runtime clock
	    difftime = trgtime - prevtime;
	    if(difftime < 0) {
	      difftime = difftime + 2147483647;
	      runtime += ((8*difftime)/1.0E6);
	      prevtime = trgtime;
	    }
	    else {
	      runtime += ((8*difftime)/1.0E6);
	      prevtime = trgtime;
	    }

	    break;

	  case 1 :
	    det1.amp = amplitude;
	    det1.l = paraL;
	    det1.s = paraS;
	    det1.cfd = CFD;
	    det1.trg = trg;
	    det1.psd = det1.s / det1.l;
      det1.pp = pposition;
	    if (det1.trg){ det1.trg = 1;}
	    else {det1.trg = 0;}

	    break;
	  }

    }

      tt->Fill();
      TEvt++;
      if (TEvt%1000==0) {cout << "\rEvent counter: " << TEvt << flush;}

    }


  for (i = 0; i < numFiles; i++)
    {
      if(fp[j].is_open()) fp[j].close();
    }
  ff->cd();
  tt->Write();
  PE->Write();
  PA->Write();
  ff->Close();

  cout << "\nFinsihed! " << TEvt << " events" << endl;

  return 0;
}

int main() {
  return Scan();
}
