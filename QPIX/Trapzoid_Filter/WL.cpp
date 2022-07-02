


#include <iostream>
#include <fstream>
#include <string>
#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TMath.h"
#include "TSpectrum.h"
#include "TF1.h"

typedef struct
{
  Float_t l;               // Long integral sqrt(ll*lr)
  Float_t s;               // Short integral sqrt(sl*sr)
  Float_t amp;             // Amplitude
  Float_t cfd;             // Trigger time
  Float_t psd;             // PSD parameter s/l
  Float_t trg;             // Det trigger
  Float_t tac;             // Zero-crossing timing used for PSD
  Float_t pp;              // Position of max amplitude

} Can;


void WL (){

  Can d0, d1;
  char 	filename[250];

  float wavelength, meas, sum, wls[200], average[200];

  int N, Nsamples, Nsteps;

  bool flag;


  cout << "Root file name to be opened: ";
  cin >> filename;

  TFile *f = new TFile(filename);
  //TFile *f = new TFile("run_386.root");
  TTree *T;
  f->GetObject("T", T);
  T->SetBranchAddress("d0", &d0);
  T->SetBranchAddress("meas", &meas);

  N = T->GetEntries();

  wavelength = 100;
  flag = 1;
  Nsamples = 0;

  for (int i = 0; i < N; i++)
  {
    T->GetEntry(i);


    if (meas > 100 && flag == 1) {

      average[Nsteps] = sum / Nsamples;
      wls[Nsteps] = wavelength;
      cout << average[Nsteps] << " " << wavelength << endl;

      if (wavelength < 110 || wavelength >= 140) {wavelength += 1;}
      else {wavelength += 0.5; }
      flag = 0;
      sum = 0;
      Nsamples = 0;
      Nsteps++;
    } else if (meas < 100 && flag == 0) {
      flag = 1;
    } else if (flag) {
      sum += d0.l;
      Nsamples++;
    }

  }

  TGraph *gr = new TGraph(Nsteps, wls, average);
  gr->Draw("ACP");


  return;
  }
