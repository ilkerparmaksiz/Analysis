#define ProtonDecayAna_cxx
#include "ProtonDecayAna.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include "TGraph.h"
#include "TGraph2D.h"
#include "TH1F.h"
#include "TNtuple.h"
#include <random>
#include <iomanip>
#include <string>
#include <map>
#include <random>
#include <cmath>
#include "TH3D.h"


using namespace std;
void ProtonDecayAna::Getdedx(Int_t Pdg_Code,Int_t limit,Int_t PionHitLimit=2) {


    if (fChain == 0) return;
    Long64_t nentries = fChain->GetEntriesFast();
    Long64_t nbytes = 0, nb = 0;
    std::vector<double> range;
    std::vector<double> dedx;
    TH1F *h2 =new TH1F("h2","h2",100,0,600);

    for (Long64_t jentry = 0; jentry < nentries; jentry++) {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);
        nbytes += nb;

        int piontrackIds;
        int CountPion=0;

        float Range;
        for (int i = 0; i < number_particles; i++) {
            if (particle_pdg_code->at(i) == Pdg_Code) {
                piontrackIds = particle_track_id->at(i) ;
                CountPion++;
            }

        }
        double hitE = 0;
        double hitlength=0;
        double fde=0;

        std::vector<double> E_;
        std::vector<double> X_;
        std::vector<double> Y_;
        std::vector<double> Z_;
        std::vector<double> T_;

        std::vector<double> length_;

        int CountPionHits=0;
        Double_t TotalEnergy=0;
        for (int i = 0; i < number_hits; i++) {
            if (hit_track_id->at(i) == piontrackIds) {
                //dedx.push_back(hit_energy_deposit->at(i)/hit_length->at(i));

                /*E_.push_back(hit_energy_deposit->at(i));
                X_.push_back(hit_start_x->at(i));
                Y_.push_back(hit_start_y->at(i));
                Z_.push_back(hit_start_z->at(i));
                 */
                //cout<<"E -> " <<hit_energy_deposit->at(i)<< " L -> " << hit_length->at(i)<<endl;
                // cout<<"T "<< hit_start_t->at(i)<<" X -> " <<hit_start_x->at(i)<< " Y -> " << hit_start_y->at(i)<< " Z -> " << hit_start_z->at(i)<<endl;
                E_.push_back(hit_energy_deposit->at(i));
                X_.push_back(hit_start_x->at(i));
                Y_.push_back(hit_start_y->at(i));
                Z_.push_back(hit_start_z->at(i));

                length_.push_back(hit_length->at(i));

                CountPionHits++;
                TotalEnergy+=hit_energy_deposit->at(i);

            }


            // cout<<"Event-> " <<jentry<<endl;
            //cout<<"HitDeposit-> "<<hitE<<endl;
            //h1->Fill(hitE);

        }

        float dx;
        float tempHitLength;

        if(CountPionHits>PionHitLimit){
            if(E_.size()==0 || X_.size()==0){ std::cout<<"Empty Vector"<<std::endl; continue;}
            h2->Fill(TotalEnergy);
            for(int i=0;i<E_.size()-1;i++){
                //cout<<"de/dx -> "<<(E_.at(i)-E_.at(i+1))/(X_.at(i)-X_.at(i+1))<<endl;
                //dx=sqrt((X_.at(i)-X_.at(i+1))*(X_.at(i)-X_.at(i+1)) + (Y_.at(i)-Y_.at(i+1))*(Y_.at(i)-Y_.at(i+1))+(Z_.at(i)-Z_.at(i+1))*(Z_.at(i)-Z_.at(i+1)) );
                //dedx.push_back(abs(E_.at(i)-E_.at(i+1))/dx);
                //dedx.push_back((E_.at(i))/dx);
                dedx.push_back((E_.at(i))/length_.at(i));
                tempHitLength=0;

                for (int k=i;k<X_.size()-1;k++){
                    //tempHitLength+=sqrt((X_.at(k)-X_.at(k+1))*(X_.at(k)-X_.at(k+1)) + (Y_.at(k)-Y_.at(k+1))*(Y_.at(k)-Y_.at(k+1))+(Z_.at(k)-Z_.at(k+1))*(Z_.at(k)-Z_.at(k+1)) );
                    tempHitLength+=length_.at(k);
                }

                range.push_back(tempHitLength);
            }

        }

        if (limit!=0 && jentry==limit) break;

    }
    std::string Tittle;
    Tittle="Hits for PDgCode -> " + std::to_string(Pdg_Code);
    cout<<Tittle<<endl;

    TCanvas *c1= new TCanvas("c1",Tittle.c_str(),1024,800);
    c1->Divide(1,2);
    TGraph *gr;
    gr = new TGraph(dedx.size(), &range[0], &dedx[0]);
    gr->SetTitle("DeDx vs  ResidualRange");
    gr->GetXaxis()->SetTitle("Residual Range");
    gr->GetYaxis()->SetTitle("dedx");
    gr->SetMinimum(0);
    gr->SetMaximum(60);
    gr->GetXaxis()->SetLimits(0,40);
    c1->cd(1);
    gr->Draw("AP");
    c1->cd(2);
    h2->Draw();

}

void ProtonDecayAna::Loop() {
//   In a ROOT session, you can do:
//      root> .L ProtonDecayAna.C
//      root> ProtonDecayAna t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch

    if (fChain == 0) return;

    Long64_t nentries = fChain->GetEntriesFast();
//TH2F *h2=new TH2F("h","Kaon",100,0,100,100,0,100);
    TGraph *gr;
    Long64_t nbytes = 0, nb = 0;
    std::vector<double> range;
    std::vector<double> dedx;
    std::vector<HitMap*> hits;
    for (Long64_t jentry = 0; jentry < nentries; jentry++) {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);
        nbytes += nb;
        // if (Cut(ientry) < 0) continue;
        /*cout << "NumberofParticles -> " << number_particles << endl;
        cout << "PdgCodeSize -> " << particle_pdg_code->size() << endl;
        cout << "HitDepositSize -> " << hit_energy_deposit->size() << endl;
        cout << " NumberofHits ->" << number_hits << endl;
        */
        HitMap *hit =new HitMap();
        hit->Pdg_code=321;
        hit->EventID=jentry;
         int piontrackIds;
        int CountPion=0;

        float Range;
        for (int i = 0; i < number_particles; i++) {
            if (particle_pdg_code->at(i) == 321) {
                piontrackIds = particle_track_id->at(i) ;
                CountPion++;

            }


        }
        double hitE = 0;
        double hitlength=0;
        double fde=0;

        std::vector<double> E_;
        std::vector<double> X_;
        std::vector<double> Y_;
        std::vector<double> Z_;
        std::vector<double> T_;
        std::vector<double> length_;

        int CountPionHits=0;
        for (int i = 0; i < number_hits; i++) {

            if (hit_track_id->at(i) == piontrackIds) {
                //dedx.push_back(hit_energy_deposit->at(i)/hit_length->at(i));

                /*E_.push_back(hit_energy_deposit->at(i));
                X_.push_back(hit_start_x->at(i));
                Y_.push_back(hit_start_y->at(i));
                Z_.push_back(hit_start_z->at(i));
                length_.push_back(hit_length->at(i));
                 */
                //cout<<"E -> " <<hit_energy_deposit->at(i)<< " L -> " << hit_length->at(i)<<endl;
               // cout<<"T "<< hit_start_t->at(i)<<" X -> " <<hit_start_x->at(i)<< " Y -> " << hit_start_y->at(i)<< " Z -> " << hit_start_z->at(i)<<endl;
                hit->DepositedEE.push_back(hit_energy_deposit->at(i));
                hit->StartX.push_back(hit_start_x->at(i));
                hit->StartY.push_back(hit_start_y->at(i));
                hit->StartZ.push_back(hit_start_z->at(i));
                CountPionHits++;

            }


            // cout<<"Event-> " <<jentry<<endl;
            //cout<<"HitDeposit-> "<<hitE<<endl;
            //h1->Fill(hitE);
            hits.push_back(hit);


        }
        /*float dx;
        float tempHitLength;
        if(CountPionHits>14){

            for(int i=0;i<E_.size()-1;i++){

                //cout<<"de/dx -> "<<(E_.at(i)-E_.at(i+1))/(X_.at(i)-X_.at(i+1))<<endl;
                //dx=sqrt((X_.at(i)-X_.at(i+1))*(X_.at(i)-X_.at(i+1)) + (Y_.at(i)-Y_.at(i+1))*(Y_.at(i)-Y_.at(i+1))+(Z_.at(i)-Z_.at(i+1))*(Z_.at(i)-Z_.at(i+1)) );
                //dedx.push_back(abs(E_.at(i)-E_.at(i+1))/dx);
                //dedx.push_back((E_.at(i))/dx);
                dedx.push_back(E_.at(i)/length_.at(i));
                tempHitLength=0;

                for (int k=i;k<X_.size()-1;k++){
                    //tempHitLength+=sqrt((X_.at(k)-X_.at(k+1))*(X_.at(k)-X_.at(k+1)) + (Y_.at(k)-Y_.at(k+1))*(Y_.at(k)-Y_.at(k+1))+(Z_.at(k)-Z_.at(k+1))*(Z_.at(k)-Z_.at(k+1)) );
                    tempHitLength+=length_.at(k);
                }

                range.push_back(tempHitLength);
            }

        }*/
        /*if (jentry==2000)
            break;
        */

    }
    gr = new TGraph(dedx.size(), &range[0], &dedx[0]);
    gr->SetMinimum(0);
    gr->SetMaximum(60);
    gr->GetXaxis()->SetLimits(0,40);

    gr->Draw("AP");
}


// Voxalizing the hits

std::map<Long64_t,Voxel*>ProtonDecayAna::VoxelizetheHits(Int_t Event_ID,std::string Phase) {
    f->cd();
    double Wvalue,E_vel,DiffL,DiffT,Life_Time,Reset,ElectronCharge;
    if(Phase=="gas"){

        //Fix these values for the gas
        Wvalue=23.6;
        E_vel=260000.0;
        DiffL=6.8223;
        DiffT=13.1586;
        Life_Time=0.1;
        Reset=6250;
        ElectronCharge=1.60217662e-19;
    }else{
        Wvalue=23.6;
        E_vel=164800.0;
        DiffL=6.8223;
        DiffT=13.1586;
        Life_Time=0.1;
        Reset=6250;
        ElectronCharge=1.60217662e-19;
    }

    HitMap*  hit= GetSingleEventHits(Event_ID);
    std::map<Long64_t ,Voxel*> voxels;
    std::vector<Electron*> Electrons;
   // if (hit->StartX.size()==0) {std::cout<<"StartX is empty .."; return voxels; }

    Int_t indexer=0;
    double voxel_X_step=0.4; //mm
    double voxel_Y_step=0.4; //mm
    double voxel_Z_step=0.4; //mm
    for (int ient = 0; ient < hit->StartX.size(); ient++) {


        // from PreStepPoint
        double const start_x = hit->StartX.at(ient);      // cm
        double const start_y = hit->StartY.at(ient);      // cm
        double const start_z = hit->StartZ.at(ient);      // cm
        double const start_t = (hit->StartT.at(ient)* 1e-9); // nsec

        // from PostStepPoint
        double const end_x = hit->EndX.at(ient);      // cm
        double const end_y = hit->EndY.at(ient);      // cm
        double const end_z = hit->EndZ.at(ient);      // cm
        double const end_t = (hit->EndT.at(ient) * 1e-9); // nsec
        int Nelectron;
        // to account for recombination or not
        // calcualte the number of electrons in the hit
        /* if (Recombination)
         {
             Nelectron = round(Recombonation * (energy_deposit*1e6/Wvalue) );
         }else
         {
             Nelectron = round( (energy_deposit*1e6/Wvalue) );
         }*/
        Nelectron = round((hit->DepositedEE.at(ient) * 1e6 / Wvalue));


        // if not enough move on
        if (Nelectron == 0) { continue; }
        //cout<<"EventID -> "<<Event_ID<<" Number of Electrons -> " << Nelectron <<endl;
        // define the electrons start position
        double electron_loc_x = hit->StartX.at(ient);
        double electron_loc_y = hit->StartY.at(ient); ;
        double electron_loc_z = hit->StartZ.at(ient); ;
        double electron_loc_t = hit->StartT.at(ient); ;

        // Determin the "step" size (pre to post hit)
        double const step_x = (end_x - start_x) / Nelectron;
        double const step_y = (end_y - start_y) / Nelectron;
        double const step_z = (end_z - start_z) / Nelectron;
        double const step_t = (end_t - start_t) / Nelectron;

        double electron_x, electron_y, electron_z;
        double T_drift, sigma_L, sigma_T;


        // Loop through the electrons


        for (int i = 0; i < Nelectron; i++) {
            // calculate drift time for diffusion
            T_drift = electron_loc_z / E_vel;
            // electron lifetime

            if (drand48() >= exp(-T_drift / Life_Time)) { continue; }

            // diffuse the electrons position


            sigma_T = sqrt(2 * DiffT * T_drift);
            sigma_L = sqrt(2 * DiffL * T_drift);


            std::default_random_engine generator;

            std::normal_distribution<double> E_x(electron_loc_x, sigma_T);
            electron_x =E_x(generator);
            std::normal_distribution<double>E_y(electron_loc_y, sigma_T);
            electron_y =E_y(generator);

            std::normal_distribution<double>E_z(electron_loc_z, sigma_L);
            electron_z =E_z(generator);

            // add the electron to the vector.
            Electrons.push_back(new Electron());

            // convert the electrons x,y to a pixel index
            /* int Pix_Xloc, Pix_Yloc;
             Pix_Xloc = (int) ceil(electron_x / Pix_Size);
             Pix_Yloc = (int) ceil(electron_y / Pix_Size);
             */

            Electrons[indexer]->Pix_ID = 1; //When adding sensors use some pixel identification as above
            Electrons[indexer]->time = electron_loc_t + (electron_z / E_vel);

            // Move to the next electron
            electron_loc_x += step_x;
            electron_loc_y += step_y;
            electron_loc_z += step_z;
            electron_loc_t += step_t;
            double VoxelX,VoxelY,VoxelZ;
            Long64_t VoxelID;

            // Create Voxels
            VoxelX=floor(electron_x/voxel_X_step);
            VoxelY=floor(electron_y/voxel_Y_step);
            VoxelZ=floor(electron_z/voxel_Z_step);


            VoxelID=(1e6*VoxelX)+(1e3*VoxelY)+(VoxelZ);
            // Generate Voxels..
            //VoxelFloor(electron_x,electron_y,electron_z,voxel_X_step,voxel_Y_step,voxel_Z_step,1,voxels );

            //cout<< voxels.size() <<" are created .. " <<endl;

            if(voxels.count(VoxelID)) {

                voxels[VoxelID]->Q+=1;

                //Find the Index of the Current PDG
                //std::vector<int>::iterator itr = std::find(voxels[VoxelID]->PDGCode.begin(), voxels[VoxelID]->PDGCode.end(), hit->Pdg_code);
                //voxels[VoxelID]->QPDg[std::distance(voxels[VoxelID]->PDGCode.begin(), itr)]+=1;
            }
            else {
                Voxel * v1=new Voxel();
                v1->VoxelID=VoxelID;
                v1->Q=1;
                v1->X= VoxelX;
                v1->Y=VoxelY;
                v1->Z=VoxelZ;
                v1->PDGCode.push_back(hit->Pdg_code);
                v1->QPDg.push_back(1);
                voxels.insert(pair<Long64_t,Voxel*>(VoxelID,v1));
            }
            indexer += 1;

        }

    }

    return voxels;
}

std::vector<HitMap*> ProtonDecayAna::GetAllHits(Int_t Pdg_code,Long64_t NumOfEvents=1){

    std::vector<HitMap*> hits;
    if (fChain == 0) {cout<<"Cant Read the Root File.. "<<endl; return hits;}
    Long64_t nentries = fChain->GetEntriesFast();

    if(NumOfEvents<1 or nentries>NumOfEvents) NumOfEvents=nentries;
//TH2F *h2=new TH2F("h","Kaon",100,0,100,100,0,100);
    Long64_t nbytes = 0, nb = 0;

    for (Long64_t jentry = 0; jentry < NumOfEvents; jentry++) {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);
        nbytes += nb;
        // if (Cut(ientry) < 0) continue;
        /*cout << "NumberofParticles -> " << number_particles << endl;
        cout << "PdgCodeSize -> " << particle_pdg_code->size() << endl;
        cout << "HitDepositSize -> " << energy_deposit->size() << endl;
        cout << " NumberofHits ->" << number_hits << endl;
        */
        HitMap *hit =new HitMap();
        hit->Pdg_code=Pdg_code;
        hit->EventID=jentry;
        int piontrackIds;
        int CountPion=0;

        float Range;
        for (int i = 0; i < number_particles; i++) {
            if (particle_pdg_code->at(i) == Pdg_code) {
                piontrackIds = particle_track_id->at(i) ;
                CountPion++;

            }


        }
        double hitE = 0;
        double hitlength=0;
        double fde=0;

        std::vector<double> E_;
        std::vector<double> X_;
        std::vector<double> Y_;
        std::vector<double> Z_;
        std::vector<double> T_;
        std::vector<double> length_;

        int CountPionHits=0;
        for (int i = 0; i < number_hits; i++) {

            if (hit_track_id->at(i) == piontrackIds) {
                hit->EventID=jentry;
                hit->TrackID=piontrackIds;
                hit->Pdg_code=Pdg_code;
                hit->DepositedEE.push_back(hit_energy_deposit->at(i));
                hit->StartX.push_back(hit_start_x->at(i));
                hit->StartY.push_back(hit_start_y->at(i));
                hit->StartZ.push_back(hit_start_z->at(i));
                hit->StartT.push_back(hit_start_t->at(i));
                hit->EndX.push_back(hit_end_x->at(i));
                hit->EndY.push_back(hit_end_y->at(i));
                hit->EndZ.push_back(hit_end_z->at(i));
                hit->EndT.push_back(hit_end_t->at(i));
                CountPionHits++;
            }

            hits.push_back(hit);


        }


    }
    return hits;

}

HitMap* ProtonDecayAna::GetSingleEventHits(Long64_t jentry){
    std::map<int,HitMap*> Hitmaplist;
    std::vector<HitMap*> hits;
    std::vector<int> ParticleTrackIDs;
    if (fChain == 0) {cout<<"Cant Read the Root File.. "<<endl; return 0;}

//TH2F *h2=new TH2F("h","Kaon",100,0,100,100,0,100);
    Long64_t nbytes = 0, nb = 0;

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) return 0;
    nb = fChain->GetEntry(jentry);
    nbytes += nb;

    HitMap *hit =new HitMap();
    hit->EventID=jentry;
    int piontrackIds;
    int CountPion=0;

    float Range;
    for (int i = 0; i < number_particles; i++) {
        if (std::find(ParticleTrackIDs.begin(), ParticleTrackIDs.end(), particle_track_id->at(i)) == ParticleTrackIDs.end()) {
            // particle_track_id->at(i) not in ParticleTrackIDs, add it
            ParticleTrackIDs.push_back(particle_track_id->at(i));
        }

    }
    double hitE = 0;
    double hitlength=0;
    double fde=0;

    std::vector<double> E_;
    std::vector<double> X_;
    std::vector<double> Y_;
    std::vector<double> Z_;
    std::vector<double> T_;
    std::vector<double> length_;

    int CountPionHits=0;
    hit->EventID=jentry;

    for (int i = 0; i < number_hits; i++) {


            hit->TrackID.push_back(particle_track_id->at(i));
            hit->Pdg_code.push_back(particle_pdg_code->at(i));
            hit->DepositedEE.push_back(hit_energy_deposit->at(i));
            hit->StartX.push_back(hit_start_x->at(i));
            hit->StartY.push_back(hit_start_y->at(i));
            hit->StartZ.push_back(hit_start_z->at(i));
            hit->StartT.push_back(hit_start_t->at(i));
            hit->EndX.push_back(hit_end_x->at(i));
            hit->EndY.push_back(hit_end_y->at(i));
            hit->EndZ.push_back(hit_end_z->at(i));
            hit->EndT.push_back(hit_end_t->at(i));
            CountPionHits++;


    }
    return hit;

}


void ProtonDecayAna::PlotHits(Int_t EventID,Int_t Pdg_Code) {

    if (fChain == 0) {cout<<"Cant Read the Root File.. "<<endl; return;};
    TGraph *gr;
    std::vector<double> range;
    std::vector<double> dedx;
    std::vector<HitMap*> hits;
    LoadTree(EventID);
    fChain->GetEntry(EventID);

    int piontrackIds;
    int CountPion=0;
    TNtuple *np =new TNtuple("Hits","Hits","x:y:z:E");
    for (int i = 0; i < number_particles; i++) {
        if (particle_pdg_code->at(i) == Pdg_Code)
        {
            piontrackIds = particle_track_id->at(i) ;
            CountPion++;
        }
    }

    std::vector<double> E_;
    std::vector<double> X_;
    std::vector<double> Y_;
    std::vector<double> Z_;
    std::vector<double> T_;
    std::vector<double> length_;

    int CountPionHits=0;
    for (int i = 0; i < number_hits; i++) {
        if (hit_track_id->at(i) == piontrackIds) {
            E_.push_back(hit_energy_deposit->at(i));
            X_.push_back(hit_start_x->at(i));
            Y_.push_back(hit_start_y->at(i));
            Z_.push_back(hit_start_z->at(i));
            CountPionHits++;
            np->Fill((float)hit_start_x->at(i),(float)hit_start_y->at(i),(float)hit_start_z->at(i),(float)(hit_energy_deposit->at(i)*10));
        }

    }
    if(CountPionHits<4){
        std::cout<<"Please choose an other event, this one has less hits"<<std::endl;
        return;
    }
    std::string Tittle;
    Tittle="Hits for EventID-> " + std::to_string(EventID) + " PDgCode ->" + std::to_string(Pdg_Code);
    TCanvas *c1= new TCanvas("c1",Tittle.c_str(),1024,800);
    c1->Divide(1,2);
    c1->cd(1);
    gStyle->SetPalette(1);
    np->SetNameTitle("np1",Tittle.c_str());

    np->Draw("x:y:z:E","1","COL");
    c1->cd(2);
    TGraph2D *gr2d=new TGraph2D(X_.size(),&X_[0],&Y_[0],&Z_[0]);
    gStyle->SetPalette(1);
    gr2d->SetMarkerStyle(20);
    gr2d->Draw("pcol");

}

void ProtonDecayAna::PlotVoxels(Int_t EventID,std::string fileName="/home/argon/Projects/Ilker/QPIX_Develop/Analysis/Voxels.root") {

    TTree *Voxels;
    Voxels=GetVoxelsFromFile(EventID,fileName);

    Voxels->GetEntry(EventID);
    std::string Title="Voxels_"+std::to_string(EventID);
    gStyle->SetCanvasPreferGL(kTRUE);

    double maxx,maxy,maxz,minx,miny,minz;
    minx = *min_element(fVoxel_x->begin(), fVoxel_x->end())-2;
    miny = *min_element(fVoxel_y->begin(), fVoxel_y->end())-2;
    minz = *min_element(fVoxel_z->begin(), fVoxel_z->end())-2;
    maxx = *max_element(fVoxel_x->begin(), fVoxel_x->end())+2;
    maxy = *max_element(fVoxel_y->begin(), fVoxel_y->end())+2;
    maxz = *max_element(fVoxel_z->begin(), fVoxel_z->end())+2;
    TH3D *h3=new TH3D("h3",Title.c_str(),maxx/0.4,minx,   maxx,maxy/0.4,miny,maxy,maxz/0.4,minz,maxz);



    gStyle->SetPalette(1);
    if(!fVoxel_x->empty()){
        for (int i=0;i<fVoxel_x->size();i++)
            h3->Fill(fVoxel_x->at(i),fVoxel_y->at(i),fVoxel_z->at(i),fVoxel_Q->at(i));
        h3->Draw("glbox4");
    }else {cout<<"Vector Information is empty for "<< fEvent_ID<<endl; return;}



}
TTree* ProtonDecayAna::GetVoxelsFromFile(Int_t EventID, std::string fileName) {
        TFile *tf;
        TTree *tt=0;
        if (tt == 0) {
            tf = (TFile*)gROOT->GetListOfFiles()->FindObject(fileName.c_str());
            if (!tf || !tf->IsOpen()) {
                tf = new TFile(fileName.c_str());
            }
            tf->GetObject("voxels",tt);
        }
        tt->SetBranchAddress("EventID",&fEvent_ID,&bEvent_ID);
        tt->SetBranchAddress("VoxelID",&fVoxelID,&bVoxelID);
        tt->SetBranchAddress("x",&fVoxel_x,&bVoxel_x);
        tt->SetBranchAddress("y",&fVoxel_y,&bVoxel_y);
        tt->SetBranchAddress("z",&fVoxel_z,&bEvent_ID);
        tt->SetBranchAddress("q",&fVoxel_Q,&bVoxel_Q);
    return tt;
}

void ProtonDecayAna::VoxelFloor(double XTrue,double YTrue,double ZTrue,double XStep,double YStep,double ZStep,double StartingSize,std::map<Long64_t,Voxel*> &Voxels ){
    //cout<<"Creating the Voxels ..."<<endl;
    // some Parameters for Voxels

    double XStart,YStart,ZStart,XSize,YSize,ZSize;
    Long64_t VoxelID;

    XStart=(XTrue-StartingSize);
    YStart=(YTrue-StartingSize);
    ZStart=(ZTrue-StartingSize);

    //VoxelID
    for (double vx=XStart;vx<= (XTrue+StartingSize);vx+=XStep){

        for (double vy=YStart;vy<= (YTrue+StartingSize);vy+=YStep){

            for (double vz=ZStart;vz<=(ZTrue+StartingSize);vz+=ZStep){
                VoxelID=ceil((1e4*vx)+(1e3*vy)+(1e2*vz));
                if(Voxels.count(VoxelID)){
                    continue;

                }else{

                    //cout<<vx <<"  "<<vy<<" "<<vz<<" " <<ZSize<<endl;
                    Voxel * v1=new Voxel();
                    v1->VoxelID=VoxelID;
                    v1->Q=0;
                    v1->X=vx;
                    v1->Y=vy;
                    v1->Z=vz;
                    Voxels.insert(pair<Long64_t,Voxel*>(VoxelID,v1));
                    //cout<<VoxelID << " is created " <<endl;
                }

            }
        }
    }



}



// Saving to File
void ProtonDecayAna::InitializeFile(string Name){
    vxfile=TFile::Open(  Name.c_str(),"RECREATE");
    vxtree =new TTree("voxels","Voxels for ProtonDecay");
    vxtree->Branch("EventID",&Event_ID,"EventID/I");
    vxtree->Branch("PdgCode",&Voxel_PdgCode);
    vxtree->Branch("QPDG",&Voxel_QPDG);
    vxtree->Branch("VoxelID",&VoxelID);
    vxtree->Branch("x",&Voxel_x);
    vxtree->Branch("y",&Voxel_y);
    vxtree->Branch("z",&Voxel_z);
    vxtree->Branch("q",&Voxel_Q);
}
void ProtonDecayAna::Clear(){
    Event_ID=-1;
    VoxelID.clear();
    Voxel_x.clear();
    Voxel_y.clear();
    Voxel_z.clear();
    Voxel_Q.clear();
    Voxel_PdgCode.clear();
    Voxel_QPDG.clear();
}
void ProtonDecayAna::SavetoFile(){
    vxtree->Write();

}

void ProtonDecayAna::SaveVoxelsToFile(Int_t StartEvent,Int_t EndEvent,std::string Phase="liquid"){
    //Create the file
    if(!((EndEvent-StartEvent)<=fChain->GetEntries() && (EndEvent-StartEvent)>0)) {
        cout<<"Choose a range less than what you have in the file"<<endl;
        return;
    }
    InitializeFile("Voxelsv2.root");
    Int_t count=0;
    for(int i=StartEvent;i<EndEvent;i++){
        std::map<Long64_t,Voxel*> voxels;
        voxels=VoxelizetheHits(i,Phase);

        if(voxels.size()!=0){
            Clear();
            Event_ID=i;
            for (auto &x:voxels){
                if(x.second->Q>=1){
                VoxelID.push_back(x.second->VoxelID);
                Voxel_PdgCode=x.second->PDGCode;
                Voxel_QPDG=x.second->QPDg;
                Voxel_x.push_back(x.second->X);
                Voxel_y.push_back(x.second->Y);
                Voxel_z.push_back(x.second->Z);
                Voxel_Q.push_back(x.second->Q);
                count++;
                } else continue;
            }

            if(!Voxel_x.empty() && !Voxel_y.empty() && !Voxel_z.empty() && !Voxel_Q.empty() ){
                cout<<"Voxel Size -> "<<Voxel_x.size()<<endl;

                vxfile->cd();
                vxtree->Fill();

            }else {
                std::cout<<"some of the vectors are empty"<<endl;
            }

        }

    }
    if(count>0){
        SavetoFile();
        vxfile->Close();
    }






}