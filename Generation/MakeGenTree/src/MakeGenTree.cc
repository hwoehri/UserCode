// -*- C++ -*-
//
// Package:    MakeGenTree
// Class:      MakeGenTree
// 
/**\class MakeGenTree MakeGenTree.cc Workspace/MakeGenTree/src/MakeGenTree.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Hermine Woehri,40 1-A32,+41227679747,
//         Created:  Thu Nov 25 15:20:12 CET 2010
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <SimDataFormats/GeneratorProducts/interface/HepMCProduct.h>

#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
//
// class declaration
//

using namespace std;
using namespace ROOT;

class MakeGenTree : public edm::EDAnalyzer {
   public:
      explicit MakeGenTree(const edm::ParameterSet&);
      ~MakeGenTree();


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------
private:
  edm::InputTag genParticlesLabel_;
  int QQbarPDG_;

  string outputFileName_;
  TFile *fOut_;
  TTree *tree_;
  TLorentzVector *jPsi_4mom_, *muPos_4mom_, *muNeg_4mom_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
MakeGenTree::MakeGenTree(const edm::ParameterSet& iConfig) :
  genParticlesLabel_(iConfig.getParameter<edm::InputTag>("genParticlesLabel")),
  QQbarPDG_(iConfig.getParameter<int>("QQbarPDG")),
  outputFileName_(iConfig.getParameter<string>("outputFileName")){
   //now do what ever initialization is needed

}


MakeGenTree::~MakeGenTree()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  delete jPsi_4mom_;
  delete muPos_4mom_;
  delete muNeg_4mom_;  

  fOut_->Close();
//   delete tree_;
  delete fOut_;
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
MakeGenTree::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

   using namespace edm;

   Handle<HepMCProduct> genEventHandle;
   iEvent.getByLabel(genParticlesLabel_, genEventHandle);
   const HepMC::GenEvent* genEvent = genEventHandle->GetEvent();

   for (HepMC::GenEvent::particle_const_iterator it = genEvent->particles_begin(); it != genEvent->particles_end(); ++it){
     HepMC::FourVector momP = (*it)->momentum();
     if((*it)->pdg_id() == 13){//neg. muon
       muNeg_4mom_->SetPxPyPzE(momP.px(), momP.py(), momP.pz(), momP.e());
     }
     else if((*it)->pdg_id() == -13){//pos. muon
       muPos_4mom_->SetPxPyPzE(momP.px(), momP.py(), momP.pz(), momP.e());
     }
     else if ((*it)->pdg_id() == QQbarPDG_){//quarkonium
       jPsi_4mom_->SetPxPyPzE(momP.px(), momP.py(), momP.pz(), momP.e());
//        HepMC::ThreeVector prodVtx = (*it)->production_vertex()->point3d();
//        originQQbar_->SetXYZ(prodVtx.x(), prodVtx.y(), prodVtx.z());
//        cout << "origin: " << prodVtx.x() << " " << prodVtx.y() <<  " " << prodVtx.z() << endl;
//        HepMC::ThreeVector decayVtx = (*it)->end_vertex()->point3d();
//        decayVtxQQbar_->SetXYZ(decayVtx.x(), decayVtx.y(), decayVtx.z());
//        cout << "decay of " << QQbarPDG_ << ": " << prodVtx.x() << " " << prodVtx.y() <<  " " << prodVtx.z() << endl;

//        QQbarTheta_ = (*it)->polarization().theta();
//        QQbarPhi_ = (*it)->polarization().phi();
       //       cout << QQbarPhi_ << " was generated with theta = " << QQbarTheta_ << " and phi = " << QQbarPhi_ << endl;
     }
   }

   tree_->Fill();

   //reset the variables for the next event
   jPsi_4mom_->SetPxPyPzE(-9999, -9999, -9999, -9999);
   muPos_4mom_->SetPxPyPzE(-9999, -9999, -9999, -9999);
   muNeg_4mom_->SetPxPyPzE(-9999, -9999, -9999, -9999);

}

// ------------ method called once each job just before starting event loop  ------------
void 
MakeGenTree::beginJob()
{
  fOut_ = new TFile(outputFileName_.c_str(), "RECREATE");

  jPsi_4mom_ = new TLorentzVector();
  muPos_4mom_ = new TLorentzVector();
  muNeg_4mom_ = new TLorentzVector();

  tree_ = new TTree("QQbarGenTree", "generator level info for acceptance calculation");

  //single muon variables:
  tree_->Branch("jPsi_4mom", "TLorentzVector", &jPsi_4mom_);
  tree_->Branch("muPos_4mom", "TLorentzVector", &muPos_4mom_);
  tree_->Branch("muNeg_4mom", "TLorentzVector", &muNeg_4mom_);

   //set default values
   jPsi_4mom_->SetPxPyPzE(-9999, -9999, -9999, -9999);
   muPos_4mom_->SetPxPyPzE(-9999, -9999, -9999, -9999);
   muNeg_4mom_->SetPxPyPzE(-9999, -9999, -9999, -9999);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MakeGenTree::endJob() {

  fOut_->cd();
  tree_->Write();
}

//define this as a plug-in
DEFINE_FWK_MODULE(MakeGenTree);
