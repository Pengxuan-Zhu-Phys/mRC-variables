// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/VisibleFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Tools/Cutflow.hh"
#include <tuple>
#include <iostream>
#include <fstream>

namespace Rivet
{

  /// @brief Add a short analysis description here
  class CEPC_MASS : public Analysis
  {
  public:
    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CEPC_MASS);

    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init()
    {

      // Initialise and register projections

      // The basic final-state projection:
      // all final-state particles within
      // the given eta acceptance
      std::ofstream df;
      // dfname = getOption("DFNAME");
      // dfname.append(".csv");

      // df.open(dfname, std::ios::out | std::ios::app);
      // MSG_INFO("Open file: " << dfname );
      // df << "Nmuon,Emiss,Njet,mVV,mRecoil,mMiss,mRCmin,mRCmax,mRCLSP,mLSPmax\n";
      // Initialise and register projections

      // The basic final-state projection:
      // all final-state particles within
      // the given eta acceptance
      declare(Beam(), "Beams");

      const FinalState fs(Cuts::abseta < 3.0);
      declare(VisibleFinalState(Cuts::abseta < 3.0), "vfs");

      Cut lepton_cuts = ((Cuts::abseta < 3.0) && (Cuts::E > 0.5 * GeV));

      FinalState bare_elec(Cuts::abspid == PID::ELECTRON && lepton_cuts);
      FinalState bare_muon(Cuts::abspid == PID::MUON && lepton_cuts);
      // PromptFinalState bare_elec(Cuts::abspid == PID::ELECTRON && lepton_cuts);
      // PromptFinalState bare_muon(Cuts::abspid == PID::MUON && lepton_cuts);

      // The final-state particles declared above are clustered using FastJet with
      // the anti-kT algorithm and a jet-radius parameter 0.4
      // muons and neutrinos are excluded from the clustering
      FastJets jetfs(fs, FastJets::ANTIKT, 0.4, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      declare(jetfs, "jets");

      // FinalState of prompt photons and bare muons and electrons in the event
      // FinalState photons(Cuts::abspid == PID::PHOTON);
      PromptFinalState photons(Cuts::abspid == PID::PHOTON);

      // Dress the prompt bare leptons with prompt photons within dR < 0.1,
      // and apply some fiducial cuts on the dressed leptons
      // DressedLeptons elec_drs(photons, bare_elec, 0.1, lepton_cuts);
      // DressedLeptons muon_drs(photons, bare_muon, 0.1, lepton_cuts);
      DressedLeptons elec_drs(photons, bare_elec, 0.1, lepton_cuts);
      DressedLeptons muon_drs(photons, bare_muon, 0.1, lepton_cuts);
      declare(muon_drs, "muons");
      declare(elec_drs, "elecs");
      // declare(bare_elec, "elecs");
      // declare(bare_muon, "muons");

      // Missing momentum
      declare(MissingMomentum(fs), "MET");

      // Book histograms
      // specify custom binning
      book(_hist_Wmin, "hist_Wmass_mRCmin", 50, 75, 85);
      book(_hist_Wmax, "hist_Wmass_mRCmax", 50, 75, 85);
      book(_Norm_hist_Wmin, "Nhist_Wmass_mRCmin", 50, 75, 85);
      book(_Norm_hist_Wmax, "Nhist_Wmass_mRCmax", 50, 75, 85);
      book(_hist_mRCmin, "hist_mRCmin", 250, 0., 125.);
      book(_hist_mRCmax, "hist_mRCmax", 250, 0., 125.);
      book(_Norm_hist_mRCmin, "Nhist_mRCmin", 250, 0., 125.);
      book(_Norm_hist_mRCmax, "Nhist_mRCmax", 250, 0., 125.);
    }

    /// Perform the per-event analysis
    void analyze(const Event &event)
    {

      // Retrieve dressed leptons, sorted by pT
      FourMomentum pTmiss;
      for (const Particle &p : apply<VisibleFinalState>(event, "vfs").particles())
      {
        pTmiss -= p.momentum();
      }
      const Particles &muonFS = apply<FinalState>(event, "muons").particlesByPt();
      const int nmuon = muonFS.size();

      const Particles &elecFS = apply<FinalState>(event, "elecs").particlesByPt();
      const int nelec = elecFS.size();

      // const FastJets &jetsAntiKt4 = apply<FastJets>(event, "jets");
      // const Jets &jets = jetsAntiKt4.jetsByPt(10.0 * GeV);
      // const int njet = jets.size();

      const ParticlePair &beams = apply<Beam>(event, "Beams").beams();
      const double eC1 = beams.first.E();
      const double eC2 = beams.second.E();
      const Vector3 axis = (beams.first.charge() > 0) ? beams.first.momentum().p3().unit() : beams.second.momentum().p3().unit();

      FourMomentum P_Sum;
      P_Sum.setE(eC1 + eC2);
      P_Sum.setPx(0.);
      P_Sum.setPy(0.);
      P_Sum.setPz(0.);

      pTmiss.setE(eC1 + eC2 + pTmiss.E());
      // const double Emiss = pTmiss.E();

      bool wmassPre = false;
      if (nmuon == 1 && nelec == 1)
      {
        if (muonFS[0].charge() * elecFS[0].charge() < 0)
        {
          wmassPre = true;
        }
      }
      if (wmassPre)
      {
        FourMomentum P_ISR = P_Sum - pTmiss - muonFS[0].momentum() - elecFS[0].momentum();
        FourMomentum P_recoil = P_Sum - muonFS[0].momentum() - elecFS[0].momentum();
        const double mRecoil = P_recoil.mass();
        // MSG_INFO("OSSF Pair Found");
        if (mRecoil > 1.0)
        {
          vector<double> mrc = mRC0(pTmiss, muonFS[0].mom(), elecFS[0].mom(), P_ISR);
          const double mRCmin = mrc[0];
          const double mRCmax = mrc[1];

          _hist_Wmin->fill(mRCmin);
          _hist_Wmax->fill(mRCmax);
          _Norm_hist_Wmin->fill(mRCmin);
          _Norm_hist_Wmax->fill(mRCmax);
          _hist_mRCmin->fill(mRCmin);
          _hist_mRCmax->fill(mRCmax);
          _Norm_hist_mRCmin->fill(mRCmin);
          _Norm_hist_mRCmax->fill(mRCmax);
        }
      }
      else
        return;
    }

    /// Normalise histograms etc., after the run
    void finalize()
    {
      const double sf = crossSection() / femtobarn;
      const double Lint = 5.05e3;
      double norm = sf * Lint;

      MSG_INFO("Norm is " << norm);
      MSG_INFO("Total Cross section is " << crossSection() / femtobarn << " femtobarn!");
      scale(_hist_mRCmin, sf * Lint / sumW()); // norm to generated cross-section in pb (after cuts)
      scale(_hist_mRCmax, sf * Lint / sumW()); // norm to generated cross-section in pb (after cuts)
      scale(_hist_Wmin, sf * Lint / sumW());   // norm to generated cross-section in pb (after cuts)
      scale(_hist_Wmax, sf * Lint / sumW());   // norm to generated cross-section in pb (after cuts)

      normalize(_Norm_hist_Wmin);
      normalize(_Norm_hist_Wmax);
      normalize(_Norm_hist_mRCmin);
      normalize(_Norm_hist_mRCmax);
    }
    vector<double> mRC0(const FourMomentum &met, const FourMomentum &l1, const FourMomentum &l2, const FourMomentum &ISR) const
    {
      // MSG_INFO("Tag 1");
      FourMomentum CM = met + l1 + l2;

      FourMomentum bmet;
      FourMomentum bl1;
      FourMomentum bl2;
      FourMomentum pISR;

      LorentzTransform LT = LorentzTransform::mkFrameTransformFromBeta(CM.betaVec());

      pISR = LT.transform(ISR);
      bmet = LT.transform(met);
      bl1 = LT.transform(l1);
      bl2 = LT.transform(l2);

      const double ss = (bmet.E() + bl1.E() + bl2.E()) / 2.;
      const double pMiss = bmet.p3().mod();
      const double pL1 = bl1.p3().mod();
      const double pL2 = bl2.p3().mod();
      const double EL1 = bl1.E();
      const double EL2 = bl2.E();
      const double EI1 = ss - EL1;
      const double EI2 = ss - EL2;

      if (EI1 > 0. && EI2 > 0. && EI1 + EI2 > pMiss && pL1 + pL2 > pMiss && abs(EI1 - EI2) < pMiss && abs(pL1 - pL2) < pMiss)
      {
        const vector<double> pos_C = solveXY(EI1, EI2, pMiss);
        const vector<double> pos_B = solveXY(pL1, pL2, pMiss);

        const double pMax2 = pow(pos_B[0] - pos_C[0], 2) + pow(pos_B[1] + pos_C[1], 2);
        const double pMin2 = pow(pos_B[0] - pos_C[0], 2) + pow(pos_B[1] - pos_C[1], 2);

        const double mYmax = sqrt(ss * ss - pMin2);
        const double mYmin = sqrt(ss * ss - pMax2);

        // return mYmin;
        const vector<double> mrc = {mYmin, mYmax};
        return mrc;
      }
      else
      {
        const vector<double> mYmax = {-1.0, -1.0};
        return mYmax;
      }
    }
    /// @}
    vector<double> solveXY(const double &p1, const double &p2, const double &pMiss) const
    {
      const double x = 0.5 / pMiss * (pow(p1, 2) - pow(p2, 2) + pow(pMiss, 2));
      const double numinator = 2.0 * (pow(p1, 2) * pow(p2, 2) + pow(p1, 2) * pow(pMiss, 2) + pow(p2, 2) * pow(pMiss, 2)) - pow(p1, 4) - pow(p2, 4) - pow(pMiss, 4);
      const double y = 0.5 / pMiss * sqrt(numinator);
      const vector<double> pos = {x, y};
      return pos;
    }

    /// @name Histograms
    /// @{
    map<string, Histo1DPtr> _h;
    map<string, Profile1DPtr> _p;
    map<string, CounterPtr> _c;
    /// @}
    Histo1DPtr _hist_Wmin;
    Histo1DPtr _hist_Wmax;
    Histo1DPtr _Norm_hist_Wmin;
    Histo1DPtr _Norm_hist_Wmax;
    Histo1DPtr _hist_mRCmin;
    Histo1DPtr _hist_mRCmax;
    Histo1DPtr _Norm_hist_mRCmin;
    Histo1DPtr _Norm_hist_mRCmax;
  };

  RIVET_DECLARE_PLUGIN(CEPC_MASS);

}
