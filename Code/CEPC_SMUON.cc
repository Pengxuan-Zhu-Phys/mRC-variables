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

namespace Rivet {


  /// @brief Add a short analysis description here
  class CEPC_SMUON : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CEPC_SMUON);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

// Initialise and register projections

      // The basic final-state projection:
      // all final-state particles within
      // the given eta acceptance
      declare(Beam(), "Beams");

      const FinalState fs(Cuts::abseta < 3.0);
      declare(VisibleFinalState(Cuts::abseta < 3.0), "vfs");

      Cut lepton_cuts = ((Cuts::abseta < 3.0) && (Cuts::pT > 5 * GeV));

      // FinalState bare_elec(Cuts::abspid == PID::ELECTRON);
      // FinalState bare_muon(Cuts::abspid == PID::MUON);
      PromptFinalState bare_elec(Cuts::abspid == PID::ELECTRON && lepton_cuts);
      PromptFinalState bare_muon(Cuts::abspid == PID::MUON && lepton_cuts);

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
      // declare(muon_drs, "muons");
      // declare(elec_drs, "elecs");
      declare(bare_elec, "elecs");
      declare(bare_muon, "muons");

      // Missing momentum
      declare(MissingMomentum(fs), "MET");

      // Book histograms
      // specify custom binning

      // take binning from reference data using HEPData ID (digits in "d01-x01-y01" etc.)

      book(_hist_mT2, "hist_mT2", 25, 0., 125.);
      book(_hist_mll, "hist_mll", 40, 0., 200.);
      book(_hist_ETmiss, "hist_ETmiss", 20, 0.0, 100.);
      book(_hist_cThetalp, "hist_CosThetalp", 20, -1., 1.);
      book(_hist_cThetalm, "hist_CosThetalm", 20, -1., 1.);

      book(_NH_mT2, "normlizedhist_mT2", 25, 0., 125.);
      book(_NH_mll, "normlizedhist_mll", 40, 0., 200.);
      book(_NH_ETmiss, "normlizedhist_ETmiss", 20, 0.0, 100.);
      book(_NH_cThetalp, "normlizedhist_CosThetalp", 20, -1., 1.);
      book(_NH_cThetalm, "normlizedhist_CosThetalm", 20, -1., 1.);

      // Book Cut-flows
      const strings cfnames = {
          "No-Cut",
          "Selection",
          "MT2",
          "cosTheta"};

      _Cutflows_ee.addCutflow("ee", cfnames);
      _Cutflows_em.addCutflow("em", cfnames);
      _Cutflows_mm.addCutflow("mm", cfnames);

      // Book SR counters
      // for (size_t i = 0; i < 3; ++i)
      // {
      //   book(_srcounts[i], "SR_" + toString(i));
      // }
    }

    /// Perform the per-event analysis
    void analyze(const Event &event)
    {

      _Cutflows_ee.fillinit();
      _Cutflows_em.fillinit();
      _Cutflows_mm.fillinit();

      _Cutflows_ee.fill(1);
      _Cutflows_em.fill(1);
      _Cutflows_mm.fill(1);

      const double weight = 1.0;
      const MissingMomentum met = apply<MissingMomentum>(event, "MET");
      double ETmiss = met.missingEt();

      FourMomentum pTmiss;
      for (const Particle &p : apply<VisibleFinalState>(event, "vfs").particles())
      {
        pTmiss -= p.momentum();
      }
      // double ETmiss = pTmiss.pT();
      // Retrieve dressed leptons, sorted by pT
      // vector<DressedLepton> leptons = apply<DressedLeptons>(event, "leptons").dressedLeptons();
      const Particles &elecFS = apply<FinalState>(event, "elecs").particlesByPt();
      const Particles &muonFS = apply<FinalState>(event, "muons").particlesByPt();
      const Particles &leptFS = elecFS + muonFS;
      sortByPt(leptFS);

      bool eeMET = false;
      bool emMET = false;
      bool mmMET = false;

      const int nelec = elecFS.size();
      const int nmuon = muonFS.size();

      const ParticlePair &beams = apply<Beam>(event, "Beams").beams();
      Vector3 axis;
      if (beams.first.charge() > 0)
      {
        axis = beams.first.momentum().p3().unit();
      }
      else
      {
        axis = beams.second.momentum().p3().unit();
      }

      if (nelec == 2 && nmuon == 0)
      {
        if (elecFS[0].charge() * elecFS[1].charge() < 0)
        {
          eeMET = true;
        }
      }
      if (nelec == 0 && nmuon == 2)
      {
        if (muonFS[0].charge() * muonFS[1].charge() < 0)
        {
          mmMET = true;
        }
      }
      if (nelec == 1 && nmuon == 1)
      {
        if (elecFS[0].charge() * muonFS[0].charge() < 0)
        {
          emMET = true;
        }
      }

      if (ETmiss < 5. * GeV)
        vetoEvent;

      if (eeMET)
      {
        double mll = (elecFS[0].momentum() + elecFS[1].momentum()).mass();
        if ((mll < 80. * GeV) || (mll > 100. * GeV))
        {
          _Cutflows_ee.fill(2);
          const double m_T2 = mT2(elecFS[0].momentum(), elecFS[1].momentum(), pTmiss, 0.0);
          // _hist_mT2->fill(m_T2, weight);

          if (m_T2 > 50. * GeV)
          {
            _Cutflows_ee.fill(3);
          }

          double ctheta_lp;
          double ctheta_lm;
          if (elecFS[0].charge() > 0)
          {
            ctheta_lp = elecFS[0].momentum().p3().unit().dot(axis);
            ctheta_lm = elecFS[1].momentum().p3().unit().dot(axis);
          }
          else
          {
            ctheta_lp = elecFS[1].momentum().p3().unit().dot(axis);
            ctheta_lm = elecFS[0].momentum().p3().unit().dot(axis);
          }
          // _hist_cThetalp->fill(ctheta_lp);
          // _hist_cThetalm->fill(ctheta_lm);

          if (ctheta_lp < 0.3 && ctheta_lm > -0.3)
          {
            _Cutflows_ee.fill(4);
          }
        }
      }
      else if (mmMET)
      {
        nPass += 1.;

        double mll = (muonFS[0].momentum() + muonFS[1].momentum()).mass();
        if ((mll < 80. * GeV) || (mll > 100. * GeV))
        {
          // {

          const double m_T2 = mT2(muonFS[0].momentum(), muonFS[1].momentum(), pTmiss, 0.0);
          _hist_mT2->fill(m_T2, weight);
          _NH_mT2->fill(m_T2, weight);

          // else
          //   vetoEvent;

          double ctheta_lp;
          double ctheta_lm;
          if (muonFS[0].charge() > 0)
          {
            ctheta_lp = muonFS[0].momentum().p3().unit().dot(axis);
            ctheta_lm = muonFS[1].momentum().p3().unit().dot(axis);
          }
          else
          {
            ctheta_lp = muonFS[1].momentum().p3().unit().dot(axis);
            ctheta_lm = muonFS[0].momentum().p3().unit().dot(axis);
          }
          _hist_cThetalp->fill(ctheta_lp);
          _hist_cThetalm->fill(ctheta_lm);
          _NH_cThetalp->fill(ctheta_lp);
          _NH_cThetalm->fill(ctheta_lm);
          _Cutflows_mm.fill(2);
          if (m_T2 > 20. * GeV)
          {
          _hist_ETmiss->fill(ETmiss, weight);
          _NH_ETmiss->fill(ETmiss, weight);
          _hist_mll->fill(mll, weight);
          _NH_mll->fill(mll, weight);
            _Cutflows_mm.fill(3);
            if (ctheta_lp < 0.3 && ctheta_lm > -0.3)
            {
              _Cutflows_mm.fill(4);
            }
          }
        }
        // }
      }
      else if (emMET)
      {
        double mll = (elecFS[0].momentum() + muonFS[0].momentum()).mass();
        // _hist_mll->fill(mll, weight);
        _Cutflows_em.fill(2);
        const double m_T2 = mT2(elecFS[0].momentum(), muonFS[0].momentum(), pTmiss, 0.0);
        // _hist_mT2->fill(m_T2, weight);

        if (m_T2 > 20. * GeV)
        {
          _Cutflows_em.fill(3);
        }

        double ctheta_lp;
        double ctheta_lm;
        if (elecFS[0].charge() > 0 && muonFS[0].charge() < 0)
        {
          ctheta_lp = elecFS[0].momentum().p3().unit().dot(axis);
          ctheta_lm = muonFS[0].momentum().p3().unit().dot(axis);
        }
        else if (elecFS[0].charge() < 0 && muonFS[0].charge() > 0)
        {
          ctheta_lp = muonFS[0].momentum().p3().unit().dot(axis);
          ctheta_lm = elecFS[0].momentum().p3().unit().dot(axis);
        }
        else
          vetoEvent;

        // _hist_cThetalp->fill(ctheta_lp);
        // _hist_cThetalm->fill(ctheta_lm);

        if (ctheta_lp < 0.3 && ctheta_lm > -0.3)
        {
          _Cutflows_em.fill(4);
        }
      }
      else
        vetoEvent;
    }

    /// Normalise histograms etc., after the run
    void finalize()
    {

      const double sf = crossSection() / (numEvents() * femtobarn);
      const double Lint = 5.e3;
      double norm = sf * Lint;

      MSG_INFO("Norm is " << norm);

      normalize(_hist_mT2, norm);
      normalize(_hist_mll, norm);
      normalize(_hist_ETmiss, norm);
      normalize(_hist_cThetalp, norm);
      normalize(_hist_cThetalm, norm);

      normalize(_NH_mT2);
      normalize(_NH_mll);
      normalize(_NH_ETmiss);
      normalize(_NH_cThetalp);
      normalize(_NH_cThetalm);

      MSG_INFO("Total Cross section is " << crossSection() / femtobarn << " femtobarn!");
      _Cutflows_ee.scale(sf);
      _Cutflows_em.scale(sf);
      _Cutflows_mm.scale(sf);
      // MSG_INFO("CUTFLOWS ee+ETmiss Case:\n\n"
      //          << _Cutflows_ee);
      // MSG_INFO("CUTFLOWS em+ETmiss Case:\n\n"
      //          << _Cutflows_em);
      MSG_INFO("CUTFLOWS mm+ETmiss Case:\n\n"
               << _Cutflows_mm);

      // normalize(_h["XXXX"]);                                 // normalize to unity
      // normalize(_h["YYYY"], crossSection() / picobarn);      // normalize to generated cross-section in pb (no cuts)
      // scale(_h["ZZZZ"], crossSection() / picobarn / sumW()); // norm to generated cross-section in pb (after cuts)
    }

    ///@}

    /// @name Histograms
    ///@{
    map<string, Histo1DPtr> _h;
    map<string, Profile1DPtr> _p;
    map<string, CounterPtr> _c;

    double nPass = 0.;
    ///@}

  private:
    Cutflow _flow;
    Cutflows _Cutflows_em;
    Cutflows _Cutflows_mm;
    Cutflows _Cutflows_ee;

    Histo1DPtr _NH_mT2;
    Histo1DPtr _NH_mll;
    Histo1DPtr _NH_ETmiss;
    Histo1DPtr _NH_cThetalp;
    Histo1DPtr _NH_cThetalm;

    Histo1DPtr _hist_mT2;
    Histo1DPtr _hist_mll;
    Histo1DPtr _hist_ETmiss;
    Histo1DPtr _hist_cThetalp;
    Histo1DPtr _hist_cThetalm;

    // CounterPtr _srcounts[4];
  };

  RIVET_DECLARE_PLUGIN(CEPC_SMUON);

}
