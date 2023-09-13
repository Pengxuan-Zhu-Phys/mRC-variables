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
  class CEPC_SMUON : public Analysis
  {
  public:
    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CEPC_SMUON);

    /// @name Analysis methods
    ///@{
    /// Book histograms and initialise projections before the run
    void init()
    {

      std::ofstream df;
      // dfname = getOption("DFNAME");
      // dfname.append(".csv");

      // df.open(dfname, std::ios::out | std::ios::app);
      // MSG_INFO("MSG_INFO Open file: " << dfname);
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
      // take binning from reference data using HEPData ID (digits in "d01-x01-y01" etc.)
      book(_SR_lowDELTAM, "SR_LOWDELTAM", 6, 0, 6);
      book(_SR_midDELTAM, "SR_MIDDELTAM", 6, 0, 6);
      book(_SR_highDELTAM, "SR_HIGHDELTAM", 3, 0, 3);

      // Book Cut-flows
      const strings cfnames = {
          "No Cuts",
          "SR-High DeltaM",
          "SR-Mid DeltaM",
          "SR-Low DeltaM"};
      _Cutflows.addCutflow("cf", cfnames);

      _CF_HighDM.addCutflow("HighDM", {"No Cuts",
                                       "Basic Cuts",
                                       "mRCmax > 110",
                                       "mRCmin > 85",
                                       "High Mass:  mRCmax > 117",
                                       "High Mass:  mRCmin > 95",
                                       "High Mass:  Emuon < 70."});
      _CF_MidDM.addCutflow("MidDM", {"No Cuts",
                                     "Basic Cuts",
                                     "mRCmax > 110",
                                     "SR01: mLSPmax in [40, 60]",
                                     "SR01: mRCmin > 85 && mRCmax(40) > 110",
                                     "SR02: mLSPmax in [50, 70]",
                                     "SR02: mRCmax(40) > 110",
                                     "SR02: mRCmin(40) > 100",
                                     "SR03: mLSPmax in [60, 80]",
                                     "SR03: mRCmin(40) > 95",
                                     "SR04: mLSPmax in [70, 85]",
                                     "SR04: mRCmin(70) > 100.",
                                     "SR05: mLSPmax in [80, 95]",
                                     "SR05: mRCmin(80) > 105"});
      _CF_LowDM.addCutflow("LowDM", {"No Cuts",
                                     "Basic Cuts",
                                     "mRCmax > 110",
                                     "DeltaM in [35, 50]",
                                     "DeltaM in [25, 40]",
                                     "DeltaM in [15, 30]",
                                     "DeltaM < 20",
                                     "DeltaM < 10"});

      book(_hist_Elm, "hist_ELm", 250, 0., 125.);
      book(_hist_Elp, "hist_ELp", 250, 0., 125.);
      book(_hist_dRlmRec, "hist_dR_lm_Recoil", 100, 0., 5.);
      book(_hist_dRlpRec, "hist_dR_lp_Recoil", 100, 0., 5.);
      book(_hist_dRll, "hist_dR_ll", 100, 0., 5.);
      book(_hist_mll, "hist_mll", 250, 0., 250.);
      book(_hist_mRecoil, "hist_mRecoil", 250, 0., 250.);
      book(_hist_mMiss, "hist_mMiss", 250, 0., 250.);

      book(_hist_mRCmin, "hist_mRCmin", 250, 0., 125.);
      book(_hist_mRCmax, "hist_mRCmax", 250, 0., 125.);
      book(_hist_mRCLSP, "hist_mRCLSP", 250, 0., 125.);
      book(_hist_mLSPmax, "hist_mLSPmax", 250, 0., 125.);
      book(_hist_mDM_RC, "hist_dM_RC", 250, 0., 125.);
      book(_hist_mDM_LSP, "hist_dM_LSP", 250, 0., 125.);

      book(_hist_Wmin, "hist_Wmass_mRCmin", 4000, 78, 82);
      book(_hist_Wmax, "hist_Wmass_mRCmax", 4000, 78, 82);
    }

    /// Perform the per-event analysis
    void analyze(const Event &event)
    {
            // MSG_INFO("MSG_INFO Open file: " << dfname );

      // const double weight = 1.0;
      // const MissingMomentum met = apply<MissingMomentum>(event, "MET");
      // double ETmiss = met.missingEt();
      // const double weight = 1.0;
      _Cutflows.fillinit();
      _Cutflows.fill(1);

      _CF_HighDM.fillinit();
      _CF_HighDM.fill(1);

      _CF_LowDM.fillinit();
      _CF_LowDM.fill(1);

      _CF_MidDM.fillinit();
      _CF_MidDM.fill(1);

      FourMomentum pTmiss;
      for (const Particle &p : apply<VisibleFinalState>(event, "vfs").particles())
      {
        pTmiss -= p.momentum();
      }
      // double ETmiss = pTmiss.pT();
      // Retrieve dressed leptons, sorted by pT
      // vector<DressedLepton> leptons = apply<DressedLeptons>(event, "leptons").dressedLeptons();

      const Particles &muonFS = apply<FinalState>(event, "muons").particlesByPt();
      const int nmuon = muonFS.size();

      const Particles &elecFS = apply<FinalState>(event, "elecs").particlesByPt();
      const int nelec = elecFS.size();

      const FastJets &jetsAntiKt4 = apply<FastJets>(event, "jets");
      const Jets &jets = jetsAntiKt4.jetsByPt(10.0 * GeV);
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

      bool smuonPre = false;
      // if (nmuon == 2 && njet == 0)
      // {
      //   if (muonFS[0].charge() * muonFS[1].charge() < 0 && Emiss > 1. * GeV)
      //   {
      //     smuonPre = true;
      //   }
      // }
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
        FourMomentum P_ISR = P_Sum - pTmiss - muonFS[0].momentum() - muonFS[1].momentum();
        FourMomentum P_recoil = P_Sum - muonFS[0].momentum() - muonFS[1].momentum();
        const double mRecoil = P_recoil.mass();
        if (mRecoil > 1.0) {
          vector<double> mrc = mRC0(pTmiss, muonFS[0].mom(), muonFS[1].mom(), P_ISR);
          const double mRCmin = mrc[0];
          const double mRCmax = mrc[1];

          _hist_Wmin -> fill(mRCmin);
          _hist_Wmax -> fill(mRCmax);
          _hist_mRCmin->fill(mRCmin);
          _hist_mRCmax->fill(mRCmax);
        }
      }
      else 
        return;

      if (smuonPre)
      {
        FourMomentum P_ISR = P_Sum - pTmiss - muonFS[0].momentum() - muonFS[1].momentum();
        FourMomentum P_recoil = P_Sum - muonFS[0].momentum() - muonFS[1].momentum();

        const double EmuonM = (muonFS[1].charge() > 0) ? muonFS[0].E() : muonFS[1].E();
        const double EmuonP = (muonFS[1].charge() > 0) ? muonFS[1].E() : muonFS[0].E();

        const double dRlmRec = (muonFS[1].charge() > 0) ? deltaR(muonFS[0].mom(), pTmiss) : deltaR(muonFS[1].mom(), pTmiss);
        const double dRlpRec = (muonFS[1].charge() > 0) ? deltaR(muonFS[1].mom(), pTmiss) : deltaR(muonFS[0].mom(), pTmiss);

        const double dRll = deltaR(muonFS[0].mom(), muonFS[1].mom());

        const double mll = (muonFS[0].momentum() + muonFS[1].momentum()).mass();
        const double mRecoil = P_recoil.mass();
        const double mMiss = pTmiss.mass();

        if (EmuonM > 40. && EmuonP > 40. && dRlmRec < 2.9 && dRlpRec < 2.9 && mll < 60. && mRecoil > 40.)
        {
          _SR_highDELTAM->fill(2.5);
        }
        if (EmuonM > 9. && EmuonP > 9. && EmuonM < 48. && EmuonP < 48. && dRlmRec < 2.8 && dRlpRec < 2.8 && dRlmRec > 1.5 && dRlpRec > 1.5 && mll < 80.)
        {
          _SR_midDELTAM->fill(5.5);
        }
        if (dRlmRec < 2.8 && dRlpRec < 2.8 && dRlmRec > 1.5 && dRlpRec > 1.5 && mRecoil > 220.)
        {
          _SR_lowDELTAM->fill(5.5);
        }

        if (mRecoil > 1.0 && (mll < 80 || mll > 100) && (mRecoil < 85 || mRecoil > 95))
        // if (mRecoil > 1.0)
        {
          vector<double> mrc = mRC(pTmiss, muonFS[0].mom(), muonFS[1].mom(), P_ISR);
          const double mRCmin = mrc[0];
          const double mRCmax = mrc[1];
          const double mLSPmax = mrc[2];
          const double mRCLSP = mrc[3];
          const double dm_RC = mRCmax - mLSPmax;
          const double dm_LSP = mRCLSP - mLSPmax;

          _CF_HighDM.fill(2);
          _CF_LowDM.fill(2);
          _CF_MidDM.fill(2);

          _hist_mRCmin->fill(mRCmin);
          _hist_mRCmax->fill(mRCmax);
          _hist_mLSPmax->fill(mLSPmax);
          _hist_mRCLSP->fill(mRCLSP);
          _hist_mDM_RC->fill(dm_RC);
          _hist_mDM_LSP->fill(dm_LSP);
          _hist_mll->fill(mll);
          _hist_Elm->fill(EmuonM);
          _hist_Elp->fill(EmuonP);
          _hist_dRlpRec->fill(dRlpRec);
          _hist_dRlmRec->fill(dRlmRec);
          _hist_dRll->fill(dRll);
          _hist_mRecoil->fill(mRecoil);
          _hist_mMiss->fill(mMiss);

          if (mRCmax > 110)
          {
            _CF_HighDM.fill(3);
            _CF_LowDM.fill(3);
            _CF_MidDM.fill(3);
            // Define SR-HighDeltaM For DeltaM in [80, 120]
            if (mRCmin > 85)
            {
              _CF_HighDM.fill(4);
              _SR_highDELTAM->fill(0.5);
            }
            if (mRCmax > 117)
            {
              _CF_HighDM.fill(5);
              if (mRCmin > 95)
              {
                _CF_HighDM.fill(6);
                if (EmuonM < 70. && EmuonP < 70. && EmuonM > 50. && EmuonP > 50.)
                {
                  _SR_highDELTAM->fill(1.5);
                  _CF_HighDM.fill(7);
                }
              }
            }
            // Define SR-LowDeltaM For DeltaM in [0, 50]
            if (dm_RC < 50. && dm_RC > 35. && EmuonP < 44. && EmuonP > 34. && EmuonM < 44. && EmuonM > 34.)
            {
              _CF_LowDM.fill(4);
              _SR_lowDELTAM->fill(0.5);
            }
            if (dm_RC < 40. && dm_RC > 25. && EmuonP < 37. && EmuonP > 28. && EmuonM < 37. && EmuonM > 28.)
            {
              _CF_LowDM.fill(5);
              _SR_lowDELTAM->fill(1.5);
            }
            if (dm_RC < 30. && dm_RC > 15. && EmuonP < 28. && EmuonP > 22. && EmuonM < 28. && EmuonM > 22.)
            {
              _SR_lowDELTAM->fill(2.5);
              _CF_LowDM.fill(6);
            }
            if (dm_RC < 20. && EmuonP < 18. && EmuonP > 15. && EmuonM < 18. && EmuonM > 15.)
            {
              _SR_lowDELTAM->fill(3.5);
              _CF_LowDM.fill(7);
            }
            if (dm_RC < 10.)
            {
              _SR_lowDELTAM->fill(4.5);
              _CF_LowDM.fill(8);
            }

            // Define SR-MedDELTAM For DeltaM in [40, 80]
            if (mLSPmax > 40. && mLSPmax < 60.)
            {
              _CF_MidDM.fill(4);
              vector<double> mrc40 = mRC(pTmiss, muonFS[0].mom(), muonFS[1].mom(), P_ISR, 40.);
              if (mrc40[1] > 110 && mRCmin > 85.)
              {
                _CF_MidDM.fill(5);
                _SR_midDELTAM->fill(0.5);
              }
            }
            if (mLSPmax > 50. && mLSPmax < 70.)
            {
              _CF_MidDM.fill(6);
              vector<double> mrc40 = mRC(pTmiss, muonFS[0].mom(), muonFS[1].mom(), P_ISR, 40.);
              if (mrc40[1] > 110.)
              {
                _CF_MidDM.fill(7);
                if (mrc40[0] > 100.)
                {
                  _CF_MidDM.fill(8);
                  _SR_midDELTAM->fill(1.5);
                }
              }
            }
            if (mLSPmax > 60. && mLSPmax < 80.)
            {
              _CF_MidDM.fill(9);
              vector<double> mrc40 = mRC(pTmiss, muonFS[0].mom(), muonFS[1].mom(), P_ISR, 40.);
              if (mrc40[0] > 95.)
              {
                _CF_MidDM.fill(10);
                _SR_midDELTAM->fill(2.5);
              }
            }
            if (mLSPmax > 70. && mLSPmax < 85.)
            {
              _CF_MidDM.fill(11);
              vector<double> mrc70 = mRC(pTmiss, muonFS[0].mom(), muonFS[1].mom(), P_ISR, 70.);
              if (mrc70[0] > 100.)
              {
                _CF_MidDM.fill(12);
                _SR_midDELTAM->fill(3.5);
              }
            }
            if (mLSPmax > 80. && mLSPmax < 95)
            {
              _CF_MidDM.fill(13);
              vector<double> mrc80 = mRC(pTmiss, muonFS[0].mom(), muonFS[1].mom(), P_ISR, 80.);
              if (mrc80[0] > 105.)
              {
                _CF_MidDM.fill(14);
                _SR_midDELTAM->fill(4.5);
              }

              // _hist_mRCmin->fill(mrc80[0]);
              // _hist_mRCmax->fill(mrc80[1]);
              // _hist_mLSPmax->fill(mrc80[2]);
              // _hist_mRCLSP->fill(mrc80[3]);
              // _hist_mDM_RC->fill(mrc80[1] - mrc80[2]);
              // _hist_mDM_LSP->fill(mrc80[3] - mrc80[2]);
              // _hist_mll->fill(mll);
              // _hist_Elm->fill(EmuonM);
              // _hist_Elp->fill(EmuonP);
              // _hist_dRlpRec->fill(dRlpRec);
              // _hist_dRlmRec->fill(dRlmRec);
              // _hist_dRll->fill(dRll);
              // _hist_mRecoil->fill(mRecoil);
              // _hist_mMiss->fill(mMiss);
            }
          }
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
      // normalize(_h["XXXX"]);                                 // normalize to unity
      // normalize(_h["YYYY"], crossSection() / picobarn);      // normalize to generated cross-section in pb (no cuts)
      scale(_hist_mRCmin, crossSection() / femtobarn / sumW());  // norm to generated cross-section in pb (after cuts)
      scale(_hist_mRCmax, crossSection() / femtobarn / sumW());  // norm to generated cross-section in pb (after cuts)

      scale(_hist_Wmin, crossSection() / femtobarn / sumW());
      scale(_hist_Wmax, crossSection() / femtobarn / sumW());

      // scale(_hist_mLSPmax, crossSection() / femtobarn / sumW()); // norm to generated cross-section in pb (after cuts)
      // scale(_hist_mRCLSP, crossSection() / femtobarn / sumW());  // norm to generated cross-section in pb (after cuts)
      // scale(_hist_mDM_RC, crossSection() / femtobarn / sumW());  // norm to generated cross-section in pb (after cuts)
      // scale(_hist_mDM_LSP, crossSection() / femtobarn / sumW()); // norm to generated cross-section in pb (after cuts)
      // scale(_hist_mll, crossSection() / femtobarn / sumW());     // norm to generated cross-section in pb (after cuts)
      // scale(_hist_Elm, crossSection() / femtobarn / sumW());     // norm to generated cross-section in pb (after cuts)
      // scale(_hist_Elp, crossSection() / femtobarn / sumW());     // norm to generated cross-section in pb (after cuts)
      // scale(_hist_dRlpRec, crossSection() / femtobarn / sumW()); // norm to generated cross-section in pb (after cuts)
      // scale(_hist_dRlmRec, crossSection() / femtobarn / sumW()); // norm to generated cross-section in pb (after cuts)
      // scale(_hist_dRll, crossSection() / femtobarn / sumW());    // norm to generated cross-section in pb (after cuts)
      // scale(_hist_mRecoil, crossSection() / femtobarn / sumW()); // norm to generated cross-section in pb (after cuts)
      // scale(_hist_mMiss, crossSection() / femtobarn / sumW());   // norm to generated cross-section in pb (after cuts)

      // scale(_SR_highDELTAM, sf * Lint / sumW());
      // scale(_SR_midDELTAM, sf * Lint / sumW());
      // scale(_SR_lowDELTAM, sf * Lint / sumW());

      // _Cutflows.scale(sf * Lint / numEvents());
      // MSG_INFO("CUTFLOWS Smuon + ETmiss Case:\n\n"
      //          << _Cutflows);

      // _CF_HighDM.scale(sf * Lint / numEvents());
      // MSG_INFO("CUTFLOWS High DeltaM Case:\n\n"
      //  << _CF_HighDM);

      _CF_LowDM.scale(sf * Lint / numEvents());
      MSG_INFO("CUTFLOWS Low DeltaM Case:\n\n"
               << _CF_LowDM);

      _CF_MidDM.scale(sf * Lint / numEvents());
      // MSG_INFO("CUTFLOWS Low DeltaM Case:\n\n"
      //          << _CF_MidDM);
      df.close();
      MSG_INFO("Close file" << dfname);
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

    vector<double> mRC(const FourMomentum &met, const FourMomentum &l1, const FourMomentum &l2, const FourMomentum &ISR, const double &mI = 0.0) const
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

      const double pI1max = sqrt(EI1 * EI1 - mI * mI);
      const double pI2max = sqrt(EI2 * EI2 - mI * mI);

      if (pI1max > 0. && pI2max > 0. && pI1max + pI2max > pMiss && pL1 + pL2 > pMiss && abs(pI1max - pI2max) < pMiss && abs(pL1 - pL2) < pMiss)
      {
        const vector<double> pos_C = solveXY(pI1max, pI2max, pMiss);
        const vector<double> pos_B = solveXY(pL1, pL2, pMiss);

        const double pMax2 = pow(pos_B[0] - pos_C[0], 2) + pow(pos_B[1] + pos_C[1], 2);
        const double pMinX2 = pow(pos_B[0] - pos_C[0], 2);
        const double pMinY2 = pos_B[1] > pos_C[1] ? pow(pos_B[1] - pos_C[1], 2) : 0.0;
        const double pLSP2 = pow(pos_B[0] - pos_C[0], 2) + pow(pos_B[1], 2);

        const double mYmax = sqrt(ss * ss - pMinX2 - pMinY2);
        const double mYmin = sqrt(ss * ss - pMax2);
        const double mImax = sqrt(EI1 * EI1 - pos_C[0] * pos_C[0]);
        const double mYLSP = sqrt(ss * ss - pLSP2);

        // return mYmin;
        const vector<double> mrc = {mYmin, mYmax, mImax, mYLSP};
        return mrc;
      }
      else
      {
        const vector<double> mYmax = {-1.0, -1.0, -1.0, -1.0};
        // const double mYmax = -1.0;
        // MSG_INFO("Something wrong in Reconstructed Mass variables, " << bmet << bl1 << bl2 << pCM);
        return mYmax;
      }
    }

    vector<double> solveXY(const double &p1, const double &p2, const double &pMiss) const
    {
      const double x = 0.5 / pMiss * (pow(p1, 2) - pow(p2, 2) + pow(pMiss, 2));
      const double numinator = 2.0 * (pow(p1, 2) * pow(p2, 2) + pow(p1, 2) * pow(pMiss, 2) + pow(p2, 2) * pow(pMiss, 2)) - pow(p1, 4) - pow(p2, 4) - pow(pMiss, 4);
      const double y = 0.5 / pMiss * sqrt(numinator);
      const vector<double> pos = {x, y};
      return pos;
    }

    ///@}

    /// @name Histograms

    double nPass = 0.;
    ///@}

  private:
    std::ofstream df;
    std::string dfname;

    Cutflow _flow;
    Cutflows _Cutflows_em;
    Cutflows _Cutflows_mm;
    Cutflows _Cutflows_ee;

    // CounterPtr _srcounts[4];

    Histo1DPtr _hist_Elm;
    Histo1DPtr _hist_Elp;
    Histo1DPtr _hist_dRlmRec;
    Histo1DPtr _hist_dRlpRec;
    Histo1DPtr _hist_dRll;
    Histo1DPtr _hist_mll;
    Histo1DPtr _hist_mRecoil;
    Histo1DPtr _hist_mMiss;

    Histo1DPtr _hist_mRCmin;
    Histo1DPtr _hist_mRCmax;
    Histo1DPtr _hist_mLSPmax;
    Histo1DPtr _hist_mRCLSP;
    Histo1DPtr _hist_mDM_RC;
    Histo1DPtr _hist_mDM_LSP;

    Histo1DPtr _hist_Wmin;
    Histo1DPtr _hist_Wmax;

    ///@}
    Cutflows _Cutflows;
    Cutflows _CF_HighDM;
    Cutflows _CF_MidDM;
    Cutflows _CF_LowDM;

    Histo1DPtr _SR_highDELTAM;
    Histo1DPtr _SR_midDELTAM;
    Histo1DPtr _SR_lowDELTAM;
  };

  RIVET_DECLARE_PLUGIN(CEPC_SMUON);
}