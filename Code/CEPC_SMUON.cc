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
#include <math.h>

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

      // Initialise and register projections

      // The basic final-state projection:
      // all final-state particles within
      // the given eta acceptance
      declare(Beam(), "Beams");

      const FinalState fs(Cuts::abseta < 3.0);
      declare(VisibleFinalState(Cuts::abseta < 3.0), "vfs");

      Cut lepton_cuts = ((Cuts::abseta < 3.0) && (Cuts::E > 0.5 * GeV));

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
      DressedLeptons elec_drs(photons, bare_elec, 0.1, lepton_cuts);
      DressedLeptons muon_drs(photons, bare_muon, 0.1, lepton_cuts);
      // declare(muon_drs, "muons");
      // declare(elec_drs, "elecs");
      declare(bare_elec, "elecs");
      declare(bare_muon, "muons");

      // Missing momentum
      declare(MissingMomentum(fs), "MET");

      // Book histograms
      // specify custom binning
      book(_SR_lowDELTAM, "SR_LOWDELTAM", 5, 0, 5);
      book(_SR_midDELTAM, "SR_MIDDELTAM", 6, 0, 6);
      book(_SR_highDELTAM, "SR_HIGHDELTAM", 4, 0, 4);

      const strings cfnames = {
          "No Cuts",
          "SR-High DeltaM",
          "SR-Mid DeltaM",
          "SR-Low DeltaM"};
      _Cutflows.addCutflow("cf", cfnames);

      _CF_HighDM.addCutflow("HighDM", {"No Cuts",
                                       "Basic Cuts",
                                       "mRCmax > 110",
                                       "Emuon in [45, 75]",
                                       "mRCmin > 80",
                                       "dR(mu, Recoil) > 2.0",
                                       "mll < 120",
                                       "Off_Shell Cuts mll < 60",
                                       "Off_Shell Cuts mRCmin > 95"});
      _CF_MidDM.addCutflow("MidDM", {"No Cuts",
                                     "Basic Cuts",
                                     "mRCmax > 110",
                                     "DeltaM in [20, 80]",
                                     "SR01 DeltaM in [40, 80]",
                                     "SR01 Emu in [40, 60]",
                                     "SR01 mll < 80",
                                     "SR01 dR(mu, Recoil) > 2.5",
                                     "SR01 mRCmin > 85",
                                     "SR02 DeltaM in [35, 70]",
                                     "SR02 Emu in [30, 55]",
                                     "SR02 mRecoil < 140",
                                     "SR02 mRCmin > 85",
                                     "SR03 DeltaM in [30, 60]",
                                     "SR03 Emu in [35, 50]",
                                     "SR03 dR(mu, Recoil) > 2.5",
                                     "SR03 mRCmin > 60",
                                     "SR03 mll < 50",
                                     "SR04 DeltaM in [25, 50]",
                                     "SR04 Emu in [35, 45]",
                                     "SR04 mRCmin > 60",
                                     "SR05 DeltaM in [20, 40]",
                                     "SR05 Emu < 40."});
      _CF_LowDM.addCutflow("LowDM", {"No Cuts",
                                     "Basic Cuts",
                                     "mRCmax > 110",
                                     "DeltaM < 40",
                                     "SR01 DeltaM < 10",
                                     "SR02 DeltaM < 20",
                                     "SR03 DeltaM < 30",
                                     "SR04 DeltaM < 40"});

      book(_hist_Elm, "hist_ELm", 1500, 0., 1500.);
      book(_hist_Elp, "hist_ELp", 1500, 0., 1500.);
      book(_hist_dRlmRec, "hist_dR_lm_Recoil", 100, 0., 20.);
      book(_hist_dRlpRec, "hist_dR_lp_Recoil", 100, 0., 20.);
      book(_hist_mll, "hist_mll", 3000, 0., 3000.);
      book(_hist_mRecoil, "hist_mRecoil", 3000, 0., 3000.);

      book(_hist_mRCmin, "hist_mRCmin", 1500, 0., 1500.);
      book(_hist_mRCmax, "hist_mRCmax", 1500, 0., 1500.);
      book(_hist_mLSPmax, "hist_mLSPmax", 1500, 0., 1500.);
      book(_hist_mRCdM, "hist_mRCdM", 1500, 0., 1500.);
    }

    /// Perform the per-event analysis
    void analyze(const Event &event)
    {
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

      const Particles &muonFS = apply<FinalState>(event, "muons").particlesByPt();
      const int nmuon = muonFS.size();

      const FastJets &jetsAntiKt4 = apply<FastJets>(event, "jets");
      const Jets &jets = jetsAntiKt4.jetsByPt(10.0 * GeV);
      const int njet = jets.size();

      const ParticlePair &beams = apply<Beam>(event, "Beams").beams();
      const double eC = beams.first.E();
      const Vector3 axis = (beams.first.charge() > 0) ? beams.first.momentum().p3().unit() : beams.second.momentum().p3().unit();

      FourMomentum P_Sum;
      P_Sum.setE(2 * eC);
      P_Sum.setPx(0.);
      P_Sum.setPy(0.);
      P_Sum.setPz(0.);

      pTmiss.setE(2 * eC + pTmiss.E());

      bool smuonPre = false;
      if (nmuon == 2 && njet == 0)
      {
        if (muonFS[0].charge() * muonFS[1].charge() < 0)
        {
          smuonPre = true;
        }
      }

      if (smuonPre)
      {
        FourMomentum P_ISR = P_Sum - pTmiss - muonFS[0].momentum() - muonFS[1].momentum();

        const double EmuonM = (muonFS[1].charge() > 0) ? muonFS[0].E() : muonFS[1].E();
        const double EmuonP = (muonFS[1].charge() > 0) ? muonFS[1].E() : muonFS[0].E();

        const double dRlmRec = (muonFS[1].charge() > 0) ? deltaR(muonFS[0].mom(), pTmiss) : deltaR(muonFS[1].mom(), pTmiss);
        const double dRlpRec = (muonFS[1].charge() > 0) ? deltaR(muonFS[1].mom(), pTmiss) : deltaR(muonFS[0].mom(), pTmiss);

        const double mll = (muonFS[0].momentum() + muonFS[1].momentum()).mass();
        const double mRecoil = pTmiss.mass();

        // if (mRecoil > 1.0 && (mll < 85 || mll > 95) && (mRecoil < 85 || mRecoil > 95))
        if (mRecoil > 1.0)
        {
          vector<double> mrc = mRC(pTmiss, muonFS[0].mom(), muonFS[1].mom(), P_ISR);
          const double mRCmin = mrc[0];
          const double mRCmax = mrc[1];
          const double mLSPmax = mrc[2];
          const double DeltaM = mRCmax - mLSPmax;

          _CF_HighDM.fill(2);
          _CF_LowDM.fill(2);
          _CF_MidDM.fill(2);

          _hist_mRCmin->fill(mRCmin);
          _hist_mRCmax->fill(mRCmax);
          _hist_mLSPmax->fill(mLSPmax);
          _hist_mRCdM->fill(DeltaM);
          _hist_mll->fill(mll);
          _hist_Elm->fill(EmuonM);
          _hist_Elp->fill(EmuonP);
          _hist_dRlpRec->fill(dRlpRec);
          _hist_dRlmRec->fill(dRlmRec);
          _hist_mRecoil->fill(mRecoil);

          // if (mRCmax > 110. * GeV)
          if (false)
          {

            _CF_HighDM.fill(3);
            _CF_LowDM.fill(3);
            _CF_MidDM.fill(3);

            // Define the SRs for High delta M, i.e. mLSP ~ 0 GeV
            if (EmuonM > 45. * GeV && EmuonP > 45. * GeV && EmuonM < 75. * GeV && EmuonP < 75. * GeV)
            {

              _CF_HighDM.fill(4);
              if (mRCmin > 85. * GeV)
              {

                _CF_HighDM.fill(5);
                if (dRlpRec > 2.0 && dRlmRec > 2.0)
                {

                  _CF_HighDM.fill(6);
                  if (mll < 120. * GeV)
                  {

                    _CF_HighDM.fill(7);
                    _SR_highDELTAM->fill(1.5);
                    if (mll < 60. * GeV)
                    {
                      _SR_highDELTAM->fill(2.5);
                      _CF_HighDM.fill(8);
                    }
                    if (mRCmin > 95. * GeV)
                    {
                      _SR_highDELTAM->fill(3.5);
                      _CF_HighDM.fill(9);
                    }
                  }
                }
              }
            }
            // Define the SRs for Low delta M, i.e. mY - mLSP ~ 0 GeV
            if (DeltaM < 40. * GeV)
            {

              _CF_LowDM.fill(4);
              if (DeltaM < 10. * GeV)
              {
                _CF_LowDM.fill(5);
                _SR_lowDELTAM->fill(1.5);
              }
              if (DeltaM < 20. * GeV)
              {
                _CF_LowDM.fill(6);
                _SR_lowDELTAM->fill(2.5);
              }
              if (DeltaM < 30. * GeV && EmuonM < 30. * GeV && EmuonM > 20. * GeV && EmuonP < 30. * GeV && EmuonP > 20. * GeV)
              {
                _CF_LowDM.fill(7);
                _SR_lowDELTAM->fill(3.5);
              }
              if (DeltaM < 40. * GeV && EmuonM < 40. * GeV && EmuonM > 25. * GeV && EmuonP < 40. * GeV && EmuonP > 25. * GeV)
              {
                _CF_LowDM.fill(8);
                _SR_lowDELTAM->fill(4.5);
              }
            }

            // Define the SRs for Mid delta M, i.e. mY - mLSP ~ 1/3 mY -> 2/3 mY
            if (DeltaM < 80. * GeV && DeltaM > 20. * GeV)
            {

              _CF_MidDM.fill(4);
              // Define the SR for Mid delta M: mY - mLSP = 1/3 mY -> 5/12 mY
              if (DeltaM < 80. * GeV && DeltaM > 40. * GeV)
              {

                _CF_MidDM.fill(5);
                if (EmuonM > 40. * GeV && EmuonP > 40. * GeV && EmuonM < 60. * GeV && EmuonP < 60. * GeV)
                {
                  _CF_MidDM.fill(6);
                  if (mll < 80. * GeV)
                  {
                    _CF_MidDM.fill(7);
                    if (dRlpRec > 2.5 && dRlmRec > 2.5)
                    {
                      _CF_MidDM.fill(8);

                      if (mRCmin > 85. * GeV)
                      {

                        _CF_MidDM.fill(9);
                        _SR_midDELTAM->fill(1.5);
                      }
                    }
                  }
                }
              }
              // Define the SR for Mid delta M: mY-mLSP = 5/12 mY -> 1/2 mY
              if (DeltaM < 70. * GeV && DeltaM > 35. * GeV)
              {

                _CF_MidDM.fill(10);
                if (EmuonM > 30. * GeV && EmuonM < 55. * GeV && EmuonP > 30. * GeV && EmuonP < 55. * GeV)
                {
                  _CF_MidDM.fill(11);
                  if (mRecoil < 140. * GeV)
                  {

                    _CF_MidDM.fill(12);
                    if (mRCmin > 85 * GeV)
                    {

                      _CF_MidDM.fill(13);
                      _SR_midDELTAM->fill(2.5);
                    }
                  }
                }
              }

              // Define the SR for Mid delta M: mY - mLSP = 1/2 mY -> 7/12 mY
              if (DeltaM < 60. * GeV && DeltaM > 30. * GeV)
              {
                _CF_MidDM.fill(14);
                if (EmuonM > 35. * GeV && EmuonM < 50. * GeV && EmuonP > 35. * GeV && EmuonP < 50. * GeV)
                {
                  _CF_MidDM.fill(15);
                  if (dRlmRec > 2.5 && dRlpRec > 2.5)
                  {
                    _CF_MidDM.fill(16);
                    if (mRCmin > 60. * GeV)
                    {
                      _CF_MidDM.fill(17);
                      if (mll < 50. * GeV)
                      {
                        _CF_MidDM.fill(18);
                        _SR_midDELTAM->fill(3.5);
                      }
                    }
                  }
                }
              }

              // Define the SR for Mid delta M: mY - mLSP = 7/12 mY -> 2/3 mY
              if (DeltaM < 50. * GeV && DeltaM > 25. * GeV)
              {
                _CF_MidDM.fill(19);
                if (EmuonM > 35. * GeV && EmuonM < 45. * GeV && EmuonP > 35. * GeV && EmuonP < 45. * GeV)
                {
                  _CF_MidDM.fill(20);
                  if (mRCmin > 60. * GeV)
                  {
                    _CF_MidDM.fill(21);
                    _SR_midDELTAM->fill(4.5);
                  }
                }
              }

              // Define the SR for Mid delta M: mY - mLSP = 2/3 mY -> 3/4 mY
              if (DeltaM < 40. * GeV && DeltaM > 20. * GeV)
              {
                _CF_MidDM.fill(22);
                if (EmuonM < 40. * GeV && EmuonP < 40. * GeV)
                {
                  _CF_MidDM.fill(23);
                  _SR_midDELTAM->fill(5.5);
                }
              }
            }
          }
        }

        // Define the SRs
        if (dRlmRec < 2.9 && dRlpRec < 2.9 && EmuonM > 40 * GeV && EmuonP > 40 * GeV && mll < 60 * GeV && mRecoil > 40 * GeV)
        {
          _SR_highDELTAM->fill(0.5);
          _Cutflows.fill(2);
        }
        if (dRlmRec > 1.5 && dRlmRec < 2.8 && dRlpRec > 1.5 && dRlpRec < 2.8 && EmuonM > 9 * GeV && EmuonM < 48 * GeV && EmuonP > 9 * GeV && EmuonP < 48 * GeV && mll < 80 * GeV)
        {
          _SR_midDELTAM->fill(0.5);
          _Cutflows.fill(3);
        }
        if (dRlmRec > 1.5 && dRlmRec < 2.8 && dRlpRec > 1.5 && dRlpRec < 2.8 && mRecoil > 220 * GeV)
        {
          _SR_lowDELTAM->fill(0.5);
          _Cutflows.fill(4);
        }
      }
      else
        return;
    }

    /// Normalise histograms etc., after the run
    void finalize()
    {
      // normalize(_h["XXXX"]); // normalize to unity
      // normalize(_h["YYYY"], crossSection()/picobarn); // normalize to generated cross-section in pb (no cuts)

      const double sf = crossSection() / (numEvents() * femtobarn);
      const double XSect = crossSection() / femtobarn;
      const double Luminosity = 5050.0;
      const double norm = Luminosity * XSect;

      scale(_SR_lowDELTAM, XSect / sumW()); // norm to generated cross-section in pb (after cuts)
      scale(_SR_midDELTAM, XSect / sumW());
      scale(_SR_highDELTAM, XSect / sumW());

      scale(_hist_Elm, norm / sumW());
      scale(_hist_Elp, norm / sumW());
      scale(_hist_dRlmRec, norm / sumW());
      scale(_hist_dRlpRec, norm / sumW());
      scale(_hist_mll, norm / sumW());
      scale(_hist_mRecoil, norm / sumW());
      scale(_hist_mRCmin, norm / sumW());
      scale(_hist_mRCmax, norm / sumW());
      scale(_hist_mLSPmax, norm / sumW());
      scale(_hist_mRCdM, norm / sumW());

      _Cutflows.scale(sf * Luminosity);
      // MSG_INFO("CUTFLOWS Smuon + ETmiss Case:\n\n"
      //          << _Cutflows);

      _CF_HighDM.scale(sf * Luminosity);
      // MSG_INFO("CUTFLOWS High DeltaM Case:\n\n"
      //          << _CF_HighDM);

      _CF_LowDM.scale(sf * Luminosity);
      // MSG_INFO("CUTFLOWS Low DeltaM Case:\n\n"
      //          << _CF_LowDM);

      _CF_MidDM.scale(sf * Luminosity);
      MSG_INFO("CUTFLOWS Mid DeltaM Case:\n\n"
               << _CF_MidDM);
    }

    vector<double> mRC(const FourMomentum &met, const FourMomentum &l1, const FourMomentum &l2, const FourMomentum &ISR) const
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

      if (EI1 > 0. && EI2 > 0. && EI1 + EI2 > pMiss && pL1 + pL2 > pMiss)
      {
        const vector<double> pos_C = solveXY(EI1, EI2, pMiss);
        const vector<double> pos_B = solveXY(pL1, pL2, pMiss);

        const double pMax2 = pow(pos_B[0] - pos_C[0], 2) + pow(pos_B[1] + pos_C[1], 2);
        const double pMinX2 = pow(pos_B[0] - pos_C[0], 2);
        const double pMinY2 = pos_B[1] > pos_C[1] ? pow(pos_B[1] - pos_C[1], 2) : 0.0;

        const double mYmax = sqrt(ss * ss - pMinX2 - pMinY2);
        const double mYmin = sqrt(ss * ss - pMax2);
        const double mImax = sqrt(EI1 * EI1 - pos_C[0] * pos_C[0]);

        // return mYmin;
        const vector<double> mrc = {mYmin, mYmax, mImax};
        return mrc;
      }
      else
      {
        const vector<double> mYmax = {-1.0, -1.0, -1.0};
        // const double mYmax = -1.0;
        // MSG_INFO("Something wrong in Reconstructed Mass variables, " << bmet << bl1 << bl2 << pCM);
        return mYmax;
      }
    }

    vector<double> solveXY(const double &p1, const double &p2, const double &pMiss) const
    {
      const double x = 0.5 / pMiss * (pow(p1, 2) - pow(p2, 2) + pow(pMiss, 2));
      const double numinator = 2.0 * (pow(p1, 2) * pow(p2, 2) + pow(p1, 2) * pow(pMiss, 2) + pow(p2, 2) * pow(pMiss, 2)) -pow(p1, 4) - pow(p2, 4) -pow(pMiss, 4);
      const double y = 0.5 / pMiss * sqrt( numinator );
      const vector<double> pos = {x, y};
      return pos;
    }

    ///@}

    /// @name Histograms
    ///@{
    map<string, Histo1DPtr> _h;
    map<string, Profile1DPtr> _p;
    map<string, CounterPtr> _c;

    Histo1DPtr _hist_Elm;
    Histo1DPtr _hist_Elp;
    Histo1DPtr _hist_dRlmRec;
    Histo1DPtr _hist_dRlpRec;
    Histo1DPtr _hist_mll;
    Histo1DPtr _hist_mRecoil;

    Histo1DPtr _hist_mRCmin;
    Histo1DPtr _hist_mRCmax;
    Histo1DPtr _hist_mLSPmax;
    Histo1DPtr _hist_mRCdM;

    ///@}
  private:
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
