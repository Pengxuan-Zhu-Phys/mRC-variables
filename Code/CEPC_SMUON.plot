BEGIN PLOT /CEPC_SLEPTON/hist_CosThetalm
XLabel=$\cos{\theta_{\ell^-}}$ 
LogY=0
YLabel=Events / 0.05
# + any additional plot settings you might like, see make-plots documentation
END PLOT

BEGIN PLOT /CEPC_SLEPTON/hist_CosThetalp
XLabel=$\cos{\theta_{\ell^+}}$ 
LogY=0
YLabel=Events / 0.05
# + any additional plot settings you might like, see make-plots documentation
END PLOT

BEGIN PLOT /CEPC_SLEPTON/hist_ETmiss
XLabel=$E_{\rm T}^{\rm miss}$ [GeV]
LogY=0
YLabel=Events / 5 GeV
# + any additional plot settings you might like, see make-plots documentation
END PLOT

BEGIN PLOT /CEPC_SLEPTON/hist_mll
XLabel=$m_{\ell \ell}$ [GeV]
LogY=0
YLabel=Events / 5 GeV
# + any additional plot settings you might like, see make-plots documentation
END PLOT

BEGIN PLOT /CEPC_SLEPTON/hist_mT2
XLabel=$m_{T2}$ [GeV]
LogY=0
YLabel=Events / 5 GeV
# + any additional plot settings you might like, see make-plots documentation
END PLOT

BEGIN PLOT /CEPC_SLEPTON/normlizedhist_CosThetalm
XLabel=$\cos{\theta_{\ell^-}}$
YLabel=Fraction of Events / 0.05
LogY=1
# + any additional plot settings you might like, see make-plots documentation
END PLOT

BEGIN PLOT /CEPC_SLEPTON/normlizedhist_CosThetalp
XLabel=$\cos{\theta_{\ell^+}}$ 
YLabel=Fraction of Events / 0.05
# + any additional plot settings you might like, see make-plots documentation
END PLOT

BEGIN PLOT /CEPC_SLEPTON/normlizedhist_ETmiss
XLabel=$E_{\rm T}^{\rm miss}$ [GeV]
YLabel=Fraction of Events / 5 GeV
RatioPlot=0
# + any additional plot settings you might like, see make-plots documentation
END PLOT

BEGIN PLOT /CEPC_SLEPTON/normlizedhist_mll
XLabel=$m_{\ell \ell}$ [GeV]
YLabel=Fraction of Events / 5 GeV
# + any additional plot settings you might like, see make-plots documentation
END PLOT

BEGIN PLOT /CEPC_SLEPTON/normlizedhist_mT2
XLabel=$m_{T2}$ [GeV]
YLabel=Fraction of Events / 5 GeV
# + any additional plot settings you might like, see make-plots documentation
END PLOT
# ... add more histograms as you need them ...