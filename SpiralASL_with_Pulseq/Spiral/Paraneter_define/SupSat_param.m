SupSat.m_dSupSatThk = VOL_THICK_SUPSAT;    %%// mm
SupSat.m_dSupSatGap = VOL_GAP_SUPSAT;      %%// mm
SupSat.m_dIRSlabThk = VOL_THICK_INV;       %%// mm
%SupSat.m_eLabelState(ASLSTATE::ASL_CONTROL)
SupSat.sus_dRFscale = 1.0;
SupSat.sus_dFlipAngle = RF_SUPSAT_FLIPANGLE;
SupSat.sus_lRampTime = 1000;
SupSat.sus_lSpoilDuration = 5000;
SupSat.sus_dSpoilAmplitude = 6.0;
%SupSat.sus_RFSet("sus_RFSet")
%SupSat.sus_RFNeg("sus_RFNeg")

SupSat.sus_RF.FlipAngle = SupSat.sus_dFlipAngle * SupSat.sus_dRFscale;
SupSat.sus_RF.InitialPhase = 180.0;
SupSat.sus_RF.Thickness = SupSat.m_dSupSatThk;
SupSat.sus_RF.SampleSize = 512;
%SupSat.sus_RF.setFamilyName  ("SAT2560A.SAT_16A2_1");  //external pulse
SupSat.sus_RF.Duration = RF_SUPSAT_DURATION;     %//pulse duration
%sus_RF.setIdent (ptIdent);

SupSat.sus_RF.LarmorConst = 42.5775;
    
SupSat.REF_SUPSAT_GSEL = REF_SUPSAT_GSEL;
%prepare the external pulse
SupSat.sus_dAmpGSel = REF_SUPSAT_GSEL * 100.0 /SupSat.m_dSupSatThk;
SupSat.sus_RF.GSAmplitude = SupSat.sus_dAmpGSel;

SupSat.m_dSupSatGap	  = MrProt.spcl.g2.Inferior_Sat_Gap ;
%SupSat.m_dIRSlabThk = MrProt.spcl.g2.Inversion_Slab_Thk ;
SupSat.m_dIRSlabThk = PASL.m_dFairIRSlabThk;