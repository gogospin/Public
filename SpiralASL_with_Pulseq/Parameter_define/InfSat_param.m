InfSat.m_dInfSatThk = VOL_THICK_INFSAT;    %% mm
InfSat.m_dInfSatGap = VOL_GAP_INFSAT;     %% mm
InfSat.m_dIRSlabThk = VOL_THICK_INV;       %% mm
%InfSat.m_eLabelState(ASLSTATE::ASL_CONTROL)
InfSat.ifs_dRFscale = 1.0;
InfSat.ifs_dFlipAngle = RF_INFSAT_FLIPANGLE;
InfSat.ifs_lRampTime = 1000;
InfSat.ifs_lSpoilDuration = 5000;
InfSat.ifs_dSpoilAmplitude = 6.0;
InfSat.REF_InfSat_GSEL = REF_INFSAT_GSEL;

%InfSat.ifs_RFSet("ifs_RFSet")
%InfSat.ifs_RFNeg("ifs_RFNeg")

InfSat.ifs_RF.FlipAngle = InfSat.ifs_dFlipAngle * InfSat.ifs_dRFscale;
InfSat.ifs_RF.InitialPhase = 180.0;
InfSat.ifs_RF.Thickness = InfSat.m_dInfSatThk; %%.maolin - sat thickness for perf quantitation
InfSat.ifs_RF.SampleSize = 512;
%InfSat.ifs_RF.setFamilyName    ("SAT2560A.SAT_36A2_1");  //external pulse
InfSat.ifs_RF.Duration = RF_INFSAT_DURATION;
if InfSat.m_dInfSatThk ~= 0
    InfSat.dAmpGSel = REF_INFSAT_GSEL * 100.0  / InfSat.m_dInfSatThk;
    InfSat.ifs_RF.GSAmplitude = InfSat.dAmpGSel;
end

InfSat.dTI1 = MrProt.ti(2)/1000; % in ms
InfSat.m_dTI2 = MrProt.ti(1)/1000 - InfSat.dTI1;

InfSat.m_dInfSatThk	 = MrProt.spcl.g2.Inferior_Sat_Thk ;
InfSat.m_dInfSatGap = MrProt.spcl.g2.Inferior_Sat_Gap ;
%InfSat.m_dIRSlabThk = MrProt.spcl.g2.Inversion_Slab_Thk ;
InfSat.m_dIRSlabThk = PASL.m_dFairIRSlabThk;