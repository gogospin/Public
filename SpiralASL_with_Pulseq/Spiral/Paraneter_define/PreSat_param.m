PreSat.m_dIRSlabThk = 100;    % get this value from running poet,PulSeq: PreSat: 
                              % m_dIRSlabThk=pMrProt->wipMemBlock().adFree[UILink::WIP_DoubleIRSlabThk] mm
PreSat.m_dIRSlabThk = PASL.m_dFairIRSlabThk;

%m_eLabelState(ASLSTATE::ASL_CONTROL)
PreSat.pre_dRFscale = 1.0;
PreSat.pre_dFlipAngle = RF_PRESAT_FLIPANGLE;
PreSat.pre_lRampTime = 1000;
PreSat.pre_lSpoilDuration = 5000;
PreSat.pre_dSpoilAmplitude = 6.0;

PreSat.pre_RF.FlipAngle = PreSat.pre_dFlipAngle * PreSat.pre_dRFscale;
PreSat.pre_RF.InitialPhase = 180.0;
PreSat.pre_RF.LarmorConst = 42.5756;

PreSat.pre_RF.setThickness = PreSat.m_dIRSlabThk;
PreSat.pre_RF.setSamples = 512;
	
PreSat.REF_PRESAT_GSEL = REF_PRESAT_GSEL;
PreSat.pre_RF.Duration = RF_PRESAT_DURATION; 

  if( PreSat.m_dIRSlabThk ~= 0)
  	
	PreSat.pre_RF.GSAmplitude = REF_PRESAT_GSEL * 100. / PreSat.m_dIRSlabThk;
  end