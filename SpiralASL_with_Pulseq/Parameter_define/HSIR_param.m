HSIR.m_dThickness = 50.0; %% default 50.0mm
HSIR.m_dFlipAngle = 240.0; %% 240 for inversion
HSIR.m_lGSSRampTime_us = 400;	%% default ramp time of slice select gradient
HSIR.m_dShift = 0.0;
HSIR.m_dShiftNonSelective = 0.0; %% 'shift' for non-selective control pulse, default 0.0
HSIR.m_dIRSlabThk = VOL_THICK_INV;       %% mm
HSIR.m_dIRSlabThkNS = VOL_THICK_CON;     %% mm
%HSIR.m_eLabelState(ASLSTATE::ASL_CONTROL)
HSIR.ss_dRFscale = 1.0;
HSIR.ss_dFlipAngle = RF_FLIPANGLE;
HSIR.ss_lRampTime = 1000;
HSIR.ss_lSpoilDuration = 5000;
HSIR.ss_dSpoilAmplitude = 3.0;

HSIR.m_sSRF_FOCI_LBL.Thickness = HSIR.m_dThickness;
HSIR.m_sSRF_FOCI_LBL.Duration = 10240; %% Renzo variable Inv Pulse Duration ######## Attention: this pulse is 10 ms long !!!

% Renzo samples festsetzen uf 1/2 duration
HSIR.m_sSRF_FOCI_LBL.SampleSize = 5120;  % sRF_PULSE_ARB.setSamples  %% Achtung, musst du vielleicht verdoppeln 30.10.13 // von 2560 auf 5120
% double sampleCount = m_sSRF_FOCI_LBL.getDuration() / 2.0;
%m_sSRF_FOCI_LBL.setSamples        (sampleCount);
HSIR.m_sSRF_FOCI_LBL.FlipAngle = HSIR.m_dFlipAngle;
HSIR.m_sSRF_FOCI_LBL.InitialPhase = 0.0;
HSIR.m_sSRF_FOCI_LBL.RequiredGSPolarity = 1.0;
HSIR.m_sSRF_FOCI_LBL.LarmorConst = 42.5756;

HSIR.m_sSRF_FOCI_CTRL = HSIR.m_sSRF_FOCI_LBL;

HSIR.m_sFociData.m_dB1Cutoff = 0.004;
HSIR.m_sFociData.m_dFOCI_Cutoff = 20.0;
HSIR.m_sFociData.m_dMu = 5.0;
HSIR.m_sFociData.m_lGradRampTime = HSIR.m_lGSSRampTime_us;

HSIR.m_sGSS_ARB = [];
