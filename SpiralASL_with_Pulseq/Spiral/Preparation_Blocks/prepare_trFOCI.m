%% ##############   prepare trFOCI pulse      #############################
%###########################################
%     MAKE FOCI pulse: m_sSRF_FOCI_LBL
%     m_sSRF_FOCI_CTRL, m_sFociData
[HSIR.m_sFociData,HSIR.m_sSRF_FOCI_LBL,HSIR.m_sSRF_FOCI_CTRL,HSIR.m_sGSS_ARB] ...
    = FPaslHSIR.foci_gen(MrProt, HSIR.m_sSRF_FOCI_LBL, HSIR.m_sSRF_FOCI_CTRL, HSIR.m_sGSS_ARB,HSIR.m_sFociData);

HSIR.m_sGSS.axis = 'z'; 
HSIR.m_sGSS.StartTime = 0;
HSIR.m_sGSS_ARB.axis = 'z'; 
HSIR.m_sGSS_ARB.StartTime = 0;

HSIR.m_sGSS_ARB.RampUpTime = HSIR.m_sFociData.m_lGradRampTime;
HSIR.m_sGSS.RampUpTime = HSIR.m_sFociData.m_lGradRampTime;


rep_points = HSIR.m_sSRF_FOCI_LBL.Duration/HSIR.m_sSRF_FOCI_LBL.SampleSize;
rf_samples = zeros(HSIR.m_sSRF_FOCI_LBL.Duration/rfras_us+20,1);% changed to correct for gap
for i = 1:rep_points
    rf_samples(10+i:rep_points:end-10) = HSIR.m_sSRF_FOCI_LBL.Samples;
end
 
% I don't understand why the frequency offset is zeros (see cpp code
% offcen_preset, Stefan made it zero...)
%lims_hsir = lims;lims_hsir.rfRasterTime = HSIR.m_sSRF_FOCI_LBL.Duration/HSIR.m_sSRF_FOCI_LBL.SampleSize*1e-6;
rf_hsir_label = mr.makeArbitraryRf(rf_samples*HSIR.m_sFociData.B1max,HSIR.m_dFlipAngle/HSIR.m_sFociData.FA_scale*pi()/180,'system',lims, ...
    'delay',(HSIR.m_sFociData.m_lGradRampTime-10)*1e-6);

gradlabel = zeros(length(HSIR.m_sGSS_ARB.Samples),1);
gradlabel(1:end) = HSIR.m_sGSS_ARB.Samples; 
%gz_hsir_label = mr.makeArbitraryGrad(HSIR.m_sGSS_ARB.axis,HSIR.m_sGSS_ARB.Samples*larmor_freq,'system',lims);
gz_hsir_label = mr.makeArbitraryGrad(HSIR.m_sGSS_ARB.axis,gradlabel*larmor_freq,'system',lims);
gz_hsir_control = mr.makeDelay(HSIR.m_sGSS_ARB.Duration*1e-6);
%HSIR.m_sGSS_ARB.Duration = HSIR.m_sGSS_ARB.Duration+20;

HSIR.ss_GPSpoil.rampup = HSIR.ss_lRampTime*1e-6 ; 
HSIR.ss_GPSpoil.rampdown = HSIR.ss_lRampTime*1e-6 ; 
HSIR.ss_GPSpoil.plateau = (HSIR.ss_lSpoilDuration - HSIR.ss_lRampTime) * 1e-6 ; 
HSIR.ss_GPSpoil.duration = HSIR.ss_GPSpoil.rampup + HSIR.ss_GPSpoil.rampdown + HSIR.ss_GPSpoil.plateau ;
HSIR.ss_GPSpoil.amplitude = HSIR.ss_dSpoilAmplitude * larmor_freq; %mT/m 
HSIR.ss_GPSpoil.grad = mr.makeTrapezoid('x',lims,'duration', HSIR.ss_GPSpoil.duration, ...
    'flatTime',HSIR.ss_GPSpoil.plateau, 'riseTime', HSIR.ss_GPSpoil.rampup,'amplitude', HSIR.ss_GPSpoil.amplitude);

% calculate timing
HSIR.m_lSBBDurationPerRequest_us = ...
    HSIR.m_sGSS_ARB.Duration + HSIR.ss_GPSpoil.duration + 260 * MrProt.sliceSeries.size; 