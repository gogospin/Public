%% ##############   prepare PreSat pulse      #############################
% prepare PreSat pulse (slice selection and spoiling
% pre_GPSSel.prepAmplitude(PreSat.pre_lRampTime, PreSat.pre_lRampTime + PreSat.pre_RF.getDuration(), PreSat.pre_lRampTime, PreSat.pre_dAmpGSel))
% pre_GPSpoil.prepAmplitude(pre_lRampTime, pre_lSpoilDuration, pre_lRampTime, pre_dSpoilAmplitude)
%
% PreSat.pre_GPSSel: 
%     rampup time       (rampup)       -> PreSat.pre_lRampTime ; 
%     plateau duration  (plateau)      -> PreSat.pre_RF.Duration ; 
%     rampdown time     (rampdown)     -> PreSat.pre_lRampTime ;
%     duration          (duration)     -> PreSat.pre_GPSSel.rampup +
%     PreSat.pre_GPSSel.rampdown + PreSat.pre_GPSSel.plateau;
%     amplitude         (amplitude)    -> PreSat.pre_dAmpGSel = 0.49;%REF_PRESAT_GSEL
%     *100. / PreSat.m_dIRSlabThk; !!!! pas d'accord, les deux
%     axis                -> z
% gradient unit: 
%      time: seconds; 
%      amplitude: kHz/m; (the result of 'MHz/T times mT/m' has the same value as 'kHz/m';
PreSat.pre_GPSSel.rampup = PreSat.pre_lRampTime*1e-6 ; 
PreSat.pre_GPSSel.rampdown = PreSat.pre_lRampTime*1e-6 ; 
PreSat.pre_GPSSel.plateau = PreSat.pre_RF.Duration*1e-6 ; 
PreSat.pre_GPSSel.duration = PreSat.pre_GPSSel.rampup + PreSat.pre_GPSSel.rampdown + PreSat.pre_GPSSel.plateau ;
PreSat.pre_GPSSel.amplitude = REF_PRESAT_GSEL*100 / PreSat.m_dIRSlabThk * larmor_freq; %mT/m  
PreSat.pre_GPSSel.RampUpTime = PreSat.pre_lRampTime;
PreSat.pre_GPSSel.RampDownTime = PreSat.pre_lRampTime;

PreSat.pre_GPSSel.axix = 'z' ;
PreSat.pre_GPSSel.grad = mr.makeTrapezoid(PreSat.pre_GPSSel.axix,lims,'duration', PreSat.pre_GPSSel.duration, ...
    'flatTime',PreSat.pre_GPSSel.plateau,'riseTime', PreSat.pre_GPSSel.rampup, 'amplitude',PreSat.pre_GPSSel.amplitude);

% read presat RF pulse:
%PreSat.rf = read_RF_pta('C:\Users\P70071217\WMWare_share\RF\SAT2560A.SAT_36A2_1.pta');
PreSat.rf = read_RF_pta(fullfile(buildpath,'RF_design\RFdta_file\SAT2560A.SAT_36A2_1.pta'));
PreSat.pre_RF.FlipAngle = PreSat.pre_dFlipAngle * PreSat.pre_dRFscale;
PreSat.pre_RF.InitialPhase = 180.0;

PreSat.pre_RF.Thickness = PreSat.m_dIRSlabThk;
PreSat.pre_RF.SampleSize = 512;
PreSat.pre_RF.FamilyName = "SAT2560A.SAT_36A2_1";  %%external pulse
PreSat.pre_RF.Duration = RF_PRESAT_DURATION;     %%pulse duration

rep_points = PreSat.pre_RF.Duration/PreSat.pre_RF.SampleSize;
rf_samples = zeros(PreSat.pre_RF.Duration/rfras_us,1);
for i = 1:rep_points
    rf_samples(i:rep_points:end) = PreSat.rf.rf;
end

PreSat = FPaslPreSat.offcen_preset(MrProt, PreSat);

PreSat.pre_RF.rf = mr.makeArbitraryRf(rf_samples,PreSat.pre_RF.FlipAngle*pi()/180,'system',lims, ...
    'delay',PreSat.pre_GPSSel.rampup,'freqOffset',PreSat.pre_RFSet.Frequency,'phaseOffset',deg2rad(PreSat.pre_RF.InitialPhase));


% pre_GPSpoil:
%     rampup time         -> PreSat.pre_lRampTime ; 
%     plateau duration    -> PreSat.pre_lSpoilDuration - PreSat.pre_lRampTime;  
%     rampdown time       -> PreSat.pre_lRampTime ;
%     amplitude           -> 6 ; 
%     axis                -> y


PreSat.pre_GPSpoil.rampup = PreSat.pre_lRampTime*1e-6 ; 
PreSat.pre_GPSpoil.rampdown = PreSat.pre_lRampTime*1e-6 ; 
PreSat.pre_GPSpoil.plateau = (PreSat.pre_lSpoilDuration - PreSat.pre_lRampTime) *1e-6 ; 
PreSat.pre_GPSpoil.duration = PreSat.pre_GPSpoil.rampup + PreSat.pre_GPSpoil.rampdown + PreSat.pre_GPSpoil.plateau ;
PreSat.pre_GPSpoil.amplitude = PreSat.pre_dSpoilAmplitude * larmor_freq; %mT/m  
PreSat.pre_GPSpoil.TotalTime = PreSat.pre_lRampTime + PreSat.pre_lSpoilDuration;

PreSat.pre_GPSpoil.axix = 'x' ;
PreSat.pre_GPSpoil.grad = mr.makeTrapezoid(PreSat.pre_GPSpoil.axix,lims,'duration', PreSat.pre_GPSpoil.duration, ...
    'flatTime',PreSat.pre_GPSpoil.plateau, 'riseTime',PreSat.pre_GPSpoil.rampup, 'amplitude',PreSat.pre_GPSpoil.amplitude);