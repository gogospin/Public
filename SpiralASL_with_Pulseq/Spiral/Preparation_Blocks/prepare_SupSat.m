%% ##############   prepare SupSat pulse      #############################
% prepare InfSat pulse (slice selection and spoiling
% sus_GPSSel.prepAmplitude(SupSat.sus_lRampTime, SupSat.sus_lRampTime + SupSat.sus_RF.getDuration(), SupSat.sus_lRampTime, SupSat.sus_dAmpGSel))
% sus_GPSpoil.prepAmplitude(sus_lRampTime, sus_lSpoilDuration, sus_lRampTime, sus_dSpoilAmplitude)
%
% SupSat.sus_GPSSel: 
%     rampup time       (rampup)       -> SupSat.sus_lRampTime ; 
%     plateau duration  (plateau)      -> SupSat.sus_RF.Duration ; 
%     rampdown time     (rampdown)     -> SupSat.sus_lRampTime ;
%     duration          (duration)     -> SupSat.sus_GPSSel.rampup +
%     SupSat.sus_GPSSel.rampdown + SupSat.sus_GPSSel.plateau;
%     amplitude         (amplitude)    -> SupSat.sus_dAmpGSel = 0.49;%REF_InfSat_GSEL
%     *100. / SupSat.m_dIRSlabThk; !!!! pas d'accord, les deux
%     axis                -> z
% gradient unit: 
%      time: seconds; 
%      amplitude: kHz/m; (the result of 'MHz/T times mT/m' has the same value as 'kHz/m';
%  sus_GPSSel.prepAmplitude(sus_lRampTime, sus_lRampTime + sus_RF.getDuration(), sus_lRampTime, sus_dAmpGSel))
SupSat.sus_GPSSel.rampup = SupSat.sus_lRampTime*1e-6 ; 
SupSat.sus_GPSSel.rampdown = SupSat.sus_lRampTime*1e-6 ; 
SupSat.sus_GPSSel.plateau = SupSat.sus_RF.Duration*1e-6 ; 
SupSat.sus_GPSSel.duration = SupSat.sus_GPSSel.rampup + SupSat.sus_GPSSel.rampdown + SupSat.sus_GPSSel.plateau ;
SupSat.sus_GPSSel.amplitude = REF_SUPSAT_GSEL * 100.0  /SupSat.m_dSupSatThk * larmor_freq; %mT/m  
SupSat.sus_GPSSel.RampUpTime = SupSat.sus_lRampTime;
SupSat.sus_GPSSel.RampDownTime = SupSat.sus_lRampTime;

SupSat.sus_GPSSel.axix = 'z' ;
SupSat.sus_GPSSel.grad = mr.makeTrapezoid(SupSat.sus_GPSSel.axix,lims,'duration', SupSat.sus_GPSSel.duration, ...
    'flatTime',SupSat.sus_GPSSel.plateau,'riseTime',SupSat.sus_GPSSel.rampup, 'amplitude', SupSat.sus_GPSSel.amplitude);

% read InfSat RF pulse:
%SupSat.rf = read_RF_pta('C:\Users\P70071217\WMWare_share\RF\SAT2560A.SAT_16A2_1.pta');
SupSat.rf = read_RF_pta(fullfile(buildpath,'RF_design\RFdta_file\SAT2560A.SAT_16A2_1.pta'));
SupSat.sus_RF.FlipAngle = SupSat.sus_dFlipAngle * SupSat.sus_dRFscale;
SupSat.sus_RF.InitialPhase = 180.0;

SupSat.sus_RF.Thickness = SupSat.m_dSupSatThk;
SupSat.sus_RF.SampleSize = 512;
SupSat.sus_RF.FamilyName = "SAT2560A.SAT_16A2_1";  %%external pulse
SupSat.sus_RF.Duration = RF_INFSAT_DURATION;     %%pulse duration

SupSat = FPaslSupSat.offcen_preset(MrProt,SupSat);

rep_points = SupSat.sus_RF.Duration/SupSat.sus_RF.SampleSize;
rf_samples = zeros(SupSat.sus_RF.Duration/rfras_us,1);
for i = 1:rep_points
    rf_samples(i:rep_points:end) = SupSat.rf.rf;
end
SupSat.sus_RF.rf = mr.makeArbitraryRf(rf_samples,SupSat.sus_RF.FlipAngle*pi()/180,'system',lims, ...
    'delay',SupSat.sus_GPSSel.rampup,'freqOffset',SupSat.sus_RFSet.Frequency,'phaseOffset',deg2rad(SupSat.sus_RF.InitialPhase));


% sus_GPSpoil:
%     rampup time         -> SupSat.sus_lRampTime ; 
%     plateau duration    -> SupSat.sus_lSpoilDuration - SupSat.sus_lRampTime;  
%     rampdown time       -> SupSat.sus_lRampTime ;
%     amplitude           -> 6 ; 
%     axis                -> y

%sus_GPSpoil.prepAmplitude(sus_lRampTime, sus_lSpoilDuration, sus_lRampTime, sus_dSpoilAmplitude))


SupSat.sus_GPSpoil.rampup = SupSat.sus_lRampTime*1e-6 ; 
SupSat.sus_GPSpoil.rampdown = SupSat.sus_lRampTime*1e-6 ; 
SupSat.sus_GPSpoil.plateau = (SupSat.sus_lSpoilDuration - SupSat.sus_lRampTime) *1e-6 ; 
SupSat.sus_GPSpoil.duration = SupSat.sus_GPSpoil.rampup + SupSat.sus_GPSpoil.rampdown + SupSat.sus_GPSpoil.plateau ;
SupSat.sus_GPSpoil.amplitude = SupSat.sus_dSpoilAmplitude * larmor_freq; %mT/m  
SupSat.sus_GPSpoil.TotalTime = SupSat.sus_lSpoilDuration + SupSat.sus_lRampTime;

SupSat.sus_GPSpoil.axix = 'x' ;
SupSat.sus_GPSpoil.grad = mr.makeTrapezoid(SupSat.sus_GPSpoil.axix,lims,'duration', SupSat.sus_GPSpoil.duration, ...
    'flatTime',SupSat.sus_GPSpoil.plateau,'riseTime',SupSat.sus_GPSpoil.rampup, 'amplitude', SupSat.sus_GPSpoil.amplitude);