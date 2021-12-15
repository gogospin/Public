%% ##############   prepare InfSat pulse      #############################
% prepare InfSat pulse (slice selection and spoiling
% ifs_GPSSel.prepAmplitude(InfSat.ifs_lRampTime, InfSat.ifs_lRampTime + InfSat.ifs_RF.getDuration(), InfSat.ifs_lRampTime, InfSat.ifs_dAmpGSel))
% ifs_GPSpoil.prepAmplitude(ifs_lRampTime, ifs_lSpoilDuration, ifs_lRampTime, ifs_dSpoilAmplitude)
%
% InfSat.ifs_GPSSel: 
%     rampup time       (rampup)       -> InfSat.ifs_lRampTime ; 
%     plateau duration  (plateau)      -> InfSat.ifs_RF.Duration ; 
%     rampdown time     (rampdown)     -> InfSat.ifs_lRampTime ;
%     duration          (duration)     -> InfSat.ifs_GPSSel.rampup +
%     InfSat.ifs_GPSSel.rampdown + InfSat.ifs_GPSSel.plateau;
%     amplitude         (amplitude)    -> InfSat.ifs_dAmpGSel = 0.49;%REF_InfSat_GSEL
%     *100. / InfSat.m_dIRSlabThk; !!!! pas d'accord, les deux
%     axis                -> z
% gradient unit: 
%      time: seconds; 
%      amplitude: kHz/m; (the result of 'MHz/T times mT/m' has the same value as 'kHz/m';
%  ifs_GPSSel.prepAmplitude(ifs_lRampTime, ifs_lRampTime + ifs_RF.getDuration(), ifs_lRampTime, dAmpGSel))
InfSat.ifs_GPSSel.rampup = InfSat.ifs_lRampTime*1e-6 ; 
InfSat.ifs_GPSSel.rampdown = InfSat.ifs_lRampTime*1e-6 ; 
InfSat.ifs_GPSSel.plateau = InfSat.ifs_RF.Duration*1e-6 ; 
InfSat.ifs_GPSSel.duration = InfSat.ifs_GPSSel.rampup + InfSat.ifs_GPSSel.rampdown + InfSat.ifs_GPSSel.plateau ;
InfSat.ifs_GPSSel.amplitude = REF_INFSAT_GSEL * 100.0  /InfSat.m_dInfSatThk * larmor_freq; %mT/m  
InfSat.ifs_GPSSel.RampUpTime = InfSat.ifs_lRampTime ; 
InfSat.ifs_GPSSel.RampDownTime = InfSat.ifs_lRampTime ; 

InfSat.ifs_GPSSel.axix = 'z' ;
InfSat.ifs_GPSSel.grad = mr.makeTrapezoid(InfSat.ifs_GPSSel.axix,lims,'duration', InfSat.ifs_GPSSel.duration, ...
    'flatTime',InfSat.ifs_GPSSel.plateau,'riseTime',InfSat.ifs_GPSSel.rampup, 'amplitude', InfSat.ifs_GPSSel.amplitude);

% read InfSat RF pulse:
InfSat.rf = read_RF_pta(fullfile(buildpath,'RF_design\RFdta_file\SAT2560A.SAT_36A2_1.pta'));
%InfSat.rf = read_RF_pta('C:\Users\P70071217\WMWare_share\RF\SAT2560A.SAT_36A2_1.pta');
InfSat.ifs_RF.FlipAngle = InfSat.ifs_dFlipAngle * InfSat.ifs_dRFscale;
InfSat.ifs_RF.InitialPhase = 180.0;
InfSat.ifs_RF.Asymmetry = 0.5;

InfSat.ifs_RF.Thickness = InfSat.m_dInfSatThk;
InfSat.ifs_RF.SampleSize = 512;
InfSat.ifs_RF.FamilyName = "SAT2560A.SAT_36A2_1";  %%external pulse
InfSat.ifs_RF.Duration = RF_INFSAT_DURATION;     %%pulse duration
InfSat.ifs_RF.LarmorConst = 42.5775;

InfSat = FPaslInfSat.offcen_preset(MrProt,InfSat);

rep_points = InfSat.ifs_RF.Duration/InfSat.ifs_RF.SampleSize;
rf_samples = zeros(PreSat.pre_RF.Duration/rfras_us,1);
for i = 1:rep_points
    rf_samples(i:rep_points:end) = InfSat.rf.rf;
end


InfSat.ifs_RF.rf = mr.makeArbitraryRf(rf_samples,InfSat.ifs_RF.FlipAngle*pi()/180,'system',lims, ...
    'delay',InfSat.ifs_GPSSel.rampup,'freqOffset',InfSat.ifs_RFSet.Frequency,'phaseOffset',deg2rad(InfSat.ifs_RF.InitialPhase));


% ifs_GPSpoil:
%     rampup time         -> InfSat.ifs_lRampTime ; 
%     plateau duration    -> InfSat.ifs_lSpoilDuration - InfSat.ifs_lRampTime;  
%     rampdown time       -> InfSat.ifs_lRampTime ;
%     amplitude           -> 6 ; 
%     axis                -> y

%ifs_GPSpoil.prepAmplitude(ifs_lRampTime, ifs_lSpoilDuration, ifs_lRampTime, ifs_dSpoilAmplitude))


InfSat.ifs_GPSpoil.rampup = InfSat.ifs_lRampTime*1e-6 ; 
InfSat.ifs_GPSpoil.rampdown = InfSat.ifs_lRampTime*1e-6 ; 
InfSat.ifs_GPSpoil.plateau = (InfSat.ifs_lSpoilDuration - InfSat.ifs_lRampTime) *1e-6 ; 
InfSat.ifs_GPSpoil.duration = InfSat.ifs_GPSpoil.rampup + InfSat.ifs_GPSpoil.rampdown + InfSat.ifs_GPSpoil.plateau ;
InfSat.ifs_GPSpoil.amplitude = InfSat.ifs_dSpoilAmplitude * larmor_freq; %mT/m  
InfSat.ifs_GPSpoil.TotalTime = InfSat.ifs_lSpoilDuration + InfSat.ifs_lRampTime;
InfSat.ifs_GPSpoil.RampUpTime = InfSat.ifs_lRampTime;


InfSat.ifs_GPSpoil.axix = 'x' ;
InfSat.ifs_GPSpoil.grad = mr.makeTrapezoid(InfSat.ifs_GPSpoil.axix,lims,'duration', InfSat.ifs_GPSpoil.duration, ...
    'flatTime',InfSat.ifs_GPSpoil.plateau,'riseTime',InfSat.ifs_GPSpoil.rampup, 'amplitude', InfSat.ifs_GPSpoil.amplitude);