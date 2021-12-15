%% ##############   prepare CSatFat pulse     #############################
% prepare CSatFat pulse ainsi le gradient precedant and suivant
% Je lis les chifres de POET mais il faut verifier chaque fois que les
% parametres sont changes.

% CSatFat.csat_grad: 
%     rampup time       (rampup)       -> CSatFat.csat_lRampTime ; 
%     plateau duration  (plateau)      -> CSatFat.csat_RF.Duration ; 
%     rampdown time     (rampdown)     -> CSatFat.csat_lRampTime ;
% 
% gradient unit: 
%      time: seconds; 
%      amplitude: kHz/m; (the result of 'MHz/T times mT/m' has the same value as 'kHz/m';
%  csat_grad.prepAmplitude(sus_lRampTime, sus_lRampTime + csat_RF.getDuration(), sus_lRampTime, sus_dAmpGSel))
CSatFat.csat_grad.rampup = CSatFat.csat_lRampTime*1e-6 ; 
CSatFat.csat_grad.rampdown = CSatFat.csat_lRampTime*1e-6 ; 
CSatFat.csat_grad.plateau = CSatFat.csat_RF.Duration*1e-6 ; 
CSatFat.csat_grad.duration = CSatFat.csat_grad.rampup + CSatFat.csat_grad.rampdown + CSatFat.csat_grad.plateau ;
CSatFat.csat_grad.amplitude = CSatFat.csat_Amplitude * larmor_freq; %mT/m  
CSatFat.csat_grad.RampUpTime = CSatFat.csat_lRampTime;
CSatFat.csat_grad.RampDownTime = CSatFat.csat_lRampTime;


CSatFat.csat_grad.grad_x1 = mr.makeTrapezoid('x',lims,'duration', CSatFat.csat_grad.duration, ...
    'flatTime',CSatFat.csat_grad.plateau,'riseTime',CSatFat.csat_grad.rampup, 'amplitude', CSatFat.csat_grad.amplitude*CSatFat.csat_polarity(1));
CSatFat.csat_grad.grad_x2 = mr.makeTrapezoid('x',lims,'duration', CSatFat.csat_grad.duration, ...
    'flatTime',CSatFat.csat_grad.plateau,'riseTime',CSatFat.csat_grad.rampup, 'amplitude', -CSatFat.csat_grad.amplitude*CSatFat.csat_polarity(1));
CSatFat.csat_grad.grad_y1 = mr.makeTrapezoid('y',lims,'duration', CSatFat.csat_grad.duration, ...
    'flatTime',CSatFat.csat_grad.plateau, 'riseTime',CSatFat.csat_grad.rampup,'amplitude', CSatFat.csat_grad.amplitude*CSatFat.csat_polarity(2));
CSatFat.csat_grad.grad_y2 = mr.makeTrapezoid('y',lims,'duration', CSatFat.csat_grad.duration, ...
    'flatTime',CSatFat.csat_grad.plateau,'riseTime',CSatFat.csat_grad.rampup, 'amplitude', -CSatFat.csat_grad.amplitude*CSatFat.csat_polarity(2));
CSatFat.csat_grad.grad_z1 = mr.makeTrapezoid('z',lims,'duration', CSatFat.csat_grad.duration, ...
    'flatTime',CSatFat.csat_grad.plateau,'riseTime',CSatFat.csat_grad.rampup, 'amplitude', CSatFat.csat_grad.amplitude*CSatFat.csat_polarity(3));
CSatFat.csat_grad.grad_z2 = mr.makeTrapezoid('z',lims,'duration', CSatFat.csat_grad.duration, ...
    'flatTime',CSatFat.csat_grad.plateau,'riseTime',CSatFat.csat_grad.rampup, 'amplitude', -CSatFat.csat_grad.amplitude*CSatFat.csat_polarity(3));

% read InfSat RF pulse:

CSatFat.rf = read_RF_pta(fullfile(buildpath,'RF_design\RFdta_file\GAUSS5120.B375.pta'));
CSatFat.csat_RF.FlipAngle = MrProt.spcl.g5.FatSat_flip_angle;

%CSatFat.csat_RF.InitialPhase = 180.0;

%CSatFat.csat_RF.Thickness = CSatFat.m_dInfSatThk;
CSatFat.csat_RF.SampleSize = 512;
%CSatFat.csat_RF.FamilyName = "SAT2560A.SAT_16A2_1";  %%external pulse

rep_points = CSatFat.csat_RF.Duration/CSatFat.csat_RF.SampleSize;
rf_samples = zeros(CSatFat.csat_RF.Duration/rfras_us,1);
for i = 1:rep_points
    rf_samples(i:rep_points:end) = CSatFat.rf.rf;
end
CSatFat.csat_RF.rf = mr.makeArbitraryRf(rf_samples,CSatFat.csat_RF.FlipAngle*pi()/180, ...
    'system',lims,'freqOffset',CSatFat.csat_RFSet.Frequency,'phaseOffset',deg2rad(CSatFat.csat_RFSet.Phase));


CSatFat.delay1 = mr.makeDelay(CSatFat.gap_RF*1e-6);
CSatFat.delay2 = mr.makeDelay(CSatFat.gap_RF*1e-6);

CSatFat.TotalTime = 2*CSatFat.csat_lRampTime + CSatFat.csat_RF.Duration + 2*CSatFat.csat_Duration + ...
    2*CSatFat.gap_RF;
