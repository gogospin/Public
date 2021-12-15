%% ##############   prepare excitation pulse  #############################
% design slice selective gradient and rf
slice_thickness = MrProt.sliceGroupList(1).thickness * MrProt.sliceGroupList(1).sliceperslab / 1000;
[rf0, gz0, gzr0] = mr.makeSincPulse(MrProt.flipAngle/180*pi,'SliceThickness',slice_thickness, ...
    'Duration',MrProt.spcl.g5.RF_pulse_duration/1e6,'system',lims,'timeBwProduct',MrProt.spcl.g5.RF_BWTP,'apodization',0.5);
clear slice_thickness

RFExcite.exc_GPSSel.grad = gz0; 

RFExcite.exc_GPSSel.rampup = gz0.riseTime;% RFExcite.GPSSel_lRampTime*1e-6 ;
RFExcite.exc_GPSSel.rampdown = gz0.fallTime;%RFExcite.GPSSel_lRampTime*1e-6 ;
RFExcite.exc_GPSSel.plateau = gz0.flatTime;%RFExcite.Duration*1e-6 ;
RFExcite.exc_GPSSel.duration = RFExcite.exc_GPSSel.rampup + RFExcite.exc_GPSSel.rampdown + RFExcite.exc_GPSSel.plateau ;
RFExcite.GPSSel_Amplitude = RFExcite.exc_GPSSel.grad.amplitude / larmor_freq ; %mT/m
RFExcite.RF.GPSSel_Amplitude = RFExcite.GPSSel_Amplitude;
RFExcite.exc_GPSSel.RampUpTime = RFExcite.exc_GPSSel.rampup*1e6;%RFExcite.GPSSel_lRampTime*1e6;
RFExcite.exc_GPSSel.RampDownTime = RFExcite.exc_GPSSel.rampdown*1e6;%RFExcite.GPSSel_lRampTime*1e6;

RFExcite.RF.Asym = 0.5;% symmetric phase
RFExcite.RF.InitialPhase = 90;
RFExcite.RF.Duration = RFExcite.Duration;


% refocusing grad

RFExcite.refoc.rampup = gzr0.riseTime;% RFExcite.refoc_lRampTime*1e-6 ; 
RFExcite.refoc.rampdown = gzr0.fallTime;% RFExcite.refoc_lRampTime*1e-6 ; 
RFExcite.refoc.plateau = gzr0.flatTime;%(RFExcite.refoc_Duration - RFExcite.refoc_lRampTime)*1e-6 ; 
RFExcite.refoc.duration = RFExcite.refoc.rampup + RFExcite.refoc.rampdown + RFExcite.refoc.plateau ;
RFExcite.refoc_Amplitude = gzr0.amplitude /larmor_freq ; %mT/m  
RFExcite.refoc.RampUpTime = RFExcite.refoc.rampup*1e6;%RFExcite.refoc_lRampTime;
RFExcite.refoc.RampDownTime = RFExcite.refoc.rampdown*1e6;%RFExcite.refoc_lRampTime;

RFExcite.refoc.grad = gzr0; 

clear gz0; clear gzr0;