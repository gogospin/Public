%% ##############   prepare excitation pulse  #############################
% design slice selective gradient and rf

rf0 = mr.makeBlockPulse(MrProt.flipAngle*pi/180,lims,'Duration',MrProt.spcl.g5.RF_pulse_duration/1e6);



RFExcite.RF.Asym = 0.5;% symmetric phase
RFExcite.RF.InitialPhase = 90;
RFExcite.RF.Duration = RFExcite.Duration;
% the following two gradients do not exist in the case of rec RF,
% the duratiuon is set here for the calculation of sequence timing
RFExcite.exc_GPSSel.duration = RFExcite.RF.Duration*1e-6;
RFExcite.refoc.duration = 0;

% refocusing grad
