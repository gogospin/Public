% the parameters are read from POET simulation, pay attention to the
% polarity: it's recorded as RO, PE, 3D (X,Y,Z) 
CSatFat.csat_lRampTime = 500;
% RFSET.Frequency and Phase hard copied from IDEA simulation,
% The values do not depend on slice numbers or position
CSatFat.csat_RFSet.Frequency = -981;
CSatFat.csat_RFSet.Phase = 451.9;
CSatFat.csat_RF.Duration = MrProt.spcl.g5.RF_pulse_duration ;
CSatFat.csat_Duration = 3000; 
CSatFat.csat_Amplitude = 8;
CSatFat.csat_polarity = [-1 1 -1];
CSatFat.csat_RF.FlipAngle = MrProt.spcl.g5.FatSat_flip_angle;
CSatFat.gap_RF = 20; % delay between the end of RF and the beginning of gradient, the beginning of RF and the end of the gradient