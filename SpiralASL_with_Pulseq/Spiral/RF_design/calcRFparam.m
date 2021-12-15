function RF = calcRFparam(RF, MrProt)
% calculate the amplitude and duration ...
%  of Slice slective and refocusing gradient 
% of the RF pulses from Siemens scanner

RF.LarmorConst = 42.5756;
if RF.FamilyName == 'SE2560A90.SE90_12A2_2'
    my_scale_GPSel = 24;
    my_maxslew = 0.24; %mT/m/us
    RF.GPSSel_Amplitude = my_scale_GPSel/MrProt.sliceSeries.aFront.thickness; % in mT/m
    RF.GPSSel_lRampTime = ceil(RF.GPSSel_Amplitude/10/my_maxslew)*10;
    RF.GPSSel_Duration = RF.GPSSel_lRampTime + RF.Duration;
    
    RF.refoc_lRampTime = 160;
    
end
    
    
    