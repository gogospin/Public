RFExcite.predelay = 110;

RFExcite.Thickness = MrProt.sliceSeries.aFront.thickness;
RFExcite.flipAngle = MrProt.flipAngle;
RFExcite.InitialPhase = 90.0;
RFExcite.FamilyName = 'SINC';%'SE2560A90.SE90_12A2_2';
RFExcite.Duration = 2560;


slab_thickness = MrProt.sliceGroupList(1).thickness * MrProt.sliceGroupList(1).sliceperslab;
RFExcite.GPSSel_lRampTime = 10;%50;
RFExcite.GPSSel_Duration = 2570;%2610;
RFExcite.GPSSel_Amplitude = 48/slab_thickness;%1.49;%12.0;  

RFExcite.refoc_lRampTime = 150;%160;
RFExcite.refoc_Duration = 150;%410;%330;
RFExcite.refoc_Amplitude = -RFExcite.GPSSel_Amplitude*RFExcite.GPSSel_Duration/2/RFExcite.refoc_Duration;%47.45;

RFExcite.postdelay = 400;