RFExcite.predelay = 110;

RFExcite.Thickness = MrProt.sliceSeries.aFront.thickness;
RFExcite.flipAngle = MrProt.flipAngle;
RFExcite.InitialPhase = 90.0;
RFExcite.FamilyName = 'SE2560A90.SE90_12A2_2';
RFExcite.Duration = 2560;

RFExcite.GPSSel_lRampTime = 50;
RFExcite.GPSSel_Duration = 2610;
RFExcite.GPSSel_Amplitude = 12.0;  

RFExcite.refoc_lRampTime = 230;%160;
RFExcite.refoc_Duration = 790;%330;
RFExcite.refoc_Amplitude = -RFExcite.GPSSel_Amplitude*RFExcite.GPSSel_Duration/2/RFExcite.refoc_Duration;%47.45;

RFExcite.postdelay = 400;