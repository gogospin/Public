
%% set parameters

MrProt.sliceGroupList(1).position.dTra = 0;% 14; % H>>F (H positive)
MrProt.sliceGroupList(1).position.dCor = 0;%-20; % P>>A (P positive)
MrProt.sliceGroupList(1).position.dSag = 0;%5; % L>>R (L positive)
MrProt.sliceGroupList(1).thickness = 2.00; %mm
MrProt.sliceGroupList(1).shift = 0;%14 ;
MrProt.sliceGroupList(1).distanceFactor = 0.5;

MrProt.flipAngle = 65;
MrProt.sliceSeries.aFront.thickness = 2;

MrProt.sliceSeries.size = 24;
MrProt.sliceSeries.Slices = 'Ascending';%'Descending','Interleaved'

MrProt.tr = 5000;
MrProt.ti(1) = 700;
MrProt.ti(2) = 1800;
%m_dTI1 = 700;
%m_dTI2 = 1800;

MrProt.asl.FlowLimit = 100;

% the following 2 parameters are actually not in the scan exam card,
% but are defined (but not registered) in the SBBPaslUILink. I don't know
% why. I just put the hard code here
MrProt.spcl.g2.Inferior_Sat_Thk = 100 ;
MrProt.spcl.g2.Inferior_Sat_Gap = 0;
MrProt.spcl.g2.Inversion_Slab_Thk = 121;

MrProt.spcl.g1.trfociAmpl = 100;
MrProt.spcl.g1.trfociBWDTH = 300;
MrProt.spcl.g1.trfociD = 100;

MrProt.spcl.g5.FatSat_flip_angle = 110;
MrProt.spcl.g5.RF_pulse_duration = 2560;

% sequence -> special card

MrProt.spcl.g3.RO_samples = 10240;
MrProt.spcl.g3.Spiral_Type = -2;
MrProt.spcl.g3.Spiral_Peak_Grad = 35;
MrProt.spcl.g3.Spiral_Slew_Rate = 153;
MrProt.spcl.g3.Interleaves = 1;
MrProt.spcl.g3.Flat_First_ZC = 0;
MrProt.spcl.g3.Spiral_BW = 400000;
MrProt.spcl.g3.VD = 1.3;
MrProt.spcl.g3.Undersampling_fac = 1.6;
MrProt.spcl.g3.paramLongInterleaves = 1;

%m_dFatSatFlipangle = 110;

MrProt.spcl.g4.SMS_factor  = 2;
MrProt.spcl.g4.CAIPI_shift = 1;
MrProt.spcl.g4.SMS_online_recon = 1; 
MrProt.spcl.g4.SMS_RF_phase_optim = 1;

MrProt.spcl.g4.CAIPI_shift_mm = 36;%36; % &paramDoubleEchoShift
MrProt.spcl.g4.CAIPI_period_us = 200; %&paramDoubleCAIPIPeriod	
MrProt.spcl.g4.CAIPI_delay_us = 200;	%&paramDoubleCAIPIDelay
	
	
%% to be filled
MrProt.sliceSeries.aFront.readoutFOV = 200;
MrProt.rxSpec.effDwellTime = 1800;

%% skope
%% process the position informaiton
slab_thickness = (MrProt.sliceSeries.size-1)*MrProt.sliceSeries.aFront.thickness*MrProt.sliceGroupList(1).distanceFactor + ...
    MrProt.sliceSeries.size*MrProt.sliceSeries.aFront.thickness;
for i = 1:MrProt.sliceSeries.size
    MrProt.sliceGroup(1).pslc(i).SliceShift = -(slab_thickness/2-MrProt.sliceSeries.aFront.thickness/2 - ...
        (i-1)*(MrProt.sliceGroupList(1).distanceFactor*MrProt.sliceSeries.aFront.thickness + MrProt.sliceSeries.aFront.thickness)) + ...
        MrProt.sliceGroupList(1).shift;
end
clear i slab_thickness

%% process the gap info

for i = 1:length(MrProt.sliceGroupList)
    MrProt.sliceGroupList(i).distance = (MrProt.sliceGroupList(i).distanceFactor+1)*MrProt.sliceGroupList(i).thickness;
end