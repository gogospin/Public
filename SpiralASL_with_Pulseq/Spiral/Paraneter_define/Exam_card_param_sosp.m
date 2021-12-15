
%% set parameters

MrProt.sliceGroupList(1).position.dTra = 0;% 14; % H>>F (H positive)
MrProt.sliceGroupList(1).position.dCor = 0;%-20; % P>>A (P positive)
MrProt.sliceGroupList(1).position.dSag = 0;%5; % L>>R (L positive)
MrProt.sliceGroupList(1).thickness = 1.5; %mm
MrProt.sliceGroupList(1).sliceperslab = 32;
MrProt.sliceGroupList(1).shift = 0;%14 ;24
MrProt.sliceGroupList(1).distanceFactor = 0.5;

MrProt.sliceOversampling = 0;%0.15;
MrProt.dynamics = 20;
MrProt.private.rftype = 'sinc'; %'rec'or ''sinc'
MrProt.private.sNormal.TC = -0; % T>C degree 
MrProt.private.pe_cos =cos(deg2rad(MrProt.private.sNormal.TC));
MrProt.private.pe_sin =sin(deg2rad(MrProt.private.sNormal.TC));

MrProt.flipAngle = 20;
l_additionalslice = 2*round(MrProt.sliceGroupList(1).sliceperslab*MrProt.sliceOversampling);
l_additionalthick = l_additionalslice*MrProt.sliceGroupList(1).thickness;
MrProt.sliceSeries.aFront.thickness = l_additionalthick + MrProt.sliceGroupList(1).sliceperslab*MrProt.sliceGroupList(1).thickness;
MrProt.private.l_additionalslice = l_additionalslice;
MrProt.private.sliceOversampling = l_additionalslice/MrProt.sliceGroupList(1).sliceperslab;
clear l_additionalthick l_additionalslice;

MrProt.sliceSeries.size = 48;
MrProt.sliceSeries.Slices = 'Ascending';%'Descending','Interleaved'

MrProt.perfusion = 1;

MrProt.tr = 3000;
MrProt.ti(1) = 700;
MrProt.ti(2) = 1800;
%m_dTI1 = 700;
%m_dTI2 = 1800;

MrProt.asl.FlowLimit = 100;

MrProt.satList.thickness = 100;
MrProt.satList.gap = 25;

MrProt.dim = 3;
if MrProt.dim == 3
    MrProt.sliceSeries.size = MrProt.sliceGroupList(1).sliceperslab; 
    MrProt.sliceGroupList(1).distanceFactor = 0.0;
end
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
MrProt.spcl.g3.Undersampling_fac = 1.6;%1.6;
MrProt.spcl.g3.paramLongInterleaves = 1;

%m_dFatSatFlipangle = 110;

MrProt.spcl.g4.SMS_factor  = 1;
MrProt.spcl.g4.CAIPI_shift = 1;
MrProt.spcl.g4.SMS_online_recon = 1; 
MrProt.spcl.g4.SMS_RF_phase_optim = 1;

MrProt.spcl.g4.CAIPI_shift_mm = 0;%36; % &paramDoubleEchoShift
MrProt.spcl.g4.CAIPI_period_us = 0; %&paramDoubleCAIPIPeriod	
MrProt.spcl.g4.CAIPI_delay_us = 0;	%&paramDoubleCAIPIDelay
	
MrProt.spcl.g5.RF_BWTP = 25;	
MrProt.spcl.private.slicepershot = 1;

MrProt.spcl.single_shot_mode = 1;
MrProt.spcl.accr_kz = 1;
MrProt.spcl.private.rotate_mode_ax = 'linear';% or 'golden_angle'
MrProt.spcl.private.increm_mode_sl = 'linear';% or 'binomial' or 'arbit'
%% to be filled
MrProt.sliceSeries.aFront.readoutFOV = 200;
%MrProt.rxSpec.effDwellTime = 1800;

%% skope
%% process the position informaiton
if MrProt.dim == 3
    MrProt.sliceGroupList(1).distanceFactor = 0;
end
 slab_thickness = (MrProt.sliceSeries.size-1)*MrProt.sliceSeries.aFront.thickness*MrProt.sliceGroupList(1).distanceFactor + ...
     MrProt.sliceSeries.size*MrProt.sliceSeries.aFront.thickness;
 for i = 1:MrProt.sliceSeries.size
     MrProt.sliceGroup(1).pslc(i).SliceShift = -(slab_thickness/2-MrProt.sliceSeries.aFront.thickness/2 - ...
         (i-1)*(MrProt.sliceGroupList(1).distanceFactor*MrProt.sliceSeries.aFront.thickness + MrProt.sliceSeries.aFront.thickness)) + ...
         MrProt.sliceGroupList(1).shift;
 end
%MrProt.sliceGroup(1).pslc(1).SliceShift = MrProt.sliceGroupList(1).shift;
clear i slab_thickness

%% process the gap info

% for i = 1:length(MrProt.sliceGroupList)
%     MrProt.sliceGroupList(i).distance = (MrProt.sliceGroupList(i).distanceFactor+1)*MrProt.sliceGroupList(i).thickness;
% end