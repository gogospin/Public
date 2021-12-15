
MAX_ADCS = 80;
MAX_TRAJ = 20481; 	   %max gradient shape length
RAMP_DOWN_POINTS = 50;   % Ramp down of the gradients at the end of a spiral gradient
MAX_DIFF_DIR = 48;
SQRT_TWO = 1.4142137;
BASE_RES = 2;

MAX_SLEW_RATE_SP = 400;
MAX_GRAD_AMP_SP = 75;

spiral.dGamma  = 42576000.0; 
spiral.dGammaRad = 267512897.638;

spiral.GradDwellTime = 10; % from gcommon.h
%gSpiral_setParams(paramLongROSamples,paramLongSpiralType, 
%  		paramLongSpGradAmp,paramLongSpSlewRate,paramLongInterleaves,0,
%  		paramLongSpiralBW,paramDoubleVD,paramDoubleAccR,paramLongZCToFlat,nIceProgram);

spiral.paramLongROSamples = MrProt.spcl.g3.RO_samples ;
spiral.paramLongSpiralType = MrProt.spcl.g3.Spiral_Type ;
spiral.paramLongSpGradAmp = MrProt.spcl.g3.Spiral_Peak_Grad ;
spiral.paramLongSpSlewRate = MrProt.spcl.g3.Spiral_Slew_Rate ;
spiral.paramLongInterleaves = MrProt.spcl.g3.Interleaves ;
spiral.paramLongZCToFlat = MrProt.spcl.g3.Flat_First_ZC ;
spiral.spBW = MrProt.spcl.g3.Spiral_BW ;
spiral.VD = MrProt.spcl.g3.VD ;
spiral.AccR = MrProt.spcl.g3.Undersampling_fac ;
spiral.paramLongShiftsToMeasure = 0;
spiral.BASE_RES = BASE_RES;

spiral.paramDoubleEchoShift = MrProt.spcl.g4.CAIPI_shift_mm;
spiral.paramDoubleCAIPIPeriod = MrProt.spcl.g4.CAIPI_period_us;	
spiral.paramDoubleCAIPIDelay = MrProt.spcl.g4.CAIPI_delay_us;

if ((MrProt.spcl.g4.SMS_factor > 1) && (MrProt.spcl.g4.CAIPI_shift_mm > 0))
    spiral.UseCAIPI = true;
else
    spiral.UseCAIPI = false;
end
% for sosp
spiral.FOVz = MrProt.sliceSeries.aFront.thickness*0.001;
spiral.nsosp_segments = 1;
spiral.nz = MrProt.private.l_additionalslice+MrProt.sliceGroupList(1).sliceperslab;
spiral.nsosp_rz = spiral.nz/MrProt.spcl.private.slicepershot;%

 spiral.single_shot_mode = MrProt.spcl.single_shot_mode;
 spiral.accr_kz = MrProt.spcl.accr_kz;
 spiral.rotate_mode_ax = MrProt.spcl.private.rotate_mode_ax;% or 'golden_angle'
 spiral.increm_mode_sl = MrProt.spcl.private.increm_mode_sl;% or 'binomial' or 'arbit'
 spiral.nsosp_rz = MrProt.sliceGroupList(1).sliceperslab;