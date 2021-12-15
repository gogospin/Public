%% PASL parameters:

% hard parameters (from pasl.h
PASL_TI1 = 700;                %% ms
PASL_TI1Stop = 1600;           %% ms
PASL_TI2 = 1800;               %% ms
PASL_NPRE_SATS = 2;
PASL_PREP_SCAN_TIME = 4000;    %% ms

FAIR_THICK_INV = 100;          %% mm
FAIR_THICK_CON = 300;          %% mm
FAIR_INFSAT_GAP = 0;           %% mm
FAIR_INFSAT_THK = 100;         %% mm

% hard parameters from PaslFair.h
% default values for initialization - are overwritten by UI input
VOL_THICK_INV = 100;  % mm
 VOL_THICK_CON = 300;    % mm
 VOL_THICK_INFSAT = 100;    % mm
 VOL_GAP_INFSAT = 0;    % mm
 VOL_GAP_SUPSAT = 0;    % mm

if isfield(MrProt,'ti')
    if length(MrProt.ti) == 2
        PASL_TI1 = MrProt.ti(1);
        PASL_TI2 = MrProt.ti(2);
    elseif length(MrProt.ti) == 1
        PASL_TI1 = MrProt.ti(1);
    end
end

% FIXED (!) RF parameters for inversion and saturation pulses
%    note: RF duration must only be changed in steps of gradient raster time
 RF_DURATION = 15360;
 RF_FLIPANGLE = 1260;

 RF_INFSAT_DURATION = 2560;   % us % Renzo probiert
 RF_PRESAT_DURATION = 2560;
 RF_SUPSAT_DURATION = 2560;
 RF_INFSAT_FLIPANGLE = 90;     % deg
 RF_PRESAT_FLIPANGLE = 90;
 RF_SUPSAT_FLIPANGLE = 90;

 VOL_THICK_SUPSAT = 40;     % mm

 SUPSAT_TOTAL = 1;       % total number of superior sat pulses
 SUPSAT_DELAY1 = 0;       % [us] 

 INFSAT_TOTAL = 4;       % total number of inferior sat pulses
 INFSAT_DELAY1 = 0;       % [us] 
% Renzo: InfSat_delay ist der delay zischen den SatPulsen, wieso er so gross ist weiss auch keiner!
% INFSAT_DELAY2         = 300000  % [us]
% INFSAT_DELAY3         = 300000  % [us]
 INFSAT_DELAY2 = 30000; % [us]
 INFSAT_DELAY3 = 30000;  % [us]
 INFSAT_DELAY4 = 0;       % [us] 

 REF_HSIR_GSEL = 0.7;     % reference gradient amplitude (mT/m) for  = 100;mm inversion
% REF_HSIR_GSEL			   8.0     % reference gradient amplitude (mT/m) for  = 100;mm inversion
 REF_INFSAT_GSEL = 0.716;   % reference gradient amplitude (mT/m) for  = 100;mm inferior saturation
 REF_SUPSAT_GSEL = 0.326;   % reference gradient amplitude (mT/m) for  = 100;mm superior saturation
 REF_PRESAT_GSEL = 0.716;   % reference gradient amplitude (mT/m) for  = 100;mm vol saturation

% soft parameters from Pasl.h.cpp
 PASL.m_dTI1 = PASL_TI1;
 PASL.m_dTI1Stop =PASL_TI1Stop;
  PASL.m_dTI2 = PASL_TI2;
  PASL.m_lNPreSats = PASL_NPRE_SATS;
  PASL.m_dASLSlabThickness = 100.0;
  PASL.m_dASLSlabOffCentre = -77.5;
  PASL.m_dASLSlabOffset = -76.25;
  PASL.m_dFairIRSlabThk = FAIR_THICK_INV;
  PASL.m_dFairIRSlabThkNS = FAIR_THICK_CON;
  PASL.m_dFairInfSatThk = FAIR_INFSAT_THK;
  PASL.m_dFairInfSatGap = FAIR_INFSAT_GAP;
  %m_lQuipssIIPulseType(ASLIR::WIP_IRBASSIa;          // ASLIR::WIP_IRBASSIa BPaslExciteType::RF_ARB_BASSI
  PASL.m_lTimeToImageSliceExc_us = 0;                % time between end of this SBB and excitation of the first image slice ...
                                     % this variable defines at which side of the image stack we put the inversion slab...
                                     % note: if you reverse this, you may want a DESCENDING slice order....
  PASL.m_dShiftDistalSign = 1;
  %m_eASLMode          (SEQ::ASL_PICOREQ2TIPS),
  %m_eLabelState       (ASLSTATE::ASL_LABEL),
  %m_dFlowLimit        (CRUSHGRAD_MAX_VELOCITY_ENC),
  PASL.m_dPrepScanTime = PASL_PREP_SCAN_TIME;  %% [ms]
  PASL.m_lPreSatTime_us = 0;                    %% time of presaturation pulses
  PASL.m_lSupSatTime_us = 0;                    %% time of superior sats
  PASL.m_lInfSatTime_us = 0;                    %% time of inferior sats
  PASL.m_lTI1FillTime_us = 0;                   %% fill time between inversion and periodic sats
  PASL.m_lTI2FillTime_us = 0;                   %% fill time after  sats
  PASL.m_lSATNumber = -1; 
  PASL.m_dFatSatFlipangle = MrProt.spcl.g5.FatSat_flip_angle;
  
  
  slab_thickness = (MrProt.sliceSeries.size)*MrProt.sliceGroupList(1).thickness *MrProt.sliceGroupList(1).distanceFactor + ...
    MrProt.sliceSeries.size*MrProt.sliceGroupList(1).thickness ;

  PASL.m_dASLSlabOffCentre = -(MrProt.satList.thickness/2 + MrProt.satList.gap+slab_thickness/2)+MrProt.sliceGroupList(1).shift;
  PASL.m_dASLSlabThickness = MrProt.satList.thickness;
  PASL.m_dSliceSlabThickness = slab_thickness;
  PASL.dASLSlabEdge = MrProt.satList.gap+PASL.m_dSliceSlabThickness/2;
  PASL.m_dFairIRSlabThk	= 2*PASL.dASLSlabEdge;
  PASL.m_dFairIRSlabThkNS = PASL.m_dFairIRSlabThk + 2 * PASL.m_dASLSlabThickness;  
  
  
  
  