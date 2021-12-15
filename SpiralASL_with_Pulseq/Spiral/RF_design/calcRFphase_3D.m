function RF = calcRFphase_3D(RF,MrProt)
%CALCRFPHASE : calculate sRFNCOSet.prepSet, sRFNCONeg.prepNeg
% calculation from IDEA manual D2-33
% rf spoiler calculated from multiple simulations from Dimo's code
dGamma = 42.5756;
%if MrProt.spcl.g4.SMS_factor == 1
    %for i =1:length(MrProt.sliceGroup.pslc)
    
        delta_x = -MrProt.sliceGroupList(1).shift;
        %MrProt.sliceGroup.pslc(i).SliceShift;
        %delta_x = MrProt.sliceGroupList(1).shift;
    %end
% elseif MrProt.spcl.g4.SMS_factor > 1
%     lNewSliceIndex = SMS_slice_reorder(MrProt);
%     for i =1:length(lNewSliceIndex)
%         delta_x(i) = MrProt.sliceGroup.pslc(lNewSliceIndex(i)).SliceShift;
%     end
% end
Gslice = 0;
if isfield(RF,'GPSSel_Amplitude')
    Gslice = RF.GPSSel_Amplitude;
elseif isfield(RF,'GSAmplitude')
    Gslice = RF.GSAmplitude;
end

Tpulse = RF.Duration;
Asym = RF.Asym;
InitPhase = RF.InitialPhase;

% I don't know how this works for 3D; I kinda know now
f = floor(dGamma*delta_x*Gslice+0.5);
    RF.freqOffset = f;
    RF.InitPhaseSet = -(f*360/1e6*Tpulse*Asym)+InitPhase;
    RF.InitPhaseNeg = -(f*360/1e6*Tpulse*(1-Asym))-InitPhase;
for i = 1:RF.nrf %MrProt.private.l_additionalslice+MrProt.sliceGroupList(1).sliceperslab
    RF.PhaseSet(i) = RF.InitPhaseSet+25*i*i+175*i+300;
    RF.PhaseNeg(i) = RF.InitPhaseNeg+25*i*i+175*i+300;
end

end

