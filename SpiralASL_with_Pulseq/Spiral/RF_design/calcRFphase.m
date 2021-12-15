function RF = calcRFphase(RF,MrProt)
%CALCRFPHASE : calculate sRFNCOSet.prepSet, sRFNCONeg.prepNeg
% calculation from IDEA manual D2-33
dGamma = 42.5756;
if MrProt.spcl.g4.SMS_factor == 1
    for i =1:length(MrProt.sliceGroup.pslc)
        if strcmp(MrProt.sliceSeries.Slices, 'Ascending')
           delta_x(i) = -MrProt.sliceGroup.pslc(i).SliceShift;
        elseif strcmp(MrProt.sliceSeries.Slices, 'Descending')
           delta_x(i) = -MrProt.sliceGroup.pslc(length(MrProt.sliceGroup.pslc)-i+1).SliceShift; 
        end
    end
elseif MrProt.spcl.g4.SMS_factor > 1
    lNewSliceIndex = SMS_slice_reorder(MrProt);
    for i =1:length(lNewSliceIndex)
        delta_x(i) = -MrProt.sliceGroup.pslc(lNewSliceIndex(i)).SliceShift;
    end
end

if isfield(RF,'GPSSel_Amplitude')
    Gslice = RF.GPSSel_Amplitude;
    RF.GSAmplitude = RF.GPSSel_Amplitude;
elseif isfield(RF,'GSAmplitude')
    Gslice = RF.GSAmplitude;
end

Tpulse = RF.Duration;
Asym = RF.Asym;
InitPhase = RF.InitialPhase;

for i = 1:length(delta_x)
    f = dGamma*delta_x(i)*Gslice+0.5;
    RF.freqOffset(i) = f;
    RF.PhaseSet(i) = -(f*360/1e6*Tpulse*Asym)+InitPhase;
    RF.PhaseNeg(i) = -(f*360/1e6*Tpulse*(1-Asym))-InitPhase;
end

end

