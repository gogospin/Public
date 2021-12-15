function SupSat = offcen_preset(MrProt,SupSat)

lSlicesToMeasure = MrProt.sliceSeries.size;       %// get number of slices

dGamma = SupSat.sus_RF.LarmorConst;
if ~isfield(SupSat,'sus_dAmpGSel')
    SupSat.sus_dAmpGSel = REF_SUPSAT_GSEL * 100.0 /SupSat.m_dSupSatThk;
end
sus_dAmpGSel = SupSat.sus_dAmpGSel;
sus_dX_offcenter = 0;
dRFAsymmetry = 0.5;

if ( mod(lSlicesToMeasure , 2) == 1 )
    sus_dX_offcenter = MrProt.sliceGroup(1).pslc(floor(lSlicesToMeasure/2 - 0.5+1)).SliceShift + ...
        (SupSat.m_dIRSlabThk + SupSat.m_dSupSatThk)/2.0 + SupSat.m_dSupSatGap;
else
    sus_dX_offcenter = (MrProt.sliceGroup(1).pslc(floor(lSlicesToMeasure/2 - 1+1)).SliceShift + ...
        MrProt.sliceGroup(1).pslc(floor(lSlicesToMeasure/2+1)).SliceShift )/2.0 + ...
        (SupSat.m_dIRSlabThk + SupSat.m_dSupSatThk)/2.0 + SupSat.m_dSupSatGap;
end

sus_lFrequency = floor(dGamma * sus_dX_offcenter * sus_dAmpGSel + 0.5);
sus_dPhase = double(sus_lFrequency * 360. / (1.0E6) * SupSat.sus_RF.Duration * dRFAsymmetry);

SupSat.sus_RFSet.Frequency = sus_lFrequency;
SupSat.sus_RFSet.Phase= 180;%sus_dPhase;
SupSat.sus_RFNeg.Frequency = 0.0;
end
