function InfSat = offcen_preset(MrProt,InfSat)

lSlicesToMeasure = MrProt.sliceSeries.size;       %// get number of slices

dGamma = InfSat.ifs_RF.LarmorConst;
if ~isfield(InfSat,'ifs_dAmpGSel')
    InfSat.ifs_dAmpGSel = InfSat.REF_InfSat_GSEL * 100.0 /InfSat.m_dInfSatThk;
end
ifs_dAmpGSel = InfSat.ifs_dAmpGSel;
ifs_dX_offcenter = 0;
dRFAsymmetry = 0.5;

if ( mod(lSlicesToMeasure , 2) == 1 )
    ifs_dX_offcenter = MrProt.sliceGroup(1).pslc(floor(lSlicesToMeasure/2 - 0.5+1)).SliceShift - ...
        (InfSat.m_dIRSlabThk + InfSat.m_dInfSatThk)/2.0 - InfSat.m_dInfSatGap;
else
    ifs_dX_offcenter = (MrProt.sliceGroup(1).pslc(floor(lSlicesToMeasure/2 )).SliceShift + ...
        MrProt.sliceGroup(1).pslc(floor(lSlicesToMeasure/2+1)).SliceShift )/2.0 - ...
        (InfSat.m_dIRSlabThk + InfSat.m_dInfSatThk)/2.0 - InfSat.m_dInfSatGap;
end

ifs_lFrequency = floor(dGamma * ifs_dX_offcenter * ifs_dAmpGSel + 0.5);
ifs_dPhase = double(ifs_lFrequency * 360. / (1.0E6) * InfSat.ifs_RF.Duration * dRFAsymmetry);

InfSat.ifs_RFSet.Frequency = ifs_lFrequency;
InfSat.ifs_RFSet.Phase = 180;%ifs_dPhase;
InfSat.ifs_RFNeg.Frequency = 0.0;
end
