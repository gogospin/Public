function  PreSat = offcen_preset(MrProt, PreSat)


lSlicesToMeasure = MrProt.sliceSeries.size;      %// get number of slices
dGamma = PreSat.pre_RF.LarmorConst;
pre_dAmpGSel = PreSat.REF_PRESAT_GSEL * 100. / PreSat.m_dIRSlabThk;
if ~isfield(PreSat,'Asymmetry')
    PreSat.Asymmetry = 0.5;
end

pre_dX_offcenter = MrProt.sliceGroup.pslc(1).SliceShift+15;
warning("PreSat offcenter freq set only implemented for Acending series");

pre_lFrequency = dGamma * pre_dX_offcenter * pre_dAmpGSel + 0.5;


PreSat.pre_RFSet.Frequency = pre_lFrequency;
PreSat.pre_RFSet.Phase = 180;
PreSat.pre_RFNeg.Frequency = 0.0;
