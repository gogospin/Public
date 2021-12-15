%% calculate the timing of sosp


if MrProt.spcl.g4.SMS_factor == 1
        sliceloop = length(MrProt.sliceGroup.pslc);
    elseif MrProt.spcl.g4.SMS_factor > 1
        sliceloop = length(SMS_slice_reorder(MrProt));
    end
    spiral.duration_spi2D = sliceloop*1e6*(RFExcite.exc_GPSSel.duration+RFExcite.refoc.duration) ...
        +sliceloop*400+sliceloop*spiral.grad_x.t(end)*1e6;
    



% calculate timing tr_fill
tr_fill = MrProt.tr - 1e-3*(PASL.lSBBPaslDurationPerRequest_us + ...
    spiral.duration_spi2D);

tr_fill_delay = mr.makeDelay(tr_fill*1e-3);
tr_fill_skope = 3;
tr_fill_skope_delay = mr.makeDelay(tr_fill_skope);

