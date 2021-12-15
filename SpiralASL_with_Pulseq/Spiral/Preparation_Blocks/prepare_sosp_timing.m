%% calculate the timing of sosp
count_rf = 1;

spiral.delay_pre_spiral = RFExcite.postdelay;
spiral.total_rf_duration = 0;
spiral.total_delay_pre_spiral = 0;
spiral.total_spoil_duration = 0;

%for idyn = 1:MrProt.dynamics
    for ll = 1:spiral.nsosp_segments
        for i = 1:spiral.nsosp_rz
            count_rf = count_rf+1;
            spiral.total_rf_duration = spiral.total_rf_duration + 1e6*(RFExcite.exc_GPSSel.duration+RFExcite.refoc.duration);
            spiral.total_delay_pre_spiral = spiral.total_delay_pre_spiral + spiral.delay_pre_spiral;
            spiral.total_spoil_duration = spiral. spoil_duration + spiral.total_spoil_duration;
        end
    end
%end
spiral.duration_sosp = spiral.total_rf_duration + spiral.total_delay_pre_spiral + ...
    spiral.total_grad_duration + spiral.total_spoil_duration;

% calculate timing tr_fill
tr_fill = MrProt.tr - 1e-3*(PASL.lSBBPaslDurationPerRequest_us + ...
    spiral.duration_sosp);

tr_fill_delay = mr.makeDelay(tr_fill*1e-3);
tr_fill_skope = 3;
tr_fill_skope_delay = mr.makeDelay(tr_fill_skope);

