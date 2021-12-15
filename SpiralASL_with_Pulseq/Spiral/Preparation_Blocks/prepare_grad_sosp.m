%% ##############   design spiral grad        #############################

[spiral,MrProt,SeqLim] = sprial_prep_sosp(spiral,MrProt,SeqLim);
RFExcite.RF.nrf = spiral.nsosp_segments*spiral.nsosp_rz*MrProt.dynamics;
RFExcite.RF = calcRFphase_3D(RFExcite.RF,MrProt);
spiral.total_grad_duration = 0;

%count_rf = 1;
for ll = 1:spiral.nsosp_segments
    for i = 1:spiral.nsosp_rz
        
        %         RFExcite.RF.rf_out(count_rf) = rf0;
        %         RFExcite.RF.rf_out(count_rf).freqOffset = RFExcite.RF.freqOffset;
        %         phase_offset_input = (RFExcite.RF.PhaseSet(count_rf)+4336)/9.1189; %strange inaccordance between what I put in here and how pulseq compiles and generates the NCO phase from poet simulation
        %         RFExcite.RF.rf_out(count_rf).phaseOffset = phase_offset_input/2/pi;
        %         count_rf = count_rf+1;
        if spiral.single_shot_mode == 1
            loopnr = 1;
        else
            loopnr = spiral.nz/spiral.nsosp_rz;
        end
        for j = 1:loopnr
            mylength = length(spiral.gx_seg(ll).shot(i).kslice(j).grad);
            spiral.grad_x.seg(ll).shot(i).kslice(j).grad = mr.makeArbitraryGrad('x',spiral.gx_seg(ll).shot(i).kslice(j).grad*larmor_freq,'system',lims_sp);
            spiral.grad_y.seg(ll).shot(i).kslice(j).grad = mr.makeArbitraryGrad('y',spiral.gy_seg(ll).shot(i).kslice(j).grad*larmor_freq,'system',lims_sp);
            spiral.grad_z.seg(ll).shot(i).kslice(j).grad = mr.makeArbitraryGrad('z',spiral.gz_seg(ll).shot(i).kslice(j).grad*larmor_freq,'system',lims_sp);
            spiral.total_grad_duration = spiral.total_grad_duration + mylength * 10;
        end
        spiral.grad_x.seg(ll).shot(i).rampdown = mr.makeArbitraryGrad('x',spiral.gx_seg(ll).shot(i).rampdown*larmor_freq,'system',lims_sp);
        spiral.grad_y.seg(ll).shot(i).rampdown = mr.makeArbitraryGrad('y',spiral.gy_seg(ll).shot(i).rampdown*larmor_freq,'system',lims_sp);
        spiral.grad_z.seg(ll).shot(i).rampdown = mr.makeArbitraryGrad('z',spiral.gz_seg(ll).shot(i).rampdown*larmor_freq,'system',lims_sp);
        spiral.total_grad_duration = spiral.total_grad_duration + spiral.RAMP_DOWN_POINTS * 10;
    end
end