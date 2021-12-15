%% ##############   design spiral grad        #############################

%[spiral,MrProt,SeqLim] = sprial_prep_sosp(spiral,MrProt,SeqLim);
RFExcite.RF.nrf = spiral.nsosp_segments*spiral.nsosp_rz*MrProt.dynamics;
RFExcite.RF = calcRFphase_3D(RFExcite.RF,MrProt);


count_rf = 1;

for idyn = 1:MrProt.dynamics
    for ll = 1:spiral.nsosp_segments
        for i = 1:spiral.nsosp_rz
            RFExcite.RF.rf_out(count_rf) = rf0;
            RFExcite.RF.rf_out(count_rf).freqOffset = RFExcite.RF.freqOffset;
            %phase_offset_input = (RFExcite.RF.PhaseSet(count_rf)+4336)/9.1189; %strange inaccordance between what I put in here and how pulseq compiles and generates the NCO phase from poet simulation
            phase_offset_input = (RFExcite.RF.PhaseSet(count_rf));
            RFExcite.RF.rf_out(count_rf).phaseOffset = phase_offset_input;% mod(phase_offset_input,360)*pi/180;
            %phase_offset_input_adc = (RFExcite.RF.PhaseSet(count_rf)-RFExcite.RF.InitPhaseSet+4336)/9.1189;
            phase_offset_input_adc = RFExcite.RF.PhaseSet(count_rf)-RFExcite.RF.InitPhaseSet;
            %phase_offset_input_adc = mod(phase_offset_input_adc,360)*pi/180;
            for j = 1:spiral.nz/spiral.nsosp_rz
                 adcDwell = 1e-5*(spiral.adc_shot(i).end(j)-spiral.adc_shot(i).start(j)+1)/spiral.paramLongROSamples;
                 spiral.adc.seg(ll).shot(i).kslice(j).dynamic(idyn).adc = mr.makeAdc(spiral.paramLongROSamples ,'Dwell',adcDwell, ...
                    'Delay',(spiral.adc_shot(i).start(j)-1)*1e-5,'phaseOffset',phase_offset_input_adc);
            end
            count_rf = count_rf+1;
        end
    end
end