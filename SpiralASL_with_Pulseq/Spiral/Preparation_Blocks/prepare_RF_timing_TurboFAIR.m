%% ##############   calculate RF pulse timing (in us)     #################
if MrProt.perfusion == 1
%     PASL.m_lTI1FillTime_us = PASL.m_dTI1*1000 - FPaslHSIR.getDurationPerRequest(HSIR,MrProt) + ...
%         FPaslHSIR.getTimeToCentreRF_us(HSIR);
    
    
    if PASL.m_lTI1FillTime_us < 0
        error ('timing error m_lTI1FillTime_us');
    end
    
    %PASL.m_lPreSatTime_us = FPaslPreSat.getDurationPerRequest(PreSat) * PASL.m_lNPreSats;
    
    %PASL.m_lSupSatTime_us = FPaslSupSat.getDurationPerRequest(SupSat) * SUPSAT_TOTAL + SUPSAT_DELAY1;
    
    if MrProt.ti(1)*1000 - MrProt.tr*1000 - (HSIR.m_sGSS_ARB.Duration + HSIR.ss_GPSpoil.duration*1e6) > 300000
        tr_orig = MrProt.tr;
        MrProt.tr = MrProt.ti(1)*1000 - (HSIR.m_sGSS_ARB.Duration + HSIR.ss_GPSpoil.duration*1e6) - 300000;

        MrProt.tr = MrProt.tr/1000;
        warnstr = strcat('TR changes from',num2str(tr_orig),' ms to ',num2str(MrProt.tr));
        warning(warnstr);
        clear tr_orig warnstr
    end
    if MrProt.ti(1)*1000 - MrProt.tr*1000 < FPaslHSIR.getDurationPerRequest(HSIR,MrProt)
        tr_orig = MrProt.tr;
        MrProt.tr = MrProt.ti(1)*1000 - FPaslHSIR.getDurationPerRequest(HSIR,MrProt) - 300000;
        MrProt.tr = MrProt.tr/1000;
        warnstr = strcat('TR changes from',num2str(tr_orig),' ms to ',num2str(MrProt.tr));
        warning(warnstr);
        clear tr_orig warnstr
    end
    
    
    % PASL.m_lTI1FillTime_us: between HSIR block and next Inversion block
    PASL.m_lTI1FillTime_us = MrProt.tr*1000 - FPaslHSIR.getDurationPerRequest(HSIR,MrProt) + ...
        FPaslHSIR.getTimeToCentreRF_us(HSIR);
    % TI2 fill time, TimeToImageSliceExc is subtracted from TI2
    PASL.m_lTI2FillTime_us = PASL.m_dTI1*1000 - PASL.m_lTimeToImageSliceExc_us ...
        - FPaslHSIR.getDurationPerRequest(HSIR,MrProt) + FPaslHSIR.getTimeToCentreRF_us(HSIR);% ...
        %- CSatFat.TotalTime - SpoilGrad.TotalTime;
    if( PASL.m_lTI2FillTime_us < 0 )
        error ('timing error m_lTI1FillTime_us');
    end
    
    PASL.lSBBPaslDurationPerRequest_us  = PASL.m_lPreSatTime_us + PASL.m_lSupSatTime_us + ...
        FPaslHSIR.getDurationPerRequest(HSIR,MrProt) + PASL.m_lTI1FillTime_us +  ...
        CSatFat.TotalTime + SpoilGrad.TotalTime();% + PASL.m_lTI2FillTime_us
elseif MrProt.perfusion == 0
    PASL.lSBBPaslDurationPerRequest_us  = CSatFat.TotalTime + SpoilGrad.TotalTime();
end