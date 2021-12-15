%% ##############   calculate RF pulse timing (in us)     #################
if MrProt.perfusion == 1
    PASL.m_lTI1FillTime_us = PASL.m_dTI1*1000 - FPaslHSIR.getDurationPerRequest(HSIR,MrProt) + ...
        FPaslHSIR.getTimeToCentreRF_us(HSIR) - FPaslInfSat.getTimeToCentreRF_us(InfSat);
    
    if PASL.m_lTI1FillTime_us < 0
        error ('timing error m_lTI1FillTime_us');
    end
    
    PASL.m_lPreSatTime_us = FPaslPreSat.getDurationPerRequest(PreSat) * PASL.m_lNPreSats;
    
    PASL.m_lSupSatTime_us = FPaslSupSat.getDurationPerRequest(SupSat) * SUPSAT_TOTAL + SUPSAT_DELAY1;
    
    PASL.m_lInfSatTime_us = FPaslInfSat.getDurationPerRequest(InfSat) * INFSAT_TOTAL ...
        + INFSAT_DELAY1 + INFSAT_DELAY2 + INFSAT_DELAY3 + INFSAT_DELAY4;
    
    % TI2 fill time, TimeToImageSliceExc is subtracted from TI2
    PASL.m_lTI2FillTime_us = PASL.m_dTI2*1000 - PASL.m_lTimeToImageSliceExc_us ...
        - FPaslHSIR.getDurationPerRequest(HSIR,MrProt) + FPaslHSIR.getTimeToCentreRF_us(HSIR) ...
        - PASL.m_lTI1FillTime_us - CSatFat.TotalTime ...
        - SpoilGrad.TotalTime- PASL.m_lInfSatTime_us;
    if( PASL.m_lTI2FillTime_us < 0 )
        error ('timing error m_lTI1FillTime_us');
    end
    
    PASL.lSBBPaslDurationPerRequest_us  = PASL.m_lPreSatTime_us + PASL.m_lSupSatTime_us + ...
        FPaslHSIR.getDurationPerRequest(HSIR,MrProt) + PASL.m_lTI1FillTime_us + PASL.m_lInfSatTime_us + ...
        PASL.m_lTI2FillTime_us + CSatFat.TotalTime + SpoilGrad.TotalTime();
elseif MrProt.perfusion == 0
    PASL.lSBBPaslDurationPerRequest_us  = CSatFat.TotalTime + SpoilGrad.TotalTime();
end