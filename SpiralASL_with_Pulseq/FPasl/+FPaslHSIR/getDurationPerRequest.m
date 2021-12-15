function m_lSBBDurationPerRequest_us = getDurationPerRequest(HSIR, MrProt)

m_lSBBDurationPerRequest_us = ...
    HSIR.m_sGSS_ARB.Duration + HSIR.ss_GPSpoil.duration*1e6 + 260 * MrProt.sliceSeries.size; 

end
