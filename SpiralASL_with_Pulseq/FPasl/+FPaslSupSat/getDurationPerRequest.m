function m_lSBBDurationPerRequest_us = getDurationPerRequest(SupSat)


m_lSBBDurationPerRequest_us = ...
    SupSat.sus_GPSSel.RampUpTime + SupSat.sus_RF.Duration + SupSat.sus_GPSSel.RampDownTime + SupSat.sus_GPSpoil.TotalTime;
