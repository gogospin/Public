function m_lSBBDurationPerRequest_us = getDurationPerRequest(PreSat)

m_lSBBDurationPerRequest_us = ...
    PreSat.pre_GPSSel.RampUpTime + PreSat.pre_RF.Duration + PreSat.pre_GPSSel.RampDownTime + PreSat.pre_GPSpoil.TotalTime;
    
 
