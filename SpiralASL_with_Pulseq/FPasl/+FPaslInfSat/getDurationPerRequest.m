function m_lSBBDurationPerRequest_us = getDurationPerRequest(InfSat)

m_lSBBDurationPerRequest_us = ...
    InfSat.ifs_GPSSel.RampUpTime + InfSat.ifs_RF.Duration + InfSat.ifs_GPSSel.RampDownTime + InfSat.ifs_GPSpoil.TotalTime;
end
