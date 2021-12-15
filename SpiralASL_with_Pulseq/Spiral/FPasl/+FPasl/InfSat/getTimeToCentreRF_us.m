function m_lSBBTimeToCentreRF_us = getTimeToCentreRF_us(InfSat)

m_lSBBTimeToCentreRF_us = InfSat.ifs_GPSSel.RampUpTime + floor(InfSat.ifs_RF.Asymmetry * InfSat.ifs_RF.Duration + 0.5);
