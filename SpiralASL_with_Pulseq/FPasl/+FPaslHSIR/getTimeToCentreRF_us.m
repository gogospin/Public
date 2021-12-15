function m_lSBBTimeToCentreRF_us = getTimeToCentreRF_us(HSIR)

m_lSBBTimeToCentreRF_us = HSIR.m_sGSS.RampUpTime + floor(HSIR.m_sSRF_FOCI_LBL.Asymmetry *HSIR. m_sSRF_FOCI_LBL.Duration + 0.5);
end
