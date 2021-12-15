%% ##############   prepare spoil gradient    #############################
SpoilGrad.rampup = SpoilGrad.m_lRampTime*1e-6 ; 
SpoilGrad.rampdown = SpoilGrad.m_lRampTime*1e-6 ; 
SpoilGrad.plateau = (SpoilGrad.m_lSpoilDuration - SpoilGrad.m_lRampTime) *1e-6 ; 

SpoilGrad.duration = SpoilGrad.rampup + SpoilGrad.rampdown + SpoilGrad.plateau ;
SpoilGrad.amplitude = SpoilGrad.m_dSpoilAmplitude * larmor_freq ; %mT/m  
SpoilGrad.TotalTime = SpoilGrad.m_lSpoilDuration + SpoilGrad.m_lRampTime;

SpoilGrad.axix = 'y' ;
SpoilGrad.grad_y = mr.makeTrapezoid(SpoilGrad.axix,lims,'duration', SpoilGrad.duration, ...
    'flatTime',SpoilGrad.plateau,'riseTime',SpoilGrad.rampup, 'amplitude', SpoilGrad.amplitude);
SpoilGrad.axix = 'z' ;
SpoilGrad.grad_z = mr.makeTrapezoid(SpoilGrad.axix,lims,'duration', SpoilGrad.duration, ...
    'flatTime',SpoilGrad.plateau,'riseTime',SpoilGrad.rampup, 'amplitude', SpoilGrad.amplitude);