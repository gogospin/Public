%% make spoil gradient
spoil_area = 45000*larmor_freq*1e-6; %mT/m*us*Hz/mT*1e-6 -> 1/m
spoil_dur = (1860+270)*1e-6; %s
spoil_rise = 270*1e-6;
spoil_flat = (1860-270)*1e-6;
spoil_amp = 24.19*larmor_freq;

spiral.spoil_duration = spoil_dur*1e6;

spiral.gxSpoil = mr.makeTrapezoid('x',lims,'FlatTime',spoil_flat,'RiseTime',spoil_rise,'amplitude',spoil_amp);
spiral.gySpoil = mr.makeTrapezoid('y',lims,'FlatTime',spoil_flat,'RiseTime',spoil_rise,'amplitude',spoil_amp);
spiral.gzSpoil = mr.makeTrapezoid('z',lims,'FlatTime',spoil_flat,'RiseTime',spoil_rise,'amplitude',spoil_amp);