%% PART 1: TEST kz sampling
spiral.accr_kz = 2; nz = 36; spiral.single_shot_mode = 1;FOVz = 0.048;spiral.increm_mode_sl = 'binomial';
%%
nz_eff = ceil(nz/spiral.accr_kz);
if mod(nz_eff,2) == 1
    nz_eff = nz_eff + 1;
end
clear kz_targ;
if spiral.single_shot_mode == 1
    if strcmp(spiral.increm_mode_sl,'linear') == 1
        kzfov_half = 1/(FOVz/nz)/2/1000; %rad/2pi/mm
        delta_kz = 1/FOVz/1000; %rad/2pi/mmm
        for i = 1:nz_eff
            kz_targ(i,1) = kzfov_half - delta_kz/2 - (i-1)*delta_kz;
        end
    elseif strcmp(spiral.increm_mode_sl,'arbit') == 1
        nz_eff_full = round(nz/3);
        % nz_eff_full: number of kz sampling points subjected to nyquist
        % sampling criteria
        if mod(nz_eff_full,2) == 1
            nz_eff_full = nz_eff_full+1;
        end
        nz_eff_und = nz_eff - nz_eff_full;
        % nz_eff_und: number of kz sampling points that are undersampled
        if mod(nz_eff_und,2) == 1
            nz_eff_und = nz_eff_und+1;
        end
        delta_kz_full = 1/FOVz/1000; %rad/2pi/mmm
        delta_kz_und = 1/(FOVz/nz)/1000 - delta_kz_full - delta_kz_full*(nz_eff_full-1);
        delta_kz_und = delta_kz_und/nz_eff_und;
        for i = 1:nz_eff_und/2+1
            kz_targ(i,1) = 1/(FOVz/nz)/2/1000 - delta_kz_full/2-(i-1)*delta_kz_und;
        end
        for i = 2:nz_eff_full/2
            kz_targ(i+nz_eff_und/2,1) = kz_targ(i+nz_eff_und/2-1,1)-delta_kz_full;
        end
        for i = 1:(nz_eff_full+nz_eff_und)/2
            kz_targ(i+(nz_eff_full+nz_eff_und)/2,1) = -kz_targ((nz_eff_full+nz_eff_und)/2-i+1,1);
        end
    elseif strcmp(spiral.increm_mode_sl, 'binomial')==1
        delta_kz_full = 1/FOVz/1000; %rad/2pi/mmm
        kzfov_half = 1/(FOVz/nz)/2/1000; %rad/2pi/mm
        % fitting the binomial curve, with the intercept temporaly set to
        % be delta_kz_full/3;
        % y = ax^2 + bx + delta_kz_full/3
        y1 = delta_kz_full/2;
        y2 = kzfov_half-delta_kz_full/2;
        x1 = 1;
        x2 = nz_eff/2;
        c = delta_kz_full/3;
        a = (y1*x2-y2*x1-c*x2+c*x1)/x1/x2/(x1-x2);
        b = (y1*x2*x2 - y2*x1*x1 - c*x2*x2 + c*x1*x1)/x1/x2/(x2-x1);
        for i = 1:nz_eff/2
            idx = nz_eff/2-i+1;
            kz_targ(i,1) = a*idx*idx+b*idx+c;
            kz_targ(i+nz_eff/2,1) = -(a*i*i+b*i+c);
        end
    end
end