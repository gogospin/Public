function Tend = VDgetTend(FOV, res, AccR, nInterleaves, alpha, MaxGradAmp, MaxSlewRate, spiral)

dGammaRad = spiral.dGammaRad;

lambda = .5/res; % in m^(-1)
n = 1/(1-pow((1-nInterleaves*AccR/(FOV*lambda)),(1/alpha)));
w = 2*pi()*n;

Tea = lambda*w/(dGammaRad*MaxGradAmp*(alpha+1)); %/ in s
Tes = sqrt(lambda*(w*w)/(MaxSlewRate*dGammaRad))/(alpha/2+1); % in s
Ts2a = pow(pow(Tes,((alpha+1)/(alpha/2+1)))*(alpha/2+1)/(Tea*(alpha+1)),(1+2/alpha)); % in s

if(Ts2a<Tes)
    Tend = Tea;
else
    Tend = Tes;
end

end