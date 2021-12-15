function firstPart = tau2( t, Tes, alpha, Ts2a, Tea,  tautrans)

firstPart = 0;
if(t>=0 && t<=Tes)
    firstPart=pow(t/Tes,1/(alpha/2+1));
end

end