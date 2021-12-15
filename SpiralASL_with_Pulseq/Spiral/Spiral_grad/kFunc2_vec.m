function temp = kFunc2_vec(lambda, w, t, Tes, alpha, Ts2a, Tea, tautrans)

% original: taut=tau2(t,Tes,alpha,Ts2a,Tea,tautrans);

taut = 0;
if(t>=0 && t<=Tes)
    taut = pow(t/Tes,1/(alpha/2+1));
end
Mag=lambda*  pow(taut,alpha);

temp.real = Mag*cos(w*taut);
temp.imag = Mag*sin(w*taut);

end