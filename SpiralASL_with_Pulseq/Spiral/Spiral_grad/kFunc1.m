function temp = kFunc1(lambda, w, t, Tes, alpha, Ts2a, Tea, tautrans)


taut = tau1(t,Tes,alpha,Ts2a,Tea,tautrans);
Mag=lambda * pow(taut,alpha);
temp.real = Mag*cos(w*taut);
temp.imag = Mag*sin(w*taut);

end
