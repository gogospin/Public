function temp = kFunc1_vec(lambda, w, t, Tes, alpha, Ts2a, Tea, tautrans)


% original : taut = tau1(t,Tes,alpha,Ts2a,Tea,tautrans);
firstPart = 0;
if(t>=0 && Ts2a>=t)
    firstPart=pow(t/Tes,1/(alpha/2+1));
end
secondPart = 0;
if(t>Ts2a && t<=Tea && Tes>=Ts2a)
    secondPart=pow((t-Ts2a)/Tea + pow(tautrans,alpha+1),1/(alpha+1));
end
taut = firstPart + secondPart;

Mag=lambda * pow(taut,alpha);
temp.real = Mag*cos(w*taut);
temp.imag = Mag*sin(w*taut);

end
