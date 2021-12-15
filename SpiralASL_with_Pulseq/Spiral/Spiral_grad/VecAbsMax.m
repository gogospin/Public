function Mx = VecAbsMax(Vec, nSamples)
Mx=-999999999999999.9;
for i=1:nSamples
    Mx=max(Mx,abs(Vec(i)));
end

end