function MaxS = GetMaxVecDiff(Vec, N)

MaxS=-100;
for i=1:N-1
    MaxS=max(MaxS,abs(Vec(1+i)-Vec(i)));
end