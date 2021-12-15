function Vec1 = VecMul(Vec, nSamples, c) 
if isfield(c,'real')
    for i=1:nSamples
        Vec1(i)=cmul(Vec(i),c);
    end
else
    for i=1:nSamples
        Vec1(i)=cmuld(Vec(i),c);
    end
end

end

