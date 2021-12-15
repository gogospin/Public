function temp = cmuld( n1, d)

if ~isfield(n1,'real')
    temp = n1*d;
else
    temp.real = n1.real*d;
    temp.imag = n1.imag*d;
end


end

