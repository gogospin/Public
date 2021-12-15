function temp = cmul(n1, n2)

temp.real = n1.real*n2.real-n1.imag*n2.imag;
temp.imag = n1.real*n2.imag + n1.imag*n2.real;

end