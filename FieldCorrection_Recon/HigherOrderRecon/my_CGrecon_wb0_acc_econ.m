function [de,out] = my_CGrecon_wb0_acc_econ(u,k,s,p,b0,t,thresh)
nCoil = size(u,2);
b = 0;

%EH = EHmtx_wb0_econ(k,p,b0,t);
E = Emtx_wb0_econ(k,p,b0,t);
s=permute(s,[2,1]);
for isens = 1:nCoil
    srep((isens-1)*size(u,1)+1:isens*size(u,1),:) = repmat(s(isens,:),size(u,1),1);
end

% method1: svd
[uu,ss,vv] = svd(repmat(E,nCoil,1).*srep,'econ');
idxss = find(diag(ss)<ss(1)/thresh);
idxss = idxss(1)-1;
out = vv(:,1:idxss)*inv(ss(1:idxss,1:idxss))*uu(:,1:idxss)';
out = out*u(:);
% method 2: pinv

% temp= repmat(E,nCoil,1).*srep;
% out = (temp'*temp)\temp'*u(:);
% de=0;
de = diag(ss);

end
