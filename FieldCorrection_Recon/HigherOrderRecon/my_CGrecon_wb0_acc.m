function [de,out] = my_CGrecon_wb0_acc(u,k,s,p,b0,t,maxit,beta)
nCoil = size(u,2);
b = 0;

tstart = tic;
EH = EHmtx_wb0(k,p,b0,t);
E = Emtx_wb0(k,p,b0,t);
A = @(z) reshape(E*(s.*repmat(z,[1 nCoil])),size(E,1)*nCoil,1);
AH = @(z) sum(conj(s).*(EH*reshape(z,[size(EH,2),nCoil])),2);

ps = AH(u);

tt = toc(tstart);
%fprintf('initialization takes %d seconds\n',tt);
%ps = squeeze(sum(pl,2));
delta = ps;
%pt= squeeze(sum(pl,2));
de=zeros(maxit,1);
pt = ps;

for nit = 1:maxit
    
    %tstart = tic;
    qt = A(ps);
    qt = reshape(qt,[],nCoil);
    
%     tt = toc(tstart);
%     fprintf('iteration %d E takes %d seconds\n',nit,tt);
%     clear tt tstart
    
    qs = [];
   % tstart = tic;
    qs = AH(qt);
%     tt = toc(tstart);
%     fprintf('iteration %d EH takes %d seconds\n',nit,tt);
%     clear tt tstart
    
    %qs = squeeze(sum(q,2));
    qs = qs+beta*ps;
    b = b+((delta'*delta)/(ps'*qs))*ps;
    delta_old = delta;
    delta = delta-((delta'*delta)/(ps'*qs))*qs;
    ps = delta+((delta'*delta)/(delta_old'*delta_old))*ps;
    out(:,nit) = b;
    de(nit) = (delta'*delta)/(pt'*pt);
    
end
