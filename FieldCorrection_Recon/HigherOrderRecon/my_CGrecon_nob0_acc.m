function [de,out] = my_CGrecon_nob0_acc(u,k,s,p,maxit,beta)
% u: k-space data (MRI), with the size of [nrRO x nrPE, nrCoils]
% k: field data (k-coefficients), with the size of [9, nrRO x nrPE] (field
% terms to the 2nd order) or [4, nrRO x nrPE] (1st order)
% s: sensitivity maps, with the size of [ xres x yres, nrCoils];
% p: position vectors, with the size of [3,xres x yres]
% [xres  yres] -> matrix size in image domain 
% maxit: number of iterations you want to run
% beta: serves for trading off noise and artefact (Wilm BJ, 2011 MRM), we use 0.000001 here,
% but this parameter is not very useful, may remove in future.
nCoil = size(u,2);
b = 0;

tstart = tic;
EH = EHmtx(k,p);
E = Emtx(k,p);
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

