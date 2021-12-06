function [de,out] = my_CGrecon_nob0_acc(u,k,s,p,b0,t,jump_rate,maxit,beta)
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
fbar = waitbar(0,'Please wait for the initialization');

tstart = tic;
% EH = EHmtx(k,p);
% E = Emtx(k,p);
% A = @(z) reshape(E*(s.*repmat(z,[1 nCoil])),size(E,1)*nCoil,1);
% AH = @(z) sum(conj(s).*(EH*reshape(z,[size(EH,2),nCoil])),2);

ps = func_AH_EH_wb0(k,p,s,u,b0,t,jump_rate);
tt = toc(tstart);
fprintf('initialization takes %d seconds\n',tt);
close(fbar);

%ps = squeeze(sum(pl,2));
delta = ps;
%pt= squeeze(sum(pl,2));
de=gpuArray(zeros(maxit,1));
pt = ps;

fbar = waitbar(0,'CG SENSE in progress');
for nit = 1:maxit
    
    %tstart = tic;
    
    %%qt = A(ps);
    qt = func_A_E_wb0(k,p,s,ps,b0,t,nCoil,jump_rate);
    qt = reshape(qt,[],nCoil);
    

    qs = [];
   % tstart = tic;
    %%qs = AH(qt);
    qs = func_AH_EH_wb0(k,p,s,qt,b0,t,jump_rate);
%     tt = toc(tstart);
%     fprintf('iteration %d EH takes %d seconds\n',nit,tt);
%     clear tt tstart
    waitbar(nit/maxit,fbar,'CG SENSE in progress');
    %qs = squeeze(sum(q,2));
    qs = qs+beta*ps;
    b = b+((delta'*delta)/(ps'*qs))*ps;
    delta_old = delta;
    delta = delta-((delta'*delta)/(ps'*qs))*qs;
    ps = delta+((delta'*delta)/(delta_old'*delta_old))*ps;
    out(:,nit) = b;
    de(nit) = gpuArray((delta'*delta)/(pt'*pt));
  
end

close(fbar);
