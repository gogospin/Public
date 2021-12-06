function out = func_A_E_cpu(k,p,s,z,nCoil,jump_rate)
% use GPU and for loop to perform the calculation rather than demanding
% size
%E = Emtx(k,p);
%A = @(z) reshape(E*(s.*repmat(z,[1 nCoil])),size(E,1)*nCoil,1);

gamma = 267.513; % gyromagnetic ratio for proton in rad/us/T
field_strength = 6.98;


% bulding spherical harmonics
if size(k,1) == 9
    hs = (zeros(9,size(p,2)));
    hs(1,:) = 1;
    hs(2,:) = squeeze(p(1,:));
    hs(3,:) = squeeze(p(2,:));
    hs(4,:) = squeeze(p(3,:));
    hs(5,:) = squeeze(p(1,:)).*squeeze(p(2,:));
    hs(6,:) = squeeze(p(2,:)).*squeeze(p(3,:));
    hs(7,:) = 2*((squeeze(p(3,:))).^2)-(squeeze(p(1,:))).^2-(squeeze(p(2,:))).^2;
    hs(8,:) = squeeze(p(1,:)).*squeeze(p(3,:));
    hs(9,:) = (squeeze(p(1,:))).^2-(squeeze(p(2,:))).^2;
elseif size(k,1) == 4
    hs = (zeros(4,size(p,2)));
    hs(1,:) = 1;
    hs(2,:) = squeeze(p(1,:));
    hs(3,:) = squeeze(p(2,:));
    hs(4,:) = squeeze(p(3,:));
    elseif size(k,1) == 16
    hs = (zeros(16,size(p,2)));
    hs(1,:) = 1;
    hs(2,:) = squeeze(p(1,:));
    hs(3,:) = squeeze(p(2,:));
    hs(4,:) = squeeze(p(3,:));
    hs(5,:) = squeeze(p(1,:)).*squeeze(p(2,:));
    hs(6,:) = squeeze(p(2,:)).*squeeze(p(3,:));
    hs(7,:) = 2*((squeeze(p(3,:))).^2)-(squeeze(p(1,:))).^2-(squeeze(p(2,:))).^2;
    hs(8,:) = squeeze(p(1,:)).*squeeze(p(3,:));
    hs(9,:) = (squeeze(p(1,:))).^2-(squeeze(p(2,:))).^2;
    hs(10,:) = 3*squeeze(p(2,:)).*(squeeze(p(1,:))).^2-(squeeze(p(3,:))).^3;
    hs(11,:) = squeeze(p(1,:)).*squeeze(p(2,:)).*squeeze(p(3,:));
    hs(12,:) = (5*squeeze(p(3,:)).^2-(squeeze(p(1,:)).^2+ squeeze(p(2,:)).^2 + squeeze(p(3,:)).^2)).*squeeze(p(2,:));
    hs(13,:) = 5*squeeze(p(3,:)).^3 - 3*squeeze(p(3,:)).*(squeeze(p(1,:)).^2+squeeze(p(2,:)).^2+squeeze(p(3,:)).^2);
    hs(14,:) = (5*squeeze(p(3,:)).^2-(squeeze(p(1,:)).^2+ squeeze(p(2,:)).^2 + squeeze(p(3,:)).^2)).*squeeze(p(1,:));
    hs(15,:) = squeeze(p(1,:)).^2.*squeeze(p(3,:)) - squeeze(p(2,:)).^2.*squeeze(p(3,:));
    hs(16,:) = squeeze(p(1,:)).^3 - 3*squeeze(p(1,:)).*squeeze(p(2,:)).^2;
end

%tic;E=single(exp(sqrt(-1)*(permute(k,[2,1])*hs)));toc;
k=permute(k,[2 1]);

for ii = 1:jump_rate:size(k,1)
    %tic;
    %parfor jj = 1:size(hs,2)
        E_temp = exp(sqrt(-1)*k(ii:ii+jump_rate-1,:)*hs(:,:));
    %end
    %tooc = toc;
    out(ii:ii+jump_rate-1,:) = E_temp*(s.*repmat(z,[1 nCoil]));
    clear E_temp;
end

