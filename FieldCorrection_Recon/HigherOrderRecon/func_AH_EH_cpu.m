function out = func_AH_EH_cpu(k,p,s,z,jump_rate)
%EH = EHmtx(k,p);

%AH = @(z) sum(conj(s).*(EH*reshape(z,[size(EH,2),nCoil])),2);

gamma = 267.513; % gyromagnetic ratio for proton in rad/s/T
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

%tic;E=single(exp(-sqrt(-1)*permute(hs,[2,1])*k));toc;
EH_temp2 = 0;

for ii = 1:jump_rate:size(k,2)
    %tic;
    %parfor jj = 1:size(hs,2)
    EH_temp = exp(-sqrt(-1)*permute(hs,[2,1])*k(:,ii:ii+jump_rate-1));
    
    %end
    %tooc = toc;
    EH_temp2 = EH_temp2 + (EH_temp*z(ii:ii+jump_rate-1,:));
    
    clear EH_temp;
end

out = sum(conj(s).*EH_temp2,2); 
clear EH_temp2;
