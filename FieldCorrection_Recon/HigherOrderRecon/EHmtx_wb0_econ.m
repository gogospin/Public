function E = EHmtx_wb0_econ(k,p,b0,t)
% E: Nrkpoints x Nrpositions
% input : 
%       k: 7 x Nrkpoints (for each x)
%       p: 3 x Nrpositions (for each x)
%       t: time in seconds
%       b0: b0map, arranged in the same way as sense map
%gamma = 267.513; % gyromagnetic ratio for proton in rad/s/T
%field_strength = 6.98;

% bulding spherical harmonics
if size(k,1) == 7
    hs = zeros(8,size(p,2));
    hs(1,:) = squeeze(p(2,:));
    hs(2,:) = squeeze(p(3,:));
    hs(3,:) = squeeze(p(1,:)).*squeeze(p(2,:));
    hs(4,:) = squeeze(p(2,:)).*squeeze(p(3,:));
    hs(5,:) = 2*((squeeze(p(3,:))).^2)-(squeeze(p(1,:))).^2-(squeeze(p(2,:))).^2;
    hs(6,:) = squeeze(p(1,:)).*squeeze(p(3,:));
    hs(7,:) = (squeeze(p(1,:))).^2-(squeeze(p(2,:))).^2;
    hs(8,:) = b0(:);
elseif size(k,1) == 4
    hs = zeros(3,size(p,2));
    hs(1,:) = squeeze(p(2,:));
    hs(2,:) = squeeze(p(3,:));
    hs(3,:) = b0(:);
elseif size(k,1) == 16
    hs = zeros(14,size(p,2));
    hs(1,:) = squeeze(p(2,:));
    hs(2,:) = squeeze(p(3,:));
    hs(3,:) = squeeze(p(1,:)).*squeeze(p(2,:));
    hs(4,:) = squeeze(p(2,:)).*squeeze(p(3,:));
    hs(5,:) = 2*((squeeze(p(3,:))).^2)-(squeeze(p(1,:))).^2-(squeeze(p(2,:))).^2;
    hs(6,:) = squeeze(p(1,:)).*squeeze(p(3,:));
    hs(7,:) = (squeeze(p(1,:))).^2-(squeeze(p(2,:))).^2;
    hs(8,:) = 3*squeeze(p(2,:)).*(squeeze(p(1,:))).^2-(squeeze(p(3,:))).^3;
    hs(9,:) = squeeze(p(1,:)).*squeeze(p(2,:)).*squeeze(p(3,:));
    hs(10,:) = (5*squeeze(p(3,:)).^2-(squeeze(p(1,:)).^2+ squeeze(p(2,:)).^2 + squeeze(p(3,:)).^2)).*squeeze(p(2,:));
    hs(11,:) = 5*squeeze(p(3,:)).^3 - 3*squeeze(p(3,:)).*(squeeze(p(1,:)).^2+squeeze(p(2,:)).^2+squeeze(p(3,:)).^2);
    hs(12,:) = (5*squeeze(p(3,:)).^2-(squeeze(p(1,:)).^2+ squeeze(p(2,:)).^2 + squeeze(p(3,:)).^2)).*squeeze(p(1,:));
    hs(13,:) = squeeze(p(1,:)).^2.*squeeze(p(3,:)) - squeeze(p(2,:)).^2.*squeeze(p(3,:));
    hs(14,:) = squeeze(p(1,:)).^3 - 3*squeeze(p(1,:)).*squeeze(p(2,:)).^2;
    hs(15,:) = b0(:);
end

%s=conj(reshape(s,length(s),1));
% t = reshape(t,1,length(t));
% t = t*gamma*field_strength;
% 
% hs = cat(1,ones(1,size(hs,2)),hs);
% if mod (size(k,2),length(t))~=0
%     error('wrong size of time variable and k traj')
%     return
% else
     trep = size(k,2)/length(t);
% end    
 t = reshape(t,1,[]);
 k = cat(1,k,repmat(t,1,trep));

%tic;E=single(zeros(size(k,2),size(p,2)));toc;
tic;E=single(exp(-sqrt(-1)*permute(hs,[2,1])*k));toc;



