function E = EHmtx(k,p)
% E: Nrkpoints x Nrpositions
% input : 
%       k: 9 x Nrkpoints
%       p: 3 x Nrpositions
%       
gamma = 267.513; % gyromagnetic ratio for proton in rad/s/T
field_strength = 6.98;

% bulding spherical harmonics
if size(k,1) == 9
    hs = zeros(9,size(p,2));
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
    hs = zeros(4,size(p,2));
    hs(1,:) = 1;
    hs(2,:) = squeeze(p(1,:));
    hs(3,:) = squeeze(p(2,:));
    hs(4,:) = squeeze(p(3,:));
elseif size(k,1) == 16
    hs = zeros(16,size(p,2));
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

%s=conj(reshape(s,length(s),1));
% t = reshape(t,1,length(t));
% t = t*gamma*field_strength;
% 
% hs = cat(1,ones(1,size(hs,2)),hs);
% if mod (size(k,2),length(t))~=0
%     error('wrong size of time variable and k traj')
%     return
% else
%     trep = size(k,2)/length(t);
% end    
% k = cat(1,repmat(t,1,trep),k);

%tic;E=single(zeros(size(k,2),size(p,2)));toc;
E=single(exp(-sqrt(-1)*permute(hs,[2,1])*k));



