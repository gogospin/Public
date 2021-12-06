clear img_recon;
for i = 7
    
    data = sq(data0(:,:,:,i,:));% my pulseq data
    data = fft1Dshift(data,3);
    
    

  
    %data = data.*repmat(reshape(exp(1i*phi.phi(:)),[],size(data,2),size(data,3)),1,1,1,size(data,4));
    % end of phase corr
    fprintf('%s %d\n','dynamic',i);
    % CG Alg recon --- needs more RAM
    %mbgap = size(data,3);
    %clear img_recon;
    fbar = waitbar(0,'CG SENSE in progress');
    for idx_sli = 1:size(data,3)
        % organizing MR data
        waitbar(idx_sli/size(data,3),fbar,'CG SENSE in progress');
        %rawdata = sq(data(:,1:end/2,idx_sli,:));
        %
        rawdata = sq(data(:,:,idx_sli,:));
        [nFE,nSpokes,nCh]=size(rawdata);
        % prepare CG matrix AND run recon
        maxit = 20;
        u = single(rawdata);
        u = reshape(u,[],nCh);
        s = single(sq(maps(:,:,:,idx_sli)));
        s = single(reshape(s,[],nCh));
        p = pos(:,:,:,idx_sli);
        p = single(reshape(p,3,[]));
        
        
        k = kdata0(:,:,idx_sli,i-1); k = k(1:4,:);
        k = reshape(k,size(k,1),[]);
        
        k(2,:) = -k(2,:);
        k(3,:) = -k(3,:);
        if recon.b0corr==1
            jump_rate = size(k,2);
            b0 = reshape(b0data(:,:,idx_sli),[],1)*2*pi;
            [de,out] = my_CGrecon_wb0_acc_loopi_gpu(u,k,s,p,b0,t_vec(:),jump_rate,maxit,0.00000);
        else
            jump_rate = size(k,2);
            [de,out] = my_CGrecon_nob0_acc_loopi_gpu(u,k,s,p,jump_rate,maxit,0.00000);
        end
        % run CG
        
        out = reshape(out,mri.xres,mri.yres,[]);
        img_recon(:,:,idx_sli,i) = out(:,:,end);
    end
    close (fbar);
end

img_recon=gather(img_recon);
