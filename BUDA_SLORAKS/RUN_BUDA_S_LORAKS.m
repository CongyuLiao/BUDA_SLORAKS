
clc
clear
close all

addpath mtimesx_20110223
addpath codes

fn  = 'raw_data/'

load([fn 'sens_gre_gcc.mat'])  % load sensitivity maps
load([fn 'esp.mat'])   % load echo-spacing of circular-EPI


Niter = 20;   % iterations of BUDA
tol = 1e-4;
%alpha = 0.0035;
Display = 1;

for dif = 1 %% 5:16
    load([fn 'gcc_k_ap_pa_dif' num2str(dif) 'rf1.mat'])
    
    sc = 0.5*(max(abs(kspace_cor_ap(:)))+max(abs(kspace_cor_pa(:))));
    
    kspace_cor_ap = kspace_cor_ap./sc;
    kspace_cor_pa = kspace_cor_pa./sc;
    
    %%
    for slice = 8
        
        kdata(:,:,:,1) = kspace_cor_ap(:,:,:,slice);
        kdata(:,:,:,2) = kspace_cor_pa(:,:,:,slice);
        %calculate Circular_mask
        tmp = squeeze(sum(abs(kdata),3));
        mask = zeros(size(tmp));
        mask(tmp>0) = 1;
        
        mask = repmat(mask,[1 1 1 size(kdata,3)]);
        mask = permute(mask,[1 2 4 3]);
        
        
        csm = squeeze(sens_gre_gcc(:,:,slice,:));  % coil sensitivity maps
        
        wmap0 = 2*pi*wmap(:,:,slice);  % fieldmap
        %%
        [A_ap, A_pa] =  FT_BUDA(csm,wmap0,esp,ky_idx_ap,ky_idx_pa);
        
        %%
        nx = 300;  % nx of cEPI image
        ny = 300;  % ny of cEPI image
        kx = 300;  % kx of cEPI kspace
        ky = 50;  % ky of cEPI kspace (R4 w/ p.f.)
        num_chan = 12;  % number of channel
        
        At=@(x) transpose_DFT_BUDA_V0(x,A_ap,A_pa,mask,kx,ky,num_chan);
        A = @(x) forward_DFT_BUDA_V0(x,A_ap,A_pa,mask,nx,ny,kx,ky,num_chan);
        
        AtA = @(x) At(A(x));
        
        Ahd = At(kdata);
        
        BUDA = CG_Recon(AtA,Ahd,Niter,tol,Display);   % initialization with BUDA recon

        
        %% BUDA-LORAKS (multi-channel)
        disp('BUDA-LORAKS Reconstruction');
        ns = size(kdata,4);  %ns: num shots
        LORAKS_type = 1; % S-matrix (Support + Phase + parallell imaging)
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$$$$$$$$$$$$$$$$$$$$$$$####################
        %% S-LORAKS parameters
        z = BUDA;
        R = 3;
        rank = 80;
        lambda = 0.0008;%0.003;%0.018;
        tol = 1e-4;
        max_iter = 8;  % change this based on necessary # iterations of S-loraks
        
        [in1,in2] = meshgrid(-R:R,-R:R);
        idx = find(in1.^2+in2.^2<=R^2);
        patchSize = numel(idx);
        
        coil_sens = cat(4,csm,csm);
        coil_sens = permute(coil_sens,[1 2 4 3]);
        
        %%
        B = @(x) ft2(coil_sens.*repmat(x,[1 1 1 num_chan])); % SENSE encoding
        Bh = @(x)  sum(conj(coil_sens).*ift2(x),4); % Adjoint of SENSE encoding
        
        P_M = @(x) LORAKS_operators(x,nx,ny,ns*num_chan,R,LORAKS_type,[]);
        Ph_M = @(x) LORAKS_operators(x,nx,ny,ns*num_chan,R,-LORAKS_type,[]);
        
        %%
        ns = 2;
        N1 = nx;
        N2 = ny;
        Nc = num_chan*ns;
        
        ZD = @(x) padarray(reshape(x,[N1 N2 Nc]),[2*R, 2*R], 'post');
        ZD_H = @(x) x(1:N1,1:N2,:,:);
        
        %%
        tic
        for iter = 1:max_iter
            
            z_cur = z;
            pz = B(z);
            MM = P_M(pz);
            l = sum(MM,1);
            l = find(abs(l)>0);
            MM = MM(:,l);
            
            Um = svd_left(MM);
            nmm = Um(:,rank+1:end)'; % null space
            Bhr = 0;
            
            if LORAKS_type == 1 % S
                
                nf = size(nmm,1);
                nmm = reshape(nmm,[nf, patchSize, 2*Nc]);
                nss_h = reshape(nmm(:,:,1:2:end)+1j*nmm(:,:,2:2:end),[nf, patchSize*Nc]);
                Nis = filtfilt(nss_h,'C',N1,N2,Nc,R);
                Nis2 = filtfilt(nss_h,'S',N1,N2,Nc,R);
                
                L1 = @(x) ZD_H(ifft2(squeeze(sum(Nis.*repmat(fft2(ZD(B(x))),[1 1 1 Nc]),3))));
                L2 = @(x) ZD_H(ifft2(squeeze(sum(Nis2.*repmat(conj(fft2(ZD(B(x)))),[1 1 1 Nc]),3))));
                LhL = @(x) 2*Bh(reshape(L1(x)-L2(x),[nx ny ns num_chan]));
                
            end
            
            % data fitting
            M = @(x) AtA(x) + lambda*LhL(x);
            % [z,~] = pcg(M, Ahd + lambda*Bhr);
            z = CG_Recon(M,Ahd+lambda*Bhr,Niter,tol,Display);
            
            t = (norm(z_cur(:)-z(:))/norm(z(:)));
            
            % display the status
            if ~rem(iter,1)
                disp(['iter ' num2str(iter) ', relative change in solution: ' num2str(t)]);
            end
            
            % check for convergence
            if t < tol
                disp('Convergence tolerance met: change in solution is small');
                break;
            end
            
        end
        toc
        
        
        BUDA_S_LORAKS = z;
        BUDA_SENSE = BUDA;
        
        
    end
end

kk_sense = fft2c(BUDA_SENSE);
figure,imagesc(log(abs(kk_sense(:,:)))), colormap gray    % kspace of BUDA w/ sense recon
title('BUDA-SENSE kspace (AP/PA)')% AP+PA shot
figure,imagesc(cat(2,rot90(abs(BUDA_SENSE(:,:,1)),3),rot90(abs(BUDA_SENSE(:,:,2)),3))), colormap gray
title('BUDA-SENSE recon (AP/PA)')% AP+PA shot


kk_sloraks = fft2c(BUDA_S_LORAKS);
figure,imagesc(log(abs(kk_sloraks(:,:)))), colormap gray    % kspace of BUDA w/ S-loraks recon
title('BUDA S-LORAKS kspace (AP/PA)')% AP+PA shot
figure,imagesc(cat(2,rot90(abs(BUDA_S_LORAKS(:,:,1)),3),rot90(abs(BUDA_S_LORAKS(:,:,2)),3))), colormap gray
title('BUDA S-LORAKS recon (AP/PA)')% AP+PA shot







