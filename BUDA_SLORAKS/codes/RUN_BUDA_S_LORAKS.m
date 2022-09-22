%function x = BUDA_recon_V2
addpath BUDA_LORAKS_UY
clear;clc;
close all

nRF = 1;  % --1 conventional EPI
%--5 gSlider factor
AccZ = 1; %% SMS-3
AccY = 6;  %% inpalne -2 

load('esp.mat')

%function x = BUDA_recon
load(['sens_gre_gcc.mat'])
[nx,ny,nslice,nc] = size(sens_gre_gcc);
       

ky = length(esp);
t_pa = [0:AccY:AccY*ky-1]'.*esp*1e-4;
t_ap = t_pa(end:-1:1);

for dif = 18   

    load(['gcc_k_ap_pa_dif' num2str(dif) 'rf1.mat'])

for  slice = 1:16
   
        csm = squeeze(sens_gre_gcc(:,:,slice,:));
        wmap_0 =  2*pi*wmap(:,:,slice);
        
        %clear wmap_buda
        for t=1:length(ky_idx_ap)
        wmap_buda(:,:,t,1) = exp(-1i*wmap_0.*t_ap(t));
        wmap_buda(:,:,t,2) = exp(-1i*wmap_0.*t_pa(t));
        end
         
        k_ap = zeros(nx,ny,nc);
        k_ap(:,ky_idx_ap,:) = kspace_cor_ap(:,:,:,slice);
       
        k_pa = zeros(nx,ny,nc);
        k_pa(:,ky_idx_pa,:) = kspace_cor_pa(:,:,:,slice);
        
        kdata(:,:,:,1) = k_ap;
        kdata(:,:,:,2) = k_pa;
        
        kmask = zeros(size(kdata));
        kmask(abs(kdata)>0) = 1;
        
        %% RECON
        Niter = 20;
        tol = 1e-4;
        alpha = 0.0035;
        Display = 1;
        
        A = @(x) forward_buda(x,csm,wmap_buda,kmask,ky_idx_ap,ky_idx_pa);
        Ah = @(x) transpose_buda(x,csm,wmap_buda,kmask,ky_idx_ap,ky_idx_pa);
        AhA = @(x) Ah(A(x));
        AhA0 = @(x) Ah(A(x))+alpha.*x;
        
        tic
        Ahd = Ah(kdata);
        toc
        
        BUDA = CG_Recon(AhA0,Ahd,Niter,tol,Display);

        
 %% BUDA-LORAKS (multi-channel)
disp('BUDA-LORAKS Reconstruction');
[nx,ny, nc, ns] = size(kdata);  %ns: num shots
LORAKS_type = 1; % S-matrix (Support + Phase + parallell imaging)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$$$$$$$$$$$$$$$$$$$$$$$####################
%% S-LORAKS parameters
z = BUDA;
R = 3;
rank = 80;
lambda = 0.005;%0.003;%0.018;
tol = 1e-4;
max_iter = 6;  % change this based on necessary # iterations

[in1,in2] = meshgrid(-R:R,-R:R);
idx = find(in1.^2+in2.^2<=R^2);
patchSize = numel(idx);

coil_sens = cat(4,csm,csm);
coil_sens = permute(coil_sens,[1 2 4 3]);
%%
%%
%coil_sens(177:176*2,:,:,:) = circshift(coil_sens(177:176*2,:,:,:),[-44 0 0 0]);

B = @(x) ft2(coil_sens.*repmat(x,[1 1 1 nc]));      % SENSE encoding
Bh = @(x)  sum(conj(coil_sens).*ift2(x),4);         % Adjoint of SENSE encoding

P_M = @(x) LORAKS_operators(x,nx,ny,ns*nc,R,LORAKS_type,[]);
Ph_M = @(x) LORAKS_operators(x,nx,ny,ns*nc,R,-LORAKS_type,[]);

N1 = nx; N2 = ny; Nc = nc*ns;

ZD = @(x) padarray(reshape(x,[N1 N2 Nc]),[2*R, 2*R], 'post');
ZD_H = @(x) x(1:N1,1:N2,:,:);

%% 
tic
for iter = 1:max_iter
    
    z_cur = z;
    pz = B(z);
    MM = P_M(pz);           

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
        LhL = @(x) 2*Bh(reshape(L1(x)-L2(x),[nx ny ns nc]));
                
    end
    
    % data fitting
    M = @(x) AhA(x) + lambda*LhL(x);
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


BUDA_S_LORAKS(:,:,:,slice) = z;
BUDA_SENSE(:,:,:,slice) = BUDA;

 end
fn = (['Recon_for_CNN_Dif_' num2str(dif) '.mat'])

save(fn,'BUDA_S_LORAKS','BUDA_SENSE')
end
