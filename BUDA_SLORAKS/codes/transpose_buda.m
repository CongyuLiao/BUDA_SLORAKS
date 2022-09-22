function IM = transpose_buda(kdata,csm,wmap,kmask,ky_idx_ap,ky_idx_pa)

[nx,ny,nc] = size(csm);

        mask_ap = kmask(:,:,:,1);
        mask_pa = kmask(:,:,:,2);
        
        k_ap = kdata(:,:,:,1).*mask_ap;
        k_pa = kdata(:,:,:,2).*mask_pa;
        
        
%         ky_idx_ap = mask_ap(1,:,1);
%         ky_idx_ap = find(ky_idx_ap>0);
%         ky_idx_pa = mask_pa(1,:,1);
%         ky_idx_pa = find(ky_idx_pa>0);
%         
        wmap_ap = wmap(:,:,:,1);
        wmap_pa = wmap(:,:,:,2);
        
        I_ap = zeros(nx,ny);
        I_pa = zeros(nx,ny);
        
        parfor k=1:length(ky_idx_ap)
            
            tmp = zeros(size(k_ap));
            tmp(:,ky_idx_ap(k),:) = k_ap(:,ky_idx_ap(k),:);
            
            I_ap = I_ap+sum(ifft2c(tmp).*conj(csm),3).*conj(wmap_ap(:,:,k));
                     
            tmp = zeros(size(k_pa));
            tmp(:,ky_idx_pa(k),:) = k_pa(:,ky_idx_pa(k),:);
            
            I_pa = I_pa+sum(ifft2c(tmp).*conj(csm),3).*conj(wmap_pa(:,:,k));
            
        end
        
        IM(:,:,1) = I_ap;
        IM(:,:,2) = I_pa;
        
        
            
  