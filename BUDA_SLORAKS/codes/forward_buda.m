
function kdata = forward_buda(IM,csm,wmap,kmask,ky_idx_ap,ky_idx_pa)


[nx,ny,nc] = size(csm);


        mask_ap = kmask(:,:,:,1);
        mask_pa = kmask(:,:,:,2);
        
        I_ap = IM(:,:,1);
        I_pa = IM(:,:,2);
                
%         ky_idx_ap = mask_ap(1,:,1);
%         ky_idx_ap = find(ky_idx_ap>0);
%         ky_idx_pa = mask_pa(1,:,1);
%         ky_idx_pa = find(ky_idx_pa>0);
%         
        wmap_ap = wmap(:,:,:,1);
        wmap_pa = wmap(:,:,:,2);
                 
        mask_ap4d = zeros(nx,ny,nc,length(ky_idx_ap));
        mask_pa4d = zeros(nx,ny,nc,length(ky_idx_pa));              
        for t=1:length(ky_idx_ap)
        mask_ap4d(:,ky_idx_ap(t),:,t) = mask_ap(:,ky_idx_ap(t),:);
        mask_pa4d(:,ky_idx_pa(t),:,t) = mask_pa(:,ky_idx_pa(t),:);
        end

        k_ap = zeros(nx,ny,nc);
        k_pa = zeros(nx,ny,nc);

        parfor k=1:length(ky_idx_ap)
            
            tmp = repmat(I_ap,[1 1 nc]);           
            tmp = fft2c(tmp.*csm.*repmat(wmap_ap(:,:,k),[1 1 nc]));
            k_ap = k_ap + tmp.*mask_ap4d(:,:,:,k);
            
            tmp = repmat(I_pa,[1 1 nc]);           
            tmp = fft2c(tmp.*csm.*repmat(wmap_pa(:,:,k),[1 1 nc]));
            k_pa = k_pa + tmp.*mask_pa4d(:,:,:,k);
                   
        end
        
        kdata(:,:,:,1) = k_ap;
        kdata(:,:,:,2) = k_pa;
 