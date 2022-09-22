function IMG = transpose_DFT_BUDA(kdata,A_ap,A_pa,mask,kx,ky,num_chan)

%% ky: time segmentation 
signal_ap = ifftc(kdata(:,:,:,1).*mask(:,:,:,1),1);
signal_pa = ifftc(kdata(:,:,:,2).*mask(:,:,:,2),1);

signal_apZ = ifftc(kdata(:,:,:,3).*(1-mask(:,:,:,1)),1);
signal_paZ = ifftc(kdata(:,:,:,4).*(1-mask(:,:,:,2)),1);

    %%    
    rhs_ap = permute(signal_ap,[2 3 1]);
    rhs_ap = reshape(rhs_ap,ky*num_chan,kx); %% rhs_ap size 300 x 1 x 136
    rhs_apZ = permute(signal_apZ,[2 3 1]);
    rhs_apZ = reshape(rhs_apZ,ky*num_chan,kx); %% rhs_ap size 300 x 1 x 136
     
    rhs_pa = permute(signal_pa,[2 3 1]);
    rhs_pa = reshape(rhs_pa,ky*num_chan,kx); %% rhs_ap size 300 x 1 x 136
    rhs_paZ = permute(signal_paZ,[2 3 1]);
    rhs_paZ = reshape(rhs_paZ,ky*num_chan,kx); %% rhs_ap size 300 x 1 x 136
   
    
clear Im
    for q=1:kx     
        Im(:,q,1) = A_ap(:,:,q)'*rhs_ap(:,q);
        Im(:,q,2) = A_pa(:,:,q)'*rhs_pa(:,q); 
        Im(:,q,3) = A_ap(:,:,q)'*rhs_apZ(:,q);
        Im(:,q,4) = A_pa(:,:,q)'*rhs_paZ(:,q); 
    end
    
    IMG(:,:,1) = Im(:,:,1)+Im(:,:,3);
    IMG(:,:,2) = Im(:,:,2)+Im(:,:,4);
end