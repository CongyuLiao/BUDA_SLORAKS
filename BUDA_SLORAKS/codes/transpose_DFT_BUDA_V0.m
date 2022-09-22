function Im = transpose_DFT_BUDA_V0(kdata,A_ap,A_pa,mask,kx,ky,num_chan)

%% ky: time segmentation 
signal_ap = ifftc(kdata(:,:,:,1).*mask(:,:,:,1),1);
signal_pa = ifftc(kdata(:,:,:,2).*mask(:,:,:,2),1);

    %%    
    rhs_ap = permute(signal_ap,[2 3 1]);
    rhs_ap = reshape(rhs_ap,ky*num_chan,kx); %% rhs_ap size 300 x 1 x 136
 
    rhs_pa = permute(signal_pa,[2 3 1]);
    rhs_pa = reshape(rhs_pa,ky*num_chan,kx); %% rhs_ap size 300 x 1 x 136 
    
clear Im
    for q=1:kx     
        Im(:,q,1) = A_ap(:,:,q)'*rhs_ap(:,q);        
        Im(:,q,2) = A_pa(:,:,q)'*rhs_pa(:,q);        
    end
    
Im = permute(Im,[2 1 3]);  
end