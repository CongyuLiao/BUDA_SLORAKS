function kdata = forward_DFT_BUDA_V0(Im,A_ap,A_pa,mask,nx,ny,kx,ky,num_chan)

    %% Forward
    Im = permute(Im,[2 1 3]);
    Iap = reshape(Im(:,:,1),[nx 1 ny]);
    Ipa = reshape(Im(:,:,2),[nx 1 ny]);
    
    clear Kap Kpa 
    for q=1:nx     
        Kap(:,q) = A_ap(:,:,q)*squeeze(Iap(:,1,q)); 
        Kpa(:,q) = A_pa(:,:,q)*squeeze(Ipa(:,1,q)); 
    end
    
    Kap = permute(Kap,[2 1]);
    Kap = reshape(Kap,[kx ky num_chan]);
    Kap = fftc(Kap,1);
    kdata(:,:,:,1) = Kap.*mask(:,:,:,1);
    
    Kpa = permute(Kpa,[2 1]);
    Kpa = reshape(Kpa,[kx ky num_chan]);
    Kpa = fftc(Kpa,1);
    kdata(:,:,:,2) = Kpa.*mask(:,:,:,2);
    
end