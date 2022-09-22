
function kdata_ud = ft2_coil(img,coil_sens,Circular_mask)

ns = size(coil_sens,3);
num_chan = size(coil_sens,4);

mask = repmat(Circular_mask,[1 1 ns num_chan]);

kdata = ft2(coil_sens.*repmat(img,[1 1 1 num_chan]));

kdata_ud(:,:,:,:,1) = kdata.*mask;
kdata_ud(:,:,:,:,2) = kdata.*(1-mask);

end