
function img = ift2_coil(kdata_ud,coil_sens)
% ns = size(coil_sens,3);
% nc = size(coil_sens,4);

img0 = sum(ift2(kdata_ud(:,:,:,:,1)).*conj(coil_sens),4);
imgC = sum(ift2(kdata_ud(:,:,:,:,2)).*conj(coil_sens),4);

img = img0+imgC;
end

