function result = isumsq( kdata )
%SUMSQ Summary of this function goes here
%   Detailed explanation goes here
    
    N = ndims(kdata);
    if N == 3
        result = sqrt(sum(abs(fftshift(ifft2(ifftshift(kdata)))).^2, 3));
    elseif N == 4
        result = sqrt( sum(abs(ift3(kdata)).^2,4));
    end
%     result = sqrt(sum(abs(fftshift(ifft2(ifftshift(kdata(:, :,:))))).^2, 3));
    

end

