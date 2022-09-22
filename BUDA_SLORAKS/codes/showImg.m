function [] = showImg( image, varargin )
%argument: kdata, name
%SHOWIMAGE Summary of this function goes here
%   Detailed explanation goes here

    if ~isreal(image)
        image = abs(image);
    end
    
    if ~ismatrix(image)
        image = squeeze(image);
    end
        
        
    if length(varargin) < 1
        name = '';
        imagesc(image);
    elseif length(varargin) < 2
        name = varargin{1};       
        imagesc(image);
        title(name);

    else 
        name = varargin{1};
        scale = varargin{2};
        imagesc(image,scale);
        title(name);
    end
    
% 
%     [N1, N2, N3, numImg] = size(kdata);
%             
%     for k = 1:numImg
%     subplot(numImg, 1 , k/2 + 1 
%     
    
    axis equal; axis tight; colormap gray; axis off; 
end

