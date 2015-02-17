function h = superimposeImage(im, map, varargin)
% SUPERIMPOSEIMAGE superimposes a color map MAP to a grayscale image IM
% using IM as an alpha map (transparency index), so that the final result
% is fake coloured. 
% 
%   Version 0.1 September 2014
%   Francisco J. Martinez. 
%   
%   SUPERIMPOSEIMAGE(im, map) depicts the grayscale image IM fake coloured
%   using the color information in MAP. If MAP is a grayscale image, it
%   uses the 'jet' default colormap. If it is a true color image (a matrix 
%   NxMx3), it uses its colour. 
%
%   SUPERIMPOSEIMAGE(..., 'colormap', cm) uses a user-defined colormap
%   instead of the default 'jet'. E.g. 'bone', 'hot', etc. 
% 
%   SUPERIMPOSEIMAGE(..., 'alpha', alpha) controls the transparency of IM
%   by adjusting its brightness. An alpha < 1 will reduce the transparency
%   of the fake color map, whereas an alpha > 1 will increase the
%   transparency of this maps. 
% 
%   h = SUPERIMPOSEIMAGE(...) returns the handle of the current figure. 
% 
%   NOTES: 
%       - Now a colorbar can be added, and it will represent exactly the
%       range of values in MAP. 

% Defaults
cm='jet';
alpha = 1;

% Check some basic requirements of the data
if (nargin == 0)||(nargin==1),
    error ('You must supply the base image and the image to be applied. Usage: superimposeImage( image , map )');
end
if ~all(size(im)==size(map))
    error ('The two images are different size.');
end
% Read options. 
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        if ~ischar (varargin{i}),
            error (['Unknown type of optional parameter name (parameter' ...
                ' names must be strings).']);
        end
        % change the value of parameter
        switch lower (varargin{i})
            case 'alpha'
                alpha = varargin{i+1};
            case 'colormap'
                cm = varargin{i+1};
            otherwise
                % Hmmm, something wrong with the parameter string
                error(['Unrecognized parameter: ''' varargin{i} '''']);
        end;
    end;
end

% Displays the grayscale base image.
imagesc(repmat(entre0y1(im-min(im(:))), [1,1,3]));
hold on;
if(min(map(:))<max(map(:)))
    h=imagesc(map,[min(map(:)), max(map(:))]);colormap(cm);
else
    
    h=imagesc(map,[min(map(:)), min(map(:))+1]);colormap(cm);
end
hold off;

% Sets the transparency for the map layer.
set(h, 'AlphaData', alpha*entre0y1(im)+(1-alpha));
end

function img = entre0y1(img)
img(isnan(img))=0;
aux=min(img(:));
tam = size(img);
img = img - aux;
aux = max(img(:));
maximum = ones(tam,'single')*single(aux);
img = single(img) ./ maximum;
end