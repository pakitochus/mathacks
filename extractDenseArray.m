function [X,Y,Z] = extractDenseArray(imSize, cubeSize)

if length(imSize)~=length(cubeSize)
    error('Dimension of image must be the same as the dimension of the cube')
end
centro = floor(imSize/2); % floor para que sea exacto o el mínimo

vectX = [centro(1):-cubeSize(1):ceil(cubeSize(1)/2), centro(1):cubeSize(1):(imSize(1)-ceil(cubeSize(1)/2))];
vectX = sort(vectX(2:end));

vectY = [centro(2):-cubeSize(2):ceil(cubeSize(2)/2), centro(2):cubeSize(2):(imSize(2)-ceil(cubeSize(2)/2))];
vectY = sort(vectY(2:end));

vectZ = [centro(3):-cubeSize(3):ceil(cubeSize(3)/2), centro(3):cubeSize(3):(imSize(3)-ceil(cubeSize(3)/2))];
vectZ = sort(vectZ(2:end));

[X,Y,Z] = meshgrid(vectX, vectY, vectZ);

