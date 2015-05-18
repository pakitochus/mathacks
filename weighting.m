function y = weighting(p,d,k)
sumatoria = 1;
for jarl=1:length(d),
    if ndims(p)==1
        sumatoria = sumatoria + d(jarl)*exp(-squeeze(abs(p(jarl)))/k);
    elseif ismatrix(p)
        sumatoria = sumatoria + d(jarl)*exp(-squeeze(abs(p(jarl,:)))/k);
    elseif ndims(p)==3
        sumatoria = sumatoria + d(jarl)*exp(-squeeze(abs(p(jarl,:,:)))/k);
    else
        error('P matrices with dimensions > 3 not supported');
    end
end
y = sumatoria/length(d);
end