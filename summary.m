function y = summary(var)
% summarizes a matrix telling how many elements of each value there are
valores=unique(var);
fprintf('\tVal\t#\n');
for i=1:length(valores)
    fprintf('\t%.3f\t%d\n',valores(i),sum(var==valores(i)));
end