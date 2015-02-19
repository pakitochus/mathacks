function y = fromStructToLatex(perfStruct, name)
% Creates a latex table from perfStruct, which can be a single strucutre or
% an array of structures, and present the data in a LaTeX table containing
% values of accuracy, sensitivty and specificity, along with their standard
% errors, and includes the name "name" in the first column.

fprintf(' & Accuracy & Sensitivity & Specificity \\\\\n\\hline\n');
for j=1:numel(perfStruct)
    fprintf('%s & $%.3f \\pm %.3f$ & $%.3f \\pm %.3f$ & $%.3f \\pm %.3f$ \\\\\n',...
        name, perfStruct(j).CorrectRate(1),perfStruct(j).CorrectRate(2), ...
        perfStruct(j).Sensitivity(1),perfStruct(j).CorrectRate(2),...
        perfStruct(j).Specificity(1),perfStruct(j).Specificity(2));
end