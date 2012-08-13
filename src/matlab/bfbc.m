function [Patterns, Genes_in_Pattern, Mcount] = bfbc(Ac)
[M, N] = size(Ac);

% Re-rank conditions across all genes
for k = 1:M
    z = Ac(k, :);
    [Y, I] = sort(z);
    [I1, I2] = sort(I);
    Acc(k, :) = I2;
    Bcc(k, :) = I;
end

% Sort according to patterns
[B, index] = sortrows(Bcc);
Li = size(index);
Mindex = index';

% Diff to find boundaries between similar patterns
DB = abs(diff(B));
DB = sum(DB, 2);

% Finds all distinct patterns and add a one to capture first one
m1 = find(DB ~= 0) + 1;
m1= [1 ; m1];

% Find how many of each pattern
m11 = [m1 ; Li(1) + 1];
Mcount = diff(m11);
Patterns = B(m1, :);
LMcount = size(Mcount);
Genes_in_Pattern = zeros(LMcount(1), max(Mcount));

for k = 1:LMcount(1)
    Mm = Mindex(1, m11(k) : m11(k + 1) - 1);
    [Mm1, Mm2] = size(Mm);
    Genes_in_Pattern(k, 1 : Mm2) = Mm;
end