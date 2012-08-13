function total = width(A, width)

[genes, conditions] = size(A);

conditionCombos = choosenk(conditions, width);
[combos, width] = size(conditionCombos); 

total = 0;
for i = 1 : combos
    [Patterns, Genes, Counts] = bfbc(A(:, conditionCombos(i, :)));
    Patterns = conditionCombos(i, Patterns);
    Genes = Genes;
    total = total + size(Patterns, 2) / width;
end