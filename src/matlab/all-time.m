Bc = Ac(:, 1:10);
for i=2:10
    i
    tic, total = width(Bc, i), toc
end