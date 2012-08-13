function [Index]=TotalSig(A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q);
a=sigrow(A);
b=sigrow(B);
c=sigrow(C);
d=sigrow(D);
e=sigrow(E);
f=sigrow(F);
g=sigrow(G);
h=sigrow(H);
i=sigrow(I);
j=sigrow(J);
k=sigrow(K);
l=sigrow(L);
m=sigrow(M);
n=sigrow(N);
o=sigrow(O);
p=sigrow(P);
q=sigrow(Q);
%All significant rows in each functional category
%now concatenate, sort, and delete redundancy
%'put together'
X=([a;b;c;d;e;f;g;h;i;j;k;l;m;n;o;p;q]);
%'and sorted'
Y=sortrows(X);
%'which rows are different'
rows=find(diff(Y)~=0);
%'Total Sig Rows'
Index=([1;Y([rows])]);
