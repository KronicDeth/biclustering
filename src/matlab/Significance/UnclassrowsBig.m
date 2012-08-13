function [Index]=UnclassrowsBig(R,S);
[r,x]=find(R==1);
[s,y]=find(S==1);
%r=r';
%s=s';
X=sortrows(([r;s]));
Y=diff(X);
[rows,cols]=find(Y~=0);
Index=([1;X([rows+1])]);