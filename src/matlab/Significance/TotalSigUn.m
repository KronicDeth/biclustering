function [Index]=TotalSigUn(R,S);
r=sigrow(R);
s=sigrow(S);
%All significant rows in each functional category
%now concatenate, sort, and delete redundancy
%'put together'
X=([r;s]);
%'and sorted'
Y=sortrows(X);
%'which rows are different'
rows=find(diff(Y)~=0);
%'Total Sig Rows'
Index=([1;Y([rows])]);
