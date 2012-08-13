
clear;
load PreBioListT14;
labelsig;
%FcnSig=TotalSig(A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q);
FcnSig=0;
UnSig=TotalSigUn(R,S);
%UnSig=0;
Rows_with_Unclass=UnclassrowsBig(R,S);
load PReBioList14TX;
labelsigcombo;
MixedSig=TotalSig(A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q);
[M,n]=size(MixedSig);
[X,x]=size(FcnSig);
[Y,y]=size(UnSig);
MixedSigNum=M-X-Y;
save PreBioList14TX FcnSig UnSig Rows_with_Unclass MixedSig MixedSigNum -append
'14'

clear;
load PreBioListT13;
labelsig;
%FcnSig=TotalSig(A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q);
FcnSig=0;
UnSig=TotalSigUn(R,S);
%UnSig=0;
Rows_with_Unclass=UnclassrowsBig(R,S);
load PReBioList13TX;
labelsigcombo;
MixedSig=TotalSig(A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q);
[M,n]=size(MixedSig);
[X,x]=size(FcnSig);
[Y,y]=size(UnSig);
MixedSigNum=M-X-Y;
save PreBioList13TX FcnSig UnSig Rows_with_Unclass MixedSig MixedSigNum -append
'13'

clear;
load PreBioListT12;
labelsig;
FcnSig=TotalSig(A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q);
%FcnSig=0;
UnSig=TotalSigUn(R,S);
%UnSig=0;
Rows_with_Unclass=UnclassrowsBig(R,S);
load PReBioList12TX;
labelsigcombo;
MixedSig=TotalSig(A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q);
[M,n]=size(MixedSig);
[X,x]=size(FcnSig);
[Y,y]=size(UnSig);
MixedSigNum=M-X-Y;
save PreBioList12TX FcnSig UnSig Rows_with_Unclass MixedSig MixedSigNum -append
'12'

clear;
load PreBioListT11;
labelsig;
FcnSig=TotalSig(A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q);
%FcnSig=0;
UnSig=TotalSigUn(R,S);
%UnSig=0;
Rows_with_Unclass=UnclassrowsBig(R,S);
load PReBioList11TX;
labelsigcombo;
MixedSig=TotalSig(A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q);
[M,n]=size(MixedSig);
[X,x]=size(FcnSig);
[Y,y]=size(UnSig);
MixedSigNum=M-X-Y;
save PreBioList11TX FcnSig UnSig Rows_with_Unclass MixedSig MixedSigNum -append
'11'