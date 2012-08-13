clear;
load PreBioList14;
labelsig;
%FcnSig=TotalSig(A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q);
FcnSig=0;
UnSig=TotalSigUn(R,S);
%UnSig=0;
Rows_with_Unclass=Unclassrows(R,S);
load PreBioList14X;
labelsigcombo;
MixedSig=TotalSig(A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q);
[M,n]=size(MixedSig);
[X,x]=size(FcnSig);
[Y,y]=size(UnSig);
MixedSigNum=M-X-Y;
save PreBioList14X FcnSig UnSig Rows_with_Unclass MixedSig MixedSigNum -append

clear;
load PreBioList13;
labelsig;
%FcnSig=TotalSig(A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q);
FcnSig=0;
UnSig=TotalSigUn(R,S);
%UnSig=0;
Rows_with_Unclass=Unclassrows(R,S);
load PreBioList13X;
labelsigcombo;
MixedSig=TotalSig(A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q);
[M,n]=size(MixedSig);
[X,x]=size(FcnSig);
[Y,y]=size(UnSig);
MixedSigNum=M-X-Y;
save PreBioList13X FcnSig UnSig Rows_with_Unclass MixedSig MixedSigNum -append

clear;
load PreBioList12;
labelsig;
FcnSig=TotalSig(A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q);
%FcnSig=0;
UnSig=TotalSigUn(R,S);
%UnSig=0;
Rows_with_Unclass=Unclassrows(R,S);
load PreBioList12X;
labelsigcombo;
MixedSig=TotalSig(A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q);
[M,n]=size(MixedSig);
[X,x]=size(FcnSig);
[Y,y]=size(UnSig);
MixedSigNum=M-X-Y;
save PreBioList12X FcnSig UnSig Rows_with_Unclass MixedSig MixedSigNum -append

clear;
load PreBioList11;
labelsig;
FcnSig=TotalSig(A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q);
%FcnSig=0;
UnSig=TotalSigUn(R,S);
%UnSig=0;
Rows_with_Unclass=UnclassrowsBig(R,S);
load PreBioList11X;
labelsigcombo;
MixedSig=TotalSig(A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q);
[M,n]=size(MixedSig);
[X,x]=size(FcnSig);
[Y,y]=size(UnSig);
MixedSigNum=M-X-Y;
save PreBioList11X FcnSig UnSig Rows_with_Unclass MixedSig MixedSigNum -append