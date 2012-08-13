function [Index]=sigrow(biosig);
%Takes a significance matrix, 
%and returns index of significant rows
%use this to find total significance

Index=[];
i=1;
[M,N]=size(biosig);
while i<=M;
   j=1;
   while j<=N;
      if biosig(i,j)==0;
         %NOT a significant row
         break;
      else
         %SIGNIFICANT ROW
         j=j+1;
         if j==N+1;
            Index=([Index;i]);
         end;
      end
   end
   i=i+1;
end;