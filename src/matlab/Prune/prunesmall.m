function [TC, TG, I]=prunesmall(smallC, smallG, bigC, bigG);
%Get sizes of matrices
[m,n]=size(smallC)
[m,g]=size(smallG);

I=[];

for (x=100001:m)
   bc=bigC;
   bg=bigG;
   for (y=1:n);
      [index,j]=find(bc==smallC(x,y));
      bc=bc([index],:);
      if isempty(bc);
         %contain=0;
         break;
      end;
   end;
   
   if bc;
      bg=bg([index],:);
      for (z=1:g);
         if (smallG(x,z)==0);
            break;
         else
            [index,k]=find(bg==smallG(x,z));
            bg=bg([index],:);
            if isempty(bg);
               contain=0;
               break;
            else
               contain=1;
            end;
         end;
      end;
      
      if contain;
         I=([I;x])
      end;
   end;
end;
TC=smallC;
TG=smallG;
         
      