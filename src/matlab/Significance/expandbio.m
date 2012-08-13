function [adjust_biores]=expandbio(biores, index);
%adjusts an aggresively trimmed list to the regular list numbers

[R,x]=size(biores);
[I,y]=size(index);
i=1;
while i<=I;
   r=1;
   while r<=R;
      if index(i)>biores(r);
         r=r+1;
         if r==R+1;
            i=i+1;
            break;
         end;
         %the gene was not deleted before r, so let the number remain
         %check next item on biores list
      elseif index(i)==biores(r);
         %the gene WAS deleted RIGHT in front of r's value
         for j=r:R;
            biores(j)=biores(j)+1;
         end;
         %the numbers were all adjusted, 
         %so "room" has been made for the deleted gene number
         %now check next item on index
         i=i+1;
         break;
      else index(i)<biores(r);
         for j=r:R;
            biores(j)=biores(j)+1;
         end;
         %same thing, adjust gene numbers
         i=i+1;
         break;
      end;
   end;
end;

adjust_biores=biores;

         