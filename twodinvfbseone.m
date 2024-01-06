function imgnew1=twodinvfbseone(Cofcol,dim)
% inverse
alfa=besselzero(1,224);
Besinv=besselj(1,(alfa/224)'*(1:224));
if dim==1
imgnew1=Cofcol*Besinv;
imshow(uint8(imgnew1),[]);

elseif dim==2
imgnew1=Cofcol'*Besinv;
imshow(uint8(imgnew1'),[]);

elseif dim==3
 imgnew=Cofcol'*Besinv;    
 imgnew1=imgnew'*Besinv;
 
end