
function Cofcol=twodfbseone(img,dim)
img=double(imresize(img,[224 224]));
A=1:224;
A=repmat(A,224,1);
alfa=besselzero(1,224);
K=2/(224*224);

Bes=besselj(1,([1:224]'.*repmat(alfa/224,224,1)))./(besselj(0,repmat(alfa,224,1))).^2;
if dim==1
Aimg=A.*img;
Cofcol=K*Aimg*Bes;

% col
elseif dim==2
imgcol=img';
Aimgcol=A.*imgcol;
Cofcol=K*Aimgcol*Bes;
Cofcol=Cofcol';

% 2d
elseif dim==3
Aimg=A.*img;
Cof=K*Aimg*Bes;
imgcol=Cof';
Aimgcol=A.*imgcol;
Cofcol=K*Aimgcol*Bes;
Cofcol=Cofcol';
end

