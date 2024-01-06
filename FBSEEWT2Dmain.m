params.globtrend = 'none';
params.degree=5; % degree for the polynomial interpolation

% Choose the wanted regularization (none,gaussian,avaerage,closing)
params.reg = 'none';
params.lengthFilter = 10;
params.sigmaFilter = 1.5;

% Choose the wanted detection method (locmax,locmaxmin,ftc,scalespace)
params.detect = 'scalespace';
params.typeDetect='otsu'; %for scalespace:otsu,halfnormal,empiricallaw,mean,kmeans

params.N = 5; % maximum number of band for the locmaxmin method
params.completion=1;%% 0 if all possible subband are required

% Perform the detection on the log spectrum instead the spectrum
params.log=0;

% Choose the results you want to display (Show=1, Not Show=0)
Bound=1;   % Display the detected boundaries on the spectrum
Comp=1;    % Display the EWT components
Rec=0;     % Display the reconstructed image

i= imread('optic.png');% TEST IMAGE

ir=i(:,:,2); %green channel

 ir=imresize(ir,[224,224]);
 ir=adapthisteq(ir);
%% row
%tic
%cof=twodfbse(ig,1);
cof=twodfbseone(ir,1);
[M,N]=size(cof);
rowspe=(sum(abs(cof)))/M;
boundariesrow = EWT_Boundaries_Detect((rowspe),params);%

% boundary complition
if(params.completion==1)
if (length(boundariesrow)<(params.N-1))
    boundariesrow=EWT_Boundaries_Completion(boundariesrow,params.N-1);
end

if (length(boundariesrow)>(params.N-1))
    boundariesrow=boundariesrow(1:params.N-1);
end
end


%% column
%cofr=twodfbse(ig,2);
cofr=twodfbseone(ir,2);
[M,N]=size(cofr);
colspe=(sum(abs(cofr),2))/N;
boundariescol = EWT_Boundaries_Detect((colspe),params);%

% boundary complition
if(params.completion==1)
if (length(boundariescol)<(params.N-1))
    boundariescol=EWT_Boundaries_Completion(boundariescol,params.N-1);
end

if (length(boundariescol)>(params.N-1))
    boundariescol=boundariescol(1:params.N-1);
end
end

%% row filter design
boundaries=boundariesrow;
Npic=length(boundaries);
% We compute gamma accordingly to the theory
gamma=1;
for k=1:Npic-1
    r=(boundaries(k+1)-boundaries(k))/(boundaries(k+1)+boundaries(k));
    if r<gamma 
       gamma=r;
    end
end


gamma=(1-1/N)*gamma; %this ensure that gamma is chosen as strictly less than the min

mfb=cell(Npic+1,1);
% We start by generating the scaling function
mfb{1}=fbseewt_Meyer_Scaling(boundaries(1),gamma,N);
% We generate the wavelets
for k=1:Npic-1
   mfb{k+1}=fbseewt_Meyer_Wavelet(boundaries(k),boundaries(k+1),gamma,N); 
end
%mfb{Npic+1}=fbseewt_Meyer_Wavelet(boundaries(end),N,gamma,N); 
mfb{Npic+1}=fbseewtlast_Meyer_Wavelet(boundaries(end),N,gamma,N); 
mfbR=mfb;
%% col filter design
boundaries=boundariescol;
Npic=length(boundaries);

gamma=1;
for k=1:Npic-1
    r=(boundaries(k+1)-boundaries(k))/(boundaries(k+1)+boundaries(k));
    if r<gamma 
       gamma=r;
    end
end
gamma=(1-1/M)*gamma;

mfb=cell(Npic+1,1);
% We start by generating the scaling function
mfb{1}=fbseewt_Meyer_Scaling(boundaries(1),gamma,M);

% We generate the wavelets
for k=1:Npic-1
   mfb{k+1}=fbseewt_Meyer_Wavelet(boundaries(k),boundaries(k+1),gamma,M); 
end
mfb{Npic+1}=fbseewtlast_Meyer_Wavelet(boundaries(end),M,gamma,M); 
mfbC=mfb;




%% Filtering
 ewtR=cell(length(mfbR),1);
 imcof=twodfbseone(ir,3);
 for k=1:length(mfbR)
     filter=[repmat(mfbR{k},1,M)]';
     filtered=filter.*imcof; 
      ewtR{k}=twodinvfbseone(filtered,3);
 end

 
 
 ewtC=cell(length(mfbR),length(mfbC));

for k=1:length(mfbC)
    filtercol=repmat(mfbC{k},1,N);
    for s=1:length(mfbR)
       imC= ewtR{s} ;%
        fftim=twodfbseone(imC,3);
        filteredcol=filtercol.*fftim;
        ewtcoffbse{s}{k}=filteredcol;
        ewtC{s}{k}=twodinvfbseone(filteredcol,3);  
    end
end
%toc


%% EWTC is decomposed images adding all will give orignal reconstructed image
%for sample sample showing few image output

imshow(uint8(ewtC{1}{1}),[]);
figure,imshow(uint8(ewtC{1}{2}),[]);
figure,imshow(uint8(ewtC{1}{3}),[]);
figure,imshow(uint8(ewtC{1}{4}),[]);
figure,imshow(uint8(ewtC{1}{5}),[]);
figure,imshow(uint8(ewtC{2}{1}),[]);
figure,imshow(uint8(ewtC{2}{1}),[]);
%% CAn apply on diffrent can concatine to get decompsed RGB image
