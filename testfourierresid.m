%function []=fourierresid(npix, nn, mm, fourierphaseplus, fourierphaseminus, triresidual, quadresidual)
%end

%NOTES ON VARIABLE NAMES:
    %varname_plus/minus_plus/minus: first plus/minus refers to the positive
    %or negative spatial frequency. The second plus/minus refers to the
    %sign of the variable 'm' which definies the spatial frequency number.

    
  %% 3PWFS
  
triresidual=kscreen;
sz=size(triresidual);
szp=size(fourierphaseminusminus);
sz0=size(fourierphaseminuszero);

for i=1:sz(3)
    for j=1:szp(3)

%% averaging for positive m
     trimodeplusplus= triresidual(:,:,i).*fourierphaseplusplus(:,:,j);
     trimodeminusplus=triresidual(:,:,i).*fourierphaseminusplus(:,:,j);
     triplussumplus= sum(sum(trimodeplusplus));
     triminussumplus= sum(sum(trimodeminusplus));
     
     triminplus(:,j,i)=triminussumplus;
     triplusplus(:,j,i)=triplussumplus;
     trisumplus(:,j,i)=sqrt(triplussumplus.^2+triminussumplus.^2);
%% averaging for ngegative m
      trimodeplusminus= triresidual(:,:,i).*fourierphaseplusminus(:,:,j);
     trimodeminusminus=triresidual(:,:,i).*fourierphaseminusminus(:,:,j);
     triplussumminus= sum(sum(trimodeplusminus));
     triminussumminus= sum(sum(trimodeminusminus));
     
     triminminus(:,j,i)=triminussumminus;
     triplusminus(:,j,i)=triplussumminus;
     trisumminus(:,j,i)=sqrt(triplussumminus.^2+triminussumminus.^2);
    end
end

%% m=zero
for i=1:sz(3)
    for j=1:sz0(3)
trimodepluszero= triresidual(:,:,i).*fourierphasepluszero(:,:,j);
     trimodeminuszero=triresidual(:,:,i).*fourierphaseminuszero(:,:,j);
     triplussumzero= sum(sum(trimodepluszero));
     triminussumzero= sum(sum(trimodeminuszero));
     triminzero(:,j,i)=triminussumzero;
     tripluszero(:,j,i)=triplussumzero;
     trisumzero(:,j,i)=sqrt(triplussumzero.^2+triminussumzero.^2);      
    end
end
%% averaging the +/- m into a single spatial frequency
triaverage=sqrt(trisumminus.^2+trisumplus.^2);

%%
mcount=[];
ncount=[];
success=1;

for n=0:nn
    for m=0:mm
       mcount(:,success)=m;
       ncount(:,success)=m;
       success=success+1;   
    end
end
success1=1;
success2=1;
sz=size(mcount);

for i=1:sz(2)
    if mcount(i)==0
        x=trisumzero(:,success2,:);
        fouriertotal(:,i,:)=x;
        success2=success2+1;
    else
        x=triaverage(:,success1,:);
        fouriertotal(:,i,:)=x;
        success1=success1+1;
    end
end

