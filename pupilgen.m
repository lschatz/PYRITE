%Pupil Gen

npix= 256; %number of pixels across the pupil diameter
lambda= 700.*10^-9;
error= 10*10^-9;
pupil=[];
%ncount=[];
%mcount=[];
%counter=1;
%G_mat=[];
success=1;
%Scanning through Z -35 to Z 35
%%
for n=0:5
    for m=-5:5
%ncount(counter)=n;
%mcount(counter)=m;
%counter=counter+1;
ma = abs(m);
if m==0 & n==0
    continue
elseif mod(n-ma,2)~=0
        continue
    elseif n<ma
        continue
    else  
% Generate the pupil with WF error        
    pupil= ((2*pi)/lambda)*error*zernike(n,m, npix);        
    %pupil = complex(zeros(Npix));
    %pupil(Npix/2-npix/2:Npix/2+npix/2-1,Npix/2-npix/2:Npix/2+npix/2-1) =ef;
    orgpupil(:,:,success)=pupil;
    success=success+1;
    
    end
    end
end

