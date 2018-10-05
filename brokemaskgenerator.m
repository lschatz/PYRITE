%function [pyramidmask]=brokemaskgenerator(Npix, tripyramid, tabletoperror ,tabletop, roof, rooferror)
%Generates pyramid focal plane mask.

%%
%Saves the pyramid mask generated in a .mat file
function [pyramidmask]=brokemaskgenerator(Npix, tripyramid, tabletoperror ,tabletop, roof, rooferror)
x= (-Npix/2:Npix/2-1);
[X,Y]=meshgrid(x,x);
[phi,r]=cart2pol(X,Y); % phi from -pi to pi


tilt=.295;
zernikez={r.*cos(phi),r.*sin(phi)};


%Masks;
%%

if tripyramid==true


W1=-tilt*(cos(pi/3).*zernikez{1}+sin(pi/3).*zernikez{2});
W2= -tilt*(cos(pi/3).*zernikez{1}-sin(pi/3)*zernikez{2});
W3= tilt*(zernikez{1});
triaperture1=phi<=-pi/3;
triaperture2=phi>=pi/3;
triaperture3=triaperture1+triaperture2<1;



focal1= exp(1i*2*pi*W1).*triaperture1;%Complex wavefront
focal2= exp(1i*2*pi*W2).*triaperture2;%Complex wavefront
focal3= exp(1i*2*pi*W3).*triaperture3;%Complex wavefront

pyramidmask=focal1+focal2+focal3;


figure
imagesc(angle(pyramidmask));
axis equal  
  
save('broketripyramidmask.mat','pyramidmask')
%%
if tabletop==true
flat=pyramidmask(Npix/2, Npix/2);
pyramidmask(Npix/2-tabletoperror/2:Npix/2+tabletoperror/2,Npix/2-tabletoperror/2:Npix/2+tabletoperror/2)=flat;
%figure; imagesc(angle(pyramidmask));

save('broketabletripyramidmask.mat','pyramidmask')
end

end
    
%%
if tripyramid==false
mask=ones(Npix,Npix);
sz=size(mask);
mask2=mask;
mask3=mask;
mask4=mask;
mask(1:sz(1),sz(2)/2+1:sz(2))=0;
mask2(sz(1)/2+1:sz(1),1:sz(2))=0;
mask3=mask2==0;
mask4=mask==0;

% figure
% subplot(2,2,1)
% imagesc(mask)
% subplot(2,2,2)
% imagesc(mask2)
% subplot(2,2,3)
% imagesc(mask3)
% subplot(2,2,4)
% imagesc(mask4)
%%

mask11=mask.*mask2;
mask22=mask2.*mask4;
mask33=mask3.*mask;
mask44=mask4.*mask3;

% figure
% subplot(2,2,1)
% imagesc(mask11);
% subplot(2,2,2)
% imagesc(mask22);
% subplot(2,2,3)
% imagesc(mask33);
% subplot(2,2,4);
% imagesc(mask44);
figure; imagesc(mask11+2*mask22+3*mask33+4*mask44);

%%
W1=tilt*(cos(-3*pi/4).*zernikez{1}+sin(-3*pi/4).*zernikez{2});

W2= tilt*(cos(-pi/4).*zernikez{1}+sin(-pi/4).*zernikez{2});
W3= -tilt*(cos(pi/4).*zernikez{1}-sin(pi/4).*zernikez{2});
W4= -tilt*(cos(3*pi/4).*zernikez{1}-sin(3*pi/4).*zernikez{2});

focalPP= exp(1i*2*pi*W1).*mask11;%Complex wavefront
focalPM= exp(1i*2*pi*W2).*mask22;%Complex wavefront
focalMP= exp(1i*2*pi*W3).*mask33;%Complex wavefront
focalMM= exp(1i*2*pi*W4).*mask44;%Complex wavefront

pyramidmask=focalPP+focalPM+focalMP+focalMM;

figure; imagesc(angle(pyramidmask)); axis equal
save('brokequadpyramidmask.mat','pyramidmask');

%%   NOT CURRENTLY WORKING                 
if roof==true
%% roof error
quad2=imtranslate(quad2, [rooferror/2,0]);
    quad4=imtranslate(quad4,[-rooferror/2,0]);
        quad1=quad4+quad2<1;
        quad3=quad1;
            quad1(Npix/2+1:Npix,:)=0;
                quad3(1:Npix/2,:)=0;
            %figure; imagesc(quad1+2*quad2+3*quad3+4*quad4)
end
  
if tabletop==true 
flat=pyramidmask(Npix/2, Npix/2);
pyramidmask(Npix/2-tabletoperror/2:Npix/2+tabletoperror/2,Npix/2-tabletoperror/2:Npix/2+tabletoperror/2)=flat;
%figure; imagesc(angle(pyramidmask));

save('broketablequadpyramidmask.mat','pyramidmask')
end
end

% %% Continuing the debugging code
% pyramid=FTpupil.*pyramidmask;
% Pupilpyramid=abs(ifftshift(fft2(ifftshift(pyramid)))).^2;
% figure; imagesc(Pupilpyramid); axis equal












end