%PYRITE: Pyramid Residual Wavefront Experiment 

%Simulator Goals:
    %Simulates a correct model of atmospheric turbulence
    %Simulates a refractive 3PWFS and 4PWFS
    %Accurately reconstructs the wavefront
    %Calculates residual wavefront error
    %Fits residual wavefront error with Fourier Modes.
    
% Lauren Schatz
% May 15, 2018

clear all, close all;

%% Simulation Setup:
% In this section are the system parameters for the simulator.


error= 10*10^-9; %Amount of wavefront error in nanometers applied in the pupil.
    npix= 256; %number of pixels across the pupil diameter
        Npix= 2560; % number of pixels across the total roster. Represents the amount of zero padding.
            lambda= 700.*10^-9; %wavelength in nanometers

nphot=2.5*10^6; % The total number of photons on the detector (from summing all the photons in detection), before noise.
    rnphot=10^15; % The total number of photons on the detector used in building the reconstructor matrix.
        rdns= 3; % Read noise applied

sampling=32; % 128; 32 %Number of pixels across each pupil. Can be 16, 32, 64
    nmodes=6; %Number of fourier modes used to fit the residual WF error
        nn=11; %Zernike order Z_nm used in reconstructor matrix calibration.
            mm=11; % Zernike order. (Will be +/-m)
nscreen=1;

%rooferror=3; % in pixels. If you have a roof error, number of pixels in diameter creating the roof effect. 
%tabletoperror=10; %in pixels. If you have a tabletop error, number of pixels in diameter creating the flattened tabletop of the pyramid tip. 




%% Set Up

%Simulator setup
runsim=true; %toggle true/false to run complete sim. If false will only generate reconstruction matrix.
duelrun=true; %toggle true/false to run 3PWFS and 4PWFS at the same time. 


tripyramid=true; %toggle true/false for 3PWFS simulation. If false, will run 4PWFS simulation
have_reconstructor=false; %toggle true/false for reconstructor matrix
have_pyramidmask=true;  %toggle true/fasle for pyramid focal plane maske.

roof=false; %toggle true/false for roof error
tabletop=false; % toggle true/false for tabletop error
    
MVM=true; %Toggle true/false. True: Will run Matrix-Vector-Multiply (intensity Centroids) for reconstructor matrix
                             %Flase: Will use the full frame (Olivier
                             %Guyon's method) for the reconstructor matrix

kolmogorov=true;
broke=false;
%% Run Simulation

%% Load in or generate the pyramid focal plane masks
if duelrun==false
    
if have_pyramidmask==false
    fprintf('generating pyramid mask')
    pyramidmask=maskgenerator(Npix, tripyramid, tabletoperror, tabletop, roof, rooferror);
    
end
    
if have_pyramidmask==true && tripyramid==true
    fprintf('loading pyramid mask')
    load 'tripyramidmask.mat';
   
end

if have_pyramidmask==true && tripyramid==false
    fprintf('loading pyramid mask')
    load 'quadpyramidmask.mat';
end



%% Load in or generate the reconstructor matrix
if have_reconstructor==false
    fprintf('Generating Reconstructor matrix')
    [rmatrix,success, Sig, ppupil]=reconstructorgenerator(nn,mm,rdns, rnphot, npix,Npix, pyramidmask, lambda, error, sampling, tripyramid, MVM);
end

if have_reconstructor==true
    if tripyramid==true && MVM ==true
        fprintf('3PWFS MVM')
        r=load('MVMtrireconstructormatrix.mat');
        rmatrix=r.rmatrix;
    end
    if tripyramid==true && MVM ==false
        fprintf('3PWFS full frame')
        r=load('trireconstructormatrix.mat');
        rmatrix=r.rmatrix;
    end
    if tripyramid==false && MVM ==true
        fprintf('4PWFS MVM')
        r=load('MVMquadreconstructormatrix.mat');
        rmatrix=r.rmatrix;
    end
    if tripyramid==false && MVM ==false
        fprintf('4PWFS full frame')
        r=load('quadreconstructormatrix.mat');
        rmatrix=r.rmatrix;
    end
end
        
%%   Run the WFS
for   error= 10*10^-9
if runsim==true
    
%Returns variable rpuil, that is a data cube of the reconstructed wavefronts
%from PWFS. Returns Rwavefront which is the zernike polynomials used in the
%reconstruction.
rpupil=[];
%Calls on function pyramidsim

if kolmogorov==false
[rpupil, Rwavefront, Pupilpyramid, phase, mcount, ncount]=zernpupil(nn,mm, nphot,rdns, error,npix, Npix, sampling, rmatrix, lambda, tripyramid, MVM, pyramidmask);
end




if kolmogorov == true
[kscreen,kphase, kpupil]=kolpupilgen(nscreen,error,npix, Npix, lambda);

for k=1:nscreen  
[rpupil(:,:,k),Rwavefront(:,:,k), phase(:,:,k)]= pyramidsim(nn,mm,nphot,rdns,error,npix, Npix, sampling, rmatrix, lambda, tripyramid, MVM, kpupil(:,:,k), pyramidmask);
residual(:,:,k)=phase(:,:,k)-kphase(:,:,k);
end


end

%rpupil= pyramidsim(npix, Npix, sampling, rmatrix, lambda, sampling, tripyramid, MVM);
end
end
end
%% Run 3PWFS versus 4 PWfs script.

%Must have prebuilt the quadpyramidmask and tripyramidmask you want. Same
%with the reconstructor matrix.
if duelrun==true
    error= 10*10^-9;
%% Load in the masks and reconstructor matrixes
if broke==false
    load 'quadpyramidmask.mat';
    quadpyramidmask=pyramidmask;
    load 'tripyramidmask.mat';
    tripyramidmask=pyramidmask;
end

if broke==true
  load 'brokequadpyramidmask.mat';
    quadpyramidmask=pyramidmask;
    load 'broketripyramidmask.mat';
    tripyramidmask=pyramidmask;  
end



if have_reconstructor==false
tripyramid=true;
[trirmatrix,success, triSig, trippupil]=reconstructorgenerator(nn,mm,rdns, rnphot, npix,Npix, tripyramidmask, lambda, error, sampling, tripyramid, MVM, broke);
tripyramid=false; 
[quadrmatrix,success, quadSig, quadppupil]=reconstructorgenerator(nn,mm,rdns, rnphot, npix,Npix, quadpyramidmask, lambda, error, sampling, tripyramid, MVM, broke);
end

if have_reconstructor==true
 r=load('MVMquadreconstructormatrix.mat');
        quadrmatrix=r.rmatrix;
 r=load('MVMtrireconstructormatrix.mat');
        trirmatrix=r.rmatrix;    
end      
        
 %% Generate phase screens
for trial=1:10000
[kscreen,kphase, kpupil]=kolpupilgen(nscreen,error,npix, Npix, lambda,nphot);
mask=zernike(0,0,npix);
for k=nscreen
tripyramid=true;  
[trirpupil(:,:,k),triRwavefront(:,:,k), triphase(:,:,k),triPupilpyramid(:,:,k),trisums(:,:,k)]= pyramidsim(nn,mm,nphot,rdns,error,npix, Npix, sampling, trirmatrix, lambda, tripyramid, MVM, kpupil(:,:,k), tripyramidmask, broke);
triresidual(:,:,k)=triphase(:,:,k)-kscreen(:,:,k);
tripyramid=false; 
[quadrpupil(:,:,k),quadRwavefront(:,:,k), quadphase(:,:,k),quadPupilpyramid(:,:,k),quadsums(:,:,k)]= pyramidsim(nn,mm,nphot,rdns,error,npix, Npix, sampling, quadrmatrix, lambda, tripyramid, MVM, kpupil(:,:,k), quadpyramidmask, broke);




quadresidual(:,:,k)=quadphase(:,:,k)-kscreen(:,:,k);
triresid=triresidual(:,:,k);
quadresid=quadresidual(:,:,k);
trivariance(k)=var(triresid(mask>0));
quadvariance(k)=var(quadresid(mask>0));
end




%%
%figure; plot(1:nscreen, trivariance)
%hold on
%plot(1:nscreen,quadvariance); title('Varience vs Trail'); legend('3PWFS','4PWFS')

%% Calculate residuals

%Fourier residual.
[trispatialmode1, quadspatialmode1,xp,yp,cp]= fourierbasis(npix,nmodes,triresidual,quadresidual);
trispatialmode(:,:,trial)=trispatialmode1;
quadspatialmode(:,:,trial)=quadspatialmode1;
end

%%
quadspatialmode(isnan(quadspatialmode))=0;
trispatialmode(isnan(trispatialmode))=0;

%% Averaged residuals
trivarspatialmode=var(trispatialmode,0,3);
quadvarspatialmode=var(quadspatialmode,0,3);

%% Fixing m and n spatial frequency indexing

m=xp-cp;
n=yp-cp;
%% Calculate k
for i=1:size(m)
    k(i)=sqrt((m(i).^2+n(i).^2)/npix.^2);
end
%%
m=m+abs(min(m))+1;
n=n+abs(min(n))+1;


%% Quad Sine and Cosine matrix
for current=1:size(m)
    quadsinematrix(m(current),n(current))=quadvarspatialmode(current*2);
    quadcosinematrix(m(current),n(current))=quadvarspatialmode(current*2-1);
    totalquadvarspatialfreq(current)=sqrt(quadsinematrix(m(current),n(current)).^2+quadcosinematrix(m(current),n(current)).^2);
end

%figure; imagesc(-nmodes:0,-nmodes:nmodes,quadsinematrix); title('sine matrix 4PWFS'); colorbar; caxis([0, max(max(quadvarspatialmode))]); 
%figure; imagesc(-nmodes:0,-nmodes:nmodes,quadcosinematrix); title('cosine matrix 4PWFS'); colorbar; caxis([0, max(max(quadvarspatialmode))]); 

totalquadspatialfreqmodes=sqrt(quadsinematrix.^2+quadcosinematrix.^2);
figure; imagesc(-nmodes:0,-nmodes:nmodes,totalquadspatialfreqmodes);title('4PWFS Variance of fitted Fourier Amplitude'); colorbar; axis equal; caxis([0, max(max(quadvarspatialmode))]); xlim([-nmodes-0.5 0.5]); xlabel('n Fourier index'); ylabel('m Fourier index')

%% Sine and Cosine matrix
for current=1:size(m)
    sinematrix(m(current),n(current))=trivarspatialmode(current*2);
    cosinematrix(m(current),n(current))=trivarspatialmode(current*2-1);
    totaltrivarspatialfreq(current)=sqrt(sinematrix(m(current),n(current)).^2+cosinematrix(m(current),n(current)).^2);
end

%figure; imagesc(sinematrix); title('sine matrix'); colorbar; caxis([0, max(max(trivarspatialmode))]); axis equal
%figure; imagesc(cosinematrix); title('cosine matrix'); colorbar; caxis([0, max(max(trivarspatialmode))]); axis equal

totaltrispatialfreqmodes=sqrt(sinematrix.^2+cosinematrix.^2);
figure; imagesc(-nmodes:0,-nmodes:nmodes, totaltrispatialfreqmodes);title('3PWFS Variance of fitted Fourier Amplitude'); colorbar; axis equal; caxis([0, max(max(trivarspatialmode))]); xlim([-nmodes-0.5 0.5]); xlabel('n Fourier index'); ylabel('m Fourier index')


%% Variance vs Spatial Frequency Plots

uniqk=unique(k);
ssz=size(uniqk);
for i=1:ssz(2)
    I=find(k==uniqk(i));
    sz=size(I);
    for j=1:sz(2)
        values(j)=totaltrivarspatialfreq(I(j));
        quadvalues(j)=totalquadvarspatialfreq(I(j));
    end
    uniqtotaltrivarspatialfreq(i)=sqrt(sum(values.^2));
    uniqtotalquadvarspatialfreq(i)=sqrt(sum(quadvalues.^2));
end

%%
figure; plot(uniqk,abs(uniqtotaltrivarspatialfreq))%, 'LineWidth', 4);
hold on
plot(uniqk, abs(uniqtotalquadvarspatialfreq))%, 'LineWidth',4); %title('Variance of fitted Amplitude vs Spatial Frequency','fontsize',56); ; xlabel('Spatial Frequency (1/m)', 'fontsize',56), ylabel('Varience of the fitted Fourier Mode Amplitudes (rms radians)','fontsize',56);
legend({'3PWFS','4PWFS'},'fontsize',56)
xt = get(gca, 'XTick');
set(gca, 'FontSize', 50)
set(gca,'LineWidth',5)



sum(sum(abs(uniqtotaltrivarspatialfreq)))
sum(sum(abs(uniqtotalquadvarspatialfreq)))
sum(sum(abs(uniqtotaltrivarspatialfreq)/abs(uniqtotalquadvarspatialfreq)))



end 
 

%% Saving
basename='/Users/laurenschatz/Documents/MATLAB/Final Figures/3 Read Noise/2.5*10^6 photons/32';
trialnumbers='/10000trial';
%%
% filename=[basename,trialnumbers, 'triphase', num2str(error)];
% filename2=[basename,trialnumbers,  'triresidual', num2str(error)];
% filename3=[basename,trialnumbers,  'kscreen', num2str(error)];
% filename4=[basename,trialnumbers,  'quadphase', num2str(error)];
% filename5=[basename,trialnumbers,  'quadresidual', num2str(error)];
% filename6=[basename,trialnumbers,  'quadvarience', num2str(error)];
% filename7=[basename,trialnumbers,  'trivarience', num2str(error)];
%%
filename8=[basename,trialnumbers,  'trivarspatialmode', num2str(error)];
filename9=[basename,trialnumbers,  'quadvarspatialmode', num2str(error)];
filename10=[basename,trialnumbers,  'totalquadspatialfreqmodes', num2str(error)];
filename11=[basename,trialnumbers,  'totaltrispatialfreqmodes', num2str(error)];
filename12=[basename,trialnumbers,  'uniqtotaltrivarspatialfreq', num2str(error)];
filename13=[basename,trialnumbers,  'uniqtotalquadvarspatialfreq', num2str(error)];
filename14=[basename,trialnumbers,  'uniqk', num2str(error)];
%%

% save(filename,'triphase', '-v7.3');
% save(filename2,'triresidual', '-v7.3');
% save(filename3,'kscreen', '-v7.3');
% save(filename4,'quadphase', '-v7.3');
% save(filename5,'quadresidual', '-v7.3');
% save(filename6,'quadvarience', '-v7.3');
% save(filename7,'trivarience', '-v7.3');
%%
save(filename8,'trivarspatialmode', '-v7.3');
save(filename9,'quadvarspatialmode', '-v7.3');
save(filename10,'totalquadspatialfreqmodes', '-v7.3');
save(filename11,'totaltrispatialfreqmodes', '-v7.3');
save(filename12,'uniqtotaltrivarspatialfreq', '-v7.3');
save(filename13,'uniqtotalquadvarspatialfreq', '-v7.3');
save(filename14,'uniqk', '-v7.3');

% %%
% filenamefits=[filename5,'.fits'];
% filenamefits2=[filename2,'.fits'];
% filenamefits3=[filename3,'.fits'];
% 
% fitswrite(quadresidual,filenamefits);
% fitswrite(triresidual,filenamefits2);
% fitswrite(kscreen,filenamefits3);

%%
% filename=[basename, '/trivarspatialfreq1200'];
% filename2=[basename, '/quadvarspatialfreq1200'];
% filename3=[basename, '/tricosinematrix1200'];
% filename4=[basename, '/trisinematrix1200'];
% filename5=[basename, '/quadcosinematrix1200'];
% filename6=[basename, '/quadsinematrix1200'];
% 
% save(filename,'trispatialfreq', '-v7.3');
% save(filename2,'quadvarspatialfreq', '-v7.3');
% save(filename,'trivarspatialfreq', '-v7.3');
% save(filename3,'cosinematrix', '-v7.3');
% save(filename4,'sinematrix', '-v7.3');
% save(filename5,'quadcosinematrix', '-v7.3');
% save(filename6,'quadsinematrix', '-v7.3');
    

    
    
    