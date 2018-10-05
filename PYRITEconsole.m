% PYRITE CONSOLE


clear all, close all;

%% Simulation Setup:
% In this section are the system parameters for the simulator.


error= 10*10^-9; %Amount of wavefront error in nanometers applied in the pupil.
    npix= 256; %number of pixels across the pupil diameter
        Npix= 2560; % number of pixels across the total roster. Represents the amount of zero padding.
            lambda= 700.*10^-9; %wavelength in nanometers

nphot=15*10^8; % The total number of photons on the detector (from summing all the photons in detection), before noise.
    rnphot=10^15; % The total number of photons on the detector used in building the reconstructor matrix.
        rdns= 0; % Read noise applied

sampling=32; % 128; 32 %Number of pixels across each pupil. Can be 16, 32, 64
    nmodes=6; %Number of fourier modes used to fit the residual WF error
        nn=6; %Zernike order Z_nm used in reconstructor matrix calibration.
            mm=6; % Zernike order. (Will be +/-m)

trial=100; %number of trials of the simultion to run
nscreen=trial; %number of unique kolmogorov phase screens to generate. 


%%

runsim=true; %toggle true/false to run complete sim. If false will only generate reconstruction matrix.

tripyramid=true; %toggle true/false for 3PWFS simulation. If false, will run 4PWFS simulation
    have_reconstructor=true; %toggle true/false for reconstructor matrix
        have_pyramidmask=true;  %toggle true/fasle for pyramid focal plane mask.
            have_fourierbasis=true; %toggle true/false for fourier basis set
                add_noise=false; %toggle true/false for read noise and photon noise.

MVM=true; %Toggle true/false. True: Will run Matrix-Vector-Multiply (intensity Centroids) for reconstructor matrix
                             %Flase: Will use the full frame (Olivier
                             %Guyon's method) for the reconstructor matrix
kolmogorov=true;


roof=false; %toggle true/false for roof error
    rooferror=0;
tabletop=false; % toggle true/false for tabletop error
    tabletoperror=0;
broke=false;
%% Start Loading/ Generating

% Generate Fourier Basis set and pupils
if have_fourierbasis==false
    fprintf('Generating Fourier Basis Set');
[fbpupil,fbs,xp,yp,cp]= fbspupil(npix, nmodes, lambda, error, rnphot); %Generate the Fourier Basis set
%[fbpupil,fbs]=zbasis(nn, mm, npix, Npix, nmodes,lambda,error);
end
%%

if have_fourierbasis==true
    fprintf('loading Fourier Basis Set')
    
    load('fbs.mat');
    load('fbpupil.mat');
    load('xp.mat');
    load('yp.mat');
    load('cp.mat');

end

%%
z=fitsread('modf_256_15_ortho.fits');
son=size(z);
onfbs=zeros(npix,npix,son(3)+2);
mask=zernike(0,0,npix);
[onfbs(:,:,1),onfbs(:,:,2)] = meshgrid([-npix/2+1:npix/2]/(npix/2),[-npix/2+1:npix/2]/(npix/2));
onfbs(:,:,3:son(3)+2)=z;
onfbs=onfbs.*mask;

for i=1:son(3)
onfbs(:,:,i)=onfbs(:,:,i).*sqrt(1/sum(sum(onfbs(:,:,i).*onfbs(:,:,i))));% normalize
end


onfbpupil=zernike(0,0,npix).*exp(((2*pi*1i)/lambda)*error.*onfbs);
    
    
%%
    
%Generate Kologorov Phase Screens.
fprintf('Calculating Kolmogorov Phase Screens')
for i= 1:nscreen
[kscreen(:,:,i), kpupil(:,:,i)]=kolpupilgen(error,npix, Npix, lambda,nphot);
end


%kpupil=fbpupil;

% Calculate the RMS of the residual Fourier Mode Amplitudes of the Kolmogorov phase screen





%Generate the pyramid focal plane phase masks
if have_pyramidmask==false
    fprintf('generating pyramid masks')
        tripyramid=true;
        tripyramidmask=maskgenerator(Npix, tripyramid, tabletoperror, tabletop, roof, rooferror);
            tripyramid=false;
            quadpyramidmask=maskgenerator(Npix, tripyramid, tabletoperror, tabletop, roof, rooferror);  
end
    
%Load in saved pyramid focal plane phase masks
if have_pyramidmask==true 
    fprintf('Loading pyramid masks')
    load 'quadpyramidmask.mat';
    quadpyramidmask=pyramidmask;
        load 'tripyramidmask.mat';
        tripyramidmask=pyramidmask;
end


%Generate reconstructor martix


if have_reconstructor==false
    fprintf('Generating Reconstructor Matrices')
        tripyramid=true;
        [trirmatrix, triSig, trippupil]=reconstructorgenerator(add_noise,rnphot,onfbpupil, rdns, npix, Npix, tripyramidmask, sampling, tripyramid, MVM, broke);
            tripyramid=false; 
            [quadrmatrix, quadSig, quadppupil]=reconstructorgenerator(add_noise,rnphot,onfbpupil, rdns, npix, Npix, quadpyramidmask, sampling, tripyramid, MVM, broke);
end

%Load in saved reconstructor matrix
fprintf('Loading Reconstructor Matrices')
if have_reconstructor==true
    r=load('MVMquadreconstructormatrix.mat');
    quadrmatrix=r.rmatrix;
        r=load('MVMtrireconstructormatrix.mat');
        trirmatrix=r.rmatrix;    
end 



%% Run Sims

%Variable Outputs:
    %tri/quadpupil: the pupil with the reconstructed phase.
        %quadpupil=exp(i*2p/lamda * reconstructed phase)
    %tri/quadRwavefront= reconsturcted wavefront (unmasked)
    %tri/quadphase= reconstructed phase 



fprintf('Running 3PWFS')
% Three Sided
tripyramid=true;                                    
    [trirpupil,triRwavefront, triphase]= pyramidsim(add_noise,nphot,lambda, error, onfbs, trirmatrix, kpupil, rdns, npix, Npix, tripyramidmask, sampling, tripyramid, MVM, broke);


tripyramid=false;
fprintf('Running 4PWFS')
% Four sided
%% Run 4PWFS
    [quadrpupil,quadRwavefront, quadphase]= pyramidsim(add_noise,nphot,lambda, error, onfbs, quadrmatrix,kpupil, rdns, npix, Npix, quadpyramidmask, sampling, tripyramid, MVM, broke);


%% Run Analysis

%Step One: Calculate variance of residual wavefront across pupil.

[trivariance, quadvariance, quadresidual, triresidual]=residualcalc(triphase, quadphase, kscreen,npix);

%%Step Two: Fit Fourier modes 
[trival]=fourierfit(fbs,triphase);
[quadval]=fourierfit(fbs,quadphase);
[Rtrival]=fourierfit(fbs,triresidual);
[Rquadval]=fourierfit(fbs,quadresidual);
[kval]=fourierfit(fbs,kscreen); 




%% Step 3: RMS across trials

%RMS

[tri]=rms(trival,2);
[quad]=rms(quadval,2);
[kspat]=rms(kval,2);
[Rtri]=rms(Rtrival,2);
[Rquad]=rms(Rquadval,2);


%% Step 4: Fit Spatial Frequencies
%Fit spatial Frequencies
[tri]=spatialfreqcalc(tri,xp,yp,cp,nmodes);
[quad]=spatialfreqcalc(quad,xp,yp,cp,nmodes);
[kspat]=spatialfreqcalc(kspat,xp,yp,cp,nmodes);
[Rtri]=spatialfreqcalc(Rtri,xp,yp,cp,nmodes);
[Rquad]=spatialfreqcalc(Rquad,xp,yp,cp,nmodes);


%% Step 5: Normalize the Residuals to get percent error
[nRtrival]=valnorm(Rtri,kspat);
[nRquadval]=valnorm(Rquad,kspat);



%% Step 6: Plot Power Spectrums: Orginal, Reconstructed, Residual vs Cycles/Aperture
figure; plot(tri)
hold on
plot(quad); plot(kspat);plot(Rtri); plot(Rquad); %plot(nRquadval), plot(nRtrival);
legend({'3PWFS','4PWFS','Original','Residual 3', 'Residual 4', 'Normalized Resid 4', 'Normalized Resid 3'})

% log scale
figure; loglog(tri)
hold on
loglog(quad); loglog(kspat);loglog(Rtri); loglog(Rquad); %plot(nRquadval), plot(nRtrival);
legend({'3PWFS','4PWFS','Original','Residual 3', 'Residual 4', 'Normalized Resid 4', 'Normalized Resid 3'})
%% Step 7: Plot Percent Error vs Cycles Per aperture
figure; plot(nRquadval); hold on
 plot(nRtrival);
legend({'4PWFS, 3PWFS'});

%log scale
figure; loglog(nRquadval); hold on
 loglog(nRtrival);
legend({'4PWFS, 3PWFS'});



%%
% %% Step Four: That the RMS across trial numbers.
% 
% totnRtrival=rms(nRtrival,2);
% totnRquadval=rms(nRquadval,2);
% 
% %% Step Five: Calculate spatial frequencies, combine cosine/sine terms
%  
% [trispatialfreq,uniqk,tritotalvar]=spatialfreqcalc(totnRtrival,xp,yp,cp,nmodes);
% [quadspatialfreq,uniqk, quadtotalvar]=spatialfreqcalc(totnRquadval,xp,yp,cp,nmodes);
% 
% %% Step Six: Plot percent error vs cycles per aperture
% 
% figure; loglog(uniqk, trispatialfreq)%, 'LineWidth', 4);
% hold on
% loglog(uniqk, quadspatialfreq)%, 'LineWidth',4); %title('Variance of fitted Amplitude vs Spatial Frequency','fontsize',56); ; xlabel('Spatial Frequency (1/m)', 'fontsize',56), ylabel('Varience of the fitted Fourier Mode Amplitudes (rms radians)','fontsize',56);
% legend({'3PWFS','4PWFS'},'fontsize',56)
% xt = get(gca, 'XTick');
% set(gca, 'FontSize', 50)
% set(gca,'LineWidth',5)
% 
% 
% figure; plot(uniqk, trispatialfreq)%, 'LineWidth', 4);
% hold on
% plot(uniqk, quadspatialfreq)%, 'LineWidth',4); %title('Variance of fitted Amplitude vs Spatial Frequency','fontsize',56); ; xlabel('Spatial Frequency (1/m)', 'fontsize',56), ylabel('Varience of the fitted Fourier Mode Amplitudes (rms radians)','fontsize',56);
% legend({'3PWFS','4PWFS'},'fontsize',56)
% xt = get(gca, 'XTick');
% set(gca, 'FontSize', 50)
% set(gca,'LineWidth',5)
% 
% 
% sum(sum(abs(trispatialfreq)))
% sum(sum(abs(quadspatialfreq)))
% sum(sum(abs(trispatialfreq)/abs(quadspatialfreq)))
% 
% 
