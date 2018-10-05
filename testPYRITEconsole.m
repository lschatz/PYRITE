% PYRITE CONSOLE


clear all, close all;

%% Simulation Setup:
% In this section are the system parameters for the simulator.


error= 10*10^-9; %Amount of wavefront error in nanometers applied in the pupil.
    npix= 256; %number of pixels across the pupil diameter
        Npix= 2560; % number of pixels across the total roster. Represents the amount of zero padding.
            lambda= 700.*10^-9; %wavelength in nanometers

nphot=2.5*10^6; % The total number of photons on the detector (from summing all the photons in detection), before noise.
    rnphot=10^15; % The total number of photons on the detector used in building the reconstructor matrix.
        rdns= 0; % Read noise applied

sampling=32; % 128; 32 %Number of pixels across each pupil. Can be 16, 32, 64
    nmodes=4; %Number of fourier modes used to fit the residual WF error
        %nn=11; %Zernike order Z_nm used in reconstructor matrix calibration.
            %mm=11; % Zernike order. (Will be +/-m)

trial=2; %number of trials of the simultion to run
nscreen=trial; %number of unique kolmogorov phase screens to generate. 


%%

runsim=true; %toggle true/false to run complete sim. If false will only generate reconstruction matrix.

tripyramid=true; %toggle true/false for 3PWFS simulation. If false, will run 4PWFS simulation
    have_reconstructor=false; %toggle true/false for reconstructor matrix
        have_pyramidmask=true;  %toggle true/fasle for pyramid focal plane mask.
            have_fourierbasis=false; %toggle true/false for fourier basis set
  

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
[fbpupil,fbs]= testfbspupil(npix, nmodes, Npix, lambda, error, rnphot); %Generate the Fourier Basis set
end

if have_fourierbasis==true
    fprintf('loading Fourier Basis Set')
    
    fbs= load('fbs.mat');
    fbpupil=load('fbpupil.mat');
end
    
    
    
    
%Generate Kologorov Phase Screens.
fprintf('Calculating Kolmogorov Phase Screens')
for i= 1:nscreen
[kscreen(:,:,i), kpupil(:,:,i)]=testkolpupilgen(error,npix, Npix, lambda,nphot);
end

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
        [trirmatrix, triSig, trippupil]=testreconstructorgeneratortest(fbpupil, rdns, npix, Npix, tripyramidmask, sampling, tripyramid, MVM, broke);
            tripyramid=false; 
            [quadrmatrix, quadSig, quadppupil]=testreconstructorgeneratortest(fbpupil, rdns, npix, Npix, quadpyramidmask, sampling, tripyramid, MVM, broke);
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
    [trirpupil,triRwavefront, triphase]= testpyramidsim(lambda, error, fbs, trirmatrix,kpupil,kscreen, rdns, npix, Npix, tripyramidmask, sampling, tripyramid, MVM, broke);


tripyramid=false;
fprintf('Running 4PWFS')
% Four sided
%% Run 4PWFS
    [quadrpupil,quadRwavefront, quadphase]= testpyramidsim(lambda, error, fbs, quadrmatrix,kpupil,kscreen, rdns, npix, Npix, quadpyramidmask, sampling, tripyramid, MVM, broke);


% %% Run Analysis
% 
% %Step One: Calculate variance of residual wavefront across pupil.
% 
% [trivariance, quadvariance, quadresidual, triresidual]=residualcalc(triphase, quadphase, kscreen,npix);
% 
% %Step Two: Fit Fourier Modes to the residual wavefront.
% 
% %[]=fourierresidfit(quadresidual, triresidual, fbs);

