% This program generates a 2D US-image for a Transesophageal Echocardiography (TEE)
% in 3 steps based on a 3D-phantom as input (preferably in Nifti format)
% and correspondingly set parameters



% STEP 1: Preprocessing - Generating the scatter Coordinates and Amplitudes using
%         the input-phantom as a map


%   Load the input data

t = 2;      % Timepoint of the valve closure (1-30)
rot = 7;    % Rotation angle of the TEE (1-17)

input_nifti = ['226950-timeseries/226950-TorsoMask_Rot',num2str(rot),'.nii.gz']
% input_nifti = ['293182-timeseries/293182-TorsoMask_Rot',num2str(rot),'.nii.gz']
% input_nifti = ['771083-timeseries-NoPM/771083-TorsoMask_Rot',num2str(rot),'.nii.gz']

img.vol     = niftiread(input_nifti);   % Loading the data from the file

%   Cut unnecessary dimensions, crop and rotate phantoms if needed

img.vol     = img.vol(:, :, :, 1, t);             % 5 dim to 3 dim image
img.vol     = img.vol(81:535, 206:547, 12:28);    % For valve 226950
% img.vol     = img.vol(91:545, 206:547, 12:28);  % For valve 293182
% img.vol     = img.vol(81:535, 206:547, 12:28);  % For valve 771083
img.vol     = imrotate(img.vol, -90);

%   Set the phantom properties

img.px_size = [size(img.vol,2) size(img.vol,1)];       % Image size in Pixels
img.mm_size = [189.6 , 16.4 , 142.6];                  % Image size in mm
img.n_sc    = 1e5/2;                                   % Number of scatters

%   Use the make_sc() function to generate coordinates and amplitudes for the scatters

[phantom_positions, phantom_amplitudes] = make_sc(img);

%   Save the data

new_pht = 'pht_data_test.mat'

save(new_pht, 'phantom_positions', 'phantom_amplitudes')



% STEP 2: Simulation - Simulating the US-data with Field II based on the generated scatters


%   Initialize the Field II system

field_init()

%   Set the properties of the transducer you want to simulate with

trans.f0       = 7e6;                  %  Transducer center frequency [Hz]
trans.fs       = 40e6;                 %  Sampling frequency [Hz] 
trans.c        = 1540;                 %  Speed of sound [m/s]
trans.lambda   = trans.c/trans.f0;     %  Wavelength [m]
trans.width    = trans.lambda/2;       %  Width of element
trans.el_h     = 5/1000;               %  Height of element [m]
trans.kerf     = trans.lambda/10;      %  Kerf [m]
trans.focus    = [0 0 70]/1000;        %  Fixed focal point [m]
trans.N_el     = 64;                   %  Number of physical elements

%   Choose the name of the folder you want the US data to be saved in

new_dir = ['Test/'];

%   Use the run_sim() function with Field II functions to run the simulation

run_sim(trans, new_pht, new_dir);



% STEP 3: Postprocessing - Adding Attenuation, Noise and Speckle reduction
%         (Putting together the US data to plot the desired 2D TEE image)


%   Set the properties for the US image

params.D           = 4;                                    %  Sampling frequency decimation factor
params.fs          = 40e6/params.D;                        %  Sampling frequency  [Hz]
params.c           = 1540;                                 %  Speed of sound [m/s]
params.no_lines    = 64;                                   %  Number of lines in image
params.image_width = 90/180*pi;                            %  Size of image sector [rad]
params.dtheta      = params.image_width/params.no_lines;   %  Increment for image
params.radius      = 0.8;                                    %  normalised US cone radius (1 for full radius)

%   Enable Gaussian noise in the image and set the parameters

params.noise       = 1;                      %  Enable/ disable noise (1/0) 
params.mu          = 7*1.15;
params.sigma       = 1e-27;                  %  Noise increase with radius

%   Enable speckle reduction

params.speckle     = 1;

%   Use the interp_n_plot() to generate the image

interp_n_plot(params, new_dir)

