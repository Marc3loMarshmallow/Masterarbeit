%  Make the image interpolation for the polar scan

function interp_n_plot(params, new_dir)

%  Set initial parameters

addpath(genpath('../kMDEV'));
addpath(new_dir);

%  Read the data and adjust it in time 

min_sample=0;
if params.noise == 1
    for i=1:params.no_lines
        %  Load the result

        %cmd=['load ', new_dir, 'rf_ln', num2str(i), '.mat']
        cmd=['load 226950_Series/Sim2_', num2str(t), '_', num2str(ind-1), '/rf_ln',num2str(i),'.mat']

        eval(cmd)

        res_rf = params.c / (2*params.fs);
   
        d = res_rf*(1:length(rf_data));
        d(1)=d(2)/2;
        amplific = exp(params.mu*d);
        noise = params.sigma * randn(1,length(rf_data)) .* amplific;
        new_rf = rf_data + noise';

        %  Find the envelope
  
        rf_env=abs(hilbert([zeros(round(tstart*params.fs-min_sample),1); new_rf]));

        env(1:max(size(rf_env)),i)=rf_env;

    end
else
    for i=1:params.no_lines

        %  Load the result

        cmd=['load ', new_dir, 'rf_ln', num2str(i), '.mat']

        eval(cmd)

        %  Find the envelope
  
        rf_env=abs(hilbert([zeros(round(tstart*params.fs-min_sample),1); rf_data]));

        env(1:max(size(rf_env)),i)=rf_env;
    end
end

%  Do logarithmic compression to 40 dB

dB_Range=50;
env=env-min(min(env));

log_env=20*log10(env(1:params.D:max(size(env)),:)/max(max(env)));
log_env=255/dB_Range*(log_env+dB_Range);

%  Get the data into the proper format

start_depth=0.02;   % Depth for start of image in meters
image_size=0.105;   % Size of image in meters
skipped_samples=0;
samples=max(size(log_env));

start_of_data=(skipped_samples+1)/params.fs*params.c/2;           % Depth for start of data in meters
end_of_data=(skipped_samples+samples+1)/params.fs*params.c/2;     % Depth for end of data in meters
delta_r=params.c/2*1/params.fs;                                   % Sampling interval for data in meters

theta_start= -params.no_lines/2*params.dtheta;     % Angle for first line in image

Nz=512;                         % Size of image in pixels
Nx=512;                         % Size of image in pixels
scaling=1;                      % Scaling factor form envelope to image

[N,M]=size(log_env);
D=floor(N/1024);

env_disp=uint8(255*log_env(1:D:floor(params.radius*N),:)/max(max(log_env)));
%env_disp=uint8(255*log_env(1:D:N,:)/max(max(log_env)));

img_data=env_disp;

values.probe=28;
values.sectorPercent = 100;       
values.resolutionDisplacement = params.c/2/params.fs; %in m
values.depthOffset = 0; %in m
parameters.resolution = [0.5 0.5]*10^-4;
[imageOut,mask, values] = scanconvert(double(env_disp), values, parameters);

% Set accurate frame

New_img=imresize(imageOut, (342/size(imageOut,1)));
frame = zeros(158, size(New_img, 2));
img_frame = cat(1, frame, New_img);
frame = zeros(100, size(New_img, 2));
img_frame = cat(1, img_frame, frame);
frame = zeros(size(img_frame, 1), 58);
img_frame = cat(2, frame, img_frame);
img_frame = cat(2, img_frame, frame);

% Speckle reduction filter

if params.speckle == 1
    img_frame = speckle(bild_mit_rand, [3 3], 2);

iw(img_frame)


