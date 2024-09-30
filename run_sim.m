% This function uses Field II to simulate the US-Data based on
% given scatter Coordinates and Amplitudes with set Transducer

function run_sim(trans, new_pht, new_dir)

%  Initialize the Field II system

% addpath('Field_II_ver_3_30_mac')  % Add Field II folder if needed
field_init()

%  Generate the transducer apertures for send and receive

%  Set the sampling frequency

set_sampling(trans.fs);
set_field ('show_times', 5)

%  Generate aperture for emission

xmit_aperture = xdc_linear_array (trans.N_el, trans.width, trans.el_h, trans.kerf, 1, 5, trans.focus);

%  Set the impulse response and excitation of the xmit aperture

impulse_response=sin(2*pi*trans.f0*(0:1/trans.fs:2/trans.f0));
impulse_response=impulse_response.*hanning(max(size(impulse_response)))';
xdc_impulse (xmit_aperture, impulse_response);

excitation=sin(2*pi*trans.f0*(0:1/trans.fs:2/trans.f0));
xdc_excitation (xmit_aperture, excitation);

%  Generate aperture for reception

receive_aperture = xdc_linear_array (trans.N_el, trans.width, trans.el_h, trans.kerf, 1, 5, trans.focus);

%  Set the impulse response for the receive aperture

xdc_impulse (receive_aperture, impulse_response);

%  Load the scatter data and make a directory for the rf data

load(new_pht)
mkdir(new_dir)

%  Set the different focal zones for reception

focal_zones=[5:1:150]'/1000;
Nf=max(size(focal_zones));
focus_times=(focal_zones-10/1000)/1540;
z_focus=60/1000;          %  Transmit focus

%  Set the apodization

apo=hanning(trans.N_el)';
xdc_apodization (xmit_aperture, 0, apo);
xdc_apodization (receive_aperture, 0, apo);

%  Do phased array imaging

no_lines=64;                   %  Number of lines in image
image_width=90/180*pi;         %  Size of image sector [rad]
dtheta=image_width/no_lines;   %  Increment for image

%  Do imaging line by line

for i=1:64

  if ~exist([new_dir, 'rf_ln',num2str(i),'.mat'])
    
    cmd=['save ', new_dir, 'rf_ln',num2str(i),'.mat i']
    eval(cmd)

    %   Set the focus for this direction

    theta= (i-1-no_lines/2)*dtheta;
    xdc_focus (xmit_aperture, 0, [z_focus*sin(theta) 0 z_focus*cos(theta)]);
    xdc_focus (receive_aperture, focus_times, [focal_zones*sin(theta) zeros(max(size(focal_zones)),1) focal_zones*cos(theta)]);
  
    %   Calculate the received response

    [rf_data, tstart]=calc_scat(xmit_aperture, receive_aperture, phantom_positions, phantom_amplitudes);
    
    %   Store the result
    
    cmd=['save ', new_dir, 'rf_ln',num2str(i),'.mat rf_data tstart'];
    eval(cmd)
    end

  end

%   Free space for apertures

xdc_free (xmit_aperture)
xdc_free (receive_aperture)

