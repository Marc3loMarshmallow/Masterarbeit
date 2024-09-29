% Creates a phantom for a Heart from a CT scan of the Heart.
% The phantom starts 2 mm from the transducer surface.

function [positions, amp] = make_sc(img)

% Define image coordinates

Nl = img.px_size(1) 
Ml = img.px_size(2);

x_size = img.mm_size(1) / 1000 ;
y_size = img.mm_size(2) / 1000 ;
z_size = img.mm_size(3) / 1000 ;

dx=x_size/Nl;          %  Sampling interval in x direction [m]
dz=z_size/Ml;          %  Sampling interval in z direction [m]
dy=y_size/17;          %  Sampling interval in y direction [m]
theta = 35/180*pi;     %   Rotation of the scatterers [rad]
theta = 0;
z_start = 2/1000;

% Calculate position data

z0 = rand(img.n_sc,1);

x0 = rand(img.n_sc,1);
x = (x0-0.5)* x_size;

z = z0*z_size+z_start;

y0 = rand(img.n_sc,1);
y = (y0-0.5)* y_size; 

amp = zeros(img.n_sc,1);

% If you want to increase the number of scatters in the valve

% [z1, x1, y1] = ind2sub(size(V(:,:,:)), find(V(:,:,:) == 107));
% x1=x1/Nl;
% z1=z1/Ml;
% y1=y1/17;
% numPixels = length(x1);
% randomIndices = randperm(numPixels, 3000);
% x0(1:3000) = x1(randomIndices);
% z0(1:3000) = z1(randomIndices);
% y0(1:3000) = y1(randomIndices);

for i=1:17
    map = img.vol(:,:,i,:,:);
    
    mask = map;                     % create mask for reflections
    mask=my_changem(mask, 1 , 0);

    % Set intensities outside the cone to 0
    
    for j=1:floor(Nl/2)
        for k=1:j
            map(floor(Nl/2)+1-j,k)=0;
            map(floor(Nl/2)+1-j,Nl+2-k)=0;
            mask(floor(Nl/2)+1-j,k)=0;
            mask(floor(Nl/2)+1-j,Nl+2-k)=0;
        end
    end

    % Create mask for reflections

    mask=my_changem(mask, 0 , 74);
    mask=my_changem(mask, 0 , 75);
    mask=my_changem(mask, 0 , 76);
    mask=my_changem(mask, 0 , 107);
    mask(mask~=0)=1;
    mask=logical(mask);

    % Set intensities / gray values

    map=my_changem(map, 100 , 106);
    map=my_changem(map, 100 , 105);
    map=my_changem(map, 0, 76);
    map=my_changem(map, 0 , 56);
    map=my_changem(map, 200 , 74);
    map=my_changem(map, 5 , 75);
    map=my_changem(map, 210 , 107);     % set to 30 if scatters increased
    map=my_changem(map, 100 , 61);
    map=my_changem(map, 120 , 73);
    map=my_changem(map, 100 , 65);
    map=my_changem(map, 100 , 67);
    map=my_changem(map, 40 , 59);
    map=my_changem(map, 100 , 55);
    map=my_changem(map, 100 , 57);        
    map=my_changem(map, 100 , 83);          
    map=my_changem(map, 100 , 81);          
    map=my_changem(map, 120 , 84);
    map=double(map);

    % Gradients and Impedance in axial direction

    kernel=[-1 0 1;-2 0 2;-1 0 1]/4;
    Gx = conv2(map,kernel, 'same');
    Gz = conv2(map,kernel', 'same');
    x_g=1:size(map, 2);
    [z_k,x_k] = meshgrid(x_g-mean(x_g), 1:size(map, 1));

    gtheta = cart2pol(z_k,x_k);

    unit_z=sin(gtheta);
    unit_x=cos(gtheta);

    g=unit_x.*Gx+unit_z.*Gz;
    g(:,1)=0;
    g(1,:)=0;
    g(:,end)=0;
    g(end,:)=0;
    edges=abs(g);

    Gx = conv2(mask,kernel, 'same');
    Gz = conv2(mask,kernel', 'same');
    g=unit_x.*Gx+unit_z.*Gz;
    g(:,1)=0;
    g(1,:)=0;
    g(:,end)=0;
    g(end,:)=0;
    edges_mask = abs(g);

    r=sqrt(z_k.^2 + x_k.^2);
    r=r/max(r(:));
    edge_kid = imdilate(edges, strel('disk',1));
    edge_mask = imdilate(edges_mask, strel('disk',1));
    edge_kid2 = imdilate(edge_mask.*mask.*r, strel('disk',20));
    edge_kid = edge_kid + 400*mask.*edge_kid2;

    map=map+edge_kid;

    map=map';

    %  Find the index for the amplitude value
    
    xindex = floor((x + 0.5*x_size)/dx + 1);
    zindex = floor((z - z_start)/dz + 1);
    yindex = floor((y + 0.5*y_size)/dy + 1);
    inside = (0 < xindex)  & (xindex <= Nl) & (0 < zindex)  & (zindex <= Ml) & (yindex==i);
    index = (xindex + (zindex-1)*Nl).*inside + 1*(1-inside);

    % TEMP
    linear_index = zeros(img.n_sc, 1);
    for l = 1:img.n_sc
        linear_index(l) = sub2ind(size(map), xindex(l), zindex(l));
    end
    amp_sl = map(linear_index);
    amp_sl(~inside) = 0;
    amp_sl(amp_sl~=0) = exp(amp_sl(amp_sl~=0)/100);

    % Amplitudes with different variance must be generated according to the 
    % input map.
    % The amplitude of the map is used to scale the variance

    amp_sl=amp_sl-min(min(amp_sl));
    amp_sl=amp_sl/max(max(amp_sl));

    amp = amp + amp_sl;           % The amplitudes of the scatters in different slices get stacked
    
end

% Generate the rotated and offset block of sample

xnew=x*cos(theta)+z*sin(theta);
znew=z*cos(theta)-x*sin(theta);
znew=znew-min(min(znew)) + z_start;

positions=[(xnew-40/1000) y znew];

% remove scatters wie amplitude zero for runtime purposes

xnew(amp==0)=[];
y(amp==0)=[];
znew(amp==0)=[];
amp(amp==0)=[];

positions=[xnew y znew];
