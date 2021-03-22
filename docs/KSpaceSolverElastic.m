clearvars;

% =========================================================================
% SIMULATION
% =========================================================================

% create the computational grid
Nx = 128;  % number of grid points in the x (row) direction
Ny = 256;  % number of grid points in the y (column) direction
dx = 0.1e-3;            % grid point spacing in the x direction [m]
dy = 0.1e-3;            % grid point spacing in the y direction [m]


% define the properties of the upper layer of the propagation medium
% cs1_magnitue = 1500;
% medium.sound_speed_compression = cs1_magnitue * ones(Nx, Ny);   % [m/s]
% medium.sound_speed_shear       = 500 * ones(Nx, Ny);    % [m/s]
% medium.density                 = 1000 * ones(Nx, Ny);   % [kg/m^3]
% 
% % define the properties of the lower layer of the propagation medium
% cs2_magnitue = 2000;
% medium.sound_speed_compression(Nx/2:end, :) = cs2_magnitue;     % [m/s]
% medium.sound_speed_shear(Nx/2:end, :)       = 800;      % [m/s]
% medium.density(Nx/2:end, :)                 = 1200;     % [kg/m^3]

cs1_magnitue = 1500;
medium.sound_speed_compression = cs1_magnitue * ones(Nx, Ny);   % [m/s]
medium.sound_speed_shear       = 500 * ones(Nx, Ny);    % [m/s]
medium.density                 = 1000 * ones(Nx, Ny);   % [kg/m^3]


% create initial pressure distribution using makeDisc
disc_magnitude = 5;         % [Pa]
disc_x_pos = 50;            % [grid points]
disc_y_pos = 50;            % [grid points]
disc_radius = 8;            % [grid points]
source.p0 = disc_magnitude * makeDisc(Nx, Ny, disc_x_pos, disc_y_pos, disc_radius);

% define a binary line sensor
sensor.mask = zeros(Nx, Ny);
sensor.mask(21, :) = 1;

kgrid.dx = dx;
kgrid.dy = dy;
kgrid.Nx = Nx;
kgrid.Ny = Ny;
kgrid.dim = 2;

% assign the grid parameters for the x and y spatial directions
kgrid.kx_vec = makeDim(kgrid.Nx, kgrid.dx);
kgrid.ky_vec = makeDim(kgrid.Ny, kgrid.dy);

% define the wavenumber based on the wavenumber components
% (using bsxfun saves memory by avoiding ndgrid)
kgrid.k = zeros(kgrid.Nx, kgrid.Ny);
kgrid.k = bsxfun(@plus, (reshape(kgrid.kx_vec, [], 1, 1).^2), kgrid.k);
kgrid.k = bsxfun(@plus, (reshape(kgrid.ky_vec, 1, [], 1).^2), kgrid.k);
kgrid.k = sqrt(kgrid.k);
                    
% define maximum supported frequency
kgrid.kx_max = max(abs(kgrid.kx_vec(:)));
kgrid.ky_max = max(abs(kgrid.ky_vec(:)));
kgrid.k_max  = min([kgrid.kx_max, kgrid.ky_max]); 

kgrid.x_size = kgrid.dx * kgrid.Nx;
kgrid.y_size = kgrid.dy * kgrid.Ny;

kgrid.kx = repmat(kgrid.kx_vec, [1, kgrid.Ny]);
kgrid.ky = repmat(kgrid.ky_vec.', [kgrid.Nx, 1]);

kgrid.x = kgrid.x_size .* kgrid.kx .* kgrid.dx ./ (2 .* pi);
kgrid.y = kgrid.y_size .* kgrid.ky .* kgrid.dy ./ (2 .* pi);

kgrid.x_vec = kgrid.x_size .* kgrid.kx_vec .* kgrid.dx ./ (2 .* pi);
kgrid.y_vec = kgrid.y_size .* kgrid.ky_vec .* kgrid.dy ./ (2 .* pi);


% if c is a matrix, find the minimum and maximum values
KSPACE_CFL = 0.3;                  % default CFL value used if kgrid.t_array is set to 'auto'
c = medium.sound_speed_compression;
c_max = max(c(:));
c_min = min(c(:));
t_end = sqrt(kgrid.x_size.^2 + kgrid.y_size.^2) ./ c_min;

min_grid_dim = min([kgrid.dx, kgrid.dy]);

% assign time step based on CFL stability criterion
cfl = KSPACE_CFL;
kgrid.dt = cfl .* min_grid_dim ./ c_max;
            
% assign number of time steps based on t_end
kgrid.Nt = floor(t_end / kgrid.dt) + 1;

% catch case were dt is a recurring number
if (floor(t_end / kgrid.dt) ~= ceil(t_end / kgrid.dt)) && (rem(t_end, kgrid.dt) == 0)
    kgrid.Nt = kgrid.Nt + 1; 
end

rho0 = medium.density;


c_ref = max(medium.sound_speed_compression(:));
pml_x_size = 20;
pml_y_size = 20;
pml_x_alpha = 2;
pml_y_alpha = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% assign the lame parameters
mu     = medium.sound_speed_shear.^2       .* medium.density;
lambda = medium.sound_speed_compression.^2 .* medium.density - 2*mu;


% rho0 is heterogeneous and staggered grids are used
%rho0_sgx = interpn(kgrid.x, kgrid.y, rho0, kgrid.x + kgrid.dx/2, kgrid.y,              '*linear');
%rho0_sgy = interpn(kgrid.x, kgrid.y, rho0, kgrid.x,              kgrid.y + kgrid.dy/2, '*linear');

% set values outside of the interpolation range to original values 
%rho0_sgx(isnan(rho0_sgx)) = rho0(isnan(rho0_sgx));
%rho0_sgy(isnan(rho0_sgy)) = rho0(isnan(rho0_sgy));
rho0_sgx = rho0;
rho0_sgy = rho0;



% invert rho0 so it doesn't have to be done each time step
rho0_sgx_inv = 1./rho0_sgx;
rho0_sgy_inv = 1./rho0_sgy;

% clear unused variables
clear rho0_sgx rho0_sgy

% mu is heterogeneous and staggered grids are used
% mu_sgxy  = 1./interpn(kgrid.x, kgrid.y, 1./mu, kgrid.x + kgrid.dx/2, kgrid.y + kgrid.dy/2, '*linear');
mu_sgxy = mu;

% set values outside of the interpolation range to original values 
mu_sgxy(isnan(mu_sgxy)) = mu(isnan(mu_sgxy));


% =========================================================================
% PREPARE DERIVATIVE AND PML OPERATORS
% =========================================================================

% get the regular PML operators based on the reference sound speed and PML settings
pml_x     = getPML(kgrid.Nx, kgrid.dx, kgrid.dt, c_ref, pml_x_size, pml_x_alpha, false, 1);
pml_x_sgx = getPML(kgrid.Nx, kgrid.dx, kgrid.dt, c_ref, pml_x_size, pml_x_alpha, true,  1);
pml_y     = getPML(kgrid.Ny, kgrid.dy, kgrid.dt, c_ref, pml_y_size, pml_y_alpha, false, 2);
pml_y_sgy = getPML(kgrid.Ny, kgrid.dy, kgrid.dt, c_ref, pml_y_size, pml_y_alpha, true,  2);

multi_axial_PML_ratio = 0.1;

% get the multi-axial PML operators
mpml_x     = getPML(kgrid.Nx, kgrid.dx, kgrid.dt, c_ref, pml_x_size, multi_axial_PML_ratio * pml_x_alpha, false, 1);
mpml_x_sgx = getPML(kgrid.Nx, kgrid.dx, kgrid.dt, c_ref, pml_x_size, multi_axial_PML_ratio * pml_x_alpha, true,  1);
mpml_y     = getPML(kgrid.Ny, kgrid.dy, kgrid.dt, c_ref, pml_y_size, multi_axial_PML_ratio * pml_y_alpha, false, 2);
mpml_y_sgy = getPML(kgrid.Ny, kgrid.dy, kgrid.dt, c_ref, pml_y_size, multi_axial_PML_ratio * pml_y_alpha, true,  2);

% define the k-space derivative operators, multiply by the staggered
% grid shift operators, and then re-order using ifftshift (the optio)
ddx_k_shift_pos = ifftshift( 1i * kgrid.kx_vec .* exp( 1i * kgrid.kx_vec * kgrid.dx / 2) );
ddx_k_shift_neg = ifftshift( 1i * kgrid.kx_vec .* exp(-1i * kgrid.kx_vec * kgrid.dx / 2) );
ddy_k_shift_pos = ifftshift( 1i * kgrid.ky_vec .* exp( 1i * kgrid.ky_vec * kgrid.dy / 2) );
ddy_k_shift_neg = ifftshift( 1i * kgrid.ky_vec .* exp(-1i * kgrid.ky_vec * kgrid.dy / 2) );

% force the derivative and shift operators to be in the correct direction
% for use with BSXFUN
ddy_k_shift_pos = ddy_k_shift_pos.'; 
ddy_k_shift_neg = ddy_k_shift_neg.';

% =========================================================================
% DATA CASTING
% =========================================================================

% preallocate the loop variables using the castZeros anonymous function
% (this creates a matrix of zeros in the data type specified by data_cast)
ux_split_x   = zeros([kgrid.Nx, kgrid.Ny]);
ux_split_y   = zeros([kgrid.Nx, kgrid.Ny]);
ux_sgx       = zeros([kgrid.Nx, kgrid.Ny]);  % **
uy_split_x   = zeros([kgrid.Nx, kgrid.Ny]);
uy_split_y   = zeros([kgrid.Nx, kgrid.Ny]);
uy_sgy       = zeros([kgrid.Nx, kgrid.Ny]);  % **

sxx_split_x  = zeros([kgrid.Nx, kgrid.Ny]);
sxx_split_y  = zeros([kgrid.Nx, kgrid.Ny]);
syy_split_x  = zeros([kgrid.Nx, kgrid.Ny]);
syy_split_y  = zeros([kgrid.Nx, kgrid.Ny]);
sxy_split_x  = zeros([kgrid.Nx, kgrid.Ny]);
sxy_split_y  = zeros([kgrid.Nx, kgrid.Ny]);

duxdx        = zeros([kgrid.Nx, kgrid.Ny]);  % **
duxdy        = zeros([kgrid.Nx, kgrid.Ny]);  % **
duydy        = zeros([kgrid.Nx, kgrid.Ny]);  % **
duydx        = zeros([kgrid.Nx, kgrid.Ny]);  % **

dsxxdx       = zeros([kgrid.Nx, kgrid.Ny]);  % **
dsxydy       = zeros([kgrid.Nx, kgrid.Ny]);  % **
dsxydx       = zeros([kgrid.Nx, kgrid.Ny]);  % **
dsyydy       = zeros([kgrid.Nx, kgrid.Ny]);  % **

p            = zeros([kgrid.Nx, kgrid.Ny]);  % **

index_start = 1;
index_step = 1;
index_end = kgrid.Nt; 

source.s_mask = ones(size(kgrid.k));
source.sxx(:, 1) = -reshape(source.p0, 1, []) / 2;
source.sxx(:, 2) = source.sxx(:, 1);
source.syy(:, 1:2) = repmat(source.sxx(:, 1), [1, 2]);
s_source_pos_index = find(source.s_mask ~= 0);

sensor.record_start_index = 1;
sensor_mask_index = find(sensor.mask ~= 0);
s_source_sig_index = ':';


% start time loop
for t_index = index_start:index_step:index_end

    % compute the gradients of the stress tensor (these variables do not
    % necessaily need to be stored, they could be computed as needed)
     tep_raw = syy_split_x + syy_split_y;
     tep_fft = fft(syy_split_x + syy_split_y, [], 2);
     tep_ifft = real(ifft(tep_fft, [], 2));
    
    dsxxdx = real( ifft( bsxfun(@times, ddx_k_shift_pos, fft(sxx_split_x + sxx_split_y, [], 1)), [], 1) );
    dsyydy = real( ifft( bsxfun(@times, ddy_k_shift_pos, fft(syy_split_x + syy_split_y, [], 2)), [], 2) );
    dsxydx = real( ifft( bsxfun(@times, ddx_k_shift_neg, fft(sxy_split_x + sxy_split_y, [], 1)), [], 1) );
    dsxydy = real( ifft( bsxfun(@times, ddy_k_shift_neg, fft(sxy_split_x + sxy_split_y, [], 2)), [], 2) );

    % calculate the split-field components of ux_sgx and uy_sgy at the next
    % time step using the components of the stress at the current time step
    ux_split_x = bsxfun(@times, mpml_y, ...
                    bsxfun(@times, pml_x_sgx, ...
                        bsxfun(@times, mpml_y, ...
                            bsxfun(@times, pml_x_sgx, ux_split_x) ...
                        ) + kgrid.dt .* rho0_sgx_inv .* dsxxdx ... 
                     )...
                  );
             
    ux_split_y = bsxfun(@times, mpml_x_sgx, bsxfun(@times, pml_y, ...
                 bsxfun(@times, mpml_x_sgx, bsxfun(@times, pml_y, ux_split_y)) ...
                 + kgrid.dt .* rho0_sgx_inv .* dsxydy));
             
    uy_split_x = bsxfun(@times, mpml_y_sgy, bsxfun(@times, pml_x, ...
                 bsxfun(@times, mpml_y_sgy, bsxfun(@times, pml_x, uy_split_x)) ...
                 + kgrid.dt .* rho0_sgy_inv .* dsxydx));
    uy_split_y = bsxfun(@times, mpml_x,     bsxfun(@times, pml_y_sgy, ...
                 bsxfun(@times, mpml_x,     bsxfun(@times, pml_y_sgy, uy_split_y)) ...
                 + kgrid.dt .* rho0_sgy_inv .* dsyydy)); 
    
    % combine split field components (these variables do not necessarily
    % need to be stored, they could be computed when needed)
    ux_sgx = ux_split_x + ux_split_y;
    uy_sgy = uy_split_x + uy_split_y;
        
    % calculate the velocity gradients (these variables do not necessarily
    % need to be stored, they could be computed when needed)
    duxdx = real( ifft( bsxfun(@times, ddx_k_shift_neg, fft(ux_sgx, [], 1)), [], 1));      
    duxdy = real( ifft( bsxfun(@times, ddy_k_shift_pos, fft(ux_sgx, [], 2)), [], 2));
    
    duydx = real( ifft( bsxfun(@times, ddx_k_shift_pos, fft(uy_sgy, [], 1)), [], 1));
    duydy = real( ifft( bsxfun(@times, ddy_k_shift_neg, fft(uy_sgy, [], 2)), [], 2));   
        
    % update the normal and shear components of the stress tensor using
    % a lossless elastic model with a split-field multi-axial pml
    sxx_split_x = bsxfun(@times, mpml_y, ...
                        bsxfun(@times, pml_x, ...
                            bsxfun(@times, mpml_y, bsxfun(@times, pml_x, sxx_split_x)) ...
                  + kgrid.dt .* (2 .* mu + lambda) .* duxdx ));
              
    sxx_split_y = bsxfun(@times, mpml_x, bsxfun(@times, pml_y, ...
                  bsxfun(@times, mpml_x, bsxfun(@times, pml_y, sxx_split_y)) ...
                  + kgrid.dt .* lambda .* duydy ));

    syy_split_x = bsxfun(@times, mpml_y, bsxfun(@times, pml_x, ...
                  bsxfun(@times, mpml_y, bsxfun(@times, pml_x, syy_split_x)) ...
                  + kgrid.dt .* lambda .* duxdx));
    syy_split_y = bsxfun(@times, mpml_x, bsxfun(@times, pml_y, ...
                  bsxfun(@times, mpml_x, bsxfun(@times, pml_y, syy_split_y)) ...
                  + kgrid.dt .* (2 .* mu + lambda) .* duydy ));

    sxy_split_x = bsxfun(@times, mpml_y_sgy, bsxfun(@times, pml_x_sgx, ...
                  bsxfun(@times, mpml_y_sgy, bsxfun(@times, pml_x_sgx, sxy_split_x)) ...
                  + kgrid.dt .* mu_sgxy .* duydx));
    sxy_split_y = bsxfun(@times, mpml_x_sgx, bsxfun(@times, pml_y_sgy, ...
                  bsxfun(@times, mpml_x_sgx, bsxfun(@times, pml_y_sgy, sxy_split_y)) ...
                  + kgrid.dt .* mu_sgxy .* duxdy));
        
        
    % add in the pre-scaled stress source terms
    if t_index == 1        
       % add the source values to the existing field values 
        sxx_split_x(s_source_pos_index) = sxx_split_x(s_source_pos_index) + source.sxx(s_source_sig_index, t_index);
        sxx_split_y(s_source_pos_index) = sxx_split_y(s_source_pos_index) + source.sxx(s_source_sig_index, t_index);
    
        % add the source values to the existing field values 
        syy_split_x(s_source_pos_index) = syy_split_x(s_source_pos_index) + source.syy(s_source_sig_index, t_index);
        syy_split_y(s_source_pos_index) = syy_split_y(s_source_pos_index) + source.syy(s_source_sig_index, t_index);                 
    end

    % compute pressure from normal components of the stress
    p = -(sxx_split_x + sxx_split_y + syy_split_x + syy_split_y) / 2;     
    
    % extract required sensor data from the pressure and particle velocity
    % fields if the number of time steps elapsed is greater than
    % sensor.record_start_index (defaults to 1) 
    if t_index >= sensor.record_start_index
    
        % update index for data storage
        file_index = t_index - sensor.record_start_index + 1;
        
        % run sub-function to extract the required data
        sensor_data.p(:, file_index) = p(sensor_mask_index);
    end
end

% =========================================================================
% CLEAN UP
% =========================================================================

sensor_data = sensor_data.p;

% reorder the simulation data
sensor_data_reordered = reorderSensorData(kgrid, sensor, sensor_data);

% =========================================================================
% VISUALISATION
% =========================================================================

% plot the re-ordered sensor data
figure;
imagesc(sensor_data_reordered, [-1, 1]);
colormap(getColorMap);
ylabel('Sensor Position');
xlabel('Time Step');
colorbar;;
colorbar;

function kx_vec = makeDim(Nx, dx)

    % define the discretisation of the spatial dimension such that
    % there is always a DC component
    if rem(Nx, 2) == 0
    	% grid dimension has an even number of points
        nx = ((-Nx/2:Nx/2-1)/Nx).';
    else
        % grid dimension has an odd number of points
        nx = ((-(Nx-1)/2:(Nx-1)/2)/Nx).';
    end

    % force middle value to be zero in case 1/Nx is a recurring
    % number and the series doesn't give exactly zero
    nx(floor(Nx/2) + 1) = 0;
            
    % define the wavenumber vector components
    kx_vec = (2*pi/dx) .* nx;       
end

function pml = getPML(Nx, dx, dt, c, pml_size, pml_alpha, staggered, dimension, axisymmetric) 

    % check for axisymmetric input
    if nargin < 9 || isempty(axisymmetric)
        axisymmetric = false;
    end

    % define x-axis
    x = 1:pml_size;

    % create absorption profile
    if staggered

        % calculate the varying components of the pml using a staggered grid
        pml_left  = pml_alpha * (c / dx) * ( ((x + 0.5) - pml_size - 1) ./ (0 - pml_size) ).^4; 
        pml_right = pml_alpha * (c / dx) * ( (x + 0.5) ./ pml_size ).^4;

    else

        % calculate the varying components of the pml using a regular grid
        pml_left  = pml_alpha * (c / dx) * ( (x - pml_size - 1) ./ (0 - pml_size) ).^4;
        pml_right = pml_alpha * (c / dx) * ( x ./ pml_size ).^4;

    end

    % exponentiation
    pml_left  = exp(-pml_left  * dt / 2);
    pml_right = exp(-pml_right * dt / 2);

    % add the components of the pml to the total function, not adding the axial
    % side of the radial PML if axisymmetric 
    pml = ones(1, Nx);
    if ~axisymmetric
        pml(1:pml_size) = pml_left;
    end
    pml(end - pml_size + 1:end) = pml_right;

    % reshape the pml vector to be in the desired direction
    switch dimension
        case 1
            pml = pml.';
        case 3
            pml = reshape(pml, [1, 1, Nx]);
    end

end

% function y = sinck(x)
%     zero_vals = (x == 0);
%     y = sin(x + pi * zero_vals) ./ (x + zero_vals) + zero_vals;
% end