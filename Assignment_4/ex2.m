%% EXERCISE 2
clc
clear 
close all

% Conventional Beamforming
Ny = [8, 32];
Nz = [8, 4];
theta_UE = [90, 105, 70, 100];
phi_UE = [0, 30, -45, 60];

UPA_BeamForming(Ny(1), Nz(1), theta_UE(1), phi_UE(1), 0);
UPA_BeamForming(Ny(1), Nz(1), theta_UE(1), phi_UE(1), 1);
 
UPA_BeamForming(Ny(1), Nz(1), theta_UE(2), phi_UE(2), 0);
UPA_BeamForming(Ny(1), Nz(1), theta_UE(2), phi_UE(2), 1);

UPA_BeamForming(Ny(1), Nz(1), theta_UE(3), phi_UE(3), 0);
UPA_BeamForming(Ny(1), Nz(1), theta_UE(3), phi_UE(3), 1);
 
UPA_BeamForming(Ny(2), Nz(2), theta_UE(4), phi_UE(4), 0);
UPA_BeamForming(Ny(2), Nz(2), theta_UE(4), phi_UE(4), 1);

% MVDR Beamforming
theta_intf1 = [86, 88];
phi_intf1 = [4, 2];
sigma_n2 = [1e-5, 1e3];

UPA_MVDR(theta_intf1(1), phi_intf1(1), sigma_n2(1), 0)
UPA_MVDR(theta_intf1(1), phi_intf1(1), sigma_n2(1), 1)

UPA_MVDR(theta_intf1(2), phi_intf1(2), sigma_n2(1), 0)
UPA_MVDR(theta_intf1(2), phi_intf1(2), sigma_n2(1), 1)

UPA_MVDR(theta_intf1(1), phi_intf1(1), sigma_n2(2), 0)
UPA_MVDR(theta_intf1(1), phi_intf1(1), sigma_n2(2), 1)

UPA_MVDR(theta_intf1(2), phi_intf1(2), sigma_n2(2), 0)
UPA_MVDR(theta_intf1(2), phi_intf1(2), sigma_n2(2), 1)

function [] = UPA_BeamForming(Ny, Nz, theta_UE, phi_UE, directive)
    % Nz: number of sensors along z
    % Ny: number of sensors along y
    
    c = 3e8;                 % speed of light
    fc = 6e9;                % frequency
    lambda = c/fc;           % wave length
    Ntot = Nz*Ny;
    d_z = lambda/2;          % sensor spacing along z [m]
    d_y = lambda/2;          % sensor spacing along y [m]
    %% ARRAY LIES IN THE Z–Y PLANE
    Tz=(2*[-Nz/2+1:Nz/2]-1)/3;
    Ty=2*[-Ny/2+1:Ny/2]-1;
    T = Tz'*ones(1,Ny)./Nz + 1j*ones(Nz,1)*Ty/Ny;
    T = T(:);
    
    %% Antenna element directivity function
    %%% 1) Define as anonymous function the directivity
    % function of the antenna element if directive ==1, otherwise define it equal to 1
    
    if(directive==1)
        d = @(theta, phi) 0.25.*(1-cos(2*theta)).*(1+cos(phi));
    else
        d = @(theta, phi) ones(size(theta));
    end
    
    %% Compute the beamforming filter for Conventional BF
    %% 2) Compute the BF vector as in slide 11 of Lec13
    
    
    m = [0: Ny-1]'; 
    a_y = @(theta, phi) exp(1j * 2*pi*d_y/lambda * sind(theta) * sind(phi) * m);
    n = [0:Nz-1]';
    a_z = @(theta) exp(1j * 2*pi*d_z/lambda * cosd(theta) * n);
    a_1 =  kron(a_y(theta_UE, phi_UE), a_z(theta_UE)); % dobbiamo farlo anche per l'altro pair
    w = 1/(sqrt(Nz*Ny))*a_1;
    
    %% Compute UPA pattern for Conventional BF
    angles_1 = 361;
    theta = linspace(0, pi, angles_1);
    angles_2 = 721;
    phi = linspace(-pi, pi, angles_2);
    
    %% 1) Use “meshgrid” to get matrices of all phi and theta
    [PHI, THETA] = meshgrid(phi, theta);
    
    %% 2) Generate the UPA pattern or AF (3D array) with isotropic antennas
    V = zeros([Ntot size(theta)]);
    theta_deg = rad2deg(theta);
    phi_deg = rad2deg(phi);
    for i = 1:length(theta_deg)
            az = a_z(theta_deg(i));  
        for j = 1:length(phi_deg)
            ay = a_y(theta_deg(i), phi_deg(j));  
            V(:, i, j)  = kron(ay, az);           
        end
    end
    
    
    %% Perform pattern multiplication
    %%% 3) Call you directivity function with input 
    % parameters the matrices PHI and THETA and 
    % perform elementwise multiplication with AF
    % Esegui pattern multiplication su tutta la griglia THETA, PHI
    
    
    % TotalPattern = d(THETA, PHI) .* (abs(V).^2).';  
    Pattern = zeros(size(THETA));
    
    for i = 1:length(theta)
        Pattern(i, :) = w' * squeeze(V(:, i, :));
    end
     
    TotalPattern = d(THETA, PHI) .* Pattern;
    Gain = abs(TotalPattern).^2;    
    
    %% Plot the pattern
    %%% 4) Convert the pattern into cartesian coordinates
    % (use matrices PHI and THETA), save into the 
    % variable r the square root of the array gain 
    % towards the user of interest (UE).
    
    % r = max(abs(TotalPattern));  % boh
    r = sqrt(Gain);
    x = r .* sin(THETA) .* cos(PHI);
    y = r .* sin(THETA) .* sin(PHI);
    z = r .* cos(THETA);
    r_C = sqrt(x.^2 + y.^2 + z.^2);  
    
    %%% 5) Plot the pattern with "mesh" with color scaling 
    % r_C and options 'FaceAlpha' and 'EdgeAlpha' 
    % at 0.5, set the axes equal and the "view" at 
    % (140, 15).
    
    figure
    mesh(x, y, z, r_C, 'FaceAlpha', 0.5, 'EdgeAlpha', 0.5);
    axis equal;
    view(140, 15);
    hold on;
    hidden on;
    plot3(zeros(Nz*Ny,1), imag(T), real(T), 'ks', 'MarkerFaceColor', 'r');
    r_broadside = max(r(:));
    x_broadside = r_broadside * sind(90) * cosd(0);
    y_broadside = r_broadside * sind(90) * sind(0);
    z_broadside = r_broadside * cosd(90);
    plot3([0, x_broadside], [0, y_broadside], [0, z_broadside], 'r--', 'LineWidth', 1.5);
    text(x_broadside*1.05, y_broadside*1.05, z_broadside*1.05, 'Broadside', 'FontSize', 12, 'Color', 'r', 'HorizontalAlignment', 'left');
    
    if theta_UE   ~= 90 || phi_UE ~= 0
        x_ue = r_broadside * sind(theta_UE) * cosd(phi_UE);
        y_ue = r_broadside * sind(theta_UE) * sind(phi_UE);
        z_ue = r_broadside * cosd(theta_UE);
        plot3([0, x_ue], [0, y_ue], [0, z_ue], 'r--', 'LineWidth', 1.5);
        label_ue = sprintf('UE (%.0f%s, %.0f%s)', ...
            theta_UE, char(176), phi_UE, char(176));
        text(x_ue*1.1, y_ue*1.1, z_ue*1.1, label_ue, 'FontSize', 12, 'Color', 'r');
    end
    
    xlabel('x'); ylabel('y'); zlabel('z');
    
    title(sprintf(['Radiation Pattern - %s\n', ...
        'Array: %s (%dx%d) | UE: el=%.0f%s, az=%.0f%s | directive=%d'], ...
        'Conventional BF', 'UPA', Nz, Ny, ...
        theta_UE, char(176), phi_UE, char(176), directive), ...
        'FontSize', 12);

    hold off;
end

function [] = UPA_MVDR(theta_intf1, phi_intf1, sigma_n2, directive)
    %% Additional parameters for MVDR
    n_intf = 5;              
    theta_intfdeg = [theta_intf1 85 80 100 105];
    theta_intf = deg2rad(theta_intfdeg);
    phi_intfdeg = [phi_intf1 20 5 -15 15];      
    phi_intf = deg2rad(phi_intfdeg);
    theta_1deg = 90;
    phi_1deg = 0;
    
    c = 3e8;                 % speed of light
    fc = 6e9;                % frequency
    lambda = c/fc;           % wave length
    d_z = lambda/2;          % sensor spacing along z [m]
    d_y = lambda/2;          % sensor spacing along y [m]
    Ny = 16;          
    Nz = 16;          
    Ntot = Ny * Nz;
    m = (0:(Ny-1)).'; 
    n = (0:(Nz-1)).';    
    a_y = @(theta, phi) exp(1j * (2*pi * d_y / lambda) .* sin(theta) .* sin(phi) .* m);
    a_z = @(theta) exp(1j * (2*pi * d_z / lambda) .* cos(theta) .* n);
    
    if(directive==1)
        d = @(theta, phi) 0.25.*(1-cos(2*theta)).*(1+cos(phi));
    else
        d = @(theta, phi) ones(size(theta));
    end
    theta_1 = deg2rad(theta_1deg);  
    phi_1 = deg2rad(phi_1deg);    
    
    %% Compute MVDR beamforming filter
    %%% 9) Compute the BF vector as in slide 13 of Lec13 with different normalization 
    % User’s steering vector
    az_u = a_z(theta_1);        
    ay_u = a_y(theta_1, phi_1);     
    a_1 = kron(az_u, ay_u);       
    
    A = zeros(Ntot, 1 + n_intf);
    A(:,1) = d(theta_1, phi_1) * a_1; 
    
    for k = 1:n_intf
        th_i = theta_intf(k);     
        ph_i = phi_intf(k);       
        az_i = a_z(th_i);          
        ay_i = a_y(th_i, ph_i);    
        a_intf = kron(az_i, ay_i);  
        % Interferer’s steering by element directivity:
        A(:, 1 + k) = d(th_i, ph_i) * a_intf;
    end
   
    % Covariance (P = 1)
    R_y = A * A' + sigma_n2 * eye(Ntot);
    
    % Compute the MVDR weight (normalized):
    w_MVDR = (R_y \ a_1) / norm(R_y \ a_1, 2);
    
    
    %% Compute UPA pattern for MVDR
    %%% 10) Perform elementwise multiplication between two 3−D arrays: 
    % the 3−D array containing all steering vector and the one with the
    % directivity function
    N_theta = 361;               
    N_phi = 721;     
    K = N_theta * N_phi;
    theta_grid = linspace(0, pi, N_theta);    
    phi_grid = linspace(-pi, pi, N_phi);    
    
    [PHI, THETA] = meshgrid(phi_grid, theta_grid);
    directivity = d(THETA, PHI);
    directivity_3D = repmat(directivity, [1, 1, Ntot]);
    directivity_3D = permute(directivity_3D, [3 1 2]);
    
    vectorized_THETA = THETA(:).';
    vectorized_PHI = PHI(:).';

    az = a_z(vectorized_THETA);              
    ay = a_y(vectorized_THETA, vectorized_PHI);
    
    az_kron = reshape(az, [Nz, 1, K]);   
    ay_kron = reshape(ay, [1, Ny, K]);   
    
    V_temp = az_kron .* ay_kron;        
    V = reshape(V_temp, [Ntot, N_theta, N_phi]);

    pattern_dir = V .* directivity_3D;

    %% 11) Multiply the BF vector with each steering vector (only use one for loop)
    TotalPattern_MVDR = zeros(N_theta, N_phi);
    for i = 1:N_theta
        V_slice = squeeze(pattern_dir(:, i, :));  
        TotalPattern_MVDR(i, :) = w_MVDR' * V_slice;          
    end
    Gain_MDVR = abs(TotalPattern_MVDR).^2; 

    %% Plot the pattern for MVDR
    %%% 12) Convert the pattern into cartesian coordinates (use matrices 
    % PHI and THETA), save into the variable r the square root of the array gain 
    % towards the user of interest (UE).
    r = sqrt(Gain_MDVR);              
    X = r .* sin(THETA) .* cos(PHI);  
    Y = r .* sin(THETA) .* sin(PHI); 
    Z = r .* cos(THETA);             
    r_C = sqrt(X.^2+Y.^2+Z.^2);                     
    
    %%% 13) Plot the pattern with "mesh" with color scaling r_C and options 
    % 'FaceAlpha' and 'EdgeAlpha' at 0.5, set the axes equal and the "view" at (96, 3).
    figure;
    mesh(X, Y, Z, r_C, 'FaceAlpha', 0.5, 'EdgeAlpha', 0.5);
    axis equal;
    view(96, 3);
    hold on;
    hidden on;
    
    % Recompute T for 16×16 UPA (to plot element positions):
    Tz = (2*[-Nz/2+1 : Nz/2] - 1)/3;   
    Ty = (2*[-Ny/2+1 : Ny/2] - 1);   
    T = (Tz' * ones(1,Ny)) / Nz + 1j*(ones(Nz,1)*(Ty/ Ny));  
    T = T(:); 
    
    plot3( zeros(Ntot,1), imag(T), real(T), 'ks', 'MarkerFaceColor','r');  % draws the array
    
    
    %% 14) Draw a red dashed line at UE DoA (90°, 0°) of length = r_user
    %%% 14) Draw a red dashed line at UE DoA (90, 0) with length equal to r with 
    % "plot3", add a red dot with 'MarkerSize' 8 at the end of the red line and 
    % add a text with writing 'User' and UE azimuth and el. angle (use 
    % char(176) for the degree symbol).
   
    r_user = max(r(:));
    x_user = r_user * sin(pi/2) * cos(0);
    y_user = r_user * sin(pi/2) * sin(0);
    z_user = r_user * cos(pi/2);
    
    plot3([0, x_user], [0, y_user], [0, z_user], 'r--', 'LineWidth', 1.5);
    plot3(x_user, y_user, z_user, 'r.', 'MarkerSize', 12);
    txt1 = sprintf('User (%.0f%s, %.0f%s)', theta_1deg, char(176), phi_1deg, char(176));
    text(x_user*1.1, y_user*1.1, z_user*1.1, txt1, 'FontSize', 12, 'Color', 'r');
    
    %% 15) Draw black dashed lines for each interferer at length r_user
    %%% 15) For each interferer intf_i: add a black dashed line with length r, 
    % add a text with el. angle and azimuth of intf_i at the end of each line
    for k = 1:n_intf
        th_i = theta_intf(k);  
        ph_i = phi_intf(k); 
    
        x_i = r_user * sin(th_i) * cos(ph_i);
        y_i = r_user * sin(th_i) * sin(ph_i);
        z_i = r_user * cos(th_i);
    
        plot3([0, x_i], [0, y_i], [0, z_i], 'k--', 'LineWidth', 1.2);
        txt_i = sprintf('intf_%d (%.0f%s, %.0f%s)', ...                     
            k, theta_intfdeg(k), char(176), phi_intfdeg(k), char(176));
        text(x_i*1.05, y_i*1.05, z_i*1.05, txt_i, 'FontSize', 12, 'Color', 'k');
        hold on;
    end
    
    %% 16) Labels & Title
    %%% 16) Add labels and the title with BF technique, array type and size, 
    % and elevation angle and azimuth of UE, set 'FontSize' to 12.
    
    xlabel('x'); ylabel('y'); zlabel('z');
    
    title(sprintf(['Radiation Pattern - %s\n', ...
        'Array: %s (%dx%d) | UE: el=%d%s, az=%d%s\n', ...
        'Intf1: el=%d%s, az=%d%s | \\sigma_{n}^{2}=%.2e | directive=%d' ], ...
        'MVDR BF', 'UPA', Nz, Ny, theta_1deg, char(176), phi_1deg, ...
        char(176), theta_intf1, char(176), phi_intf1, char(176), sigma_n2, ...
        directive ), 'FontSize', 12 );
    hold off;
end
