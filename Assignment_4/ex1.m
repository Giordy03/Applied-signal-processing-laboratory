clear 
clc
close all

M = 8;
fc = 9e8;
c = 3e8;
theta_1 = 0;

lambda = c/fc;
d = lambda/2;
psi_1 = 2*pi*d/lambda*sin(theta_1);

m = 0:(M-1);
m = m';
a_theta1 = exp(1j*2*pi*d/lambda*sind(theta_1)*m); 
w = 1/sqrt(M)*a_theta1;

S = 3601;
theta = linspace(-90, 90, S);
V = zeros(M, S);

for i = 1:S
    a_theta = exp(1j*2*pi*d/lambda*m*sind(theta(i)));
    V(:,i) = a_theta;
end

w_hermitian = w'; % hermitian -> complex conjugate transpose
pattern = w_hermitian*V;

% plot
figure
plot(theta, abs(pattern))
grid on
title('Beam pattern obtained with w^H a(\theta), \theta_1 = 0^\circ')
xlabel('\theta [°]')
ylabel('|pattern(\theta)|')


% diric
psi_diric = 2*pi*d/lambda*(sind(theta)-sind(theta_1));
dir = sqrt(M)*diric(psi_diric, M); 

figure
plot(theta, abs(dir))
grid on
title('Dirichlet function idealized, M = 8, \theta_1 = 0^\circ')
xlabel('\theta [°]')
ylabel('|dir(\theta)|')


figure
polarplot(deg2rad(theta), abs(dir), 'r', 'LineWidth', 2);
ax = gca;
ax.ThetaZeroLocation = 'top';
ax.ThetaDir = 'clockwise';
ax.ThetaLim = [-90 90];
title(['Dirichlet pattern: M = ', num2str(M), ', \theta_1 = ', num2str(theta_1), '^o'])

% polar
figure
polarplot(deg2rad(theta), abs(pattern), 'r', 'LineWidth', 2);
ax = gca;
ax.ThetaZeroLocation = 'top';
ax.ThetaDir = 'clockwise';
ax.ThetaLim = [-90 90];
title(['Conventional beamforming: M=',num2str(M),', \theta_1 = ', num2str(theta_1), '^o'])

%% point 5
M_vector = [2, 15, 30, 50];

figure
for i = 1:length(M_vector)
    M = M_vector(i);
    m = 0:(M-1);
    m = m';
    a_theta1 = exp(1j*2*pi*d/lambda*sind(theta_1)*m); 
    w = 1/sqrt(M)*a_theta1;
    V = zeros(M, S);
    
    for j = 1:S
        a_theta = exp(1j*2*pi*d/lambda*m*sind(theta(j)));
        V(:,j) = a_theta;
    end
    
    w_hermitian = w'; 
    pattern = w_hermitian*V;
    
    % polar
    subplot(2,2,i)
    polarplot(deg2rad(theta), abs(pattern), 'r', 'LineWidth', 1.5);
    ax = gca;
    ax.ThetaZeroLocation = 'top';
    ax.ThetaDir = 'clockwise';
    ax.ThetaLim = [-90 90];
    title(['M = ', num2str(M_vector(i))])
    fprintf('Max for M = %d is %f.\n', M, max(abs(pattern)));
end
sgtitle('Conventional Beamforming: Polar Patterns for Different M at \theta_1 = 0^\circ', 'FontSize',11)

% change d over lambda
d_over_lambda = [0.25, 0.5, 1, 2];
M = 8;
m = 0:(M-1);
m = m';
m_pm = [-(M-1):(M-1)]';

figure
for i = 1:length(d_over_lambda)
    a_theta1 = exp(1j*2*pi*d_over_lambda(i)*sind(theta_1)*m); 
    w = 1/sqrt(M)*a_theta1;
    V = zeros(M, S);

    for j = 1:S
        a_theta = exp(1j*2*pi*d_over_lambda(i)*m*sind(theta(j)));
        V(:,j) = a_theta;
    end
    
    w_hermitian = w'; % hermitian -> complex conjugate transpose
    pattern = w_hermitian*V;
    
    % polar
    subplot(2,2,i)
    polarplot(deg2rad(theta), abs(pattern), 'r', 'LineWidth', 2);
    ax = gca;
    ax.ThetaZeroLocation = 'top';
    ax.ThetaDir = 'clockwise';
    ax.ThetaLim = [-90 90];
    title(['d/\lambda = ', num2str(d_over_lambda(i))])

    valid_ind = abs(sind(theta_1) + m_pm/(d_over_lambda(i))) <= 1;
    theta_max = asind(sind(theta_1) + m_pm/(d_over_lambda(i)));
    theta_max = theta_max(valid_ind);

    fprintf('Theta of the max for d/lambda = %d are: %s\n', d_over_lambda(i), mat2str(theta_max, 4));

end
sgtitle('Conventional Beamforming: Polar Patterns for Different d/lambda at \theta_1 = 0^\circ', 'FontSize',11)

theta_null_max = 90;
m_pm = m_pm(m_pm ~= 0);

lambda_over_d_max_direct = (sind(theta_null_max) - sind(theta_1)) ./ (m_pm/M);
valid = lambda_over_d_max_direct > 0;
d_over_lambda = 1 ./ lambda_over_d_max_direct(valid);
disp('Possible d/lambda values to avoid grating lobes:')
disp(d_over_lambda)


%% point 7
M = 8;
d_lambda = 0.5;
theta_1 = 30;

theta1_null = asind(sind(theta_1) + 1/(d_lambda)/M);
thetaminus1_null = asind(sind(theta_1) - 1/(d_lambda)/M);
FNBW = abs(theta1_null - thetaminus1_null);
disp(FNBW)

theta_1 = 60;
M_options = 1:20;
for M = 1:length(M_options)
    theta1_null = asind(sind(theta_1) + 1/(d_lambda)/M);
    if theta1_null < 90
        fprintf('Minimum M: %d, Theta_null_1: %f \n', M, theta1_null)
        break
    end
end

%% point 8
M = 8;
m = 0:(M-1)';
fc = 9e8;
c = 3e8;

theta_1 = 0;
theta1_possible = [20, 10, 5];
for i = 1:length(theta1_possible)
    theta = [theta_1, theta1_possible(i), -40, 60, -75, 80];
    
    lambda = c/fc;
    d = lambda/2;
    A = zeros(M, length(theta));
    
    for p = 1:length(theta)
        a_theta = exp(1j*2*pi*d/lambda*m*sind(theta(p)));
        A(:,p) = a_theta;
    end 
    
    % without for loop
    A = exp(1j*2*pi*d/lambda*m'*sind(theta));
  
    sigma = 1e-5;
    R_y = A*A'+sigma*eye(M);
    w_MDVR = (R_y \ A(:,1)) / (A(:,1)' * (R_y \ A(:,1)));
    
    
    S = 3601;
    theta_val = linspace(-90, 90, S);
    V = zeros(M, S);
    
    for k = 1:S
        a_theta = exp(1j*2*pi*d/lambda*m*sind(theta_val(k)));
        V(:,k) = a_theta;
    end
   
    % without for loop
    V = exp(1j*2*pi*d/lambda*m'*sind(theta_val));
    
    w_hermitian = w_MDVR';
    pattern = w_hermitian*V;
    
    % polar
    figure
    polarplot(deg2rad(theta_val), abs(pattern), 'r', 'LineWidth', 2);
    ax = gca;
    ax.ThetaZeroLocation = 'top';
    ax.ThetaDir = 'clockwise';
    ax.ThetaLim = [-90 90];
    title(['MVDR Beamforming: \theta_1 = ', num2str(theta_1), '^\circ, first interferer at \theta = ', num2str(theta1_possible(i)), '^\circ'])
    hold on
    for c= 2:length(theta)
        polarplot(deg2rad([theta(c), theta(c)]), [0, max(abs(pattern))], '--b')
    end
end

%% point 14 
interferers = 31;
theta = linspace(-90, 90, interferers);
theta1 = 0;
theta = [theta1, theta];
lambda = c/fc;
d = lambda/2; 
A = zeros(M, length(theta));

% Si può fare senza for 
A = exp(1j*2*pi*d/lambda*m'*sind(theta));

sigma = 1e-5;
R_y = A*A'+sigma*eye(M);
w_MDVR = (R_y \ A(:,1)) / norm(R_y \ A(:,1), 2);

S = 3601;
theta_val = linspace(-90, 90, S);
V = zeros(M, S);

for k = 1:S
    a_theta = exp(1j*2*pi*d/lambda*m*sind(theta_val(k)));
    V(:,k) = a_theta; 
end
   
% without for loop 
V = exp(1j*2*pi*d/lambda*m'*sind(theta_val));

w_hermitian = w_MDVR';
pattern = w_hermitian*V;

A = exp(1j*2*pi*d/lambda*m*sind(theta1));
w_conv = 1/sqrt(M)*A';
w_hermitian_conv = w_conv';
pattern_conv = w_hermitian_conv*V;

% polar
figure
polarplot(deg2rad(theta_val), abs(pattern), 'r', 'LineWidth', 2);
ax = gca;
ax.ThetaZeroLocation = 'top';
ax.ThetaDir = 'clockwise';
ax.ThetaLim = [-90 90];
hold on;
polarplot(deg2rad(theta_val), abs(pattern_conv), '--b', 'LineWidth', 2);
title(['MVDR vs Conventional Beamforming (M = ', num2str(M), ', DoA = ', num2str(theta1), '°)'], 'FontSize', 14);
legend('MVDR', 'Conventional', 'Location', 'southoutside');
