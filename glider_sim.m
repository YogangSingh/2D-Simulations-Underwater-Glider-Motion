%% Set Glider Parameters

clear;                      % Clear workspace
clc;                        % Clear screen              

%% Set glider trajectory and other inputs


Xi_d = [deg2rad(-25),deg2rad(25),deg2rad(-25),deg2rad(25)];     % Glider trajectory by specifying desired glide path angles
V_d = 0.3;                  % Glider speed in the vertical plane
rp3_d = 0.05;               % Fixed position of primary moving mass in body coordinate along e3 axis
ballast_rate = 0.0250;      % Ballast rate input (kg/s) (Magnitude only, sign is assigned later)

%% Constants

g = 9.816;                  % Acceleration due to gravity
I3 = eye(3);                % Identity matrix (3x3)
Z3 = zeros(3);              % Zero matric (3x3)

i = [1 0 0]';               % Unit vector alog earth frame x-axis
j = [0 1 0]';               % Unit vector alog earth frame y-axis
k = [0 0 1]';               % Unit vector alog earth frame z-axis

%% Mass Properties

mh = 40;                    % Hull mass
mbar = 9;                   % Primary moving inernal mass
mb = 1.0;                   % Variable ballast mass     
mw = 0;                     % Fixed point mass
mv = mh + mw + mb + mbar;   % Total vehicle mass
m = 50;                     % Displacement
m0 = mv - m;                % Buoyancy

mf1 = 5;                    % Added mass term
mf2 = 60;                   % Added mass term
mf3 = 70;                   % Added mass term     
Mf = diag([mf1, mf2, mf3]); % Added mass matrix

M = mh*I3 + Mf;             % Total mass
m1 = M(1,1);
m2 = M(2,2);
m3 = M(3,3);

j1 = 4;                     % Inertia term  
j2 = 12;                    % Inertia term
j3 = 11;                    % Inertia term
J = diag([j1, j2, j3]);     % Total intertia

%% Force and Moment Coefficients 

KL = 132.5;                 % Lift coefficient
KL0 = 0;                    % Lift coefficient

KD = 25;                    % Drag Coefficient
KD0 = 2.15;                 % Drag Coefficient

KM = -100;                  % Moment Coefficient
KM0 = 0;                    % Moment Coefficient

KOmega1_2 = -50;            % Rotational Damping Coefficient
KOmega2_2 = -50;            % Rotational Damping Coefficient

%% Admissible values of xi_d

lim1 = rad2deg(atan(2*(KD/KL)*((KL0/KL) + nthroot(((KL0/KL)^2) + (KD0/KD), 2))));
lim2 = rad2deg(atan(2*(KD/KL)*((KL0/KL) - nthroot(((KL0/KL)^2) + (KD0/KD), 2))));
fprintf('Admissible values of xi_d     = (-90, %f) U (%f, 90)', lim2, lim1);

%% Calculate glider trajectory

len = length(Xi_d);
t = [];
r = [];

for i=1:len
    fprintf('\n\nIteration %f\n',i);

    %% Desired value of glide angle

    xi_d = Xi_d(i);
    fprintf('Chosen value of xi_d          = %f degrees\n', rad2deg(xi_d));
 
    %% Set glide direction & ballast rate
    % Glide direction (U = Upwards, D = Downwards)

    if xi_d > 0
        glide_dir = 'U';
        disp('Glider direction: Upward');
        ballast_rate = -abs(ballast_rate);          % Negative ballast rate for upward glide
    elseif xi_d < 0
        glide_dir = 'D';
        disp('Glider direction: Downward');
        ballast_rate = abs(ballast_rate);           % Positive ballast rate for downward glide
    end
  
    %% Desired angle of attack

    alpha_d = (1/2)*(KL/KD)*(tan(xi_d))*(-1 + nthroot(1 - 4*(KD/(KL^2))*(cot(xi_d))*(KD0*cot(xi_d) + KL0),2));
    fprintf('Desired valued of alpha_d     = %f degrees\n', rad2deg(alpha_d));

    %% Desired velocity in body frame

    v1_d = V_d*cos(alpha_d);
    v3_d = V_d*sin(alpha_d);

    %% Desired ballast mass

    mb_d = (m - mbar - mh) + (1/g)*( (-sin(xi_d))*(KD0 + KD*(alpha_d^2)) + (cos(xi_d))*(KL0 + KL*(alpha_d)))*(V_d^2);
    fprintf('Desired valued of mb_d        = %f kg\n', mb_d);

    %% Desired position of longitudinal moving mass

    theta_d = alpha_d + xi_d;
    rp1_d = -rp3_d*tan(theta_d) + (1/(mbar*g*cos(theta_d)))*((mf3 - mf1)*v1_d*v3_d + (KM0 + KM*alpha_d)*(V_d^2));
    fprintf('Desired valued of rp1_d       = %f m\n', rp1_d);
    
    %% Save workspace

    save('glider_vars.mat');

    %% Set initial conditions for ODE solver

    % Initial params in first in the first glide
    if i == 1
        y0 = [v1_d                  % v1_d
              v3_d                  % v3_d
              0.4                   % theta
              0                     % Omega
              rp1_d                 % rp1
              rp3_d                 % rp3
              0                     % Pp1
              0                     % Pp3
              0                     % x
              0                     % z
              mb_d]';               % mb
    else
        l_idx = length(r);
        y0 = [r(l_idx,1)            % v1_d
              r(l_idx,2)            % v3_d
              r(l_idx,3)            % theta
              r(l_idx,4)            % Omega
              r(l_idx,5)            % rp1
              r(l_idx,6)            % rp3
              r(l_idx,7)            % Pp1
              r(l_idx,8)            % Pp3
              r(l_idx,9)            % x
              r(l_idx,10)           % z
              r(l_idx,11)]';        % mb
    end
    
    %% Solve ODE

    tspan = [(i-1)*375 i*375];
    [T, R] = ode45(@vertical_plane_model,tspan,y0);
    t = [t;T];
    r = [r;R];
end

%% Plot results 

h = figure;
set(h, 'Position', [100 0 600 650]);

% Depth
subplot(5,1,1);
plot(t, r(:,10))
set(gca,'YDir','reverse');
ylabel('Depth (m)','interpreter','latex');

% Pitch angle
subplot(5,1,2);
plot(t, r(:,3)*180/pi);
ylabel('$\theta$ (deg)','interpreter','latex');

% Omega_2
subplot(5,1,3);
plot(t,r(:,4)*180/pi);
ylabel('$\Omega_2$ (deg/s)','interpreter','latex');

% Velocity
V = (r(:,1).^2 + r(:,2).^2).^0.5;
subplot(5,1,4);
plot(t,V);
ylabel('Velocity (m/s)','interpreter','latex');

% Angle of attack
subplot(5,1,5);
alpha = atan(r(:,2)./r(:,1))*180/pi;
plot(t,alpha);
xlabel('Time (Seconds)');
ylabel('$\alpha$ (deg)','interpreter','latex');

% Title
set(gcf,'NextPlot','add');
axes;
h = title(sprintf('Model Validation Study on SLOCUM Glider\n Vertical Plane Simulation'));
set(gca,'Visible','off');
set(h,'Visible','on');

% Save to disk
hgexport(gcf, 'Glider_VerticalPlane_Results.jpg', hgexport('factorystyle'), 'Format', 'jpeg');
saveas(gcf, 'Glider_VerticalPlane_Results','fig');
save('Glider_VerticalPlane_Results.mat');