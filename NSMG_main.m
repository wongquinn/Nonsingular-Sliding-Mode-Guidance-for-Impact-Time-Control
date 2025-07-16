%% Nonsingular Sliding Mode Guidance for Impact Time Control - MATLAB Code
% This script reproduces the numerical simulation results of the following paper:
%   Cho D, Kim H J, Tahk M J. "Nonsingular Sliding Mode Guidance for Impact Time Control,"
%   *Journal of Guidance, Control, and Dynamics*, Vol. 39, No. 1, 2016, pp. 61–68.
% Author: Dr_ikun
% The video explaining the code can be found at: 
% https://www.bilibili.com/video/BV11G4y1V7kt/?spm_id_from=333.1387.upload.video_card.click&vd_source=ef8d042607d313515e1983dfa42e807d

%% === Global Parameters ===
global V N epsilon1 epsilon2 k M p tf target_x target_y data;

V = 250;               % Constant speed [m/s]
N = 4;                 % Navigation gain
epsilon1 = 1e-3;       % constant design parameter
epsilon2 = 0.015;      % constant design parameter
k = 50;                % constant design parameter
M = 10;                % constant design parameter
p = 2;                 % constant design parameter
tf = 60;               % Desired impact time [s]
target_x = 10000;      % Target x-position [m]
target_y = 0;          % Target y-position [m]
data = [];             % Data storage for analysis

%% === Initial State ===
initial_x = 0;
initial_y = 0;
initial_theta = deg2rad(120);   % Convert initial heading to radians

%% === ODE Solver Setup ===
options = odeset('events', @(t, z) event2(t, z), ...
                 'AbsTol', 1e-12, 'RelTol', 1e-12);

[t, trajectory] = ode45(@(t, z) dynamics(t, z), ...
                        [0, 0.01:100], ...
                        [initial_x, initial_y, initial_theta], ...
                        options);

%% === Plotting ===
figure(1);
plot(trajectory(:,1), trajectory(:,2), '-', 'LineWidth', 5);
axis equal;
xlabel('x [m]'); ylabel('y [m]');
title('Missile Trajectory');
hold on;

figure(2);
plot(data(:,1), data(:,2)/9.81, '--', 'LineWidth', 5);
xlabel('Time [s]'); ylabel('Lateral Acceleration [g]');
title('Lateral Acceleration History');
hold on;

%% === Cost Function (control effort of lateral acceleration) ===
cost = 0.5 * trapz(linspace(0, tf, length(data)), data(:,2).^2);
disp(['control effort: ', num2str(cost)]);

%%%%%%%%%%%%% Function definitions
function dz = dynamics(t, z)
    % Dynamics function for missile guidance
    global V N k M p tf target_x target_y data;

    dz = zeros(3,1);
    x = z(1);
    y = z(2);
    gamma = z(3);  % Heading angle

    % Kinematics
    dz(1) = V * cos(gamma);  % x_dot
    dz(2) = V * sin(gamma);  % y_dot

    % Compute range and look angle
    r = sqrt((x - target_x)^2 + (y - target_y)^2);
    theta = gamma - atan2(target_y - y, target_x - x);

    % Rate of line-of-sight angle
    lambda_dot = -V * sin(theta) / r;

    % Estimated time-to-go (TGO)
    tgo = r / V * (1 + 0.5 * theta^2 / (2 * N - 1));

    % Sliding surface
    s = t + tgo - tf;

    % Equivalent control + linear control + switching control
    if theta == 0
        am_eq = 0;
        amMcon = 0;
    else
        am_eq = V^2 * (2*N - 1) / r * (cos(theta) - 1) / theta + ...
                0.5 * V^2 / r * cos(theta) * theta + ...
                V * lambda_dot;
        amMcon = -hfunction(theta) * k * s / theta;
    end

    am_sw = -M * (p * sgmf(theta) + 1) * sgmf(s);
    am = am_eq + amMcon + am_sw;

    % Limit lateral acceleration
    am = max(min(am, 100), -100);

    dz(3) = am / V;  % Heading rate

    % Store data: [time, acceleration, s, theta]
    data(end+1, :) = [t, am, s, theta];
end


function out = sgmf(in)
    % Smooth sign function for sliding mode
    alpha = 20;
    out = 2 * (1 / (1 + exp(-alpha * in)) - 0.5);
end


function h_fun = hfunction(theta)
    global epsilon1 epsilon2;

    if abs(theta) < epsilon1
        h_fun = sgmf(theta) * theta;
    elseif abs(theta) <= epsilon2
        h_fun = (abs(theta) * (1 - epsilon1) + epsilon1 * (epsilon2 - 1)) / (epsilon2 - epsilon1);
    else
        h_fun = 1;
    end
end


function [value, isterminal, direction] = event2(t, z)
    % Event function to stop integration near the target
    global target_x target_y;
    x_position = z(1);
    y_position = z(2);
    value = sqrt((x_position - target_x)^2 + (y_position - target_y)^2) - 0.5477;  % ≈ sqrt(0.3)
    isterminal = 1;   % Stop integration
    direction = 0;    % Detect all zero-crossings
end
