clc
clear all
% Define beam parameters
L = 1;          % beam length
H = 0.2;        % beam height
b = 0.02;       % beam thickness
E = 500000;        % Young's modulus
nu = 0.3;       % Poisson's ratio
P = 2000;       % applied load
% Define grid parameters
nx = 10;        % number of grid points in x-direction
ny = 10;        % number of grid points in y-direction
dx = L/nx;      % grid spacing in x-direction
dy = H/ny;      % grid spacing in y-direction
% Define force distribution
f = zeros(nx,ny);
for i = 1:nx
    for j = 1:ny
        if i*dx < L/2
            f(i,j) = P/(b*H);
        end
    end
end
% Set up initial guess for displacement function
u = zeros(nx,ny);
% Define tolerance for convergence
tol = 1e-6;
% Iteratively solve for displacement function using finite difference method
max_iter = 1000;
for iter = 1:max_iter
    u_old = u;
    for i = 2:nx-1
        for j = 2:ny-1
            u(i,j) = (dy^2*(u(i+1,j) + u(i-1,j)) + dx^2*(u(i,j+1) + u(i,j-1)) - dx^2*dy^2*f(i,j))/(2*(dx^2 + dy^2));
        end
    end
    % Check for convergence
    err = norm(u - u_old)/norm(u);
    if err < tol
        break
    end
end
% Plot results
x = linspace(0,L,nx);
y = linspace(0,H,ny);
[X,Y] = meshgrid(x,y);
figure;
surf(X,Y,u);
xlabel('x');
ylabel('y');
zlabel('u');
title('Displacement of Beam');