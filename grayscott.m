%% ============================================================
%  CONSISTENT REACTION-DIFFUSION SOLVER (GRAY-SCOTT)
%  Aligned with Schnakenberg numerical style
%% ============================================================
clear; clc; close all;

% --- Simulation Domain ---
L2D    = 100;           % Physical length
N2D    = 150;           % Grid resolution
dx2    = L2D / N2D;
dt2    = 0.2;           % Time step
Tend2  = 4000;          % Total time
nsteps2 = round(Tend2 / dt2);

% --- Diffusion Constants ---
Du = 0.16; 
Dv = 0.08;

% --- Laplacian Operator (Vectorized) ---
Lap2D = @(Z) (circshift(Z,-1,1) + circshift(Z,1,1) + ...
              circshift(Z,-1,2) + circshift(Z,1,2) - 4*Z) / dx2^2;

% --- Parameter Sets (F, k) ---
psets(1) = struct('F',0.024, 'k',0.060, 'label','Spots (F=0.024, k=0.060)');
psets(2) = struct('F',0.027, 'k',0.060, 'label','Denser Spots (F=0.027, k=0.060)');
psets(3) = struct('F',0.040, 'k',0.060, 'label','Labyrinthine (F=0.040, k=0.060)');
results2D = cell(1,3);
noise_amp = 0.02;

fprintf('=== Running Gray-Scott Simulation ===\n');

for p = 1:3
    Fp = psets(p).F;
    kp = psets(p).k;
    
    % Initial Conditions: Perturbed Steady State (u=1, v=0)
    rng(13); 
    u = ones(N2D, N2D);
    v = zeros(N2D, N2D);
    
    % Seed a small random area in the center to trigger instability
    % (Standard for Gray-Scott to ensure pattern birth)
    mid = floor(N2D/2);
    range = mid-10:mid+10;
    u(range, range) = 0.5 + noise_amp*randn(21,21);
    v(range, range) = 0.25 + noise_amp*randn(21,21);

    % --- Time Stepping ---
    for n = 1:nsteps2
        uvv = u .* (v.^2);
        
        % Kinetics
        Ru = -uvv + Fp*(1 - u);
        Rv =  uvv - (Fp + kp)*v;
        
        % Update
        u = u + dt2 * (Du * Lap2D(u) + Ru);
        v = v + dt2 * (Dv * Lap2D(v) + Rv);
        
        % Numerical Stability Clamp
        u = max(0, min(1, u));
        v = max(0, min(1, v));
    end
    results2D{p} = v; % Storing v for Gray-Scott visual contrast
    fprintf('  Set %d (%s) done.\n', p, psets(p).label);
end

%% --- Visualization ---
figure('Name','Unified Gray-Scott Results','Color','white','Position',[50 50 1000 300]);
for p = 1:3
    subplot(1,3,p);
    imagesc(results2D{p});
    colormap(gca, hot);
    colorbar;
    axis square off;
    title(psets(p).label, 'FontSize', 14,'FontWeight','bold');
end
sgtitle('2D Turing Patterns — Gray-Scott Model', 'FontSize', 15, 'FontWeight', 'bold');