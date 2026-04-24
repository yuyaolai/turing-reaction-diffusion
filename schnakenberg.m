%% ============================================================
%  SCHNAKENBERG
%% ============================================================

clear; clc; close all;

%% ============================================================
%  PARAMETERS  (rescaled for visible pattern emergence)
%% ============================================================

a  = 0.1;
b  = 0.9;
Du = 0.05;      % FIXED: smaller Du so diffusion length matches domain
Dv = 2.0;       % FIXED: ratio Dv/Du = 40 maintained, but absolute values smaller

u_ss = a + b;
v_ss = b / (a + b)^2;

fprintf('Steady state: u* = %.3f, v* = %.3f\n', u_ss, v_ss);

%% ============================================================
%  SECTION 1: DISPERSION RELATION (re-plot with fixed params)
%% ============================================================

fu =  2*u_ss*v_ss - 1;
fv =  u_ss^2;
gu = -2*u_ss*v_ss;
gv = -u_ss^2;
detJ = fu*gv - fv*gu;

k2 = linspace(0, 60, 1000);
Bk = -(fu + gv) + (Du + Dv)*k2;
Ck = Du*Dv*k2.^2 - (Du*gv + Dv*fu)*k2 + detJ;
disc = Bk.^2 - 4*Ck;
sigma = 0.5*(-Bk + sqrt(abs(disc)));
sigma(disc < 0) = real(0.5*(-Bk(disc<0) + sqrt(complex(disc(disc<0)))));
sigma = real(sigma);

figure('Name','Dispersion Relation (Fixed)','Color','white','Position',[50 500 650 380]);
plot(sqrt(k2), sigma, 'b-', 'LineWidth', 2.5); hold on;
yline(0,'k--','LineWidth',1.5);
unstable = sigma > 0;
if any(unstable)
    k_vals = sqrt(k2(unstable));
    patch([k_vals(1) k_vals(end) k_vals(end) k_vals(1)], ...
          [min(sigma) min(sigma) max(sigma(unstable)) max(sigma(unstable))], ...
          [0.8 0.9 1.0], 'FaceAlpha',0.35,'EdgeColor','none');
    text(mean(k_vals), max(sigma(unstable))*0.6, 'Unstable modes', ...
         'HorizontalAlignment','center','FontSize',12,'Color',[0 0 0.7]);
end
xlabel('Wavenumber k','FontSize',13);
ylabel('Growth rate \sigma(k)','FontSize',13);
title('Dispersion Relation','FontSize',14,'FontWeight','bold');
legend('Growth rate \sigma(k)','\sigma = 0 threshold','Location','northeast','FontSize',11);
grid on; box on; set(gca,'FontSize',11);

% Find dominant wavenumber (peak of sigma)
[~, idx_peak] = max(sigma);
k_peak = sqrt(k2(idx_peak));
lambda_pred = 2*pi / k_peak;
fprintf('Dominant wavenumber k* = %.2f  →  predicted wavelength λ* = %.2f\n', k_peak, lambda_pred);

%% ============================================================
%  SECTION 2: 1D SIMULATION — FIXED
%% ============================================================

fprintf('\n=== 1D Simulation (fixed) ===\n');

L1D  = 20;          % FIXED: smaller domain (~3-4 wavelengths)
N1D  = 200;
dx   = L1D / N1D;
x    = linspace(0, L1D, N1D);

dt     = 0.001;     % FIXED: smaller dt for stability with larger Dv
Tend   = 80;        % FIXED: run long enough for patterns to saturate
nsteps = round(Tend / dt);
save_every = round(nsteps / 50);  % 50 snapshots total

rng(42);
noise_amp = 0.02;
u1 = u_ss + noise_amp * randn(1, N1D);
v1 = v_ss + noise_amp * randn(1, N1D);

snap_t = [];
snap_u = [];

Lap1 = @(f) (circshift(f,-1) - 2*f + circshift(f,1)) / dx^2;

for n = 1:nsteps
    Ru = a - u1 + u1.^2 .* v1;
    Rv = b - u1.^2 .* v1;
    u1 = u1 + dt*(Du*Lap1(u1) + Ru);
    v1 = v1 + dt*(Dv*Lap1(v1) + Rv);

    if mod(n, save_every) == 0
        snap_t(end+1) = n*dt;
        snap_u(end+1,:) = u1;
    end
    if mod(n, round(nsteps/5)) == 0
        fprintf('  t = %.1f / %.1f\n', n*dt, Tend);
    end
end

figure('Name','1D Turing Patterns (Fixed)','Color','white','Position',[50 50 900 420]);

subplot(1,2,1);
imagesc(x, snap_t, snap_u);
colormap(gca, parula); colorbar;
xlabel('Space x','FontSize',13); ylabel('Time t','FontSize',13);
title('Activator u(x,t) — Space-Time Evolution','FontSize',13,'FontWeight','bold');
set(gca,'YDir','normal','FontSize',11);

subplot(1,2,2);
cols = lines(4);
idx4 = round(linspace(1, length(snap_t), 4));
hold on;
leg_str = {};
for i = 1:4
    plot(x, snap_u(idx4(i),:), 'Color', cols(i,:), 'LineWidth', 2);
    leg_str{i} = sprintf('t = %.1f', snap_t(idx4(i)));
end
xlabel('Space x','FontSize',13); ylabel('Activator u(x)','FontSize',13);
title('Pattern Snapshots at 4 Time Points','FontSize',13,'FontWeight','bold');
legend(leg_str,'Location','best','FontSize',10);
grid on; box on; set(gca,'FontSize',11);

sgtitle('1D Turing Reaction-Diffusion', ...
        'FontSize',14,'FontWeight','bold');

%% ============================================================
%  SECTION 3: 2D SIMULATION — FIXED
%  Spots vs Stripes side by side
%% ============================================================

fprintf('\n=== 2D Simulation (fixed) ===\n');
fprintf('Estimated time: 2-4 minutes total...\n');

L2D  = 20;
N2D  = 100;         % resolution
dx2  = L2D / N2D;
dt2    = 0.001;
% Tend2  = 100;
Tend2  = 200;
nsteps2 = round(Tend2 / dt2);

Lap2D = @(F) (circshift(F,-1,1) + circshift(F,1,1) + ...
              circshift(F,-1,2) + circshift(F,1,2) - 4*F) / dx2^2;

% Two parameter sets giving spots vs stripes
psets(1) = struct('a',0.1,'b',0.9,'label','Spots  (a=0.1, b=0.9)');
psets(2) = struct('a', 0.1, 'b', 1.5, 'label', 'Denser Spots (a=0.1, b=1.5)');
psets(3) = struct('a', 0.1, 'b', 2.2, 'label', 'Labyrinthine (a=0.1, b=2.2)');
results2D = cell(1,3);
noise_amp = 0.02;

fprintf('=== Running Schnakenberg Simulation ===\n');

for p = 1:3
    ap = psets(p).a;
    bp = psets(p).b;
    uss2 = ap + bp;
    vss2 = bp / (ap+bp)^2;

    rng(13);
    u2 = uss2 + noise_amp * randn(N2D, N2D);
    v2 = vss2 + noise_amp * randn(N2D, N2D);

    % --- Time Stepping ---
    for n = 1:nsteps2

        % Kinetics
        Ru2 = ap - u2 + u2.^2.*v2;
        Rv2 = bp - u2.^2.*v2;
        
        % Update
        u2 = u2 + dt2*(Du*Lap2D(u2) + Ru2);
        v2 = v2 + dt2*(Dv*Lap2D(v2) + Rv2);
    end
    results2D{p} = u2;
    fprintf('  Set %d (%s) done.\n', p, psets(p).label);
end

% --- Visualization ---
figure('Name','2D Turing Patterns (Fixed)','Color','white','Position',[50 50 1000 300]);
for p = 1:3
    subplot(1,3,p);
    imagesc(results2D{p});
    colormap(gca, hot);
    colorbar;
    axis square off;
    title(psets(p).label,'FontSize',14,'FontWeight','bold');
end
% sgtitle('2D Turing Patterns — Spots vs. Stripes','FontSize',15,'FontWeight','bold');
sgtitle('2D Turing Patterns — Schnakenberg Model','FontSize',15,'FontWeight','bold');

%% ============================================================
%  SECTION 4: SENSITIVITY — Effect of Dv/Du ratio
%% ============================================================

figure('Name','Sensitivity: Diffusion Ratio','Color','white','Position',[50 50 750 430]);

ratios = [10, 20, 40, 60, 80];   % Dv/Du ratio
cols_s = cool(length(ratios));
k2s = linspace(0, 60, 800);
hold on;
leg_s = {};

for i = 1:length(ratios)
    Dvi = Du * ratios(i);
    Bki = -(fu+gv) + (Du+Dvi)*k2s;
    Cki = Du*Dvi*k2s.^2 - (Du*gv+Dvi*fu)*k2s + detJ;
    di  = Bki.^2 - 4*Cki;
    si  = real(0.5*(-Bki + sqrt(complex(di))));
    plot(sqrt(k2s), si, 'Color', cols_s(i,:), 'LineWidth', 2);
    leg_s{i} = sprintf('D_v/D_u = %d', ratios(i));
end
yline(0,'k--','LineWidth',1.5);
xlabel('Wavenumber k','FontSize',13);
ylabel('Growth rate \sigma(k)','FontSize',13);
title('Effect of Diffusion Ratio D_v/D_u on Pattern Scale','FontSize',13,'FontWeight','bold');
legend(leg_s,'Location','northeast','FontSize',11);
grid on; box on; set(gca,'FontSize',11);
annotation('textbox',[0.13 0.13 0.4 0.09],'String', ...
    'Higher D_v/D_u  \rightarrow  finer patterns (larger k)', ...
    'FontSize',11,'EdgeColor','none','Color',[0.1 0.1 0.6]);

fprintf('\n=== ALL DONE — 4 figures ready for slides ===\n');
