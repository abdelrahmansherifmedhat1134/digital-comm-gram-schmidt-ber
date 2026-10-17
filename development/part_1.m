clc; clear; close all;
%% ----------------------------------------------------------------
%  STEP 1: User inputs – amplitude matrix & time setup
Tb = 1;     
dt = 0.001; 
fprintf('Enter the signal matrix S.\n');
fprintf('  - Each ROW    = one signal\n');
fprintf('  - Each COLUMN = amplitude for one time slot of duration Tb = %g s\n', Tb);
fprintf('  - Signals can have different lengths (different number of columns)\n');
fprintf('  - Use 0 to pad slots where the signal is zero\n');
fprintf('Example:  [2 -1 -1 -1 ; -2 1 1 0 ; 0 2 2 -2 ; 1 -2 -2 2]\n\n');

S_amp = input('S = ');       

M         = size(S_amp, 1);    % number of signals
num_slots = size(S_amp, 2);    % number of time slots (matrix width)
T_total   = num_slots * Tb;    % total duration

t = 0 : dt : T_total - dt;    % common time axis
L = length(t);

% Build continuous signals from amplitude rows
signals      = zeros(M, L);
signal_names = cell(1, M);

for m = 1 : M
    for k = 1 : num_slots
        t_start = (k-1) * Tb;
        t_end   =  k    * Tb;
        signals(m, :) = signals(m, :) + S_amp(m, k) * (t >= t_start & t < t_end);
    end
    signal_names{m} = sprintf('s_%d(t)', m);
end
%% ----------------------------------------------------------------
%  STEP 2: Gram-Schmidt – produce orthonormal basis functions
basis_raw = zeros(M, L);
N         = 0;

for m = 1 : M

    v = signals(m, :);

    for n = 1 : N
        inner = trapz(v .* basis_raw(n,:)) * dt;
        v     = v - inner * basis_raw(n,:);
    end

    E_v = trapz(v .* v) * dt;

    if E_v > 1e-8
        N = N + 1;
        basis_raw(N, :) = v / sqrt(E_v);
    end

end

basis = basis_raw(1:N, :);

fprintf('\nM = %d signals   ->   N = %d basis functions\n', M, N);
%% ----------------------------------------------------------------
%  STEP 3: Constellation coordinates  c_mn

coords = zeros(M, N);

for m = 1 : M
    for n = 1 : N
        coords(m, n) = trapz(signals(m,:) .* basis(n,:)) * dt;
    end
end
%% ----------------------------------------------------------------
%  STEP 4: Symbol energies  (Parseval verification)

E_symbol = zeros(1, M);
E_direct = zeros(1, M);

for m = 1 : M
    E_symbol(m) = sum( coords(m,:).^2 );
    E_direct(m) = trapz(signals(m,:).^2) * dt;
end

fprintf('\n  Symbol energies (Parseval = direct integration):\n');
for m = 1 : M
    fprintf('    %s :  E_parseval = %7.4f   |   E_direct = %7.4f\n', ...
            signal_names{m}, E_symbol(m), E_direct(m));
end
%% ----------------------------------------------------------------
%  STEP 5: PLOT – Input Signals
figure('Name', 'Input Signals', 'NumberTitle','off', 'Color',[1 1 1]);

rows_plot = ceil(M / 2);
cols_plot = 2;

for m = 1 : M
    subplot(rows_plot, cols_plot, m);
    plot(t, signals(m,:), 'b', 'LineWidth', 2);
    xlabel('t  (s)');
    ylabel(signal_names{m});
    title(signal_names{m}, 'Interpreter','tex');
    grid on;
    xlim([t(1)  t(end)]);
    yline(0, 'k--', 'LineWidth', 0.8);
end

sgtitle(sprintf('Input Signals  (M = %d)', M), ...
        'FontSize', 13, 'FontWeight','bold');
%% ----------------------------------------------------------------
%  STEP 6: PLOT – Basis Functions

figure('Name', 'Basis Functions', 'NumberTitle','off', 'Color',[1 1 1]);

rows_b = ceil(N / 2);
cols_b = min(N, 2);

for n = 1 : N
    subplot(rows_b, cols_b, n);
    plot(t, basis(n,:), 'r', 'LineWidth', 2);
    xlabel('t  (s)');
    ylabel(sprintf('\\phi_%d(t)', n), 'Interpreter','tex');
    title(sprintf('Basis function  \\phi_%d(t)', n), 'Interpreter','tex');
    grid on;
    xlim([t(1)  t(end)]);
    yline(0, 'k--', 'LineWidth', 0.8);
end

sgtitle(sprintf('Gram-Schmidt Basis Functions  (N = %d)', N),'FontSize', 13, 'FontWeight','bold');
%% ----------------------------------------------------------------
%  STEP 7: PLOT – Constellation Diagram

figure('Name', 'Constellation Diagram', 'NumberTitle','off', 'Color',[1 1 1]);

colors_c = lines(M);

if N == 1
    hold on;
    for m = 1 : M
        scatter(coords(m,1), 0, 120, colors_c(m,:), 'filled');
        text(coords(m,1), 0.08, signal_names{m}, ...
             'HorizontalAlignment','center', 'FontSize', 10,'Interpreter','tex', 'Color', colors_c(m,:));
    end
    xlabel('\phi_1  coordinate', 'Interpreter','tex');
    title('Constellation Diagram  (1-D)', 'FontSize', 12);
    ylim([-0.5  0.5]);
    yline(0,'k-');
    grid on;

elseif N == 2
    hold on;
    for m = 1 : M
        scatter(coords(m,1), coords(m,2), 120, colors_c(m,:), 'filled');
        text(coords(m,1) + 0.05, coords(m,2) + 0.05, signal_names{m}, ...
             'FontSize', 10, 'Interpreter','tex', 'Color', colors_c(m,:));
    end
    xlabel('\phi_1  coordinate', 'Interpreter','tex');
    ylabel('\phi_2  coordinate', 'Interpreter','tex');
    title('Constellation Diagram  (2-D)', 'FontSize', 12);
    grid on;
    xline(0,'k--'); yline(0,'k--');

else
    hold on;
    for m = 1 : M
        scatter3(coords(m,1), coords(m,2), coords(m,3),120, colors_c(m,:), 'filled');
        text(coords(m,1), coords(m,2), coords(m,3), ...
             ['  ' signal_names{m}], 'FontSize', 10,'Interpreter','tex', 'Color', colors_c(m,:));
    end
    xlabel('\phi_1', 'Interpreter','tex');
    ylabel('\phi_2', 'Interpreter','tex');
    zlabel('\phi_3', 'Interpreter','tex');
    title(sprintf('Constellation Diagram  (3-D, first 3 of %d dims)', N), ...
          'FontSize', 12);
    grid on;
    view(35, 25);
end

e_str = '';
for m = 1 : M
    e_str = [e_str  sprintf('%s: E = %.3f    ', signal_names{m}, E_symbol(m))]; %#ok<AGROW>
end
text(0.02, -0.09, e_str, 'Units','normalized', 'FontSize', 8,'Interpreter','tex', 'Color',[0.3 0.3 0.3]);