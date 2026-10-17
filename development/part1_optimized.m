% ECE352_Project_Phase1.m
% ---------------------------------------------------------
% Main Script Execution
% ---------------------------------------------------------
clc; clear; close all;

% 1. Get inputs from the user via the command window
[S_amp, Tb, dt] = get_user_inputs();

% 2. Process Data
[signals, t, signal_names] = generate_signals(S_amp, Tb, dt);
[basis, N]                 = apply_gram_schmidt(signals, dt);
[coords, E_symbol, E_dir]  = calc_constellation(signals, basis, dt);

% 3. Output Verification to Console
fprintf('\nGenerated %d basis functions from %d signals.\n', N, size(S_amp,1));

% 4. Plot the results
plot_input_signals(t, signals, signal_names);
plot_basis_functions(t, basis);
plot_constellation(coords, signal_names, E_symbol, N);

% ---------------------------------------------------------
% Local Functions Definitions (Must be at the bottom)
% ---------------------------------------------------------

function [S_amp, Tb, dt] = get_user_inputs()
    % Defines time constants and prompts the user for the signal matrix
    Tb = 1;     
    dt = 0.001; 
    
    fprintf('--- ECE 352: Phase 1 Input ---\n');
    fprintf('Enter the signal matrix S.\n');
    fprintf('  - Each ROW    = one signal\n');
    fprintf('  - Each COLUMN = amplitude for one time slot of duration Tb = %g s\n', Tb);
    fprintf('  - Use 0 to pad slots where the signal is zero\n');
    fprintf('Example for Set 1:  [2 2 -1 -1 ; 1 1 -1 -1 ; 2 0 0 -2 ; 1 -2 -2 2]\n\n');
    
    S_amp = input('S = ');       
    
    % Basic error check
    if isempty(S_amp)
        error('Input matrix cannot be empty. Please run the script again.');
    end
end

function [signals, t, signal_names] = generate_signals(S_amp, Tb, dt)
    % Generates continuous time arrays from amplitude slots
    M = size(S_amp, 1);
    num_slots = size(S_amp, 2);
    T_total = num_slots * Tb;
    t = 0 : dt : T_total - dt;
    L = length(t);

    signals = zeros(M, L);
    signal_names = cell(1, M);
    
    for m = 1 : M
        for k = 1 : num_slots
            t_start = (k-1) * Tb;
            t_end   =  k    * Tb;
            % Build the piecewise signal using logical indexing
            signals(m, :) = signals(m, :) + S_amp(m, k) * (t >= t_start & t < t_end);
        end
        signal_names{m} = sprintf('s_%d(t)', m);
    end
end

function [basis, N] = apply_gram_schmidt(signals, dt)
    % Takes continuous signals and outputs N orthonormal basis functions
    [M, L] = size(signals);
    basis_raw = zeros(M, L);
    N = 0;
    
    for m = 1 : M
        v = signals(m, :);
        for n = 1 : N
            % Calculate projection using trapezoidal integration
            inner = trapz(v .* basis_raw(n,:)) * dt;
            v = v - inner * basis_raw(n,:);
        end
        E_v = trapz(v .* v) * dt;
        % Check if remaining energy is significant
        if E_v > 1e-8
            N = N + 1;
            basis_raw(N, :) = v / sqrt(E_v);
        end
    end
    basis = basis_raw(1:N, :); % Trim unused rows
end

function [coords, E_symbol, E_direct] = calc_constellation(signals, basis, dt)
    % Projects signals onto basis functions to get coordinates and energies
    [M, ~] = size(signals);
    N = size(basis, 1);
    coords = zeros(M, N);
    E_symbol = zeros(1, M);
    E_direct = zeros(1, M);

    for m = 1 : M
        for n = 1 : N
            coords(m, n) = trapz(signals(m,:) .* basis(n,:)) * dt;
        end
        % Energy calculated via Parseval's theorem (Constellation distance)
        E_symbol(m) = sum(coords(m,:).^2);
        % Energy calculated directly from the time-domain waveform
        E_direct(m) = trapz(signals(m,:).^2) * dt;
    end
end

function plot_input_signals(t, signals, signal_names)
    M = size(signals, 1);
    figure('Name', 'Input Signals', 'NumberTitle','off', 'Color',[1 1 1]);
    rows_plot = ceil(M / 2);
    cols_plot = 2;
    for m = 1 : M
        subplot(rows_plot, cols_plot, m);
        plot(t, signals(m,:), 'b', 'LineWidth', 2);
        xlabel('t  (s)'); ylabel(signal_names{m});
        title(signal_names{m}, 'Interpreter','tex');
        grid on; xlim([t(1)  t(end)]); yline(0, 'k--', 'LineWidth', 0.8);
    end
    sgtitle(sprintf('Input Signals  (M = %d)', M), 'FontSize', 13, 'FontWeight','bold');
end

function plot_basis_functions(t, basis)
    N = size(basis, 1);
    figure('Name', 'Basis Functions', 'NumberTitle','off', 'Color',[1 1 1]);
    rows_b = ceil(N / 2);
    cols_b = min(N, 2);
    for n = 1 : N
        subplot(rows_b, cols_b, n);
        plot(t, basis(n,:), 'r', 'LineWidth', 2);
        xlabel('t  (s)'); ylabel(sprintf('\\phi_%d(t)', n), 'Interpreter','tex');
        title(sprintf('Basis function  \\phi_%d(t)', n), 'Interpreter','tex');
        grid on; xlim([t(1)  t(end)]); yline(0, 'k--', 'LineWidth', 0.8);
    end
    sgtitle(sprintf('Gram-Schmidt Basis Functions  (N = %d)', N),'FontSize', 13, 'FontWeight','bold');
end

function plot_constellation(coords, signal_names, E_symbol, N)
    M = size(coords, 1);
    figure('Name', 'Constellation Diagram', 'NumberTitle','off', 'Color',[1 1 1]);
    colors_c = lines(M);
    
    hold on;
    if N == 1
        for m = 1 : M
            scatter(coords(m,1), 0, 120, colors_c(m,:), 'filled');
            text(coords(m,1), 0.08, signal_names{m}, 'HorizontalAlignment','center', 'Color', colors_c(m,:));
        end
        xlabel('\phi_1  coordinate'); ylim([-0.5  0.5]);
    elseif N == 2
        for m = 1 : M
            scatter(coords(m,1), coords(m,2), 120, colors_c(m,:), 'filled');
            text(coords(m,1) + 0.05, coords(m,2) + 0.05, signal_names{m}, 'Color', colors_c(m,:));
        end
        xlabel('\phi_1  coordinate'); ylabel('\phi_2  coordinate');
        xline(0,'k--'); yline(0,'k--');
    else
        for m = 1 : M
            scatter3(coords(m,1), coords(m,2), coords(m,3), 120, colors_c(m,:), 'filled');
            text(coords(m,1), coords(m,2), coords(m,3), ['  ' signal_names{m}], 'Color', colors_c(m,:));
        end
        xlabel('\phi_1'); ylabel('\phi_2'); zlabel('\phi_3'); view(35, 25);
    end
    
    title(sprintf('Constellation Diagram (%d-D)', min(N, 3)), 'FontSize', 12); grid on;
    
    % Print Energies at the bottom
    e_str = '';
    for m = 1 : M
        e_str = [e_str sprintf('%s: E = %.3f    ', signal_names{m}, E_symbol(m))]; 
    end
    text(0.02, -0.09, e_str, 'Units','normalized', 'FontSize', 8, 'Color',[0.3 0.3 0.3]);
end