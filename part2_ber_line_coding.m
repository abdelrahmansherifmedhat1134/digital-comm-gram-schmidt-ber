clc; clear; close all;

%%MINI GUI

prompt = {'Number of generated bits:', 'Eb/No Min (dB):', 'Eb/No Max (dB):', 'Eb/No Step (dB):'};
dlgtitle = 'Simulation Parameters';
dims = [1 55];
definput = {'100000', '-10', '10', '2'}; % Default values
answer = inputdlg(prompt, dlgtitle, dims, definput);

% Exit if user closes the window
if isempty(answer)
    disp('Simulation canceled by user.');
    return; 
end

% Assign inputs to variables
num_bits = str2double(answer{1});
ebno_min = str2double(answer{2});
ebno_max = str2double(answer{3});
ebno_step = str2double(answer{4});

% Generate the Eb/No array from GUI inputs
EbNo_dB = ebno_min:ebno_step:ebno_max;
EbNo_lin = 10.^(EbNo_dB / 10);

% Prompt for Line Coding Selection
list = {'Unipolar NRZ', 'Manchester', 'Both'};
[indx, tf] = listdlg('PromptString', 'Select Line Coding:', ...
                     'SelectionMode', 'single', ...
                     'ListSize', [200, 80], ...
                     'ListString', list);

% Exit if user clicks cancel
if tf == 0
    disp('Simulation canceled by user.');
    return;
end

% Set boolean flags based on selection
do_nrz = (indx == 1 || indx == 3);
do_man = (indx == 2 || indx == 3);

%==========================================
% Setup & Generation of Base Waveforms

ts = 0.001;
T_sym = 1;
n_per_sym = T_sym / ts;
t = 0:ts:T_sym - ts;
A = 1;

if do_nrz
    s0_unrz = zeros(1, n_per_sym);
    s1_unrz = A * ones(1, n_per_sym);
    matrix_unrz = [s1_unrz; s0_unrz];
    [basis_unrz, coeff_unrz] = gram_schmidt(matrix_unrz, ts);
    const_diag(coeff_unrz);
    title('Constellation Diagram: Unipolar NRZ');
end

if do_man
    half_n = n_per_sym / 2;
    s0_man = [-A*ones(1, half_n), A*ones(1, half_n)];
    s1_man = [A*ones(1, half_n), -A*ones(1, half_n)];
    matrix_man = [s1_man; s0_man];
    [basis_man, coeff_man] = gram_schmidt(matrix_man, ts);
    const_diag(coeff_man);
    title('Constellation Diagram: Manchester');
end

%==========================================
% Generate Random Data & Encode

bits = randi([0, 1], 1, num_bits);

if do_nrz
    tx_unrz = zeros(num_bits, size(coeff_unrz, 2));
    tx_unrz(bits == 0, :) = repmat(coeff_unrz(2,:), sum(bits == 0), 1);
    tx_unrz(bits == 1, :) = repmat(coeff_unrz(1,:), sum(bits == 1), 1);
    Eb_unrz = mean(sum(coeff_unrz.^2,2));
    ber_unrz_sim = zeros(1, length(EbNo_dB));
    fig_rx_unrz = figure('Name', 'Received Constellation Diagram- Unipolar NRZ');
end

if do_man
    tx_man = zeros(num_bits, size(coeff_man, 2));
    tx_man(bits == 0, :) = repmat(coeff_man(2,:), sum(bits == 0), 1);
    tx_man(bits == 1, :) = repmat(coeff_man(1,:), sum(bits == 1), 1);
    Eb_man = mean(sum(coeff_man.^2,2));
    ber_man_sim = zeros(1, length(EbNo_dB));
    fig_rx_man = figure('Name', 'Received Constellation Diagram - Manchester');
end
% =========================================================================
%               AWGN Channel & Decoding

% Calculate dynamic subplot grid dimensions to prevent index errors
num_plots = length(EbNo_dB);
n_cols = ceil(sqrt(num_plots)); 
n_rows = ceil(num_plots / n_cols);

for k = 1:num_plots
    
    if do_nrz
        % Unipolar NRZ Channel
        N0_unrz = Eb_unrz / EbNo_lin(k);
        noise_unrz = sqrt(N0_unrz/2) * randn(size(tx_unrz));
        rx_unrz = tx_unrz + noise_unrz;
        
        % Plot Rx Constellation with DYNAMIC grid
        figure(fig_rx_unrz); 
        subplot(n_rows, n_cols, k); 
        scatter(rx_unrz, zeros(size(rx_unrz)), 1, 'filled', 'MarkerFaceAlpha', 0.5); hold on;
        % Plot the ideal points
        plot(coeff_unrz(:,1), [0 0], 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 5);
        xline(0, 'k--', 'LineWidth', 1);
        title(sprintf('Eb/No = %d dB', EbNo_dB(k)));
        
        % Decode Data (Minimum Euclidean Distance)
        dist0_unrz = sum((rx_unrz - coeff_unrz(2,:)).^2, 2);
        dist1_unrz = sum((rx_unrz - coeff_unrz(1,:)).^2, 2);
        rx_bits_unrz = zeros(1, num_bits);
        rx_bits_unrz(dist1_unrz < dist0_unrz) = 1;
        
        % Calculate Simulated BER
        ber_unrz_sim(k) = sum(rx_bits_unrz ~= bits) / num_bits;
    end
    
    if do_man
        % Manchester Channel
        N0_man = Eb_man / EbNo_lin(k);
        noise_man = sqrt(N0_man/2) * randn(size(tx_man));
        rx_man = tx_man + noise_man;
        
        % Plot Rx Constellation with DYNAMIC grid
        figure(fig_rx_man); 
        subplot(n_rows, n_cols, k); 
        scatter(rx_man, zeros(size(rx_man)), 1, 'filled', 'MarkerFaceAlpha', 0.5);  hold on;
        % Plot the ideal points
        plot(coeff_man(:,1), [0 0], 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 5);
        yticks([]);
        xline(0, 'k--', 'LineWidth', 1);
        title(sprintf('Eb/No = %d dB', EbNo_dB(k)));
        
        % Decode Data (Minimum Euclidean Distance)
        dist0_man = sum((rx_man - coeff_man(2,:)).^2, 2);
        dist1_man = sum((rx_man - coeff_man(1,:)).^2, 2);
        rx_bits_man = zeros(1, num_bits);
        rx_bits_man(dist1_man < dist0_man) = 1;
        
        % Calculate Simulated BER
        ber_man_sim(k) = sum(rx_bits_man ~= bits) / num_bits;
    end
end

% =========================================================================
%               Theoretical BER Calculation & Plotting

figure('Name', 'BER Performance Curves');

if do_nrz
    ber_unrz_th = qfunc(sqrt(EbNo_lin));
    semilogy(EbNo_dB, ber_unrz_sim, 'bo', 'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', 'Unipolar NRZ (Simulated)');
    hold on;
    semilogy(EbNo_dB, ber_unrz_th, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Unipolar NRZ (Theoretical)');
end

if do_man
    ber_man_th = qfunc(sqrt(2 * EbNo_lin));   
    semilogy(EbNo_dB, ber_man_sim, 'rx', 'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', 'Manchester (Simulated)');
    hold on;
    semilogy(EbNo_dB, ber_man_th, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Manchester (Theoretical)');
end

% Force the Y-axis to be logarithmic
set(gca, 'YScale', 'log');

grid on;
xlabel('E_b/N_o (dB)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Bit Error Rate (BER)', 'FontSize', 12, 'FontWeight', 'bold');
title('Theoretical vs Practical BER', 'FontSize', 14);
legend('Location', 'southwest');

% =========================================================================
%                 ENERGY USING TIME DOMAIN FUNCTION
function E = energy(signal,ts)
    E = sum(abs(signal).^2)*ts;
end
% =========================================================================
%                  GRAM SCHMIDT FUNCTION
function [basis_functions,basis_coeff] = gram_schmidt(signals_matrix,ts)
    M = length(signals_matrix(:,1));
    L = length(signals_matrix);
    basis_functions = zeros(M,L);
    basis_coeff = zeros(M);
    E_s1 = energy(signals_matrix(1,:),ts);
    basis_functions(1,:) = signals_matrix(1,:)/sqrt(E_s1);
    basis_coeff(1) = sqrt(E_s1); 
    N = 1;
    for i = 2:M
        for j = 1:N
            basis_coeff(i,j) = sum(basis_functions(j,:).*signals_matrix(i,:))*ts;
        end
        g_i = signals_matrix(i,:) - basis_coeff(i,:)*basis_functions;
        if(max(abs(g_i)) < 1e-10)
            continue;
        else
            E_gi = energy(g_i,ts);
            basis_functions(i,:) = g_i/sqrt(E_gi); 
            N = N+1;
       end
    end
    basis_functions = basis_functions(1:N,:);
    basis_coeff = basis_coeff(:,1:N);
end
% =========================================================================
%                  PLOTTING CONSTELLATION
function const_diag(basis_coeff)
    [M, N] = size(basis_coeff);
    figure('Name', 'Constellation Diagram', 'NumberTitle','off', 'Color',[1 1 1]);
    colors_c = lines(M);
    hold on; grid on;
    
    if N == 1
        for m = 1 : M
            scatter(basis_coeff(m,1), 0, 120, colors_c(m,:), 'filled');
            text(basis_coeff(m,1), 0.08, sprintf('s_{%d}',m), 'HorizontalAlignment', 'center', 'Color', colors_c(m,:), 'FontWeight', 'bold');
        end
        xlabel('\phi_1'); ylim([-1 1]);
        xline(0,'k--');
    elseif N == 2
        for m = 1 : M
            scatter(basis_coeff(m,1), basis_coeff(m,2), 120, colors_c(m,:), 'filled');
            text(basis_coeff(m,1) + 0.05, basis_coeff(m,2) + 0.05, sprintf('s_{%d}',m), 'Color', colors_c(m,:), 'FontWeight', 'bold');
        end
        xlabel('\phi_1'); ylabel('\phi_2');
        xline(0,'k--'); yline(0,'k--');
    else
        for m = 1 : M
            scatter3(basis_coeff(m,1), basis_coeff(m,2), basis_coeff(m,3), 120, colors_c(m,:), 'filled');
            text(basis_coeff(m,1), basis_coeff(m,2), basis_coeff(m,3), sprintf('  s_{%d}',m), 'Color', colors_c(m,:), 'FontWeight', 'bold');
        end
        xlabel('\phi_1'); ylabel('\phi_2'); zlabel('\phi_3');
        view(35, 25);
    end
    title(sprintf('Constellation Diagram (%d-D)', min(N, 3)), 'FontSize', 12);
end