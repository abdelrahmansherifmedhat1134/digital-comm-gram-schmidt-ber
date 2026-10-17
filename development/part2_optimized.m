% ECE352_Project_Phase2.m
% ---------------------------------------------------------
% Main Script Execution for System Simulation & BER
% ---------------------------------------------------------
clc; clear; close all;

% 1. System Parameters (GUI Inputs)
A = 1;              % Amplitude
Tb = 1;             % Bit duration
dt = 0.001;         % Time step
num_bits = 100000;  % Number of bits to simulate
EbNo_dB = -10:2:10; % Eb/No sweep range

% 2. Generate Line Codes & Link to Phase 1 Gram-Schmidt
[t, sigs_unrz, sigs_man] = generate_line_codes(A, Tb, dt);

% Phase 1 Link: Unipolar NRZ
[basis_unrz, N_unrz] = apply_gram_schmidt(sigs_unrz, dt);
[coords_unrz, ~, ~]  = calc_constellation(sigs_unrz, basis_unrz, dt);

% Phase 1 Link: Manchester
[basis_man, N_man]   = apply_gram_schmidt(sigs_man, dt);
[coords_man, ~, ~]   = calc_constellation(sigs_man, basis_man, dt);

% 3. Generate Random Bits & Map to Constellation
bits = randi([0, 1], 1, num_bits);
tx_coords_unrz = modulate_data(bits, coords_unrz);
tx_coords_man  = modulate_data(bits, coords_man);

% 4. Simulate AWGN Channel & Receiver (Signal Space)
[ber_sim_unrz, rx_clouds_unrz] = simulate_awgn_channel(bits, tx_coords_unrz, coords_unrz, EbNo_dB);
[ber_sim_man, rx_clouds_man]   = simulate_awgn_channel(bits, tx_coords_man, coords_man, EbNo_dB);

% 5. Calculate Theoretical BER
ber_th_unrz = calculate_theoretical_ber(EbNo_dB, 'unrz');
ber_th_man  = calculate_theoretical_ber(EbNo_dB, 'manchester');

% 6. Plot Results
plot_noisy_constellation(rx_clouds_unrz, EbNo_dB, 'Unipolar NRZ');
plot_noisy_constellation(rx_clouds_man, EbNo_dB, 'Manchester');
plot_ber_curves(EbNo_dB, ber_sim_unrz, ber_th_unrz, ber_sim_man, ber_th_man);

% ---------------------------------------------------------
% Include Phase 1 Functions Here (apply_gram_schmidt, calc_constellation)
% ---------------------------------------------------------


%% ----------------------------------------------------------------
%  function applying gram schmidt
function [basis, N] = apply_gram_schmidt(signals, dt)
   
    [M, L] = size(signals);
    basis_raw = zeros(M, L);
    N = 0;
    
    for m = 1 : M
        v = signals(m, :);
        for n = 1 : N
           
            inner = trapz(v .* basis_raw(n,:)) * dt;
            v = v - inner * basis_raw(n,:);
        end
        E_v = trapz(v .* v) * dt;
       
        if E_v > 1e-8
            N = N + 1;
            basis_raw(N, :) = v / sqrt(E_v);
        end
    end
    basis = basis_raw(1:N, :); 
end

%% ----------------------------------------------------------------
%  function calculating constellation values 
function [coords, E_symbol, E_direct] = calc_constellation(signals, basis, dt)
   
    [M, ~] = size(signals);
    N = size(basis, 1);
    coords = zeros(M, N);
    E_symbol = zeros(1, M);
    E_direct = zeros(1, M);

    for m = 1 : M
        for n = 1 : N
            coords(m, n) = trapz(signals(m,:) .* basis(n,:)) * dt;
        end
       
        E_symbol(m) = sum(coords(m,:).^2);
       
        E_direct(m) = trapz(signals(m,:).^2) * dt;
    end
end

%% ----------------------------------------------------------------
%  function Line Code Generator
function [t, sigs_unrz, sigs_man] = generate_line_codes(A, Tb, dt)
    t = 0 : dt : Tb - dt;
    n_per_sym = length(t);
    half_n = floor(n_per_sym / 2);
    
    % Unipolar NRZ: Row 1 = Logic 1, Row 2 = Logic 0
    s1_unrz = A * ones(1, n_per_sym);
    s0_unrz = zeros(1, n_per_sym);
    sigs_unrz = [s1_unrz; s0_unrz];
    
    % Manchester: Row 1 = Logic 1, Row 2 = Logic 0
    s1_man = [A * ones(1, half_n), -A * ones(1, n_per_sym - half_n)];
    s0_man = [-A * ones(1, half_n), A * ones(1, n_per_sym - half_n)];
    sigs_man = [s1_man; s0_man];
end

%% ----------------------------------------------------------------
%  function Data Modulator
function tx_coords = modulate_data(bits, ideal_coords)
    num_bits = length(bits);
    num_dims = size(ideal_coords, 2);
    tx_coords = zeros(num_bits, num_dims);
    
    % In our matrices, Row 1 is Logic 1, Row 2 is Logic 0
    tx_coords(bits == 1, :) = repmat(ideal_coords(1, :), sum(bits == 1), 1);
    tx_coords(bits == 0, :) = repmat(ideal_coords(2, :), sum(bits == 0), 1);
end

%% ----------------------------------------------------------------
%  function AWGN Simulator & Optimal Receiver
function [ber_sim, rx_clouds] = simulate_awgn_channel(bits, tx_coords, ideal_coords, EbNo_dB)
    num_bits = length(bits);
    num_dims = size(ideal_coords, 2);
    EbNo_lin = 10.^(EbNo_dB / 10);
    
    % Calculate Average Energy per Bit
    Eb = mean(sum(ideal_coords.^2, 2)); 
    
    ber_sim = zeros(1, length(EbNo_dB));
    rx_clouds = cell(1, length(EbNo_dB)); % Store noisy data for plotting later
    
    for k = 1:length(EbNo_dB)
        % Add AWGN
        N0 = Eb / EbNo_lin(k);
        noise = sqrt(N0 / 2) * randn(size(tx_coords));
        rx_coords = tx_coords + noise;
        
        % Save the first 1000 noisy points for the constellation scatter plot (GUI efficiency)
        rx_clouds{k} = rx_coords(1:min(1000, num_bits), :);
        
        % Minimum Distance Decoder (Optimal Receiver)
        dist1 = sum((rx_coords - ideal_coords(1, :)).^2, 2); % Dist to Logic 1
        dist0 = sum((rx_coords - ideal_coords(2, :)).^2, 2); % Dist to Logic 0
        
        rx_bits = zeros(1, num_bits);
        rx_bits(dist1 < dist0) = 1;
        
        % Calculate Errors
        ber_sim(k) = sum(rx_bits ~= bits) / num_bits;
    end
end

%% ----------------------------------------------------------------
%  function Theoretical Math
function ber_th = calculate_theoretical_ber(EbNo_dB, code_type)
    EbNo_lin = 10.^(EbNo_dB / 10);
    
    if strcmp(code_type, 'unrz')
        % Unipolar NRZ Theoretical: Q(sqrt(Eb/N0)) -- Note: Formula adjusted from your draft
        ber_th = qfunc(sqrt(EbNo_lin)); 
    elseif strcmp(code_type, 'manchester')
        % Manchester (Antipodal) Theoretical: Q(sqrt(2*Eb/N0))
        ber_th = qfunc(sqrt(2 * EbNo_lin));
    else
        error('Unknown line code type.');
    end
end

%% ----------------------------------------------------------------
%  function plotting 
function plot_noisy_constellation(rx_clouds, EbNo_dB, code_name)
    figure('Name', ['Received Constellation - ' code_name], 'Color', 'w');
    for k = 1:length(EbNo_dB)
        subplot(3, 4, k);
        
        % Check if 1D or 2D (Usually Gram-Schmidt makes these 1D)
        dims = size(rx_clouds{k}, 2);
        if dims == 1
            scatter(rx_clouds{k}(:,1), zeros(size(rx_clouds{k},1), 1), 10, 'filled', 'MarkerFaceAlpha', 0.5);
        else
            scatter(rx_clouds{k}(:,1), rx_clouds{k}(:,2), 10, 'filled', 'MarkerFaceAlpha', 0.5);
        end
        
        xline(0, 'k--', 'LineWidth', 1);
        title(sprintf('Eb/No = %d dB', EbNo_dB(k)));
        grid on; yticks([]);
    end
    sgtitle(sprintf('Noisy Received Constellation Clouds: %s', code_name), 'FontWeight', 'bold');
end

function plot_ber_curves(EbNo_dB, ber_sim_unrz, ber_th_unrz, ber_sim_man, ber_th_man)
    figure('Name', 'BER Performance Curves', 'Color', 'w');
    semilogy(EbNo_dB, ber_sim_unrz, 'bo', 'LineWidth', 1.5, 'MarkerSize', 6); hold on;
    semilogy(EbNo_dB, ber_th_unrz, 'b-', 'LineWidth', 1.5);
    semilogy(EbNo_dB, ber_sim_man, 'rx', 'LineWidth', 1.5, 'MarkerSize', 6);
    semilogy(EbNo_dB, ber_th_man, 'r-', 'LineWidth', 1.5);
    
    grid on;
    xlabel('E_b/N_o (dB)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Bit Error Rate (BER)', 'FontSize', 12, 'FontWeight', 'bold');
    title('Theoretical vs Simulated BER Performance', 'FontSize', 14);
    legend('Unipolar NRZ (Simulated)', 'Unipolar NRZ (Theoretical)', ...
           'Manchester (Simulated)', 'Manchester (Theoretical)', 'Location', 'southwest');
    ylim([1e-5, 1]); % Standardize graph bounds
end