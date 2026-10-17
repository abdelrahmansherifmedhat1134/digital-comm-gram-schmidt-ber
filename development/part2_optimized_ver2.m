% ECE352_Project_Phase2_GUI.m
% ---------------------------------------------------------
% Main Script Execution for System Simulation & BER
% ---------------------------------------------------------
clc; clear; close all;

% 1. System Parameters & GUI Menus
A = 1; Tb = 1; dt = 0.001;
[num_bits, EbNo_dB, selected_codes] = get_phase2_inputs();

% 2. Generate Random Bits ONCE (So all line codes are tested fairly)
bits = randi([0, 1], 1, num_bits);

% 3. Prepare the Master BER Plot
fig_ber = figure('Name', 'Master BER Performance Curves', 'Color', 'w');
hold on; grid on;
colors = lines(length(selected_codes)); % Generate distinct colors
legend_entries = {};

% 4. Loop through every selected Line Code
for c = 1:length(selected_codes)
    code_name = selected_codes{c};
    
    % A. Generate Waveforms
    [t, sigs] = generate_line_codes(A, Tb, dt, code_name);

    % B. Phase 1 Link (Gram-Schmidt)
    [basis, N]     = apply_gram_schmidt(sigs, dt);
    [coords, ~, ~] = calc_constellation(sigs, basis, dt);

    % C. Modulate Data
    tx_coords = modulate_data(bits, coords);

    % D. Simulate AWGN Channel & Receiver
    [ber_sim, rx_clouds] = simulate_awgn_channel(bits, tx_coords, coords, EbNo_dB);

    % E. Calculate Theoretical BER
    ber_th = calculate_theoretical_ber(EbNo_dB, code_name);

    % F. Plot the Noisy Constellations for this specific code
    plot_noisy_constellation(rx_clouds, EbNo_dB, code_name);

    % G. Add the curves to the Master BER Plot
    figure(fig_ber); % Focus on the BER plot
    semilogy(EbNo_dB, ber_sim, 'o', 'Color', colors(c,:), 'LineWidth', 1.5);
    semilogy(EbNo_dB, ber_th, '-', 'Color', colors(c,:), 'LineWidth', 1.5);
    
    legend_entries{end+1} = sprintf('%s (Sim)', code_name);
    legend_entries{end+1} = sprintf('%s (Theory)', code_name);
end

% 5. Finalize Master BER Plot Graphics
figure(fig_ber);
xlabel('E_b/N_o (dB)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Bit Error Rate (BER)', 'FontSize', 12, 'FontWeight', 'bold');
title('Theoretical vs Simulated BER Performance', 'FontSize', 14);
legend(legend_entries, 'Location', 'southwest');
ylim([1e-5, 1]);



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
function [t, sigs] = generate_line_codes(A, Tb, dt, code_name)
    t = 0 : dt : Tb - dt;
    n_per_sym = length(t);
    half_n = floor(n_per_sym / 2);
    
    switch code_name
        case 'Unipolar NRZ'
            s1 = A * ones(1, n_per_sym);
            s0 = zeros(1, n_per_sym);
        case 'Polar NRZ'
            s1 = A * ones(1, n_per_sym);
            s0 = -A * ones(1, n_per_sym);
        case 'Unipolar RZ'
            s1 = [A * ones(1, half_n), zeros(1, n_per_sym - half_n)];
            s0 = zeros(1, n_per_sym);
        case 'Polar RZ'
            s1 = [A * ones(1, half_n), zeros(1, n_per_sym - half_n)];
            s0 = [-A * ones(1, half_n), zeros(1, n_per_sym - half_n)];
        case 'Manchester'
            s1 = [A * ones(1, half_n), -A * ones(1, n_per_sym - half_n)];
            s0 = [-A * ones(1, half_n), A * ones(1, n_per_sym - half_n)];
        case 'ASK (OOK)'
            fc = 5 / Tb; % 5 Carrier cycles per bit
            s1 = A * sin(2*pi*fc*t);
            s0 = zeros(1, n_per_sym);
        case 'BPSK'
            fc = 5 / Tb; % 5 Carrier cycles per bit
            s1 = A * sin(2*pi*fc*t);
            s0 = -A * sin(2*pi*fc*t);
        otherwise
            error('Unknown line code selected.');
    end
    sigs = [s1; s0]; % Package into the Phase 1 Matrix
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
function ber_th = calculate_theoretical_ber(EbNo_dB, code_name)
    EbNo_lin = 10.^(EbNo_dB / 10);
    
    % Group 1: Orthogonal Signals (Distance = sqrt(Eb))
    if ismember(code_name, {'Unipolar NRZ', 'Unipolar RZ', 'ASK (OOK)'})
        ber_th = qfunc(sqrt(EbNo_lin)); 
        
    % Group 2: Antipodal Signals (Distance = 2*sqrt(Eb)) -> Much better performance!
    elseif ismember(code_name, {'Polar NRZ', 'Polar RZ', 'Manchester', 'BPSK'})
        ber_th = qfunc(sqrt(2 * EbNo_lin));
        
    else
        error('Unknown line code for theoretical calculation.');
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


%% ----------------------------------------------------------------
%  function drop menu mini gui 
function [num_bits, EbNo_dB, selected_codes] = get_phase2_inputs()
    % 1. Numeric Inputs Form
    prompt = {'Number of bits to simulate:', 'Eb/No Start (dB):', 'Eb/No Step (dB):', 'Eb/No End (dB):'};
    dlgtitle = 'Simulation Parameters';
    definput = {'100000', '-10', '2', '10'};
    ans_num = inputdlg(prompt, dlgtitle, 1, definput);
    
    if isempty(ans_num), error('Simulation cancelled by user.'); end
    
    num_bits = str2double(ans_num{1});
    start_db = str2double(ans_num{2});
    step_db  = str2double(ans_num{3});
    end_db   = str2double(ans_num{4});
    EbNo_dB  = start_db : step_db : end_db;
    
    % 2. Multi-Select Line Code Menu
    list_codes = {'Unipolar NRZ', 'Polar NRZ', 'Unipolar RZ', 'Polar RZ', 'Manchester', 'ASK (OOK)', 'BPSK'};
    [indx, tf] = listdlg('PromptString', 'Select Line Codes (Hold Ctrl to select multiple):', ...
                         'SelectionMode', 'multiple', ...
                         'ListSize', [250, 150], ...
                         'ListString', list_codes);
                     
    if ~tf, error('Simulation cancelled by user.'); end
    selected_codes = list_codes(indx);
end