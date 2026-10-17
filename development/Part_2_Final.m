clc; clear; close all;

%% Setup & Generation of Base Waveforms [0, Tb]
ts = 0.001;
T_sym = 1;
n_per_sym = T_sym / ts;
t = 0:ts:T_sym - ts;

% Unipolar NRZ 
A = 1;
s0_unrz = zeros(1, n_per_sym);
s1_unrz = A * ones(1, n_per_sym);
matrix_unrz = [s1_unrz; s0_unrz];

% Manchester
half_n = n_per_sym / 2;
s0_man = [-A*ones(1, half_n), A*ones(1, half_n)];
s1_man = [A*ones(1, half_n), -A*ones(1, half_n)];
matrix_man = [s1_man; s0_man];

%% Gram-Schmidt & Constellation Diagram 
[basis_unrz, coeff_unrz] = gram_schmidt(matrix_unrz, ts);
[basis_man, coeff_man] = gram_schmidt(matrix_man, ts);

const_diag(coeff_unrz);
title('Constellation Diagram: Unipolar NRZ');

const_diag(coeff_man);
title('Constellation Diagram: Manchester');

%% Generate Random Data & Encode
num_bits = 100000;
bits = randi([0, 1], 1, num_bits);

tx_unrz = zeros(num_bits, size(coeff_unrz, 2));
tx_unrz(bits == 0, :) = repmat(coeff_unrz(2,:), sum(bits == 0), 1);
tx_unrz(bits == 1, :) = repmat(coeff_unrz(1,:), sum(bits == 1), 1);

tx_man = zeros(num_bits, size(coeff_man, 2));
tx_man(bits == 0, :) = repmat(coeff_man(2,:), sum(bits == 0), 1);
tx_man(bits == 1, :) = repmat(coeff_man(1,:), sum(bits == 1), 1);

%% AWGN Channel
EbNo_dB = -10:2:10;
EbNo_lin = 10.^(EbNo_dB / 10);

Eb_unrz = mean(sum(coeff_unrz.^2,2));
Eb_man = mean(sum(coeff_man.^2,2));

ber_unrz_sim = zeros(1, length(EbNo_dB));
ber_man_sim = zeros(1, length(EbNo_dB));

fig_rx_unrz = figure('Name', 'Received Constellation Diagram- Unipolar NRZ');
fig_rx_man = figure('Name', 'Received Constellation Diagram - Manchester');

for k = 1:length(EbNo_dB)
    % Unipolar NRZ Channel
    N0_unrz = Eb_unrz / EbNo_lin(k);
    noise_unrz = sqrt(N0_unrz/2) * randn(size(tx_unrz));
    rx_unrz = tx_unrz + noise_unrz;
    
    % Manchester Channel
    N0_man = Eb_man / EbNo_lin(k);
    noise_man = sqrt(N0_man/2) * randn(size(tx_man));
    rx_man = tx_man + noise_man;
    
    %% Sketch the Constellation Diagram of Received Signals
        figure(fig_rx_unrz); subplot(3, 4, k);
        scatter(rx_unrz, zeros(size(rx_unrz)), 1, 'filled', 'MarkerFaceAlpha', 0.5);
        
        xline(0, 'k--', 'LineWidth', 1);
        title(sprintf('Eb/No = %d dB', EbNo_dB(k)));
        
        figure(fig_rx_man); subplot(3, 4, k);
        scatter(rx_man, zeros(size(rx_man)), 1, 'filled', 'MarkerFaceAlpha', 0.5);
        yticks([]);
        xline(0, 'k--', 'LineWidth', 1);
        title(sprintf('Eb/No = %d dB', EbNo_dB(k)));

    %% Decode Data (Minimum Euclidean Distance / Correlation Receiver)
    % Find distance between received symbol and Ideal Basis Coefficients
    dist0_unrz = sum((rx_unrz - coeff_unrz(2,:)).^2, 2);
    dist1_unrz = sum((rx_unrz - coeff_unrz(1,:)).^2, 2);
    rx_bits_unrz = zeros(1, num_bits);
    rx_bits_unrz(dist1_unrz < dist0_unrz) = 1;
    
    dist0_man = sum((rx_man - coeff_man(2,:)).^2, 2);
    dist1_man = sum((rx_man - coeff_man(1,:)).^2, 2);
    rx_bits_man = zeros(1, num_bits);
    rx_bits_man(dist1_man < dist0_man) = 1;

    %% Calculate Simulated BER
    ber_unrz_sim(k) = sum(rx_bits_unrz ~= bits) / num_bits;
    ber_man_sim(k) = sum(rx_bits_man ~= bits) / num_bits;
end
%% Theoretical BER Calculation

% Unipolar NRZ theoretical BER is Q(sqrt(Eb/(2*N0)))
ber_unrz_th = qfunc(sqrt(EbNo_lin));

% Manchester theoretical BER is Q(sqrt(2*Eb/No))
ber_man_th = qfunc(sqrt(2 * EbNo_lin));   

%% 6. Plot the Final Results (BER vs Eb/No)
figure('Name', 'BER Performance Curves');
semilogy(EbNo_dB, ber_unrz_sim, 'bo', 'LineWidth', 1.5, 'MarkerSize', 6); hold on;
semilogy(EbNo_dB, ber_unrz_th, 'b-', 'LineWidth', 1.5);
semilogy(EbNo_dB, ber_man_sim, 'rx', 'LineWidth', 1.5, 'MarkerSize', 6);
semilogy(EbNo_dB, ber_man_th, 'r-', 'LineWidth', 1.5);

grid on;
xlabel('E_b/N_o (dB)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Bit Error Rate (BER)', 'FontSize', 12, 'FontWeight', 'bold');
title('Theoretical vs Practical BER: Unipolar NRZ and Manchester', 'FontSize', 14);
legend('Unipolar NRZ (Simulated)', 'Unipolar NRZ (Theoretical)', ...
       'Manchester (Simulated)', 'Manchester (Theoretical)', 'Location', 'southwest');