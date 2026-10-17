% ECE352_Project_Phase1_Interactive.m
% ---------------------------------------------------------
% Main Script Execution
% ---------------------------------------------------------
clc; clear; close all;

% 1. Open the interactive forms to get custom signals
[sig_defs, T_total, dt] = get_user_inputs();

% 2. Process Data using the new generator
[signals, t, signal_names] = generate_signals(sig_defs, T_total, dt);
[basis, N]                 = apply_gram_schmidt(signals, dt);
[coords, E_symbol, E_dir]  = calc_constellation(signals, basis, dt);

% 3. Output Verification to Console
fprintf('\nGenerated %d basis functions from %d signals.\n', N, length(sig_defs));

% 4. Plot the results
plot_input_signals(t, signals, signal_names);
plot_basis_functions(t, basis);
plot_constellation(coords, signal_names, E_symbol, N);


%% ----------------------------------------------------------------
%  gui taking user inputs 
function [sig_defs, T_total, dt] = get_user_inputs()

    ans_time = inputdlg({'Total Simulation Time (s):', 'Time Step (dt):'}, ...
                        'Global Settings', 1, {'4', '0.001'});
    
    if isempty(ans_time)
        error('User canceled the input.');
    end
    
    T_total = str2double(ans_time{1});
    dt = str2double(ans_time{2});

    sig_defs = {}; 
    M = 0;         
    
   
    while true
       
        choice = menu(sprintf('Build Signal %d', M+1), ...
                      '1. Time-Slot Array (Original)', ...
                      '2. Sine / Cosine Wave', ...
                      '3. Exponential Curve', ...
                      '4. FINISH AND RUN');

        if choice == 1
            
            ans_slot = inputdlg({'Enter Amplitude Array (e.g., [2 2 -1 -1]):', 'Duration of ONE slot (Tb):'}, ...
                                'Time-Slot Input', 1, {'[2 2 -1 -1]', '1'});
            if ~isempty(ans_slot)
                M = M + 1;
                sig_defs{M}.type = 'slot';
                sig_defs{M}.arr = str2num(ans_slot{1});
                sig_defs{M}.Tb = str2double(ans_slot{2});
                sig_defs{M}.name = sprintf('s_%d (Slot)', M);
            end
            
        elseif choice == 2
           
            ans_sin = inputdlg({'Amplitude (A):', 'Frequency in Hz (f):', ...
                                'Time Shift in secs (t0):', 'Amplitude/DC Shift:'}, ...
                                'Sine Wave Input', 1, {'2', '1', '0', '0'});
            if ~isempty(ans_sin)
                M = M + 1;
                sig_defs{M}.type = 'sin';
                sig_defs{M}.A = str2double(ans_sin{1});
                sig_defs{M}.f = str2double(ans_sin{2});
                sig_defs{M}.t0 = str2double(ans_sin{3});
                sig_defs{M}.dc = str2double(ans_sin{4});
                sig_defs{M}.name = sprintf('s_%d (Sin)', M);
            end
            
        elseif choice == 3
            % Form for Exponential (Decaying or Rising)
            ans_exp = inputdlg({'Amplitude (A):', 'Exponent rate (alpha) [Negative=Decay, Positive=Rise]:'}, ...
                               'Exponential Input', 1, {'1', '-2'});
            if ~isempty(ans_exp)
                M = M + 1;
                sig_defs{M}.type = 'exp';
                sig_defs{M}.A = str2double(ans_exp{1});
                sig_defs{M}.alpha = str2double(ans_exp{2});
                sig_defs{M}.name = sprintf('s_%d (Exp)', M);
            end
            
        elseif choice == 4 || choice == 0
            % If they click Finish or close the menu
            if M == 0
                error('You must add at least one signal before running!');
            end
            break; % Exit the while loop and move on to generation
        end
    end
end

%% ----------------------------------------------------------------
%  function generating signals 
function [signals, t, signal_names] = generate_signals(sig_defs, T_total, dt)
  
    t = 0 : dt : T_total - dt;
    L = length(t);
    M = length(sig_defs);
    
   
    signals = zeros(M, L);
    signal_names = cell(1, M);
    
    
    for m = 1 : M
        def = sig_defs{m};          
        signal_names{m} = def.name; 
        
        switch def.type
            case 'slot'
                
                num_slots = length(def.arr);
                for k = 1 : num_slots
                    t_start = (k-1) * def.Tb;
                    t_end   =  k    * def.Tb;
                    signals(m, :) = signals(m, :) + def.arr(k) * (t >= t_start & t < t_end);
                end
                
            case 'sin'
                
                signals(m, :) = def.A * sin(2 * pi * def.f * (t - def.t0)) + def.dc;
                
            case 'exp'
                
                signals(m, :) = def.A * exp(def.alpha * t);
        end
    end
end


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
%  function plotting input signals
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

%% ----------------------------------------------------------------
%  function plotting basis functions 
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

%% ----------------------------------------------------------------
%  function plotting constellation diagram
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
    
    
    e_str = '';
    for m = 1 : M
        e_str = [e_str sprintf('%s: E = %.3f    ', signal_names{m}, E_symbol(m))]; 
    end
    text(0.02, -0.09, e_str, 'Units','normalized', 'FontSize', 8, 'Color',[0.3 0.3 0.3]);
end