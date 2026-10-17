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
    % Global Time Settings
    ans_time = inputdlg({'Total Simulation Time (s):', 'Time Step (dt):'}, ...
                        'Global Settings', 1, {'4', '0.001'});
    if isempty(ans_time), error('User canceled the input.'); end
    T_total = str2double(ans_time{1});
    dt = str2double(ans_time{2});

    sig_defs = {}; 
    M = 0;         
    
    while true
        main_choice = menu(sprintf('Main Menu (Current Signals: %d)', M), '1. Add a New Signal', '2. FINISH AND RUN');
        if main_choice == 2 || main_choice == 0
            if M == 0, error('You must add at least one signal!'); end
            break; 
        end

        M = M + 1;
        sig_defs{M}.name = sprintf('s_%d(t)', M);
        sig_defs{M}.segments = {}; 
        seg_count = 0;
        
        while true
            seg_choice = menu(sprintf('Building Signal %d (Segments: %d)', M, seg_count), ...
                              '1. Add a Flat / DC Segment', ...
                              '2. Add a Sine / Cosine Segment', ...
                              '3. Add an Exponential Segment', ...
                              '4. Add a Custom Math Equation (eval)', ... % <--- THE BRIDGE
                              '5. DONE WITH THIS SIGNAL');
                          
            if seg_choice == 5 || seg_choice == 0, break; end
            
            % Ask WHERE this segment should live in time
            ans_bounds = inputdlg({'Start Time (s):', 'End Time (s):'}, ...
                                  sprintf('Time Bounds for Seg %d', seg_count+1), 1, {'0', num2str(T_total)});
            if isempty(ans_bounds), continue; end
            
            seg_count = seg_count + 1;
            seg = struct(); 
            seg.t_start = str2double(ans_bounds{1});
            seg.t_end   = str2double(ans_bounds{2});
            
            % Ask for parameters based on choice
            if seg_choice == 1
                ans_dc = inputdlg({'Amplitude (V):'}, 'Flat/DC Input', 1, {'1'});
                if ~isempty(ans_dc)
                    seg.type = 'dc'; seg.A = str2double(ans_dc{1});
                end
            elseif seg_choice == 2
                ans_sin = inputdlg({'Amplitude (A):', 'Frequency in Hz (f):', 'Time Shift (t0):', 'DC Shift:'}, ...
                                   'Sine Wave Input', 1, {'2', '1', '0', '0'});
                if ~isempty(ans_sin)
                    seg.type = 'sin'; seg.A = str2double(ans_sin{1}); seg.f = str2double(ans_sin{2});
                    seg.t0 = str2double(ans_sin{3});  seg.dc = str2double(ans_sin{4});
                end
            elseif seg_choice == 3
                ans_exp = inputdlg({'Amplitude (A):', 'Exponent rate (alpha):'}, 'Exponential Input', 1, {'1', '-2'});
                if ~isempty(ans_exp)
                    seg.type = 'exp'; seg.A = str2double(ans_exp{1}); seg.alpha = str2double(ans_exp{2});
                end
            elseif seg_choice == 4
                % --- THE BRIDGE INPUT ---
                ans_eq = inputdlg({'Enter MATLAB equation (use "t" as the variable):', 'Example: t.^2 + 2*t'}, ...
                                  'Custom Math Input', 1, {'sawtooth(2*pi*1*t)'});
                if ~isempty(ans_eq)
                    seg.type = 'custom'; seg.eq = ans_eq{1};
                end
            end
            
            % Save the segment if one was successfully created
            if isfield(seg, 'type')
                sig_defs{M}.segments{seg_count} = seg;
            end
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
        signal_names{m} = sig_defs{m}.name;
        segs = sig_defs{m}.segments; 
        
        for k = 1 : length(segs)
            seg = segs{k};
            
            % Create the logical mask for the specific time boundaries
            mask = (t >= seg.t_start & t < seg.t_end);
            t_window = t(mask); 
            
            if strcmp(seg.type, 'dc')
                signals(m, mask) = signals(m, mask) + seg.A;
                
            elseif strcmp(seg.type, 'sin')
                signals(m, mask) = signals(m, mask) + (seg.A * sin(2 * pi * seg.f * (t_window - seg.t0)) + seg.dc);
                
            elseif strcmp(seg.type, 'exp')
                signals(m, mask) = signals(m, mask) + (seg.A * exp(seg.alpha * t_window));
                
            elseif strcmp(seg.type, 'custom')
                % --- THE BRIDGE EXECUTION ---
                % 1. Evaluate the user's string over the ENTIRE 't' vector
                full_wave = eval(seg.eq); 
                
                % 2. Catch if the user typed a static number (like "5") instead of an equation with 't'
                if length(full_wave) == 1 
                    full_wave = full_wave * ones(1, L);
                end
                
                % 3. Apply the mask so it only cuts out the piece they requested!
                signals(m, mask) = signals(m, mask) + full_wave(mask);
            end
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