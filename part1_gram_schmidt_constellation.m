clc; clear; close all;
GramSchmidt_GUI();  % launch the GUI (call was missing in the submitted version)

%%MINI GUI
function GramSchmidt_GUI()
  
    fig = uifigure('Name', 'Gram-Schmidt & Constellation GUI', 'Position', [100 100 500 550]);

    
    uilabel(fig, 'Position', [20 500 250 22], 'Text', 'Max time duration of signals (seconds):');
    time_edit = uieditfield(fig, 'numeric', 'Position', [280 500 100 22], 'Value', 4);

    uilabel(fig, 'Position', [20 460 250 22], 'Text', 'Number of signals (M):');
    num_edit = uieditfield(fig, 'numeric', 'Position', [280 460 100 22], 'Value', 4);

    uilabel(fig, 'Position', [20 420 450 22], 'Text', 'Signal expressions (one per line):');
    expr_area = uitextarea(fig, 'Position', [20 280 450 130], ...
        'Value', { ...
            '2.*( t<1) - 1.*(t>=1 & t<4)'; ... 
            'sin(2*pi*t).*(t<2.5) + cos(2*pi*t).*(t>=2.5)'; ...
            '4.*(t<1) - 3.*(t>=1 & t<2) + 1.*(t>=2)'; ...
            'exp(t)' ...
        });

 
    btn = uibutton(fig, 'push', 'Position', [175 230 150 35], ...
        'Text', 'Calculate & Plot', 'FontWeight', 'bold');

  
    uilabel(fig, 'Position', [20 190 200 22], 'Text', 'Energy Calculations:');
    result_area = uitextarea(fig, 'Position', [20 20 450 160], 'Editable', 'off');

    btn.ButtonPushedFcn = @(btn,event) run_calculations(time_edit, num_edit, expr_area, result_area);
end

function run_calculations(time_edit, num_edit, expr_area, result_area)
    No_sym = time_edit.Value;
    M = num_edit.Value;
    exprs = expr_area.Value;

    if length(exprs) < M, uialert(time_edit.Parent, 'Provide enough signal expressions.', 'Error'); return; end

    ts = 0.001;
    t = 0: ts: No_sym - ts;
    L = length(t);
    signals_matrix = zeros(M, L);

    for i = 1:M
        try
            expr = exprs{i};
            signals_matrix(i,:) = eval(expr);
        catch ME
            uialert(time_edit.Parent, sprintf('Error in signal %d: %s', i, ME.message), 'Error'); return;
        end
    end

    my_plot(signals_matrix,'s');
    [basis_func,basis_coeff] = gram_schmidt(signals_matrix,ts);
    my_plot(basis_func,'phi');
    
    result_text = {};
    for i = 1:M
        E_time = energy(signals_matrix(i,:), ts);
        E_const = energy_constellation(basis_coeff(i,:));
        
        result_text{end+1} = sprintf('--- Signal %d ---', i);
        result_text{end+1} = sprintf('Energy (Time Domain): %f', E_time);
        result_text{end+1} = sprintf('Energy (Constellation): %f', E_const);
        result_text{end+1} = ''; 
    end
    result_area.Value = result_text;

    const_diag(basis_coeff); 
end

% =========================================================================
%                CONSTELLATION FUNCTION 

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

% =========================================================================
%                  GRAM SCHMIDT FUNCTION


function [basis_functions,basis_coeff] = gram_schmidt(signals_matrix,ts)
    M = size(signals_matrix, 1);  
    L = size(signals_matrix, 2);
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
        g_i = signals_matrix(i,:) - basis_coeff(i,1:N)*basis_functions(1:N,:);
        g_i(abs(g_i)<1e-10) = 0;
        if(max(abs(g_i)) < 1e-10)
           continue;
        else
            N = N+1;
            E_gi = energy(g_i,ts);
            basis_functions(N,:) = (g_i)/sqrt(E_gi);
            basis_coeff(i,N) = sqrt(E_gi);
       end
    end
    basis_functions = basis_functions(1:N,:); 
     basis_coeff = basis_coeff(:,1:N);
end
% =========================================================================
%                  MY PLOT FUNCTION

function my_plot(signals_matrix, symbol_str)
    M = size(signals_matrix, 1);
     ts = 0.001;
    t = 0:ts:(size(signals_matrix, 2)-1)*ts;
    num_rows = ceil(sqrt(M)); num_cols = ceil(M / num_rows);
    main_fig = figure('Name', symbol_str); set(main_fig, 'color', 'w'); 
    for i = 1:M
        ax = subplot(num_rows, num_cols, i);
        plot(ax, t, signals_matrix(i, :), 'b', 'LineWidth', 1.5);
        grid(ax, 'on'); ax.Box = 'on'; ax.XColor = [0 0 0]; ax.YColor = [0 0 0];
        if(symbol_str == "phi"), title_str = sprintf('$\\%s_{%d}(t)$', symbol_str, i);
        else, title_str = sprintf('$%s_{%d}(t)$', symbol_str, i); end
        title(ax, title_str, 'Interpreter', 'latex', 'FontSize', 12);
        xlabel(ax, 'Time (s)'); ylabel(ax, 'Amplitude');
        ylim(ax, [min(signals_matrix(:))-0.5, max(signals_matrix(:))+0.5]);
    end
end
% =========================================================================
%                 ENERGY USING TIME DOMAIN FUNCTION
function E = energy(signal,ts), E = sum(abs(signal).^2)*ts; end

% =========================================================================
%                 ENERGY USING CONSTELLATION FUNCTION

function E = energy_constellation(basis_coeff), E = sum(basis_coeff.^2); end