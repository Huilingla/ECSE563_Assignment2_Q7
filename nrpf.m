function [V, delta, Ps1, Qgv, N, time] = nrpf(Y, is, ipq, ipv, Pg, Qg, Pd, Qd, V0, Sbase, toler, maxiter)
%NRPF Newton-Raphson Power Flow Solver
%
% Inputs:
%   Y       - System admittance matrix
%   is      - Slack bus index
%   ipq     - PQ bus indices
%   ipv     - PV bus indices
%   Pg      - Generator active power injections (MW)
%   Qg      - Generator reactive power injections (MVAR)
%   Pd      - Load active power demands (MW)
%   Qd      - Load reactive power demands (MVAR)
%   V0      - Initial voltage magnitudes (pu)
%   Sbase   - System base power (MVA)
%   toler   - Convergence tolerance (pu)
%   maxiter - Maximum number of iterations
%
% Outputs:
%   V       - Final voltage magnitudes (pu)
%   delta   - Final voltage angles (radians)
%   Ps1     - Slack bus active power (MW)
%   Qgv     - Reactive power generations at PV buses (MVAR)
%   N       - Number of iterations to convergence
%   time    - CPU time (seconds)

    tstart = tic;
    
    % Input validation and initialization
    nbus = size(Y, 1);
    Pg = Pg(:); Qg = Qg(:); Pd = Pd(:); Qd = Qd(:); V0 = V0(:);
    ipq = ipq(:); ipv = ipv(:);
    
    % Convert to per unit
    Pg_pu = Pg / Sbase;
    Qg_pu = Qg / Sbase;
    Pd_pu = Pd / Sbase;
    Qd_pu = Qd / Sbase;
    
    % Initialize variables
    V = V0;
    delta = zeros(nbus, 1);
    
    % Remove slack from PV buses
    ipv_noslack = setdiff(ipv, is);
    
    % Identify bus types
    npq = length(ipq);
    npv = length(ipv_noslack);
    non_slack = sort([ipv_noslack; ipq]);
    
    fprintf('Newton-Raphson Power Flow:\n');
    fprintf('  System: %d buses (%d PQ, %d PV, 1 Slack)\n', nbus, npq, npv);
    fprintf('  Tolerance: %.2e, Max iterations: %d\n', toler, maxiter);
    fprintf('  Unknown variables: %d\n\n', npq + npv + npq);
    
    % Main Newton-Raphson iteration loop
    for N = 1:maxiter
        % Calculate power injections
        [P_calc, Q_calc] = calculate_power(Y, V, delta);
        
        % Calculate power mismatches
        dP = (Pg_pu - Pd_pu) - P_calc;
        dQ = (Qg_pu - Qd_pu) - Q_calc;
        
        % Zero mismatches for specified quantities
        dP(is) = 0;
        dQ(is) = 0;
        dQ(ipv) = 0;
        
        % Create mismatch vector
        dP_vec = dP(non_slack);
        dQ_vec = dQ(ipq);
        mismatch = [dP_vec; dQ_vec];
        
        max_mismatch = max(abs(mismatch));
        fprintf('  Iter %2d: Max mismatch = %12.6f pu\n', N, max_mismatch);
        
        % Check convergence
        if max_mismatch < toler
            fprintf('  *** CONVERGED in %d iterations ***\n', N);
            break;
        end
        
        % Build Jacobian matrix
        J = build_jacobian(Y, V, delta, non_slack, ipq);
        
        % Solve for corrections
        corrections = linsolve(J, mismatch);
        
        % Extract corrections
        n_angles = length(non_slack);
        d_delta = corrections(1:n_angles);
        dV = corrections(n_angles+1:end);
        
        % Apply corrections
        delta(non_slack) = delta(non_slack) + d_delta;
        V(ipq) = V(ipq) + dV;
        
        % Maintain specified voltages
        V(ipv) = V0(ipv);
        V(is) = V0(is);
    end
    
    % Handle convergence
    if N == maxiter
        [P_calc, Q_calc] = calculate_power(Y, V, delta);
        dP = (Pg_pu - Pd_pu) - P_calc;
        dQ = (Qg_pu - Qd_pu) - Q_calc;
        dP(is) = 0; dQ(is) = 0; dQ(ipv) = 0;
        final_mismatch = max(max(abs(dP)), max(abs(dQ)));
        
        if final_mismatch > toler
            warning('Did not converge within %d iterations. Final mismatch: %.6f', maxiter, final_mismatch);
        end
    end
    
    % Calculate final outputs
    [P_calc, Q_calc] = calculate_power(Y, V, delta);
    Ps1 = P_calc(is) * Sbase;
    Qgv = Q_calc(ipv_noslack) * Sbase;
    
    time = toc(tstart);
    
    % Display results
    display_results(V, delta, Ps1, Qgv, N, time, is, ipv_noslack, ipq);
end

function [P_calc, Q_calc] = calculate_power(Y, V, delta)
% Calculate power injections using admittance matrix
    nbus = length(V);
    P_calc = zeros(nbus, 1);
    Q_calc = zeros(nbus, 1);
    
    for i = 1:nbus
        for j = 1:nbus
            theta_ij = delta(i) - delta(j);
            P_calc(i) = P_calc(i) + V(i) * V(j) * ...
                       (real(Y(i,j)) * cos(theta_ij) + imag(Y(i,j)) * sin(theta_ij));
            Q_calc(i) = Q_calc(i) + V(i) * V(j) * ...
                       (real(Y(i,j)) * sin(theta_ij) - imag(Y(i,j)) * cos(theta_ij));
        end
    end
end

function J = build_jacobian(Y, V, delta, non_slack, ipq)
% Build Jacobian matrix for Newton-Raphson method
    nbus = length(V);
    G = real(Y);
    B = imag(Y);
    
    n_angles = length(non_slack);
    n_voltages = length(ipq);
    J = zeros(n_angles + n_voltages, n_angles + n_voltages);
    
    % J11 = dP/dDelta
    for i = 1:n_angles
        m = non_slack(i);
        for j = 1:n_angles
            n = non_slack(j);
            if m == n
                % Diagonal element
                sum_term = 0;
                for k = 1:nbus
                    if k ~= m
                        theta_mk = delta(m) - delta(k);
                        sum_term = sum_term + V(m) * V(k) * ...
                                  (-G(m,k) * sin(theta_mk) + B(m,k) * cos(theta_mk));
                    end
                end
                J(i,j) = sum_term;
            else
                % Off-diagonal element
                theta_mn = delta(m) - delta(n);
                J(i,j) = V(m) * V(n) * (G(m,n) * sin(theta_mn) - B(m,n) * cos(theta_mn));
            end
        end
    end
    
    % J12 = dP/dV
    for i = 1:n_angles
        m = non_slack(i);
        for j = 1:n_voltages
            n = ipq(j);
            if m == n
                % Diagonal element
                sum_term = V(m) * G(m,m);
                for k = 1:nbus
                    theta_mk = delta(m) - delta(k);
                    sum_term = sum_term + V(k) * (G(m,k) * cos(theta_mk) + B(m,k) * sin(theta_mk));
                end
                J(i, n_angles + j) = sum_term;
            else
                % Off-diagonal element
                theta_mn = delta(m) - delta(n);
                J(i, n_angles + j) = V(m) * (G(m,n) * cos(theta_mn) + B(m,n) * sin(theta_mn));
            end
        end
    end
    
    % J21 = dQ/dDelta
    for i = 1:n_voltages
        m = ipq(i);
        for j = 1:n_angles
            n = non_slack(j);
            if m == n
                % Diagonal element
                sum_term = 0;
                for k = 1:nbus
                    if k ~= m
                        theta_mk = delta(m) - delta(k);
                        sum_term = sum_term + V(m) * V(k) * ...
                                  (G(m,k) * cos(theta_mk) + B(m,k) * sin(theta_mk));
                    end
                end
                J(n_angles + i, j) = sum_term;
            else
                % Off-diagonal element
                theta_mn = delta(m) - delta(n);
                J(n_angles + i, j) = -V(m) * V(n) * (G(m,n) * cos(theta_mn) + B(m,n) * sin(theta_mn));
            end
        end
    end
    
    % J22 = dQ/dV
    for i = 1:n_voltages
        m = ipq(i);
        for j = 1:n_voltages
            n = ipq(j);
            if m == n
                % Diagonal element
                sum_term = -V(m) * B(m,m);
                for k = 1:nbus
                    theta_mk = delta(m) - delta(k);
                    sum_term = sum_term + V(k) * (G(m,k) * sin(theta_mk) - B(m,k) * cos(theta_mk));
                end
                J(n_angles + i, n_angles + j) = sum_term;
            else
                % Off-diagonal element
                theta_mn = delta(m) - delta(n);
                J(n_angles + i, n_angles + j) = V(m) * (G(m,n) * sin(theta_mn) - B(m,n) * cos(theta_mn));
            end
        end
    end
end

function display_results(V, delta, Ps1, Qgv, N, time, is, ipv, ipq)
% Display final results
    fprintf('\n=== NEWTON-RAPHSON POWER FLOW RESULTS ===\n');
    
    fprintf('\nREQUIRED OUTPUTS:\n');
    fprintf('=================\n');
    
    % 1. Node Voltage Magnitudes and Angles
    fprintf('\n1. NODE VOLTAGE SOLUTION:\n');
    fprintf('Bus  Type   V(pu)      delta(rad)   delta(deg)\n');
    fprintf('---- -----  ----------  -----------  ----------\n');
    for i = 1:length(V)
        if i == is
            bus_type = 'SL';
        elseif ismember(i, ipv)
            bus_type = 'PV';
        else
            bus_type = 'PQ';
        end
        fprintf('%2d   %-2s    %10.6f  %11.6f  %10.3f\n', ...
            i, bus_type, V(i), delta(i), rad2deg(delta(i)));
    end
    
    % 2. Slack Bus Active Power
    fprintf('\n2. SLACK BUS ACTIVE POWER:\n');
    fprintf('   Ps1 = %.6f MW\n', Ps1);
    
    % 3. Reactive Power at Voltage-Controlled Nodes
    fprintf('\n3. REACTIVE POWER AT VOLTAGE-CONTROLLED NODES:\n');
    fprintf('Bus   Qgv(MVAR)\n');
    fprintf('----  ----------\n');
    for i = 1:length(ipv)
        fprintf('%2d    %10.6f\n', ipv(i), Qgv(i));
    end
    
    % 4. Iterations and CPU Time
    fprintf('\n4. CONVERGENCE INFORMATION:\n');
    fprintf('   Iterations (N): %d\n', N);
    fprintf('   CPU Time: %.6f seconds\n', time);
end

function deg = rad2deg(rad)
    deg = rad * 180 / pi;
end
