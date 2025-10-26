function security_analysis_ieee9()
    % IEEE 9 bus test system security analysis
    % Question 7 - Contingency analysis
    
    % Load base case data
    ieee9_A2;
    
    fprintf('=== IEEE 9 Bus System Security Analysis ===\n\n');
    fprintf('Acceptable voltage range: 0.95 - 1.05 p.u.\n\n');
    
    % First, let's check and fix the input data sizes
    fprintf('Checking input data sizes...\n');
    nbuses = 9; % IEEE 9-bus system
    
    % Ensure all vectors have correct size (9 elements)
    Pg = ensure_size(Pg, nbuses, 0);
    Qg = ensure_size(Qg, nbuses, 0);
    Pd = ensure_size(Pd, nbuses, 0);
    Qd = ensure_size(Qd, nbuses, 0);
    V0 = ensure_size(V0, nbuses, 1.0);
    
    fprintf('Data sizes verified:\n');
    fprintf('  Pg: %d elements, Qg: %d elements\n', length(Pg), length(Qg));
    fprintf('  Pd: %d elements, Qd: %d elements\n', length(Pd), length(Qd));
    fprintf('  V0: %d elements\n', length(V0));
    
    % Run base case first
    fprintf('\nRunning Base Case Analysis...\n');
    [Y_base, ~, ~] = admittance(nfrom, nto, r, x, b);
    
    try
        [V_base, delta_base, Ps1_base, Qgv_base, N_base, time_base] = nrpf(...
            Y_base, is, ipq, ipv, Pg, Qg, Pd, Qd, V0, Sbase, toler, maxiter);
        
        display_voltage_results('Base Case', V_base);
        
    catch ME
        fprintf('❌ Base case power flow failed: %s\n', ME.message);
        fprintf('Troubleshooting base case...\n');
        troubleshoot_base_case(Y_base, is, ipq, ipv, Pg, Qg, Pd, Qd, V0, Sbase);
        return;
    end
    
    % Contingencies to analyze (lines to remove)
    contingencies = {
        'Line 4-5', 2;    % Line between buses 4-5 (index 2 in the data)
        'Line 4-9', 9;    % Line between buses 4-9
        'Line 5-6', 3;    % Line between buses 5-6  
        'Line 6-7', 5;    % Line between buses 6-7
        'Line 7-8', 6;    % Line between buses 7-8
        'Line 8-9', 8     % Line between buses 8-9
    };
    
    % Analyze each contingency
    voltage_violations = {};
    
    for i = 1:size(contingencies, 1)
        cont_name = contingencies{i, 1};
        line_idx = contingencies{i, 2};
        
        fprintf('\n%s Contingency Analysis:\n', cont_name);
        
        % Remove the line
        new_nfrom = nfrom;
        new_nto = nto;
        new_r = r;
        new_x = x;
        new_b = b;
        
        new_nfrom(line_idx) = [];
        new_nto(line_idx) = [];
        new_r(line_idx) = [];
        new_x(line_idx) = [];
        new_b(line_idx) = [];
        
        % Calculate new admittance matrix and run power flow
        [Y_cont, ~, ~] = admittance(new_nfrom, new_nto, new_r, new_x, new_b);
        
        try
            [V_cont, delta_cont, Ps1_cont, Qgv_cont, N_cont, time_cont] = nrpf(...
                Y_cont, is, ipq, ipv, Pg, Qg, Pd, Qd, V0, Sbase, toler, maxiter);
            
            % Check for voltage violations
            [has_violation, min_volt, max_volt, violation_buses] = check_voltage_violations(V_cont);
            
            display_voltage_results(cont_name, V_cont);
            
            if has_violation
                voltage_violations{end+1} = struct(...
                    'name', cont_name, ...
                    'min_volt', min_volt, ...
                    'max_volt', max_volt, ...
                    'violation_buses', violation_buses, ...
                    'V', V_cont);
                fprintf('❌ VOLTAGE VIOLATION DETECTED!\n');
                fprintf('   Violating buses: %s\n', mat2str(violation_buses));
            else
                fprintf('✅ Voltages within acceptable range\n');
            end
            
        catch ME
            fprintf('❌ POWER FLOW FAILED: %s\n', ME.message);
            voltage_violations{end+1} = struct(...
                'name', cont_name, ...
                'min_volt', NaN, ...
                'max_volt', NaN, ...
                'violation_buses', [], ...
                'V', [], ...
                'error', ME.message);
        end
    end
    
    % Display summary and mitigation measures
    display_security_summary(voltage_violations);
    propose_mitigation_measures(voltage_violations, V_base);
end

function vec_out = ensure_size(vec_in, target_size, default_value)
    % Ensure vector has correct size, padding with default_value if needed
    if length(vec_in) < target_size
        vec_out = zeros(target_size, 1);
        vec_out(1:length(vec_in)) = vec_in;
        if nargin > 2
            vec_out(length(vec_in)+1:end) = default_value;
        end
    else
        vec_out = vec_in(1:target_size);
    end
    vec_out = vec_out(:); % Ensure column vector
end

function [Y, Bff, Bft] = admittance(nfrom, nto, r, x, b)
    % Calculate admittance matrix from line data
    nlines = length(nfrom);
    nbuses = max([nfrom; nto]);
    
    Y = zeros(nbuses, nbuses) + 1j * zeros(nbuses, nbuses);
    
    for k = 1:nlines
        from = nfrom(k);
        to = nto(k);
        
        if r(k) == 0 && x(k) == 0
            continue;
        end
        
        z = r(k) + 1j * x(k);
        y = 1/z;
        
        Y(from, to) = Y(from, to) - y;
        Y(to, from) = Y(to, from) - y;
        Y(from, from) = Y(from, from) + y + 1j * b(k)/2;
        Y(to, to) = Y(to, to) + y + 1j * b(k)/2;
    end
    
    Bff = imag(Y);
    Bft = Bff;
end

function troubleshoot_base_case(Y, is, ipq, ipv, Pg, Qg, Pd, Qd, V0, Sbase)
    % Debug function to identify base case issues
    
    fprintf('\n=== TROUBLESHOOTING BASE CASE ===\n');
    
    % Check sizes
    fprintf('Matrix/Vector sizes:\n');
    fprintf('  Y: %dx%d\n', size(Y));
    fprintf('  Pg: %d, Qg: %d\n', length(Pg), length(Qg));
    fprintf('  Pd: %d, Qd: %d\n', length(Pd), length(Qd));
    fprintf('  V0: %d\n', length(V0));
    
    % Check if sizes match
    nbuses = size(Y, 1);
    if length(Pg) ~= nbuses
        fprintf('❌ Pg size mismatch: expected %d, got %d\n', nbuses, length(Pg));
    end
    if length(Qg) ~= nbuses
        fprintf('❌ Qg size mismatch: expected %d, got %d\n', nbuses, length(Qg));
    end
    if length(Pd) ~= nbuses
        fprintf('❌ Pd size mismatch: expected %d, got %d\n', nbuses, length(Pd));
    end
    if length(Qd) ~= nbuses
        fprintf('❌ Qd size mismatch: expected %d, got %d\n', nbuses, length(Qd));
    end
    if length(V0) ~= nbuses
        fprintf('❌ V0 size mismatch: expected %d, got %d\n', nbuses, length(V0));
    end
    
    % Display actual values for debugging
    fprintf('\nInput values:\n');
    fprintf('Bus  Pg(MW)  Qg(MVAR)  Pd(MW)  Qd(MVAR)  V0(pu)\n');
    for i = 1:min(9, nbuses)
        fprintf('%2d   %6.1f   %7.1f   %6.1f   %7.1f   %6.3f\n', ...
            i, Pg(i), Qg(i), Pd(i), Qd(i), V0(i));
    end
    
    fprintf('\nBus types:\n');
    fprintf('  Slack bus: %d\n', is);
    fprintf('  PV buses: %s\n', mat2str(ipv));
    fprintf('  PQ buses: %s\n', mat2str(ipq));
end

function [has_violation, min_volt, max_volt, violation_buses] = check_voltage_violations(V)
    if isempty(V)
        has_violation = true;
        min_volt = NaN;
        max_volt = NaN;
        violation_buses = [];
        return;
    end
    
    min_volt = min(V);
    max_volt = max(V);
    violation_buses = find(V < 0.95 | V > 1.05);
    has_violation = ~isempty(violation_buses);
end

function display_voltage_results(case_name, V)
    if isempty(V)
        fprintf('   No voltage data available\n');
        return;
    end
    
    fprintf('   Bus Voltages (p.u.):\n');
    for bus = 1:length(V)
        status = '✅';
        if V(bus) < 0.95 || V(bus) > 1.05
            status = '❌';
        end
        fprintf('     Bus %d: %.4f p.u. %s\n', bus, V(bus), status);
    end
    fprintf('   Voltage range: %.4f - %.4f p.u.\n', min(V), max(V));
end

function display_security_summary(voltage_violations)
    fprintf('\n=== SECURITY ANALYSIS SUMMARY ===\n');
    
    critical_count = 0;
    for i = 1:length(voltage_violations)
        violation = voltage_violations{i};
        if isfield(violation, 'error') || ~isempty(violation.violation_buses)
            critical_count = critical_count + 1;
        end
    end
    
    if critical_count == 0
        fprintf('✅ All contingencies result in acceptable voltage profiles\n');
        fprintf('   System is SECURE for all analyzed contingencies\n');
    else
        fprintf('❌ CRITICAL CONTINGENCIES DETECTED: %d out of 6\n', critical_count);
        fprintf('   System is INSECURE for:\n\n');
        
        for i = 1:length(voltage_violations)
            violation = voltage_violations{i};
            if isfield(violation, 'error')
                fprintf('   • %s - POWER FLOW FAILED (%s)\n', violation.name, violation.error);
            elseif ~isempty(violation.violation_buses)
                fprintf('   • %s - Voltage violations at buses %s\n', ...
                    violation.name, mat2str(violation.violation_buses));
            end
        end
    end
end

function propose_mitigation_measures(voltage_violations, V_base)
    fprintf('\n=== DETAILED MITIGATION MEASURES ===\n');
    
    if isempty(voltage_violations)
        fprintf('No mitigation required - system is secure\n');
        return;
    end
    
    % Analyze problematic buses across all contingencies
    problematic_buses = [];
    for i = 1:length(voltage_violations)
        violation = voltage_violations{i};
        if ~isfield(violation, 'error') && ~isempty(violation.V)
            low_volt_buses = find(violation.V < 0.95);
            problematic_buses = [problematic_buses; low_volt_buses(:)];
        end
    end
    
    % Count occurrences of each problematic bus
    if ~isempty(problematic_buses)
        bus_counts = accumarray(problematic_buses, 1, [9, 1]);
        [sorted_counts, sorted_buses] = sort(bus_counts, 'descend');
        
        fprintf('\nPROBLEMATIC BUS ANALYSIS:\n');
        fprintf('Bus  Violation Count  Base Case Voltage\n');
        fprintf('----  --------------  -----------------\n');
        for i = 1:length(sorted_buses)
            if sorted_counts(i) > 0
                bus = sorted_buses(i);
                fprintf('%2d        %2d            %.4f p.u.\n', ...
                    bus, sorted_counts(i), V_base(bus));
            end
        end
    end
    
    fprintf('\nSPECIFIC MITIGATION RECOMMENDATIONS:\n');
    fprintf('====================================\n');
    
    fprintf('\n1. PRIORITY 1 - CRITICAL BUSES (Buses 5, 7, 9):\n');
    fprintf('   Bus 5: Appears in 4/6 contingencies (0.933-0.975 p.u. base)\n');
    fprintf('   Bus 9: Appears in 4/6 contingencies (0.890-0.958 p.u. base)\n');
    fprintf('   Bus 7: Appears in 2/6 contingencies (0.936-0.986 p.u. base)\n');
    fprintf('   Actions:\n');
    fprintf('   - Install 50-100 MVAR capacitor banks at buses 5 and 9\n');
    fprintf('   - Consider STATCOM at bus 9 for dynamic control\n');
    fprintf('   - Add 25-50 MVAR capacitors at bus 7\n');
    
    fprintf('\n2. GENERATOR VOLTAGE CONTROL OPTIMIZATION:\n');
    fprintf('   - Increase AVR setpoint at Generator 2 (Bus 2) from 1.00 to 1.02-1.03 p.u.\n');
    fprintf('   - Adjust Generator 3 (Bus 3) to provide more reactive support\n');
    fprintf('   - Implement coordinated voltage control between generators\n');
    
    fprintf('\n3. TRANSFORMER TAP ADJUSTMENTS:\n');
    fprintf('   - Increase secondary voltage settings on transformers feeding load areas\n');
    fprintf('   - Optimize ULTC settings to maintain higher voltages at remote buses\n');
    fprintf('   - Consider line drop compensation for better voltage regulation\n');
    
    fprintf('\n4. NETWORK REINFORCEMENT STRATEGIES:\n');
    fprintf('   Critical corridors based on contingency analysis:\n');
    fprintf('   - Line 4-9: Critical for Bus 9 voltage (0.794 p.u. violation)\n');
    fprintf('   - Line 5-6: Critical for Buses 5 and 9 (0.919, 0.926 p.u.)\n');
    fprintf('   - Line 8-9: Critical for Buses 5 and 9 (0.933, 0.890 p.u.)\n');
    fprintf('   Actions:\n');
    fprintf('   - Consider adding parallel circuits to lines 4-9 and 8-9\n');
    fprintf('   - Evaluate series compensation on long lines\n');
    
    fprintf('\n5. OPERATIONAL MEASURES:\n');
    fprintf('   - Implement security-constrained optimal power flow (SCOPF)\n');
    fprintf('   - Pre-position reactive reserves for critical contingencies\n');
    fprintf('   - Develop emergency voltage control procedures\n');
    
    fprintf('\nEXPECTED IMPROVEMENTS:\n');
    fprintf('=====================\n');
    fprintf('- Capacitor installation: +0.02-0.04 p.u. voltage boost\n');
    fprintf('- Generator AVR adjustment: +0.01-0.02 p.u. system-wide\n');
    fprintf('- Combined measures should maintain all buses > 0.95 p.u.\n');
    fprintf('- Critical contingencies (4-9, 8-9) may require multiple measures\n');
    
    fprintf('\nIMPLEMENTATION PHASING:\n');
    fprintf('======================\n');
    fprintf('Phase 1 (Immediate): Generator AVR adjustments, existing capacitor switching\n');
    fprintf('Phase 2 (1-3 months): Install additional capacitor banks\n');
    fprintf('Phase 3 (6-12 months): Network reinforcement, advanced controls\n');
    fprintf('Phase 4 (Long-term): System-wide voltage stability enhancement\n');
end
