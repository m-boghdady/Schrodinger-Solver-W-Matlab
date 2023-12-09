
function schrodingerSolver()
    % Constants
    hbar = 1.0545718e-34; % Reduced Planck's constant (J s)
    e = 1.60217662e-19;   % Elementary charge (C)
    me = 9.10938356e-30 ;  % Electron mass (kg)

    % User inputs in standard units
    disp('Select the potential type:');
    disp('1. Infinite Square Well');
    disp('2. Finite Square Well');
    disp('3. Harmonic Oscillator');
    V_type = input('Enter your choice (1/2/3): ');
    L_nm = input('Enter the width of the well in nm: '); % Width in nanometers
    L = L_nm * 1e-9; % Convert nm to m for calculations
    m_eff = input('Enter the effective mass of the particle as a ratio to me (electron mass): '); % Effective mass as a ratio of electron mass
    m = m_eff * me; % Effective mass in kg
    N = input('Enter the number of grid points: ');
    numEigenvaluesToPlot = input('Enter the number of eigenenergies to plot: ');

    % Define spatial grid
    x = linspace(-L/2, L/2, N)'; % x goes from -L to L to center the harmonic oscillator
    dx = x(2) - x(1); % Grid spacing

    % Initialize the potential V for each type of well
    V = zeros(N, 1);
    V_inf = 1e9 * e;

    %theoritical energies from the analytic formula
    theoreticalEnergies = zeros(numEigenvaluesToPlot, 1);


    if V_type == 1 % Infinite Square Well
        V0_eV=V_inf;
        V(:) = 0; 
        V(1) = V0_eV * e; % Set the first and last point to the barrier height
        V(end) = V0_eV * e;
    elseif V_type == 2 % Finite Square Well
        V0_eV = input('Enter the height of the finite well in eV: ');  % Height in electron volts
        V(1) = V0_eV * e;
        V(end) = V0_eV * e;
    elseif V_type == 3 % Harmonic Oscillator
        k_eV_per_nm2 = input('Enter the spring constant for the harmonic oscillator in eV/nm^2: '); % Spring constant in eV/nm^2
        k = k_eV_per_nm2 * e * (1e9)^2; % Convert eV/nm^2 to J/m^2 for calculations
        V = 0.5 * k * (x).^2; % Parabolic potential centered at L/2
    end

    % Construct the Hamiltonian matrix H using finite differences
    t0 = hbar^2 / (2 * m * dx^2); % Kinetic energy coefficient
    main_diag = 2 * t0 + V; % Main diagonal with kinetic and potential energy
    off_diag = -t0 * ones(N-1, 1); % Off-diagonal with kinetic energy
    H = diag(main_diag) + diag(off_diag, 1) + diag(off_diag, -1);

    % Solve the eigenvalue problem to find energy levels and wavefunctions
    [psi, E] = eig(H);
    energies_eV = diag(E) / e; % Convert the energies from J to eV for output

    % Calculate theoretical energies for Infinite Square Well and Harmonic Oscillator
    if V_type == 1 % Infinite Square Well
        for n = 1:numEigenvaluesToPlot
            theoreticalEnergies(n) = (n^2 * pi^2 * hbar^2) / (2 * m * L^2) / e; % Convert to eV
        end
    elseif V_type == 3 % Harmonic Oscillator
        omega = sqrt(k / m); % Angular frequency
        for n = 1:numEigenvaluesToPlot
            theoreticalEnergies(n) = (n - 0.5) * hbar * omega / e; % Convert to eV
        end
    end

    % Plot the eigenenergies
    figure;
    hold on;
    title('Eigenenergies');
    xlabel('Eigenvalue number');
    ylabel('Energy (eV)');
    
    % Plot eigenenergies as points
    for i = 1:numEigenvaluesToPlot
        plot(i, energies_eV(i), 'bo', 'MarkerFaceColor', 'b');
    end
    hold off;

    % Display calculated and theoretical energies
    if V_type == 1 || V_type == 3 || V_type == 2
        disp('Quantum Well Specifications:');
        fprintf('Well Type: %s\n', getWellType(V_type));
        fprintf('Width of the Well (nm): %.2f\n', L_nm);
        if V_type == 2 % Finite Square Well
            fprintf('Height of the Well (eV): %.2f\n', V0_eV);
        elseif V_type == 3 % Harmonic Oscillator
            fprintf('Spring Constant (eV/nm^2): %.2f\n', k_eV_per_nm2);
        end
        fprintf('Effective Mass (ratio to electron mass): %.2f\n', m_eff);
        fprintf('Number of Grid Points: %d\n', N);

        disp('Calculated and Theoretical Energies (eV):');
        for i = 1:numEigenvaluesToPlot
            fprintf('n = %d: Calculated = %.4f eV, Theoretical = %.4f eV\n', ...
                    i, energies_eV(i), theoreticalEnergies(i));
        end
    end

    % Helper function to get the well type as a string
    function wellType = getWellType(typeNum)
        switch typeNum
            case 1
                wellType = 'Infinite Square Well';
            case 2
                wellType = 'Finite Square Well';
            case 3
                wellType = 'Harmonic Oscillator';
            otherwise
                wellType = 'Unknown';
        end
    end


    % Plot the potential V(x) and eigenfunctions psi(x)
    figure;
    hold on;
    title('Potential and Eigenfunctions');
    xlabel('Position (nm)');
    ylabel('Energy (eV)');
    
    % max_energy = max(energies_eV(1:numEigenvaluesToPlot));
       
    % Plot the potential for the harmonic oscillator
    if V_type == 3
        psi = abs(psi);
        plot(x * 1e9, V / e, 'k-', 'LineWidth', 2, 'DisplayName', 'Potential');
    end

    % Plot eigenfunctions offset by their eigenenergies
    for i = 1:numEigenvaluesToPlot
        % Offset each eigenfunction by its eigenenergy
        offset = energies_eV(i);
        plot(x * 1e9, psi(:,i) + offset, 'DisplayName', ['Eigenfunction ' num2str(i)]);
    end

    legend show;
    hold off;
end
