%% Team 3 – Steady-State Calculation with fsolve
% (1) CA*, T*, Tc*, Tline* 계산 + 에너지 밸런스 확인

clear; clc;

% -----------------------------
% [1] Parameters (Team #3)
% -----------------------------
V      = 27.0;          % L, reactor volume
Vc     = 6.0;           % L, jacket volume
q      = 8.0;           % L/min, feed flow rate
qc0    = 18.0;          % L/min, coolant flow rate (nominal SS)
Caf0   = 1.3;           % mol/L, feed concentration of A
Tf     = 325.0;         % K, feed temperature
Tcin   = 315.0;         % K, coolant inlet temperature
Vline  = 9.0;           % L, line holdup volume
Tline0 = 355.0;         % K, initial line temperature (동특성 초기값용)
k0     = 9.0e10;        % 1/min, pre-exponential factor
EoverR = 9250.0;        % K, activation term E/R
dH     = -70000;        % J/mol, heat of reaction (exothermic)
UA     = 9000000.0;       % J/(min·K), heat transfer coefficient × area
rhoCp  = 18000.0;       % J/(L·K), reacting mixture heat capacity
rhoCpc = 25000.0;       % J/(L·K), coolant heat capacity

% Steady-state에서 사용할 냉각수 유량
qc = qc0;               % L/min

% -----------------------------
% [2] fsolve 설정
%     x = [CA; T; Tc; Tline]
% -----------------------------
x0 = [0.05; 380.0; 330.0; 320.0];   % initial guess

fun = @(x) cstr_ss_eq(x, V, Vc, q, qc, Caf0, Tf, Tcin, ...
                      Vline, k0, EoverR, dH, UA, rhoCp, rhoCpc);

options = optimoptions('fsolve', ...
    'Display', 'iter', ...          % 반복 과정 보고
    'TolFun', 1e-12, ...
    'TolX',   1e-12);

[x_ss, fval, exitflag, output] = fsolve(fun, x0, options);

CA_ss    = x_ss(1);
T_ss     = x_ss(2);
Tc_ss    = x_ss(3);
Tline_ss = x_ss(4);

fprintf('\n=== Steady State (Team 3, fsolve) ===\n');
fprintf('CA*    = %.6f mol/L\n', CA_ss);
fprintf('T*     = %.3f K\n',     T_ss);
fprintf('Tc*    = %.3f K\n',     Tc_ss);
fprintf('Tline* = %.3f K\n\n',   Tline_ss);

% -----------------------------
% [3] 에너지 밸런스(reactor) 확인
% -----------------------------
k_ss  = k0 * exp(-EoverR / T_ss);   % 1/min
rA_ss = k_ss * CA_ss;              % mol/(L·min)

Q_flow    = q  * rhoCp  * (Tf - T_ss);       % J/min
Q_rxn     = (-dH) * rA_ss * V;              % J/min
Q_cooling = - UA * (T_ss - Tc_ss);          % J/min

fprintf('Energy balance in reactor (J/min):\n');
fprintf(' Q_flow    = %+12.4e\n', Q_flow);
fprintf(' Q_rxn     = %+12.4e\n', Q_rxn);
fprintf(' Q_cooling = %+12.4e\n', Q_cooling);
fprintf(' Sum       = %+12.4e\n', Q_flow + Q_rxn + Q_cooling);

% =========================================================
%  Local function: steady-state equations
% =========================================================
function F = cstr_ss_eq(x, V, Vc, q, qc, Caf0, Tf, Tcin, ...
                         Vline, k0, EoverR, dH, UA, rhoCp, rhoCpc)
    % x = [CA; T; Tc; Tline]
    CA    = x(1);   % mol/L
    T     = x(2);   % K
    Tc    = x(3);   % K
    Tline = x(4);   % K

    % 반응 속도
    k  = k0 * exp(-EoverR / T);   % 1/min
    rA = k * CA;                  % mol/(L·min)

    % (1) A에 대한 물질수지
    eq1 = q/V * (Caf0 - CA) - rA;

    % (2) 반응기 에너지 수지
    eq2 = q/V * (Tf - T) ...
        + (-dH / rhoCp) * rA ...
        - UA/(rhoCp * V) * (T - Tc);

    % (3) 재킷 에너지 수지
    eq3 = qc/Vc * (Tline - Tc) ...
        + UA/(rhoCpc * Vc) * (T - Tc);

    % (4) 라인(홀드업 탱크) 에너지 수지
    eq4 = qc/Vline * (Tcin - Tline);

    F = [eq1; eq2; eq3; eq4];
end
