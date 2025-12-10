%% prac2 / prac3 – Nonlinear Step Simulations (q0 = 8 L/min 기준)
%  이제는 qc는 12로 고정, q를 ±50% step (8 -> 12, 8 -> 4)

clear; clc;

%% [1] Parameters (Team #3 그대로 사용)
V      = 27.0;          % L, reactor volume
Vc     = 6.0;           % L, jacket volume
q0     = 8.0;           % L/min, nominal feed flow rate (기준 q)
qc0    = 12.0;          % L/min, nominal coolant flow rate (항상 이 값으로 고정)
Caf0   = 1.3;           % mol/L, feed concentration of A
Tf     = 325.0;         % K, feed temperature
Tcin   = 315.0;         % K, coolant inlet temperature
Vline  = 9.0;           % L, line holdup volume
Tline0 = 355.0;         % K, initial line temp
k0     = 9.0e10;        % 1/min, pre-exponential factor
EoverR = 9250.0;        % K, activation term E/R
dH     = -1.5e6;        % J/mol, heat of reaction (exothermic)
UA     = 55000.0;       % J/(min·K), heat transfer coefficient × area
rhoCp  = 18000.0;       % J/(L·K), heat capacity of reacting mixture
rhoCpc = 25000.0;       % J/(L·K), heat capacity of coolant

%% [2] (예전처럼) q = q0, qc = qc0일 때의 정상상태를 초기조건으로 사용
CA_ss    = 0.041097;    % mol/L
T_ss     = 401.870;     % K
Tc_ss    = 328.459;     % K
Tline_ss = 315.000;     % K

x0_step = [CA_ss; T_ss; Tc_ss; Tline_ss];

%% [3] 공통 파라미터 struct
p.V      = V;
p.Vc     = Vc;
p.q0     = q0;      % 참고용
p.q      = q0;      % 이 run에서 쓸 q (나중에 up/down에서 바뀜)
p.qc     = qc0;     % qc는 모든 run에서 12로 고정
p.Caf0   = Caf0;
p.Tf     = Tf;
p.Tcin   = Tcin;
p.Vline  = Vline;
p.Tline0 = Tline0;
p.k0     = k0;
p.EoverR = EoverR;
p.dH     = dH;
p.UA     = UA;
p.rhoCp  = rhoCp;
p.rhoCpc = rhoCpc;

%% [4] 시뮬레이션 설정
tspan = [0 60];          % min

q_up   = 1.5 * q0;       % 8 -> 12 L/min  ( +50% )
q_down = 0.5 * q0;       % 8 ->  4 L/min  ( -50% )

opts = odeset('RelTol',1e-6, 'AbsTol',1e-8);

%% [5] Up step 시뮬레이션 (q : 8 -> 12, qc는 12로 고정)
p_up    = p;
p_up.q  = q_up;          % 이번 run에서의 feed flow rate
% qc는 p_up.qc = 12.0 그대로

odefun_up = @(t,x) cstr_team3_odes(t, x, p_up);
[t_up, x_up] = ode15s(odefun_up, tspan, x0_step, opts);

CA_up    = x_up(:,1);
T_up     = x_up(:,2);
Tc_up    = x_up(:,3);
Tline_up = x_up(:,4);

%% [6] Down step 시뮬레이션 (q : 8 -> 4, qc는 12로 고정)
p_dn    = p;
p_dn.q  = q_down;

odefun_dn = @(t,x) cstr_team3_odes(t, x, p_dn);
[t_dn, x_dn] = ode15s(odefun_dn, tspan, x0_step, opts);

CA_dn    = x_dn(:,1);
T_dn     = x_dn(:,2);
Tc_dn    = x_dn(:,3);
Tline_dn = x_dn(:,4);

%% [7] 그래프 플로팅

figure;

% (a) Reactor temperature T(t)
subplot(3,1,1);
plot(t_up, T_up, 'b-', 'LineWidth',1.5); hold on;
plot(t_dn, T_dn, 'r--', 'LineWidth',1.5);
yline(T_ss, 'k:');
xlabel('Time (min)');
ylabel('T (K)');
title('Reactor temperature T(t)');
legend('Up step: q 8 \rightarrow 12 (qc=12 고정)', ...
       'Down step: q 8 \rightarrow 4 (qc=12 고정)', ...
       'T^* at q = 8, qc = 12', ...
       'Location','best');

% (b) Jacket temperature Tc(t)
subplot(3,1,2);
plot(t_up, Tc_up, 'b-', 'LineWidth',1.5); hold on;
plot(t_dn, Tc_dn, 'r--', 'LineWidth',1.5);
yline(Tc_ss, 'k:');
xlabel('Time (min)');
ylabel('T_c (K)');
title('Coolant (jacket) temperature T_c(t)');
legend('Up step', 'Down step', 'T_c^* at q = 8, qc = 12', 'Location','best');

% (c) Line temperature Tline(t)
subplot(3,1,3);
plot(t_up, Tline_up, 'b-', 'LineWidth',1.5); hold on;
plot(t_dn, Tline_dn, 'r--', 'LineWidth',1.5);
yline(Tline_ss, 'k:');
xlabel('Time (min)');
ylabel('T_{line} (K)');
title('Line temperature T_{line}(t)');
legend('Up step', 'Down step', 'T_{line}^* at q = 8, qc = 12', 'Location','best');

sgtitle('Team 3: Nonlinear step responses for \pm50% change in q (qc0 = 12 L/min 고정)');

%% [8] 역응답 구간 강조 – T(t)에서만 표시 (원래 코드 그대로 유지)
dT_up_final = T_up(end) - T_ss;
sign_up     = sign(dT_up_final);
idx_inv_up  = find( (T_up - T_ss) .* sign_up < 0 );

dT_dn_final = T_dn(end) - T_ss;
sign_dn     = sign(dT_dn_final);
idx_inv_dn  = find( (T_dn - T_ss) .* sign_dn < 0 );

subplot(3,1,1);
hold on;
plot(t_up(idx_inv_up), T_up(idx_inv_up), 'ko', 'MarkerSize',4, 'MarkerFaceColor','c');
plot(t_dn(idx_inv_dn), T_dn(idx_inv_dn), 'ko', 'MarkerSize',4, 'MarkerFaceColor','m');
text(0.5, T_ss, '  inverse-response region', 'VerticalAlignment','bottom');

%% ==== CSTR ODE 함수 (qc는 struct 안에서 항상 12로 유지) ====
function dxdt = cstr_team3_odes(t, x, p)
% x = [CA; T; Tc; Tline]

CA    = x(1);   % mol/L
T     = x(2);   % K
Tc    = x(3);   % K
Tline = x(4);   % K

% ----- parameters -----
V      = p.V;
Vc     = p.Vc;
q      = p.q;      % 이번 run에서의 feed flow rate (step 후 값)
qc     = p.qc;     % coolant flow rate (모든 run에서 12로 고정)
Caf0   = p.Caf0;
Tf     = p.Tf;
Tcin   = p.Tcin;
Vline  = p.Vline;
k0     = p.k0;
EoverR = p.EoverR;
dH     = p.dH;
UA     = p.UA;
rhoCp  = p.rhoCp;
rhoCpc = p.rhoCpc;

% ----- reaction rate -----
k  = k0 * exp(-EoverR / T);   % 1/min
rA = k * CA;                  % mol/(L·min)

% ----- ODEs -----
dCA = q/V * (Caf0 - CA) - rA;

dT  = q/V * (Tf - T) ...
      + (-dH / rhoCp) * rA ...
      - UA/(rhoCp * V) * (T - Tc);

dTc = qc/Vc * (Tline - Tc) ...
      + UA/(rhoCpc * Vc) * (T - Tc);

dTline = qc/Vline * (Tcin - Tline);

dxdt = [dCA; dT; dTc; dTline];
end
