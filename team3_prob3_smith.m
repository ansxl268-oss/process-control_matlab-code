%% Team 3 – Smith method identification for T/qc using qc = 18 L/min
%  - 비선형 CSTR + 재킷 + 라인 모델(qc=18) 동특성 시뮬레이션
%  - Smith method (t20/t60) 로 2차 + dead time 모델 식별
%  - 비선형 응답 vs Smith 2차 모델 응답 비교

clear; clc; close all;

%% [1] Parameters (Team #3)
V      = 27.0;          % L, reactor volume
Vc     = 6.0;           % L, jacket volume
q      = 8.0;           % L/min, feed flow rate
qc_18  = 18.0;          % L/min, coolant flow rate (이번엔 이 케이스만 사용)
Caf0   = 1.3;           % mol/L, feed concentration of A
Tf     = 325.0;         % K, feed temperature
Tcin   = 315.0;         % K, coolant inlet temperature
Vline  = 9.0;           % L, line holdup volume
Tline0 = 355.0;         % K, initial line temp (hot line)
k0     = 9.0e10;        % 1/min, pre-exponential factor
EoverR = 9250.0;        % K, activation term E/R
dH     = -70000.0;      % J/mol, heat of reaction (exothermic)
UA     = 9.0e6;         % J/(min·K), heat transfer coefficient × area
rhoCp  = 18000.0;       % J/(L·K), heat capacity of reacting mixture
rhoCpc = 25000.0;       % J/(L·K), heat capacity of coolant

%% [2] Initial conditions (from problem statement)
CA0   = Caf0;   % mol/L
T0    = Tf;     % K
Tc0   = Tcin;   % K
% Tline0는 위에서 이미 정의됨 (355 K)

x0 = [CA0; T0; Tc0; Tline0];   % [CA; T; Tc; Tline]

%% [3] Parameter struct (for ODE function)
p.V      = V;
p.Vc     = Vc;
p.q      = q;
p.qc     = qc_18;     % 이번 run에서는 qc=18 고정
p.Caf0   = Caf0;
p.Tf     = Tf;
p.Tcin   = Tcin;
p.Vline  = Vline;
p.k0     = k0;
p.EoverR = EoverR;
p.dH     = dH;
p.UA     = UA;
p.rhoCp  = rhoCp;
p.rhoCpc = rhoCpc;

%% [4] Nonlinear simulation for qc = 18 (approach to steady state)
tspan = [0 60];                      % min
opts  = odeset('RelTol',1e-6,'AbsTol',1e-8);

odefun_18 = @(t,x) cstr_team3_odes(t, x, p);
[t18, x18] = ode15s(odefun_18, tspan, x0, opts);

CA_18    = x18(:,1);
T_18     = x18(:,2);   % 우리가 식별에 사용할 출력 (reactor temperature)
Tc_18    = x18(:,3);
Tline_18 = x18(:,4);

%% [5] Smith method – 1) Gain K 계산
T_init  = T_18(1);
T_final = T_18(end);
dT      = T_final - T_init;      % 출력 변화량 ΔT

du      = qc_18;                 % qc step 크기 (0 -> 18 L/min 이라 가정)
K       = dT / du;               % 시간영역 이득 [K / (L/min)]

fprintf('Step gain K (from qc to T) = %.5f [K / (L/min)]\n', K);

%% [6] Smith method – 2) 정규화 응답 & dead time(theta) 추정

% 정규화 응답 (0 → 1)
y_norm = (T_18 - T_init) / dT;

% dead time 추정:
%  - y_norm이 아주 작다가 일정 threshold(예: 0.02) 이상이 되는 첫 시점을 θ로 잡음
threshold = 0.02;
idx_theta = find(y_norm > threshold, 1, 'first');

if isempty(idx_theta)
    theta = 0.0;   % 만약 응답이 너무 느리면 fallback
else
    theta = t18(idx_theta);
end

fprintf('Estimated dead time theta ≈ %.6f min\n', theta);

% θ 이후의 시간축으로 shift (Smith method에서 "delay 제거" 개념)
t_shift = t18 - theta;
mask    = t_shift >= 0;      % delay 이후 구간만 사용

t_s = t_shift(mask);         % shifted time (t')
y_s = y_norm(mask);          % 그에 해당하는 정규화 응답

%% [7] Smith method – 3) t20, t60 계산

% y_s(t')가 0.2, 0.6이 되는 시간 찾기
t20 = interp1(y_s, t_s, 0.2);
t60 = interp1(y_s, t_s, 0.6);

fprintf('t20 (y=0.2) = %.3f min (after removing delay)\n', t20);
fprintf('t60 (y=0.6) = %.3f min (after removing delay)\n', t60);

R = t20 / t60;
fprintf('R = t20/t60 = %.4f  --> 이 값을 Smith chart(Fig.7.7)에 사용\n', R);

%% [8] Smith chart에서 읽은 zeta, tau 반영 (네가 이미 읽은 값 사용)

% ======= 네가 Smith chart에서 읽은 값 =======
zeta = 2.3;   % Identified zeta (from Smith chart)
tau  = 0.28;   % Identified tau  (from Smith chart)
% ==========================================

fprintf('Identified zeta (from Smith chart) = %.3f\n', zeta);
fprintf('Identified tau  (from Smith chart) = %.3f min\n', tau);

%% [9] Smith 2nd-order + dead time 모델 구성

% G(s) = K * e^(-theta s) / (tau^2 s^2 + 2 zeta tau s + 1)
num = K;
den = [tau^2, 2*zeta*tau, 1];

G_smith = tf(num, den, 'InputDelay', theta);

disp('Smith-identified 2nd-order + dead-time model G(s) from qc to T:');
G_smith

%% [10] Nonlinear vs Smith model 응답 비교 (step 응답)

% === 1) Smith 모델 step 응답: "등간격" 시간축으로 계산 ===
t_lin = linspace(0, max(t18), 2000);               % 등간격 시간축
[y_lin_step, t_lin] = step(G_smith, t_lin);        % unit step 응답

% === 2) 비선형 t18 시간축에 맞게 보간 ===
y_lin_interp = interp1(t_lin, y_lin_step, t18, 'pchip');

% === 3) 실제 입력 step 크기 du 반영해서 T_smith 계산 ===
T_smith = T_init + du * y_lin_interp;              % T(t) 예측

figure;

% (a) T(t) 절대값 비교
subplot(2,1,1);
plot(t18, T_18,     'k',  'LineWidth',1.5); hold on;
plot(t18, T_smith,  'r--','LineWidth',1.5);
xlabel('Time (min)');
ylabel('T (K)');
title('Nonlinear CSTR vs Smith 2nd-order model (absolute T)');
legend('Nonlinear (qc=18)', 'Smith model', 'Location','best');
grid on;

%% [11] Normalized 응답 비교 (delay 제거 후)

% --- 비선형: 이미 t_s, y_s (shifted & normalized) 존재 ---
% --- Smith 모델: 동일한 정상상태 변화량 사용해서 normalized 만들기 ---

% Smith step 응답 전체를 normalized로 변환
y_lin_norm_all = (du * y_lin_step) / dT;          % (T - T_init)/dT 와 비교할 수 있도록

% 시간축도 delay 제거: t_lin_shift = t - theta
t_lin_shift = t_lin - theta;
mask_lin    = t_lin_shift >= 0;

t_ls  = t_lin_shift(mask_lin);      % shifted time for Smith model
y_ls  = y_lin_norm_all(mask_lin);   % 그에 대응하는 normalized 응답

% y_ls를 비선형 t_s 시간축으로 보간
y_lin_on_ts = interp1(t_ls, y_ls, t_s, 'pchip');

subplot(2,1,2);
plot(t_s, y_s,        'k',  'LineWidth',1.5); hold on;
plot(t_s, y_lin_on_ts,'r--','LineWidth',1.5);
yline(0.2,'k:'); yline(0.6,'k:');
xline(t20,'b:','t_{20}'); xline(t60,'b:','t_{60}');
xlabel('Shifted time t - \theta (min)');
ylabel('Normalized T');
title('Normalized response after removing delay');
legend('Nonlinear (normalized)','Smith model (normalized)', ...
       'Location','best');
grid on;

%% [12] t20, t60 위치를 보여주는 별도 그래프 (비선형만 보기용)

figure;
plot(t_s, y_s, 'k','LineWidth',1.5); hold on;
yline(0.2,'r--'); yline(0.6,'b--');
xline(t20,'r:','t_{20}'); xline(t60,'b:','t_{60}');
xlabel('Shifted time t - \theta (min)');
ylabel('Normalized T');
title('Smith method: t_{20} and t_{60} on normalized response (nonlinear)');
grid on;

%% ===== Local function: CSTR + jacket + line ODE =====
function dxdt = cstr_team3_odes(~, x, p)
    % x = [CA; T; Tc; Tline]

    CA    = x(1);   % mol/L
    T     = x(2);   % K
    Tc    = x(3);   % K
    Tline = x(4);   % K

    % parameters
    V      = p.V;
    Vc     = p.Vc;
    q      = p.q;
    qc     = p.qc;
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

    % reaction rate
    k  = k0 * exp(-EoverR / T);   % 1/min
    rA = k * CA;                  % mol/(L·min), consumption of A

    % ODEs
    dCA    = q/V * (Caf0 - CA) - rA;

    dT     = q/V * (Tf - T) ...
             + (-dH / rhoCp) * rA ...
             - UA/(rhoCp * V) * (T - Tc);

    dTc    = qc/Vc * (Tline - Tc) ...
             + UA/(rhoCpc * Vc) * (T - Tc);

    dTline = qc/Vline * (Tcin - Tline);

    dxdt = [dCA; dT; dTc; dTline];
end
