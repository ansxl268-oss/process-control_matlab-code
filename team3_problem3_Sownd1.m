% Team 3 – Dynamic simulation + SOWND identification using qc = 18 L/min
clear; clc; close all;

%% [1] Parameters (Team #3)
V      = 27.0;          % L, reactor volume
Vc     = 6.0;           % L, jacket volume
q      = 8.0;           % L/min, feed flow rate
qc0    = 12.0;          % L/min, nominal coolant flow rate (reference)
Caf0   = 1.3;           % mol/L, feed concentration of A
Tf     = 325.0;         % K, feed temperature
Tcin   = 315.0;         % K, coolant inlet temperature
Vline  = 9.0;           % L, line holdup volume
Tline0 = 355.0;         % K, initial line temperature
k0     = 9.0e10;        % 1/min, pre-exponential factor
EoverR = 9250.0;        % K, activation term E/R
dH     = -70000;        % J/mol, heat of reaction (exothermic)
UA     = 9000000.0;     % J/(min·K), heat transfer coefficient × area
rhoCp  = 18000.0;       % J/(L·K), heat capacity of reacting mixture
rhoCpc = 25000.0;       % J/(L·K), heat capacity of coolant

%% [2] Initial conditions (common for all qc)
CA0 = Caf0;
T0  = Tf;
Tc0 = Tcin;
x0  = [CA0; T0; Tc0; Tline0];

%% [3] Parameter struct
p.V      = V;
p.Vc     = Vc;
p.q      = q;
p.qc0    = qc0;
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

%% [4] Simulation settings
tspan = [0 60];               % min
qc_6   = 6.0;
qc_12  = 12.0;
qc_18  = 18.0;
opts   = odeset('RelTol',1e-6, 'AbsTol',1e-8);

%% [5] Case qc = 6 L/min
p6 = p; p6.qc = qc_6;
[t6, x6] = ode15s(@(t,x) cstr_team3_odes(t,x,p6), tspan, x0, opts);
CA_6    = x6(:,1);
T_6     = x6(:,2);
Tc_6    = x6(:,3);
Tline_6 = x6(:,4);

%% [6] Case qc = 12 L/min
p12 = p; p12.qc = qc_12;
[t12, x12] = ode15s(@(t,x) cstr_team3_odes(t,x,p12), tspan, x0, opts);
CA_12    = x12(:,1);
T_12     = x12(:,2);
Tc_12    = x12(:,3);
Tline_12 = x12(:,4);

%% [7] Case qc = 18 L/min  (이 케이스로 SOWND 식별)
p18 = p; p18.qc = qc_18;
[t18, x18] = ode15s(@(t,x) cstr_team3_odes(t,x,p18), tspan, x0, opts);
CA_18    = x18(:,1);
T_18     = x18(:,2);
Tc_18    = x18(:,3);
Tline_18 = x18(:,4);

%% [8] Plot – original approach-to-SS curves (확인용)
figure;

% (a) Reactor temperature
subplot(3,1,1);
plot(t6,  T_6,  'r--', 'LineWidth',1.5); hold on;
plot(t12, T_12, 'k-',  'LineWidth',1.5);
plot(t18, T_18, 'b-',  'LineWidth',1.5);
xlabel('Time (min)');
ylabel('T (K)');
title('Reactor temperature T(t): approach to steady state');
legend('qc = 6','qc = 12','qc = 18','Location','best');

% (b) Jacket temperature
subplot(3,1,2);
plot(t6,  Tc_6,  'r--', 'LineWidth',1.5); hold on;
plot(t12, Tc_12, 'k-',  'LineWidth',1.5);
plot(t18, Tc_18, 'b-',  'LineWidth',1.5);
xlabel('Time (min)');
ylabel('T_c (K)');
title('Coolant (jacket) temperature T_c(t)');
legend('qc = 6','qc = 12','qc = 18','Location','best');

% (c) Line temperature
subplot(3,1,3);
plot(t6,  Tline_6,  'r--', 'LineWidth',1.5); hold on;
plot(t12, Tline_12, 'k-',  'LineWidth',1.5);
plot(t18, Tline_18, 'b-',  'LineWidth',1.5);
yline(Tcin, 'k:', 'LineWidth',1.0);
xlabel('Time (min)');
ylabel('T_{line} (K)');
title('Line temperature T_{line}(t)');
legend('qc = 6','qc = 12','qc = 18','T_{cin}=315K','Location','best');

sgtitle('Team 3: Approach to steady state from initial condition for different q_c');

%% [9] SOWND model identification using qc = 18 trajectory

% ---- 9-1) ode15s 결과(비등간격)를 등간격 시간축으로 보간 ----
t_raw = t18;          % 원래 시간 벡터 (비등간격)
T_raw = T_18;         % 원래 온도

T0    = T_raw(1);     % t=0 온도
y_raw = T_raw - T0;   % 편차 출력 y(t) = T(t) - T(0)

N    = length(t_raw);
t_id = linspace(t_raw(1), t_raw(end), N).';          % 등간격 column vector
y_id = interp1(t_raw, y_raw, t_id, 'pchip');         % 보간
y_id = y_id(:);                                      % 확실하게 column으로

du   = qc_18;        % 입력 스텝 크기 (0 -> 18 L/min)

% ---- 9-2) 초기 추정값 ----
K0 = y_id(end) / du;   % 최종값 기반 이득 추정

% theta_log = [K, log(tau1), log(tau2), log(|tauz|), log(theta_d)]
theta0 = [K0, log(1), log(10), log(1), log(0.5)];

% ---- 9-3) fminsearch로 제곱오차 최소화 ----
theta_hat = fminsearch(@(th) sownd_obj(th, t_id, y_id, du), theta0);

% ---- 9-4) 물리 파라미터로 변환 ----
K_hat      = theta_hat(1);
tau1_hat   = exp(theta_hat(2));
tau2_hat   = exp(theta_hat(3));
tauz_hat   = -exp(theta_hat(4));   % 항상 음수 (RHP zero)
theta_dhat = exp(theta_hat(5));    % dead time > 0

% ---- 9-5) 피팅 응답 계산 ----
y_fit = sownd_step(theta_hat, t_id, du);
T_fit = T0 + y_fit;

% ---- 9-6) 결과 플롯 (qc = 18만 비교) ----
figure;
plot(t_raw, T_raw, 'k', 'LineWidth',1.5); hold on;
plot(t_id,  T_fit, 'r--', 'LineWidth',1.5);
xlabel('Time (min)');
ylabel('T (K)');
legend('Nonlinear CSTR (qc = 18)','SOWND fit','Location','best');
title('SOWND identification from qc = 18 trajectory');
grid on;

% ---- 9-7) 파라미터 출력 ----
fprintf('\nSOWND parameters identified from qc = 18 trajectory:\n');
fprintf('  K       = %g   [K / (L/min)]\n', K_hat);
fprintf('  tau1    = %g   [min]\n', tau1_hat);
fprintf('  tau2    = %g   [min]\n', tau2_hat);
fprintf('  tau_z   = %g   [min]  (RHP zero if < 0)\n', tauz_hat);
fprintf('  theta_d = %g   [min]  (dead time)\n\n', theta_dhat);

%% [10] 동일 SOWND 모델로 qc = 6, 12, 18 응답 예측해서 비교

% 비교용 시간축: SOWND 식별에 썼던 등간격 t_id 사용
t_cmp = t_id;

% 비선형 모델 응답을 같은 시간축으로 보간
T6_cmp  = interp1(t6,  T_6,  t_cmp, 'pchip');
T12_cmp = interp1(t12, T_12, t_cmp, 'pchip');
T18_cmp = interp1(t18, T_18, t_cmp, 'pchip');

% SOWND 모델 예측 (동일 파라미터 theta_hat, 스텝 크기만 qc 값으로)
y6_lin  = sownd_step(theta_hat, t_cmp, qc_6);
y12_lin = sownd_step(theta_hat, t_cmp, qc_12);
y18_lin = sownd_step(theta_hat, t_cmp, qc_18);

T6_lin  = T0 + y6_lin;
T12_lin = T0 + y12_lin;
T18_lin = T0 + y18_lin;

% 플롯: qc별로 Nonlinear vs SOWND 비교
figure;

subplot(3,1,1);
plot(t_cmp, T6_cmp,  'k',  'LineWidth', 1.5); hold on;
plot(t_cmp, T6_lin,  'r--','LineWidth', 1.5);
xlabel('Time (min)'); ylabel('T (K)');
title('qc = 6 L/min: Nonlinear CSTR vs SOWND TF');
legend('Nonlinear','SOWND','Location','best'); grid on;

subplot(3,1,2);
plot(t_cmp, T12_cmp, 'k',  'LineWidth', 1.5); hold on;
plot(t_cmp, T12_lin, 'r--','LineWidth', 1.5);
xlabel('Time (min)'); ylabel('T (K)');
title('qc = 12 L/min: Nonlinear CSTR vs SOWND TF');
legend('Nonlinear','SOWND','Location','best'); grid on;

subplot(3,1,3);
plot(t_cmp, T18_cmp, 'k',  'LineWidth', 1.5); hold on;
plot(t_cmp, T18_lin, 'r--','LineWidth', 1.5);
xlabel('Time (min)'); ylabel('T (K)');
title('qc = 18 L/min: Nonlinear CSTR vs SOWND TF');
legend('Nonlinear','SOWND','Location','best'); grid on;

sgtitle('T/q_c SOWND model (identified at qc=18) vs nonlinear CSTR for qc = 6, 12, 18');

%% -------- Local functions ----------------------------------------------

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
    k  = k0 * exp(-EoverR / T);
    rA = k * CA;

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

function y_model = sownd_step(theta_log, t, du)
    % theta_log = [K, log(tau1), log(tau2), log(|tauz|), log(theta_d)]
    K       = theta_log(1);
    tau1    = exp(theta_log(2));
    tau2    = exp(theta_log(3));
    tauz    = -exp(theta_log(4));   % 항상 음수 (RHP zero)
    theta_d = exp(theta_log(5));    % dead time > 0

    num = K * [tauz 1];
    den = conv([tau1 1], [tau2 1]);
    G   = tf(num, den, 'InputDelay', theta_d);

    [y_step, ~] = step(G, t);   % unit step 응답
    y_step  = y_step(:);       % column vector
    y_model = du * y_step;     % 실제 스텝 크기로 스케일
end

function J = sownd_obj(theta_log, t, y_data, du)
    y_model = sownd_step(theta_log, t, du);
    y_model = y_model(:);
    y_data  = y_data(:);

    e = y_model - y_data;
    J = sum(e.^2);             % 스칼라
end
