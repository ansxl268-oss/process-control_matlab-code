%% Team 3 – Inverse-response (Zhang et al.) identification for T/qc using qc = 18 L/min
%  - 비선형 CSTR + 재킷 + 라인 모델(qc=18) 동특성 시뮬레이션
%  - 논문(Zhang et al.)에서 사용하는 inverse-response 모델
%       G(s) = K1/(tau1 s + 1) - K2/(tau2 s + 1)
%    을 비선형 step 응답에 nonlinear regression으로 피팅
%  - (옵션) 식별된 모델로 논문식 controller R(s) 골격까지 생성

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
p.qc     = qc_18;     % 이번 run에서는 qc=18 고정 (0 -> 18 step으로 해석)
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

% --- 입력 step 크기: qc 0 -> 18 L/min 으로 해석 ---
du      = qc_18;          % L/min
T_init  = T_18(1);
T_final = T_18(end);
dT      = T_final - T_init;

fprintf('Observed ΔT = %.4f K for qc step 0 -> %.1f L/min\n', dT, du);

%% [5] Inverse-response 모델 구조 (논문식)
%     G_inv(s) = K1/(tau1 s + 1) - K2/(tau2 s + 1)
%  - K1, K2, tau1, tau2 를 비선형 최소자승으로 피팅

% (1) 측정된 "편차" 출력 (ΔT) 정의
y_meas = T_18 - T_init;   % ΔT(t), [K]

% (2) inverse-response 모델 step 응답을 계산하는 helper 함수 핸들
cost_fun = @(param) invresp_cost(param, t18, y_meas, du);

% param = [K1, K2, log_tau1, log_tau2] 로 잡아서 tau>0 보장
% 초기 추정:
K_total = dT / du;         % 전체 steady-state gain ≈ K1 - K2
K1_0    = 1.5 * K_total;   % 대충 K1 > K2 > 0 가 되도록
K2_0    = 0.5 * K_total;
tau1_0  = 10;              % 느린 time constant (slower dynamics)
tau2_0  = 1;               % 빠른 time constant (faster dynamics)

param0  = [K1_0, K2_0, log(tau1_0), log(tau2_0)];

% (3) fminsearch로 최소자승 피팅
options_opt = optimset('Display','iter','TolX',1e-8,'TolFun',1e-8);
[param_hat, J_min] = fminsearch(cost_fun, param0, options_opt);

K1_hat   = param_hat(1);
K2_hat   = param_hat(2);
tau1_hat = exp(param_hat(3));
tau2_hat = exp(param_hat(4));

fprintf('\n=== Identified inverse-response model parameters ===\n');
fprintf('K1   = %.5f\n', K1_hat);
fprintf('K2   = %.5f\n', K2_hat);
fprintf('tau1 = %.5f min\n', tau1_hat);
fprintf('tau2 = %.5f min\n', tau2_hat);
fprintf('Sum of squared error J = %.6e\n', J_min);

% (4) 식별된 inverse-response 모델 전달함수 구성
s = tf('s');
G1_hat = K1_hat / (tau1_hat*s + 1);
G2_hat = K2_hat / (tau2_hat*s + 1);
G_inv_hat = G1_hat - G2_hat;    % 논문 구조의 inverse-response 모델

% zero 위치 (inverse-response 여부 확인)
[num_G, den_G] = tfdata(G_inv_hat, 'v');
% G_inv_hat(s) = [(K1*tau2 - K2*tau1) s + (K1 - K2)] / [(tau1 s+1)(tau2 s+1)]
zero_hat = roots(num_G);
fprintf('\nEstimated zeros of G_inv_hat(s):\n');
disp(zero_hat);

%% [6] 비선형 응답 vs inverse-response 모델 응답 비교

% 1) 식별된 G_inv_hat 에 qc step (크기 du=18) 적용
%    → step()에는 균일한 시간축을 사용하고, 나중에 t18에 보간
t_final = t18(end);
N_sim   = 2000;                              % step 시뮬레이션용 샘플 수
t_step  = linspace(0, t_final, N_sim);       % 균일한 시간축

[y_inv_step, t_step] = step(G_inv_hat * du, t_step);   % 여기서 y_inv_step = ΔT_model(t_step)

% 2) 비선형 시뮬레이션 시간축 t18에 맞게 보간
y_inv_interp = interp1(t_step, y_inv_step, t18, 'pchip');   % ΔT_model(t18)

% 3) 절대 온도 T(t)로 복원
T_model = T_init + y_inv_interp;

figure;
subplot(2,1,1);
plot(t18, T_18,   'k',  'LineWidth',1.5); hold on;
plot(t18, T_model,'r--','LineWidth',1.5);
xlabel('Time (min)');
ylabel('T (K)');
title('Nonlinear CSTR vs inverse-response model (absolute T)');
legend('Nonlinear (qc=18)', 'Inverse-response model', 'Location','best');
grid on;

% 4) 정규화 응답 비교 (ΔT/dT)
y_meas      = T_18 - T_init;          % 이미 위에서 썼던 측정 ΔT
y_meas_norm = y_meas      / dT;
y_model_norm= y_inv_interp / dT;

subplot(2,1,2);
plot(t18, y_meas_norm,   'k',  'LineWidth',1.5); hold on;
plot(t18, y_model_norm,  'r--','LineWidth',1.5);
xlabel('Time (min)');
ylabel('Normalized ΔT');
title('Normalized response: Nonlinear vs inverse-response model');
legend('Nonlinear (normalized)', 'Inverse-response model (normalized)', ...
       'Location','best');
grid on;


%% [7] (옵션) 논문식 controller R(s), C(s) 골격 만들기
% Zhang et al. 방식:
%   - 보상기 구조:
%       R(s) = (tau1 s + 1)(tau2 s + 1) / ((K1 - K2)(lambda s + 1)^2)
%   - G_c(s) = G2(s) - G1(s)
%   - unity feedback 관점에서 등가 C(s)는:
%       C(s) = R(s) / (1 + G_c(s) R(s))
%   여기서는 lambda를 사용자가 조정 (성능-강인성 trade-off)

lambda_ctrl = (tau1_hat + tau2_hat)/2;    % 예시: 중간 정도 time constant

R_num = conv([tau1_hat 1], [tau2_hat 1]);                     % (tau1 s+1)(tau2 s+1)
R_den = (K1_hat - K2_hat) * conv([lambda_ctrl 1],[lambda_ctrl 1]);  % (K1-K2)(lambda s+1)^2
R_tf  = tf(R_num, R_den);

G_c_hat = G2_hat - G1_hat;        % 논문에서의 G_c(s)
C_tf    = minreal( R_tf / (1 + G_c_hat * R_tf) );   % 등가 unity-feedback controller

disp(' ');
disp('=== Controller R(s) and equivalent C(s) (Zhang et al. style) ===');
R_tf
C_tf

% 참고용: 식별 모델 + C_tf를 쓴 폐루프 setpoint step 응답
CL_tf = feedback(C_tf * G_inv_hat, 1);   % unity feedback
figure;
step(CL_tf);
title('Closed-loop step response using Zhang-style controller C(s)');
xlabel('Time (min)');
ylabel('T (deviation)');
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

%% ===== Local function: inverse-response 모델 cost 함수 =====
function J = invresp_cost(param, t, y_meas, du)
    % param = [K1, K2, log_tau1, log_tau2]
    K1   = param(1);
    K2   = param(2);
    tau1 = exp(param(3));   % tau1 > 0 보장
    tau2 = exp(param(4));   % tau2 > 0 보장

    s = tf('s');
    G1 = K1 / (tau1*s + 1);
    G2 = K2 / (tau2*s + 1);
    G  = G1 - G2;           % inverse-response 모델

    % === 핵심 수정 부분 ===
    % step()은 "균일한 시간 벡터"를 요구하므로,
    % t(end)까지 균일한 시간축 t_sim을 새로 만든다.
    t_final = t(end);
    N_sim   = numel(t);                          % 샘플 개수를 데이터와 맞춤
    t_sim   = linspace(0, t_final, N_sim);       % 0~t_final 균일 시간축

    % 균일 시간축 t_sim 에서 step 응답 계산
    [y_step, t_step] = step(G * du, t_sim);      % 여기선 t_step == t_sim

    % 이제 y_step(t_sim)을, 실제 데이터 시간축 t로 보간
    y_model = interp1(t_step, y_step, t, 'pchip');

    % 오차 및 SSE 계산
    e = y_meas - y_model;
    J = sum(e.^2);
end
