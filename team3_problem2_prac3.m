% Team 3 – Dynamic approach to steady state for different qc (6, 12, 18 L/min)
% 초기상태에서 시작해서, qc를 각각 6 / 12 / 18 L/min으로 고정했을 때
% CSTR + 재킷 + 라인 온도가 정상상태로 어떻게 접근하는지 보는 코드

clear; clc;

%% [1] Parameters (Team #3)
V      = 27.0;          % L, reactor volume
Vc     = 6.0;           % L, jacket volume
q      = 8.0;           % L/min, feed flow rate
qc0    = 12.0;          % L/min, nominal coolant flow rate (기준점)
Caf0   = 1.3;           % mol/L, feed concentration of A
Tf     = 325.0;         % K, feed temperature
Tcin   = 315.0;         % K, coolant inlet temperature
Vline  = 9.0;           % L, line holdup volume
Tline0 = 355.0;         % K, initial line temp (문제에서 제시)
k0     = 9.0e10;        % 1/min, pre-exponential factor
EoverR = 9250.0;        % K, activation term E/R
dH     = -70000;        % J/mol, heat of reaction (exothermic) -수정함
UA     = 9000000.0;       % J/(min·K), heat transfer coefficient × area - 수정함 
rhoCp  = 18000.0;       % J/(L·K), heat capacity of reacting mixture
rhoCpc = 25000.0;       % J/(L·K), heat capacity of coolant

%% [2] (참고) qc = 12 L/min에서의 정상상태 값 (fsolve 결과)
CA_ss_12    = 1.210892;    % mol/L
T_ss_12     = 318.430;     % K 
Tc_ss_12    = 318.320;     % K
Tline_ss_12 = 315.000;     % K

%% [3] 초기 상태 설정
% "초기상태에서 정상상태로 접근하는" 걸 보고 싶은 거니까
%   - 반응기 농도: feed 농도와 같다고 가정 (CA0 = Caf0)
%   - 반응기 온도: feed 온도 (T0 = Tf)
%   - 재킷 온도: 냉각수 입구 온도 (Tc0 = Tcin)
%   - 라인 온도: 문제에서 제시된 뜨거운 체류 온도 355 K (Tline0)
CA0    = Caf0;
T0     = Tf;
Tc0    = Tcin;
% Tline0는 이미 위에서 355로 정의함.

x0 = [CA0; T0; Tc0; Tline0];

%% [4] 공통 파라미터 struct 만들기
p.V      = V;
p.Vc     = Vc;
p.q      = q;
p.qc0    = qc0;        % 참고용
p.Caf0   = Caf0;
p.Tf     = Tf;
p.Tcin   = Tcin;
p.Vline  = Vline;
p.Tline0 = Tline0;     % 참고용
p.k0     = k0;
p.EoverR = EoverR;
p.dH     = dH;
p.UA     = UA;
p.rhoCp  = rhoCp;
p.rhoCpc = rhoCpc;

%% [5] 시뮬레이션 설정
tspan = [0 60];        % min (필요하면 100, 200으로 늘려도 됨)

qc_6   = 6.0;          % L/min
qc_12  = 12.0;         % L/min
qc_18  = 18.0;         % L/min

opts = odeset('RelTol',1e-6, 'AbsTol',1e-8);

%% [6] Case 1: qc = 6 L/min (고정)

p6 = p;
p6.qc = qc_6;

odefun_6 = @(t,x) cstr_team3_odes(t, x, p6);
[t6, x6] = ode15s(odefun_6, tspan, x0, opts);

CA_6    = x6(:,1);
T_6     = x6(:,2);
Tc_6    = x6(:,3);
Tline_6 = x6(:,4);

%% [7] Case 2: qc = 12 L/min (고정, nominal)

p12 = p;
p12.qc = qc_12;

odefun_12 = @(t,x) cstr_team3_odes(t, x, p12);
[t12, x12] = ode15s(odefun_12, tspan, x0, opts);

CA_12    = x12(:,1);
T_12     = x12(:,2);
Tc_12    = x12(:,3);
Tline_12 = x12(:,4);

%% [8] Case 3: qc = 18 L/min (고정)

p18 = p;
p18.qc = qc_18;

odefun_18 = @(t,x) cstr_team3_odes(t, x, p18);
[t18, x18] = ode15s(odefun_18, tspan, x0, opts);

CA_18    = x18(:,1);
T_18     = x18(:,2);
Tc_18    = x18(:,3);
Tline_18 = x18(:,4);

%% [9] 그래프 플로팅 – 각 qc에서 정상상태로 접근하는 모양 비교

figure;

% (a) Reactor temperature T(t)
subplot(3,1,1);
plot(t6,  T_6,  'r--', 'LineWidth',1.5); hold on;   % qc=6
plot(t12, T_12, 'k-',  'LineWidth',1.5);            % qc=12
plot(t18, T_18, 'b-',  'LineWidth',1.5);            % qc=18
yline(T_ss_12, 'k:', 'LineWidth',1.0);              % 참고: qc=12 정상상태 T*
xlabel('Time (min)');
ylabel('T (K)');
title('Reactor temperature T(t): approach to steady state');
legend('qc = 6', 'qc = 12', 'qc = 18', 'T^* at qc=12', ...
       'Location','best');

% (b) Jacket temperature Tc(t)
subplot(3,1,2);
plot(t6,  Tc_6,  'r--', 'LineWidth',1.5); hold on;
plot(t12, Tc_12, 'k-',  'LineWidth',1.5);
plot(t18, Tc_18, 'b-',  'LineWidth',1.5);
yline(Tc_ss_12, 'k:', 'LineWidth',1.0);             % 참고: qc=12 정상상태 Tc*

xlabel('Time (min)');
ylabel('T_c (K)');
title('Coolant (jacket) temperature T_c(t): approach to steady state');
legend('qc = 6', 'qc = 12', 'qc = 18', 'T_c^* at qc=12', ...
       'Location','best');

% (c) Line temperature Tline(t)
subplot(3,1,3);
plot(t6,  Tline_6,  'r--', 'LineWidth',1.5); hold on;
plot(t12, Tline_12, 'k-',  'LineWidth',1.5);
plot(t18, Tline_18, 'b-',  'LineWidth',1.5);
yline(Tcin, 'k:', 'LineWidth',1.0);                 % 라인 정상상태는 결국 Tcin(=315K) 근처

xlabel('Time (min)');
ylabel('T_{line} (K)');
title('Line temperature T_{line}(t): approach to steady state');
legend('qc = 6', 'qc = 12', 'qc = 18', 'T_{cin}=315K', ...
       'Location','best');

sgtitle('Team 3: Approach to steady state from initial condition for different q_c');

%% ---- 시스템 ODE 정의 함수 ----
function dxdt = cstr_team3_odes(t, x, p)
% x = [CA; T; Tc; Tline]

CA    = x(1);   % mol/L
T     = x(2);   % K
Tc    = x(3);   % K
Tline = x(4);   % K

% ----- parameters -----
V      = p.V;
Vc     = p.Vc;
q      = p.q;
qc     = p.qc;       % 이 run에서의 냉각수 유량 (고정 값)
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
rA = k * CA;                  % mol/(L·min), A consumption rate

% ----- ODEs -----
% 1) Component balance for A
dCA = q/V * (Caf0 - CA) - rA;

% 2) Reactor energy balance
dT  = q/V * (Tf - T) ...
      + (-dH / rhoCp) * rA ...
      - UA/(rhoCp * V) * (T - Tc);

% 3) Jacket energy balance
dTc = qc/Vc * (Tline - Tc) ...
      + UA/(rhoCpc * Vc) * (T - Tc);

% 4) Coolant line energy balance
dTline = qc/Vline * (Tcin - Tline);

dxdt = [dCA; dT; dTc; dTline];
end
