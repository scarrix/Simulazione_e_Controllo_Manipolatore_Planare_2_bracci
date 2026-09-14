clc;
clear;

% PD compensazione gravita + Feedforward


%% Definizione parametri robot
a1 = 1.0;   a2 = 0.6;
l1 = a1/2;  l2 = a2/2;
m1 = 5;     m2 = 3;
I1 = 14;    I2 = 8;
g  = 9.81;

%% Parametri di controllo
% Legge PD + compensazione gravità + feedforward nello spazio dei giunti:
%   tau = g(q) + B(q)*q_ddot_d + Kp*(qd - q) - Kd*q_dot
%
% Il termine B(q)*q_ddot_d (feedforward) anticipa la coppia necessaria
% per seguire la traiettoria, riducendo drasticamente l'errore di
% inseguimento anche con guadagni moderati e realistici.
% Kp = 190 * eye(2);
Kp = 300 * eye(2);
Kd =  40 * eye(2);

%% Parametri della traiettoria
tf = 10;
N  = 10000;
t  = linspace(0, tf, N);
dt = tf / (N-1);

q0 = [pi/4;  pi/3];
qf = [pi/2; -pi/2];

% -------------------------------------------------------
%% 1) PIANIFICAZIONE TRAIETTORIA (punto 1 traccia)
% -------------------------------------------------------

% --- Giunto 1: polinomio cubico ---
% Condizioni al contorno: q(0)=q0, q(tf)=qf, q_dot(0)=0, q_dot(tf)=0
a0_1 =  q0(1);
a1_1 =  0;
a2_1 =  3*(qf(1)-q0(1)) / tf^2;
a3_1 = -2*(qf(1)-q0(1)) / tf^3;

q1_d      = a0_1 + a1_1*t    + a2_1*t.^2 + a3_1*t.^3;
q1_dot_d  =        a1_1      + 2*a2_1*t   + 3*a3_1*t.^2;
q1_ddot_d =                    2*a2_1     + 6*a3_1*t;

% --- Giunto 2: profilo trapezoidale (Dq2 < 0) ---
% ta = periodo di accelerazione = tf/4
ta2 = tf / 4;
Dq2 = qf(2) - q0(2);          % escursione totale (negativa)
vm2 = Dq2 / (tf - ta2);       % velocità di crociera
ac2 = vm2 / ta2;               % accelerazione (negativa)

q2_d      = zeros(1, N);
q2_dot_d  = zeros(1, N);
q2_ddot_d = zeros(1, N);

for i = 1:N
    ti = t(i);
    if ti <= ta2
        % Fase di accelerazione
        q2_d(i)      = q0(2) + 0.5*ac2*ti^2;
        q2_dot_d(i)  = ac2*ti;
        q2_ddot_d(i) = ac2;
    elseif ti <= tf - ta2
        % Fase a velocità costante
        q2_d(i)      = q0(2) + 0.5*ac2*ta2^2 + vm2*(ti - ta2);
        q2_dot_d(i)  = vm2;
        q2_ddot_d(i) = 0;
    else
        % Fase di decelerazione
        q2_d(i)      = qf(2) - 0.5*ac2*(tf - ti)^2;
        q2_dot_d(i)  = ac2*(tf - ti);
        q2_ddot_d(i) = -ac2;
    end
end
q2_d(end) = qf(2);   % correzione numerica al punto finale

% -------------------------------------------------------
%% 2+3) SIMULAZIONE (dinamica diretta) + CONTROLLO PD+g+feedforward
%       NELLO SPAZIO DEI GIUNTI (punti 2 e 3 traccia)
% -------------------------------------------------------
q     = q0;
q_dot = zeros(2,1);

salva_q     = zeros(2, N);
salva_q_dot = zeros(2, N);
salva_tau   = zeros(2, N);
salva_e     = zeros(2, N);
salva_x     = zeros(2, N);

for i = 1:N
    % Traiettoria desiderata al passo corrente
    qd = [q1_d(i); q2_d(i)];

    % Errore di posizione nello spazio dei giunti
    e = qd - q;

    % Modello dinamico
    B  = calcolaB(q, m1, m2, l1, l2, a1, I1, I2);
    C  = calcolaC(q, q_dot, m2, l2, a1);
    gq = calcolaG(q, m1, m2, l1, l2, a1, g);

    % Accelerazione desiderata al passo corrente
    q_ddot_d = [q1_ddot_d(i); q2_ddot_d(i)];

    % Legge di controllo PD + compensazione gravità + feedforward:
    %   tau = g(q) + B(q)*q_ddot_d + Kp*(qd - q) - Kd*q_dot
    % Il termine B(q)*q_ddot_d anticipa la coppia necessaria per
    % seguire la traiettoria (feedforward sull'accelerazione desiderata).
    tau = gq + B*q_ddot_d + Kp*e - Kd*q_dot;

    % Dinamica diretta: q_ddot = B^{-1} * (tau - C*q_dot - g(q))
    q_ddot = B \ (tau - C*q_dot - gq);

    % Salvataggio stato
    salva_x(:, i)     = cinematica_diretta(q, a1, a2);
    salva_q(:, i)     = q;
    salva_q_dot(:, i) = q_dot;
    salva_tau(:, i)   = tau;
    salva_e(:, i)     = e;

    % Integrazione con metodo di Eulero esplicito
    q_dot = q_dot + q_ddot * dt;
    q     = q     + q_dot  * dt;
end

% Stampa errori finali
fprintf('\n=== Errori finali a t=%.1f s ===\n', tf);
fprintf('e1(tf) = %+.2e rad  (%+.4f deg)\n', salva_e(1,end), rad2deg(salva_e(1,end)));
fprintf('e2(tf) = %+.2e rad  (%+.4f deg)\n', salva_e(2,end), rad2deg(salva_e(2,end)));
fprintf('IAE1   = %.6f rad*s\n', sum(abs(salva_e(1,:)))*dt);
fprintf('IAE2   = %.6f rad*s\n', sum(abs(salva_e(2,:)))*dt);

% -------------------------------------------------------
%% FIGURE
% -------------------------------------------------------

fs_legend = 13;

% --- Figure 1: Profilo cubico giunto 1 ---
figure(1);
subplot(2,1,1);
plot(t, q1_dot_d, 'b', 'LineWidth', 2);
title('Profilo di velocità polinomiale cubica - Giunto 1', 'FontSize', 13);
xlabel('Tempo (s)', 'FontSize', 12);
ylabel('Velocità (rad/s)', 'FontSize', 12);
grid on;

subplot(2,1,2);
plot(t, q1_ddot_d, 'b', 'LineWidth', 2);
title('Profilo di accelerazione - Giunto 1', 'FontSize', 13);
xlabel('Tempo (s)', 'FontSize', 12);
ylabel('Accelerazione (rad/s^2)', 'FontSize', 12);
grid on;

% --- Figure 2: Profilo trapezoidale giunto 2 ---
figure(2);
subplot(2,1,1);
plot(t, q2_dot_d, 'r', 'LineWidth', 2);
title('Profilo di velocità trapezoidale - Giunto 2', 'FontSize', 13);
xlabel('Tempo (s)', 'FontSize', 12);
ylabel('Velocità (rad/s)', 'FontSize', 12);
grid on;

subplot(2,1,2);
plot(t, q2_ddot_d, 'r', 'LineWidth', 2);
title('Profilo di accelerazione - Giunto 2', 'FontSize', 13);
xlabel('Tempo (s)', 'FontSize', 12);
ylabel('Accelerazione (rad/s^2)', 'FontSize', 12);
grid on;

% --- Figure 3: Velocità dei giunti ---
figure(3);
plot(t, salva_q_dot(1,:), 'b',   'LineWidth', 2); hold on;
plot(t, salva_q_dot(2,:), 'r',   'LineWidth', 2);
plot(t, q1_dot_d,          'b--', 'LineWidth', 1.5);
plot(t, q2_dot_d,          'r--', 'LineWidth', 1.5);
leg = legend('\dot{q}_1', '\dot{q}_2', '\dot{q}_{1d}', '\dot{q}_{2d}', ...
             'Location', 'best');
set(leg, 'FontSize', fs_legend);
title('Velocità dei giunti', 'FontSize', 13);
xlabel('Tempo (s)', 'FontSize', 12);
ylabel('Velocità (rad/s)', 'FontSize', 12);
grid on; hold off;

% --- Figure 4: Errore di posizione ---
figure(4);
plot(t, salva_e(1,:), 'b', t, salva_e(2,:), 'r', 'LineWidth', 2);
yline(0, 'k--', 'LineWidth', 1);
leg = legend('e_1 = q_{1d} - q_1', 'e_2 = q_{2d} - q_2', 'Location', 'best');
set(leg, 'FontSize', fs_legend);
title('Errore di posizione dei giunti', 'FontSize', 13);
xlabel('Tempo (s)', 'FontSize', 12);
ylabel('Errore (rad)', 'FontSize', 12);
e_max = max(max(abs(salva_e)));
ylim([-1.2*e_max, 1.2*e_max]);
grid on;

% --- Figure 5: Coppie di controllo ---
figure(5);
plot(t, salva_tau(1,:), 'b', t, salva_tau(2,:), 'r', 'LineWidth', 2);
leg = legend('\tau_1', '\tau_2', 'Location', 'best');
set(leg, 'FontSize', fs_legend);
title('Coppie di controllo', 'FontSize', 13);
xlabel('Tempo (s)', 'FontSize', 12);
ylabel('Coppia (Nm)', 'FontSize', 12);
grid on;

% --- Figure 6: Traiettoria spazio operativo ---
x_des = zeros(2, N);
for i = 1:N
    x_des(:,i) = cinematica_diretta([q1_d(i); q2_d(i)], a1, a2);
end
figure(6);
plot(x_des(1,:), x_des(2,:), 'r', 'LineWidth', 2); hold on;
plot(salva_x(1,:), salva_x(2,:), 'b--', 'LineWidth', 2);
scatter(x_des(1,1),   x_des(2,1),   80, 'go', 'filled');
scatter(x_des(1,end), x_des(2,end), 80, 'ro', 'filled');
leg = legend('Traiettoria desiderata', 'Traiettoria effettiva', ...
             'Punto di partenza', 'Punto di arrivo', 'Location', 'best');
set(leg, 'FontSize', fs_legend);
title('Traiettoria end-effector nello spazio operativo', 'FontSize', 13);
xlabel('X (m)', 'FontSize', 12);
ylabel('Y (m)', 'FontSize', 12);
axis equal; grid on; hold off;

% --- Figure 7: Posizione end-effector nel tempo ---
figure(7);
plot(t, x_des(1,:),   'r',   'LineWidth', 2); hold on;
plot(t, x_des(2,:),   'r--', 'LineWidth', 2);
plot(t, salva_x(1,:), 'b',   'LineWidth', 2);
plot(t, salva_x(2,:), 'b--', 'LineWidth', 2);
leg = legend('x_{des}', 'y_{des}', 'x', 'y', 'Location', 'best');
set(leg, 'FontSize', fs_legend);
title('Posizione end-effector nel tempo', 'FontSize', 13);
xlabel('Tempo (s)', 'FontSize', 12);
ylabel('Posizione (m)', 'FontSize', 12);
grid on; hold off;

% --- Figure 8: Angoli dei giunti ---
figure(8);
plot(t, salva_q(1,:), 'b',   'LineWidth', 2); hold on;
plot(t, salva_q(2,:), 'r',   'LineWidth', 2);
plot(t, q1_d,         'b--', 'LineWidth', 1.5);
plot(t, q2_d,         'r--', 'LineWidth', 1.5);
leg = legend('q_1', 'q_2', 'q_{1d}', 'q_{2d}', 'Location', 'best');
set(leg, 'FontSize', fs_legend);
title(sprintf('Angoli dei giunti  [Kp=%d  Kd=%d  + feedforward]', Kp(1,1), Kd(1,1)), 'FontSize', 13);
xlabel('Tempo (s)', 'FontSize', 12);
ylabel('Angolo (rad)', 'FontSize', 12);
grid on; hold off;

% --- Figure 9: Confronto q vs qd per giunto ---
figure(9);
subplot(2,1,1);
plot(t, q1_d, 'b', t, salva_q(1,:), 'r--', 'LineWidth', 2);
leg = legend('q_{1d}', 'q_1', 'Location', 'best');
set(leg, 'FontSize', fs_legend);
title('Giunto 1: posizione desiderata vs simulata', 'FontSize', 13);
xlabel('Tempo (s)', 'FontSize', 12);
ylabel('Angolo (rad)', 'FontSize', 12);
grid on;

subplot(2,1,2);
plot(t, q2_d, 'b', t, salva_q(2,:), 'r--', 'LineWidth', 2);
leg = legend('q_{2d}', 'q_2', 'Location', 'best');
set(leg, 'FontSize', fs_legend);
title('Giunto 2: posizione desiderata vs simulata', 'FontSize', 13);
xlabel('Tempo (s)', 'FontSize', 12);
ylabel('Angolo (rad)', 'FontSize', 12);
grid on;


% =========================================================
%% FUNZIONI LOCALI
% =========================================================

function B = calcolaB(q, m1, m2, l1, l2, a1, I1, I2)
% Matrice di inerzia B(q) del manipolatore planare 2-DOF
    th2 = q(2);
    B11 = I1 + m1*l1^2 + I2 + m2*(a1^2 + l2^2 + 2*a1*l2*cos(th2));
    B12 = I2 + m2*(l2^2 + a1*l2*cos(th2));
    B22 = I2 + m2*l2^2;
    B = [B11, B12; B12, B22];
end

function C = calcolaC(q, q_dot, m2, l2, a1)
% Matrice di Coriolis/centrifuga C(q, q_dot)
    th2  = q(2);
    dth1 = q_dot(1);
    dth2 = q_dot(2);
    h = m2*a1*l2*sin(th2);
    C = [-h*dth2,  -h*(dth1+dth2);
          h*dth1,   0];
end

function gq = calcolaG(q, m1, m2, l1, l2, a1, g)
% Vettore delle forze gravitazionali g(q)
    th1 = q(1); th2 = q(2);
    g1 = (m1*l1 + m2*a1)*g*cos(th1) + m2*l2*g*cos(th1 + th2);
    g2 = m2*l2*g*cos(th1 + th2);
    gq = [g1; g2];
end

function pos = cinematica_diretta(q, a1, a2)
% Cinematica diretta: restituisce [x; y] dell'end-effector
    th1 = q(1); th2 = q(2);
    x = a1*cos(th1) + a2*cos(th1 + th2);
    y = a1*sin(th1) + a2*sin(th1 + th2);
    pos = [x; y];
end