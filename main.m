clc;
clear;
close all;

%% ---------------- LOAD DATA ----------------
mpc = case14;

bus = mpc.bus;
gen = mpc.gen;
branch = mpc.branch;

nb = size(bus,1);

%% ---------------- BUS TYPES ----------------
SLACK = find(bus(:,2) == 3);
PV    = find(bus(:,2) == 2);
PQ    = find(bus(:,2) == 1);

%% ---------------- INITIAL VALUES ----------------
V = bus(:,8);
delta = deg2rad(bus(:,9));

Pd = bus(:,3)/100;
Qd = bus(:,4)/100;

Pg = zeros(nb,1);
Qg = zeros(nb,1);

for i = 1:size(gen,1)
    bus_id = gen(i,1);
    Pg(bus_id) = gen(i,2)/100;
    Qg(bus_id) = gen(i,3)/100;
end

%% ---------------- YBUS ----------------
Y = zeros(nb);

for k = 1:size(branch,1)
    i = branch(k,1);
    j = branch(k,2);
    r = branch(k,3);
    x = branch(k,4);
    b = branch(k,5);

    z = r + 1i*x;
    y = 1/z;

    Y(i,i) = Y(i,i) + y + 1i*b/2;
    Y(j,j) = Y(j,j) + y + 1i*b/2;
    Y(i,j) = Y(i,j) - y;
    Y(j,i) = Y(i,j);
end

G = real(Y);
B = imag(Y);

%% ---------------- ITERATION ----------------
max_iter = 20;
tol = 1e-6;

for iter = 1:max_iter

    P = zeros(nb,1);
    Q = zeros(nb,1);

    for i = 1:nb
        for j = 1:nb
            P(i) = P(i) + V(i)*V(j)*( G(i,j)*cos(delta(i)-delta(j)) + B(i,j)*sin(delta(i)-delta(j)) );
            Q(i) = Q(i) + V(i)*V(j)*( G(i,j)*sin(delta(i)-delta(j)) - B(i,j)*cos(delta(i)-delta(j)) );
        end
    end

    %% ---------------- MISMATCH ----------------
    dP = Pg - Pd - P;
    dQ = Qg - Qd - Q;

    % Only include required equations
    mismatch = [dP([PV; PQ]); dQ(PQ)];

    if max(abs(mismatch)) < tol
        disp(['Converged in ', num2str(iter), ' iterations']);
        break;
    end

    %% ---------------- JACOBIAN ----------------
    npv = length(PV);
    npq = length(PQ);

    J1 = zeros(npv+npq);
    J2 = zeros(npv+npq, npq);
    J3 = zeros(npq, npv+npq);
    J4 = zeros(npq);

    buses = [PV; PQ];

    % J1 (dP/dδ)
    for i = 1:length(buses)
        m = buses(i);
        for j = 1:length(buses)
            n = buses(j);
            if m == n
                for k = 1:nb
                    if k ~= m
                        J1(i,j) = J1(i,j) + V(m)*V(k)*(-G(m,k)*sin(delta(m)-delta(k)) + B(m,k)*cos(delta(m)-delta(k)));
                    end
                end
                J1(i,j) = J1(i,j) - Q(m);
            else
                J1(i,j) = V(m)*V(n)*(G(m,n)*sin(delta(m)-delta(n)) - B(m,n)*cos(delta(m)-delta(n)));
            end
        end
    end

    % J2 (dP/dV)
    for i = 1:length(buses)
        m = buses(i);
        for j = 1:npq
            n = PQ(j);
            if m == n
                J2(i,j) = (P(m)/V(m)) + G(m,m)*V(m);
            else
                J2(i,j) = V(m)*(G(m,n)*cos(delta(m)-delta(n)) + B(m,n)*sin(delta(m)-delta(n)));
            end
        end
    end

    % J3 (dQ/dδ)
    for i = 1:npq
        m = PQ(i);
        for j = 1:length(buses)
            n = buses(j);
            if m == n
                for k = 1:nb
                    if k ~= m
                        J3(i,j) = J3(i,j) + V(m)*V(k)*(G(m,k)*cos(delta(m)-delta(k)) + B(m,k)*sin(delta(m)-delta(k)));
                    end
                end
                J3(i,j) = J3(i,j) + P(m);
            else
                J3(i,j) = -V(m)*V(n)*(G(m,n)*cos(delta(m)-delta(n)) + B(m,n)*sin(delta(m)-delta(n)));
            end
        end
    end

    % J4 (dQ/dV)
    for i = 1:npq
        m = PQ(i);
        for j = 1:npq
            n = PQ(j);
            if m == n
                J4(i,j) = (Q(m)/V(m)) - B(m,m)*V(m);
            else
                J4(i,j) = V(m)*(G(m,n)*sin(delta(m)-delta(n)) - B(m,n)*cos(delta(m)-delta(n)));
            end
        end
    end

    J = [J1 J2; J3 J4];

    %% ---------------- UPDATE ----------------
    dx = J \ mismatch;

    d_delta = dx(1:(npv+npq));
    d_V     = dx((npv+npq)+1:end);

    delta([PV; PQ]) = delta([PV; PQ]) + d_delta;
    V(PQ) = V(PQ) + d_V;
end

%% ---------------- RESULTS ----------------
disp('Final Bus Voltages:');
disp(V);

disp('Final Voltage Angles (deg):');
disp(rad2deg(delta));

%% ---------------- PLOTS ----------------
figure;
plot(V,'o-','LineWidth',2);
title('Bus Voltage Profile');
xlabel('Bus Number');
ylabel('Voltage (p.u)');
grid on;

figure;
plot(rad2deg(delta),'*-','LineWidth',2);
title('Voltage Angle');
xlabel('Bus Number');
ylabel('Angle (deg)');
grid on;