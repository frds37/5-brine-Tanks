clear
clc

A = [-15/500 5/1000 0 0 0
    15/500 -17/1000 2/400 0 0
    0 12/1000 -17/400 5/700 0
    0 0 15/400 -20/700 5/500
    0 0 0 15/700 -15/500];

REs = eig(A); % 실제 e-value

Evs = zeros(5,1); % e-value 저장
Exs = zeros(5,5); % e-vector 저장
step = zeros(5,1); % 수렴하는데까지 걸린 step수 저장
trEv = zeros(40,5); % e-value 추적
trEx = zeros(40,5,5); % e-vector 추적

% Power Method를 사용해 가장 큰 e-value 계산

q = rand(5,1); 
q = q/norm(q); % random unit 2-norm q 선택
Ev_prev = 1e-4; % 수렴확인 조건문 오류방지
j = 0; % step 저장에 사용

while 1
    j = j + 1;
    z = A*q;
    q = z/norm(z);
    Ev = q'*A*q;
    trEv(j,1) = Ev;
    trEx(j,:,1) = q';
    if abs(Ev-Ev_prev)/abs(Ev_prev) < 1e-5
        break
    end
    Ev_prev = Ev;
end
Evs(1) = Ev;
Exs(:,1) = q;
step(1) = j;

% Inversed power method 사용해 가장 작은 e-value 계산

q = rand(5,1);
q = q/norm(q);
Ev_prev = 1e-4;
j = 0;

while 1
    j = j + 1;
    z = A \ q;
    q = z/norm(z);
    Ev = q'*A*q;
    trEv(j,5) = Ev;
    trEx(j,:,5) = q';
    if abs(Ev-Ev_prev)/abs(Ev_prev) < 1e-5
        break
    end
    Ev_prev = Ev;
end
Evs(5) = Ev;
Exs(:,5) = q;
step(5) = j;

% Rayleigh Quotient Iteration 사용해 e-value 계산
% 
% q = rand(5,1);
% q = q/norm(q);
% Ev_prev = 1e-4;
% Ev = q'*A*q; % 초기 shift를 지정
% j = 0;
% 
% while 1
%     j = j + 1;
%     z = (A- Ev*eye(5)) \ q;
%     q = z/norm(z);
%     Ev = q'*A*q;
%     trEv(j,5) = Ev;
%     trEx(j,:,5) = q';
%     if abs(Ev-Ev_prev)/abs(Ev_prev) < 1e-5
%         break
%     end
%     Ev_prev = Ev;
% end
% Evs(5) = Ev;
% Exs(:,5) = q;
% step(5) = j;

% 위에서 최대와 최소의 e-value를 구했으므로 
% Shifted Inversed Power Method를 사용해 나머지 e-value들을 계산
% Gershgorin Theorem에 의해 모든 e-value는 diagonal element값이 중심인 
% 원 안에 존재하므로 diagonal element의 값들을 shift로 정하고 
% SIPM을 사용하는 것은 reasonable 하다고 말할 수 있음

s = [0 -0.03 -0.017 -0.0425 0]; % shift값

for i = 2:4
    j = 0;
    q = rand(5,1);
    q = q/norm(q);
    Ev_prev = 1e-4;
    while 1
        j = j + 1;
        z = (A-s(i)*eye(5)) \ q;
        q = z/norm(z);
        Ev = q'*A*q;
        trEv(j,i) = Ev;
        trEx(j,:,i) = q';
        if abs(Ev-Ev_prev)/abs(Ev_prev) < 1e-5
            break
        end
        Ev_prev = Ev;
    end
    Evs(i) = Ev;
    Exs(:,i) = q;
    step(i) = j;
end        


X_0 = [10 0 0 0 0]'; % 초기값
f = [50 0 0 0 0]'; % 상수항
v = [500 1000 400 700 500]'; % 부피
Sb = A \ (-f); % 특이해
C = Exs \ (X_0 - Sb); % 초기값을 통해 연립미분방정식의 해의 계수를 계산

t = (0:0.1:1000)'; % t=20일때
X = zeros(length(t),5); % t시간일때의 소금양을 저장
D = zeros(length(t),5); % t시간일때의 농도

for i = 1:5
    X(:,i) = (Exs(i,:) * (C .* exp(Evs*t')) + Sb(i))'; % 연립미분방정식의 해
end

for i = 1:5
    D(:,i) = X(:,i) / v(i);
end

subplot(3,2,1)
plot(t,X(:,1),t,X(:,2),t,X(:,3),t,X(:,4),t,X(:,5))
title('PM-variant t-salt plot')
xlabel('t(min)')
ylabel('salt(kg)')
legend({'x1','x2','x3','x4','x5'},'Fontsize',12)

subplot(3,2,2)
plot(t,D(:,1),t,D(:,2),t,D(:,3),t,D(:,4),t,D(:,5))
title('PM-variant t-Concentration plot')
xlabel('t(min)')
ylabel('Concentration(kg/l)')
legend({'C1','C2','C3','C4','C5'},'Fontsize',12)

tspan = [0 1000];
[t, x] = ode45(@(t,x) DE(A, x, f), tspan, X_0);
d = zeros(size(x));

for i = 1:5
    d(:,i) = x(:,i) / v(i);
end

subplot(3,2,3)
plot(t,x(:,1),t,x(:,2),t,x(:,3),t,x(:,4),t,x(:,5))
title('ode45 t-salt plot')
xlabel('t(min)')
ylabel('salt(kg)')
legend({'x1','x2','x3','x4','x5'},'Fontsize',12)

subplot(3,2,4)
plot(t,d(:,1),t,d(:,2),t,d(:,3),t,d(:,4),t,d(:,5))
title('ode45 t-Concentration plot')
xlabel('t(min)')
ylabel('Concentration(kg/l)')
legend({'C1','C2','C3','C4','C5'},'Fontsize',12)

H = A;
n = 5;
Es = zeros(5,1);
step = zeros(5,1);

for i = 1:4
    j = 0;
    while 1
        j = j + 1;
        s = H(n,n);
        [Q, R] = qr(H - s*eye(n));
        H = R * Q + s*eye(n);
        if abs(H(n,n-1)) < 1e-25
            break
        end
    end
    Es(i) = H(n,n);
    step(i) = j;
    n = n-1;
    H = H(1:n,1:n);
end
Es(5) = H(1,1);

Ex = zeros(5);
 
for i = 1:5
    B = A - Es(i)*eye(5);
    Ex(:,i) = null(B);
end

X_0 = [10 0 0 0 0]'; % 초기값
f = [50 0 0 0 0]'; % 상수항
v = [500 1000 400 700 500]'; % 부피
Sb = A \ (f); % 특이해
C = Ex \ (X_0 + Sb); % 초기값을 통해 연립미분방정식의 해의 계수를 계산

t = (0:0.1:1000)'; % t=20일때
X = zeros(length(t),5); % t시간일때의 소금양을 저장
D = zeros(length(t),5); % t시간일때의 농도

for i = 1:5
    X(:,i) = (Ex(i,:) * (C .* exp(Es*t')) - Sb(i))'; % 연립미분방정식의 해
end

for i = 1:5
    D(:,i) = X(:,i) / v(i);
end

subplot(3,2,5)
plot(t,X(:,1),t,X(:,2),t,X(:,3),t,X(:,4),t,X(:,5))
title('Hessen QR t-salt plot')
xlabel('t(min)')
ylabel('salt(kg)')
legend({'x1','x2','x3','x4','x5'},'Fontsize',12)

subplot(3,2,6)
plot(t,D(:,1),t,D(:,2),t,D(:,3),t,D(:,4),t,D(:,5))
title('Hessen QR t-Concentration plot')
xlabel('t(min)')
ylabel('Concentration(kg/l)')
legend({'C1','C2','C3','C4','C5'},'Fontsize',12)

function dxdt = DE(A, x, f)
dxdt = A * x + f;
end
