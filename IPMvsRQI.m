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

% Inversed power method 사용해 가장 작은 e-value 계산
% Rayleigh Quotient Iteration 보다 느려 여기서 사용하지는 않음

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

% Rayleigh Quotient Iteration 사용해 가장 작은 e-value 계산

q0 = rand(5,1);
q = q0/norm(q0);
Ev_prev = 1e-4;
Ev = q'*A*q; % 초기 shift를 지정
j = 0;

while 1
    j = j + 1;
    z = (A - Ev*eye(5)) \ q;
    q = z/norm(z);
    Ev = q'*A*q;
    trEv(j,1) = Ev;
    trEx(j,:,1) = q';
    if abs(Ev-Ev_prev)/abs(Ev_prev) < 1e-6
        break
    end
    Ev_prev = Ev;
end
Evs(1) = Ev;
Exs(:,1) = q;
step(1) = j;
