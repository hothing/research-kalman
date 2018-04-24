// Testa kalman filtration for a simple movement model

T = 1 // sampling time [sec]

// state space model

A = [1, T; 0, 1]

B = [T^2; T]

C = [1, 0; 0, 1]

// check a controlablility

R = [B, A*B]

rank(R)

// check an observability

O = [C; C*A]

rank(O)

// simulation a liner movement with constant speed

n = 500

x = zeros(2, n)
y = zeros(2, n)
u = zeros(1, n)

x(:,1) = [1500; 0]
u(1) = -0.33
u(2) = u(1)

for i = 2:n do
    x(:,i) = A * x(:, i - 1) + B*u(:, i)
    y(:,i) = C * x(:, i)  
end

// simulate a measurement with noise
sigma_pos = 3.5
sigma_vel = 0.1

n1 = grand(1, n, "nor", 0, sigma_pos)
n2 = grand(1, n, "nor", 0, sigma_vel)

ns = [n1; n2]

z = y + ns


// Kalman filter
Q = [0.1, 0; 0, 0.01]
S = [4^2, 0; 0, 0.2^2]

xhat = zeros(2,n)
P = zeros(2, 2)

for i = 2:n do
    xhatminus = A * xhat(:, i - 1) + B * u(:, i - 1)
    Pminus = A * P * A' + Q
    Kp = (C * Pminus * C' + S)^-1
    K = Pminus * C' * Kp
    xhat(:, i) = xhatminus + K * (z(:,i) - C * xhatminus)
    P = (eye(K) - K*C) * Pminus
end


    
