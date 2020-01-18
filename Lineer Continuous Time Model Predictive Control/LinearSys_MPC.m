function [x] = LinearSys_MPC(x,u,Ts,A,B)
xini = x;
K1X = Ts * (A * x + B * u);

% 2. Order
x = xini + K1X / 2;
K2X = Ts * (A * x + B * u);

% 3. Order
x = xini + K1X / 2;
K3X = Ts * (A * x + B * u);

% 4. Order
x = xini + K3X;
K4X = Ts * (A * x + B * u);

x = (xini) + (K1X / 6) + (K2X / 3) + (K3X / 3) + (K4X / 6);
end

