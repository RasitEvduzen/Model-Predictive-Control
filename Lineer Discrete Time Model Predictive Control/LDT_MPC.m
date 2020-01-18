clc,clear all,close all,warning off;
%% Lineer Discrete Time Model Predictive Control
% 18-Jan-2017
%% System Modeling
A = [1.1 0.6;0.2 0.1];
b = [1;0];
c = [1 1]';

%% Algorithm Parameters
N = 60;
umin = -1;       % input Constrain
umax = 1;        % input Constrain
ymin = 0;        % Output Constrain
ymax = 1;        % Output Constrain

Ku = 1;          % Control Horizon
Ky = 5;          % Prediction Horizon
deltaumax = 0.2; % input Change Speed Constrain
lamda = 1;     % Penalty Parameter For input 
% Optimum Lamda
% lamda = 0.1*(Ky*(ymax-ymin)^2)/(Ku*(umax-umin)^2);

x(:,1) = [0;0];         % State initial Parameters
uini = zeros(Ku+1,1);   % input initial Parameters
u(:,1) = uini;


%% Referance input 
% yref = 1*ones(N+Ky-1,1);        % Step Referance Signal
yref = sin((1:N+4)'*0.1*pi);  % Sin Referance Signal
% yref = [0.5*ones(N/(4),1);ones(N/(4),1);0.6*ones(N/(4),1);0.25*ones(N/(4)+4,1)];


%% Linear MPC Matrices
[Haxx,L,M,Z] = LDTS_MPC_matrices(A,b,c,Ku,Ky,lamda);

% Hessian and inverse 
Hes = (M'*M + lamda * L);
Hters = inv(Hes);
tic
%% Algorithm Start
for i=1:N
    g = (M' * M + lamda * L) * u(:,i) - (M' * (yref(i:i+Ky-1) - Z * x(:,i)) + [lamda * u(1,i);zeros(Ku,1)]);   % Gradient Calculation
    deltau = -Hters * g;
    mu(i) = deltaumax / max([max(abs(deltau)),max(abs(diff(u(:,i) + deltau)))]);
    if (mu(i) > 1)
        mu(i) = 1;
    end
    
    deltau = mu(i) * deltau;
    u(:,i+1) = u(:,i) + deltau;
    if (u(1,i+1) > umax)
        u(1,i+1) = umax;
    end
    if (u(1,i+1) < umin)
        u(1,i+1) = umin;
    end
    
    x(:,i+1) = A * x(:,i) + b * u(1,i+1);
    y(i+1) = c' * x(:,i+1);
end
Elapsed_Time = toc;
%% Plot Output
subplot(221)
grid on, hold on;
plot(yref,'r','LineWidth',1)
plot(y,'b','LineWidth',1)
xlabel('n');
ylabel('Y_{ref}[n],Y[n]');
title(['K_{y} = ',num2str(Ky),' | ','K_{u} =',num2str(Ku),' | ','lamda = ',num2str(lamda),' | ','u_{min} = ',num2str(umin),' | ','u_{max} = ',num2str(umax),' | ','du_{max} = ',num2str(deltaumax),' | ','Elapsed Time = ',num2str(Elapsed_Time)]);
legend('= Yref[n]','= Y[n]');
% axis([0 inf 0 2])

subplot(222)
Error = yref(1:length(y),:) - y';
plot(Error,'b','LineWidth',1)
grid on
title('Referance Tracking Error')
% axis([0 inf -1.5 1.5])

subplot(223)
grid on, hold on;
plot(u(1,:),'b','LineWidth',1)
xlabel('n');
ylabel('u[n]');
title('Ýnput')
% axis([0 inf -1.5 1.5])

subplot(224)
grid on, hold on;
plot(x(1,:),'k--','LineWidth',1)
plot(x(2,:),'r--','LineWidth',1)
xlabel('n');
ylabel('x_{i}[n]');
title('State')
legend('= State 1','= State 2')
% axis([0 inf 0 1.3])




