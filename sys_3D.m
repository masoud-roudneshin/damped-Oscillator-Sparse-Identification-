clc 
clear 
 
% This is a MATLAB code for examples of nonlinear system identification.
% We use measurements form the system state to recover the governing
% equations. For more info see: 
% https://www.pnas.org/action/downloadSupplement?doi=10.1073%2Fpnas.1517384113&file=pnas.1517384113.sapp.pdf
% This code implements identification for example 1.b 
% Objective: identify a system with the following dynamics from pure state
% measurements
% dX/dt = [-0.1 -2 0;2 -0.1 0;0 0 -0.3] * [x_1; x_2; x_3]
% For system recovery we assume that the base functions are polynomials up
% to degree 3
%% Simulation Parameters
delT = 0.001; % sampling time for the original data
delT_Measurement = 0.002; % sampling time for the state measurement
simTime = delT: delT: 25; % simulation time
 
X(:,1) = [2; 0; 1]; % initial state 
XX = []; 
dXdt = [];  
polyDeg = 3; 
 
XOld = X(:,1); 
 
for i = 2 : length(simTime) 
    delX1 = [-0.1 -2 0;2 -0.1 0;0 0 -0.3] * X(:,i-1) * delT; 
    X(:,i) = X(:,i-1) + delX1; 

    if mod(simTime(i),delT_Measurement) == 0 
        tempPolyBase = polynomialBase_3States(X(1,i), X(2,i), X(3,i), polyDeg); 
        XX = [XX;tempPolyBase]; 
        tempdXdt = (X(:,i) - XOld)/delT_Measurement; 
        dXdt = [dXdt;tempdXdt']; 
        XOld = X(:,i); 

    end 
end 
 
  
 
%% compute Sparse regression: sequential least squares 
% Despite the large size of the measurement matrix
% XX, we assume that the coefficients matrix is sparse and can be recovered
% with some optimization that imposes sparsity. 
% One solution is L1 norm minimization (LASSO) which 
% is quite expensive given the size of the measrement matrix. 
% We use a smart method through succesive least-squares by making zero 
% the solutions that are smaller than a threshold lambda. 

Xi = XX\dXdt; % initial guess: Least-squares 
lambda = 0.01; 
stateDim = 2; 
% lambda is our sparsification knob. 
for k = 1:10 
    smallinds = (abs(Xi) < lambda); % find small coefficients 
    Xi(smallinds) = 0; % and threshold 
    for ind = 1: stateDim % n is state dimension 
        biginds = ~smallinds(:,ind); 
        % Regress dynamics onto remaining terms to find sparse Xi 
        Xi(biginds,ind) = XX(:,biginds)\dXdt(:,ind); 
    end 
end 
 
 
%% Signal Recovery 
XRecovered(:,1) = X(:,1); 

for i = 2: size(XX,1) 
    tempPolyBase = polynomialBase_3States(XRecovered(1,i-1), XRecovered(2,i-1),  XRecovered(3,i-1), polyDeg); 
    dXRecovered = Xi' * tempPolyBase'; 
    XRecovered(:,i) = XRecovered(:,i-1) + dXRecovered * delT_Measurement; 
 
end 
 
simTime_Measurement = delT_Measurement: delT_Measurement: delT_Measurement* size(XX,1); 
 
 
%% Plots 
subplot(2,1,1) 
hold on 
plot(simTime,X(1,:), 'LineWidth',2, "LineStyle","--") 
plot(simTime,X(2,:), 'LineWidth',2, "LineStyle","--") 
plot(simTime,X(3,:), 'LineWidth',2, "LineStyle","--") 

plot(simTime_Measurement,XRecovered(1,:), 'LineWidth',2) 
plot(simTime_Measurement,XRecovered(2,:), 'LineWidth',2) 
plot(simTime_Measurement,XRecovered(3,:), 'LineWidth',2) 

grid on
xlabel("Time")
ylabel("x_k")
legend("x_1", "x_2","x_3", "recovered x_1", "recovered x_2", "recovered x_3")
 
subplot(2,1,2) 
hold on 
hold on
plot3(X(1,:),X(2,:),X(3,:), 'LineWidth',2, "LineStyle","--") 
plot3(XRecovered(1,:),XRecovered(2,:),XRecovered(3,:), 'LineWidth',2) 
grid on
xlabel("x_1")
ylabel("x_2")
ylabel("x_3")
legend("Original Data", "Recovered Data")