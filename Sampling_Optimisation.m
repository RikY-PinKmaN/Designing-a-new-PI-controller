clear all;
%% Sampling Plan
%
%
%
%% Fullfactorial Hypercube
x = fullfactorial([10 10], 1);
% Calculation of phi
Phi_f =  mmphi(x, 5, 2);
% Plot
subplot(1,3,1)
scatter(x(:,1), x(:,2),'red')
xlabel('Kp')
ylabel('Ki')
title('Full Factorial')

%% Latin Hypercube
l = rlh(100, 2, 1);
% Calculation of phi
Phi_l = mmphi(l, 5, 2);
% Plot
subplot(1,3,2)
scatter(l(:,1), l(:,2),'green')
xlabel('Kp')
ylabel('Ki')
title('Random Latin Hyper Cube')

%% Sobol
so = sobolset(2);
so = net(so,100);
% Calculation of phi
Phi_s = mmphi(so, 5, 2);
% Plot
subplot(1,3,3)
scatter(so(:,1), so(:,2),'blue')
xlabel('Kp')
ylabel('Ki')
title('Sobel')
% Scaling
S = rescale(so, 0.08, 0.4);
% Evaluation
Z_S = evaluateControlSystem(S);

%% Knowledge Discovery
%
%
%
%% Plot
z = [S,Z_S];
% Scatterplot
figure;
[~,a] =plotmatrix(z);

a(11,1).XLabel.String='Kp';
a(11,2).XLabel.String='Ki';
a(11,3).XLabel.String='Closed loop pole';
a(11,4).XLabel.String='Gain margin';
a(11,5).XLabel.String='Phase margin';
a(11,6).XLabel.String='Rise time';
a(11,7).XLabel.String='Peak time';
a(11,8).XLabel.String='Max overshoot';
a(11,9).XLabel.String='Max undershoot';
a(11,10).XLabel.String='Settling time';
a(11,11).XLabel.String='Steady State error';

ylabel(a(1,1),'Kp','Rotation',0,'HorizontalAlignment','right')
a(1,1).YLabel.String='Kp';
ylabel(a(2,1),'Y Axis Label','Rotation',0,'HorizontalAlignment','right')
a(2,1).YLabel.String='Ki';
ylabel(a(3,1),'Y Axis Label','Rotation',0,'HorizontalAlignment','right')
a(3,1).YLabel.String='Closed loop pole';
ylabel(a(4,1),'Y Axis Label','Rotation',0,'HorizontalAlignment','right')
a(4,1).YLabel.String='Gain margin';
ylabel(a(5,1),'Y Axis Label','Rotation',0,'HorizontalAlignment','right')
a(5,1).YLabel.String='Phase margin';
ylabel(a(6,1),'Y Axis Label','Rotation',0,'HorizontalAlignment','right')
a(6,1).YLabel.String='Rise time';
ylabel(a(7,1),'Y Axis Label','Rotation',0,'HorizontalAlignment','right')
a(7,1).YLabel.String='Peak time';
ylabel(a(8,1),'Y Axis Label','Rotation',0,'HorizontalAlignment','right')
a(8,1).YLabel.String='Max overshoot';
ylabel(a(9,1),'Y Axis Label','Rotation',0,'HorizontalAlignment','right')
a(9,1).YLabel.String='Max undershoot';
ylabel(a(10,1),'Y Axis Label','Rotation',0,'HorizontalAlignment','right')
a(10,1).YLabel.String='Settling time';
ylabel(a(11,1),'Y Axis Label','Rotation',0,'HorizontalAlignment','right')
a(11,1).YLabel.String='Steady state error';

% Parallelplot
% figure;
% lab = {'Kp','ki ','Closed loop pole','Gain margin','Phase margin','Rise time','Peak time','Max overshoot','Max undershoot','Settling time','Steady state error'};
% parallelplot(z,'CoordinateTickLabels',lab)
%% Optimization
%
%
%
%% Loop
nP = S;
nZ = Z_S;
it = 1;
while it<=250
    Z_S = optimizeControlSystem(nZ);
    %% Fitness
    % Non Dominating Sorting
    r = rank_nds(Z_S);
    % Crodwing distance
    c_d = crowding(Z_S, r);
    % Fitness
    i_r =max(r)-r(:,1);
    f = [i_r,c_d];
    %% Variation
    % Selection
    s = btwr(f);
    % Parent
    for i = 1:length(s)
        if s(i)>100
            s(i) = 100;
        end
    end
    P = nP(s,:);
    % Bounds
    u_b = [max(P(:,1)), max(P(:,2))];
    l_b = [min(P(:,1)), min(P(:,2))];
    b = [l_b; u_b];
    % Child
    C = polymut(P,b);
    %% Performance
    % 2n population
    Pop = cat(1,P,C);
    % Evaluation
    Z = evaluateControlSystem(Pop);
    Z = optimizeControlSystem(Z);
    % Non Dominating Sorting
    r_p = rank_nds(Z);
    % Crodwing distance
    c_d_p = crowding(Z, r_p);
    % Reducer
    S = reducerNSGA_II(Pop, r_p, c_d_p);
    U_S = unique(S);
    %% Convergence
    nP = Pop(U_S,:);
    nZ = evaluateControlSystem(nP);
    nZ1 = optimizeControlSystem(nZ);
    conv(it) = Hypervolume_MEX(nZ1);
    if conv == 0
        break
    end
    it = it + 1;
end
plot(conv)

%% Preference Optimization
%
%
%
%% Loop
nP = S;
nZ = Z_S;
it = 1;
% Goals
goal = [1,-6,8,2,10,10,8,20,1];
% Priority
priority = [4,3,3,2,1,2,1,1,2];
while it<=150
    Z_l = optimizeControlSystem(nZ);
    %% Fitness
    % Non Dominating Sorting
    r = rank_prf(Z_S, goal, priority);
    % Crodwing distance
    c_d = crowding(Z_S, r);
    % Fitness
    i_r =max(r)-r(:,1);
    f = [i_r,c_d];
    %% Variation
    % Selection
    s = btwr(f);
    % Parent
    P = nP(s,:);
    % Bounds
    u_b = [max(P(:,1)), max(P(:,2))];
    l_b = [min(P(:,1)), min(P(:,2))];
    b = [l_b; u_b];
    % Child
    C = polymut(P,b);
    %% Performance
    % 2n population
    Pop = cat(1,P,C);
    % Evaluation
    Z = evaluateControlSystem(Pop);
    Z = optimizeControlSystem(Z);
    % Non Dominating Sorting
    r_p = rank_prf(Z, goal, priority);
    % Crodwing distance
    c_d_p = crowding(Z, r_p);
    % Reducer
    Re = reducerNSGA_II(Pop, r_p, c_d_p);
    %% Convergence
    nP = Pop(Re,:);
    nZ = evaluateControlSystem(nP);
    nZ1 = optimizeControlSystem(nZ);
    conv(it) = Hypervolume_MEX(nZ1);
    if conv == 0
        break
    end
    it = it + 1;
end
plot(conv)