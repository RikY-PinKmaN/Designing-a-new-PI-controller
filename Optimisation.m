clear all;
%% Sampling Plan
%
%
%
%% Fullfactorial Hypercube
x = fullfactorial([10 10], 1);
% Scaling
X = rescale(x, 0.06, 0.5);
% Calculation of phi
Phi_f =  mmphi(X, 5, 2);
% Evaluation
Z_f = evaluateControlSystem(X);

%% Latin Hypercube
l = rlh(100, 2, 1);
% Scaling
L = rescale(l, 0.06, 0.5);
% Calculation of phi
Phi_l = mmphi(L, 5, 2);
% Evaluation
Z_l = evaluateControlSystem(L);

%% Optimization
%
%
%
%% Loop
for it = 1:250
    %% Fitness
    % Non Dominating Sorting
    r = rank_nds(Z_l);
    % Crodwing distance
    c_d = crowding(Z_l,r);
    %% Variation
    % Selection
    s = btwr(c_d);
    u_s = unique(s);
    % Performance
    for i = 1:length(u_s)
        row = i;
        P(i,:) = l(row,:); % Parent
    end
    % Bounds
    u_b = [max(P(:,1)), max(P(:,2))];
    l_b = [min(P(:,1)), min(P(:,2))];
    b = [l_b; u_b];
    % Child
    c = polymut(P,b);
    C = rescale(c, 0.06, 0.5);
    %% Selection for Survival
    % 2n population
    p = cat(1,P,C);
    % Evaluation
    Z = evaluateControlSystem(p);
    % Non Dominating Sorting
    r_p = rank_nds(Z);
    % Crodwing distance
    c_d_p = crowding(Z,r_p);
    % Reducer
    S = reducerNSGA_II(p,r_p,c_d_p);
    U_S = unique(S);
    %% Convergence
    for i = 1:length(U_S)
        row = i;
        nP(i,:) = p(row,:); % Parent
    end
    
end

