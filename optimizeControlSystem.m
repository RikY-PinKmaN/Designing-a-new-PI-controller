function Ev_Z= optimizeControlSystem(Z1)
% Phase margin

Z1(:,3) = abs(45 - Z1(:,3));

% Gain margin
Z1(:,2) = (-1).*(20*log10(Z1(:,2)));

Ev_Z = Z1;

end


