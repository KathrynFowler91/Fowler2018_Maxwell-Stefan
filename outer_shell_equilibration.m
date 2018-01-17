function [DeltaV, shell_flag] = outer_shell_equilibration(u,ind2,vol,Xw_shell,n,rho_s,m_weights,r0,rN,ip)
% function outer_shell_equilibration is used by the MS diffucion model to
% ensure that the outer shell of the aerosol particle is always
% equilibrated with the external RH or air water content.

num_shells = 1;
DeltaMin = 0.1;

% shell indicies
% nend = ind2(end);
% ni = ind2(end+1-num_shells);

% number of moles in the shells
t = u(ind2,:).*repmat(vol(ind2)',[1 n]);
t2 = t(end,:);

% change in the number of moles
%Deltat = ((1-Xw_shell).*t(end,1)-Xw_shell.*t(end,2))./(Xw_shell-1);
Deltat = ((1-Xw_shell).*t2(1)-Xw_shell.*t2(2))./(Xw_shell-1);
DeltaV = Deltat./rho_s(1).*m_weights(1);
shell_flag = 1;

% change in volume check it is greater than defined amount
RN = linspace(r0,rN,ip+1);
vollayer = 4./3.*pi.*(RN(ind2(end))-RN(ind2(end)-1));


if (abs(DeltaV)<DeltaMin.*vollayer && DeltaV<0)
    % calculate change in water volume based on the final two shells
    t2 = t(end,:)+t(end-1,:);
    Deltat = ((1-Xw_shell).*t2(1)-Xw_shell.*t2(2))./(Xw_shell-1);
    DeltaV = Deltat./rho_s(1).*m_weights(1);
    shell_flag = 2;
end


end