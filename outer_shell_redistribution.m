function t = outer_shell_redistribution(u,vol,ind2,ind3,n,Xw_init,Xw_shell,m_weights,rho_s) %,deltaV)
% function outer_shell_equilibration is used by the MS diffucion model to
% ensure that the outer shell of the aerosol particle is always
% equilibrated with the external RH or air water content. if the outer
% shell is not equilibrated, the moles are averaged over the outer 2
% shells.

    ndig = 4;
    t = u(ind2,:).*repmat(vol(ind2)',[1 n]);
    Xw = t(ind2(end),1)./sum(t(ind2(end),:),2);
    Xw = round(Xw.*(10^ndig))./(10^ndig);
    rat = Xw_shell./(1-Xw_shell);
    
%     if (ind2(end)-ind3>2)
%         indd = ind2(end)-ind3;
%         u(ind2(end-indd:end),1)./sum(u(ind2(end-indd:end),:),2);
%         tot = sum(t(ind2(end-indd:end),:),1);
%         for i=(indd:-1:0)
%             t(ind2(end-i),1) = vol(ind2(end-i))./(m_weights(1)./rho_s(1)+m_weights(2)./rho_s(2)./rat);
%             t(ind2(end-i),2) = t(ind2(end-i),1)./rat;            
%         end
%         disp({'entered'})
        
    if any(t(ind2(end),:)<1e-40)
        tot = t(ind2(end),:)+t(ind2(end-1),:);
        t(ind2(end),1) = vol(ind2(end))./(m_weights(1)./rho_s(1)+m_weights(2)./rho_s(2)./rat);
        t(ind2(end),2) = t(ind2(end),1)./rat;
        t(ind2(end-1),1) = tot(1,1) - t(ind2(end),1);
        t(ind2(end-1),2) = tot(1,2) - t(ind2(end),2);
    
        
    elseif (Xw<Xw_shell && Xw_init>Xw_shell) %in the case of water removal
       tot = t(ind2(end),:)+t(ind2(end)-1,:);
       ind2c = ind2(end)-1;
       t(ind2(end)-1,1) = vol(ind2c)./(m_weights(1)./rho_s(1)+m_weights(2)./rho_s(2)./rat);
       t(ind2(end)-1,2) = t(ind2(end)-1,1)./rat;
       t(ind2(end),1) = tot(1,1) - t(ind2(end)-1,1);
       t(ind2(end),2) = tot(1,2) - t(ind2(end)-1,2);
 
    end
    
end