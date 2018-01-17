%--UNIFAC--
function [gamma_i] = UNIFAC_gamma(molecules,x_i,sub_groups,molecule_group_flag,molecule_group_stoich,T)
%*****************ambient conditions*******************
%T=318;

%***************component specification****************
%molecules=int8(molecules);

%--acetone--
%molecule_group_flag(1,1)=1; %ch3
%molecule_group_flag(1,2)=19; %ch3co
%sub_groups(1)=2;

%molecule_group_stoich(1,1)=1;
%molecule_group_stoich(1,2)=1;

%--chloroform--
%molecule_group_flag(2,1)=45; %ch3cl
%sub_groups(2)=1;

%molecule_group_stoich(2,1)=1;

%--n-Hexane--
%molecule_group_flag(3,1)=1;  %ch3
%molecule_group_flag(3,2)=2; %ch2
%sub_groups(3)=2;

%molecule_group_stoich(3,1)=2;
%molecule_group_stoich(3,2)=4;

%*******define mole fractions of each compound********
%x_i(1)=0.599;
%x_i(2)=0.102;
%x_i(3)=0.299;


%*********************************************************
%pull out the unifac parameters and save to the workspace:
save_step=0;
if save_step==0
UFCData_ri_save = UFC_datari_v1();
UFCData_qi_save = UFC_dataqi_v1();
UFCData_Main_save = UFC_datamain();
UFCData2_save = UFC_interact_params_v1();
save_step=save_step+1;
end

%*********************************************************
%pull out the volume, size and stoichiometric information
vk_i(:,:)=0;
group_step=1;
func_groups=0;
for mol=1:molecules
%order the functional groups in some way
%[r,c]=size(molecule_group_flag(mol,:));
 for i=1:sub_groups(mol)
     %here look at each group identifier then pull out
     %the group number and stoichiometry seperately
     %dosnt matter if a group is repeated as the interactions
     %between each are set to zero anyway (which you can do manually)
     group_flag_array(group_step)=molecule_group_flag(mol,i);
     Q_k(group_step)=UFCData_qi_save(molecule_group_flag(mol,i));
     R_k(group_step)=UFCData_ri_save(molecule_group_flag(mol,i));
     vk_i(mol,group_step)=molecule_group_stoich(mol,i);
     %mol_func_groups_rec(i)=molecule_group_flag(mol,i)
     group_step=group_step+1;
     func_groups=func_groups+1;
  end
end
%*********************************************************
%pull out interaction parameters based on the groups set out above
a_m_n(:,:)=0.0;
[r,c]=size(Q_k);
for i=1:c
    for j=1:c
    a_m_n(i,j)=UFCData2_save(UFCData_Main_save(group_flag_array(i)),UFCData_Main_save(group_flag_array(j)));
    end
end


%--------------defining groups-------------------------
%-------groups involved with the biogenic compound-----
%--ch3--
%Q_k(1)=0.848;
%R_k(1)=0.9011;
%--ch2--
%Q_k(2)=0.54;
%R_k(2)=0.6744; 
%--ch--
%Q_k(3)=0.228;
%R_k(3)=0.4469;
%--c--
%Q_k(4)=0; 
%R_k(4)=0.2195;
%--oh--
%Q_k(5)=1.2;
%R_k(5)=1;
%--ch3co--
%Q_k(6)=1.488;
%R_k(6)=1.6724;
%--cho--
%Q_k(7)=0.948;
%R_k(7)=0.998;
%--cooh--
%Q_k(8)=1.224;
%R_k(8)=1.3013;
%nk=1 (CH3)
%nk=2 (CH2)
%nk=3 (CH3CO)
%nk=4 (CHCL3)



%------------interaction parameters--------------------
% this loop pulls out the relevant interaction parameters for the groups
% you are considering

%a_m_n(1,3)=476.4;
%a_m_n(1,4)=24.90;

%a_m_n(2,3)=476.4;
%a_m_n(2,4)=24.90;

%a_m_n(3,1)=26.76;
%a_m_n(3,2)=26.76;
%a_m_n(3,4)=-354.6;

%a_m_n(4,1)=36.7;
%a_m_n(4,2)=36.7;
%a_m_n(4,3)=552.1;

%--------------define compounds being considered--------
%functional group stoichiometric coefficients
%first set all matrices to zero 
%vk_i(:,:)=0;
%acetone
%vk_i(1,1)=1;
%vk_i(1,3)=1;
%chloroform
%vk_i(2,4)=1;
%n-Hexane
%vk_i(3,1)=2;
%vk_i(3,2)=4;

%compound flags
%flag(:,:)=0;
%flag(1,1)=1;
%flag(1,3)=1;
%flag(2,4)=1;
%flag(3,1)=1;
%flag(3,2)=1;


%********************************************************************
%A)--calculate the combinatorial contribution to activity coefficient
%--calculate volume parameters--
for mol=1:molecules
    for k=1:func_groups
      Q_k_i(mol,k)=Q_k(k)*vk_i(mol,k);
    end
    q_i(mol)=sum(Q_k_i(mol,:));
end

Q=sum(x_i.*q_i);
for mol=1:molecules
    w_i(mol)=q_i(mol)/Q;
end
omega_i=log(w_i);

%--calculate surface parameters--
for mol=1:molecules
    for k=1:func_groups
        R_k_i(mol,k)=R_k(k)*vk_i(mol,k);
    end
    r_i(mol)=sum(R_k_i(mol,:));
end

R=sum(x_i.*r_i);
for mol=1:molecules
    row_i(mol)=r_i(mol)/R;
end
P_i=log(row_i);

%--combine the two--
for mol=1:molecules
    delta_i(mol)=row_i(mol)*(5*Q-1);
    cross_i(mol)=5*q_i(mol);
end


%--calculate combinatorial part--
for mol=1:molecules
    Ln_gamma_i_C(mol)=1+delta_i(mol)+P_i(mol)+cross_i(mol)*(omega_i(mol)-P_i(mol)-1);
end

%******************************************************************
%B)--calculate the residual contribution to activity coefficient

%calculate interaction energies between appropriate functional groups
for m=1:func_groups
    for n=1:func_groups
        alpha_m_n(m,n)=a_m_n(m,n)/T;
        tau_m_n(m,n)=exp(-1*alpha_m_n(m,n));
    end
end

Q_m_i=Q_k_i;

for mol=1:molecules
    for m=1:func_groups
        for n=1:func_groups
            t_i_m_n(mol,m,n)=Q_m_i(mol,m).*tau_m_n(m,n);
            u_i_m_n(mol,m,n)=Q_m_i(mol,m).*tau_m_n(n,m);
        end
    end
end

for mol=1:molecules
    for m=1:func_groups
            s_i_n(mol,m)=sum(t_i_m_n(mol,:,m));   %sum over all functional groups in molecule i
    end
end

for m=1:func_groups
        temp=s_i_n(:,m)';
        ss_n(m)=sum(x_i.*temp);  %sum over all molecules
end

for m=1:func_groups
    for n=1:func_groups
        temp=u_i_m_n(:,m,n)';
        uu_m_n(m,n)=sum(x_i.*temp);  %sum over all molecules
    end
end

for mol=1:molecules
    for n=1:func_groups
        x_weird_i_n(mol,n)=s_i_n(mol,n)/ss_n(n); 
        xx_weird_i_n(mol,n)=log(x_weird_i_n(mol,n));
    end
end

%--calculate the residual part--

for mol=1:molecules
    summation1=0;
    for n=1:func_groups
        summation2=0;
        brac1=0;
        brac2=0;
        for m=1:func_groups
            brac1=brac1+(u_i_m_n(mol,m,n))/s_i_n(mol,m);
            brac2=brac2+(uu_m_n(m,n))/ss_n(m);             
        end
        summation2=xx_weird_i_n(mol,n)-omega_i(mol)+brac1-brac2;
        summation1=summation1+Q_k_i(mol,n)*summation2;
    end
Ln_gamma_i_R(mol)=summation1;
end

%-------------------------------
%--final activity coefficient---

for mol=1:molecules
Ln_gamma_i(mol)=Ln_gamma_i_R(mol)+Ln_gamma_i_C(mol);
gamma_i(mol)=exp(Ln_gamma_i(mol));
end
%-------------------------------



