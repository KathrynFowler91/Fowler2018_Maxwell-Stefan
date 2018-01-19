function [u_save, RN] = implicit_maxwell_stefan03(varargin)

%--------------------------------------------------------------------------
%
% A function to solve the diffusion equation on a spherical core shell 
% model with neumann boundary conditions. Fick's second law is solved using
% the backward Euler method of finite differences and the non-ideal effects 
% of diffusion are included using the Maxwell-Stefan framework, relating 
% diffusion flux to gradient in activity coefficients. The UNIFAC group
% contribution model is used to estimate the activity coefficients.
%
% Call the model by typing [u_save, RN] = implicit_maxwell_stefan03().
% Input parameters maybe changed below in the 'setting up parameters 
% section' or using 'implicit_maxwell_stefan_runscript'.
%
% Outputs variables are:
% u_save = concentration (or molar density) for a given shell, timestep and 
%          component.
% RN = inital shell boundaries as determined from the number of shells and
%      size of grid.
%
% Other functions required to run the script include:
%   outer_shell_equilibration.m
%   outer_shell_redistribution.m
%   UFC_datamain.m
%   UFC_dataqi_v1.m
%   UFC_interact_params_v1.m
%   UNIFAC_gamma.m
%
% Details of the calculations contained within this script can be found: 
% https://www.atmos-chem-phys-discuss.net/acp-2017-424/#discussion
%
% Copyright (C) 2018  Kathryn Fowler

warning off;

%--------------------------------------------------------------------------
% setting up parameters

fickian = false;            % true==fickian / false==MS
diffusion = 'darken';       % 'darken' / 2 component 'vignes'

n_components = 2;
dt = 10.^(-9);             % time-step
ip = 1e2;                   % number of grid points
ntm = 1e3;                  % number of time steps

DSelf = 1e-18;
aerosol = 'Monocarboxylic'; % see MS variable switch (only for 2 components)
CNum = 4;

% initial particle conditions - to change initial concentration profile see below
Xw_init = 0.1;              % initial water mole fraction
Xw_shell = zeros(ntm,1);
Xw_shell(:) = 0.8;
% Xw_shell(1:0.3e3) = 0.1;
% Xw_shell(0.3e3+1:end) = 0.9; % step up to equilibrium RH of the outer shell 
 

r0 = 0;                     % lowest radius edge
rN = 1.2e-7;                % upper radius edge
R = 1e-7;                   % initial radius of aerosol

T = 293.15;                 % temperature (Kelvin)

%--------------------------------------------------------------------------
% user input - to use in batch function

if nargin > 0 fickian=varargin{1}; end
if nargin > 1 diffusion=varargin{2}; end
if nargin > 2 n_components=varargin{3}; end
if nargin > 3 dt=varargin{4}; end
if nargin > 4 ip=varargin{5}; end
if nargin > 5 ntm=varargin{6}; end
if nargin > 6 rN=varargin{7}; end
if nargin > 7 R=varargin{8}; end
if nargin > 8 aerosol=varargin{9}; end
if nargin > 9 CNum=varargin{10}; end
if nargin > 10 DSelf=varargin{11}; end
if nargin > 11 Xw_init=varargin{12}; end
if nargin > 12 Xw_shell=varargin{13}; end
if nargin > 13 T=varargin{14}; end

%--------------------------------------------------------------------------
% initialising matrices

n = n_components;
Dself=[2e-9 DSelf];         % m^2/s
D=zeros(ip,n);
RN=linspace(r0,rN,ip+1);    % radius bin-edges
dR=diff(RN);                % distance between grid faces
RM=RN(1:end-1)+dR./2;
dR2=[diff(RM) 100];         % distance between grid cells

ind=find(RN<R);

RN(ind(end)+1)=R;           % radius bin-edges
dR=diff(RN);                % distance between grid faces        
RM=RN(1:end-1)+dR./2;
dR2=[diff(RM) 100];         % distance between grid cells
vol=4.*pi./3.*(RN(2:end).^3-RN(1:end-1).^3);

ni=zeros(ip,n);             % molefraction
mi=zeros(ip,n);             % massfraction
u=zeros(ip,n);              % solution - concentration
A=sparse(ip,ip);            % matrix A

u_save=zeros(ip,ntm+1,n);

B=zeros(ip,1);              % solution vector
act=zeros(ip,n);            % activity
act05=zeros(ip,n);          % activity
x05=zeros(ip,n);            % molefraction

%-------------------------------------------------------------------------
% maxwell-stefan variables

% DT (Dave Topping) notes
% Given we have not got any ionic compounds im using the UNIFAC model with
% the Hansel et al parameters updated with new EDB variants where
% appropriate.

% Set the UNIFAc parameters here:

m_weights=[400.*ones(n,1)]./1000;
m_weights(1)=18./1000;
rho_s=[1500.*ones(n,1)];    % density of pure component
rho_s(1)=1000;
J=zeros(ip,n);
lhs=zeros(ip,n+1);
lhs=zeros(ip,n);
v=J;

molecules=n;

% Water
molecule_group_flag(1,1)=17; %h2o
sub_groups(1)=1;
molecule_group_stoich(1,1)=1;

% 2nd component switch
switch aerosol
    case 'Monocarboxylic'
        if CNum==1 % Formic acid
            molecule_group_flag(2,1)=43; %HCOOH
            sub_groups(2)=1;
            molecule_group_stoich(2,1)=1;
        elseif CNum==2 % Acetic acid
            molecule_group_flag(2,1)=1;  %ch3
            molecule_group_flag(2,2)=43; %cooh
            sub_groups(2)=2;
            molecule_group_stoich(2,1)=1;
            molecule_group_stoich(2,2)=1;
        else % +2 Carbon chains  monocarboxylic acid
            molecule_group_flag(2,1)=1;  %ch3
            molecule_group_flag(2,2)=2; %ch2
            molecule_group_flag(2,3)=43; %cooh
            sub_groups(2)=3;
            molecule_group_stoich(2,1)=1;
            molecule_group_stoich(2,2)=CNum-2;
            molecule_group_stoich(2,3)=1;
        end
    case 'Dicarboxylic'
        if CNum<=1
            print('CNum too low')
        elseif CNum==2 % Oxalic acid
            molecule_group_flag(2,1)=43; %cooh        
            sub_groups(2)=1;
            molecule_group_stoich(2,1)=CNum;
        else
            molecule_group_flag(2,1)=43; %cooh
            molecule_group_flag(2,2)=2; %ch2
            molecule_group_flag(2,3)=43; %cooh 
            sub_groups(2)=3;
            molecule_group_stoich(2,1)=1;
            molecule_group_stoich(2,2)=CNum-2;
            molecule_group_stoich(2,3)=1;
        end

    case 'Ethanoic'
        % acetic acid (ethanoic) CH3COOH
        molecule_group_flag(2,1)=1;  %ch3
        molecule_group_flag(2,2)=43; %cooh
        sub_groups(2)=2;
        molecule_group_stoich(2,1)=1;
        molecule_group_stoich(2,2)=1;
    case 'Propanoic'
        % propionic acid (propanoic) C3H7COOH
        molecule_group_flag(2,1)=1;  %ch3
        molecule_group_flag(2,2)=2; %ch2
        molecule_group_flag(2,3)=43; %cooh
        sub_groups(2)=3;
        molecule_group_stoich(2,1)=1;
        molecule_group_stoich(2,2)=1;
        molecule_group_stoich(2,3)=1;
    case 'Butanoic'
        % butyric acid (butanoic) C3H7COOH
        molecule_group_flag(2,1)=1;  %ch3
        molecule_group_flag(2,2)=2; %ch2
        molecule_group_flag(2,3)=43; %cooh
        sub_groups(2)=3;
        molecule_group_stoich(2,1)=1;
        molecule_group_stoich(2,2)=2;
        molecule_group_stoich(2,3)=1;
    case 'Pentanoic'
        % valeric acid (pentanoic) C4H9COOH
        molecule_group_flag(2,1)=1;  %ch3
        molecule_group_flag(2,2)=2; %ch2
        molecule_group_flag(2,3)=43; %cooh
        sub_groups(2)=3;
        molecule_group_stoich(2,1)=1;
        molecule_group_stoich(2,2)=3;
        molecule_group_stoich(2,3)=1;
    case 'Hexanoic'
        % caprioc acid (hexanoic) C4H9COOH
        molecule_group_flag(2,1)=1;  %ch3
        molecule_group_flag(2,2)=2; %ch2
        molecule_group_flag(2,3)=43; %cooh
        sub_groups(2)=3;
        molecule_group_stoich(2,1)=1;
        molecule_group_stoich(2,2)=4;
        molecule_group_stoich(2,3)=1;
    case 'Sucrose'
        % sucrose C12H22O11
        molecule_group_flag(2,1)=18;  %ACOH
        molecule_group_flag(2,2)=15; % OH
        molecule_group_flag(2,3)=13; % ACCH2
        molecule_group_flag(2,4)=10; % ACH
        molecule_group_flag(2,5)=114; % CHOCH
        molecule_group_flag(2,6)=115; % CHOC
        sub_groups(2)=4;
        molecule_group_stoich(2,1)=5; % checked
        molecule_group_stoich(2,2)=3; % checked
        molecule_group_stoich(2,3)=3; % checked
        molecule_group_stoich(2,4)=8;
        molecule_group_stoich(2,5)=1; % checked
        molecule_group_stoich(2,6)=1; % checked
    otherwise
        disp('no aerosol found')
end

%-------------------------------------------------------------------------

% populate gamma vs molefraction
N=1e3; % the accuracy of the MS model when diffusion coefficient is close to zero
clear mole_fraction;
mole_fraction(1,:)=linspace(0,1.,N);
for j=2:n
    mole_fraction(j,:)=1-mole_fraction(1,:);
end
gamma=zeros(N,n);
gamma1=zeros(ip,n);
gamma2=zeros(ip,n);
for i=1:N
    mf=mole_fraction(:,i);
    mf(n)=max(mf(n),0.0);
    gamma(i,:)=UNIFAC_gamma(molecules,mf',sub_groups,molecule_group_flag,molecule_group_stoich,T); 
end
if fickian == true
    gamma(:)=1;
end

%-------------------------------------------------------------------------
% Set up initial concentration profile

% initial particle in equilibrium with water mole fraction Xw
ind2 = find(RM<=R);
ni(ind2,1) = Xw_init;

% Example alternative concentration profile
% ind2=find(RM<=R);
% ni(ind2,1)=0;
% ind=find(RM<0);
% % ind=find(RM>0.5e-7 & RM<0.8e-7);
% ni(ind,1)=1.0;

for j=2:n
    ni(ind2,j)=1-ni(ind2,1)./(n-1);
    ni(ind,j)=1-ni(ind,1)./(n-1);
end

rhot=sum(ni(ind2,:).*repmat(m_weights',[length(ind2) 1]),2) ...
    ./sum(ni(ind2,:).*repmat(m_weights'./rho_s',[length(ind2) 1]),2);
mt=rhot.*vol(ind2)';
nt=mt./sum(ni(ind2,:).*repmat(m_weights',[length(ind2) 1]),2);
u(ind2,:)=ni(ind2,:).*repmat(nt./vol(ind2)',[1 n]);
r=dt./dR2(1); % used in implicit scheme

%-------------------------------------------------------------------------
% Time loop

% save
u_save(ind2,1,:)=u(ind2,:); % concentration

for i=1:ntm % solve over ntm time-steps
    
    % Solve for Maxwell-Stefan fluxes here:
    % calculate the fluxes
    u05=[0.5.*(u(ind2(2:end),:)+u(ind2(1:end-1),:)); 0.5.*(u(ind2(end),:)+u(ind2(end-1),:))]; % averages the difference between the two shells
    u05(:)=max(u05(:),1e-30);
    
    for j=1:n
        gamma1(ind2,j)=min(interp1(mole_fraction(j,:)',gamma(:,j),u05(ind2,j)./sum(u05(ind2,:),2),'linear'),100000);
        gamma2(ind2,j)=min(interp1(mole_fraction(j,:)',gamma(:,j),u(ind2,j)./sum(u(ind2,:),2),'linear'),100000);
    end

    
    small=0.; 1e-60;
    act(ind2,:)=gamma2(ind2,:).*u(ind2,:)./(repmat(sum(u(ind2,:),2),[1 n]));
    
    lhs(ind2(1:end-1),1:n)=(act(ind2(2:end),:)-act(ind2(1:end-1),:)); 
    lhs(ind2(end),1:n)=(act(ind2(end),:)-act(ind2(end-1),:)); % forward difference

    act05(ind2,:)=gamma1(ind2,:).*u05(ind2,:)./(repmat(sum(u05(ind2,:),2),[1 n])+small);       
    x05(ind2,:)=u05(ind2,:)./(repmat(sum(u05(ind2,:),2),[1 n])+small);    
        
    ind=find(abs(lhs(ind2,1))>0 & abs(lhs(ind2,2))>0);
    J(:)=0.;
    
    if fickian==true
       lhs(ind2(1:end-1),1:n)= (u(ind2(2:end),:)-u(ind2(1:end-1),:))./repmat(dR2(ind2(1:end-1))',[1 n]);
       lhs(ind2(end),1:n)= lhs(ind2(end-1),1:n);
    end
    
    %----------------------------------------------------------------------
    % Diffusion Coefficients
    
    for k=ind2(ind) % for each point on the grid
                
        switch diffusion
            case 'darken'
                % Darken equation (equation 4.12 and 4.13 in http://homepage.tudelft.nl/v9k6y/thesis-xin-liu.pdf) :  
                Dij=zeros(n,n);
                for ii=1:n
                    for jj=1:n
                        Dij(ii,jj)=Dself(ii)*Dself(jj)*sum(x05(k,:)'./Dself(:));
                    end
                end
                Dij2=(1-eye(n,n))./Dij;
            case 'vignes'
                if n==2
                    % Viscous liquids - Vignes equation
                    Dij=zeros(n,n);
                    for ii=1:n
                        for jj=1:n
                            Dij(ii,jj)=Dself(ii).^x05(k,ii).*Dself(jj).^x05(k,jj);
                        end
                    end
                else
                    print('the vignes equation is not multicomponent');
                end
                Dij2=(1-eye(n,n))./Dij;
        end
        
        if fickian==false
            % 2. (ii) now find the vector fluxes (create the Matrix A, which when
            % multiplied by the vector of 1-d fluxes gives the rhs of Maxwell-Stefan  
            m1=0.0;
            m2=0.0;

            % First the off-diagonal elements (easy):
            B=(1-eye(n,n)).* ...
                (1./Dij-rho_s(n)./(m_weights(n).*repmat(Dij(:,n),[1 n]) ).* ...
                repmat(m_weights'./rho_s',[n 1])).* ...
                repmat((max(x05(k,:),m1))',[1 n]);
            B(:,n)=0;
            % Now the diagonal elements:
            mat1=sum(Dij2.*repmat(max(x05(k,:),m1),[n 1]),2)';
            mat1=mat1+x05(k,:).*(rho_s(n).*m_weights./(m_weights(n).*rho_s.*Dij(:,n)))';
            B=B-diag(mat1,0);
            B=B./(sum(max(u05(k,:),m2),2)).*repmat(gamma1(k,:)',[1 n]).*dR2(ind2(k));
            sol=B\lhs(k,1:n)';
            if(min(abs(sol.*m_weights./rho_s))<1e-40)
                sol(:)=0;
            end
            
            J(k,:)=sol(1:n);
        
        else
            J(k,1:n)=-Dij([n:n-1:n.^2-(n-1)]).*lhs(k,1:n);
        end
    end
    
    %--------------------------------------------------------------------------
    
    % Fick's first law to find diffusion coefficients
    D(ind2(1:end-1),:)=-J(ind2(1:end-1),:).*repmat(dR2(ind2(1:end-1))',[1 n]).*1./( u(ind2(2:end),:)-u(ind2(1:end-1),:) );

    ind=find(isnan(D(:,1)) |isinf(D(:,1)));
    ind1a=find(~isnan(D(:,1)) &~isinf(D(:,1)));
    D(ind,1)=interp1(ind1a,D(ind1a,1),ind,'linear','extrap');

    D(ind,:)=0.;
    ind=find(isnan(D(:,2)) |isinf(D(:,2)));
    ind1a=find(~isnan(D(:,2)) &~isinf(D(:,2)));
    D(ind,2)=interp1(ind1a,D(ind1a,2),ind,'linear','extrap');
    D(ind,:)=0.;

    indd=find(D(:)>1e-9);
    D(indd)=1e-9;    
    indd=find(D(:)<-1e-9);
    D(indd)=-1e-9;   
    indd=find(D(:,1)<0 | D(:,2)<0);
    D(indd,:)=0;   
    

    D2=D;
    for j=1:n
        % old mass
        mold=sum(u(ind2,j).*vol(ind2)');
        
        % tridiagonal terms:
        Ac=(1+...
            dt./(dR(ind2)'.*dR2(ind2)'.*RM(ind2)'.^2).*RN([ind2(2:end) ind2(end)+1])'.^2.*D(ind2,j)+...
            dt./(dR(ind2)'.*dR2([ind2(1) ind2(1:end-1)])'.*RM(ind2)'.^2).*RN(ind2(1:end))'.^2.*D([ind2(1) ind2(1:end-1)],j));
        
        Au=-dt./(dR(ind2(1:end-1))'.*dR2(ind2(1:end-1))'.*RM(ind2(1:end-1))'.^2).*(RN(ind2(2:end))'.^2.*D(ind2(1:end-1),j));
        Al=-dt./(dR(ind2(2:end))'.*dR2(ind2(1:end-1))'.*RM(ind2(2:end))'.^2).*(RN([ind2(2:end)])'.^2.*D(ind2(1:end-1),j));

        AuN=-dt./(dR(ind2(end))'.*dR2(ind2(end))'.*RM(ind2(end))'.^2).*(RN(ind2(end)+1)'.^2.*D(ind2(end),j));
        Al1=-dt./(dR(ind2(1))'.*dR2(ind2(1))'.*RM(ind2(1))'.^2).*(RN(ind2(1))'.^2.*D(ind2(1),j));
        
        % Neumann BCs:
        Au(1)=Au(1)+Al1;
        Al(end)=Al(end)+AuN;
 
        % Matrix A:
        A(1:ip+1:ip.*ind2(end))=Ac;
        A(2:ip+1:ip.*(ind2(end)-1))=Al;
        A(ip+1:ip+1:ip.*ind2(end))=Au;
        
        % Solution vector:
        B=u(ind2,j);
        u(ind2,j)=A(ind2,ind2)\B(:);
        u(:,j)=max(u(:,j),1e-30);
        mnew=sum(u(ind2,j).*vol(ind2)');   
        u(ind2,j)=u(ind2,j).*mold./(mnew+1e-60);

    end

%-------------------------------------------------------------------------
% Particle growth or shrinkage 

    % calculates amount of water to add to equilibriate outer shell
    deltaV = outer_shell_equilibration(u,ind2,vol,Xw_shell(i),n,rho_s,m_weights,r0,rN,ip);
    deltaR = (3.*deltaV./4./pi+R.^3).^(1./3)-R;

    if deltaV > 0  % case 1: put water on particle, taking care that the last previous layer fills up as particle grows

        R=R+deltaR;
        t=u(ind2,:).*repmat(vol(ind2)',[1 n]); % total moles in each layer
        RN=linspace(r0,rN,ip+1); % radius bin-edges
        vollayer=4./3.*pi.*min(R.^3-RN(ind2(end)).^3,RN(ind2(end)+1).^3-RN(ind2(end)).^3);
        voladdlayer=vollayer-sum(u(ind2(end),:).*vol(ind2(end)).*m_weights'./rho_s');
        molewateradd=voladdlayer.*rho_s(1)./m_weights(1);
        
        t(ind2(end),1)=t(ind2(end),1)+molewateradd;
        
        ind2a=find(RN<R);
        ind2b=ind2(end)+1:ind2a(end);

        u(ind2b,1)=(1-1e-30).*rho_s(1)./m_weights(1);    

        for j=2:n
            u(ind2b,j)=(1e-30./(n-1)).*rho_s(j)./m_weights(j);
        end

        RN(ind2a(end)+1)=R; % radius bin-edges
        dR=diff(RN);  % distance between grid faces        
        RM=RN(1:end-1)+dR./2;
        dR2=[diff(RM) 100]; % distance between grid cells
        vol=4.*pi./3.*(RN(2:end).^3-RN(1:end-1).^3);

        u(ind2,:)=t(ind2,:)./repmat(vol(ind2)',[1 n]);
     
        ind2=ind2a;
              
    elseif deltaV<0 % case 2: take water off particle, take the water off first, then add up the other stuff and reset particle radius

        deltaV=-deltaV; % defined above
        v=0;
        t=u(ind2,:).*repmat(vol(ind2)',[1 n]); % total moles in each layer
        volwtot=sum(t(ind2,1).*m_weights(1)./rho_s(1));
        
        RN=linspace(r0,rN,ip+1); % radius bin-edges
        dR=diff(RN);  % distance between grid faces
        RM=RN(1:end-1)+dR./2;
        dR2=[diff(RM) 100]; % distance between grid cells
        
        deltaV=min(volwtot.*0.99,deltaV);
        
        if volwtot > 0. % only do this if there is any water on the particle

            % calculate how many layers we need to rid the water of
            k=ind2(end);
            while v< deltaV
                v=v+t(ind2(k),1).*m_weights(1)./rho_s(1);
                k=k-1;
            end
            istart=k+1;
            
            % calculate new radius, not including the other diffusing species
            volw_first_layer= v-deltaV;
            % set water volume
            t(ind2(istart),1)=volw_first_layer.*rho_s(1)./m_weights(1);
            % set water volume in outer layers to zero:
            t(ind2(istart)+1:ind2(end),1)=0.;            
            % volume of other to add
            volo_outer=...
                sum(t(ind2(istart:end),2:n).*m_weights(2:n)./rho_s(2:n),1);
            % now fill up shells with the other material
            k=istart;
            % then the next empty shells:
            while volo_outer > 0e-24                
                deltaVo=min(...
                    4.*pi./3.*(RN(ind2(k)+1).^3-RN(ind2(k)).^3)-t(ind2(k),1).*m_weights(1)./rho_s(1), volo_outer);

                t(ind2(k),2)=rho_s(2)./m_weights(2).*deltaVo;
                volo_outer=volo_outer-deltaVo;
                k=k+1;
            end
            iend=k-1;
            Rtemp=(3.*(sum(t(ind2(iend),1:n).*m_weights(1:n)'./rho_s(1:n)'))./(4.*pi)+...
                RN(ind2(iend)).^3 ).^(1./3);
            
            % now set indices
            ind2a=find(RN<R);
            ind2b=ind2a(end)+1:ind2(end);
            t(ind2b,:)=0.;
            
            RN=linspace(r0,rN,ip+1); % radius bin-edges
            dR=diff(RN);  % distance between grid faces
            RM=RN(1:end-1)+dR./2;
            dR2=[diff(RM) 100]; % distance between grid cells

            RN(ind2a(end)+1)=R; % radius bin-edges
            dR=diff(RN);  % distance between grid faces        
            RM=RN(1:end-1)+dR./2;
            dR2=[diff(RM) 100]; % distance between grid cells
            vol=4.*pi./3.*(RN(2:end).^3-RN(1:end-1).^3);

            u(ind2a,:)=t(ind2a,:)./repmat(vol(ind2a)',[1 n]);

            u(ind2a,:)=max(u(ind2a,:),1e-30);
            ind2=ind2a;
        end

    end

    % if number of shells increases or decreases, need to redistribute outer 2 layers
    ind3(i+1) = ind2(end);
    t = outer_shell_redistribution(u,vol,ind2,ind3(i),n,Xw_init,Xw_shell(i),m_weights,rho_s);
    u(ind2,:) = t(ind2,:)./repmat(vol(ind2)',[1 n]);
    
    % save concentration profile for current timestep
    u_save(ind2,1+i,:)=u(ind2,:);
    
end
