%%% spin transport through MTJ
%% Ref. D. Datta et al, IEDM 10-548
%% coded by JG, May 15
%% reciprocacity of torques proven
clear all
%close all

global kBT1 kBT2

%% physical parameters
kBT=0.0259;  % thermal energy at 300K
q=1.6e-19;
hbar=1.055e-34;
m0=9.11e-31;
a0=2e-10;
t0=hbar^2/(2*m0*a0^2*q);  % TB parameter in eV


%% temperature of contacts
dT=0;  T0=300; 
T1=T0-dT/2;
T2=T0+dT/2;  % in K, temperature of contact2
kBT1=(T1/300)*kBT;
kBT2=(T2/300)*kBT;

%%% the bias condition
Vd0=0.0; % applied voltage
Vd_step=0.1; 
Nd_step=10;

%%% device parameters
Ns=10;
Nd=10;
Nox=10;  % the oxide thickness is a0*Nox
NI=1; % interface
Ntot=Ns+NI+Nox+NI+Nd;
Ub=1.4;
delt=2.15;
Ef=2.25;
mfm=0.73;   % effective mass of FM contacts
mox=0.2;    % effective mass of oxide
msi=0.19;
tfm=t0/mfm;
tox=t0/mox;
tsi=t0/msi; % tight binding parameter for Si
Efc=0.1;   % Ef-Ec in Si
Ec_si=Ef-Efc; % Ec edge in Si

%%% the FM contacts
sita=0;
Mu=[0 0 1];
mu=[sin(sita) 0 cos(sita)];
I2=eye(2);
sigx=[0 1; 1 0]; sigy=[0 -1i; 1i 0]; sigz=[1 0; 0 -1];

flag_spec=0; 

%%% the transverse wave vectors
Etmax=2; % in eV
ktmax=sqrt(2*mfm*m0*Etmax*q)/hbar;
Nkt=20;
kt_grid_n=1/Nkt*(1/2:Nkt-1/2); % normalized to ktmax
dkt_n=1/Nkt; % normalized to ktmax


%%%% set up the Hamiltonian matrix
HD=cell(Ntot,1);
AUD=cell(Ntot-1,1);
ALD=cell(Ntot-1,1);
    
for ii=1:Ntot-1
    if ii<=Ns  % FM
        AUD{ii}=tfm*I2;
    elseif ii<=Ns+Nox+1 % insulator
        AUD{ii}=tox*I2;
    else  % FM or Si contacts
        AUD{ii}=tsi*I2;
    end
end
for ii=1:Ntot-1
    ALD{ii}=AUD{ii}';
end

%%% compute current

Npm=Ns+NI+Nox+NI; % the left boundary of the magnet m
NpM=Ns;  % the right boundary of the magnet M

%%% initialization
Id=zeros(Nd_step+1,1);
torquem=zeros(Nd_step+1,3); % torque to m
torqueM=zeros(Nd_step+1,3); % torque to M
Ispos=zeros(Ntot-1,4,Nd_step+1);
for ii_vd=1:Nd_step+1
    Vd=Vd0+(ii_vd-1)*Vd_step
    Ef1=Ef+Vd/2;
    Ef2=Ef-Vd/2;
    Emax=max(Ef1,Ef2)+10*kBT;
    Emin=min(Ef1,Ef2)-10*kBT;
    dE=1e-3;
    E_grid=Emin:dE:Emax;
    for ii=1:Ntot
        if ii<=Ns
            HD{ii}=(2*tfm+Vd/2)*I2+delt/2*(I2-Mu(1)*sigx-Mu(2)*sigy-Mu(3)*sigz);
        elseif ii==Ns+NI
            HD{ii}=(tfm+tox+Vd/2+Ub/2)*I2+delt/2*(I2-Mu(1)*sigx-Mu(2)*sigy-Mu(3)*sigz);
        elseif ii>Ns+NI & ii<=Ns+NI+Nox
            HD{ii}=(2*tox+delt+Ub+Vd*(1/2-(ii-Ns-NI)/(Nox+1)))*I2;
        elseif ii==Ns+NI+Nox+NI
            HD{ii}=(tox+tsi-Vd/2+(Ub+delt+Ec_si)/2)*I2;
        else  % semiconductor 
            HD{ii}=(2*tsi-Vd/2+Ec_si)*I2;
        end
    end
    
    if flag_spec==1  % uniform energy grid for 1 kt mode kt=0
        for indE=1:length(E_grid)
            [JEnorm]=func_current(E_grid(indE), HD, AUD, ALD, Ef1, Ef2);
            JE(indE,1)=JEnorm(1,1);
        end
        Id(ii_vd)=q^2/(2*pi*hbar)*trapz(E_grid,JE); % in A, one mode current
    else
        
        %% Gaussian quadrature
       
        for ii_kt=1:Nkt
            kt=kt_grid_n(ii_kt)*ktmax; % in/m
            HD_kt=cell(Ntot,1);
            for ii=1:Ntot
                if (ii<=Ns) | (ii>Ns+NI+Nox+NI)
                    HD_kt{ii}=HD{ii}+(hbar^2*kt^2/(2*mfm*m0*q))*eye(2);
                elseif ii==Ns+NI | ii==Ns+NI+Nox+NI
                    HD_kt{ii}=HD{ii}+(hbar^2*kt^2/(4*m0*q))*(1/mfm+1/mox)*eye(2);
                else
                    HD_kt{ii}=HD{ii}+(hbar^2*kt^2/(2*mox*m0*q))*eye(2);
                end
            end
            [Inorm dum]=quadv(@func_current,min(E_grid),max(E_grid),1e-7,[],HD_kt, AUD, ALD, Ef1, Ef2);
            %% torque to m
            Is_m=Inorm(Npm,[3 4 2]); % the order of Inorm is c, z, x, y; Isl is[x y z]
            torquem(ii_vd,:)=torquem(ii_vd,:)-cross(mu,cross(mu,Is_m))*2*pi*kt_grid_n(ii_kt)*dkt_n;  % in eV
            %% torque to M
            Is_M=-Inorm(NpM,[3 4 2]);
            torqueM(ii_vd,:)=torqueM(ii_vd,:)-cross(Mu,cross(Mu,Is_M))*2*pi*kt_grid_n(ii_kt)*dkt_n;  % in eV
            %% position resolved charge and spin currents
            Ispos(:,:,ii_vd)=Ispos(:,:,ii_vd)+Inorm*2*pi*kt_grid_n(ii_kt)*dkt_n;   % in eV
            Id(ii_vd)=Id(ii_vd)+Inorm(1,1)*2*pi*kt_grid_n(ii_kt)*dkt_n;  % in eV
        end
        %% integral over transverse k in the polar coordinate.

        %% end of Gaussian quadrature approach
    end
    
end
G0=q^2/(2*pi*hbar);
Id=G0*(1/4/pi^2)*ktmax^2*Id; % in A/m^2
torquem=G0*(1/4/pi^2)*ktmax^2*torquem; % in A/m^2, torque to m, times hbar/2/q to get torque 
torqueM=G0*(1/4/pi^2)*ktmax^2*torqueM; % in A/m^2, torque to M
Ispos=G0*(1/4/pi^2)*ktmax^2*Ispos; % in A/m^2


%% visualization
figure(1)
xv=1:Ntot-1;
plot(xv,Ispos(:,1,Nd_step+1),'b-','linewidth',[2]); hold on;
plot(xv,Ispos(:,2,Nd_step+1),'k--','linewidth',[2]); hold on;
plot(xv,Ispos(:,3,Nd_step+1),'r--','linewidth',[2]); hold on;
plot(xv,Ispos(:,4,Nd_step+1),'g--','linewidth',[2]); hold on;
legend('charge','z','x','y')
set(gca,'linewidth',[2],'fontsize',[20]);
xlabel('x')
ylabel('I, I_s')

Vdv=Vd0:Vd_step:Vd0+Nd_step*Vd_step;
figure(2)
plot(Vdv,Id,'b--','linewidth',[2]); hold on;
set(gca,'linewidth',[2],'fontsize',[20]);
xlabel('V_D [V]')
ylabel('I [A/m^2]')

if flag_spec==1
    figure(21)
    df=1./(1+exp((E_grid-Ef1)./kBT))-1./(1+exp((E_grid-Ef2)./kBT));
    TE=JE./df'; % transmission
    plot(TE,E_grid,'linewidth',[2]);
    set(gca,'linewidth',[2],'fontsize',[20]);
    xlabel('T_E')
    ylabel('E [eV]')
end



    
    


    
    
    








