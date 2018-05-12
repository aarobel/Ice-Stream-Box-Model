%% Ice stream box model (Robel et al., JGR, 2013)
% This is the driver for the ice stream box model described in 
% Robel et al., JGR, 2013.
% In this script, model parameters and settings are specified and the model
% is integrated using ODE45 (though other MATLAB integrators may be used as
% well)
%
% See Ice_Stream_Box_Model_RHS for implementation of model equations.
% 
% Code written by Alex Robel
% Debugging help from Christian Schoof, Eli Tziperman, Eric DeGiuli and
% Elisa Mantelli
% Last updated for MATLAB 2018 by Colin Meyer, May 2018

%% Set parameters
p.year = 3600*24*365;   %seconds in a year
t_final = 1e4;          %total time of integration

p.L=500e3;              %ice stream length
p.W=40e3;               %ice stream width
p.n=3;                  %nye-glen law exponent
p.q_g = 0.07;           %geothermal heat flux
p.htill_init = 1;       %initial till thickness

p.T_s=23;               %ice surface temperature

p.rho_i = 917;          %density of glacial ice
p.L_f = 3.335e5;        %latent heat of fusion for ice
p.K_i = 2.1;            %thermal diffusivity of ice
p.A_f = 5e-25;          %nye-glen law rate factor
p.g = 9.81;             %acceleration due to gravity

p.e_c = 0.3;            %till consolidation void ratio
p.tau0 = 9.44e8;        %empirical till coefficient
p.c = 21.7;             %empirical till exponent
p.C_ice = 1.94e6;       %specific heat capacity of ice

p.A = p.L*p.W;              %ice stream area

p.eta_b = 10;           %basal ice layer thickness
p.h_t_min = 1e-3;       %minimum till thickness (prevents unfrozen till thickness from going to zero)

p.a = 0.1./p.year;      %accumulation rate

%% Set initial conditions and intergration time

p.ic=[700 .6 p.htill_init 0];   %Initial conditions
p.tspan=[0,p.year*t_final];     %time steps

%% Run Model

options = odeset('RelTol',1e-6,'AbsTol',1e-6);      %set ode integration settings
[time,T] = ode45(@(t,X) Ice_Stream_Box_Model_RHS(t,X,p),p.tspan,p.ic,options);  %integrate box model

%% Diagnostic
h = T(:,1);
e = T(:,2);
h_till = T(:,3);
T_b = T(:,4);
deltaT = p.T_s - T_b;               %difference between basal and surface ice temperature
e(e>=p.e_c)=p.e_c;               %make sure till void ratio doesn't go below threshold
h_till(h_till<=p.h_t_min) = p.h_t_min; h_till(h_till>=p.htill_init) = p.htill_init;  %make sure unfrozen till thickness stays within bounds
T_b(T_b<=0)=0;      %make sure basal ice temperature never goes above zero

tau_d = p.rho_i*p.g*(h.^2/p.L);        %calculate driving stress
tau_f = p.tau0*exp(-p.c*e);         %calculate basal shear stress

U = (p.A_f/256)*(p.W^(p.n+1))*(((tau_d-tau_f)./h).^p.n); U = max(U,0); %calculate centerline ice stream velocity
% % Make some plots % %
figure(1);set(1,'units','pixels','position',[0 0 1002 1202])
subplot(5,1,1);
plot(time./p.year,T(:,2),'k','linewidth',2);hold on;
ylabel('Void Ratio','fontsize',16);
set(gca,'fontsize',16)

subplot(5,1,2);
plot(time./p.year,T(:,1),'k','linewidth',2);hold on;
ylabel('Ice Thickness (m)','fontsize',16);
set(gca,'fontsize',16)

subplot(5,1,3);
plot(time./p.year,p.year.*U,'k','linewidth',2);hold on;
ylabel('Sliding Velocity (m/yr)','fontsize',16);
set(gca,'fontsize',16)

subplot(5,1,4);
plot(time./p.year,T(:,3),'k','linewidth',2);hold on;
ylabel('Till Thickness (m)','fontsize',16);
set(gca,'fontsize',16)

subplot(5,1,5);
plot(time./p.year,T(:,4),'k','linewidth',2);hold on;
ylabel('Basal Temp (K dep p-m)','fontsize',16);
set(gca,'fontsize',16)
