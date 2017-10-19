%% Run Phillips et al. model of the adenosine system and PVT responses, assuming quasistatic kinetics
% AJK Phillips 6/2015

function MP = chronic_model_final(pars,times)

global sleeptime

figs = 1;

sleeptimes = [0,3.7,5.5,6.7];
Ts = [4*24,16*24,16*24,16*24];
cols = {'k','r','g','b'};
%mm = zeros(1,length(sleeptimes));

trestart = find(diff(times)<0); % find places in the times vector where a new schedule begins
trestart = [0;trestart;length(times)]; % add end values

mod_pvt = [];

tic % start timer

%% Define model parameter values

Autip = 30;

%Kd1 = 0.014991;
Kd1 = pars(1);
Kd1 = max(Kd1,1);
Kd1 = min(Kd1,10);
%Kd2 = 5;
 Kd2 = pars(2);
 Kd2 = max(Kd2,100);
 Kd2 = min(Kd2,10000);

params{1} = Kd1; % KD1
params{2} = 0.1; % k1, binding rate constant
params{3} = params{2}*params{1}; % k2, unbinding rate constant 
params{4} = Autip/(Autip + params{1}); % B0R0, target receptor occupancy
params{5} = 300/(Kd2+300); % Beta = R2u/(K2+R2u)

Atot = Autip*(Autip+Kd1+600*(1-params{5}))/((Autip+Kd1)*(1-params{5}))

%Atot = [1:0.1:2000];
%jim = ((1-params{5}).*(Atot - 600 - Kd1/(1-params{5}) + sqrt((Atot + 600 + Kd1/(1-params{5})).^2 - 2400*Atot))-2*Autip)>0;
%Atot = Atot(find(diff(jim)>0))

pars(3) = max(pars(3),0);
pars(3) = min(pars(3),1);

params{6} = pars(3)*Atot; % lower asymptote for A_tot(t) during sleep (nM)
params{7} = (Atot - 0.65*params{6})/0.36; % upper asympotote for A_tot(t) during wake (nM)

params{8} = 18.18; % chiw, time constant for exponential saturation during wake (h)
params{9} = 4.20; % chis, time constant for exponential decay during sleep (h)
%params{10} = 1/200; % lambda, time constant for slow process (/h)
params{10} = 1/300;
params{11} = 60; % PVT sigmoid max value (lapses)
params{12} = pars(4); % PVT sigmoid half max argument
%params{12} = 542.5;
params{13} = pars(5); % PVT sigmoid sigma
%params{13} = 0.1;
params{14} = pars(6); % amplitude of circadian rhythm
%params{14} = 0.04;
%pars(7) = 6.2;
params{15} = pars(7); % circadian phase
%params{15} = 56;

% Kd1 = pars(1);
% Kd1 = max(Kd1,1);
% Kd1 = min(Kd1,10);
% params{1} = Kd1; % Specify Kd1
% Kd2 = pars(2); % Specify Kd2
% Kd2 = max(Kd2,100);
% Kd2 = min(Kd2,10000);
% params{5} = 300/(300+Kd2); 

% Atot = [1:1:2000];
% jim = ((1-params{5}).*(Atot - 600 - Kd1/(1-params{5}) + sqrt((Atot + 600 + Kd1/(1-params{5})).^2 - 2400*Atot))-2*Autip)>0;
% Atot = Atot(find(diff(jim)>0));
% 
% params{4} = Autip/(Autip + Kd1);
% 
% pars(3) = max(pars(3),0);
% pars(3) = min(pars(3),1);
% params{6} = pars(3)*Atot;
% params{7} = (Atot - 0.65*params{6})/0.36;
% params{11} = pars(4);
% params{12} = pars(5);
% params{13} = pars(6);
% params{14} = pars(7);
% params{15} = pars(8);
% params{10} = pars(10);

params



%% Simulate baseline

T = 80*24; % length of baseline run
sleeptime = -1; % number of hours of sleep per night (baseline setting)

%% Initial conditions

% estimate sensible initial conditions
%Atoti = params{6} + 0.3*(params{7}-params{6});
Atoti = params{6} + 0.6237*(params{7}-params{6});
Amean = params{6} + 0.302*(params{7}-params{6});
%Rtoti = (Amean-params{5})/params{4} - params{1}/(1-params{4});
Rtoti = (Amean)/params{4} - params{1}/((1-params{4})*(1-params{5}));
%Rtoti = 80;
Abi = Rtoti*params{4}*1.05;
%xi = [Abi;Atoti;Rtoti];
xi = [Atoti;Rtoti];


%% Solve differential equations using chronic_de function (with m-factor)

mfactor = 20; % m-factor to make slow time scale faster, so as to allow a shorter first phase of the baseline run
params{10} = params{10}*mfactor;

[t,Y] = chronic_de_final(xi,params,T);


Atoti = Y(1,end);
Rtoti = Y(2,end);

%% Solve again (without m-factor)

T = 80*24; % length of baseline run
params{10} = params{10}/mfactor;
xi = [Atoti;Rtoti];

[t,Y] = chronic_de_final(xi,params,T);


Atoti = Y(1,end);
Rtoti = Y(2,end);

%% Simulate experimental conditions

for ii = 1:length(sleeptimes),
    
    sleeptime = sleeptimes(ii);
    T = Ts(ii);
    xi = [Atoti;Rtoti];

    %% Solve differential equations using chronic_de function

[t,Y] = chronic_de_final(xi,params,T);
    
%% Define model variables

%Ab = Y(1,:); % concentration of bound adenosine molecules
%Rb = Ab; % concentration of bound adenosine receptors
Atot = Y(1,:); % concentration of adenosine molecules
Rtot = Y(2,:); % concentration of adenosine receptors
%Ab =  0.5*((Atot+Rtot+params{1}-params{5}) - sqrt((Atot+Rtot+params{1}-params{5}).^2 - 4*(Atot-params{5}).*Rtot));
Ab =  0.5*((Atot+Rtot+params{1}/(1-params{5})) - sqrt((Atot+Rtot+params{1}/(1-params{5})).^2 - 4*Atot.*Rtot));
Rb = Ab;
A2 = (Atot-Ab)*params{5};
%Au = Atot - Ab - params{5}; % concentration of unbound adenosine molecules 
Au = Atot - Ab - A2;
Ru = Rtot - Rb; % concentration of unbound adenosine receptors

%% Define sleep drive and PVT performance

D = Rb + params{14}*cos((t-params{15})*2*pi/24); % sleep drive is homeostatic drive (Rb) + sinusoid (process C)
P = params{11}./(1 + exp((params{12} - D)/params{13})); %+ 0.1/((R(end)-100)*(1000-R(end)))];

%mm(ii) = Rb(end)-Rb(1);

tvector = times(trestart(ii)+1:trestart(ii+1))+6; % Add 6 to bring the times in line with our clock time (sleep from midnight to 8am)

for k = 1:length(tvector),
    indo = find(abs(t-tvector(k))==min(abs(t-tvector(k)))); % index for this time point
    %mod_pvt = [mod_pvt; PVTa*(B(indo).^PVTc) + PVTb + Camp*cos((t(indo)-Cphase)*2*pi/24)];
    mod_pvt = [mod_pvt; P(indo)]; %+ 0.1/((R(end)-100)*(1000-R(end)))];
end

%% Plot results

if figs ==1,

figure(54)
plot(t./24,Rb,cols{ii})
hold on
xlabel('Time (days)')
ylabel('Rb (nM)')
else
end

end
    
hold off

%ymetrics = [mm(2)/mm(1),mean(Au)/20,mean(Rtot)/600];

%mean(Au)
%mean(Rtot)

if figs==1,

figure(1)
subplot(4,1,1)
plot(t,Au,'b')
xlabel('Time (h)')
ylabel('Au (nM)')
subplot(4,1,2)
plot(t,Atot,'b')
xlabel('Time (h)')
ylabel('Atot (nM)')
subplot(4,1,3)
plot(t,Rtot,'b')
xlabel('Time (h)')
ylabel('Rtot (nM)')
subplot(4,1,4)
plot(t,Rb,'r')
xlabel('Time (h)')
ylabel('Rb (nM)')

figure(2)
plot(t,P)
xlabel('Time (h)')
ylabel('PVT lapses')

figure(3)
plot(t,state_final(t))
ylim([-0.1,1.1])

else
end

toc

MP = mod_pvt;

end

