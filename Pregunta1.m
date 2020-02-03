%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%% Computational codes for problem set 3 and 4 %%%
%
%%% Overlapping generation models with stochatis dynamic process %%%
%
%           CRTISTIAN ORTIZ
%           PAUL PONCE
%           KEVIN ROJAS
%           SANTIAGO SANDOVAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

%%%%% SET OF PARAMETERES

beta = 1.011;
eta = 1.5;
alpha = 0.36;
delta = 0.06;
hlabor = 0.3;
rho = 0.96;
sigmazeta = 0.045;
sigmaz1 = 0.38;
mediaz1 = 1;
tau = 0.12;
normalmedia = 0;
normalsd = 0.045;
rh = 1.081241;

%%%%% AUXILIAR VARIABLES

households = 1000;
timeworking = 40;
timeresting = 20;
guessforcapital = 5;
gridcapitalsize = 60;
salaryinitial = 0.5;
rateinitial = 0.1;
flaghunt=1;

%%%%% COMPUTATIONAL PARAMETERS
tolED =  0.01;
maxiterED=100;

%%%%% FIXED PRODUCTIVITY AND SURVIVAL PROBABILITIES

alive = [   0.9991;	0.9991;	0.999;	0.999;	0.999;	0.999;	0.9991;	0.9995;	
    0.9985;	0.999;	0.9989;	0.9988;	0.9988;	0.9987; 0.9986;	0.9985;	0.9984;	
    0.9983;	0.9982;	0.9981; 0.9979;	0.9977;	0.9976;	0.9974;	0.9972;	0.997;	
    0.9968;	0.9965; 0.9963;	0.996; 0.9957;	0.9953;	0.9949;	0.9944; 0.9939;	
    0.9934;	0.9928;	0.9921;	0.9913;	0.9905; 0.9896;	0.9886;	0.9874;	0.9861;	
    0.9846;	0.9829; 0.9812;	0.9794;	0.9778;	0.9762; 0.9743;	0.9722; 0.9699;	
    0.9674;	0.9647;	0.9616;	0.958;	0.954;	0.9494;	0.9442; 0.9385];

productivity = [1;	1.0719;	1.1438;	1.2158;	1.2842;	1.3527;	1.4212;	1.4897;
    1.5582;	1.6267; 1.6952;	1.7217;	1.7438;	1.7748;	1.8014;	1.8279;	1.8545;
    1.881;	1.9075;	1.9341; 1.9606;	1.9623;	1.964;	1.9658;	1.9675;	1.9692;
    1.9709;	1.9726;	1.9743;	1.976; 1.9777;	1.97;	1.9623;	1.9546;	1.9469;
    1.9392;	1.9315;	1.9238;	1.9161;	1.9084; 1.9007];	



%%%%% Population vector

Population = ones(60,1);
Population(1,1) = 1000;

    
for i = 2:60    
Population(i,1) = Population(i-1,1)*alive(i,1);
end

%%%% Population floor

Populationfloor = floor(Population);
Muridos = 1000 - Populationfloor;

%%%%% Proportionf of population

Muedades = Population ./ sum(Population);
Muedades(61,1) = 0;


%%%%% Create the initial assignations 

pd = makedist('Lognormal','mu',mediaz1,'sigma',sigmaz1);
rng('default');  
ziniciales = random(pd,households,1);

ztotales = zeros(households, timeworking);

%%% Stochastic process for every househols based on AR(1)

pde = makedist('Normal','mu',normalmedia,'sigma',normalsd);
eproductivos = random(pde,households,timeworking);

for j = 1:households
    
        for i = 1:timeworking-1
    ztotales(j,1) = ziniciales(j,1);
    ztotales(j,i+1) = 2.7172^( (rho*log(ztotales(j,i))) + eproductivos(j, i+1) );  
        end 
    
end


%%%%% Condition of iteration in r for general equilibria

%%% First we define an r, ((then we need to iterate it))

r = rh;
critED=1;
iterED=0;

%%% Calcaulate the excess demand of capital

while ((critED>tolED) && (iterED<maxiterED))
iterED=iterED+1;
% Compute Demand of Capital

kdemand = ((r + delta)/alpha)^(1/(alpha - 1));

% Also compute the salary

salary = (1-alpha) * (kdemand)^alpha; 

% Calculation for (constant) b 

SumRetirement =  sum(Muedades(41:60,1));
SumWorkers = zeros(40,1);
for i = 1:40
    
    
    SumWorkers(i,1) = Muedades(i,1)*productivity(i,1)*ztotales(1,mean(i))*hlabor;
    
end

SumWorkersSum = sum(SumWorkers);

b = (tau*salary*SumWorkersSum)/SumRetirement;

% Compute Supply of Capital based on the household problem

% First, we identify the guess for the value of k ( based on stationary)

kstat = ( (1/alpha)*( (1/beta) - 1 + delta) )^(1/(alpha - 1));

% Creating a vector for consumption 

Consumes = zeros(timeresting + timeworking, households);
Goverments = zeros(timeresting + timeworking, households);
Savings = zeros(timeresting + timeworking, households);
Savings(60,:) = kstat;
Goverments(60,:) = (1 - alive(60,1))*Muedades(59,1)*kstat; 
Consumes(60,:) = (1 + r)*kstat + b; 
% Finding consumption and savings for RETIRED 

for j = 1:households

for i = 59:-1:41
Consumes(i,j) =  (((1 - (1 - alive(i+1,1))*Muedades(i,1))*Consumes(i+1,j)^(eta)) / ... 
(beta*alive(i+1,1)*(1+r)*(1 - (1-alive(i+2,1))*Muedades(i+1,1))))^(1/eta); %EULER OLD = EULER WORKER
Savings(i,j) = (Consumes(i,j) + Savings(i+1,j) - b - Goverments(i+1,j)) / (1+r); %Budget constraint Retired
Goverments(i,j) = sum((1-alive(i,1))*Muedades(i-1,1)*Savings(i,j)); %Transfer for accidental death 
end

end

% Finding consumption and savings for WORKER
for j =1:households 
for i = 40:-1:1
Consumes(i,j) =  (((1 - (1 - alive(i+1,1))*Muedades(i,1))*Consumes(i+1,j)^(eta)) / ... 
(beta*alive(i+1,1)*(1+r)*(1 - (1-alive(i+2,1))*Muedades(i+1,1))))^(1/eta); %EULER WORKER
Savings(i,j) = (Consumes(i,j) + Savings(i+1,j) - (1-tau)*(salary)*(hlabor)*...
(ztotales(j,i))*(productivity(i,1)) - Goverments(i+1,j)) / (1+r); %Budget constraint Worker
Goverments(i+1,j) = sum((1-alive(i+1,1))*Muedades(i,1)*Savings(i+1,j)); %Transfer for accidental death 
end
end
Savings(1,:) = 0;   %First k=0

% Correct table to include death people

for i = 1:60
    counter = Muridos(i,1);
    while (counter > 0)
       Savings(i, counter) = NaN; 
       counter = counter - 1;
    end
    
end

% Calculate mean(savings) by age
for i = 1:60
MeanAgeSavings(i,1) = nanmean(Savings(i,:));
end
% Calculating the excess of demand for k

MeanAgeSavings(61,1) = 0;
ksupply = sum((MeanAgeSavings).*(Muedades));

demandexcess = kdemand - ksupply;

    if flaghunt==1
        if demandexcess < 0                                      % still not bracketed true R
            rh=rh-0.01;
            r=rh;
        else                                            % brackted true R for first time
            rl=r;
            rh=rh+0.01;
            r=(rh+rl)/2;
            flaghunt=0;
        end
    else                                                % bisect
        if demandexcess < 0
            rh=r;
        else
            rl=r;
        end
        r=(rh+rl)/2;
    end
    disp('new interest rate')
    disp(r)
    critED=abs(demandexcess); 

end

disp(r);
disp(salary);
disp(b);

C = sum(Consumes(:));
K = nansum(Savings(:));
N = SumWorkersSum;

%%%%% Question 1

T = table(r, salary, b, C , K , N);
T

%%%%% Question 2

meanSavings = nanmean(Savings(25,:));
sdSavings = nanstd(Savings(25,:));
mean2stdPos = meanSavings+(2*sdSavings);
mean2stdNeg = meanSavings-(2*sdSavings);

figure
plot(Savings(25,:));
hold on
line([0,1000],[meanSavings,meanSavings],'Color','red')
line([0,1000],[mean2stdPos,mean2stdPos],'Color','black')
line([0,1000],[mean2stdNeg,mean2stdNeg],'Color','black')

figure
plot(Consumes);
figure
plot(Savings);

%%%%% Question 3

laborinc = (SumWorkers*salary)./ Muedades(1:40,1);

gini_cons = ginicoeff(Population,Consumes(:,1));
gini_k = ginicoeff(Population,MeanAgeSavings(1:60,1));
gini_laborincome = ginicoeff(Population(1:40,1),laborinc);

lorenzcurve(Population,Consumes(:,1));
lorenzcurve(Population,MeanAgeSavings(1:60,1));
lorenzcurve(Population(1:40,1),laborinc);
