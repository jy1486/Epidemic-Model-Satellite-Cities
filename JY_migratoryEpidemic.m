%Jason Yang 
clf;
clear;
close all;

dt=0.5; % time interval
tmax = 100; % 100 days
clockmax=ceil(tmax/dt);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ncity = 4; %number of separate cities. Assume 1 is central city.

%Initially I wanted to see what if I varied parameters of all cities, so I
%decided not to use the non-dimensionalised simplified equation. Hard to tell what is going on
% and not very useful so I commented out this part below
%
% a = [4 , 6, 8, 10]/20; %infectivity [1/time ]
% b = [2.86, 2, 2 , 2]/20; %recovery [1/time]
% d = [1.25, 1.55 ,1.55 , 1.55]/50; %death rate constant [1/time]
% %SIR = zeros(3, Ncity); %each row is a city, col is SIR  #SIR x #CITY   
% SIR = [550, 250, 250 ,250;
%          10, 5 , 5    ,5;
%          0, 0 , 0    ,0];
% D = [0 ,0 ,0 ,0];
% %migration rate constant as a matrix
% migMat = [ 20, 35 ,35 ,35 ; 
%            10, 85 ,85 ,85; 
%            20, 35 ,35 ,35 ]/200;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Assume satellite cities have the same parameters
%modify parameters in this section of code
%%%satellite city parameters
a1 = 8/30; %infectivity 
b1 = 2/30; %recovery
d1 = 1.25/50; %death rate
SIR1 = [ 250; 5; 0]; %Starting population size, 250 susceptible, 5 infected, 0 recovered
D1 = 0; %Dead group
migMat1 = [35; 55; 35]/200; %migration rate constants, eg. 35/200 is migration rate of susceptible
%55/200 is migration rate of infected, 35/200 is migration rate of
%recovered

%%%central city parameters
a = 4.20/30 ;
b = 2.86/30;
d = 1.25/50;
SIR = [ 550; 10;0];
D = 0
migMat = [20; 10; 20]/200;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%migMat1 = zeros(3,1);
%migMat = zeros(3,1);




%%%%%%%%%%%Forming the matrices using variables above 
a = [a, ones(1,Ncity - 1) * a1];
b = [b, ones(1,Ncity - 1) * b1];
d = [d, ones(1,Ncity - 1) * d1];
SIR = [SIR, repmat(SIR1,1,Ncity - 1) ];
D = [0, repmat(D1,1,Ncity - 1)];
migMat = [migMat, repmat(migMat1,1,Ncity -1)];

       

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%total number of alive people in each city as a vector
N = sum(SIR);

     
% city 1 is the central city
%first variable in addition creates a matrix
% tempMat[0 ,1 , 1...   diag[-n+1, 1, 1 ...     
%         1, 0 , 0 ...         0 , -1, 0 ...   
%         1, 0 , 0 ...]        0 , 0, -1 ...]
tempMat = cat(2,[0,ones(1,Ncity-1)]',cat(1,ones(1,Ncity-1),zeros(Ncity-1,Ncity-1)));
flowMat = tempMat + diag([-Ncity+1 ,repmat(-1,1,Ncity-1)]); %city x city
%flow = [-2,  1 , 1
%         1, -1 , 0
%         1,  0 ,-1] 

%note use -Ncity + 1 because city doesn't migrate to itself and diag is -1
%for each neighbor migration to center city, first row fill with 1 for migration
%into neighbor city

%during time step TrueMigMat = migMat .* SIR * flowMat ;
%#SIR~3 x #city * #city x city .* #SIR~3 x #city 

%temp variables used during differential equation calculations
S2I = zeros(1,Ncity);
I2R = zeros(1,Ncity);
I2D = zeros(1,Ncity);

S2S = zeros(1,Ncity); %Net migration of suseptible
R2R = zeros(1,Ncity); %Net migration of recovered
I2I = zeros(1,Ncity); %Net migration of infected

%variables used for graphing
Ssave=zeros(Ncity,clockmax+1); 
Isave=zeros(Ncity,clockmax+1);
Rsave=zeros(Ncity,clockmax+1);
Nsave=zeros(Ncity,clockmax+1);

Dsave =zeros(Ncity,clockmax+1);
tsave=zeros(1,clockmax+1);

%setting up initial number for each city
for ii = 1:Ncity
    Ssave(ii,1)=SIR(1,ii);
    Isave(ii,1)=SIR(2,ii);
    Rsave(ii,1)=SIR(3,ii);
    Nsave(ii,1)=N(1,ii);
end
tsave(1)=0;
%Note that (migMat .* SIR) is the migration matrix in the equation sheet
%flowMat is the flow matrix in the equation sheet
%TrueMigMat is the Net migration matrix in the equation sheet
for clock=1:clockmax
    t = clock*dt;
    TrueMigMat = migMat .* SIR * flowMat;
    for cur = 1:Ncity
    %S = SIR(1,Ncity)
    %I = SIR(2,Ncity)
    %R = SIR(3,Ncity)
    %N = N(1,Ncity) 
    S2I(cur) = dt*SIR(1,cur)*(SIR(2,cur)/N(1,cur))*a(cur);
    I2R(cur) = dt*SIR(2,cur)*b(cur);
    I2D(cur) = dt*SIR(2,cur)*d(cur); 
    
    %Net Migration
    S2S(cur) = dt*TrueMigMat(1,cur); 
    I2I(cur) = dt*TrueMigMat(2,cur);
    R2R(cur) = dt*TrueMigMat(3,cur); 
    
    %update population
    SIR(1,cur) = SIR(1,cur)- S2I(cur) + S2S(cur);
    SIR(2,cur) = SIR(2,cur) + S2I(cur) - I2R(cur) + I2I(cur) - I2D(cur);
    SIR(3,cur) = SIR(3,cur) + I2R(cur) + R2R(cur);
    N(1,cur) = SIR(1,cur) + SIR(2,cur) + SIR(3,cur);
    D(cur) = D(cur) + I2D(cur);
    %save points for graphing
    Ssave(cur, clock+1) = SIR(1,cur);
    Isave(cur, clock+1) = SIR(2,cur);
    Rsave(cur, clock+1) = SIR(3,cur);
    
    Dsave(cur, clock+1) = D(cur);
    
    Nsave(cur, clock+1) = N(1,cur);
    end
    tsave(1, clock+1) = t;
end



%plotting simple(identical neighbors), central city, and sum
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%if complexplot == false
fwid = 800;
fhigh = 800;
figure('Position',[0,0,fwid,fhigh]);
 
subplot(2,2,1);

plot(tsave(1),Ssave(1,1),'ro',tsave(1),Isave(1,1),'bo',...
    tsave(1),Rsave(1,1),'go',tsave(1),Nsave(1,1),'ko'...
    ,tsave(1),Dsave(1,1),'mo','markers',3);
title('central island/city')
axis([0,clockmax*dt,0,1.5*(Nsave(1,1))])
axis manual


subplot(2,2,2);
plot(tsave(1),Ssave(2,1),'ro',tsave(1),Isave(2,1),'bo',...
    tsave(1),Rsave(2,1),'go',tsave(1),Nsave(2,1),'ko'...
    ,tsave(1),Dsave(2,1),'mo','markers',3);
title('one of the neighbor island/city')
axis([0,clockmax*dt,0,1.5*(Nsave(2,1))])
axis manual



subplot(2,2,3);
plot(tsave(1),sum(Ssave(:,1)),'ro',tsave(1),sum(Isave(:,1)),'bo',...
    tsave(1),sum(Rsave(:,1)),'go',tsave(1),sum(Nsave(:,1)),'ko'...
    ,tsave(1),sum(Dsave(:,1)),'mo','markers',3)
title('summation of all island/city')
axis manual
axis([0,clockmax*dt,0,1.2*sum(Nsave(:,1))])

legend({'Susptible','Infected','Recovered','Sum','Dead'}, [fwid/2 ,fhigh/8,0.5,0.5])


for i=2:clockmax+1
    for iter = 1:2
        subplot(2,2,iter);     
        hold on
        plot(tsave(i),Ssave(iter,i),'ro',tsave(i),Isave(iter,i),'bo',...
             tsave(i),Rsave(iter,i),'go',tsave(i),Nsave(iter,i),'ko',...
             tsave(i),Dsave(iter,i),'mo','markers',3);
        hold off
    end

    
    subplot(2,2,3);
    hold on 
    plot(tsave(i),sum(Ssave(:,i)),'ro',tsave(i),sum(Isave(:,i)),'bo',...
    tsave(i),sum(Rsave(:,i)),'go',tsave(i),sum(Nsave(:,i)),'ko',...
    tsave(i),sum(Dsave(:,i)),'mo','markers',3);
    
    
    end
hold off
