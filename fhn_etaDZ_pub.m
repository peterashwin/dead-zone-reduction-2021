% Simulations of:
% Two coupled FitzHugh Nagumo with Dead Zone coupling
%
% Peter Ashwin
%
% June 2021

% Modified file for approximate dead zones, CB Jul 2021


clc
clear variables

twopi=2*pi;
LiWi = 1.3;
Lbls = {};

% parameters: bias
p.i=0.33;
% timescale sep
mus = [0.05, 0.01, 0.005,0.001];
% pulse width
pw = 1;

% number of steps in discretization
nsteps=500;

% rotation centre to define angle
cy1=0;
cy2=0.5;
pangle=@(y1,y2) (unwrap(atan2(y1-cy1,cy2-y2)));


% plot solutions and angles for single oscillator
f1=figure(1);
clf;
f1.PaperUnits='centimeters';
f1.PaperSize=[15 8];
f1.Units='centimeters';
f1.InnerPosition=[3 3 15 8];

cmap = linspace(0.4,1,length(mus))'*ones(1,3);

for ii=1:length(mus)
    p.mu = mus(ii);
    Lbls{ii} = sprintf('$\\mu=%g$', p.mu);


% compute trajectory for uncoupled FHN to get phase
y0=[2; 0];
options = odeset('RelTol',1e-6,'AbsTol',1e-7,'MaxStep',0.1,'Events',@myEventF);
sol=ode45(@(t,y) fhn(y,p),[0 10],y0, options);

% find period
ta=sol.xe(end-1);
tb=sol.xe(end);
period=tb-ta;

% parameterize limit cycle by phase theta in (0,1)
fhnlc=@(theta) (deval(sol,ta+mod(theta,twopi)*period/twopi));

% create lookup of phase theta vs angle psi
% theta proceeds with theta' constant
% psi is geometric angle
fhntheta=[0:1/(nsteps-1):1]*twopi;
fhny=fhnlc(fhntheta);
fhnpsi=pangle(fhny(1,:),fhny(2,:));

% angle to phase
psitotheta=@(psi) (interp1(fhnpsi,fhntheta,mod(psi,twopi))+twopi*floor(psi/twopi));

% check: this should be identity
fhnthetacheck=psitotheta(fhnpsi);

% plot solutions and angles for single oscillator
f1=figure(1);
subplot(2,2,1)
box on
hold on
plot(fhntheta,fhny(1,:),'Color', 1-cmap(ii, :),'LineWidth', LiWi);
text(-1.6, 2, '{\textbf{(a)}}','Interpreter','latex');
xlim([0 2*pi]);
xlabel('$\theta$', 'Interpreter','latex');
ylabel('$v$', 'Interpreter','latex');
xticks([0 pi 2*pi]);
set(gca,'TickLabelInterpreter','latex')
xticklabels({'$0$','$\pi$','$2\pi$'});
legend(Lbls,'Interpreter','latex');



%% compute approx of branches of gamma_0(t) and hence H(chi) from Izhikevich (3.2)

% set up AVE
fhnave=@(t,dv,p) (-fhnjac(fhnlc(t*twopi/period),p)'*dv);

% solve to find periodic solution for AVE by running backwards in time
y0=[0.5 ; 0.5];
options = odeset('RelTol',1e-6,'AbsTol',1e-7,'MaxStep',0.1,'Events',@myEventF);
solave=ode45(@(t,y) fhnave(t,y,p),[period*10 0],y0, options);

% normalize using values at t=0
aveQQ= deval(solave,0)'*fhn(fhnlc(0),p);

% normalised soln of ave
avenorm=@(t) deval(solave,mod(t,period))/aveQQ;

% arrays giving soln of AVE and RHS of FHN vs T
tt=[0:1/(nsteps-1):1]*period;
avtt=avenorm(tt);
fhnlctt=fhnlc(tt*twopi/period);
fhntt=fhn(fhnlctt,p);

% should be constant
checkQQ = avtt.*fhntt;

% sketch 
f1 = figure(1);
subplot(2,2,3);
box on
hold on;
plot((tt/max(tt)*2*pi),avtt(1,:)/p.mu,'Color', 1-cmap(ii, :),'LineWidth', LiWi);
ph = linspace(0,2*pi, nsteps);
plot(ph, zeros(size(ph)), ':k')
xlabel('$\theta$','Interpreter','latex');
ylabel('$Z/\mu$','Interpreter','latex');
xlim([0 2*pi]);
if ii > 3
text(-1.6, max(avtt(1,:)/p.mu), '{\textbf{(c)}}','Interpreter','latex');
end
xticks([0 pi 2*pi]);
set(gca,'TickLabelInterpreter','latex')
xticklabels({'$0$','$\pi$','$2\pi$'});

%subplot(2,2,2);
%plot(tt,avtt(2,:));
%xlabel('$t$','Interpreter','latex');
%ylabel('$q$','Interpreter','latex');
%xlim([0 period]);
%text(-1, 3, '{\textbf{(b)}}','Interpreter','latex');

% Pulse
P = @Pc; % Cosine pulse
P = @Pd; % Exp pulse
nrm = trapz(ph, P(ph, 1, 1));
subplot(2,2,2);
%hold on
box on
plot(ph, P(ph,nrm,pw), 'Color', 1-cmap(ii, :),'LineWidth', LiWi)
xlabel('$\theta$','Interpreter','latex');
ylabel('$P$','Interpreter','latex');
xlim([0 2*pi]);
text(-1.4, 2, '{\textbf{(b)}}','Interpreter','latex');
xticks([0 pi 2*pi]);
set(gca,'TickLabelInterpreter','latex')
xticklabels({'$0$','$\pi$','$2\pi$'});
%return

H=zeros(size(tt));
for i=1:nsteps
    temp=0;
    % phase difference in [0,1)
    for j=1:nsteps
        Q=avtt(:,j);
        %gamma=fhnlctt(:,j);
        %gammashift=fhnlctt(:,mod(j+i-2,nsteps)+1);
        shft = mod(i*2*pi/nsteps+j*2*pi/nsteps, 2*pi);
        %temp=temp+Q(1)*g(gamma(1),gamma(2),gammashift(1),gammashift(2),p)/nsteps;
        %gammashift(1)
        temp=temp+Q(1)*P(shft,nrm,pw)/nsteps;
    end
    H(i)=temp;
end
% Odd part of H on times tt
Hodd=(H-fliplr(H));

subplot(2,2,4);
box on
hold on
plot((tt/max(tt)*2*pi),H/p.mu,'Color', 1-cmap(ii, :),'LineWidth', LiWi); %/max(abs(H)
ph = linspace(0,2*pi, nsteps);
plot(ph, zeros(size(ph)), ':k')
if ii>3
text(-1.4, max(H/p.mu), '{\textbf{(d)}}','Interpreter','latex');
end
xlabel('$\theta$','Interpreter','latex');
ylabel('$\mathrm{h}/\mu$','Interpreter','latex');
xlim([0 2*pi]);
%text(-1, 3, '{\textbf{(a)}}','Interpreter','latex');
xticks([0 pi 2*pi]);
set(gca,'TickLabelInterpreter','latex')
xticklabels({'$0$','$\pi$','$2\pi$'});


%subplot(2,2,4);
%plot(tt,Hodd);
%xlabel('$t$','Interpreter','latex');
%ylabel('$H_\mathrm{odd}$','Interpreter','latex');
%xlim([0 period]);


end

%return




print(f1,'fhn_dz_etaDZ_v2_f1.pdf','-dpdf');
%print(f2,'fhn_dz_v4_f2.pdf','-dpdf');
%print(f3,'fhn_dz_v4_f3.pdf','-dpdf');
%print(f4,'fhn_dz_v4_f4.pdf','-dpdf');

%keyboard;
disp('DONE!')
%
% end



%% rhs of coupled FNH oscillator system
function dydt=frhs(y,p)
dydt=zeros(size(y));

v1=y(1);
w1=y(2);
v2=y(3);
w2=y(4);

c1=g(v1,w1,v2,w2,p);
c2=g(v2,w2,v1,w1,p);

dydt(1)=(f(v1,w1)+p.i+p.eps*c1)/p.mu;
dydt(2)=v1+0.7-0.8*w1;
dydt(3)=(f(v2,w2)+p.i+p.eps*c2)/p.mu;
dydt(4)=v2+0.7-0.8*w2;
end

%nonlinearity
function temp=f(v,w)
temp=v-v.^3/3-w;
end

% df/dv
function temp=df1(v,w)
temp=1-v.^2;
end

% df/dw
function temp=df2(v,w)
temp=-1;
end


function temp=Pc(theta,nrm,pw)
temp=(((cos(theta)+1)/2).^pw)/nrm;
%pw
temp=((cos(theta)+1)/2).^3;
end


function y = uFuncBump(x, D, dlta)
% Localized bump function
y = (abs(x)<dlta).*(D*exp((-dlta^2)./(dlta^2-(x.*x))));
y(isnan(y)) = 0;
end

function temp=Pd(theta,nrm,~)
        % Bump

        % Params
        D = 1;
        %D = 0;
        dlta = 0.5;
        shft = 1;
        
        bmpx = @(x) uFuncBump(x, D, dlta);
        bmpd = @(x) bmpx(mod(x-shft, 2*pi)) + bmpx(mod(x-shft, 2*pi)-2*pi);
        temp = bmpd(theta)/nrm;
end




function temp=g(v1,w1,v2,w2,p)
% % DZ in both v1 and v2
% temp=zeros(size(v1));
% if v1>0 && v2>0
%     temp=v1*v2;
% end
% DZ in v2 only
% temp=zeros(size(v2));
% if v2>0
%     temp=v2;
% end
% coupling depends on fast derivative
%temp=(f(v2,w2)+p.i);
%
temp=(w2-1.1).*(w2>1.1);
end

%% rhs of single FNH oscillator
function dydt=fhn(y,p)
dydt=zeros(size(y));

v1=y(1,:);
w1=y(2,:);

dydt(1,:)=(f(v1,w1)+p.i)/p.mu;
dydt(2,:)=v1+0.7-0.8*w1;
end

% Used to find period of single FHN oscillator
function [v,ist,dir]=myEventF(t,y)
v=y(1); % event when v=0
ist=0; % don't stop
dir=1; % positive going
end

%% jacobian for single FNH oscillator
function jac=fhnjac(y,p)
v1=y(1);
w1=y(2);
jac=zeros(2,2);

jac(1,1)= df1(v1,w1)/p.mu;
jac(1,2)= df2(v1,w1)/p.mu;
jac(2,1)=1;
jac(2,2)=-0.8;
end
