% Simulations of:
% Two coupled FitzHugh Nagumo with Dead Zone coupling
%
% Peter Ashwin, June 2021
%
% Plotting modifications, CB 2021
%
% for paper "Dead zones and phase reduction of coupled oscillators"
% by Ashwin, Bick and Poignard 2021

clear variables

% Line width
LiWi = 1.3;

twopi=2*pi;

% parameters: bias
p.i=0.33;
% timescale sep
p.mu=0.05;

% number of steps in discretization
nsteps=200;

% rotation centre to define angle
cy1=0;
cy2=0.5;
pangle=@(y1,y2) (unwrap(atan2(y1-cy1,cy2-y2)));

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
clf;
f1.PaperUnits='centimeters';
f1.PaperSize=[14 4];
f1.Units='centimeters';
f1.InnerPosition=[0 0 f1.PaperSize];
offst = 0.2;


subplot(1,7,[1 2 3 4])
box on
hold on
plot(fhntheta,fhny(2,:),'LineWidth',LiWi,'Color',[0.7, 0.7, 0.7]);
plot(fhntheta,fhny(1,:),'LineWidth',LiWi,'Color',[0 0 0]);
text(-1.05, 2, '{\textbf{(a)}}','Interpreter','latex');
xlim([0 2*pi]);
xlabel('$\theta$', 'Interpreter','latex');
%tStr = sprintf('\\color[rgb]{%f, %f, %f}%s', titleColor, titleString);
ylabel('$v, w$', 'Interpreter','latex');
xticks([0 pi 2*pi]);
set(gca,'TickLabelInterpreter','latex')
xticklabels({'$0$','$\pi$','$2\pi$'});
pos = get(gca, 'Position');
set(gca, 'Position', [pos(1) pos(2)+offst pos(3) pos(4)-offst]);

subplot(1,7, [6 7])
plot(fhntheta,fhnpsi,'LineWidth',LiWi,'Color',[0 0 0]);
text(-2.6, max(fhnpsi), '{\textbf{(b)}}','Interpreter','latex');
xlim([0 2*pi]);
ylim([0 2*pi]);
xlabel('$\theta$', 'Interpreter','latex');
ylabel('$\psi$', 'Interpreter','latex');
set(gca,'TickLabelInterpreter','latex')
xticks([0 pi 2*pi]);
xticklabels({'$0$','$\pi$','$2\pi$'});
yticks([0 pi 2*pi]);
yticklabels({'$0$','$\pi$','$2\pi$'});
pos = get(gca, 'Position');
set(gca, 'Position', [pos(1) pos(2)+offst pos(3) pos(4)-offst]);

print(f1,'fhn_dz_v6_f1.pdf','-dpdf', '-r600');



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

% plot soln of AVE etc for single oscillator
f2=figure(2);
clf;
%f2.PaperUnits='centimeters';
%f2.PaperSize=[10 8];
%f2.Units='centimeters';
%f2.InnerPosition=[0 0 f2.PaperSize];
f2.PaperUnits='centimeters';
f2.PaperSize=[14 7];
f2.Units='centimeters';
f2.InnerPosition=[3 3 f2.PaperSize];

% arrays giving soln of AVE and RHS of FHN vs T
tt=[0:1/(nsteps-1):1]*period;
avtt=avenorm(tt);
fhnlctt=fhnlc(tt*twopi/period);
fhntt=fhn(fhnlctt,p);

% should be constant
checkQQ = avtt.*fhntt;

% sketch 
subplot(2,2,1);
box on
hold on;
plot(tt,avtt(1,:),'LineWidth',LiWi,'Color',[0, 0, 0]);
text(-1.2, max(avtt(1,:)), '{\textbf{(a)}}','Interpreter','latex');
xlim([0 period]);
xlabel('$t$', 'Interpreter','latex');
ylabel('$Z=Z_v$', 'Interpreter','latex');
set(gca,'TickLabelInterpreter','latex')

subplot(2,2,2);
plot(tt,avtt(2,:),'LineWidth',LiWi,'Color',[0 0 0]);
xlabel('$t$', 'Interpreter','latex');
ylabel('$Z_w$', 'Interpreter','latex');
xlim([0 period]);
text(-0.95, max(avtt(2,:)), '{\textbf{(b)}}','Interpreter','latex');
xlim([0 period]);
set(gca,'TickLabelInterpreter','latex')

H=zeros(size(tt));
for i=1:nsteps
    temp=0;
    % phase difference in [0,1)
    for j=1:nsteps-1
        Q=avtt(:,j);
        gamma=fhnlctt(:,j);
        gammashift=fhnlctt(:,mod(j+i-2,nsteps)+1);
        temp=temp+Q(1)*g(gamma(1),gamma(2),gammashift(1),gammashift(2),p)/nsteps;
    end
    H(i)=temp;
end
% Odd part of H on times tt
Hodd=(H-fliplr(H));

subplot(2,2,3);
box on
hold on
plot(tt,H,'LineWidth',LiWi,'Color',[0, 0, 0]);
xlabel('$t$', 'Interpreter','latex');
ylabel('$\mathrm{h}$', 'Interpreter','latex');
xlim([0 period]);
ph = linspace(0,2*pi, nsteps);
plot(ph, zeros(size(ph)), ':k')
text(-1.2, max(H), '{\textbf{(c)}}','Interpreter','latex');
set(gca,'TickLabelInterpreter','latex')

subplot(2,2,4);
box on
hold on
plot(tt,Hodd,'LineWidth',LiWi,'Color',[0, 0, 0]);
xlabel('$t$', 'Interpreter','latex');
ylabel('$\mathrm{h}_\mathrm{odd}$', 'Interpreter','latex');
xlim([0 period]);
ph = linspace(0,2*pi, nsteps);
plot(ph, zeros(size(ph)), ':k')
text(-0.95, max(Hodd), '{\textbf{(d)}}','Interpreter','latex');
set(gca,'TickLabelInterpreter','latex')

print(f2,'fhn_dz_v6_f2.pdf','-dpdf', '-r600');

return


%% simulate ensemble of trajectories for two coupled FHNs
ntraj=100;
tspan=60;

% set up figures
f3=figure(3);
clf;
f3.PaperUnits='centimeters';
f3.PaperSize=[14 8];
f3.Units='centimeters';
f3.InnerPosition=[1 0 f3.PaperSize];

f4=figure(4);
clf;
f4.PaperUnits='centimeters';
f4.PaperSize=[14 8];
f4.Units='centimeters';
f4.InnerPosition=[1 0 f4.PaperSize];


for j=3:4
    figure(j);
    if j==3
        p.eps=0.05;
    else
        p.eps=-0.05;
    end
    
    cmap = 0.8*[gray(ntraj/2); 1-gray(ntraj/2)];

for i=1:ntraj
    % start second oscillator at fixed phase on limit cycle
    yi1=fhnlc(0);
    yi2=fhnlc(twopi*i/(ntraj+1));
    y0=[yi1(1);yi1(2);yi2(1);yi2(2)];
    
    options = odeset('RelTol',1e-6,'AbsTol',1e-7,'MaxStep',0.1);
    [tt,yy]=ode45(@(t,y) frhs(y,p),[0 tspan],y0, options);
    psi1=pangle(yy(:,1),yy(:,2));
    psi2=pangle(yy(:,3),yy(:,4));
    ph1=psitotheta(psi1);
    ph2=psitotheta(psi2);
    
    % plot timeseries of v_1-v_2
    subplot(2,1,1);
    plot(tt,yy(:,1)-yy(:,3),'Color', cmap(i, :), 'LineWidth', 1);
    hold on;
    
    % plot timeseries of theta_1-theta_2
    subplot(2,1,2);
    plot(tt,mod(ph1-ph2,twopi),'Color', cmap(i, :),'LineWidth', 1);
    hold on;
end
subplot(2,1,1)
hold off;
xlabel('$t$', 'Interpreter','latex');
ylabel('$v_1-v_2$', 'Interpreter','latex');
text(-7, 5, sprintf('{\\textbf{(a%i)}}',j-2),'Interpreter','latex');
set(gca,'TickLabelInterpreter','latex')

subplot(2,1,2)
hold off;
xlabel('$t$', 'Interpreter','latex');
ylabel('$\theta_1-\theta_2$', 'Interpreter','latex');
text(-7, twopi, sprintf('{\\textbf{(b%i)}}',j-2),'Interpreter','latex');
set(gca,'TickLabelInterpreter','latex')
ylim([0,twopi])
yticks([0 pi 2*pi]);
yticklabels({'$0$','$\pi$','$2\pi$'});



end


print(f2,'fhn_dz_v6_f2.pdf','-dpdf', '-r600');
print(f3,'fhn_dz_v6_f3.pdf','-dpdf', '-r600');
print(f4,'fhn_dz_v6_f4.pdf','-dpdf', '-r600');

%keyboard;
disp('DONE!')
return
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

function temp=g(v1,w1,v2,w2,p)
% DZ in both v1 and v2
temp=zeros(size(v1));
if v1>0 && v2>0
    temp=v1*v2;
end
% DZ in v2 only
% temp=zeros(size(v2));
% if v2>0
%     temp=v2;
% end
% coupling depends on fast derivative
%temp=(f(v2,w2)+p.i);
%
% temp=(w2-1.1).*(w2>1.1);
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
