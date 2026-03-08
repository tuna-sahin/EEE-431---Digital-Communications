%% pick p(t) as a triangular function
clf;
T=2;
x = linspace(-T/2,T/2,21);
x2 = linspace(0,T,21);
p = 1-abs(x);
subplot(1,3,1)
fineplot(x2,p,"s_0(t)","t","p(t)",[-0.2 2.5],[-1.1 1.1],"off",[450 600],"\Psi(t)","-");
grid on
subplot(1,3,2)
fineplot(x2,-p,"s_1(t)","t","-p(t)",[-0.2 2.5],[-1.1 1.1],"off",[450 600],"\Psi(t)","-");
grid on
subplot(1,3,3)
psi = p/norm(p);
fineplot(x2,psi,"Basis Function","t","\Psi(t)",[-0.2 2.5],[-0.75 0.75],"off",[900 300],"\Psi(t)","-");
grid on


%% Signal constellation:
clf;
fineplot([-1*norm(p) norm(p)],[0 0],"Signal Constellation","\psi_1","",[-1.2*norm(p) 1.2*norm(p)],[-1 1],"off",[600 160],"","x");
xticks([-norm(p) -2 -1 0 1 2 norm(p)]);
text([-1*norm(p)-0.3 norm(p)-0.3],[0.3 0.3],{"m=1","m=0"});

%% b)
clf;
fineplot([-1*norm(p) norm(p)],[0 0],"Signal Constellation","\psi_1","",[-1.2*norm(p) 1.2*norm(p)],[-1 1],"on",[600 160],"","x");
xticks([-norm(p) -2 -1 0 1 2 norm(p)]);
yticks([])
text([-1*norm(p)-0.3 norm(p)-0.3],[0.3 0.3],{"m=1","m=0"});
d0 = plot([0 norm(p)*1.2],[0 0],"LineWidth",3,"Color",[0 0.9 0.9 0.3]);
d1 = plot([-norm(p)*1.2 0],[0 0],"LineWidth",3,"Color",[0.9 0.9 0 0.3]);
legend([d0 d1],"Decision Reigon for 0","Decision Reigon for 1","Location","south");

%% f)
clf;
gamma_bar = 1;%for plotting purposes
threshold = (norm(p)/(4*gamma_bar))*log(1/3);
fineplot([-1*norm(p) norm(p)],[0 0],"Signal Constellation","\psi_1","",[-1.2*norm(p) 1.2*norm(p)],[-1 1],"on",[600 160],"","x");
xticks([-norm(p) -2 -1 0 1 2 norm(p)]);
yticks([]);

plot([threshold threshold],[-0.1 0.1])
text(threshold-0.5,0.2,strcat("Threshold = ",num2str(threshold,"%.2f")),"FontSize",8)
text([-1*norm(p)-0.3 norm(p)-0.3],[0.3 0.3],{"m=1","m=0"});
d0 = plot([threshold norm(p)*1.2],[0 0],"LineWidth",3,"Color",[0 0.9 0.9 0.3]);
d1 = plot([-norm(p)*1.2 threshold],[0 0],"LineWidth",3,"Color",[0.9 0.9 0 0.3]);
legend([d0 d1],"Decision Reigon for 0","Decision Reigon for 1","Location","south");

%% montecarlos)
gamma_bars = linspace((qfuncinv(1/2).^2)/2,(qfuncinv(10^(-4)).^2)/2,20);
N_0s = norm(p)^2./gamma_bars;
rng(67);
p_errors1 = zeros(like=gamma_bars);
for i = 1:length(gamma_bars)
    s = (((rand([1 10^6])-0.5) > 0) .* 2 - 1) .*  norm(p); %transmit 0 or 1 with equal probabilities
    N_0 = (norm(p)^2)/gamma_bars(i);
    n = normrnd(0,sqrt(N_0/2),[1 10^6]);
    r = s + n;
    s_hat = ((r > 0) .* 2 - 1) .* norm(p);
    errors = s ~= s_hat;
    p_errors1(i) = mean(errors);
end

rng(67);
p_errors2 = zeros(like=gamma_bars);
for i = 1:length(gamma_bars)
    s = (((rand([1 10^6])-0.25) > 0) .* 2 - 1) .*  norm(p); %transmit 0 or 1 with equal probabilities
    N_0 = (norm(p)^2)/gamma_bars(i);
    n = normrnd(0,sqrt(N_0/2),[1 10^6]);
    r = s + n;
    threshold = (N_0 / (4*(norm(p)^2))) * log(1/3);
    s_hat = ((r > threshold) .* 2 - 1) .* norm(p);
    errors = s ~= s_hat;
    p_errors2(i) = mean(errors);
end

rng(67);
p_errors3 = zeros(like=gamma_bars);
for i = 1:length(gamma_bars)
    s = (((rand([1 10^6])-0.25) > 0) .* 2 - 1) .*  norm(p); %transmit 0 or 1 with equal probabilities
    N_0 = (norm(p)^2)/gamma_bars(i);
    n = normrnd(0,sqrt(N_0/2),[1 10^6]);
    r = s + n;
    s_hat = ((r > 0) .* 2 - 1) .* norm(p);
    p_errors3(i) = mean(s ~= s_hat);
end

%% Plots
clf;
subplot(1,3,1)
semilogy(gamma_bars,p_errors1)
hold on
grid on
semilogy(gamma_bars,qfunc(sqrt(2*gamma_bars)),"x")
title("ML Rule with Equal Priors")
legend({"Simulation Results","Theoretical Results"},"location","south")
subplot(1,3,2)
semilogy(gamma_bars,p_errors2,"--")
hold on 
grid on
semilogy(gamma_bars,0.25*(qfunc(sqrt(2.*gamma_bars)+threshold./sqrt(0.5.*N_0s)))+0.75*qfunc(sqrt(2.*gamma_bars)-threshold./sqrt(0.5.*N_0s)),"x");
title("MAP Rule")
legend({"Simulation Results","Theoretical Results"},"location","south")
subplot(1,3,3)
semilogy(gamma_bars,p_errors3,"-.")
hold on
grid on
semilogy(gamma_bars,qfunc(sqrt(2*gamma_bars)),"x")
title("ML Rule with Nonequal Priors")
set(gcf,'position',[6 10^-3 2100 500]);
legend({"Simulation Results","Theoretical Results"},"location","south")

%%
figure;
hold off
semilogy(gamma_bars,p_errors1)
hold on
grid on
semilogy(gamma_bars,p_errors2,"-")
semilogy(gamma_bars,p_errors3,"--")
title("Comparison of Different Receivers")
grid on
legend({"ML Rule (Equal Priors)","MAP Rule","ML Rule (Nonequal Priors)"},"Location","south")
set(gcf,'position',[6 10^-3 600 500]);

%% Question 2
clf;
A = 5;
fineplot([0 0 0],[0 0 0],"Initial Signal Constellation","\psi_1","\psi_2",[-1.2*A 1.2*A],[-1.2*A 1.2*A],"on",[400 400],"","k");
xticks([-A -A/2 0 A/2 A])
xticklabels({"-A","-A/2","","A/2","A"})
yticks([-A -A/2 0 A/2 A])
yticklabels({"-A","-A/2","","A/2","A"})
grid on
plot([-A -A A A],[-A 0 -A A],"x","LineWidth",2)
plot([-1.2*A 0],[-A/2 -A/2],"b--","LineWidth",1)
plot([A/4 1.2*A],[0 0],"b--","LineWidth",1)
plot([0 A/4],[-A/2 0],"b--","LineWidth",1)
plot([0 0],[-A/2 -1.2*A],"b--","LineWidth",1.2)
plot([-0.35*A A/4],[1.2*A 0],"b--","LineWidth",1)
text(-A-A/20,A/10,"s_4")
text(-A-A/20,-9*A/10,"s_3")
text(A-A/20,-9*A/10,"s_2")
text(A-A/20,11*A/10,"s_1")

%% b)
Es = (norm([A A]).^2+norm([A -A]).^2+norm([-A -A]).^2+norm([-A 0]).^2)/4;
gammabars = linspace(1,40,20);
N0s = Es./gammabars;
pe = zeros(like=N0s);
samples = 10^7;
rng(69);
for i = 1:length(gammabars)
    sm1 = (((((rand([1 samples])-0.5)>0)*2)-1) .* A);
    sm2 = ((((rand([1 samples])-0.5)>0)-1+(0.5*(sm1>0)))).*((sm1>0)+1).*A;
    
    n1 = normrnd(0,sqrt(N0s(i)/2),[1 samples]);
    n2 = normrnd(0,sqrt(N0s(i)/2),[1 samples]);
    
    r1 = sm1 + n1;
    r2 = sm2 + n2;

    sm1_hat = -A * ((r1 <= 0) .* (r2 <= -A/2)) + A * ((r1>=0) .* (r2<=0) .* (r2-2*r1<=-1)) + A * ((r2>=0) .* (r2+2*r1>=1)) - A * ((r2>=-A/2) .* (r2-2*r1>=-1) .* (r2+2*r1<=1));
    sm2_hat = -A * ((r1 < 0) .* (r2 < -A/2)) - A * ((r1>0) .* (r2<0) .* (r2-2*r1<-1)) + A * ((r2>0) .* (r2+2*r1>1));
    
    pe(i) = mean(((sm1 ~= sm1_hat) + (sm2 ~= sm2_hat))>0);
end

%%
clf;
semilogy(gammabars,pe)    
hold on
grid on
semilogy(gammabars,qfunc(sqrt(gammabars*8/7)) + qfunc(sqrt(gammabars*10/7)) + 0.5 * qfunc(sqrt(gammabars*16/7)) + 0.5 * qfunc(sqrt(gammabars*2/7)),"x")
title("Error Probability of Initial Constellation")
legend({"Simulation Results","Theoretical Union Bound"})

%% part d) constellation

clf;
A = 5;
fineplot([0 0 0],[0 0 0],"Modified Signal Constellation","\psi_1","\psi_2",[-1.2*A 1.2*A],[-1.2*A 1.2*A],"on",[400 400],"","k");
xticks([-A -A/2 0 A/2 A])
xticklabels({"-A","-A/2","","A/2","A"})
yticks([-A -A/2 0 A/2 A])
yticklabels({"-A","-A/2","","A/2","A"})
grid on
plot([-A -A A A],[-A A -A A],"x","LineWidth",2)
plot([-1.2*A 1.2*A],[0 0],"b--","LineWidth",1)
plot([0 0],[-1.2*A 1.2*A],"b--","LineWidth",1)
text(-A-A/20,11*A/10,"s_4")
text(-A-A/20,-9*A/10,"s_3")
text(A-A/20,-9*A/10,"s_2")
text(A-A/20,11*A/10,"s_1")

%% part d) simulation
Es = (norm([A A]).^2+norm([A -A]).^2+norm([-A -A]).^2+norm([-A A]).^2)/4;
N0s = Es./gammabars;
pe2 = zeros(like=N0s);
rng(70);
for i = 1:length(gammabars)
    sm1 = (((((rand([1 samples])-0.5)>0)*2)-1) .* A);
    sm2 = (((((rand([1 samples])-0.5)>0)*2)-1) .* A);
    
    n1 = normrnd(0,sqrt(N0s(i)/2),[1 samples]);
    n2 = normrnd(0,sqrt(N0s(i)/2),[1 samples]);
    
    r1 = sm1 + n1;
    r2 = sm2 + n2;

    sm1_hat = -A* (r1<0) + A*(r1>0);
    sm2_hat = -A* (r2<0) + A*(r2>0);
    
    pe2(i) = mean(((sm1 ~= sm1_hat) + (sm2 ~= sm2_hat))>0);
end

%%
clf;
semilogy(gammabars,pe2)
hold on
grid on
semilogy(gammabars,pe)    
title("Comparison of Error Probabilities")
% semilogy(gammabars,qfunc(sqrt(gammabars*8/7)) + qfunc(sqrt(gammabars*10/7)) + 0.5 * qfunc(sqrt(gammabars*16/7)) + 0.5 * qfunc(sqrt(gammabars*2/7)),"x")
legend({"Updated Constellation","Initial Constellation"})

%% e)
Es = (norm([A A]).^2+norm([A -A]).^2+norm([-A -A]).^2+norm([-A A]).^2)/4;
N0s = Es./gammabars;
pe3 = zeros(like=N0s);
rng(70);
for i = 1:length(gammabars)
    sm1 = (((((rand([1 samples])-0.5)>0)*2)-1) .* A);
    sm2 = (((((rand([1 samples])-0.5)>0)*2)-1) .* A);
    
    b1 = (sm1>0);
    b2 = sm1.*sm2 < 0;

    n1 = normrnd(0,sqrt(N0s(i)/2),[1 samples]);
    n2 = normrnd(0,sqrt(N0s(i)/2),[1 samples]);
    
    r1 = sm1 + n1;
    r2 = sm2 + n2;

    b1_hat = (r1>0);
    b2_hat = r1.*r2 < 0;

    pe3(i) = mean(((b1 ~= b1_hat) + (b2 ~= b2_hat)));
end

%%
clf;

semilogy(gammabars./2,pe2)
hold on
grid on
semilogy(gammabars./2,pe3)
title("Symbol and Bit Error Probabilities")
% semilogy(gammabars,qfunc(sqrt(gammabars*8/7)) + qfunc(sqrt(gammabars*10/7)) + 0.5 * qfunc(sqrt(gammabars*16/7)) + 0.5 * qfunc(sqrt(gammabars*2/7)),"x")
legend({"Symbol Error Probability","Bit Error Probability"})

%% f)

Es = (norm([A A]).^2+norm([A -A]).^2+norm([-A -A]).^2+norm([-A A]).^2)/4;
N0s = Es./gammabars;
pe4 = zeros(like=N0s);
rng(70);
for i = 1:length(gammabars)
    sm1 = (((((rand([1 samples])-0.5)>0)*2)-1) .* A);
    sm2 = (((((rand([1 samples])-0.5)>0)*2)-1) .* A);
    
    b1 = sm1 > 0;
    b2 = sm2 > 0;

    n1 = normrnd(0,sqrt(N0s(i)/2),[1 samples]);
    n2 = normrnd(0,sqrt(N0s(i)/2),[1 samples]);
    
    r1 = sm1 + n1;
    r2 = sm2 + n2;

    b1_hat = r1 > 0;
    b2_hat = r2 > 0;

    
    pe4(i) = mean(((b1 ~= b1_hat) + (b2 ~= b2_hat)));
end
%%

%%
clf;
subplot(1,2,1)
semilogy(gammabars./2,pe4)
hold on
grid on
semilogy(gammabars./2,pe3,"-.")
semilogy(gammabars./2,pe2,"--")
title("Comparison of Bit and Symbol Errors")
legend({"With Gray Mapping","With Ordinal Coding","Symbol Error Probability"})
subplot(1,2,2)
semilogy(gammabars./2,pe4)
hold on
grid on
semilogy(gammabars./2,pe3,"-.")
semilogy(gammabars./2,pe2,"--")
title("Comparison of Bit and Symbol Errors (Zoomed)")
legend({"With Gray Coding","With Ordinal Coding","Symbol Error Probability"})
xlim([0 2])
ylim([10^-1 1])
set(gcf,'position',[6 10^-3 1400 700]);
