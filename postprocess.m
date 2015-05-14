%%% plot the tokance
%%% tqdata contains tl and th, the tokance array at two differential bias arrays
%%% Need to input Vdv, dVd, and sita
%%% plot 3 components of the tokance.
load tqdata
q=1.6e-19;
Vdv=0.1*[0:length(tl)-1];
dVd=0.01;
sita=0.05;
dtq=th-tl;
tokance=q*dtq./dVd./sin(sita);  % in A/m^2
figure()
plot(Vdv, tokance(:,1),'r--','linewidth',[2]); hold on;
plot(Vdv, tokance(:,2),'b--','linewidth',[2]); hold on;
plot(Vdv, tokance(:,3),'g--','linewidth',[2]); hold on;


