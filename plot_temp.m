% Plot different properties vs T
t = 1; u = 15.3;
m = 3;
n = -1;
dno = -1;
order0 = 2:5;
E_cut = -1;
T = .1:.01:1;

load(join(['NLCE 2D Low Precision Data\muTgrid_m' num2str(m) '_u' num2str(double(u)) '_dno' num2str(dno) '_n' num2str(n)  '_nlce.mat']));
plotrange = 1:numel(T);
for i=1
    analytic = [];
    muz = muq(i);
    %for T = Tq0
        %adensity = @(kx,ky) m/(2*pi)^2*1./(exp(1/T*((-2*t)*(cos(kx)+cos(ky))-muz))+1);
        %aenergy = @(kx,ky) m/(2*pi)^2*1./(exp(1/T*((-2*t)*(cos(kx)+cos(ky))-muz))+1)*(-2*t).*(cos(kx)+cos(ky));
        %aentropy = @(kx,ky) m/(2*pi)^2.*(log( exp( ( 2*t.*(cos(kx)+cos(ky))+muz )/T) +1) - 1/T.*( 2*t.*(cos(kx)+cos(ky))+muz )./(1 + exp(1/T.*(-2*t.*(cos(kx)+cos(ky))-muz))));
        %analytic(end+1,:) = [integral2(adensity,-pi,pi,-pi,pi) integral2(aenergy,-pi,pi,-pi,pi) integral2(aentropy,-pi,pi,-pi,pi)];
    %end
    %subplot(2,3,1,'XScale','log'); hold on
    %for order = order0
     %   txt = [num2str(order)];
      %  array = csvread(join(['Densities_m' num2str(m) '_cut' num2str(E_cut) '_u' num2str(u) '_' num2str(order) '_dno' num2str(dno) '_n' num2str(n) '_nlce.csv']));
       % plot(T,array,'DisplayName',txt);    
        %ylabel('density');
        %xlabel('T');
        %ylim([0.99 1.01]);
    %end
    %writematrix(analytic(:,1),'Densities_non interacting 2d.csv');
    %legend('Location','SouthWest');
    %legend show
    %hold off
    subplot(2,4,1,'XScale','log');  hold on
    for order = order0
        txt = [num2str(order)];
        array = csvread(join(['Doubleoccupancies_m' num2str(m) '_cut' num2str(E_cut) '_u' num2str(u) '_' num2str(order) '_dno' num2str(dno) '_n' num2str(n)  '_nlce.csv']));
        plot(T,array,'-','DisplayName',txt);
        ylabel('Doubleoccupancies');
        xlabel('T');
        %xlim([0.3,0.4]);
        ylim([0. 0.05]);
    end
    legend('Location','NorthEast');
    legend show
    hold off
    
    subplot(2,4,2,'XScale','log');  hold on
    for order = order0
        txt = [num2str(order)];
        array1 = csvread(join(['Energies_m' num2str(m) '_cut' num2str(E_cut)...
            '_u' num2str(u) '_' num2str(order) '_dno' num2str(dno) '_n' num2str(n)  '_nlce.csv']));
        array2 = csvread(join(['Doubleoccupancies_m' num2str(m) '_cut' num2str(E_cut)...
            '_u' num2str(u) '_' num2str(order) '_dno' num2str(dno) '_n' num2str(n)  '_nlce.csv']));
        array = array1 - u*array1;
        plot(T,array,'-','DisplayName',txt);
        ylabel('kinetic energy');
        xlabel('T');
        %xlim([0.3,0.9]);
        %ylim([-0.3 -.2]);
    end
   % plot(Tq0,analytic(:,2),'DisplayName','Analytic');
    legend('Location','South');
    legend show
    hold off    
    
    subplot(2,4,3,'XScale','log');  hold on
    for order = order0
        txt = [num2str(order)];
        array = csvread(join(['Energies_m' num2str(m) '_cut' num2str(E_cut) '_u' num2str(u) '_' num2str(order) '_dno' num2str(dno) '_n' num2str(n)  '_nlce.csv']));
        plot(T,array,'-','DisplayName',txt);
        ylabel('energy');
        xlabel('T');
        %xlim([0.3,0.9]);
        ylim([-0.3 -.2]);
    end
   % plot(Tq0,analytic(:,2),'DisplayName','Analytic');
    legend('Location','South');
    legend show
    hold off
    
    subplot(2,4,4,'XScale','log');  hold on
    for order = order0
        txt = [num2str(order)];
        array = csvread(join(['Entropies_m' num2str(m) '_cut' num2str(E_cut) '_u'  num2str(u) '_' num2str(order) '_dno' num2str(dno) '_n' num2str(n) '_nlce.csv']));
        plot(T,array,'-','DisplayName',txt);
        ylabel('entropy');
        xlabel('T');
        %xlim([.3,.7]);
        ylim([0.6 1.4]);
    end
   % plot(Tq0,analytic(:,3),'DisplayName','Analytic');
    legend('Location','SouthEast');
    legend show
    hold off
    
    subplot(2,4,5,'XScale','log');  hold on
    for order = order0
        txt = [num2str(order)];
        arr = csvread(join(['Doubleoccupancies_m' num2str(m) '_cut'...
            num2str(E_cut) '_u' num2str(u) '_' num2str(order) '_dno' num2str(dno) '_n' num2str(n)  '_nlce.csv']));
        array = diff(arr)./diff(T);
        plot(T(1:end-1),array,'-','DisplayName',txt);
        ylabel('dD/dT');
        xlabel('T');
        %xlim([0.3,0.4]);
        ylim([-0.025 0.02]);
    end
    legend('Location','NorthEast');
    legend show
    hold off
    
    subplot(2,4,6,'XScale','log');  hold on
    for order = order0
        txt = [num2str(order)];
        array1 = csvread(join(['Energies_m' num2str(m) '_cut' num2str(E_cut)...
            '_u' num2str(u) '_' num2str(order) '_dno' num2str(dno) '_n' num2str(n)  '_nlce.csv']));
        array2 = csvread(join(['Doubleoccupancies_m' num2str(m) '_cut' num2str(E_cut)...
            '_u' num2str(u) '_' num2str(order) '_dno' num2str(dno) '_n' num2str(n)  '_nlce.csv']));
        arr = array1 - u*array1;
        array = midpt_diff(T,arr);
        plot(T(2:end-1),array,'-','DisplayName',txt);
        ylabel('dK/dT');
        xlabel('T');
        %xlim([0.3,0.9]);
        %ylim([-0.3 -.2]);
    end
   % plot(Tq0,analytic(:,2),'DisplayName','Analytic');
    legend('Location','South');
    legend show
    hold off   
    
    subplot(2,4,7,'XScale','log');  hold on
    for order = order0
        txt = [num2str(order)];
        arr = csvread(join(['Energies_m' num2str(m) '_cut' num2str(E_cut) '_u'  num2str(u) '_' num2str(order) '_dno' num2str(dno) '_n' num2str(n) '_nlce.csv']));
        array = midpt_diff(T,arr);
        plot(T(2:end-1),array,'-','DisplayName',txt);
        ylabel('specific heat');
        xlabel('T');
        %xlim([.1,1]);
        ylim([0. .3]);
    end
   % plot(Tq0,analytic(:,3),'DisplayName','Analytic');
    legend('Location','SouthEast');
    legend show
    hold off
    
    subplot(2,4,8,'XScale','log');  hold on
    for order = order0
        txt = [num2str(order)];
        arr = csvread(join(['Entropies_m' num2str(m) '_cut' num2str(E_cut) '_u'  num2str(u) '_' num2str(order) '_dno' num2str(dno) '_n' num2str(n) '_nlce.csv']));
        array = midpt_diff(T,arr);
        plot(T(2:end-1),array,'-','DisplayName',txt);
        ylabel('dS/dT');
        xlabel('T');
        %xlim([.3,.7]);
        ylim([-1 5]);
    end
   % plot(Tq0,analytic(:,3),'DisplayName','Analytic');
    legend('Location','SouthEast');
    legend show
    hold off    

end
 
sgtitle(join(['SU(' num2str(m) ') mu =' num2str(muq) ' u= ' num2str(u)]));
%subplot(2,3,5); plot(muq,cnns); title('CNN');