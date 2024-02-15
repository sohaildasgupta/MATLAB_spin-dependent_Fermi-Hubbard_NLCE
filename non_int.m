t = 1; u = 15.3;
m = 6;
n = 6;
dno = 3;


load(join(['NLCE 2D Low Precision Data\muTgrid_m' num2str(m) '_u' num2str(double(u)) '_dno' num2str(dno) '_n' num2str(n)  '_nlce.mat']));

analytic = [];    
for i=1:numel(T)
    Tq = T(i);
    muz = mu_array(i);
    adensity = @(kx,ky) m/(2*pi)^2*1./(exp(1/Tq*((-2*t)*(cos(kx)+cos(ky))-muz))+1);
    aenergy = @(kx,ky) m/(2*pi)^2*1./(exp(1/Tq*((-2*t)*(cos(kx)+cos(ky))-muz))+1)*(-2*t).*(cos(kx)+cos(ky));
    aentropy = @(kx,ky) m/(2*pi)^2.*(log( exp( ( 2*t.*(cos(kx)+cos(ky))+muz )/Tq) +1) - 1/Tq.*( 2*t.*(cos(kx)+cos(ky))+muz )./(1 + exp(1/Tq.*(-2*t.*(cos(kx)+cos(ky))-muz))));
    analytic(end+1,:) = [integral2(adensity,-pi,pi,-pi,pi) integral2(aenergy,-pi,pi,-pi,pi) integral2(aentropy,-pi,pi,-pi,pi)];
end

hold on;
plot(T,analytic(:,3));
hold off;
