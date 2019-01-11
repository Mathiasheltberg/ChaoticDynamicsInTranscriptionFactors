clear all; close all; clc
%%%% This creates figures like the ones from Fig 3
Ms1 = zeros(10,11);
Ms2 = zeros(10,11);
Ms3 = zeros(10,11);

for NK = 2:10
    NK
    N1 = 0; N2 = 0;
%    for test = 2:2
     for test = 1:NK+1
        N1 = NK-N2;
        N2 = NK-N1;
        MsE = zeros(2,3);
        
        for i = 1:2
            i
            A = 0.10 + 0.30*(i-1);
            tspan = [0 100000];
            %y0 = 10000*[0.0886; 0.2625; 4.1076 ; 0.2167; 1.8662; 0.0];
            y0 = 10000*[0.0886; 0.4625; 4.1076 ; 0.2167; 1.8662; 0.0];
            for u = 1:2*N1; y0 = [y0; 10.0;]; end
            for u = 1:2*N2; y0 = [y0; 10.0;]; end
            y0 = [y0; 10.0;];
            [t,y] = ode45(@(t,y) NFKB_mP2(t,y,A,N1,N2), tspan, y0);
           % Ec = y(:,8)+y(:,10)+y(:,12)+y(:,14);
%             Ec = y(:,8)+y(:,10);
%              plot(t(1:10:end),y(1:10:end,6+2*(N1+N2)+1),'color',[(0.3+0.6*(i-1)) (0.4 + 0.2*i) 0.3],'LineWidth',3); hold on;
%             goodplot             
            MsE(i,1) = mean(y(round(end/2):end,6+2*(N1+N2)+1));
            MsE(i,2) = mean(y(round(end/2):end,6+2*(N1+N2)+1))/mean(y(round(end/2):end,1));
            Eco = 0;
            for ie = 8:2:(6+2*(N1+N2))
                ie
                Eco = Eco + mean(y(round(end/2):end,ie));
            end
            MsE(i,3) = mean(y(round(end/2):end,6+2*(N1+N2)+1))/Eco;
        end
        Ms1(NK,test) = MsE(2,1)/MsE(1,1);
        Ms2(NK,test) = MsE(2,2)/MsE(1,2);
        Ms3(NK,test) = MsE(2,3)/MsE(1,3);
        N2 = N2+1;
    end    
end
