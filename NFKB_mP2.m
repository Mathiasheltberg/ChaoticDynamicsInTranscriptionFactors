function dydt = NFKB_mP2(t,y,A,N1,N2)
Vol = 2.0*10^(-14); NA = 6.02*10^(23); Cal = NA*Vol*10^(-6);

kNin = 5.4; Ntot = 1*Cal; KN = 0.029*Cal; kIin = 0.018;
kt = 1.03/Cal; ym = 0.017;
ktl = 0.24; a = 1.05/Cal; KI = 0.035*Cal;
ka = 0.24; IKKtot = 2*Cal; ki = 0.18;
kp = 0.036; kA20 = 0.0018; A20 = 0.0026;
T = 50.0;
Ntot = 1*Cal; KN = 0.029*Cal; kt = 1.03/Cal;  a = 1.05/Cal; KI = 0.035*Cal;   IKKtot = 2*Cal;
h1 = 2.0; K1 = 1.0; del1 = 0.01; rhoA1 = 1.0; RhoA1 = 1.0; Del1 = 0.001;
h2 = 4.0; K2 = 4.5; del2 = 0.01; rhoA2 = 5.0; RhoA2 = 5.0; Del2 = 0.001;
LamN = 1.0/Cal^(N1+N2-1); DelN = 0.001;

 %*(RTot - y(7) - y(8))
 %(RTot - 1*y(9) - 1*y(10))^np1/((1*y(9) + 1*y(10))^np1 + (RTot - 1*y(9) - 1*y(10))^np1)
dydt = [kNin*( Ntot-y(1) ) * KI/(KI + y(3)) - kIin*y(3)*y(1)/(KN + y(1));
    kt*y(1)*y(1) - ym*y(2);
    ktl*y(2) - a*y(4)*( Ntot - y(1) )*y(3)/(y(3) + KI);
    ka*(0.5 + A*sin(2*pi/T*y(6)))*(IKKtot - y(4) - y(5)) - ki*y(4)
    ki*y(4) - kp*y(5)*kA20/(kA20 + A20*(0.5 + A*sin(2*pi/T*y(6))));
    1;    
    ];
    PP = 1;
    for i = 1:N1+N2;  PP = PP*y(6+2*i); end
    for i = 1:N1
        dydt = [dydt; rhoA1*(y(1)/(2000*6.023))^h1/((y(1)/(2000*6.023))^h1 + (K1/6.023)^h1)-del1*y(5 + 2*i);];
        dydt = [dydt; RhoA1*y(5+2*i) - LamN*PP - Del1*y(6+2*i);];
    end 
    for i = 1:N2
        dydt = [dydt;  rhoA2*(y(1)/(2000*6.023))^h2/((y(1)/(2000*6.023))^h2 + (K2/6.023)^h2) - del2*y(5 + 2*N1 + 2*i);];
        dydt = [dydt;  RhoA2*y(5 + 2*N1 + 2*i) - LamN*PP - Del2*y(6 + 2*N1 + 2*i);];
    end

    PP2 = PP*LamN - DelN*y(6+2*(N1+N2)+1);
    dydt = [dydt;  PP2];
    