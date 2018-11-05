clear all; close all; clc

A = load('saveMs1B.mat');
B = load('saveMs2B.mat');
C = load('saveMs3B.mat');
A = A.Ms1;
B = B.Ms2;
C = C.Ms3;

for i = 2:10
    for j = 1:i+1
        for u = 1:i
            if u<j
                figure(1)
                X = A(i,j)*(1+0.01*randn);
                semilogy(i+0.05*(u-0.5*i),X,'.b','MarkerSize',30); hold on
                semilogy(i+0.05*(u-0.5*i),X,'ow','MarkerSize',10); hold on
                figure(2)
                X = B(i,j)*(1+0.01*randn);
                semilogy(i+0.05*(u-0.5*i),X,'.b','MarkerSize',30); hold on
                semilogy(i+0.05*(u-0.5*i),X,'ow','MarkerSize',10); hold on
                figure(3)
                X = C(i,j)*(1+0.01*randn);
                semilogy(i+0.05*(u-0.5*i),X,'.b','MarkerSize',30); hold on
                semilogy(i+0.05*(u-0.5*i),X,'ow','MarkerSize',10); hold on
            else
                figure(1)
                X = A(i,j)*(1+0.01*randn);
                semilogy(i+0.05*(u-0.5*i),X,'.r','MarkerSize',30); hold on
                semilogy(i+0.05*(u-0.5*i),X,'ow','MarkerSize',10); hold on
                figure(2)
                X = B(i,j)*(1+0.01*randn);
                semilogy(i+0.05*(u-0.5*i),X,'.r','MarkerSize',30); hold on
                semilogy(i+0.05*(u-0.5*i),X,'ow','MarkerSize',10); hold on
                figure(3)
                X = C(i,j)*(1+0.01*randn);
                semilogy(i+0.05*(u-0.5*i),X,'.r','MarkerSize',30); hold on
                semilogy(i+0.05*(u-0.5*i),X,'ow','MarkerSize',10); hold on               
            end
        end
    end
end
figure(1)
loglog(linspace(1,11,100),ones(1,100),'--k','LineWidth',3); goodplot
set (gca,'color','none')
set(gcf,'units','points','position',[100,100,700,500])
figure(2)
loglog(linspace(1,11,100),ones(1,100),'--k','LineWidth',3); goodplot
set (gca,'color','none')
set(gcf,'units','points','position',[100,100,700,500])
figure(3)
loglog(linspace(1,11,100),ones(1,100),'--k','LineWidth',3); goodplot
set (gca,'color','none')
set(gcf,'units','points','position',[100,100,700,500])
