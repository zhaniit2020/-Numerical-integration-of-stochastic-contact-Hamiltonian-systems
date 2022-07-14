

randn('state',100)      
T=30;  N=200;  dt=0.1;  M=30;  

alpha=0.1;epsilon=0.02;

y=zeros(1,N);
t=zeros(1,N);
t(1,1)=0;
for i=2:N
    t(1,i)=t(1,i-1)+dt;
end


dW=zeros(M,N);
W=zeros(1,N);

dW=sqrt(dt)*randn(M,N);      
W(1,1)=dW(1,1);


 for i=2:N
     
        W(1,i)=W(1,i-1)+dW(1,i);
     
 end
dW1=zeros(1,N); 
dW1=mean(dW,1);

 x=zeros(1,N);
y=zeros(1,N);
z=zeros(1,N);   
x(1,1)=0.75;
y(1,1)=-0.25;
z(1,1)=0.08;



 for j=2:N     
     x(1,j)=x(1,j-1)+y(1,j-1)*dt;    
     y(1,j)=y(1,j-1)-x(1,j-1)*dt-alpha*z(1,j-1)*dt*(x(1,j)-x(1,j-1))+0.5*sin(x(1,j-1))*cos(x(1,j-1))*dt-cos(x(1,j-1))*dW1(1,j-1); 
     z(1,j)=z(1,j-1)+0.5*(x(1,j)-x(1,j-1))^2/dt-0.5*x(1,j-1)^2*dt-0.5*alpha*z(1,j-1)^2*dt-0.5*sin(x(1,j-1))*cos(x(1,j-1))*dt-sin(x(1,j-1))*dW1(1,j-1);
 end
 
x1=zeros(1,N);
y1=zeros(1,N);
z1=zeros(1,N);

x1(1,1)=0.75;
y1(1,1)=-0.25;
z1(1,1)=0.08;




 for j=2:N    
     
     x1(1,j)=x1(1,j-1)+dt*(1-0.5*dt*alpha*z1(1,j-1))*y1(1,j-1)-0.5*dt*dt*x1(1,j-1)-dt*cos(x1(1,j-1))*dW1(1,j-1)+0.5*cos(x1(1,j-1))*sin(x1(1,j-1))*dt;    
     y1(1,j)=((1-dt*alpha*z1(1,j-1))*y1(1,j-1)-0.5*x1(1,j-1)*dt-0.5*x1(1,j)*dt-cos(x1(1,j-1))*dW1(1,j-1)+0.5*dt*cos(x1(1,j-1))*sin(x1(1,j-1)))/(1+0.5*dt*alpha*z1(1,j)); 
     z1(1,j)=z1(1,j-1)+0.5*(x1(1,j)-x1(1,j-1))^2*1/dt-0.25*dt*(x1(1,j)^2+x1(1,j-1)^2)-0.25*alpha*dt*(z1(1,j)^2+z1(1,j-1)^2)-sin(x1(1,j-1))*dW1(1,j-1)-0.5*cos(x1(1,j-1))*sin(x1(1,j-1))*dt;
 end

 clear dW1;

plot(t,x,'-b');             
 xlabel('Time','FontSize',16)
 ylabel('Solutions','FontSize',16)
hold on

 plot(t,x1,'-r','LineWidth',0.8);
 xlabel('Time','FontSize',16)
 ylabel('Comparison of q','FontSize',16)
 legend('Reference:Euler-Maruyama scheme(25)','Contact scheme(24)') 
  hold on
  figure

 
plot(t,y,'-b');               
 xlabel('Time','FontSize',16)
 ylabel('Solutions','FontSize',16)
hold on

 plot(t,y1,'-r','LineWidth',0.8);
 xlabel('Time','FontSize',16)
 ylabel('Comparison of p','FontSize',16)
 legend('Reference:Euler-Maruyama scheme(25)','Contact scheme(24)') 
  hold on
  figure

plot(t,z,'-b');                
 xlabel('Time','FontSize',16)
hold on

  plot(t,z1,'-r','LineWidth',0.3);
  xlabel('Time','FontSize',16)
   ylabel('Comparison of s','FontSize',16)
 legend('Reference: Euler-Maruyama scheme (25)','By Contact scheme(24)') 
  
 hold on
 figure
% 
  plot3(x1,y1,z1,':');        
  xlabel('q','FontSize',16)
  ylabel('p','FontSize',16)
  zlabel('s','FontSize',16)
  hold on
  figure
 
 
%                            
 z2=zeros(1,N);
 z3=zeros(1,N);
 z3(1,1)=1;
 for i=1:N
     z2(1,i)=z1(1,i)-y(1,i)*x(1,i);
 end
 plot(t,z2,'-r');                  
 xlabel('Time','FontSize',16)
ylabel('ds-pdq','FontSize',16)
hold on
figure
%                            
  for i=2:N
     z3(1,i)=z2(1,i)/z(1,i-1);
 end
 plot(t,z3,'-b');                  
 xlabel('Time','FontSize',16)
ylabel('lamda_j','FontSize',16)
hold on

%
                          

z4=zeros(1,N);
for i=1:N
    z4(1,i)=exp(-alpha*dt);
end
plot(t,z4,'-r');                  
 xlabel('Time','FontSize',16)
ylabel('conformal factor','FontSize',16)
 legend('By contact scheme(24)','By continuous case') 
hold off

