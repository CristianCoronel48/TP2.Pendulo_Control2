% Caso de estudio 2. Sistema no lineal de cuatro variables de estado
% Con las variables de estado y los valores de los coeficientes.
% Determinar deltaT y el tiempo de simulación adecuados
clc;clear all;
m=0.01;Fricc=0.1;long=1.2;g=9.8;M=0.5;h=1e-5;
tF=10;tiempo=(tF/h);t=0:h:tiempo*h;
% Variables no lineales
del_pp=0;
phi_pp=0;
del_p=0:h:tiempo*h;
del=0:h:tiempo*h;
phi_p=0:h:tiempo*h;
phi=0:h:tiempo*h;
% Variables lineales
fi=0:h:tiempo*h;
fi_p=0:h:tiempo*h;
delta=0:h:tiempo*h;
delta_p=0:h:tiempo*h;
u=linspace(0,0,tiempo+2);
%Condiciones iniciales
fi(1)=-0.01; phi(1)=-0.01; color='b'; i=1;
%Versión linealizada en el equilibrio estable
Mat_A=[0 1 0 0;0 -Fricc/M -g*m/M 0; 0 0 0 1;0 Fricc/(M*long) g*(M+m)/(M*long) 0];
Mat_B=[0; 1/M; 0; -1/(long*M)];
X0=[0 0 0 0]';x=[0 0 fi(1) 0]';
while(i<(tiempo+1))
% Sistema no lineal   
 u(i)=0; 
 del_pp=(1/(M+m))*(u(i)-m*long*phi_pp*cos(phi(i))+m*long*phi_p(i)^2*sin(phi(i))- Fricc*del_p(i));
 phi_pp=(1/long)*(g*sin(phi(i))-del_pp*cos(phi(i)));
 del_p(i+1)=del_p(i)+h*del_pp;
 del(i+1)=del(i)+h*del_p(i);
 phi_p(i+1)=phi_p(i)+h*phi_pp;
 phi(i+1)=phi(i)+h*phi_p(i);
%Sistema lineal
xp=Mat_A*(x-X0)+Mat_B*u(i);
x=x+h*xp;
i=i+1;
%Variables del sistema lineal
delta(i)=x(1); delta_p(i)=x(2);fi(i)=x(3);fi_p(i)=x(4);
end
figure(1);hold on;
subplot(2,2,1);plot(t,fi_p,color);grid on;hold on;
plot(t,phi_p,'r');title('Velocidad Ángulo');
axis([0 tF -10 10]);
subplot(2,2,2);plot(t,fi,color);hold on;
plot(t,phi,'r');hold on;grid on;title('Ángulo');
axis([0 tF -8 1]);
subplot(2,2,3);plot(t,delta,color);hold on;
plot(t,del,'r');grid on;title('Posición carro');
axis([0 tF -.5 .5]);
subplot(2,2,4);plot(t,delta_p,color);hold on;
plot(t,del_p,'r');grid on;title('Velocidad carro');
axis([0 tF -2 2]);
