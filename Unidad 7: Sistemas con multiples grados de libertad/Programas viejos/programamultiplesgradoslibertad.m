function multiplesgradosdelibertad
clc
clear
%Importante:
%Las siguientes variables se introducen manualmente: M,C,K,X0,dX0,tf,Ramort y F(i) (Fza que actua sobre sist)
%Nota: Al programa le falta el caso para cuando tenemos una fuerza aplicada
%sobre el sistema


#Primero debemos ingresar la matriz M, C y K de manera manual
N=input('Ingrese el numero de grados de libertad de su sistema: ');
M=[2 0 0; 0 2 0; 0 0 2]; #Modificar manualmente segun los grados de libertad que tenga
K=[3600 -1200 0; -1200 1800 -600; 0 -600 600]; #Idem
#Lo mismo para los vectores de condiciones iniciales X(0)=X0 y X'(0)=dX0
X0=[0.3;-0.8;0.3];
dX0=[0;0;0];
#Recordando que luego de un desarrollo se obtiene [K-w^2*M]*X=0, queremos encontrar los autovalores
#y autovectores de dicha funcion, para ello creamos la matriz modal X y la matriz diagonal
#que contendra a los valores propios VP. Luego obtenemos ambos por medio de la sig. func.
[Xn,VP]=eig(K,M);
#Luego de dicha funcion obtenemos la matriz modal de vectores propios NORMALIZADOS
#y una matriz diagonal que contiene los valores propios
VP=diag(VP); #Con esta funcion convertimos la matriz diagonal de VP en un vector columna
disp('La matriz modal NORMALIZADA es: ');
disp(Xn);
#Para obtener la matriz modal sin normalizar (como me dio analiticamente):
for (i=1:1:N)
X(:,i)=Xn(:,i)/Xn(1,i);
endfor
disp('La matriz modal SIN NORMALIZAR es: ');
disp(X);
#Luego los valores propios seran:
disp('Los valores propios son: ');
disp(VP);
#Luego, las frecuencias modales a partir de los autovalores seran:
disp('Las frecuencias modales son: ');
wn=sqrt(VP); #Creo el vector columna de frecuencias modales
disp(wn);
tf=((2*pi)/wn(1)); #Ingresamos manualmente el tiempo final hasta donde queremos calcular el despl.
#Ahora graficamos los modos
for (i=1:1:N+1)
  ejey(i)=i-1;
endfor
for (i=1:1:N)
figure(i);
hold
plot([0 0],[0,N],'k'); #Para dibujar linea negra al centro
plot([0;X(:,i)],ejey);
endfor
#Ahora lo que queremos obtener es la solucion a la ecuacion de movimiento de N grados de libertad
#Eso se logra mediante el metodo de descomposicion modal
#Podemos tener 2 casos, para un sistema SIN amortiguamiento o CON amortiguamiento
pregamort=input('El sistema es amortiguado? Ingrese Si o No: ','s');
if (pregamort=='No')
  disp('El sistema es NO AMORTIGUADO');
  #Pasamos a coordenadas modales
  qinicial=transpose(Xn)*M*X0;
  dqinicial=transpose(Xn)*M*dX0;
dt=input('Ingrese el paso de tiempo dt: ');
#Nuestro sistema puede o no estar forzado
  pregforz=input('El sistema esta forzado?.Ingrese Si o No: ','s');
  if (pregforz=='No')
    comp=1;
    for(t=0:dt:tf)
     F(:,comp)=diag(zeros(N));
     comp++;
     endfor
  elseif(pregforz=='Si')
  comp=1;
   for(t=0:dt:tf)
    F(:,comp)=[0; 0; 5000*sin(1.1*wn(1)*t)];#Armo manualmente la matriz de fuerza, una fila por cada G.L.
    comp++;
   endfor
  endif
#Hallamos el vector de carga modal
 Q=transpose(Xn)*F;
 disp(Q)
Aviejo=0; #Inicializamos el valor del acumulador A (Como en la teoria, admitimos A0=0)
Bviejo=0; #Inicializamos el valor del acumulador B (Como en la teoria, adminitos A0=0)
ycviejo=0; #Inicializamos el valor ycviejo
ysviejo=0; #Inicializamos el valor ysviejo
  for(i=1:1:N)
  comp=1;
   for(t=0:dt:tf)
    #Primero va la parte de Duhamel
    ycnuevo=Q(i,comp)*cos(wn(i)*t);
    ysnuevo=Q(i,comp)*sin(wn(i)*t);
    Anuevo=Aviejo+(dt/(2*M(i,i)*wn(i)))*(ycviejo+ycnuevo);
    Bnuevo=Bviejo+(dt/(2*M(i,i)*wn(i)))*(ysviejo+ysnuevo);
    #Tambien debemos calcular las componentes del vector X(tn), el me entrega el desplazamiento en forma discreta
    D(comp)=Anuevo*sin(wn(i)*t)-Bnuevo*cos(wn(i)*t);
    #Ahora calculamos las coordenadas modales qi
    q(i,comp)=qinicial(i)*cos(wn(i)*t)+(dqinicial(i)/wn(i))*sin(wn(i)*t)+D(comp);
    #Ahora debemos intercambias las variables viejas por las nuevas
    Aviejo=Anuevo;
    ycviejo=ycnuevo;
    ysviejo=ysnuevo;

    tiempo(comp)=t;
    comp++;
   endfor
endfor
#Ahora pasamos a coordenadas geometricas nuevamente para hallar los desplazamientos x(t);
x=Xn*q;
#Ahora, para un sistema CON amortiguamiento
elseif (pregamort=='Si')
  disp('El sistema es AMORTIGUADO')
  Ramort=[0.10;0.10;0.10]; #Debemos llenar este vector manualmente de acuerdo a los grados de libertad que tenga
  #Debo calcular la matriz diagonal de amortiguamiento modal CM
  for(i=1:1:N)
   wd(i)=wn(i)*sqrt(1-(Ramort(i))^2); #Creo el vector de frecuencia natural asociada al modo "i"
   vector(i,1)=2*Ramort(i)*wn(i); #Creamos un vector columna de elementos 2*Ramorti*Wni
  endfor
  CM=diag(vector); #Creo la matriz diagonal de amortiguamiento modal
  #Ahora puedo calcular la matriz de amortiguamiento proporcional C
  C=M*X*CM*transpose(Xn)*M;
  #Pasamos a coordenadas modales
  qinicial=transpose(Xn)*M*X0;
  dqinicial=transpose(Xn)*M*dX0;
  dt=input('Ingrese el paso de tiempo dt: ');
#Nuestro sistema puede o no estar forzado
  pregforz=input('El sistema esta forzado?.Ingrese Si o No: ','s');
  if (pregforz=='No')
    comp=1;
    for(t=0:dt:tf)
     F(:,comp)=diag(zeros(N));
     comp++;
     endfor
  elseif(pregforz=='Si')
  comp=1;
   for(t=0:dt:tf)
    F(:,comp)=[0; 0; 5000*sin(1.1*wn(1)*t)];#Armo manualmente la matriz de fuerza, una fila por cada G.L.
    comp++;
   endfor
  endif
#Hallamos el vector de carga modal
 Q=transpose(Xn)*F;
Aviejo=0; #Inicializamos el valor del acumulador A (Como en la teoria, admitimos A0=0)
Bviejo=0; #Inicializamos el valor del acumulador B (Como en la teoria, adminitos A0=0)
ycviejo=0; #Inicializamos el valor ycviejo
ysviejo=0; #Inicializamos el valor ysviejo
  for(i=1:1:N)
    comp=1;
   for(t=0:dt:tf)
#Primero va la parte de Duhamel
    ycnuevo=Q(i,comp)*cos(wd(i)*t);
    ysnuevo=Q(i,comp)*sin(wd(i)*t);
    Anuevo=Aviejo*exp(-Ramort(i)*wn(i)*dt)+(dt/(2*M(i,i)*wd(i)))*(ycviejo*exp(-Ramort(i)*wn(i)*dt)+ycnuevo);
    Bnuevo=Bviejo*exp(-Ramort(i)*wn(i)*dt)+(dt/(2*M(i,i)*wd(i)))*(ysviejo*exp(-Ramort(i)*wn(i)*dt)+ysnuevo);
    #Tambien debemos calcular las componentes del vector X(tn), el me entrega el desplazamiento en forma discreta
    D(comp)=Anuevo*sin(wd(i)*t)-Bnuevo*cos(wd(i)*t);
    #Ahora calculamos las coordenadas modales qi
    q(i,comp)=exp(-Ramort(i)*wn(i)*t)*(cos(wd(i)*t)+(Ramort(i)/sqrt(1-Ramort(i)))*sin(wd(i)*t))*qinicial(i)+((1/wd(i))*exp(-Ramort(i)*wn(i)*t)*sin(wd(i)*t))*dqinicial(i)+D(comp);
    #Ahora debemos intercambias las variables viejas por las nuevas
    Aviejo=Anuevo;
    ycviejo=ycnuevo;
    ysviejo=ysnuevo;

    tiempo(comp)=t;
    comp++;
    endfor
  endfor
#Ahora pasamos a coordenadas geometricas nuevamente para hallar los desplazamientos x(t);
x=Xn*q;
endif
#Ploteamos
figure(4)
for (i=1:1:N)
plot(tiempo,x(i,:));
hold on
endfor
hold off
endfunction

