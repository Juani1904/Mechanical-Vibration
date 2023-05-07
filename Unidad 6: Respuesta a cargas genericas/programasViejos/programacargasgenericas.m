function programacargasgenericas
clc
clear
#Primero ingresamos los datos del problema
ti=input("Ingrese el tiempo inicial del que partiremos (en segundos): "); #Este seria el Tau
tf=input("Ingrese hasta que tiempo vamos a medir desplazamiento (en segundos): ");
dt=input("Ingrese el paso de muestreo de la fuerza: ");
#Creamos el vector F(Tau) que ira en el yc y el ys
t=0; #contador de tiempo
comp=1; #Contador de componentes del vector (Empieza en 2 porque la componente 1 ya la defini
for t=ti:dt:tf
  Ftau(comp)=input(['Ingrese valor de fuerza en tiempo ',num2str(t),' segundos: ']);
 comp++;  
endfor
disp(Ftau);
#La frecuencia natural puede ser dato o no. Para eso preguntamos
pregfrec=input("Posee el dato de la frecuencia natural? Ingrese Si o No: ",'s');
if(pregfrec=='Si')
 wn=input("Ingrese frecuencia natural: ");
elseif (pregfrec=='No')
 m=input("Ingrese el valor de la masa: ");
 k=input("Ingrese el valor de la constante de rigidez k: ");
 wn=sqrt(k/m);
endif
#El sistema podria o no estar amortiguado. Por eso preguntamos
pregamort=input("El sistema esta amortiguado? Ingrese Si o No: ",'s');
if (pregamort=='Si')
  disp("El sistema es amortiguado");
#El sistema podria tener como dato la constante c o directamente la relacion de amortiguamiento. Por eso preguntamos
  pregRela=input("Conoce la relacion de amortiguamiento? Ingrese Si o No: ",'s');
    if (pregRela=='Si')
      Ramort=input("Ingrese relacion de amortiguamiento Zeta: ");
    elseif (pregRela=='No')
      c=input("Ingrese el valor de la constante de amortiguamiento c: ");
      Ramort=c/(2*m*wn);
    endif
wd=sqrt(1-(Ramort)^2)*wn;    
#Ahora debemos calcular el An y Bn para el sistema amortiguado.
Aviejo=0; #Inicializamos el valor del acumulador A (Como en la teoria, admitimos A0=0)
Bviejo=0; #Inicializamos el valor del acumulador B (Como en la teoria, adminitos A0=0)
t=0; #contador de tiempo
  if (ti==0)
    ycviejo=0;
    ysviejo=0;
  else
    Fvieja=input(['Ingrese el valor de la fuerza en el instante ',num2str(ti-dt),' segundos']);
    ycviejo=Fvieja*cos(wd*(ti-dt))
    ysviejo=Fvieja*sin(wd*(ti-dt))
  endif
  comp=1;
#Para resolver el problema podemos hacerlo por regla del trapecio o formulacion incremental paso a paso
pregmetodo=input('Ingrese metodo que quiere utilizar. Escriba Regla del trapecio(RT) o Incremental Paso a Paso o Suma simple (SS): ','s');
if (pregmetodo=='RT')    
  for (t=ti:dt:tf)
    ycnuevo=Ftau(comp)*cos(wd*t);
    ysnuevo=Ftau(comp)*sin(wd*t);
    Anuevo=Aviejo*exp(-Ramort*wn*dt)+(dt/(2*m*wd))*(ycviejo*exp(-Ramort*wn*dt)+ycnuevo)
    Bnuevo=Bviejo*exp(-Ramort*wn*dt)+(dt/(2*m*wd))*(ysviejo*exp(-Ramort*wn*dt)+ysnuevo)
#Tambien debemos calcular las componentes del vector X(tn), el me entrega el desplazamiento en forma discreta
    X(comp)=Anuevo*sin(wd*t)-Bnuevo*cos(wd*t)
#Ahora debemos intercambias las variables viejas por las nuevas
    Aviejo=Anuevo;
    ycviejo=ycnuevo;
    ysviejo=ysnuevo;
    tiempo(comp)=t #Variable creada para luego plotear el eje de abscisas  
    comp++; #Incrementamos en 1 el contador de componentes   
  endfor
elseif(pregmetodo=='SS')
    for (t=ti:dt:tf)
    ycnuevo=Ftau(comp)*cos(wd*t);
    ysnuevo=Ftau(comp)*sin(wd*t);
    Anuevo=Aviejo*exp(-Ramort*wn*dt)+(dt/(2*m*wd))*(ycviejo*exp(-Ramort*wn*dt));
    Bnuevo=Bviejo*exp(-Ramort*wn*dt)+(dt/(2*m*wd))*(ysviejo*exp(-Ramort*wn*dt));
#Tambien debemos calcular las componentes del vector X(tn), el me entrega el desplazamiento en forma discreta
    X(comp)=Anuevo*sin(wd*t)-Bnuevo*cos(wd*t)
#Ahora debemos intercambias las variables viejas por las nuevas
    Aviejo=Anuevo;
    ycviejo=ycnuevo;
    ysviejo=ysnuevo;
    tiempo(comp)=t #Variable creada para luego plotear el eje de abscisas  
    comp++; #Incrementamos en 1 el contador de componentes   
  endfor
endif  
  
#Ahora lo mismo para el caso no amortiguado
elseif (pregamort=='No')
disp("El sistema no esta amortiguado "); 
  #Debemos calcular el An y Bn para el sistema NO amortiguado.
Aviejo=0; #Inicializamos el valor del acumulador A (Como en la teoria, admitimos A0=0)
Bviejo=0; #Inicializamos el valor del acumulador B (Como en la teoria, adminitos A0=0)
t=0; #contador de tiempo
  if (ti==0)
    ycviejo=0;
    ysviejo=0;
  else
    Fvieja=input(['Ingrese el valor de la fuerza en el instante ',num2str(ti-dt),' segundos']);
    ycviejo=Fvieja*cos(wn*(ti-dt));
    ysviejo=Fvieja*sin(wn*(ti-dt));
  endif
  comp=1;
#Para resolver el problema podemos hacerlo por regla del trapecio o formulacion incremental paso a paso
pregmetodo=input('Ingrese metodo que quiere utilizar. Escriba Regla del trapecio(RT) o Incremental Paso a Paso o Suma simple (SS): ','s');
if (pregmetodo=='RT') 
  for (t=ti:dt:tf)
    ycnuevo=Ftau(comp)*cos(wn*t); 
    ysnuevo=Ftau(comp)*sin(wn*t);
    Anuevo=Aviejo+(dt/(2*m*wn))*(ycviejo+ycnuevo);
    Bnuevo=Bviejo+(dt/(2*m*wn))*(ysviejo+ysnuevo);
#Tambien debemos calcular las componentes del vector X(tn), el me entrega el desplazamiento en forma discreta
    X(comp)=Anuevo*sin(wn*t)-Bnuevo*cos(wn*t)
#Ahora debemos intercambias las variables viejas por las nuevas
    Aviejo=Anuevo;
    ycviejo=ycnuevo;
    ysviejo=ysnuevo;
    tiempo(comp)=t #Variable creada para luego plotear el eje de abscisas
    comp++; #Incrementamos en 1 el contador de componentes  
  endfor
elseif (pregmetodo=='SS')
for (t=ti:dt:tf)
    ycnuevo=Ftau(comp)*cos(wn*t);
    ysnuevo=Ftau(comp)*sin(wn*t);
    Anuevo=Aviejo+(dt/(2*m*wn))*(ycviejo);
    Bnuevo=Bviejo+(dt/(2*m*wn))*(ysviejo);
#Tambien debemos calcular las componentes del vector X(tn), el me entrega el desplazamiento en forma discreta
    X(comp)=Anuevo*sin(wn*t)-Bnuevo*cos(wn*t)
#Ahora debemos intercambias las variables viejas por las nuevas
    Aviejo=Anuevo;
    ycviejo=ycnuevo;
    ysviejo=ysnuevo;
    tiempo(comp)=t #Variable creada para luego plotear el eje de abscisas
    comp++; #Incrementamos en 1 el contador de componentes  
  endfor  
endif
endif
#Finalmente graficamos la carga y la respuesta con respecto al tiempo y entregamos las respuestas
figure(1)
plot(tiempo,Ftau,'b')
figure(2)
plot(tiempo,X,'r')
grid on
endfunction