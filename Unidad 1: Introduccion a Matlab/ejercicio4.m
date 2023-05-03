function ejercicio4
  clc
  #Primero definimos el vector de coordenadas en X
  coordenadasX=-30:0.5:30;
  #Ahora definimos el vector de coordenadas en Y
  coordenadasY=[];
  index=1;
  for x=-30:0.5:30
    if x < -2
      coordenadasY(index)=2;
    elseif (x >= -2) && (x < 3)
      coordenadasY(index)=x^2;
    elseif (x >= 3) && (x <= 10)
      coordenadasY(index)=1/x;
    elseif (x > 15) && (x < 20)
      coordenadasY(index)=x;
    elseif (x==22)
      coordenadasY(index)=3-x;
    else
      coordenadasY(index)=0;
    endif
    index++;
  endfor

  #Ahora debemos plotear
  hold on
  title("Grafica 1")
  xlabel("x")
  ylabel("f(x)")
  #legend("Grafica 1")
  plot(coordenadasX, coordenadasY)
  hold off


endfunction
