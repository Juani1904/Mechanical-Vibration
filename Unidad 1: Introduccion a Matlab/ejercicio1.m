function ejercicio1
  clc
  producto=1;
  suma=0;
  for enterosPositivos=1:1:50
    producto=producto*enterosPositivos;
    suma=suma+enterosPositivos;
  endfor
  disp("El producto es: "), disp(producto);
  disp("La suma es: "), disp(suma);

endfunction
