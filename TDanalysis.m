% Time Domain analysis from VNA Frecuency Domain measurements
% Author: Emilio José Martínez Pérez (Github user: emu0)

% Script para el análisis en el dominio del tiempo a partir de las medidas
% en el dominio de la frecuencia del puerto s21 de un analizador de redes.
% Por Emilio José Martínez Pérez.
%
% El objetivo de esta práctica es la de obtener la respuesta de un filtro
% pasabanda de microondas en el dominio del tiempo, tomando como punto de
% partida los datos tomados de un analizador de redes, en el dominio de la
% frecuencia
% 
% Para ello, los siguientes pasos han de ser llevados a cabo:
%                        _______________________________________________
%                       | 2: Generación espectro en frecuencias (FFT)   |
%         ____________  |  ______________   ________    ______________  |
%        |            | | |              | |        |  |              | |
% X(f) =>|1: Pre-proc.|_|_|2a: Colocación|_|2b: Zero|__|2c: Colocación|_|_
%        |   de datos | | |espectro izd. | |padding |  | espectro drc.| | |
%        |____________| | |______________| |________|  |______________| | |
%                       |_______________________________________________| |
%         ____________                                                    |
%        |            |                                                   |
% x(t) <=|  3:  IFFT  |___________________________________________________|
%        |            |
%        |____________|
%
% Paso 1: 
% Preprocesamiento de los datos obtenidos del analizador de redes para 
% facilitar los cálculos.
%
% Paso 2: 
% Generar un vector de la respuesta en frecuencia. Este vector tendrá que 
% ser igual a la salida que ofrece el algoritmo FFT de matlab.
% 
% Paso 3: 
% Procesamos con el algoritmo IFFT el vector con la respuesta en frecuencia
% del paso 2.
%
%Estos pasos serán explicados en detalle en su correspondiente apartado.

clear all
close all

%% PASO 1: Pre-procesamiento de los datos
% Primero se debemos realizar un pre-procesamiento de los datos obtenidos
% desde el analizador de redes. 
%
% Este script ha sido desarrollado para importar los datos de un fichero 
% *.dat y generar una variable de dato tipo struct. 
%
% La variable contendrá la información de amplitud y fase para las 
% frecuencias medidas.

[~,name] = fileparts('s21.dat');
s21.(genvarname(name)) = importdata('s21.dat');

vars = fieldnames(s21);
for i = 1:length(vars)
    assignin('base', vars{i}, s21.(vars{i}));
end

% A partir del tipo de dato struct creamos un vector complejo que nos será
% de utilidad durante el procedimiento que se desarrollará
s21_complex = zeros (length(s21.data),1);
s21_complex(:,1) = s21.data(:,2)+1i.*(s21.data(:,3));

% Ahora, para asegurar la generalización del código, calculamos la
% separación en frecuencias de dos puntos de esta señal discreta. Esta
% diferencia está directamente relacionada con la frecuencia de muestreo.

Af=s21.data(2,1)-s21.data(1,1);

% A partir de ahora se trabajará con posiciones del vector relativas a esta
% delta de f y no con los valores de frecuencia medidos.

% Podremos calcular la posición f1 y f2 correspondientes con la frecuencia 
% de inicio y final en las medidas del analizador de redes.

f1=s21.data(1,1)/Af;
f2=s21.data(end,1)/Af;

% Estimamos la frecuencia de muestreo que querremos en 4 veces el último
% valor de la frecuencia medida.
Fs=4*f2;

% Tras la finalización del primer paso, podemos observar los datos del
% analizador de redes.

figure(1);
subplot(1, 3, 1); plot(s21_complex);
title('s21');
subplot(1, 3, 2); plot (s21.data(:,1), sqrt(real(s21_complex).^2 + imag(s21_complex).^2));
title('s21 modulo'); xlabel('Frecuencia'); ylabel('|s21|');
subplot(1, 3, 3); plot (s21.data(:,1), atan(imag(s21_complex)./real(s21_complex)));
title('s21 fase'); xlabel('Frecuencia'); ylabel('Radianes');


%% PASO 2: Generación espectro en frecuencias (FFT)
% En este paso se deberá crear un vector que tenga la misma forma que el
% que obtendríamos si procesamos una señal en el dominio del tiempo con 
% el algoritmo FFT.
%
% El vector con los valores de la FFT estará compuesto por el espectro en
% frecuencias obtenido en el analizador de redes, que estará situado en la
% frecuencia medida, lo que llamaremos la parte izquierda del espectro, y 
% luego los mismos datos pero reflejados en el eje horizontal y situados en
% la zona de alta frecuencia, lo que llamamos la parte derecha del
% espectro.
% Para la construcción de este vector, hemos dividido el proceso en tres
% pasos:
%

%%% 2a: Colocación de la parte izquierda del espectro
% En este paso generamos la parte izquierda del espectro en frecuencias.
% Primero sabemos que el espectro tendrá tantas muestras como frecuencia de
% muestreo hayamos escogido.
%
% Tomamos esta parte izquierda como el espectro en frecuencia que tome
% valores desde f=0 hasta la mitad de la frecuencia de muestreo, Fs/2.
vector_FFT_izquierda=zeros(Fs/2+1, 1);
% Añadimos los datos obtenidos del analizador de redes.
vector_FFT_izquierda(f1+1:f2+1, 1) = s21_complex;

%%% 2b: Zero padding
% Este paso es necesario para obtener una buena resolución en tiempo del
% vector que obtendremos. 
% Basta con incluir ceros en el vector de FFT en valores de frecuencia 
% alta, entre las dos partes que ofrecen información dentro del espectro 
% en frecuencia.
%
% Para obtener un buen muestreo de la señal, es necesario realizar un
% sobremuestreo a una frecuencia mayor que la frecuencia de Nyquist. En
% nuestro caso hemos tomado la Fs como cuatro veces la mayor frecuencia
% medida desde el analizador de redes, en nuestro caso particular, 
% 4*3GHz=12GHz
%
% Este es un truco matemático que, como se ha indicado anteriormente,
% permite obtener un resultado con más resolución en el dominio del tiempo
% y es además, un paso que está implícito en el paso 2a y 2c, pero que es 
% conveniente explicar su utilidad.

%%% 2c: Colocación de la parte derecha del espectro
% La parte derecha del espectro se trata del espectro tomado en el
% analizador de redes, salvo con la diferencia de que será el vector
% hérmítico del mismo. Es decir, la parte derecha del espectro se formará  
% a partir del vector rotado y conjugado complejo de la parte izquierda del
% espectro. No se debe incluir el valor 0 del espectro.
vector_FFT_derecha=conj(flipud(vector_FFT_izquierda(2:end)));

% El espectro en frecuencia final será la unión de ambos vectores.
vector_FFT=[vector_FFT_izquierda; vector_FFT_derecha];

% Se puede seguir el proceso seguido con las siguientes imágenes, tomando
% como ejemplo, la parte real del espectro.
figure(2);
subplot(2, 3, 1); plot(real(s21_complex)); title('Datos importados del VNA');
subplot(2, 3, 2); plot(real(vector_FFT_izquierda)); title('Parte izquierda');
subplot(2, 3, 3); plot(real(vector_FFT_derecha)); title('Parte derecha');
subplot(2, 3, 4:6); plot((0:Af:Fs*Af),real(vector_FFT)); 
title('Vector FFT'); xlabel('Frecuencia')


%% PASO 3: IFFT
% Una vez obtenido el vector del espectro en frecuencia, solamente
% tendremos que realizar la transformada de Fourier inversa para obtener la
% señal en el dominio del tiempo.
vector_t=ifft(vector_FFT);
% Para conocer si se ha realizado correctamente todos los procesos, hemos
% de observar el vector obtenido. Si es un vector real, entonces es que se
% ha procedido correctamente.

% Podemos observar el vector obtenido
figure(3);
plot((0:1/(Fs*Af):1/Af), vector_t);
title('Dominio del tiempo'); xlabel('tiempo'); ylabel('V normalizado')

% En este vector de tiempos se observa la respuesta al impulso en el
% dominio del tiempo del filtro pasabanda analizado. 