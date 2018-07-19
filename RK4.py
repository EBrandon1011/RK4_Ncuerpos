# -*- coding: utf-8 -*-
"""
Subido el 19 de julio de 2018

@author: Eddy Brandon de León
"""
#El origen se fija en el baricentro del Sistema Solar

import numpy as np
import time

AU = 149597870700    # Definiendo una UA
ao=1.2*10**(-10)         #Constante a_o

class Mond():
    nombre='ninguno'    #Será dado después
    masa=0.0
    x=y=z=0.0
    vx=vy=vz=0.0
    r=0.0
    
    def acc(self, cuerpok, kxi, kyi, kzi, kxk, kyk, kzk, h):
        Dxki = (self.x+h*kxi/2) - (cuerpok.x+h*kxk/2)
        Dyki = (self.y+h*kyi/2) - (cuerpok.y+h*kyk/2)
        Dzki = (self.z+h*kzi/2) - (cuerpok.z+h*kzk/2)
        rki = np.sqrt(Dxki**2+Dyki**2+Dzki**2)   
        Gmk = cuerpok.masa          #G*M_k
        lmk = np.sqrt(Gmk/ao)       #Longitud de masa
        Xki = lmk/rki               #Cantidad adimensional
        #f = Xki**2			    #Función Newton
        f = Xki*(1+Xki+Xki**2+Xki**3)/(1+Xki+Xki**2) 	#Función MOND 
        ax = -ao*f*Dxki/rki         #Componente x de la aceleración
        ay = -ao*f*Dyki/rki         #Componente y de la aceleración
        az = -ao*f*Dzki/rki         #Componente z de la aceleración
        return ax, ay, az
        
        
def iteraciones(astros):
    counter=0
    hh=1                #h = 1d/hh
    h=24*3600/hh        #Tamaño de paso: 1 día
    NI=20000000*hh      #Número de iteraciones        
    arx=open("RK2016R2.dat", "w")        #Archivo de texto de salida
    
    while counter<NI:       #Método de Runge-Kutta de cuarto orden para un sistema 6N ecuaciones
        k1v={}
        k1p={}
        k2v={}
        k2p={}
        k3v={}
        k3p={}
        k4v={}
        k4p={}
                                        
        for cuerpoi in astros:      #Cálculo de las k1 (6N veces)
            axi=0
            ayi=0
            azi=0
            k0xi = k0yi = k0zi = 0
            for cuerpok in astros:   
                k0xk = k0yk = k0zk = 0
                if cuerpok is cuerpoi:
                    continue
                k0xk = k0yk = k0zk = 0
                axki, ayki, azki=cuerpoi.acc(cuerpok, k0xi, k0yi, k0zi, k0xk, k0yk, k0zk, h)
                axi+=axki
                ayi+=ayki
                azi+=azki
            k1p[cuerpoi]=(cuerpoi.vx, cuerpoi.vy, cuerpoi.vz)
            k1v[cuerpoi]=(axi, ayi, azi)
            
        for cuerpoi in astros:      #Cálculo de las k2 (6N veces)
            axi=0
            ayi=0
            azi=0
            k1xi, k1yi, k1zi = k1p[cuerpoi]
            k1vx, k1vy, k1vz = k1v[cuerpoi]
            for cuerpok in astros:      
                k1xk, k1yk, k1zk = k1p[cuerpok]
                if cuerpok is cuerpoi:
                    continue
                axki, ayki, azki=cuerpoi.acc(cuerpok, k1xi, k1yi, k1zi, k1xk, k1yk, k1zk, h)
                axi+=axki
                ayi+=ayki
                azi+=azki
            k2p[cuerpoi]=(cuerpoi.vx+k1vx*h/2, cuerpoi.vy+k1vy*h/2, cuerpoi.vz+k1vz*h/2)
            k2v[cuerpoi]=(axi, ayi, azi)
            
        for cuerpoi in astros:      #Cálculo de las k3 (6N veces)
            axi=0
            ayi=0
            azi=0
            k2xi, k2yi, k2zi = k2p[cuerpoi]
            k2vx, k2vy, k2vz = k2v[cuerpoi]
            for cuerpok in astros:    
                k2xk, k2yk, k2zk = k2p[cuerpok]
                if cuerpok is cuerpoi: 
                    continue
                axki, ayki, azki=cuerpoi.acc(cuerpok, k2xi, k2yi, k2zi, k2xk, k2yk, k2zk, h)
                axi+=axki
                ayi+=ayki
                azi+=azki
            k3p[cuerpoi]=(cuerpoi.vx+k2vx*h/2, cuerpoi.vy+k2vy*h/2, cuerpoi.vz+k2vz*h/2)
            k3v[cuerpoi]=(axi, ayi, azi)
            
        for cuerpoi in astros:      #Cálculo de las k4 (6N veces)
            axi=0
            ayi=0
            azi=0
            k3xi, k3yi, k3zi = k3p[cuerpoi]
            k3vx, k3vy, k3vz = k3v[cuerpoi]
            for cuerpok in astros:     
                k3xk, k3yk, k3zk = k3p[cuerpok]
                if cuerpok is cuerpoi:
                    continue
                axki, ayki, azki=cuerpoi.acc(cuerpok, 2*k3xi, 2*k3yi, 2*k3zi, 2*k3xk, 2*k3yk, 2*k3zk, h)
                axi+=axki
                ayi+=ayki
                azi+=azki
            k4p[cuerpoi]=(cuerpoi.vx+k3vx*h, cuerpoi.vy+k3vy*h, cuerpoi.vz+k3vz*h)
            k4v[cuerpoi]=(axi, ayi, azi)
                
        for cuerpoi in astros:     #Calculando las u_(i+1), donde u es x, y, z, vx, vy, vz
            k1vx, k1vy, k1vz = k1v[cuerpoi]
            k1x, k1y, k1z = k1p[cuerpoi]
            k2vx, k2vy, k2vz = k2v[cuerpoi]
            k2x, k2y, k2z = k2p[cuerpoi]
            k3vx, k3vy, k3vz = k3v[cuerpoi]
            k3x, k3y, k3z = k3p[cuerpoi]
            k4vx, k4vy, k4vz = k4v[cuerpoi]
            k4x, k4y, k4z = k4p[cuerpoi]
            
            if cuerpoi.nombre is'16R2':     #Archivo de texto de salida. Solo se genera para el asteroide/cometa en estudio.
                escribir='{:>11.15f} {:>11.15f} {:>11.15f}'.format(cuerpoi.x/AU, cuerpoi.y/AU, cuerpoi.z/AU) #Datos x, y, z
                arx.write(escribir+"\n")

            cuerpoi.vx += (k1vx+k2vx*2+k3vx*2+k4vx)*h/6
            cuerpoi.vy += (k1vy+k2vy*2+k3vy*2+k4vy)*h/6
            cuerpoi.vz += (k1vz+k2vz*2+k3vz*2+k4vz)*h/6
            cuerpoi.x += (k1x+k2x*2+k3x*2+k4x)*h/6
            cuerpoi.y += (k1y+k2y*2+k3y*2+k4y)*h/6
            cuerpoi.z += (k1z+k2z*2+k3z*2+k4z)*h/6          

        arx.closed
        counter+=1

#Lee las condiciones iniciales de cada planeta. Estos datos están en archivos de texto.           
def condiciones(astros):
    for cuerpoi in astros:
        file=open('datos{}.dat'.format(cuerpoi.nombre), 'r')
        cuerpoi.masa = float(file.readline())*1000**3  	#GM [km**3/s**2]
        cuerpoi.x = float(file.readline())*1000            #X [km]
        cuerpoi.y = float(file.readline())*1000
        cuerpoi.z = float(file.readline())*1000
        cuerpoi.vx = float(file.readline())*1000           #Vx [km/s]
        cuerpoi.vy = float(file.readline())*1000
        cuerpoi.vz = float(file.readline())*1000
        file.close
                        
def main():
    start = time.time()     #Inicio del tiempo de ejecución
    
    sol = Mond()        #Se indican las coordenas del baricentro de cada cuerpo [metros].
    sol.nombre = 'Sol'

    mercurio = Mond()
    mercurio.nombre = 'Mercurio'

    venus = Mond()
    venus.nombre = 'Venus'

    tierra = Mond()
    tierra.nombre = 'Tierra'

    marte = Mond()
    marte.nombre = 'Marte'
  
    jupiter = Mond()
    jupiter.nombre = 'Jupiter'
 
    saturno = Mond()
    saturno.nombre = 'Saturno'

    urano = Mond()
    urano.nombre = 'Urano'
    
    neptuno = Mond()
    neptuno.nombre = 'Neptuno'
    
    pluton = Mond()
    pluton.nombre = 'Pluton'
    
    eris = Mond()
    eris.nombre = 'Eris'
    
    body = Mond()   #Cometa / Asteroide
    body.nombre = '16R2'
    
    astros=[sol, mercurio, venus, tierra, marte, jupiter, saturno, urano, neptuno, pluton, body]
   
    condiciones(astros)	#Lectura de las condiciones iniciales
    iteraciones(astros)    #Se indican las iteraciones sobre los N cuerpos
    
    end = time.time()       #Final del tiempo de ejecución
    print('Tiempo de ejecución C/2016 R2 (Newton): '+ str(end - start)+' s')
    
if __name__ == '__main__':
    main()
    
    
    
