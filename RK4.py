# -*- coding: utf-8 -*-
"""
Created on July 19 2018

@author: eddybrandon
"""
#Trabajando con el origen en el centro de masas

import numpy as np

AU = 149597870700       # Definiendo una UA
ao=1.2*10**(-10)        # Constante a_o

class Mond():
    nombre='ninguno'    #será dado después
    masa=0.0
    x=y=z=0.0
    vx=vy=vz=0.0
    r=0.0
    
    def acc(self, cuerpok, kxi, kyi, kzi, kxk, kyk, kzk, h):
        Dxki = (self.x+h*kxi/2) - (cuerpok.x+h*kxk/2)
        Dyki = (self.y+h*kyi/2) - (cuerpok.y+h*kyk/2)
        Dzki = (self.z+h*kzi/2) - (cuerpok.z+h*kzk/2)
        rki = np.sqrt(Dxki**2+Dyki**2+Dzki**2)   
        Gmk = cuerpok.masa          #GM de cada cuerpo
        lmk = np.sqrt(Gmk/ao)       #Longitud de masa
        Xki = lmk/rki               #Cantidad adimensional
        #f = Xki**2		        	#Función Newton
        f = Xki*(1+Xki+Xki**2+Xki**3)/(1+Xki+Xki**2) 			#Función MOND 
        ax = -ao*f*Dxki/rki
        ay = -ao*f*Dyki/rki
        az = -ao*f*Dzki/rki
        return ax, ay, az
        
        
def RungeKutta(astros):
    counter=0
    hh=1
    h=24*3600/hh   #Tamaño de paso: 1 día.
    NI=5000*hh     #Número de iteraciones. --> NI=8000000*hh para el cometa C/2016 R2 
    arch={}
    for cuerpoi in astros:
        arch[cuerpoi.nombre]=open("RK{}newton.dat".format(cuerpoi.nombre), "w")  #Archivos de texto donde se guardan los datos de salida.
        
    #Cálculo del Baricentro:
    mt=0    
    cx=cy=cz=0
    cvx=cvy=cvz=0
    for cuerpoi in astros:  #Masa total
        mt+=cuerpoi.masa
    for cuerpoi in astros:  #Rcm, Vcm
        cx+=(cuerpoi.masa)*(cuerpoi.x)/mt
        cy+=(cuerpoi.masa)*(cuerpoi.y)/mt
        cz+=(cuerpoi.masa)*(cuerpoi.z)/mt
        cvx+=(cuerpoi.masa)*(cuerpoi.vx)/mt
        cvy+=(cuerpoi.masa)*(cuerpoi.vy)/mt
        cvz+=(cuerpoi.masa)*(cuerpoi.vz)/mt
    for cuerpoi in astros:
        cuerpoi.x=cuerpoi.x-cx
        cuerpoi.y=cuerpoi.y-cy
        cuerpoi.z=cuerpoi.z-cz
        cuerpoi.vx=cuerpoi.vx-cvx
        cuerpoi.vy=cuerpoi.vy-cvy
        cuerpoi.vz=cuerpoi.vz-cvz        
    #Iteración del método de Runge-Kutta    
    while counter<NI:
        counter+=1
        k1v={}
        k1p={}
        k2v={}
        k2p={}
        k3v={}
        k3p={}
        k4v={}
        k4p={}
                                        
        for cuerpoi in astros:      #Cálculo de las k1 (6Np veces).
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
            
        for cuerpoi in astros:      #Cálculo de las k2 (6Np veces).
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
            
        for cuerpoi in astros:      #Cálculo de las k3 (6Np veces).
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
            
        for cuerpoi in astros:      #Cálculo de las k4 (6Np veces).
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
                
        for cuerpoi in astros:     #Calculando las u_i
            k1vx, k1vy, k1vz = k1v[cuerpoi]
            k1x, k1y, k1z = k1p[cuerpoi]
            k2vx, k2vy, k2vz = k2v[cuerpoi]
            k2x, k2y, k2z = k2p[cuerpoi]
            k3vx, k3vy, k3vz = k3v[cuerpoi]
            k3x, k3y, k3z = k3p[cuerpoi]
            k4vx, k4vy, k4vz = k4v[cuerpoi]
            k4x, k4y, k4z = k4p[cuerpoi]
        
            escribir='{:>11.15f} {:>11.15f} {:>11.15f}'.format(cuerpoi.x/AU, cuerpoi.y/AU, cuerpoi.z/AU) #Datos x, y, z en UA.
            
            if cuerpoi.nombre is not 'Sol':
                cuerpoi.vx += (k1vx+k2vx*2+k3vx*2+k4vx)*h/6
                cuerpoi.vy += (k1vy+k2vy*2+k3vy*2+k4vy)*h/6
                cuerpoi.vz += (k1vz+k2vz*2+k3vz*2+k4vz)*h/6
                cuerpoi.x += (k1x+k2x*2+k3x*2+k4x)*h/6
                cuerpoi.y += (k1y+k2y*2+k3y*2+k4y)*h/6
                cuerpoi.z += (k1z+k2z*2+k3z*2+k4z)*h/6          
                arch[cuerpoi.nombre].write(escribir+"\n")
            else:
                cuerpoi.x=cuerpoi.y=cuerpoi.z=0
                for cuerpok in astros:
                    if cuerpok is not cuerpoi:
                        cuerpoi.x-=(cuerpok.masa)*(cuerpok.x)/cuerpoi.masa      #La posición del Sol se calcula usando el hecho que el baricentro se encuentra en el origen.
                        cuerpoi.y-=(cuerpok.masa)*(cuerpok.y)/cuerpoi.masa
                        cuerpoi.z-=(cuerpok.masa)*(cuerpok.z)/cuerpoi.masa
                        arch[cuerpoi.nombre].write(escribir+"\n")

        for cuerpoi in astros:
            arch[cuerpoi.nombre].closed

#Lee las condiciones iniciales de cada objeto (Sol, planetas, asteroides y cometas)          
def condiciones(astros):
    for cuerpoi in astros:
        file=open('datos{}.dat'.format(cuerpoi.nombre), 'r')
        cuerpoi.masa = float(file.readline())*1000**3   #Valor de GM
        cuerpoi.x = float(file.readline())*1000         #Posiciones
        cuerpoi.y = float(file.readline())*1000
        cuerpoi.z = float(file.readline())*1000
        cuerpoi.vx = float(file.readline())*1000        #Velocidades
        cuerpoi.vy = float(file.readline())*1000
        cuerpoi.vz = float(file.readline())*1000
        file.close
                        
def main():
    sol = Mond()      
    sol.nombre = 'Sol'

    mercurio = Mond()
    mercurio.nombre = 'Mercurio'

    venus = Mond()
    venus.nombre = 'Venus'

    btierra = Mond()        #Baricentro Tierra-Luna.
    btierra.nombre = 'Tierra'
    
    tierra = Mond()        #Tierra
    tierra.nombre = 'TierraX'
    
    luna = Mond()
    luna.nombre = 'Luna'

    #Se utiliza la masa total del sistema planeta-lunas y la posición de su centro de masa.
    
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
    
    astro = Mond()
    astro.nombre = 'Astro'
    
    astros=[sol, mercurio, venus, btierra, marte, jupiter, saturno, urano, neptuno, pluton, eris, astro]
    condiciones(astros)	#condiciones iniciales
    RungeKutta(astros) 
    
if __name__ == '__main__':
    main()
