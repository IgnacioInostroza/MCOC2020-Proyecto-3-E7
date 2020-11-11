import numpy as np
import matplotlib.pylab as plt
from matplotlib import cm
from calor_de_hidratacion import Calor_de_hidratacion as Q_h

def Truncate(n, decimal=0):
    mult= 10**decimal
    return int (n * mult)/mult

def coords(i,j):
    return (dx*i,dy*j)

def imshowbien(u):
    plt.imshow(u.T[N_x::-1,:],cmap=plt.cm.coolwarm,interpolation='bilinear')
    cbar=plt.colorbar(extend='both', cmap=cm.coolwarm)
    ticks=np.arange(0,35,5)
    ticks_Text=["{}°".format(deg) for deg in ticks]
    cbar.set_ticks(ticks)
    #cbar.set_tickslabels(ticks_Text)
    plt.clim(0,30)
    plt.xlabel('b')
    plt.ylabel('a')
    xTicks_N=np.arange(0,N_x+1,3)
    yTicks_N=np.arange(0,N_y+1,3) 
    xTicks=[coords(i,0)[0] for i in xTicks_N]
    yTicks=[coords(0,i)[1] for i in yTicks_N]
    xTicks_Text=["{0:.2f}".format(tick) for tick in xTicks]
    yTicks_Text=["{0:.2f}".format(tick) for tick in yTicks]
    plt.xticks(xTicks_N,xTicks_Text,rotation='vertical')
    plt.yticks(yTicks_N,yTicks_Text)
    plt.margins(0.2)
    plt.subplots_adjust(bottom=0.15)
    plt.show()
    
largo= 1
N_x=  40

ancho= 1
N_y=  40

alto= 1
N_z = 40

dx =  largo/N_x
dy=  ancho/N_y
dz =  alto/ N_z


seg    = 1
minuto = 60 * seg
hora   = 60 * minuto
dia    = 24 * hora


dt     = 5
Days = 500*dt
paso=10
frame_num=1

puntos= []
for i in range(14):
    puntos.append([])

K=0.001495/1000
c= 0.023/1000
rho= 2676

alpha= K*dt / (c*rho*dx**2)

u_k= np.zeros( (N_x+1,N_y+1,N_z+1),  dtype= np.double)
u_k_1= np.zeros( (N_x+1,N_y+1,N_z+1),  dtype= np.double)

u_k[:,:,:]= 21 # condicion inicial

for t in range (np.int32(Days/dt)):
    Q_h_t= Q_h(t*100)

    #u_k[0,:,:] =  20              # izquerda fija
    u_k[0,:,:] =      u_k[1,:,:]   # izquerda derivada

    #u_k[-1,:,:] =  20                # derecha fija
    u_k[-1,:,:] =     u_k[-2,:,:]     # derecha derivada

    #u_k[:,0,:]=  20                        # atras fija
    u_k[:,0,:]=       u_k[:,1,:]            # atras derivada

    #u_k[:,-1,:] =  20                    # frente fija
    u_k[:,-1,:] =    u_k[:,-2,:]          # frente derivada

    #u_k[:,:,0]=  20                      # abajo fijo
    u_k[:,:,0] =          u_k[:,:,1]      # abajo derivada

    temp_t= 21 + 10*np.sin(t*2*np.pi*dt/Days)
    u_k[:,:,-1]=  temp_t                    # arriba fijo
    #u_k[:,:,-1] =  1*dz + u_k[:,:,-2]      # arriba derivada
    
    
    for i in range(1,N_x):
        for j in range(1,N_y):
            for k in range (1,N_z):
                nabla_u_k = (u_k[i-1,j,k] + u_k[i+1,j,k] + u_k[i,j-1,k]+ u_k[ i,j+1,k] +u_k[ i,j,k-1]+ u_k[ i,j,k+1] -6*u_k[i,j,k])
                u_k_1[i,j,k]= u_k[i,j,k]+ alpha* nabla_u_k + Q_h_t

    u_k=u_k_1
    
    L_0= 1
    L_1= int(N_x/4)
    L_2= int(N_x/2)
    L_3= int(3*N_x/4)
    L_4= int(N_x)
    

    if t % paso ==0:
        puntos[0].append(u_k[L_3,L_2,L_2])
        puntos[1].append(u_k[L_2,L_2,L_2])
        puntos[2].append(u_k[L_2,L_0,L_2])
        
        puntos[3].append(u_k[L_2,L_2,L_4])
        puntos[4].append(u_k[L_0,L_2,L_2])
        puntos[5].append(u_k[L_2,L_4,L_2])
        
        puntos[6].append(u_k[L_2,L_2,L_3])
        puntos[7].append(u_k[L_2,L_2,L_1])
        puntos[8].append(u_k[L_2,L_3,L_2])
        
        puntos[9].append(u_k[L_4,L_2,L_2])
        puntos[10].append(u_k[L_2,L_2,L_0])
        puntos[11].append(u_k[L_1,L_2,L_2])
        
        puntos[12].append(temp_t)       # ambiente = cond borde
        puntos[13].append(u_k[L_2,L_1,L_2])
        
        tiem=(t+1)*dt
        
        dias    = Truncate(tiem/dia)
        horas   = Truncate((tiem-dia*dias)/hora)
        minutos = Truncate((tiem-dia*dias- horas*hora)/minuto)
        segundos = Truncate((tiem-dia*dias- horas*hora-minutos*minuto))
        deci_seg = Truncate((tiem-dia*dias- horas*hora-minutos*minuto-segundos)*10)
        plt.figure()
        
        
        titulo_x=  f" Corte en x  tiempo =  {int(minutos)} minutos, {int(segundos)} segundos, deci_seg= {deci_seg} "
        titulo_y=  f" Corte en y  tiempo =  {int(minutos)} minutos, {int(segundos)} segundos, deci_seg= {deci_seg} "
        titulo_z=  f" Corte en z  tiempo =  {int(minutos)} minutos, {int(segundos)} segundos, deci_seg= {deci_seg} "
        
        plt.figure()
        plt.title(titulo_x)
        imshowbien(u_k[:,:,int(N_z/2)])
        plt.savefig("Ejemplo/frame_corte_x{0:04.0f}.png".format(frame_num))
        
        plt.figure()
        plt.title(titulo_y)
        imshowbien(u_k[int(N_z/2),:,:])
        plt.savefig("Ejemplo/frame_corte_y{0:04.0f}.png".format(frame_num))
        
        plt.figure()
        plt.title(titulo_z)
        imshowbien(u_k[:,int(N_z/2),:])
        plt.savefig("Ejemplo/frame_corte_z{0:04.0f}.png".format(frame_num))
        
        frame_num+=1


def graficar_avance_puntos():
    a=1
    plt.figure()
    for i in puntos:
        plt.plot(i,label= f"Sensor {a}")
        a+=1
    
    plt.legend(loc=1)
    plt.title("Evolucion de temperatura en puntos")
    plt.ylabel("Temeratura, $T$  [°C]")
    plt.xlabel("Tiempo $t$ [segundos]")
    plt.show()
