import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import time
import matplotlib.animation as animation

L1=1
L2=2
m1=5
m2=1

theta1_zero=np.pi/7 
theta1_zero_dot=0.5
theta2_zero=np.pi
theta2_zero_dot=1.7

'''
L1=float(input('Comprimendo fio primeiro fio: '))
L2=float(input('Comprimendo fio primeiro fio: '))
m1=float(input('Massa 1: '))
m2=float(input('Massa 2: '))
theta1_zero=float(input('Angulo Inicial 1: '))
theta1_zero_dot=float(input('Frequencia Inicial 1: '))
theta2_zero=float(input('Angulo 2 inicial: '))
theta2_zero_dot=float(input('Frequencia Inicial 2: '))
'''

g = 9.81 #g=pi²

def deriv(y, t, L1, L2, m1, m2):
    theta1, z1, theta2, z2 = y

    c= np.cos(theta1-theta2)
    s= np.sin(theta1-theta2) # é igual escrever c,s=np.cos(theta1-theta2), np.sin(theta1-theta2)

    theta1dot = z1
    
    z1dot = (m2*g*np.sin(theta2)*c - m2*s*(L1*z1**2*c + L2*z2**2) -
             (m1+m2)*g*np.sin(theta1)) / L1 / (m1 + m2*s**2)
    
    theta2dot = z2
    
    z2dot = ((m1+m2)*(L1*z1**2*s - g*np.sin(theta2) + g*np.sin(theta1)*c) + 
             m2*L2*z2**2*s*c) / L2 / (m1 + m2*s**2)
    
    return theta1dot, z1dot, theta2dot, z2dot


tempo_final, dt = 50, 0.01

t = np.arange(0, tempo_final+dt, dt) #valores de tempo

y0 = np.array([theta1_zero,theta1_zero_dot ,theta2_zero, theta2_zero_dot]) # lista das condições iniciais

y = odeint(deriv, y0, t, args=(L1, L2, m1, m2)) #scipy.integrate.odeint(func,y0,args=())
                                                #integra um sistema de EDO's
theta1 = y[:,0] #y=[theta1,frequencia1,theta2,frequencia2]
theta2 = y[:,2] #y[:,0]--->y=[theta1,frequencia1] e y[:,2]--->y=[theta2,frequencia2]

x1 = L1 * np.sin(theta1) #passando de coordenadas polares para cartesianas
y1 = -L1 * np.cos(theta1)
x2 = x1 + L2 * np.sin(theta2)
y2 = y1 - L2 * np.cos(theta2)

#Corte aqui####################################################################################################
#Apenas o pêndulo Duplo
'''
fig = plt.figure()
ax = fig.add_subplot(111, autoscale_on=False, xlim=(-1.2*(L1+L2), 1.2*(L1+L2)), ylim=(-1.2*(L1+L2), 1.2*(L1+L2)))
ax.grid()

line, = ax.plot([], [], 'o--',color='r', lw=3)
time_template = 'time = %.1fs'
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)


def init():
    line.set_data([], [])
    time_text.set_text('')
    return line, time_text


def animate(i):
    xx = [0, x1[i], x2[i]]
    yy = [0, y1[i], y2[i]]


    line.set_data(xx, yy)
    time_text.set_text(time_template % (i*dt))
    return line, time_text

ani = animation.FuncAnimation(fig, animate, np.arange(1, len(y)),
                              interval=25, blit=True, init_func=init,repeat=False)

plt.show()

'''

#############################################################################################################
#Pêndulo Duplo+ Desvio de Energia + Fantasma

theta1_dot=y[:,1]
theta2_dot=y[:,3]
v1=theta1_dot*L1
v2=theta2_dot*L2

fig = plt.figure()
ax = fig.add_subplot(111, autoscale_on=False, xlim=(-1.2*(L1+L2), 1.2*(L1+L2)), ylim=(-1.2*(L1+L2), 1.2*(L1+L2)))

line, = ax.plot([], [], '-o',color='b', lw=2)
line2, = ax.plot([], [], '.-',color='r', lw=0.000001)

time_template = 'Tempo = %.1fs '
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)
ener_template = 'Desvio Percentual da Energia = %.1f'
ener_text = ax.text(0.05, 0.8, '', transform=ax.transAxes)

def ener(theta1,theta1_dot,theta2,theta2_dot):
    V=-(m1+m2)*g*L1*np.cos(theta1)-m2*g*L2*np.cos(theta2)

    T=(1/2)*m1*L1**2*theta1_dot**2+(1/2)*m2*(L1**2*theta1**2+L2**2*theta2**2+2*L1*L2*theta1_dot*theta2_dot*np.cos(theta1-theta2))

    return T + V   

def init():
    line.set_data([], [])
    time_text.set_text('')
    ener_text.set_text('')
    line2.set_data([],[])
    return line, time_text,ener_text


def animate(i):
    xx = [0, x1[i], x2[i]]
    yy = [0, y1[i], y2[i]]
    xxx=[x1[0:i],x2[0:i]]
    yyy=[y1[0:i],y2[0:i]]

    line.set_data(xx, yy)
    time_text.set_text(time_template % (i*dt))
    line2.set_data(xxx,yyy)

    ener_text.set_text(ener_template % (((ener(x1[i],v1[i],x2[i],v2[i]-ener(y0[0],y0[1],y0[2],y0[3]))/(ener(y0[0],y0[1],y0[2],y0[3]))**2)**1/2)))
    return line, ener_text, time_text, line2

ani = animation.FuncAnimation(fig, animate, np.arange(1, len(y)),
                              interval=1, blit=True, init_func=init, repeat=False)                       
plt.show()

##################################################################################################################


