from math import *
import numpy as np
import scipy.optimize
import matplotlib.pyplot as plt
from tkinter import *
import os

os.system("clear")


def f(t, x, y):
    return y


def g(t, x, y):
    return ((y ** 2. * x) - ((k - 1.) * x * (Vmax ** 2. - x ** 2. - y ** 2.)) - (
            (k - 1.) / 2. * (Vmax ** 2. - x ** 2. - y ** 2.) * (y * 1. / tan(t)))) / (
                   (k - 1.) / 2. * (Vmax ** 2. - x ** 2. - y ** 2.) - y ** 2.)


def calculo(M1, thetacono):
    global Vmax
    global M3
    global T3
    global p3
    global rho3

    thetashock = 74 * pi / 180.  # Angulo de la onda de choque tentativa.

    h = -0.5 * pi / 180.  # Variable de iteracion del angulo theta

    # Es la lista que va a tener los distintos valores de theta que hacen Vtheta igual a 0 (theta surface)
    thetaiteracion = [thetashock]

    j = 0
    interval = 0.3 * pi
    while thetaiteracion[j] > thetacono:

        func = lambda delta: sin(thetashock) ** 6. + (-(M1 ** 2. + 2.) / (M1 ** 2.) - k * sin(delta) ** 2.) * sin(
            thetashock) ** 4. + ((2. * M1 ** 2. + 1.) / M1 ** 4. + sin(delta) ** 2. * (
                (k + 1.) ** 2. / 4. + (k - 1.) / M1 ** 2.)) * sin(thetashock) ** 2. + (-cos(delta) ** 2. / M1 ** 4.)

        # Angulo de la cuña 2D (direccion del flujo aguas abajo de la OCO)
        delta = scipy.optimize.fsolve(func, np.array(interval))
        M1n = M1 * sin(thetashock)
        M2n = (((k - 1.) * M1n ** 2. + 2.) / (2. * k * M1n ** 2. - (k - 1.))) ** 0.5
        M2 = M2n / (sin(thetashock - delta[0]))

        # Chequeo si no estoy calculando la onda fuerte
        if M2 > 0.999 * M1:
            break
        T0 = T1 + (M1 ** 2. * k * R * T1) / (2. * Cp)
        Vmax = (2. * Cp * T0) ** 0.5
        V = Vmax * ((2. / ((k - 1.) * M2 ** 2.)) + 1.) ** (-0.5)

        theta = []
        Vr = []
        Vtheta = []

        # Condiciones iniciales
        theta.append(thetashock)  # t
        Vr.append(V * cos(thetashock - delta[0]))  # x
        Vtheta.append((-V) * sin(thetashock - delta[0]))  # y

        i = 0

        # Runge-Kutta de 4to orden para una ecuacion diferencial de 2do orden
        while -Vtheta[i] >= 0.:
            k1 = h * f(theta[i], Vr[i], Vtheta[i])
            l1 = h * g(theta[i], Vr[i], Vtheta[i])
            k2 = h * (f(theta[i], Vr[i], Vtheta[i]) + 0.5 * l1)
            l2 = h * g(theta[i] + 0.5 * h, Vr[i] + 0.5 * k1, Vtheta[i] + 0.5 * l1)
            k3 = h * (f(theta[i], Vr[i], Vtheta[i]) + 0.5 * l2)
            l3 = h * g(theta[i] + 0.5 * h, Vr[i] + 0.5 * k2, Vtheta[i] + 0.5 * l2)
            k4 = h * (f(theta[i], Vr[i], Vtheta[i]) + l3)
            l4 = h * g(theta[i] + 0.5 * h, Vr[i] + k3, Vtheta[i] + l3)

            theta.append(theta[i] + h)
            Vr.append(Vr[i] + (k1 + 2. * k2 + 2. * k3 + k4) / 6.)
            Vtheta.append(Vtheta[i] + (l1 + 2. * l2 + 2. * l3 + l4) / 6.)

            i = i + 1

        thetaiteracion.append(theta[i])  # Guardo el valor de theta surface en la lista thetaiteracion
        thetashock = thetashock - 0.1 * pi / 180.  # Reduzco un poco el theta de la onda de choque inicial
        j = j + 1

    M3 = (2 / ((k - 1) * ((Vr[i] / Vmax) ** (-2) - 1))) ** 0.5
    T3 = T0 / (1 + (k - 1) / 2 * M3 ** 2)
    p3 = p1 * ((1 + ((k - 1) / 2 * M2 ** 2)) / (1 + ((k - 1) / 2 * M3 ** 2))) ** (k / (k - 1)) * (
            1 + (2 * k * (((M1 * sin(thetashock)) ** 2 - 1)) / (k + 1)))
    rho3 = p3 / R / T3
    return thetashock


def dibuja():
    caly = 200.  # Valores para afinar la posicion de los dibujos
    calx = 50.
    intervalo_mach = 2. * caly / 5.  # Separa las lineas que representan el flujo libre, en y. *2 para que el perfil
    # este centrado y dividido 5 porque es la cantidad de lineas que hay

    canvas.delete(ALL)  # Borra lo dibujado para que no se superpongan

    # Inicializo los vertices del cono en 0
    x = [0, 0, 0]  # lista de puntos x de los paneles, para generar el poligono
    y = [0, 0, 0]  # Idem anterior pero en y

    # Borde de ataque del cono
    x[0] = w_max / 2. - 100.
    y[0] = h_max / 2.

    x[1] = x[0] + largo_cono
    y[1] = y[0] + r_cono

    x[2] = x[1]
    y[2] = y[0] - r_cono

    points = []  # Lista a llenar con los puntos del cono

    for i in range(3):
        points.append(x[i])
        points.append(y[i])

    # Crea el cono
    canvas.create_polygon(points, outline="black", fill="cyan", width=1)  # Dibuja el cono como un poligono

    # Flujo libre, dibuja las lineas del flujo libre que representan la velocidad que tiene, como el mach es chico,
    # este se multiplica por 30, y se le suman 50 para separarlas del borde
    canvas.create_line(50, intervalo_mach * 1, M1 * 30 + 50, intervalo_mach * 1, width=2, arrow=LAST)
    canvas.create_line(50, intervalo_mach * 2, M1 * 30 + 50, intervalo_mach * 2, width=2, arrow=LAST)
    canvas.create_line(50, intervalo_mach * 3, M1 * 30 + 50, intervalo_mach * 3, width=2, arrow=LAST)
    canvas.create_line(50, intervalo_mach * 4, M1 * 30 + 50, intervalo_mach * 4, width=2, arrow=LAST)
    canvas.create_text(220, h_max / 2., font="Arial 10", text="Minf = " + str(M1))

    # Se dibujan ambas ondas de choque
    puntos_onda1 = [x[0], y[0], x[0] + 1000., y[0] + tan(thetashock) * 1000.]
    puntos_onda2 = [x[0], y[0], x[0] + 1000., y[0] - tan(thetashock) * 1000.]
    canvas.create_line(puntos_onda1, width=2, fill="gray")
    canvas.create_line(puntos_onda2, width=2, fill="gray")

    # Se crea un poligono que borra (sobreescribe) la onda de choque inferior en el lugar donde va el texto
    x_topleft = 175
    y_topleft = 450
    x_topright = w_max - x_topleft
    y_topright = y_topleft
    x_botleft = x_topleft
    y_botleft = y_topleft + 125
    x_botright = x_topright
    y_botright = y_botleft

    puntospoligono = [x_topleft, y_topleft, x_topright, y_topright, x_botright, y_botright, x_botleft, y_botleft]
    canvas.create_polygon(puntospoligono, outline="gray", fill="white", width=1)

    # Coordenadas de calibracion del texto que muestra el angulo obtenido
    xtexto = w_max / 2
    ytexto = y_topleft + (y_botleft - y_topleft) / 2

    # Se escribe el texto
    thetashock_text = "{:.1f}".format(thetashock * 180 / pi)
    M3_text = "{:.2f}".format(M3)
    T3_text = "{:.2f}".format(T3)
    p3_text = "{:.0f}".format(p3)
    rho3_text = "{:.3f}".format(rho3)
    canvas.create_text(xtexto, ytexto, font="Calibri 15",
                       text="El angulo de la onda de choque es: " + thetashock_text + " grados\n"
                            + "El Mach en la superficie del cono es: " + M3_text + "\n"
                            + "La temperatura en la superficie del cono es: " + T3_text + " K\n"
                            + "La presion en la superficie del cono es: " + p3_text + " Pa\n"
                            + "La densidad en la superficie del cono es: " + rho3_text + " kg/m³")

    canvas.update()  # Actualiza el dibujo

    return


def dibuja_3d():
    fig = plt.figure()

    ax = fig.add_subplot(111, projection='3d')
    limite_ejes = 10

    r_final = np.tan(thetacono) * limite_ejes
    r_shock = np.tan(thetashock) * limite_ejes
    # ~ print("r_final =",r_final,"thetacono = ",thetacono*180/pi,"L =",L)

    theta = np.linspace(0, 2 * np.pi, 100)
    r = np.linspace(0, r_final, 100)
    r_shock = np.linspace(0, r_shock, 100)

    T, R = np.meshgrid(theta, r)
    T, R_shock = np.meshgrid(theta, r_shock)

    ANA = np.zeros((100, 100))
    for i in range(100):
        for j in range(100):
            ANA[i][j] = i

    z = R * np.cos(T)
    y = R * np.sin(T)
    x = ANA / limite_ejes

    z_shock = R_shock * np.cos(T)
    y_shock = R_shock * np.sin(T)

    ax.plot_surface(x, y, z, rstride=100, cstride=2, color='cyan')
    ax.plot_surface(x, y_shock, z_shock, rstride=100, cstride=2, color='white', alpha=0.3)

    # Para que el dibujo sea cuadrado
    ax.set_ylim([-limite_ejes, limite_ejes])
    ax.set_zlim([-limite_ejes, limite_ejes])
    ax.set_xlim([0, limite_ejes])

    plt.show()


def boton_calcular():
    global M1
    global thetacono
    global r_cono
    global thetashock
    global T1
    global p1

    thetacono = cono.get() * pi / 180
    M1 = mach.get()
    p1 = presion.get()
    T1 = temp.get()

    # Geometria del cono
    r_cono = tan(thetacono) * largo_cono  # Radio del cono

    thetashock = calculo(M1, thetacono)
    dibuja()
    dibuja_3d()


# Parametros de la ventana principal
w_max = 800  # Ancho de ventana
h_max = 400  # Alto de ventana

# Crea la ventana del dibujo
root = Tk()
root.title("Solucion de Taylor-Maccoll para un cono en flujo supersonico")
canvas = Canvas(root, width=w_max, height=h_max + 200, bg='white')
canvas.grid(column=0, row=0, columnspan=2, rowspan=8)

# Parametros del fluido
k = 1.4
Cp = 1004.
R = (k - 1.0) * Cp / k
T1 = 200
p0 = 1000000
p1 = 80000
rho1 = p1 / (R * T1)

# Valor inicial del Mach del flujo
M1 = 1.1

# Geometria del cono inicial
thetacono = 1. * pi / 180.
largo_cono = 400.  # Largo del cono
r_cono = tan(thetacono) * largo_cono  # Radio del cono

# Onda de choque correspondiente a las condiciones iniciales
thetashock = calculo(M1, thetacono)

# Se crean los widgets
text_slider1 = Label(root, anchor=W, text="Semi-angulo      \n del cono [°]       ")
text_slider1.grid(row=0, column=4, sticky=E)
cono = Scale(root, orient="horizontal", from_=1, to=20, resolution=1)
cono.grid(row=0, column=5)

text_slider2 = Label(root, anchor=W, text="Mach de la \n corriente libre       ")
text_slider2.grid(row=1, column=4, sticky=E)
mach = Scale(root, variable=M1, orient="horizontal", from_=1.1, to=5, resolution=0.1)
mach.grid(row=1, column=5)

text_slider3 = Label(root, anchor=W, text="Temperatura de la \n corriente libre [K]   ")
text_slider3.grid(row=2, column=4, sticky=E)
temp = Scale(root, variable=T1, orient="horizontal", from_=200, to=500, resolution=10)
temp.grid(row=2, column=5)

text_slider4 = Label(root, text="Presion de la \n corriente libre [Pa]   ")
text_slider4.grid(row=3, column=4, sticky=E)
presion = Scale(root, variable=p1, orient="horizontal", from_=80000, to=500000, resolution=10000)
presion.grid(row=3, column=5)

calcular = Button(root, text="Calcular", command=boton_calcular)
calcular.grid(row=4, column=5)

dibuja()  # Dibuja el primero, antes de tocar nada
dibuja_3d()
root.mainloop()  # Mantiene el programa dentro del loop, esperando a hacer algo
