#MODELO wCDM PLANO
#CONSIDERANDO PRIOR EN "h"
#DATOS EMPLEADOS: PANTHEON
#ANÁLISIS MEDIANTE EMCEE
#FECHA ÚLTIMA MODIFICACIÓN: 17/04/2025

import numpy as np
import emcee
import matplotlib.pyplot as plt
from scipy import integrate

#calculo para distancia comovil
def distancia_comovil(p, c, z):
    exp_w = 3*(1+p[3]) #calculamos la exponencial del termino de omega_lambda
    comovil =[]
    for z in z[0]:
        intg, error = integrate.quad(lambda x: 1/np.sqrt(p[1]*(1+x)**3+(1-p[1])*(1+x)**exp_w), 0, z)
        comovil.append(intg * c/(p[0]*100))
    return comovil

#calculo para hallar la funcion D_A(z)
def Da(p, c, z):
    comovil = distancia_comovil(p, c, z)
    return comovil/(1+z[0])

#calculo para hallar el modulo de distancia del modelo 
def modulo_distanica(p, c, z):
    da = Da(p, c, z)
    return 5*np.log10((1+z[0])*(1+z[1])*da)+25

#calculo magnitud aparente del modelo
def mag_aparente(p, c, z):
    mu_modelo = modulo_distanica(p, c, z)
    return mu_modelo+p[2]

#calculo del logartimo de la likelihood de SN (-0.5 chi^2)
def likelihood_SN(p, c, z, m_b_corr, invcov):
    m_b_modelo = mag_aparente(p, c, z)
    delta = m_b_corr - m_b_modelo
    return -0.5 * np.dot(np.dot(delta,invcov), delta)

#calculo de la likelihood externa
#fijamos un valor aproximado de h
def likelihood_externa(p):
    return -0.5*((p[0]-0.7301)/0.0099)**2

#calculamos la likelihood total
def likelihood_total(p, c, z, m_b_corr, invcov):
    l_sn = likelihood_SN(p, c, z, m_b_corr, invcov)
    l_ext = likelihood_externa(p)
    print(2*(l_sn + l_ext))
    return l_sn + l_ext

def log_probability(p, c, z, m_b_corr, invcov):
    if ((0.2<p[0]<1.0) and (0.0<p[1]<1.0)):
        return likelihood_total(p, c, z, m_b_corr, invcov)
    return -np.inf 

#cargamos los datos
z_hd = np.loadtxt('Pantheon+SH0ES.dat', skiprows = 1, usecols = 2, unpack = False)
z_hel = np.loadtxt('Pantheon+SH0ES.dat', skiprows = 1, usecols = 6, unpack = False)
m_b_corr = np.loadtxt('Pantheon+SH0ES.dat', skiprows = 1, usecols = 8, unpack = False)
cov_datos = np.loadtxt('Pantheon+SH0ES_STAT+SYS.cov', skiprows = 1)

#realizamos un filtro para trabajar únicamente con redshift
#mayores a 0.01
mascara = (z_hd > 0.01) & (z_hel > 0.01)
z_hd = z_hd[mascara]
z_hel = z_hel[mascara]
m_b_corr = m_b_corr[mascara]
indices_z = np.where(~mascara)[0] #recogemos los indices elimindos
z = np.array([z_hd, z_hel]) #array con los dos tipos de redshift

#de igual manera, hay que tratar a la matriz de covarianza
#elimanando filas y columnas correspondientes a z<=0.01
n = len(cov_datos) #numero de datos cov
k = np.sqrt(n) #dimension de la matriz cuadrada
cov_matriz = cov_datos.reshape((int(k),int(k))) #transformamos en matriz

   #inicializamos bucle para eliminar las filas y las columnas de la matriz
for elemento in indices_z: 
    cov_matriz = np.delete(cov_matriz, elemento, axis=0)
    cov_matriz = np.delete(cov_matriz, elemento, axis=1)
cov = cov_matriz

#calculamos la matriz inversa
invcov = np.linalg.inv(cov)

#inicializamos las variables a estudiar en unos valores próximos
#a los esperados
h = 0.7 #cte hubble reducida
omega_m = 0.3 #parametro densidad materia no relativista
M = -19.1 #magnitud absoluta supernovas
w = -1.0 #parametro de ecuación de estado

#recogemos dichas variables en un array
p = np.array([h,omega_m, M, w])

c = 299792.458 #velocidad luz vacio (SI)

#configuramos el proceso emcee
ndim = 4
nwalkers = 170

#randomizacion de los valores iniciales
p0 = np.random.randn(nwalkers, ndim)

p0[:,0] = p0[:,0]*0.02+ 0.7
p0[:,1] = p0[:,1]*0.06 + 0.31
p0[:,2] = p0[:,2]*0.004 -19.1
p0[:,3] = p0[:,3]*0.003 - 1.0

sampler = emcee.EnsembleSampler(
    nwalkers, ndim, log_probability, args=(c, z, m_b_corr, invcov)
)
sampler.run_mcmc(p0, 12500, progress=True);

samples = sampler.get_chain(flat = True)



#calculo valores medios y errores
h_media = np.mean(samples[:,0])
h_er = np.std(samples[:,0])
print(h_media, "+-", h_er)

om_media = np.mean(samples[:,1])
om_er = np.std(samples[:,1])
print(om_media, "+-", om_er)

M_media = np.mean(samples[:,2])
M_er = np.std(samples[:,2])
print(M_media, "+-", M_er)

w_media = np.mean(samples[:,3])
w_er = np.std(samples[:,3])
print(w_media, "+-", w_er)

omega_lambda = 1 - samples[:,1]
ol_media = np.mean(omega_lambda)
ol_er = np.std(omega_lambda)
print(ol_media, "+-", ol_er)


variables = np.array([samples[:,0], samples[:,1], samples[:,3],
                      omega_lambda]).T
labels = ["h", "omega_m", "w", "omega_lambda"]

encabezado_s = "h\tomega_m\tw\tomega_lambda"
np.savetxt("wcdm_pr_param.txt", variables, fmt='%.16f', delimiter='\t',
           header=encabezado_s, comments='')

#para comprobar convergencia
print(
    "Mean acceptance fraction: {0:.3f}".format(
        np.mean(sampler.acceptance_fraction)
    )
)

print(
    "Mean autocorrelation time: {0:.3f} steps".format(
        np.mean(sampler.get_autocorr_time())
    )
)
