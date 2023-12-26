import os
import numpy as np
import math
import pyshtools
from scipy import interpolate

FACT = 180./np.pi # Если делим на это, то в радианы, если умножаем, то в градусы

class igrf: # A simple class to put the igrf file values into
  def __init__(self, time, coeffs, parameters):
     self.time = time
     self.coeffs = coeffs
     self.parameters = parameters
     
           
def myReadCOEF(date):
    g, h = np.zeros((14,14), float), np.zeros((14,14), float)
    flag = True
    with open("IGRF13.txt", "r") as file:
    
        for line in file:

            s = ' '.join(line.split()).split() # беру значения строки c нормальным пробелами        

            # Беру заголовки
            if flag == True:
                if line[0] == '#':
                    continue
   
                if line[:2] == "cs":
                    head_name = s[:]
                    continue
        
                if line[:2] == "gh":
                    head_type = s[:]
                    head_years =  [float(c) for c in s[3:]]
                
                    flag = False
                    continue
        
            k = s[:1][0]
            n, m = int(s[1]), int(s[2])

            coeffs = [float(c) for c in s[3:]]
            coeffs[-1] = coeffs[-2] + 5 * coeffs[-1] 
            f = interpolate.interp1d(head_years, coeffs, fill_value='extrapolate')
            coef = f(date)
        
            if k == 'g':
                g[n, m] = coef
            elif k == 'h':
                h[n, m] = coef
    return g, h



def load_coeffs(filename):
    """
    load igrf12 coeffs from file
    :param filename: file which save coeffs (str)
    :return: g and h list one by one (list(float))
    """
    gh = []
    gh2arr = []
    buf = []
    with open(filename) as f:
        text = f.readlines()
        for a in text:
            if a[:2] == 'g ' or a[:2] == 'h ':
                b = a.split()[3:]
                b = [float(x) for x in b]
                gh2arr.append(b)
        gh2arr = np.array(gh2arr).transpose()
        N = len(gh2arr)
        for i in range(N):
            if i < 19:
                for j in range(120):
                    gh.append(gh2arr[i][j])
            else:
                ident = -1
                for p in gh2arr[i]:
                    ident+=1
                    if i == 25:
                        gh.append(buf[ident] + p*5)
                        
                    else:
                        if i == 24:
                            buf.append(p)
                        gh.append(p)
        gh.append(0)
        return gh


def get_coeffs(filename, date):
    """
    :param filename: str in dir COEF
    :param gh: list from load_coeffs
    :param date: float
    :return: list: g, list: h
    """
    gh = load_coeffs(os.path.dirname(os.path.abspath(__file__)) + '\\' + filename)
    
    if date < 1900.0 or date > 2030.0:
        print('This subroutine will not work with a date of ' + str(date))
        print('Date must be in the range 1900.0 <= date <= 2030.0')
        print('On return [], []')
        return [], []
    elif date >= 2020.0:
        if date > 2025.0:
            # not adapt for the model but can calculate
            print('This version of the IGRF is intended for use up to 2025.0.')
            print('values for ' + str(date) + ' will be computed but may be of reduced accuracy')
        t = date - 2020.0
        tc = 1.0
        #     pointer for last coefficient in pen-ultimate set of MF coefficients...
        ll = 3060+195
        nmx = 13
        nc = nmx * (nmx + 2)
    else:
        t = 0.2 * (date - 1900.0)
        ll = int(t)
        t = t - ll
        #     SH models before 1995.0 are only to degree 10
        if date < 1995.0:
            nmx = 10
            nc = nmx * (nmx + 2)
            ll = nc * ll
        else:
            nmx = 13
            nc = nmx * (nmx + 2)
            ll = int(0.2 * (date - 1995.0))
            #     19 is the number of SH models that extend to degree 10
            ll = 120 * 19 + nc * ll
        tc = 1.0 - t

    g, h = [], []
    temp = ll-1
    default = 0
    for n in range(nmx+1):
        buf = []
        g.append([])
        h.append([])
        if n == 0:
            g[0].append(default)
        for m in range(n+1):
            if m != 0:
                g[n].append(tc*gh[temp] + t*gh[temp+nc])
                h[n].append(tc*gh[temp+1] + t*gh[temp+nc+1])
                temp += 2
                # print(n, m, g[n][m], h[n][m])
            else:
                g[n].append(tc*gh[temp] + t*gh[temp+nc])
                h[n].append(default)
                temp += 1
                # print(n, m, g[n][m], h[n][m])
    

    datamy = []
    cell = 0
    """
    for column in range(26):
        # Перебираем весь одномерный массив...
        if column <= 19:
            cell = 120*(column+1)
            val = gh[120*column:120*(column+1)] + [0] * (195-120)
            if column == 19:
                cell += (195-120)
        else:
            val = gh[cell:cell+195]
            cell += 195
        datamy.append(val)
    """
    for column in range(26):
        # Перебираем весь одномерный массив...
        if column <= 19:
            cell = 120*(column+1)
            val = gh[120*column:120*(column+1)] + [0] * (195-120)
            if column == 19:
                cell += (195-120)
        else:
            val = gh[cell:cell+195]
            cell += 195
        datamy.append(val)
    datamy = np.array(datamy).transpose()
    return datamy, gh, g, h


def gg_to_geo(h, gdcolat):
    """
    Compute geocentric colatitude and radius from geodetic colatitude and
    height.

    Parameters
    ----------
    h : ndarray, shape (...)
        Altitude in kilometers.
    gdcolat : ndarray, shape (...)
        Geodetic colatitude

    Returns
    -------
    radius : ndarray, shape (...)
        Geocentric radius in kilometers.
    theta : ndarray, shape (...)
        Geocentric colatitude in degrees.
    
    sd : ndarray shape (...) 
        rotate B_X to gd_lat 
    cd :  ndarray shape (...) 
        rotate B_Z to gd_lat 

    References
    ----------
    Equations (51)-(53) from "The main field" (chapter 4) by Langel, R. A. in:
    "Geomagnetism", Volume 1, Jacobs, J. A., Academic Press, 1987.
    
    Malin, S.R.C. and Barraclough, D.R., 1981. An algorithm for synthesizing 
    the geomagnetic field. Computers & Geosciences, 7(4), pp.401-405.

    """
    # Use WGS-84 ellipsoid parameters

    myBuf = gdcolat * FACT
    eqrad = 6378.137 # equatorial radius
    flat  = 1/298.257223563
    plrad = eqrad*(1-flat) # polar radius
    ctgd  = np.cos(gdcolat)
    stgd  = np.sin(gdcolat)
    a2    = eqrad*eqrad
    a4    = a2*a2
    b2    = plrad*plrad
    b4    = b2*b2
    c2    = ctgd*ctgd
    s2    = 1-c2
    rho   = np.sqrt(a2*s2 + b2*c2)
    
    rad   = np.sqrt(h*(h+2*rho) + (a4*s2+b4*c2)/rho**2)

    cd    = (h+rho)/rad
    sd    = (a2-b2)*ctgd*stgd/(rho*rad)
    
    cthc  = ctgd*cd - stgd*sd           # Also: sthc = stgd*cd + ctgd*sd
    thc   = np.rad2deg(np.arccos(cthc)) # arccos returns values in [0, pi]
    
    return rad, thc/FACT, sd, cd


def geo_to_gg(radius, theta):
    """
    Compute geodetic colatitude and vertical height above the ellipsoid from
    geocentric radius and colatitude.

    Parameters
    ----------
    radius : ndarray, shape (...)
        Geocentric radius in kilometers.
    theta : ndarray, shape (...)
        Geocentric colatitude in degrees.

    Returns
    -------
    height : ndarray, shape (...)
        Altitude in kilometers.
    beta : ndarray, shape (...)
        Geodetic colatitude

    Notes
    -----
    Round-off errors might lead to a failure of the algorithm especially but
    not exclusively for points close to the geographic poles. Corresponding
    geodetic coordinates are returned as NaN.

    References
    ----------
    Function uses Heikkinen's algorithm taken from:

    Zhu, J., "Conversion of Earth-centered Earth-fixed coordinates to geodetic
    coordinates", IEEE Transactions on Aerospace and Electronic Systems}, 1994,
    vol. 30, num. 3, pp. 957-961

    """
    
    # Use WGS-84 ellipsoid parameters
    a =  6378.137  # equatorial radius
    b =  6356.752  # polar radius
    
    a2 = a**2
    b2 = b**2

    e2 = (a2 - b2) / a2  # squared eccentricity
    e4 = e2*e2
    ep2 = (a2 - b2) / b2  # squared primed eccentricity

    r = radius * np.sin(theta)
    z = radius * np.cos(theta)

    r2 = r**2
    z2 = z**2

    F = 54*b2*z2

    G = r2 + (1. - e2)*z2 - e2*(a2 - b2)

    c = e4*F*r2 / G**3

    s = (1. + c + np.sqrt(c**2 + 2*c))**(1./3)

    P = F / (3*(s + 1./s + 1.)**2 * G**2)

    Q = np.sqrt(1. + 2*e4*P)

    r0 = -P*e2*r / (1. + Q) + np.sqrt(
        0.5*a2*(1. + 1./Q) - P*(1. - e2)*z2 / (Q*(1. + Q)) - 0.5*P*r2)

    U = np.sqrt((r - e2*r0)**2 + z2)

    V = np.sqrt((r - e2*r0)**2 + (1. - e2)*z2)

    z0 = b2*z/(a*V)

    height = U*(1. - b2 / (a*V))

    beta = 90 - np.arctan2(z + ep2*z0, r)

    return height, beta/FACT


def load_coef(filename, leap_year=None):
    """
    Load shc-file and return coefficient arrays.

    Parameters
    ----------
    filepath : str
        File path to spherical harmonic coefficient shc-file.
    leap_year : {True, False}, optional
        Take leap year in time conversion into account (default). Otherwise,
        use conversion factor of 365.25 days per year.

    Returns
    -------
    time : ndarray, shape (N,)
        Array containing `N` times for each model snapshot in modified
        Julian dates with origin January 1, 2000 0:00 UTC.
    coeffs : ndarray, shape (nmax(nmax+2), N)
        Coefficients of model snapshots. Each column is a snapshot up to
        spherical degree and order `nmax`.
    parameters : dict, {'SHC', 'nmin', 'nmax', 'N', 'order', 'step'}
        Dictionary containing parameters of the model snapshots and the
        following keys: ``'SHC'`` shc-file name, `nmin` minimum degree,
        ``'nmax'`` maximum degree, ``'N'`` number of snapshot models,
        ``'order'`` piecewise polynomial order and ``'step'`` number of
        snapshots until next break point. Extract break points of the
        piecewise polynomial with ``breaks = time[::step]``.

    """
    leap_year = True if leap_year is None else leap_year

    with open(filename, 'r') as f:

        data = np.array([])
        words = f.readlines()
        for i, line in enumerate(words):

            if line[0] == '#':
                continue

            read_line = np.fromstring(line, sep=' ') # беру значения строки (там пробелы в разнобой)
            size = len(read_line)
            if size == 7:
                name = os.path.split(filename)[1]  # file name string
                val_read_line = read_line.astype(np.int64).tolist()
                values = [name] + val_read_line

            else:
                data = np.append(data, read_line) # добавляем каждую строку справа

        # unpack parameter line
        keys = ['SHC', 'nmin', 'nmax', 'N', 'order', 'step', 'start_year', 'end_year']
        parameters = dict(zip(keys, values))
        
        time = data[:parameters['N']]
        coeffs = data[parameters['N']:].reshape((-1, parameters['N']+2))
        coeffs = np.squeeze(coeffs[:, 2:])  # discard columns with n and m


    return igrf(time, coeffs, parameters)

def legendre_poly(nmax, theta):
    """
    Returns associated Legendre polynomials `P(n,m)` (Schmidt quasi-normalized)
    and the derivative :math:`dP(n,m)/d\\theta` evaluated at :math:`\\theta`.

    Parameters
    ----------
    nmax : int, positive
        Maximum degree of the spherical expansion.
    theta : ndarray, shape (...)
        Colatitude in radians :math:`[0^\\circ, 180^\\circ]`
        of arbitrary shape.

    Returns
    -------
    Pnm : ndarray, shape (n, m, ...)
          Evaluated values and derivatives, grid shape is appended as trailing
          dimensions. `P(n,m)` := ``Pnm[n, m, ...]`` and `dP(n,m)` :=
          ``Pnm[m, n+1, ...]``

    """

    #costh = np.cos(theta)
    costh = theta
    sinth = np.sqrt(1-costh**2)

    Pnm = np.zeros((nmax+1, nmax+2) + costh.shape)
    Pnm[0, 0] = 1  # is copied into trailing dimenions
    Pnm[1, 1] = sinth  # write theta into trailing dimenions via broadcasting

    rootn = np.sqrt(np.arange(2 * nmax**2 + 1))

    # Recursion relations after Langel "The Main Field" (1987),
    # eq. (27) and Table 2 (p. 256)
    for m in range(nmax):
        Pnm_tmp = rootn[m+m+1] * Pnm[m, m]
        Pnm[m+1, m] = costh * Pnm_tmp

        if m > 0:
            Pnm[m+1, m+1] = sinth*Pnm_tmp / rootn[m+m+2]

        for n in np.arange(m+2, nmax+1):
            d = n * n - m * m
            e = n + n - 1
            Pnm[n, m] = ((e * costh * Pnm[n-1, m] - rootn[d-e] * Pnm[n-2, m])
                         / rootn[d])

    # dP(n,m) = Pnm(m,n+1) is the derivative of P(n,m) vrt. theta
    Pnm[0, 2] = -Pnm[1, 1]
    Pnm[1, 2] = Pnm[1, 0]
    for n in range(2, nmax+1):
        Pnm[0, n+1] = -np.sqrt((n*n + n) / 2) * Pnm[n, 1]
        Pnm[1, n+1] = ((np.sqrt(2 * (n*n + n)) * Pnm[n, 0]
                       - np.sqrt((n*n + n - 2)) * Pnm[n, 2]) / 2)

        for m in np.arange(2, n):
            Pnm[m, n+1] = (0.5*(np.sqrt((n + m) * (n - m + 1)) * Pnm[n, m-1]
                           - np.sqrt((n + m + 1) * (n - m)) * Pnm[n, m+1]))

        Pnm[n, n+1] = np.sqrt(2 * n) * Pnm[n, n-1] / 2

    return Pnm

N = 13 # Количество суммы сумм...

Rz = 6371.2

alt = 100 # radius in km (6370.20 to 6971.20)
lat = 25 # Широта Enter the decimal latitude (-90 to 90) (- for Southern hemisphere)
lon = 50 # Долгота Enter th decimal longitude (-180 to 180) (- for Western hemisphere)
date = 2023

print("Geodetic:", alt, lat, lon)

colat = (90. - lat) / FACT # Нахожу широту в сферических координатах в радианах
lon = lon / FACT # Перевожу долготу в радианы

r, theta, sd, cd  = gg_to_geo(alt, colat)
lamda = lon
height, beta = geo_to_gg(r, theta)
#r, theta, lamda = geodetic_to_geocentricWGS84(alt, colat, lon) # Преобразование эллипсоидных геодезических координат (h, φ, λ) в сферические геоцентрические координаты (r, φ', λ)
#r, theta, lamda = Rz + alt, colat, lon
print("Sphere:", r, theta * FACT, lamda * FACT)
print("GeodeticReturn:", height, beta * FACT, lamda * FACT)

# Load in the file of coefficients
IGRF_FILE = 'IGRF13.shc'
igrf = load_coef(IGRF_FILE, None)
datamy, gh, g, h = get_coeffs('IGRF13.txt', 2023)
# Interpolate the geomagnetic coefficients to the desired date(s)
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
f = interpolate.interp1d(igrf.time, igrf.coeffs, fill_value='extrapolate')
coeffs = f(date)
print(igrf.time, len(igrf.time))
print(igrf.coeffs)
print("_________")
my =np.array(datamy).reshape((-1, 26))
print(my)
fmy = interpolate.interp1d(igrf.time, my, fill_value='extrapolate')
coeffsmy = fmy(date)
print(len(coeffsmy))

g, h = myReadCOEF(date)
# Можно увидеть как коэффициенты меняются
for i in range(195):
    print(f"{i} {coeffsmy[i]} =? {coeffs[i]} ===? {my[i][-1]}")
#fnew = interpolate.interp1d(gh.time, gh, fill_value='extrapolate')




# Получив коэффициенты и координаты сферические, мне нужно уметь вычислять полиномы ассоциированные Лежандра
# Вычислять аналитически черех Maxima и реккурентные формулы - не дало успеха
# поэтому попробуем через Scipy
from scipy.special import lpmn, lpmv # Присоед. Полиномы и их производные в одной функции clpmn

#Пример получения значений
n_test = 2 # порядок полиномов
m_test = 2  # степень производных
x = np.cos(theta)
#x = 0.4203

# Необходима нормировка
def norm(m, n):
    
    #kron = 1 if m==0 else 0
    if m == 0:
        semiNormalSmidth = 1
    if m > 0:
        semiNormalSmidth = (-1)**m * math.sqrt(2*math.factorial(n-m)/math.factorial(n+m))
    return semiNormalSmidth

def PnmNorm(m, n, x):
    '''
    val = lpmn(m, n, x)
    Pnm = values[0][-1][-1] * norm(m, n)
    dPnm = values[1][-1][-1] * norm(m, n)
    '''
    #print(x*180/np.pi)
    val = legendre_poly(n, x)
    if m == 0:
        Pnm = val[n, 0]
        dPnm = val[0, n+1]
    else:
        Pnm = val[n, m]
        dPnm = val[m, n+1]

    return Pnm, dPnm

values = lpmn(m_test, n_test, x)
values2 = pyshtools.legendre.legendre(n_test, x, 'schmidt', 1, 0, False)
val, dval = pyshtools.legendre.PlmSchmidt_d1(n_test, x, 1, 0)
pdp = legendre_poly(n_test, x) 
# Вычисляются члены к которым нужно обращаться так Pnm0: pdp[n, 0] pdp[0, n+1],
# и если m > 0 Pnm:pdp[n, m], pdp[m, n+1] . У них в коде сначало член суммы с этим, а потом другой член суммы с другим..
#  не знаю насколько это верно
# но там есть еще и косинус внутри, поэтому надо угол в радианах подставлять!
print(f"Pmn:{values[0][-1][-1] * norm(m_test, n_test)} dP:{values[1][-1][-1] * norm(m_test, n_test)}")
print(f"Pmn:{values2[n_test][m_test]}")
#print(f"Pmn:{val[n_test**2+m_test-1]} dPmn:{dval[n_test**2+m_test-1]}")
print(f"n:{n_test} m:{m_test}")
print(f"Pmn:{pdp[n_test,m_test]} dP:{pdp[m_test, 1+n_test]}")
print(f"Pmn:{PnmNorm(m_test, n_test, x)[0]} dP:{PnmNorm(m_test, n_test, x)[1]}")

resultUr, resultUt, resultUl = 0, 0, 0
# Получим Ur
"""
num = 0
nmin = 1
radius = r / 6371.2
r_n = radius**(-(nmin+2))
for n in range(1, N + 1):
    resultUr += (n+1) * PnmNorm(0, n, x)[0] * r_n * coeffs[..., num]
    num += 1
    for m in range(1, n + 1):
        #s = g[n][m] * math.cos(m*lamda) + h[n][m] * math.sin(m*lamda)
        Pnm = PnmNorm(m, n, x)[0]
        resultUr += ((n+1) * Pnm * r_n
                             * (coeffs[..., num] * math.cos(m*lamda)
                                + coeffs[..., num+1] * math.sin(m*lamda)))
        num += 2
    r_n = r_n / radius # equivalent to r_n = radius**(-(n+2))
"""
# Получим Ur
radius = 6371.2 / r
resultUr = 0
for n in range(1, N + 1):
    result = 0
    for m in range(0, n + 1):
        Pnm = PnmNorm(m, n, x)[0]
        result += (Pnm * (g[n][m] * math.cos(m*lamda) + h[n][m] * math.sin(m*lamda)))
    resultUr += radius**(n+2) * (n+1) * result


# Получим Ut/Uphi
for n in range(1, N + 1):
    koef_start = (Rz/r)**(n+2)
    result = 0
    for m in range(0, n + 1):
        s = g[n][m] * math.cos(m*lamda) + h[n][m] * math.sin(m*lamda)
        
        # handle poles using L'Hopital's rule
        dPnm = PnmNorm(m, n, x)[1]
        
        result+= s * dPnm
        #result+= s * div_Pnm
    resultUt += koef_start * result
resultUt = -resultUt

# Получим Ul
koef_start1 = -1/(math.sin(theta))
for n in range(1, N + 1):
    koef_start2 = (Rz/r)**(n+2)
    result = 0
    for m in range(0, n + 1):
        s = m*(-g[n][m] * math.sin(m*lamda) + h[n][m] * math.cos(m*lamda))
        Pnm = PnmNorm(m, n, x)[0]
        result+= s * Pnm
    resultUl += koef_start2 * result
resultUl = koef_start1 * resultUl

#Pnm, dPnm = PnmdPnm[0], PnmdPnm[1]
#div_Pnm = Pnm if theta == 0. else dPnm / math.sin(theta)
#div_Pnm = -Pnm if theta == np.pi else div_Pnm
        

"""
Br = resultUr
Bp = resultUl
Bt = resultUt
"""     
Bx = -resultUt;
By = resultUl;
Bz = -resultUr
        
mag = (Bx*Bx + By*By + Bz*Bz)**0.5
#print(Br, Bp, Bt)
print(Bx, By, Bz, mag)
#print(mag)

def xyz2dhif(x, y, z):
    """Calculate D, H, I and F from (X, Y, Z)
      
    Based on code from D. Kerridge, 2019
    
    Parameters
    ---------------
    X: north component (nT) : float
    Y: east component (nT) : float
    Z: vertical component (nT) : float
    
    Returns
    ------
    A tuple: (D, H, I, F) : float
    D: declination (degrees) : float
    H: horizontal intensity (nT) : float
    I: inclination (degrees) : float
    F: total intensity (nT) : float
    


    """
    hsq = x*x + y*y
    hoz  = np.sqrt(hsq)
    eff = np.sqrt(hsq + z*z)
    dec = np.arctan2(y,x)
    inc = np.arctan2(z,hoz)
    
    return dec*FACT, inc*FACT, hoz, eff


dec, inc, hoz, eff = xyz2dhif(Bx,By,Bz)
print(xyz2dhif(Bx,By,Bz))
#Formuls 32077.226410163617 1543.5968084227234 26270.388645348514 41490.534641693004
#2.755020311988419, 39.284079590977925, 32114.344852041915, 41490.534641693004
# Наверное тут не хватает интерполяции для коэффициентов

#IGRF script 32063.576130545727, 1558.0452976369447, 26498.338270523123, 41625.2609664835
#2.781949752267602, 39.53821418740564, 32101.40837142181

#Github 32143.522520374827 1543.5968084227231 26189.22934583435
#2.7493468019415515, 39.13943103612582, 32180.56449356997, 41490.53464169301 <--- Ближе всех к коду от производителя

print(0)



