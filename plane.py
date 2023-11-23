from math import sqrt, sin, tan, cos, pi
import numpy as np
import matplotlib.pyplot as plt

# consts
dt = 0.01
lat_0 = 55.508370 * pi / 180
lon_0 = 37.518120 * pi / 180
h_0 = 170
V1_0 = 0
V2_0 = 0
V3_0 = 0

u = 366.25 / (365.25 * 24 * 3600)
a = 6378137
e_2 = 6.6943799901413 * 10 ** (-3)


def RE(lat):
    return a / sqrt(1 - e_2 * sin(lat) ** 2)


def RN(lat):
    return a * (1 - e_2) / (sqrt(1 - e_2 * sin(lat) ** 2)) ** 3


def Omega(V1, V2, lat, h):
    RE_ = RE(lat)
    RN_ = RN(lat)
    first = -V2 / (RN_ + h)
    second = V1 / (RE_ + h)
    third = V1 * tan(lat) / (RE_ + h)
    return np.array([first, second, third])


def ux(lat):
    return np.array([0, u * cos(lat), u * sin(lat)])


def roof(vec):
    answer = [[0] * 3 for _ in range(3)]
    answer[0][1] = vec[2]
    answer[0][2] = -vec[1]
    answer[1][0] = -vec[2]
    answer[1][2] = vec[0]
    answer[2][0] = vec[1]
    answer[2][1] = -vec[0]
    return np.array(answer)


def L(corse, tang, kren):
    l1 = [cos(tang) * sin(corse), cos(tang) * cos(corse), sin(tang)]

    l2 = [-sin(tang) * sin(corse) * cos(kren) + cos(corse) * sin(kren),
          -sin(tang) * cos(corse) * cos(kren) - sin(corse) * sin(kren), cos(tang) * cos(kren)]

    l3 = [sin(tang) * sin(corse) * sin(kren) + cos(corse) * cos(kren),
          sin(tang) * cos(corse) * sin(kren) - sin(corse) * cos(kren), - cos(tang) * sin(kren)]
    return np.array([l1, l2, l3])


def gx(lat, h):
    ge = 9.78030
    beta1 = 5.302e-3
    beta2 = 7e-6
    delta_g = 14e-5
    g_ = ge * (1 + beta1 * sin(lat) ** 2 - beta2 * sin(2 * lat) ** 2 - 2 * h / a) - delta_g
    return np.array([0, 0, -g_])


def geogr_to_grinvich(lat, lon, h):
    R2 = a / sqrt(1 - e_2 * sin(lat) ** 2)
    x = (R2 + h) * cos(lon) * cos(lat)
    y = (R2 + h) * sin(lon) * cos(lat)
    z = sin(lat) * ((1 - e_2) * R2 + h)
    return [x, y, z]


def main():
    # global w_z, f_z
    # подсчет L(0)
    '''
    # осреднение угловой скорости и ньютонометров
    f = open('imu.txt', 'r').readlines()[1:]
    for line in f[:17000]:
        lin = list(map(float, line.split()))
        w_z[0] += lin[1]
        w_z[1] += lin[2]
        w_z[2] += lin[3]
        f_z[0] += lin[4]
        f_z[1] += lin[5]
        f_z[2] += lin[6]
    w_z = np.array(w_z) / 17000
    f_z = np.array(f_z) / 17000
    tmp_v = np.cross(w_z, f_z)
    v1 = tmp_v / np.linalg.norm(tmp_v)
    v2 = np.cross(f_z, tmp_v) / np.linalg.norm(np.cross(f_z, tmp_v))
    v3 = f_z / np.linalg.norm(f_z)
    L_0 = np.transpose(np.array([v1, v2, v3]))
    '''

    #  интегрирование координат и скорости
    imu = open('imu.txt', 'r').readlines()[1:]
    trj = open('trj.txt', 'r').readlines()[1:]
    lat_i, lon_i, h_i, V1_i, V2_i, V3_i = lat_0, lon_0, h_0, V1_0, V2_0, V3_0
    Vx_i = np.array([V1_0, V2_0, V3_0])

    lats, lons, alts = [], [], []
    lats_real, lons_real, alts_real = [], [], []
    x, y, z = [], [], []
    x_real, y_real, z_real = [], [], []
    V1, V2, V3 = [], [], []
    V1_real, V2_real, V3_real = [], [], []
    for i in range(len(imu)):
        # читаю данные
        l_imu = list(map(float, imu[i].split()))
        fz_i = np.array(l_imu[4:])
        l_trj = list(map(float, trj[i].split()))

        # добавляю координаты и скорости в массив
        lats.append(lat_i * 180 / pi)
        lons.append(lon_i * 180 / pi)
        alts.append(h_i)
        x_, y_, z_ = geogr_to_grinvich(lat_i * 180 / pi, lon_i * 180 / pi, h_i)
        x.append(x_)
        y.append(y_)
        z.append(z_)
        lat, lon, h = l_trj[1:4]
        lats_real.append(lat)
        lons_real.append(lon)
        alts_real.append(h)
        x_r, y_r, z_r = geogr_to_grinvich(lat, lon, h)
        x_real.append(x_r)
        y_real.append(y_r)
        z_real.append(z_r)
        V1.append(V1_i)
        V2.append(V2_i)
        V3.append(V3_i)
        V = l_trj[4:7]
        V1_real.append(V[0])
        V2_real.append(V[1])
        V3_real.append(V[2])

        #  просчитываю константы
        RN_i = RN(lat_i)
        RE_i = RE(lat_i)
        Omega_roof_i = roof(Omega(V1_i, V2_i, lat_i, h_i))
        u_roof_i = roof(ux(lat_i))
        gx_i = gx(lat_i, h_i)
        course_i, tang_i, kren_i = l_trj[7] * pi / 180, l_trj[8] * pi / 180, l_trj[9] * pi / 180
        L_i_T = np.transpose(L(course_i, tang_i, kren_i))

        #  метод эйлера
        lat_i_1 = lat_i + dt * V2_i / (RN_i + h_i)
        lon_i_1 = lon_i + dt * V1_i / ((RE_i + h_i) * cos(lat_i))
        h_i_1 = h_i + dt * V3_i
        Vx_i_1 = Vx_i + dt * ((Omega_roof_i + 2 * u_roof_i) @ Vx_i + gx_i + L_i_T @ fz_i)

        # обновляю переменные
        lat_i, lon_i, h_i = lat_i_1, lon_i_1, h_i_1
        Vx_i = Vx_i_1
        V1_i, V2_i, V3_i = Vx_i_1[0], Vx_i_1[1], Vx_i_1[2]

    # траектория и график скорости
    def draw3D(v1, v2, v3, v1_real, v2_real, v3_real, text):
        plt.figure()
        ax = plt.axes(projection='3d')
        ax.plot3D(v1, v2, v3, 'r', label=f'real and integral {text}')
        ax.plot3D(v1_real, v2_real, v3_real)
        plt.legend()
        plt.savefig(text)

    draw3D(x, y, z, x_real, y_real, z_real, 'grinvich coords')
    draw3D(V1, V2, V3, V1_real, V2_real, V3_real, 'speeds')
    draw3D(lats, lons, alts, lats_real, lons_real, alts_real, 'geograph coords')

    def draw(v, v_real, label):
        plt.figure()
        plt.plot(range(len(v)), np.array(v) - np.array(v_real), label=label)
        plt.legend()
        plt.savefig(label)

    draw(x, x_real, 'x difference')
    draw(y, y_real, 'y difference')
    draw(z, z_real, 'z difference')

    draw(lats, lats_real, 'lats difference')
    draw(lons, lons_real, 'lons difference')
    draw(alts, alts_real, 'alts difference')

    draw(V1, V1_real, 'V1 difference')
    draw(V2, V2_real, 'V2 difference')
    draw(V3, V3_real, 'V3 difference')

    plt.figure()
    plt.plot(range(len(x)), x, label='x and x_real coords')
    plt.plot(range(len(x_real)), np.array(x_real), 'r')
    plt.legend()
    plt.show()


main()
