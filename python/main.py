import pandas as pd
import matplotlib.pyplot as plt
import math
import numpy as np

class Point3D:
    def __init__(self, x, y, z, name):
        self.X = x
        self.Y = y
        self.Z = z
        self.Name = name


def equation_plane(p1: Point3D, p2: Point3D, p3: Point3D):
    a1 = p2.X - p1.X
    b1 = p2.Y - p1.Y
    c1 = p2.Z - p1.Z
    a2 = p3.X - p1.X
    b2 = p3.Y - p1.Y
    c2 = p3.Z - p1.Z
    a = b1 * c2 - b2 * c1
    b = a2 * c1 - a1 * c2
    c = a1 * b2 - b1 * a2
    d = (- a * p1.X - b * p1.Y - c * p1.Z)
    return a, b, c, d


def get_distance(p1: Point3D, p2: Point3D):
    return math.sqrt((p2.X - p1.X) * (p2.X - p1.X) + (p2.Y - p1.Y) * (p2.Y - p1.Y) + (p2.Z - p1.Z) * (p2.Z - p1.Z))


def get_sum_distance(X, Y, Z, MN: Point3D, M0: Point3D):
    max_dist = get_distance(MN, M0)
    XX = X.tolist();
    YY = Y.tolist();
    ZZ = Z.tolist();
    distances = []
    for i in range(len(X)):
        dist = get_distance(MN, Point3D(XX[i], YY[i], ZZ[i], "M"))
        if dist < max_dist/2:
            distances.append(1)

    return sum(distances)


def get_sum_count(X, Y, Z, MN: Point3D, M0: Point3D):
    max_dist = get_distance(MN, M0)
    XX = X.tolist()
    YY = Y.tolist()
    ZZ = Z.tolist()
    N = 0
    for i in range(len(X)):
        dist = get_distance(MN, Point3D(XX[i], YY[i], ZZ[i], "M"))
        if dist < max_dist / 2:
            N += 1

    return N


def plot3D(X, Y, Z, M0, M1, M2, M3, M4, M5, M6, M7, M8, P0, P1):
    delta = 5
    fig = plt.figure(figsize=(4, 4), dpi=300)
    ax = plt.axes(projection='3d')
    ax.scatter3D(X, Y, Z, c='Blue')

    ax.scatter3D(M1.X, M1.Y, M1.Z, c='Red', alpha=0.5)
    ax.text(M1.X + delta, M1.Y + delta, M1.Z + delta, 'M1', c='Red', fontsize=10, alpha=0.7)

    ax.scatter3D(M2.X, M2.Y, M2.Z, c='Red', alpha=0.5)
    ax.text(M2.X + delta, M2.Y + delta, M2.Z + delta, 'M2', c='Red', fontsize=10, alpha=0.7)

    ax.scatter3D(M3.X, M3.Y, M3.Z, c='Red', alpha=0.5)
    ax.text(M3.X + delta, M3.Y + delta, M3.Z + delta, 'M3', c='Red', fontsize=10, alpha=0.7)

    ax.scatter3D(M4.X, M4.Y, M4.Z, c='Red', alpha=0.5)
    ax.text(M4.X + delta, M4.Y + delta, M4.Z + delta, 'M4', c='Red', fontsize=10, alpha=0.7)

    ax.scatter3D(M5.X, M5.Y, M5.Z, c='Red', alpha=0.5)
    ax.text(M5.X + delta, M5.Y + delta, M5.Z + delta, 'M5', c='Red', fontsize=10, alpha=0.7)

    ax.scatter3D(M6.X, M6.Y, M6.Z, c='Red', alpha=0.5)
    ax.text(M6.X + delta, M6.Y + delta, M6.Z + delta, 'M6', c='Red', fontsize=10, alpha=0.7)

    ax.scatter3D(M7.X, M7.Y, M7.Z, c='Red', alpha=0.5)
    ax.text(M7.X + delta, M7.Y + delta, M7.Z + delta, 'M7', c='Red', fontsize=10, alpha=0.7)

    ax.scatter3D(M8.X, M8.Y, M8.Z, c='Red', alpha=0.5)
    ax.text(M8.X + delta, M8.Y + delta, M8.Z + delta, 'M8', c='Red', fontsize=10, alpha=0.7)

    X_z_min = [Xmin, Xmax, Xmax, Xmin, Xmin]
    Y_z_min = [Ymin, Ymin, Ymax, Ymax, Ymin]
    Z_z_min = [Zmin, Zmin, Zmin, Zmin, Zmin]
    ax.plot3D(X_z_min, Y_z_min, Z_z_min, c='Red', linewidth=0.75, alpha=0.5)

    X_z_max = [Xmin, Xmax, Xmax, Xmin, Xmin]
    Y_z_max = [Ymin, Ymin, Ymax, Ymax, Ymin]
    Z_z_max = [Zmax, Zmax, Zmax, Zmax, Zmax]
    ax.plot3D(X_z_max, Y_z_max, Z_z_max, c='Red', linewidth=0.75, alpha=0.5)

    X_x_min = [Xmin, Xmin, Xmin, Xmin, Xmin]
    Y_x_min = [Ymin, Ymin, Ymax, Ymax, Ymin]
    Z_x_min = [Zmin, Zmax, Zmax, Zmin, Zmin]
    ax.plot3D(X_x_min, Y_x_min, Z_x_min, c='Blue', linewidth=0.75, alpha=0.5)

    X_x_max = [Xmax, Xmax, Xmax, Xmax, Xmax]
    Y_x_max = [Ymin, Ymin, Ymax, Ymax, Ymin]
    Z_x_max = [Zmin, Zmax, Zmax, Zmin, Zmin]
    ax.plot3D(X_x_max, Y_x_max, Z_x_max, c='Blue', linewidth=0.75, alpha=0.5)

    ax.scatter(Xo, Yo, Zo, s=20, c='Green', alpha=0.5)
    ax.text(M0.X + 3, M0.Y + 3, M0.Z + 3, 'M0', c='Green', fontsize=10, alpha=0.7)
    ax.set_title('3D Point Cloud')

    ax.scatter3D(P0.X, P0.Y, P0.Z, c='Magenta')
    ax.scatter3D(P1.X, P1.Y, P1.Z, c='Magenta')

    ax.plot3D([P0.X, P1.X], [P0.Y, P1.Y], [P0.Z, P1.Z], c='Green', linewidth=2, alpha=0.5)

    # Выравнивание масштабной сетки
    ax.set_box_aspect((np.ptp(X), np.ptp(Y), np.ptp(Z)))  # aspect ratio is 1:1:1 in data space
    plt.show()


if __name__ == '__main__':
    data = pd.read_csv('D:\\TEMP\\pointsClaster_1001.txt', sep='\t', header=None)
    data.columns = ['X2D', 'Y2D', 'X', 'Y', 'Z', 'R', 'G', 'B', 'ID']
    data.head()
    # claster #3 -- 3D object
    df = data[data.ID == 3]

    X = df.X
    Y = df.Y
    Z = df.Z

    Xmin = X.min()
    Xmax = X.max()

    Ymin = Y.min()
    Ymax = Y.max()

    Zmin = Z.min()
    Zmax = Z.max()

    Xo = Xmin + (Xmax - Xmin) / 2
    Yo = Ymin + (Ymax - Ymin) / 2
    Zo = Zmin + (Zmax - Zmin) / 2
    print(f'Xo = {Xo:.0f}; Yo = {Yo:.0f}; Zo = {Zo:.0f};')

    M0 = Point3D(Xo, Yo, Zo, "M0")
    M1 = Point3D(Xmin, Ymin, Zmin, "M1")
    M2 = Point3D(Xmax, Ymin, Zmin, "M2")
    M3 = Point3D(Xmax, Ymax, Zmin, "M3")
    M4 = Point3D(Xmin, Ymax, Zmin, "M4")
    M5 = Point3D(Xmin, Ymin, Zmax, "M5")
    M6 = Point3D(Xmax, Ymin, Zmax, "M6")
    M7 = Point3D(Xmax, Ymax, Zmax, "M7")
    M8 = Point3D(Xmin, Ymax, Zmax, "M8")

    m1 = get_sum_count(X, Y, Z, M1, M0)
    m2 = get_sum_count(X, Y, Z, M2, M0)
    m3 = get_sum_count(X, Y, Z, M3, M0)
    m4 = get_sum_count(X, Y, Z, M4, M0)
    m5 = get_sum_count(X, Y, Z, M5, M0)
    m6 = get_sum_count(X, Y, Z, M6, M0)
    m7 = get_sum_count(X, Y, Z, M7, M0)
    m8 = get_sum_count(X, Y, Z, M8, M0)

    MD = {M1: m1, M2: m2, M3: m3, M4: m4, M5: m5, M6: m6, M7: m7, M8: m8}
    MDD = dict(sorted(MD.items(), key=lambda item: item[1], reverse=True))

    # Суммарные расстояния до точек
    print(f"M1: {m1:.2f}")
    print(f"M2: {m2:.2f}")
    print(f"M3: {m3:.2f}")
    print(f"M4: {m4:.2f}")
    print(f"M5: {m5:.2f}")
    print(f"M6: {m6:.2f}")
    print(f"M7: {m7:.2f}")
    print(f"M8: {m8:.2f}")

    print("\n")
    for keys, value in MDD.items():
        print(f"{keys.Name}: {value:.2f}")



    # Проверяем, чтобы осевые точки не оказались близки между собой
    M = list(MDD.keys()) # список узловых отсортированных точек
    P = [] # Пустой список для осевых точек
    P.append(M[0])
    print(f"\nP0: {M[0].Name} ({M[0].X}, {M[0].Y}, {M[0].Z})")
    for i in range(1, 8):
        m0_dist = get_distance(M[i], M0)  # до центра масс
        pp_dist = get_distance(M[i], M[0])  # между осевыми точками M[0] и текущей

        same_x = math.isclose(M[0].X, M[i].X)
        same_y = math.isclose(M[0].Y, M[i].Y)
        same_z = math.isclose(M[0].Z, M[i].Z)

        same_xyz = not(same_x or same_y or same_z)

        if ((m0_dist < pp_dist) and same_xyz):
            P.append(M[i])
            print(f"P1: {M[i].Name} ({M[i].X}, {M[i].Y}, {M[i].Z})")
            break

    # Осевые точки
    P0 = P[0]
    P1 = P[1]

    plot3D(X, Y, Z, M0, M1, M2, M3, M4, M5, M6, M7, M8, P0, P1)

    # Уравнение плоскости, проходящей через три точки
    # https://www.geeksforgeeks.org/program-to-find-equation-of-a-plane-passing-through-3-points/
    # M1 = Point3D(1, -2, 0)
    # M2 = Point3D(2, 0, -1)
    # M3 = Point3D(0, -1, 2)
    # a, b, c, d = equation_plane(M1, M2, M3)
    # print(f"a = {a}; b = {b}; c = {c}; d = {d}")