#include <iostream>
#include <fstream>
#include <string>     // для std::getline
#include <cmath>
#include <vector>
#include <numeric>
#include <algorithm>
#include <map>
#include <list>

class Point3D {
public:
    // Координаты вершины
    double X;
    double Y;
    double Z;
    // Кол-во точек возле вершины
    int NumberOfPoint;
    // Имя точки
    std::string Name;
    Point3D(double x, double y, double z, std::string name) {
        X = x;
        Y = y;
        Z = z;
        Name = name;
    }
    // Метод для установки кол-ва точек
    void setNumberOfPoint(int numberofpoints)
    {
        NumberOfPoint = numberofpoints;
    }
    // Переопределяем оператор для возможности сортировки
    bool operator< (const Point3D& p) {
        return NumberOfPoint < p.NumberOfPoint;
    }
};

double get_distance(Point3D p1, Point3D p2) {
    return std::sqrt((p2.X - p1.X) * (p2.X - p1.X) + (p2.Y - p1.Y) * (p2.Y - p1.Y) + (p2.Z - p1.Z) * (p2.Z - p1.Z));
}

double get_sum_distance(std::vector<double> X,
                        std::vector<double> Y,
                        std::vector<double> Z,
                        Point3D MN, Point3D M0) {
    double max_dist = get_distance(MN, M0);
    std::vector<double> distances;
    for (size_t i = 0; i < X.size(); i++) {
        double dist = get_distance(MN, Point3D(X[i], Y[i], Z[i], "M"));
        if (dist < max_dist / 2)
            distances.push_back(1);
    }
    return std::accumulate(distances.begin(), distances.end(), 0);
}

size_t get_sum_count(std::vector<double> X,
                     std::vector<double> Y,
                     std::vector<double> Z,
                     Point3D MN, Point3D M0) {
    double max_dist = get_distance(MN, M0);
    size_t N = 0;
    for (size_t i = 0; i < X.size(); i++) {
        double dist = get_distance(MN, Point3D(X[i], Y[i], Z[i], "M"));
        if (dist < max_dist / 2)
            N += 1;
    }
    return N;
}

std::vector<double> sum(std::vector<double> a, std::vector<double> b) {
    std::vector<double> res = a;
    for (size_t i = 0; i < a.size(); i++) {
        res[i] += b[i];
    }
    return res;
}

std::vector<double> substr(std::vector<double> a, std::vector<double> b) {
    std::vector<double> res = a;
    for (size_t i = 0; i < a.size(); i++) {
        res[i] -= b[i];
    }
    return res;
}

std::vector<double> cross(std::vector<double> a, std::vector<double> b) {
    double x1 = a[1] * b[2] - a[2] * b[1];
    double y1 = a[2] * b[0] - a[0] * b[2];
    double z1 = a[0] * b[1] - a[1] * b[0];
    std::vector<double> res = {x1, y1, z1};
    return res;
}

std::vector<double> np_divide(std::vector<double> a, std::vector<double> b) {
    std::vector<double> res = a;
    for (size_t i = 0; i < a.size(); i++) {
        res[i] /= b[i];
    }
    return res;
}

double np_linalg_norm(std::vector<double> a) {
    double res = 0;
    for (size_t i = 0; i < a.size(); i++) {
        res += std::pow(a[i], 2);
    }
    return std::sqrt(std::abs(res));
}

int np_dot(std::vector<double> v1, std::vector<double> v2) {
    double product = 0;
    for (size_t i = 0; i < v1.size(); i++)
        product += v1[i] * v2[i];
    return product;
}

double lineseg_dist(std::vector<double> p, std::vector<double> a, std::vector<double> b) {
    std::vector<double> AB = substr(b, a);
    std::vector<double> AC = substr(p, a);
    double area = np_linalg_norm(cross(AB,AC));
    double CD = area / np_linalg_norm(AB);
    return CD;
}

int main() {
    std::vector<int> X2D;
    std::vector<int> Y2D;
    std::vector<double> X3D;
    std::vector<double> Y3D;
    std::vector<double> Z3D;
    std::vector<int> R;
    std::vector<int> G;
    std::vector<int> B;
    std::vector<int> C;
    int a, b, f, g, h, i;
    double c, d, e;
    
    

    //std::string line;

    std::ifstream in("points.txt"); // окрываем файл для чтения
    if (in.is_open()) {
        while (in >> a >> b >> c >> d >> e >> f >> g >> h >> i) {
            //std::cout << line << std::endl;
           
            // обработать line
            X2D.push_back(a); 
            Y2D.push_back(b);
            X3D.push_back(c);
            Y3D.push_back(d);
            Z3D.push_back(e);
            R.push_back(f);
            G.push_back(g);
            B.push_back(h);
            C.push_back(i);
        }
    }
    in.close();     // закрываем файл

    std::vector<double> X;
    std::vector<double> Y;
    std::vector<double> Z;
    const int CLUSTER_NR = 3;
    for (size_t i = 0; i < C.size(); i++) {
        if (C[i] == CLUSTER_NR) {
            X.push_back(X3D[i]);
            Y.push_back(Y3D[i]);
            Z.push_back(Z3D[i]);
        }
    }

    for (size_t i = 0; i < X.size(); i++) {
        std::cout << X[i] << '\n';
    }
    std::cout << *max_element(X.begin(), X.end());
    
    double Xmin = *min_element(X.begin(), X.end());
    double Xmax = *max_element(X.begin(), X.end());

    double Ymin = *min_element(Y.begin(), Y.end());
    double Ymax = *max_element(Y.begin(), Y.end());

    double Zmin = *min_element(Z.begin(), Z.end());
    double Zmax = *max_element(Z.begin(), Z.end());

    double X0 = Xmin + (Xmax - Xmin) / 2.0;
    double Y0 = Ymin + (Ymax - Ymin) / 2.0;
    double Z0 = Zmin + (Zmax - Zmin) / 2.0;

    Point3D M0 = Point3D(X0, Y0, Z0, "M0");
    Point3D M1 = Point3D(Xmin, Ymin, Zmin, "M1");
    Point3D M2 = Point3D(Xmax, Ymin, Zmin, "M2");
    Point3D M3 = Point3D(Xmax, Ymax, Zmin, "M3");
    Point3D M4 = Point3D(Xmin, Ymax, Zmin, "M4");
    Point3D M5 = Point3D(Xmin, Ymin, Zmax, "M5");
    Point3D M6 = Point3D(Xmax, Ymin, Zmax, "M6");
    Point3D M7 = Point3D(Xmax, Ymax, Zmax, "M7");
    Point3D M8 = Point3D(Xmin, Ymax, Zmax, "M8");

    size_t m1 = get_sum_count(X, Y, Z, M1, M0);
    size_t m2 = get_sum_count(X, Y, Z, M2, M0);
    size_t m3 = get_sum_count(X, Y, Z, M3, M0);
    size_t m4 = get_sum_count(X, Y, Z, M4, M0);
    size_t m5 = get_sum_count(X, Y, Z, M5, M0);
    size_t m6 = get_sum_count(X, Y, Z, M6, M0);
    size_t m7 = get_sum_count(X, Y, Z, M7, M0);
    size_t m8 = get_sum_count(X, Y, Z, M8, M0);

#pragma region Расчет осевых точек
    // Устанавливаем кол-во точек около габаритной точки
    M1.setNumberOfPoint(m1);
    M2.setNumberOfPoint(m2);
    M3.setNumberOfPoint(m3);
    M4.setNumberOfPoint(m4);
    M5.setNumberOfPoint(m5);
    M6.setNumberOfPoint(m6);
    M7.setNumberOfPoint(m7);
    M8.setNumberOfPoint(m8);

    // Формируем список
    std::list<Point3D> M;
    M.push_back(M1);
    M.push_back(M2);
    M.push_back(M3);
    M.push_back(M4);
    M.push_back(M5);
    M.push_back(M6);
    M.push_back(M7);
    M.push_back(M8);
    // Вывод неотсортированного списка
    std::cout << std::endl;
    std::cout << std::endl << "M[] points:" << std::endl;
    for (auto it = M.begin(); it != M.end(); it++)
        std::cout << it->Name << ": [" << it->NumberOfPoint << "]" << std::endl;
        
    // Сортировка
    M.sort();
    M.reverse();
    // Вывод отсортированного списка
    std::cout << std::endl;
    std::cout << std::endl << "M[] points (sorted):" << std::endl;
    for (auto it = M.begin(); it != M.end(); it++)
        std::cout << it->Name << ": [" << it->NumberOfPoint << "]" << std::endl;

    std::cout << std::endl;
    // Формируем список осевых точек
    std::list<Point3D> P;
    // Первая точка берется из отсортированного списка MD
    P.push_back(M.front());
    auto& P0 = M.front();

    double m0_dist;
    double p0_dist;
    bool same_x;
    bool same_y;
    bool same_z;
    bool same_xyz;

    // Начинаем проверку со второй точки
    for (auto it = std::next(M.begin()); it != M.end(); ++it)
    {
        m0_dist = get_distance(*it, M0);  // Расстояние от текущей точки до центра масс
        p0_dist = get_distance(*it, P0);  // Расстояние между точками P0 и текущей

        // Проверяем, не лежат ли точки в одной плоскости
        same_x = it->X == P0.X;
        same_y = it->Y == P0.Y;
        same_z = it->Z == P0.Z;
        same_xyz = !(same_x || same_y || same_z);
        // if (m0_dist < p0_dist)
        if ((m0_dist < p0_dist) && same_xyz)
        {
            // Добавляем найденную точку в список
            P.push_back(*it);
            // Если нужная точка найдена, прекращаем перебор
            break;
        }   
    }
    auto& P1 = P.back();

    // Выводим найденые осевые точки на экран
    std::cout << std::endl;
    std::cout << P0.Name << "[" << P0.NumberOfPoint << "] (" << P0.X << ", " << P0.Y << ", " << P0.Z << ")" << std::endl;
    std::cout << P1.Name << "[" << P1.NumberOfPoint << "] (" << P1.X << ", " << P1.Y << ", " << P1.Z << ")" << std::endl;
#pragma endregion

    // Вычисление ширины
    std::vector<double> dists;
    for (size_t i = 0; i < X.size(); i++) {
        //
    }

}
