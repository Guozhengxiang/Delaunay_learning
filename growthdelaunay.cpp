/* 
Filename: growthdelaunay.cpp
Author: Froven
Date: 2021-07-4
Description: 使用区域增长法实现了空间曲面的Delaunay三角剖分，输入空间三维点坐标，返回三角剖分后的三角形
*/

#include <math.h>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <functional>
#include <string>
#include <vector>
#include <queue>

#define eps 1e-10

class Point3d {     // 点
 public:
    double x, y, z;
    Point3d() {}
    Point3d(double x_, double y_, double z_) :x(x_), y(y_), z(z_) {}
    friend double operator* (const Point3d& p, const Point3d& q) { return p.x * q.x + p.y * q.y + p.z * q.z; }   // 向量的点乘
    friend Point3d operator* (double s, const Point3d& p) { return Point3d(p.x * s, p.y * s, p.z * s); }    // 向量乘以常数
    friend Point3d operator/ (const Point3d& p, double s) { return Point3d(p.x / s, p.y / s, p.z / s); }    // 向量除以常数
    friend Point3d operator- (const Point3d& p, const Point3d& q) { return Point3d(p.x - q.x, p.y - q.y, p.z - q.z); } // 向量的差
    friend Point3d operator+ (const Point3d& p, const Point3d& q) { return Point3d(p.x + q.x, p.y + q.y, p.z + q.z); } // 向量的和
    friend bool operator == (const Point3d& p, const Point3d& q) { return fabs(p.x - q.x) + fabs(p.y - q.y) + fabs(p.z - q.z) <= eps; } // 两个点近似相等
};
class Edge {    // 边
 public:
    Point3d st, ed, other;
    Edge() {}
    Edge(const Point3d& p0_, const Point3d& p1_, const Point3d& ot_) : st(p0_), ed(p1_), other(ot_) {}
    friend bool operator == (const Edge& e1, const Edge& e2) { return ((e1.st == e2.st) && (e1.ed == e2.ed)) || ((e1.st == e2.ed) && (e1.ed == e2.st)); }   // 判断两条边是否相等
};
class Triangle {    // 三角形
 public:
    Point3d p0, p1, p2;
    Triangle() {}
    Triangle(const Point3d& p0_, const Point3d& p1_, const Point3d& p2_) : p0(p0_), p1(p1_), p2(p2_) {}
    friend bool operator==(const Triangle& left, const Triangle& right) {      // 判断两三角形是否相等
        return (left.p0 == right.p0 || left.p0 == right.p1 || left.p0 == right.p2)
            && (left.p1 == right.p0 || left.p1 == right.p1 || left.p1 == right.p2)
            && (left.p2 == right.p0 || left.p2 == right.p1 || left.p2 == right.p2);
    }
};
class Plane {    // 平面
 public:
    double a, b, c, d;    // 满足平面方程的参数ax + bx + cz + d = 0
    Point3d norm;
    Plane() {}
    Plane(double a_, double b_, double c_, double d_) :a(a_), b(b_), c(c_), d(d_) { norm = Point3d(a, b, c); }
    explicit Plane(const Triangle& t) {
        a = t.p0.y * (t.p1.z - t.p2.z) + t.p1.y * (t.p2.z - t.p0.z) + t.p2.y * (t.p0.z - t.p1.z);
        b = t.p0.z * (t.p1.x - t.p2.x) + t.p1.z * (t.p2.x - t.p0.x) + t.p2.z * (t.p0.x - t.p1.x);
        c = t.p0.x * (t.p1.y - t.p2.y) + t.p1.x * (t.p2.y - t.p0.y) + t.p2.x * (t.p0.y - t.p1.y);
        d = -t.p0.x * (t.p1.y * t.p2.z - t.p2.y * t.p1.z) - t.p1.x * (t.p2.y * t.p0.z - t.p0.y * t.p2.z) - t.p2.x * (t.p0.y * t.p1.z - t.p1.y * t.p0.z);
        norm = Point3d(a, b, c);
    }
};

struct hashfunc_edge {
    size_t operator() (const Edge& e) const {
        return  std::hash<double>()(e.st.x) ^ std::hash<double>()(e.st.y) ^ std::hash<double>()(e.st.z)^
                std::hash<double>()(e.ed.x) ^ std::hash<double>()(e.ed.y) ^ std::hash<double>()(e.ed.z);
    }
};
struct equalkey_edge {
    bool operator () (const Edge& le, const Edge& re) const {
        return le == re;
    }
};

double norm3d2(const Point3d& p, const Point3d& q) {   // 两点距离
    return (p.x - q.x) * (p.x - q.x) + (p.y - q.y) * (p.y - q.y) + (p.z - q.z) * (p.z - q.z);
}
double norm3d2(const Point3d& p) {   // 向量的模
    return p.x * p.x + p.y * p.y + p.z * p.z;
}
double calc_cos_angle(const Point3d& p, const Point3d& q) {
    return (p * q) / (std::sqrt(norm3d2(p) * norm3d2(q)));
}

/* 参数载入 */
bool part_params_load(int& N, const std::string& k_input_filename, std::vector<Point3d>& points_3d) {
    FILE* fp;
    auto err = fopen_s(&fp, k_input_filename.c_str(), "r");
    if (err) {
        printf_s("The file was not opened!!\n");
        exit(0);
    } else {
        fscanf_s(fp, "%d\n", &N);
        points_3d.resize(N);
        for (int i = 0; i < N; i++) {
            fscanf_s(fp, "%lf %lf %lf\n", &points_3d[i].x, &points_3d[i].y, &points_3d[i].z);
        }
    }
    fclose(fp);
    return 1;
}

/* 生长法构建三角网 */
void creat_tri_net(const int N, std::vector<Point3d>& points, std::vector<Triangle>& triangles) {

    //----------------生长法算法初始化------------------//
    std::vector<Edge> edges;    // 定义外轮廓边
    Point3d p1, p2, p3;
    // 构建初始较为饱满的三角形
    p1 = points[0];
    // 确定第二个点（距离最小的点）
    double dis = 1e5;
    int flag2 = 1;
    for (int i = 1; i < N; i++) {
        double t_dis = norm3d2(p1, points[i]);
        if (t_dis < dis) {
            dis = t_dis;
            flag2 = i;
        }
    }
    p2 = points[flag2];
    // 确定第三个点（邻域内角度最大的点, cos值最小）
    double min_cos = 1;
    for (auto& p : points) {
        if (!(p == p1) && !(p == p2)) {
            double tmp = calc_cos_angle(p1 - p, p2 - p);
            if (tmp < min_cos) {
                min_cos = tmp;
                p3 = p;
            }
        }
    }

    triangles.emplace_back(Triangle(p1, p2, p3));
    std::queue<Edge> border;  // 定义边界边
    std::vector<Edge> current_border;
    // 将首个三角形的边加入
    border.emplace(Edge(p1, p2, p3));
    border.emplace(Edge(p1, p3, p2));
    border.emplace(Edge(p3, p2, p1));

    std::unordered_map<Edge, int, hashfunc_edge, equalkey_edge> tmp_edge;  // 定义下一次遍历的边界边,值为对应的三角形数量
    // ----------------------------------------------------- //
    // ---------------边的队列不空时,开始扩张--------------- //
    while (!border.empty()) {   //
        const auto& e = border.front();
        border.pop();

        // 求过该边的三角形顶点到该边的垂向量
        Point3d vertical_ver = (e.other - e.st) - ((e.other - e.st) * (e.ed - e.st)) * (e.ed - e.st) / norm3d2(e.ed - e.st);

        // 计算异侧、角度最大（cos值最小）的点
        double min_cos = 1;
        Point3d *nextpoint = nullptr;
        for (auto& p : points) {
            if ((p - e.st) * vertical_ver < 0) {   // 计算异侧的点
                double tmp = calc_cos_angle(e.st - p, e.ed - p);  // 计算角度
                if (min_cos > tmp) {
                    min_cos = tmp;
                    nextpoint = &p;
                }
            }
        }
        if (nextpoint != nullptr) {
            // 若三角形集合没有该三角形，则将新生成的三角形加入三角形集合
            Triangle tm_tr(*nextpoint, e.st, e.ed);
            bool tm_tr_repeat = false;
            for (const auto& tr : triangles) {
                if (tm_tr == tr) tm_tr_repeat = true;
            }
            // 对于不同的三角形
            if (!tm_tr_repeat) {
                triangles.emplace_back(tm_tr);
                // 计算新生成的边界边对应的三角形数量
                Edge tmp1 = Edge(*nextpoint, e.st, e.ed);
                Edge tmp2 = Edge(*nextpoint, e.ed, e.st);
                bool tmp1_repeat = false, tmp2_repeat = false;
                for (auto& tm_e : current_border) {
                    if (tmp1 == tm_e) tmp1_repeat = true;
                    if (tmp2 == tm_e) tmp2_repeat = true;
                }
                if (!tmp1_repeat) tmp_edge[tmp1] += 1;
                if (!tmp2_repeat) tmp_edge[tmp2] += 1;
            }
        }
        // 上一组边界边遍历完后，将新的边界边加入队列
        if (border.empty()) {
            current_border.clear();
            for (const auto& te : tmp_edge) {
                if (te.second == 1) {
                    current_border.emplace_back(te.first);
                    border.emplace(te.first);
                }
            }
            tmp_edge.clear();
        }
    }
}

/* 函数作用：按格式输出路径信息，输出失败则直接退出 */
bool data_to_txt(const std::vector<Triangle>& face, const std::string k_output_filename) {
    std::ofstream file;
    file.open(k_output_filename, std::ios::out);
    if (!file) return 0;
    int size = face.size();
    file << size << std::endl;
    for (int i = 0; i < size; i++) {
        file << face[i].p0.x << " " << face[i].p0.y << " " << face[i].p0.z << " "
             << face[i].p1.x << " " << face[i].p1.y << " " << face[i].p1.z << " "
             << face[i].p2.x << " " << face[i].p2.y << " " << face[i].p2.z << std::endl;
    }
    file.close();
    return 1;
}

int main(int argc, char* argv[]) {
    const std::string k_input_filename {argv[1]};   // 原始数据文件路径
    const std::string k_output_filename {argv[2]};  // 结果文件路径

    int N;   // 点集数量
    std::vector<Point3d> points_3d;   // 输入的点集
    std::vector<Triangle> triangles;  // 输出的三角形集
    if (!part_params_load(N, k_input_filename, points_3d))return 0;    // 数据载入
    creat_tri_net(N, points_3d, triangles);    // 生成局部三角网
    data_to_txt(triangles, k_output_filename);    // 写入到文件
    return 0;
}
