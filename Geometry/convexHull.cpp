#include <bits/stdc++.h>

using namespace std;

struct Point
{
    int x, y;
    Point() {}
    Point(int _x, int _y) : x(_x), y(_y) {}
    int dist(const Point& other) const
    {
        return (x - other.x) * (x - other.x) + (y - other.y) * (y - other.y);
    }
};

struct Vector
{
    int x, y;
    Vector() {}
    Vector(int _x, int _y) : x(_x), y(_y) {}
    double cross_product(const Vector& other)
    {
        return x * other.y - y * other.x;
    }
};

Vector operator - (const Point& a, const Point& b)
{
    return Vector(a.x - b.x, a.y - b.y);
}

// Xác định hướng đi từ A -> B -> C.
int orientation(const Point& a, const Point& b, const Point& c)
{
    Vector x = b - a;
    Vector y = c - b;
    int orient = x.cross_product(y);

    if (orient < 0) // Cùng chiều kim đồng hồ.
        return -1;
    else if (orient == 0) // Thẳng hàng.
        return 0;
    return 1;
}

bool collinear(const Point& a, const Point& b, const Point& c)
{
    return orientation(a, b, c) == 0;
}

bool cw(const Point& a, const Point& b, const Point& c, bool include_collinear)
{
    int orient = orientation(a, b, c);
    return (orient < 0) || (orient == 0 && include_collinear);
}

vector < Point > convex_hull(const int& n, vector < Point >& p, bool include_collinear = false)
{
    Point p0 = *min_element(p.begin() + 1, p.end(), [](Point a, Point b)
    {
        return make_pair(a.y, a.x) < make_pair(b.y, b.x);
    });

    sort(p.begin() + 1, p.end(), [&p0](const Point& a, const Point& b)
    {
        int orient = orientation(p0, a, b);

        // Góc cầu hình của A = Góc cầu hình của B thì sắp xếp theo khoảng cách tới P0.
        if (orient == 0)
            return p0.dist(a) < p0.dist(b);
        // Thứ tự mong muốn là Góc cầu hình của A < Góc cầu hình của B.
        // Nếu P0 -> A -> B theo chiều kim đồng hồ tức là B ở dưới A trên mặt phẳng, thì
        // góc cầu hình của A sẽ lớn hơn góc cầu hình của B.
        return orient < 0;
    });

    // Trường hợp bao lồi yêu cầu phải chứa các điểm thẳng hàng.
    if (include_collinear)
    {
        int i = n;
        while (i >= 1 && collinear(p0, p[i], p.back()))
            --i;

        reverse(p.begin() + i + 1, p.end());
    }

    vector < Point > hull;
    for (int i = 1; i < (int) p.size(); ++i)
    {
        while (hull.size() > 1 && !cw(hull[hull.size() - 2], hull.back(), p[i], include_collinear))
            hull.pop_back();

        hull.push_back(p[i]);
    }

    return hull;
}

int min_edge_length(const vector<Point>& hull) {
    int min_length = INT_MAX;

    for (size_t i = 0; i < hull.size(); ++i) {
        // Tính khoảng cách giữa hai điểm liên tiếp (vòng lại điểm đầu nếu i là điểm cuối)
        int dist = hull[i].dist(hull[(i + 1) % hull.size()]);
        min_length = min(min_length, dist);
    }

    return min_length;
}

double convex_hull_area(const vector<Point>& hull) {
    double area = 0.0;

    for (size_t i = 0; i < hull.size(); ++i) {
        Point p1 = hull[i];
        Point p2 = hull[(i + 1) % hull.size()]; // Điểm tiếp theo, vòng lại điểm đầu nếu là điểm cuối
        area += (p1.x * p2.y - p1.y * p2.x);
    }

    return fabs(area) / 2.0; // Trả về giá trị tuyệt đối và chia đôi
}

signed main()
{
    ios_base::sync_with_stdio(false);
    cin.tie(nullptr);

    int n;
    cout << "Enter number of points : ";
    cin >> n;

    vector < Point > p(n + 1);
    for (int i = 1; i <= n; ++i)
        cin >> p[i].x >> p[i].y;

    cout << "convexHull: " << endl;
    vector < Point > hull = convex_hull(n, p);
    for (auto p: hull)
        cout << p.x << ' ' << p.y << endl;
    cout << "Minimum edge length: " << min_edge_length(hull) << endl;
    cout << "Convex hull area: " << convex_hull_area(hull) << endl;
    return 0;
}
