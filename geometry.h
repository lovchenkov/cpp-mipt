#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>



namespace geom {
  const double EPS = 1e-6;
  const double PI = 3.14159265359;
  const double PI_DEGREE = 180;
}

bool equal(double x, double y) { return fabs(x - y) < geom::EPS; }

struct Point {
  double x;
  double y;
  Point();
  Point(double x, double y) : x(x), y(y) {}

  bool operator==(const Point& other) const { return (equal(x, other.x) && equal(y, other.y)); }
  bool operator!=(const Point& other) const {return !(*this == other); }
};

std::ostream& operator<<(std::ostream& out, const Point& p) {
  out << p.x << ' ' << p.y;
  return out;
}

double dist(Point p1, Point p2) {
  return sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y));
}

double cross_product(Point p1, Point p2, Point p3, Point p4) {  // [p1p2, p3p4]
  return (p2.x - p1.x) * (p4.y - p3.y) - (p4.x - p3.x) * (p2.y - p1.y);
}

double signed_angle(Point p1, Point p2, Point p3, Point p4) {
  return atan2(cross_product(p1, p2, p3, p4), (p2.x - p1.x) * (p4.x - p3.x) + (p2.y - p1.y) * (p4.y - p3.y));
}

double angle(Point p1, Point p2, Point p3, Point p4) {
  return fabs(signed_angle(p1, p2, p3, p4));
}


Point middle(Point p1, Point p2) { return Point((p1.x + p2.x) / 2, (p1.y + p2.y) / 2); }

Point rotation(Point p1, Point p2, double alpha) {  // rotates p2 over p1.
  return Point(p1.x + (p2.x - p1.x) * cos(alpha * geom::PI / geom::PI_DEGREE) - (p2.y - p1.y) * sin(alpha * geom::PI / geom::PI_DEGREE),
               p1.y + (p2.x - p1.x) * sin(alpha * geom::PI / geom::PI_DEGREE) + (p2.y - p1.y) * cos(alpha * geom::PI / geom::PI_DEGREE));
}

Point reflection(Point p1, Point p2) {  // reflects p2 over p1
  return Point(2 * p1.x - p2.x, 2 * p1.y - p2.y);
}

Point homothety(Point center, Point p, double coefficient) {
  return Point(center.x + (p.x - center.x) * coefficient,
               center.y + (p.y - center.y) * coefficient);
}


class Line {
 public:
  double a;  // line is ax+by+c=0
  double b;
  double c;

  Line(Point p1, Point p2)
      : a(p1.y - p2.y),
        b(p2.x - p1.x),
        c(-p1.y * b - p1.x * a) {}
  Line(double tg, double shift) : Line(Point(shift, 0), Point(shift + 1, tg)) {}
  Line(Point p, double tg) : Line(p, Point(p.x + 1, p.y + tg)) {}
  Line(double a, double b, double c) : a(a), b(b), c(c) {}
};

std::pair<Point, Point> move_point_on_line(Point p, Line l, double dist) {
  double t = dist / sqrt((l.a * l.a + l.b * l.b));
  return {Point(p.x - l.b * t, p.y + l.a * t), Point(p.x + l.b * t, p.y - l.a * t)};
}

Line perpendicular_through_point(Point p, Line l) {
  return Line(-l.b, l.a, l.b * p.x - l.a * p.y);
}

double distance_to_line(Point p, Line l) { return fabs((l.a * p.x + l.b * p.y + l.c) /
                                                       sqrt(l.a * l.a + l.b * l.b)); }

bool operator==(const Line& l1, const Line& l2) {
  double s = l1.a + l1.b + l1.c;
  double other_s = l2.a + l2.b + l2.c;
  return (equal(l1.a / s, l2.a / other_s) &&
          equal(l1.b / s, l2.b / other_s) &&
          equal(l1.c / s, l2.c / other_s));
}
bool operator!=(const Line& l1, const Line& l2) { return !(l1 == l2); }

Point intersection(Line l1, Line l2) {
  return Point(-(l1.c * l2.b - l1.b * l2.c) / (l1.a * l2.b - l1.b * l2.a),
               -(l1.a * l2.c - l1.c * l2.a) / (l1.a * l2.b - l1.b * l2.a));
}

Line bisector_line(Point p1, Point p2) {
  return perpendicular_through_point(middle(p1, p2), Line(p1, p2));
}

Point segment_division(Point p1, Point p2, double ratio) {
  return Point((p1.x - ratio * p2.x) / (1 - ratio), (p1.y - ratio * p2.y) / (1 - ratio));
}

Point reflection_over_line(Line l, Point p) {
  Point height = intersection(l, perpendicular_through_point(p, l));
  return reflection(height, p);
}












class Shape {
 public:
  virtual double area() const = 0;
  virtual double perimeter() const = 0;
  virtual bool isEqual(const Shape& other) const = 0;
  virtual bool isCongruentTo(const Shape& other) const = 0;
  virtual bool isSimilarTo(const Shape& other) const = 0;
  virtual bool containsPoint(const Point& p) const = 0;
  virtual void rotate(const Point& center, double angle) = 0;
  virtual void reflect(const Point& center) = 0;
  virtual void reflect(const Line& axis) = 0;
  virtual void scale(const Point& center, double coefficient) = 0;
  virtual ~Shape() = default;
};

bool operator==(const Shape& s1, const Shape& s2) {
  return s1.isEqual(s2);
}

class Ellipse : public Shape {
 protected:
  Point f1;
  Point f2;
  double radius;

 public:
  Ellipse(Point p1, Point p2, double radius)
      : f1(p1),
        f2(p2),
        radius(radius) {}
  std::pair<Point, Point> focuses() const { return {f1, f2}; }
  double eccentricity() const { return dist(f1, f2) / radius; }
  std::pair<Line, Line> directrices() const;
  double small_axis() const { return sqrt(radius * radius / 4 - dist(f1, f2) * dist(f1, f2) / 4); }
  double big_axis() const  {return radius / 2; }
  Point center() const { return middle(f1, f2); }

  double area() const override { return geom::PI * small_axis() * big_axis(); }
  double perimeter() const override { return geom::PI * (3 * small_axis() + 3 * big_axis() -
                                                         sqrt((3 * small_axis() + big_axis()) * (3 * big_axis() + small_axis()))); }
  bool isEqual(const Shape& other) const override {
    const Ellipse* other_ptr = dynamic_cast<const Ellipse*>(&other);
    if (other_ptr == nullptr) {
      return false;
    }
    if (!equal(radius, other_ptr->radius)) {
      return false;
    }
    if ((f1 == other_ptr->f1 && f2 == other_ptr->f2) || (f1 == other_ptr->f2 && f2 == other_ptr->f1)) {
      return true;
    }
    return false;
  }
  bool isCongruentTo(const Shape& other) const override {
    const Ellipse* other_ptr = dynamic_cast<const Ellipse*>(&other);
    if (other_ptr == nullptr) {
      return false;
    }
    return (equal(dist(f1, f2), dist(other_ptr->f1, other_ptr->f2)) && equal(radius, other_ptr->radius));
  }
  bool isSimilarTo(const Shape& other) const override {
    const Ellipse* other_ptr = dynamic_cast<const Ellipse*>(&other);
    if (other_ptr == nullptr) {
      return false;
    }
    if (equal(dist(f1, f2), dist(other_ptr->f1, other_ptr->f2))) {
      return isCongruentTo(*other_ptr);
    }
    return (equal(dist(f1, f2) / dist(other_ptr->f1, other_ptr->f2), radius / other_ptr->radius));
  }
  bool containsPoint(const Point& p) const override {
    return (dist(p, f1) + dist(p, f2) - radius < geom::EPS);
  }
  void rotate(const Point& center, double angle) override {
    f1 = rotation(center, f1, angle);
    f2 = rotation(center, f2, angle);
  }
  void reflect(const Point& center) override {
    f1 = reflection(center, f1);
    f2 = reflection(center, f2);
  }
  void reflect(const Line& axis) override {
    f1 = reflection_over_line(axis, f1);
    f2 = reflection_over_line(axis, f2);
  }
  void scale(const Point& center, double coefficient) override {
    f1 = homothety(center, f1, coefficient);
    f2 = homothety(center, f2, coefficient);
    radius *= coefficient;
  }
};


std::pair<Line, Line> Ellipse::directrices() const {
  auto points = move_point_on_line(center(), Line(f1, f2), radius / (2 * eccentricity()));
  Line dir1 = perpendicular_through_point(points.first, Line(f1, f2));
  Line dir2 = perpendicular_through_point(points.second, Line(f1, f2));
  return {dir1, dir2};
}

class Circle : public Ellipse {
 public:
  Circle(Point p, double radius) : Ellipse(p, p, radius * 2) {}
  Circle(Point p1, Point p2) : Ellipse(p1, p1, 2 * dist(p1, p2)) {}
  double radius() const { return Ellipse::radius / 2; }
};

class Polygon : public Shape {
 protected:
  std::vector<Point> vertices;
 public:
  template<typename... Points>
  Polygon (Points... points) {
    vertices = {std::forward<Points>(points)...};
  }
  Polygon(std::vector<Point> vertices) : vertices(vertices) {}
  size_t verticesCount() const { return vertices.size(); }
  std::vector<Point> getVertices() const { return vertices; }
  bool isConvex() const;
  double area() const override;
  double perimeter() const override;
  void print() const {
    for (size_t i = 0; i < verticesCount(); ++i) {
      std:: cout << vertices[i].x << ' ' << vertices[i].y << '\n';
    }
  }
  bool isEqual(const Shape& other) const override {
    const Polygon* other_ptr = dynamic_cast<const Polygon*>(&other);
    if (other_ptr == nullptr) {
      return false;
    }
    if (verticesCount() != other_ptr->verticesCount()) {
      return false;
    }
    bool is_equal = false;
    auto v1 = getVertices();
    auto v2 = other_ptr->getVertices();
    for (int i = 0; i < static_cast<int>(v2.size()); ++i) {
      int cnt = 0;
      for (int k = 0; k < static_cast<int>(v2.size()); ++k) {
        if (v1[k] == v2[(i + k) % v2.size()]) {
          ++cnt;
        }
      }
      if (cnt == static_cast<int>(v1.size())) {
        is_equal = true;
      }
      cnt = 0;
      for (int k = 0; k < static_cast<int>(v2.size()); ++k) {
        if (v1[k] == v2[(i + v2.size() - k) % v2.size()]) {
          ++cnt;
        }
      }
      if (cnt == static_cast<int>(v1.size())) {
        is_equal = true;
      }
    }
    return is_equal;
  }

  bool isPermutationsSame(const std::vector<Point>& v1, const std::vector<Point>& v2, bool order) const {
    int size = v1.size();
    for (int k = 0; k < size; ++k) {
      bool is_equal = true;
      for (int i = 0; i < size; ++i) {
        Point next1 = v1[(i + 1) % size];
        Point next11 = v1[(i + 2) % size];
        Point next2 = order ? v2[(i + k + 1) % size] : v2[(-i + 2 * size + k - 1) % size];
        Point next22 = order ? v2[(i + k + 2) % size] : v2[(-i + 2 * size + k - 2) % size];
        int i2 = (order ? i + k : -i + size + k) % size;
        if (!equal(dist(v1[i], next1), dist(v2[i2], next2))) {
          is_equal = false;
        }
        if (!equal(angle(v1[i], next1, next1, next11), angle(v2[i2], next2, next2, next22))) {
          is_equal = false;
        }
      }
      if (is_equal) {
        return true;
      }
    }
    return false;
  }

  bool isCongruentTo(const Shape& other) const override {
    const Polygon* other_ptr = dynamic_cast<const Polygon*>(&other);
    if (other_ptr == nullptr) {
      return false;
    }
    if (verticesCount() != other_ptr->verticesCount()) {
      return false;
    }
    std::vector<Point> v1 = getVertices();
    std::vector<Point> v2 = other_ptr->getVertices();
    return (isPermutationsSame(v1, v2, true) || isPermutationsSame(v1, v2, false));
  }
  bool isSimilarTo(const Shape& other) const override {
    const Polygon* other_ptr = dynamic_cast<const Polygon*>(&other);
    if (other_ptr == nullptr) {
      return false;
    }
    Polygon copy(other_ptr->getVertices());
    copy.scale(Point(0, 0), perimeter() / other_ptr->perimeter());
    return isCongruentTo(copy);
  }
  bool containsPoint(const Point& p) const override {
    double angle_sum = 0;
    for (size_t i = 0; i < verticesCount(); ++i) {
      Point next = vertices[(i + 1) % verticesCount()];
      angle_sum += signed_angle(p, vertices[i], p, next);
    }
    return (equal(fabs(angle_sum), 2 * geom::PI));
  }
  void rotate(const Point& center, double angle) override {
    for (size_t i = 0; i < verticesCount(); ++i) {
      vertices[i] = rotation(center, vertices[i], angle);
    }
  }
  void reflect(const Point& center) override {
    for (size_t i = 0; i < verticesCount(); ++i) {
      vertices[i] = reflection(center, vertices[i]);
    }
  }
  void reflect(const Line& axis) override {
    for (size_t i = 0; i < verticesCount(); ++i) {
      vertices[i] = reflection_over_line(axis, vertices[i]);
    } 
  }
  void scale(const Point& center, double coefficient) override {
    for (size_t i = 0; i < verticesCount(); ++i) {
      vertices[i] = homothety(center, vertices[i], coefficient);
    }
  }
};

double Polygon::perimeter() const {
  double perimeter = 0;
  for (size_t i = 0; i < verticesCount(); ++i) {
    int next = (i + 1) % verticesCount();
    perimeter += dist(vertices[i], vertices[next]);
  }
  return perimeter;
}

double Polygon::area() const {
  double area = 0;
  Point o(0, 0);
  for (size_t i = 0; i < verticesCount(); ++i) {
    int next = (i + 1) % verticesCount();
    area += cross_product(o, vertices[i], o, vertices[next]);
  }
  return fabs(area) / 2;
}

bool Polygon::isConvex() const {
  int positive = 0, negative = 0;
  for (size_t i = 0; i < verticesCount(); ++i) {
    int next1 = (i + 1) % verticesCount();
    int next2 = (i + 2) % verticesCount();
    cross_product(vertices[i], vertices[next1], vertices[next1], vertices[next2]) > 0 ? ++positive : ++negative;
  }
  return (positive * negative == 0);
}

class Rectangle : public Polygon {
 public:
  Rectangle(Point p1, Point p2, double ratio) {
    ratio = std::min(ratio, 1 / ratio);
    double small_edge = dist(p1, p2) / sqrt(1 + 1 / (ratio * ratio));
    Point point_to_rotate = move_point_on_line(p1, Line(p1, p2), small_edge).second;
    vertices.push_back(p1);
    vertices.push_back(rotation(p1, point_to_rotate, acos(small_edge / dist(p1, p2)) * geom::PI_DEGREE / geom::PI));
    vertices.push_back(p2);
    vertices.push_back(reflection(middle(p1, p2), vertices[1]));
  }

  Point center() const { return middle(vertices[0], vertices[2]); }
  std::pair<Line, Line> diagonals() const { return {Line(vertices[0], vertices[2]),
                                                    Line(vertices[1], vertices[3])}; }
};

class Square : public Rectangle {
 public:
  Square(Point p1, Point p2) : Rectangle(p1, p2, 1) {}
  Circle circumscribedCircle() const { return Circle(center(), dist(vertices[0], center())); }
  Circle inscribedCircle() const { return Circle(center(), dist(vertices[0], vertices[1]) / 2); }
};

class Triangle : public Polygon {
 public:
  using Polygon::Polygon;
  std::vector<double> getEdges() const { return {dist(vertices[0], vertices[1]),
                                                 dist(vertices[1], vertices[2]),
                                                 dist(vertices[2], vertices[0])}; }
  Circle inscribedCircle() const;
  Circle circumscribedCircle() const { return Circle(circumcenter(),vertices[0]); }
  Point centroid() const { return Point(segment_division(vertices[0], middle(vertices[1], vertices[2]), -2)); }
  Point circumcenter() const;
  Point orthocenter() const { return segment_division(circumcenter(), centroid(), 1.5); }
  Line EulerLine() const { return Line(orthocenter(), centroid()); }
  Circle ninePointsCircle() const { return Circle(middle(orthocenter(), circumcenter()), middle(vertices[0], vertices[1])); }
};

Circle Triangle::inscribedCircle() const {
  auto edges = getEdges();
  Point i((edges[0] * vertices[2].x + edges[1] * vertices[0].x + edges[2] * vertices[1].x) / (edges[0] + edges[1] + edges[2]),
          (edges[0] * vertices[2].y + edges[1] * vertices[0].y + edges[2] * vertices[1].y) / (edges[0] + edges[1] + edges[2]));
  return Circle(i, distance_to_line(i, Line(vertices[0], vertices[1])));
}

Point Triangle::circumcenter() const {
  Line l1 = bisector_line(vertices[0], vertices[1]);

  Line l2 = bisector_line(vertices[1], vertices[2]);
  return intersection(l1, l2);
}