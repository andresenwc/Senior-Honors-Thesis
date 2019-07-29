#include <cmath>
class Point {
	public:
		double x;
		double y;
		Point(double x, double y) { this->x = x; this->y = y; }
        double distance(Point * P) { return sqrt(pow(P->x - x, 2) + pow(P->y - y, 2)); }
};
