class Hyperbola {
    private:
        double h, k, a, b2;
        bool horizontal;
        std::shared_ptr<Point> bisection(double a, double b,
			Hyperbola * hyp, bool p1, bool p2);
		double solve(double x, bool p);
    public:
        Hyperbola(Point * M1, Point * M2, double TDOA);
        Hyperbola(double h, double k, double a, double b2, bool horizontal);
		std::vector<std::shared_ptr<Point>> intersect(Hyperbola * hyp);
};
