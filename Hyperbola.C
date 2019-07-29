#include <vector>
#include <cmath>
#include <memory>
#include <iostream>
#include "Point.h"
#include "Hyperbola.h"
#include "config.h"

using namespace std;

/* Hyperbola - Constructs a hyperbola given mic locations and TDOA.
 *  @param: Point * M1 - Coordinate location of Mic 1. (First focus)
 *  @param: Point * M2 - Coordinate location of Mic 2. (Second focus)
 *  @param: double TDOA - TDOA between M1 and M2.
 */
Hyperbola::Hyperbola(Point * M1, Point * M2, double TDOA) {

    double x1 = M1->x;
    double y1 = M1->y;
    double x2 = M2->x;
    double y2 = M2->y;

    if (y1 == y2)
        this->horizontal = true;
    else if (x1 == x2)
        this->horizontal = false;
    else {
        cerr << "Hyperbola not axis aligned." << endl;
        exit(0);
    }

    // (h,k) is the center of the hyperbola
    this->h = (x1+x2)/2;
    this->k = (y1+y2)/2;

    // time delta (TDOA) is used to build the hyperbola
    double delta = TDOA;
    // Find distance covered by delta in centimeters
    double d_delta = delta * SPEEDOFSOUND * 100;
    
    this->a = d_delta/2;

    // Distance between M1 and M2
    double d = sqrt(pow(x2-x1,2) + pow(y2-y1,2));

    this->b2 = pow(d/2,2) - pow(a,2);
}

/* Hyperbola - Constructs a hyperbola given all information.
 *  @param: double h - x coordinate of midpoint.
 *  @param: double k - y coordinate of midpoint.
 *  @param: double a - Half the constant difference. (shortest distance from
 *   hyperbola's center to a vertex)
 *  @param: double b2 - The distance, perp to the vertex, from the vertex to the
 *   hyperbola's asymptotic line, squared.
 *  @param: bool horizontal - Type of transverse axis. (horizontal or vertical)
 */
Hyperbola::Hyperbola(double h, double k, double a, double b2,
                     bool horizontal) {
    this->h = h;
    this->k = k;
    this->a = a;
    this->b2 = b2;
    this->horizontal = horizontal;
}

/* intersect - Intersects two hyperbolae using bisection, returning their
 * intersection point.
 *  @param: Hyperbola * hyp - Second hyperbola. (the other being this->)
 *  @return: All found intersection points.
 */
std::vector<std::shared_ptr<Point>> Hyperbola::intersect(Hyperbola * hyp) {
    
    std::vector<std::shared_ptr<Point>> intersections;

    // Vertical with vertical
    if (!this->horizontal && !hyp->horizontal) {
        // Intersecting two vertical hyperbolae may generate up to one
        // intersection point.
        
        std::shared_ptr<Point> p;
        bool b1 = true, b2 = true;
        // Alternate b1
        for (int i = 0; i < 2; i++) {
            // Alternate b2
            for (int j = 0; j < 2; j++) {
                p = this->bisection(-1, X_DIM+1, hyp, b1, b2);
                intersections.push_back(p);
                b2 = !b2;
            }
            b1 = !b1;
        }
    }
    // Vertical with horizontal
    else if (!this->horizontal && hyp->horizontal) {
        // Intersecting a vertical hyperbola section with a horizontal
        // hyperbola section may generate up to two intersection points
        
        // Bounds are necessary when intersecting horizontal hyperbolae
        // because the space between the two branches is undefined.

        // upperbound
        double u = hyp->h - hyp->a;
        // lowerbound
        double l = hyp->h + hyp->a;

        std::shared_ptr<Point> p;
        bool b1 = true, b2 = true;
        // Alternate b1
        for (int i = 0; i < 2; i++) {
            // Alternate b2
            for (int j = 0; j < 2; j++) {
                p = this->bisection(-1, u, hyp, b1, b2);
                intersections.push_back(p);
                p = this->bisection(l, X_DIM+1, hyp, b1, b2);
                intersections.push_back(p);
                b2 = !b2;
            }
            b1 = !b1;
        }
    }
    // Horizontal with vertical
    else if (this->horizontal && !hyp->horizontal) {
        // Intersecting a horizontal hyperbola section with a vertical
        // hyperbola section may generate up to two intersection points

        // Bounds are necessary when intersecting horizontal hyperbolae
        // because the space between the two branches is undefined.       

        // upperbound
        double u = this->h - this->a;
        // lowerbound
        double l = this->h + this->a;

        std::shared_ptr<Point> p;
        bool b1 = true, b2 = true;
        // Alternate b1
        for (int i = 0; i < 2; i++) {
            // Alternate b2
            for (int j = 0; j < 2; j++) {
                p = this->bisection(-1, u, hyp, b1, b2);
                intersections.push_back(p);
                p = this->bisection(l, X_DIM+1, hyp, b1, b2);
                intersections.push_back(p);
                b2 = !b2;
            }
            b1 = !b1;
        }
    }
    // Horizontal with horizontal
    else if (this->horizontal && hyp->horizontal) {
        // Intersecting two horizontal hyperbolae may generate up to
        // two intersection points.
        
        // Bounds are necessary when intersecting horizontal hyperbolae
        // because the space between the two branches is undefined.

        // upper bound
        double u;
        // lower bound
        double l;

        if (this->a >= hyp->a) {
            // upperbound
            u = this->h - this->a;
            // lowerbound
            l = this->h + this->a;
        }
        else {
            // upperbound
            u = hyp->h - hyp->a;
            // lowerbound
            l = hyp->h + hyp->a;
        }

        std::shared_ptr<Point> p;
        bool b1 = true, b2 = true;
        // Alternate b1
        for (int i = 0; i < 2; i++) {
            // Alternate b2
            for (int j = 0; j < 2; j++) {
                p = this->bisection(-1, u, hyp, b1, b2);
                intersections.push_back(p);
                p = this->bisection(l, 32, hyp, b1, b2);
                intersections.push_back(p);
                b2 = !b2;
            }
            b1 = !b1;
        }
    }
    
    return intersections;
}

/* Performs bisection method to find intersection point.
 *  @param: double a - left bound
 *  @param: double b - right bound
 *  @param: Hyperbola * hyp - second hyperbola
 *  @param: bool p1 - select positive or negative solution of
 *      square root for this
 *  @param: bool p2 - select positive or negative solution of
 *      square root for hyp
 *  @return: Intersection point (if there is one) or nullptr (if there is
 *   not).
 */
std::shared_ptr<Point> Hyperbola::bisection(double a, double b,
    Hyperbola * hyp, bool p1, bool p2) {
    // The "function" that we are applying bisection to is
    // <this-hyp>. By substracting the two functions, bisection
    // should find the point at which the two functions are equal,
    // i.e. where they intersect.
    double c, solvea, solveb, solvec;
    solvea = this->solve(a,p1) - hyp->solve(a,p2);
    solveb = this->solve(b,p1) - hyp->solve(b,p2);
    if ((solvea * solveb >= 0) || std::isnan(solvea*solveb)) {
        return nullptr;
    }
    while ((b-a) > EPSILON) {
        // Calculate midpoint
        c = (a+b)/2;
        solvec = this->solve(c,p1) - hyp->solve(c,p2);
        solvea = this->solve(a,p1) - hyp->solve(a,p2);
        if (solvec == 0.0) 
            break;
        else if (solvec * solvea < 0)
            b = c;
        else
            a = c;
    }
    return std::shared_ptr<Point>(new Point(c, this->solve(c,p1)));
}

/* solve - Given x value of hyperbola, solves for y value.
 * @param: double x - x value to use.
 * @param: bool p - Designates whether to return positive or negative root.
 * @return: y value.
 */
double Hyperbola::solve(double x, bool p) {
    
    double sq;
    double xh2 = pow(x-h,2);
    double a2 = pow(a,2);
    
    if (horizontal)
        sq = sqrt(-1*(1-(xh2/a2))*b2);
    else
        sq = sqrt((1+(xh2/b2))*a2);
    
    if (p)
        return sq+k;
    else
        return (-1*sq)+k;
}
