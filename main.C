#include <boost/multiprecision/cpp_int.hpp>
#include <fstream>
#include "Point.h"
#include "ImpactFinder.h"
#include "config.h"

/* main - Initializes ImpactFinder object and calls dataProcessing(). */
int main(int argc, char ** argv) {
    ImpactFinder * impactFinder = new ImpactFinder(new Point(M0X, M0Y),
        new Point(M1X, M1Y), new Point(M2X, M2Y), new Point(M3X, M3Y),
        FILENAME, false, false);
    impactFinder->dataProcessing();
}