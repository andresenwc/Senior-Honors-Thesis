#include <algorithm>
#include <cmath>
#include <stdint.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <boost/multiprecision/cpp_int.hpp>
#include "Point.h"
#include "Hyperbola.h"
#include "ImpactFinder.h"
#include "config.h"

using namespace std;
using namespace boost::multiprecision;

/* ImpactFinder - Constructor to initialize all fields.
 *  @param: Point * M0 - Coordinate location of Mic 0.
 *  @param: Point * M1 - Coordinate location of Mic 1.
 *  @param: Point * M2 - Coordinate location of Mic 2.
 *  @param: Point * M3 - Coordinate location of Mic 3.
 *  @param: string filename - Name of file to be read.
 *  @param: bool verbose - Toggles verbose printing.
 *  @param: bool debug - Toggles debug printing.
 */
ImpactFinder::ImpactFinder(Point * M0, Point * M1, Point * M2, Point * M3,
    string filename, bool verbose, bool debug) {
    this->M0 = M0;
    this->M1 = M1;
    this->M2 = M2;
    this->M3 = M3;
    this->filename = filename;
    this->verbose = verbose;
    this->debug = debug;
    fout.open("out.txt");
}

/* dataProcessing - Processes input data from file in order to determine sound locations. */
void ImpactFinder::dataProcessing() {
    int count = 0;
    ifstream file(filename.c_str());
    string line;
    
    if (verbose) {
        cout << "====================" << endl;
    } 

    while (getline(file, line)) {
        
        // Contains all 1024 samples for each microphone
        uint1024_t m0, m1, m2, m3;
        // x and y coordinates of sound
        // State of microphone and number of cycles it remained in that state
        uint64_t x, y, state, cycles;
        // Various strings for holding information
        string fluff, coords, data;
        
        // Check for start of data region
        if (!(line.find("Ready") != string::npos)) {
            if (debug) 
                cerr << "Error: Start of data region not found. Continuing." << endl;
            continue;
        }

        // Create stringstream for parsing lines
        istringstream iss(line);

        if (debug) {
            cout << "Location line: \"" << line << "\"" << endl;
        }

        // Get grid location: "Ready (x,y)"
        if ((iss >> fluff >> coords)) {
            if (!(sscanf(coords.c_str(), "(%lu,%lu)", &x, &y)))
                continue;
        }
        else {
            cerr << "Error: Failed to read grid location. Continuing." << endl;
            continue;
        }
/*
        if (x < 15 || x > 17 || y < 15 || y > 17)
            continue;
*/
/*
        if (x != 1 || y != 20)
            continue;
*/
        // Skip fluff: "tstart <time>"
        getline(file, line);
        iss.clear();
        iss.str(line);

        if (debug) {
            cout << "Fluff line: \"" << line << "\"" << endl;
        }

        // Read the data from each port:
        //  "init:    0    trace:    (state,cycles)    ...
        getline(file, line);
        iss.clear();
        iss.str(line);

        if (debug) {
            cout << "Data line: \"" << line << "\"" << endl;
        }

        // Skip fluff: "init:    0    trace:"
        iss >> fluff >> fluff >> fluff;
        // Read and use each piece of data: "(state,cycles)"
        int numCycles = 0;
        while (iss >> data) {
            sscanf(data.c_str(), "(%lu,%lu)", &state, &cycles);
            // Append data for each mic
            for (uint64_t i = 0; i < cycles; i++) {
                m0 = (m0 << 1) + (state & 1);
                m1 = (m1 << 1) + ((state & 2) >> 1);
                m2 = (m2 << 1) + ((state & 4) >> 2);
                m3 = (m3 << 1) + ((state & 8) >> 3);
                numCycles++;
            }
        }

        if (numCycles != 1024) {
            cerr << "Incorrect number of data points: <" << numCycles << ">\nContinuing." << endl;
        }

        // At this point, each microphone variable should contain
        // the data for that point. Now we can attempt to determine
        // TDOAs.
        
        double TDOA[4];
        int closestMic = findTDOAs(TDOA, m0, m1, m2, m3);
        
        if (verbose) {
            cout << "====================" << endl;
            cout << "(" << x << "," << y << ")" << endl;
            cout << "TDOA[0]: " << TDOA[0] << endl;
            cout << "TDOA[1]: " << TDOA[1] << endl;
            cout << "TDOA[2]: " << TDOA[2] << endl;
            cout << "TDOA[3]: " << TDOA[3] << endl;
        }
        fout << x << " " << y << " ";
        
        // Create and intersect hyperbolae using TDOAs.
        findPoint(TDOA, closestMic);
        
        count++;
        //if (count == 1)
        //    exit(0);
    
    }
}

/* findTDOAs - Determines TDOAs using one of three defined algorithms, placing
 * this information into passed in TDOA array. Also determines the mic closest
 * to the sound.
 *  @param: double (& TDOA)[4] - TDOA array, passed by reference.
 *  @param: uint1024_t m0 - Raw signal data of Mic 0.
 *  @param: uint1024_t m1 - Raw signal data of Mic 1.
 *  @param: uint1024_t m2 - Raw signal data of Mic 2.
 *  @param: uint1024_t m3 - Raw signal data of Mic 3.
 *  @return: Index of the closest microphone.
 */
int ImpactFinder::findTDOAs(double (& TDOA)[4], uint1024_t m0,
    uint1024_t m1, uint1024_t m2, uint1024_t m3) {
    
    int closestMic = 0;
    
    // USING FIRST SUSTAINED HIT OR FIRST HIT ALGORITHM
    #if defined(FIRST_SUSTAINED_HIT) || defined(FIRST_HIT)
        // Create an array to store signal "true start times."
        double start[4];
        // Perform cross correlation between the closest mic
        // and each other mic to determine TDOA.
        #ifdef FIRST_SUSTAINED_HIT
            start[0] = firstSustainedHit(m0);
            start[1] = firstSustainedHit(m1);
            start[2] = firstSustainedHit(m2);
            start[3] = firstSustainedHit(m3);
        #else
            start[0] = firstHit(m0);
            start[1] = firstHit(m1);
            start[2] = firstHit(m2);
            start[3] = firstHit(m3);
        #endif
        // We now have start times. The least start time is
        // considered the closest mic, and thus the mic to
        // which all TDOA will be built against.
        //
        // Find closest mic start time and index.
        double leastTime = start[0];
        for (int i = 1; i < 4; i++) {
            if (start[i] < leastTime) {
                closestMic = i;
                leastTime = start[i];
            }
        }
        // Determine TDOAs
        for (int i = 0; i < 4; i++) {
            if (closestMic == i)
                TDOA[i] = 0.0;
            else
                TDOA[i] = start[i] - start[closestMic];
        }
    #endif
    
    // USING CROSS CORRELATION ALGORITHM
    #ifdef CROSS_CORRELATION
        // Mic array to ease calculations
        uint1024_t mics[4] = {m0, m1, m2, m3};
        // Determine closest mic
        for (int i = 0; i < 4; i++) {
            if ((mics[i] >> 1023) > 0)
                closestMic = i;
        }
        TDOA[0] = crossCorrelation(mics[closestMic], mics[0]);
        TDOA[1] = crossCorrelation(mics[closestMic], mics[1]);
        TDOA[2] = crossCorrelation(mics[closestMic], mics[2]);
        TDOA[3] = crossCorrelation(mics[closestMic], mics[3]);
    #endif
    
    return closestMic;
}

/* firstSustainedHit - Determines "true start time" as first sustained
 * hit, and therefore relative to the speed of sound in air.
 *  @param: uint1024_t m - Raw signal data of the microphone.
 *  @return: True start time of mic m.
 */
double ImpactFinder::firstSustainedHit(uint1024_t m) {
    uint1024_t one = 1;
    for (int i = 1023; i > 256; i--) {
        int count = 0;
        for (int j = 0; j < 80; j++) {
            if ((m & (one << (i-j))) > 0) {
                count++;
                if (count >= 80) {
                    return (double)(1023-i) * 2.0 *.000001;
                }
            }
            else
                break;
        }
    }
    return 0;
}

/* firstHit - Determines "true start time" as time of first "on" data
 * point, and therefore relative to the speed of sound in the apparatus.
 *  @param: uint1024_t m - Raw signal data of the microphone.
 *  @return: True start time of mic m.
 */
double ImpactFinder::firstHit(uint1024_t m) {
    uint1024_t one = 1;
    for (int i = 1023; i > 0; i--) {
        if ((m & (one<<i)) > 0) {
            return (double)(1023-i) * 2.0 * .000001;
        }
    }
    return 0;
}

/* crossCorrelation - Determines TDOA between two microphone signals using
 * cross correlation.
 *  @param: uint1024_t c - Raw signal data of the "closest" microphone.
 *  @param: uint1024_t m - Raw signal data of the microphone to be compared
 *   against.
 *  @return: TDOA between mic c and mic m.
 */
double ImpactFinder::crossCorrelation(uint1024_t c, uint1024_t m) {
    // Would generally only apply to the comparison of the closest mic
    // with itself.
    if (c == m) {
        return 0;
    }
    else {
        // bestFit holds the shift at which the highest sum was achieved.
        // bestSum holds this highest sum.
        int bestFit = -1, bestSum = 0;
        // Perform sliding dot product between c and each "iteration" of
        // m. We will only shift halfway.
        for (int i = 0; i < 511; i++) {
            // Temporary sum variable.
            int tempSum = 0;
            // Shift m over
            uint1024_t shiftm = m << i;
            // Perform a bitwise and between c and m.
            uint1024_t dotprod = c & shiftm;
            // Count the number of ones in dotprod. This is the sum. We will
            // only look at the first half of dotprod.
            uint1024_t shift = 1;
            for (int j = 0; j < 1023; j++) {
                if ((shift & dotprod) > 0)
                    tempSum++;
                shift = shift << 1;
            }
            // If tempSum is higher, we have a new
            // best fit.
            if (tempSum >= bestSum) {
                bestSum = tempSum;
                bestFit = i;
            }
        }
        // bestFit contains the number of "slides" making up the best sum.
        // Multiply this number of slides by the time elapsed between each
        // data point (2 microseconds). This is the TDOA in microseconds.
        // Then, convert to seconds.
        return bestFit * 2 * .000001;
    }
}

/* findPoint - Determines point of sound.
 *  @param: double TDOA[] - Contains TDOA information.
 *  @param: int closest - Contains index of the closest mic to the sound.
 */
void ImpactFinder::findPoint(double TDOA[], int closest) {
    
    double TDOA_M0M1, TDOA_M0M2, TDOA_M1M3, TDOA_M2M3;
    
    // Pull out TDOA information between each microphone.
    if (closest == 0) {
        TDOA_M0M1 = TDOA[1];
        TDOA_M0M2 = TDOA[2];
        TDOA_M1M3 = abs(TDOA[1]-TDOA[3]);
        TDOA_M2M3 = abs(TDOA[2]-TDOA[3]);
    }
    else if (closest == 1) {
        TDOA_M0M1 = TDOA[0];
        TDOA_M0M2 = abs(TDOA[0]-TDOA[2]);
        TDOA_M1M3 = TDOA[3];
        TDOA_M2M3 = abs(TDOA[2]-TDOA[3]);
    }
    else if (closest == 2) {
        TDOA_M0M1 = abs(TDOA[0]-TDOA[1]);
        TDOA_M0M2 = TDOA[0];
        TDOA_M1M3 = abs(TDOA[1]-TDOA[3]);
        TDOA_M2M3 = TDOA[3];
    }
    else if (closest == 3) {
        TDOA_M0M1 = abs(TDOA[0]-TDOA[1]);
        TDOA_M0M2 = abs(TDOA[0]-TDOA[2]);
        TDOA_M1M3 = TDOA[1];
        TDOA_M2M3 = TDOA[2];
    }
   
    // Build hyperbolae using TDOA 
    Hyperbola * M0M1 = new Hyperbola(M0, M1, TDOA_M0M1);
    Hyperbola * M0M2 = new Hyperbola(M0, M2, TDOA_M0M2);
    Hyperbola * M1M3 = new Hyperbola(M1, M3, TDOA_M1M3);
    Hyperbola * M2M3 = new Hyperbola(M2, M3, TDOA_M2M3);
    
    // Perform all intersection combinations
    
    std::vector<std::shared_ptr<Point>> int1, int2, int3, int4, int5, int6;

    // Intersections between M0M1 and M0M2
    int1 = M0M1->intersect(M0M2);

    // Intersections between M0M1 and M1M3
    int2 = M0M1->intersect(M1M3);
    
    // Intersections between M0M1 and M2M3
    int3 = M0M1->intersect(M2M3);
    
    // Intersections between M0M2 and M1M3
    int4 = M0M2->intersect(M1M3);
    
    // Intersections between M0M2 and M2M3
    int5 = M0M2->intersect(M2M3);
    
    // Intersections between M1M3 and M2M3
    int6 = M1M3->intersect(M2M3);
    
    // Combine all intersection vectors
    std::vector<std::shared_ptr<Point>> intersections;
    intersections.reserve(int1.size() + int2.size() + int3.size() + int4.size() + int5.size() + int6.size());
    intersections.insert(intersections.end(), int1.begin(), int1.end());
    intersections.insert(intersections.end(), int2.begin(), int2.end());
    intersections.insert(intersections.end(), int3.begin(), int3.end());
    intersections.insert(intersections.end(), int4.begin(), int4.end());
    intersections.insert(intersections.end(), int5.begin(), int5.end());
    intersections.insert(intersections.end(), int6.begin(), int6.end());
    
    Point * mics[4] = {M0, M1, M2, M3};
    double cxval = mics[closest]->x;
    double cyval = mics[closest]->y;
    
    // Center of the microphone array
    Point * center = new Point(X_DIM/2, Y_DIM/2);
    // Center of the microphone's quadrant
    // (midpoint between rectangle center of microphone and the microphone itself)
    Point * quadcenter = new Point((center->x + cxval) / 2, (center->y + cyval) / 2);
    
    double minx = quadcenter->x - X_DIM/4 - 1;
    double maxx = quadcenter->x + X_DIM/4 + 1;
    double miny = quadcenter->y - Y_DIM/4 - 1;
    double maxy = quadcenter->y + Y_DIM/4 + 1;
    
    // Marks points outside the array as invalid.
    for (uint i = 0; i < intersections.size(); i++) {
        if (intersections[i]) {
            double xval = intersections[i]->x;
            double yval = intersections[i]->y;
            if (xval < minx || xval > maxx || yval < miny || yval > maxy) {
                intersections[i]->x = -1;
                intersections[i]->y = -1;
            }
        }
    }
    
    // Removes nullptrs and invalidated (-1,-1) points.
    intersections.erase(std::remove_if(intersections.begin(), intersections.end(),
        [](std::shared_ptr<Point> p){
            if (p == nullptr)
                return true;
            else
                return (p->x == -1 && p->y == -1);
        }),
    intersections.end());
    
    #ifdef CLUSTER
    // Use DBSCAN to get optimal cluster
    intersections = DBSCAN(intersections, 4, 2);
    #elif defined STAT
    // Use standard deviation (SD) to exclude outliers
    intersections = SD(intersections);
    #endif
    
    // Final removal pass.
    intersections.erase(std::remove_if(intersections.begin(), intersections.end(),
        [](std::shared_ptr<Point> p){
            if (p == nullptr)
                return true;
            else
                return (p->x == -1 && p->y == -1);
        }),
    intersections.end());
    
    // At this point the remaining points should all be viable
    // candidates. Simply average and round them.
    double finalx = 0.0;
    double finaly = 0.0;
    for (uint i = 0; i < intersections.size(); i++) {
        finalx += intersections[i]->x;
        finaly += intersections[i]->y;
    }
    finalx = round(finalx/intersections.size());
    finaly = round(finaly/intersections.size());
    
    if (verbose) {
        cout << "(" << finalx << "," << finaly << ")" << endl;
        cout << "====================" << endl;
    }
    
    fout << finalx << " " << finaly << "\n";
    fout.flush();
}

/* DBSCAN - (Density-Based Spatial Clustering of Applications with Noise)
 * A data clustering algorithm to weed out problem points.
 *  @param: std::vector<std::shared_ptr<Point>> intersections - vector of
 *   intersection points
 *  @param: double maxdist - minimum distance between points
 *  @param: int minpts - minimum number of points to constitute a cluster
 *  @return: Cluster with enough points.
 */
std::vector<std::shared_ptr<Point>> ImpactFinder::DBSCAN(
    std::vector<std::shared_ptr<Point>> intersections,
    double maxdist, uint minpts) {
    
    // Convert intersections to a vector of pairs of <Point, int>.
    // This allows us to assign each point a cluster (int) while
    // maintaining them all within the same overarching data structure.
    std::vector<std::pair<std::shared_ptr<Point>, int>> points;
    for (uint i = 0; i < intersections.size(); i++) {
        std::pair<std::shared_ptr<Point>, int> p(intersections[i], 0);
        points.push_back(p);
    };
    
    // Cluster labels
    int UNDEFINED = 0;
    int NOISE = -1;
    int C = 0;
    // Determine clusters
    for (uint i = 0; i < points.size(); i++) {
        std::pair<std::shared_ptr<Point>, int> * P = &(points[i]);
        // P is already part of a cluster.
        if (P->second != UNDEFINED)
            continue;
        // Gather points neighboring P
        std::vector<std::pair<std::shared_ptr<Point>, int>*> Pneighbors;
        Pneighbors = inRange(points, maxdist, P);
        // Too small to become a cluster, is instead noise.
        if (Pneighbors.size() < minpts)
            P->second = NOISE;
        // Big enough to be a cluster
        else {
            // Cluster label
            C++;
            // Make P part of the new cluster.
            P->second = C;
            // Expand cluster.
            uint j = 0;
            uint Psize = Pneighbors.size();
            while (j < Psize) {
                std::pair<std::shared_ptr<Point>, int> * N = Pneighbors[j];
                // N has already been processed.
                if (N->second != UNDEFINED && N->second != NOISE) {
                    j++;
                    Psize = Pneighbors.size();
                    continue;
                }
                // Add N to this cluster, as well as N's neighbors.
                else {
                    N->second = C;
                    // Gather points neighboring N
                    std::vector<std::pair<std::shared_ptr<Point>, int>*> Nneighbors;
                    Nneighbors = inRange(points, maxdist, N);
                    // Extend cluster if dense enough around N
                    if (Nneighbors.size() >= minpts) {
                        for (uint i = 0; i < Pneighbors.size(); i++) {
                        }
                        for (uint i = 0; i < Nneighbors.size(); i++) {
                            if (std::find(Pneighbors.begin(), Pneighbors.end(), Nneighbors[i]) == Pneighbors.end()) {
                                Pneighbors.push_back(Nneighbors[i]);
                            }
                        }
                        //Pneighbors.insert(Pneighbors.end(), Nneighbors.begin(), Nneighbors.end());
                    }
                }
                j++;
                Psize = Pneighbors.size();
            }
        }
    }
    
    return intersections;
}

/* inRange - Helper for DBSCAN. Returns all Points in a range of maxdist
 * of Point P.
 *  @param: std::vector<std::shared_ptr<Point>> intersections - vector of
 *   intersection points
 *  @param: double maxdist - minimum distance between points
 *  @param: std::shared_ptr<Point> P - point to compare against
 *  @return: All neighboring points.
 */
std::vector<std::pair<std::shared_ptr<Point>, int>*> ImpactFinder::inRange(
    std::vector<std::pair<std::shared_ptr<Point>, int>> & points,
    double maxdist, std::pair<std::shared_ptr<Point>, int> * P) {
    
    std::vector<std::pair<std::shared_ptr<Point>, int>*> neighbors;
    for (uint i = 0; i < points.size(); i++) {
        double x = points[i].first->x;
        double y = points[i].first->y;
        double dist = sqrt(pow(x - P->first->x, 2) + pow(y - P->first->y, 2));
        if (dist <= maxdist)
            neighbors.push_back(&(points[i]));
    }
    return neighbors;
}

/* SD - Uses standard deviation of mean distance of average point to
 * exclude outlier points.
 *  @param: std::vector<std::shared_ptr<Point>> intersections - vector of
 *   intersection points
 */
std::vector<std::shared_ptr<Point>> ImpactFinder::SD(
    std::vector<std::shared_ptr<Point>> intersections) {
        
    // Mean x and y values to determine mean point.
    double xsum = 0.0;
    double ysum = 0.0;
    for (uint i = 0; i < intersections.size(); i++) {
        xsum += intersections[i]->x;
        ysum += intersections[i]->y;
    }
    double xav = xsum/intersections.size();
    double yav = ysum/intersections.size();
    
    // Calculate mean distance from mean point.
    double avdist = 0.0;
    for (uint i = 0; i < intersections.size(); i++) {
        double xval = intersections[i]->x;
        double yval = intersections[i]->y;
        double dist = sqrt(pow(xval-xav, 2) + pow(yval-yav, 2));
        avdist += dist;
    }
    avdist = avdist/intersections.size();
    
    // Calculate standard deviation of mean distance.
    double sd = 0.0;
    for (uint i = 0; i < intersections.size(); i++) {
        double xval = intersections[i]->x;
        double yval = intersections[i]->y;
        double dist = sqrt(pow(xval-xav, 2) + pow(yval-yav, 2));
        sd += pow(dist-avdist, 2);
    }
    sd = sqrt(sd / intersections.size());
    
    // Now that we have the standard deviation, we will invalidate
    // all points greater than two standard deviations from the mean
    // distance. These points should be outliers.
    for (uint i = 0; i < intersections.size(); i++) {
        double xval = intersections[i]->x;
        double yval = intersections[i]->y;
        double dist = sqrt(pow(xval-xav, 2) + pow(yval-yav, 2));
        if (dist > (sd * 2)) {
            intersections[i]->x = -1;
            intersections[i]->y = -1;
        }
    }
    
    return intersections;
}
