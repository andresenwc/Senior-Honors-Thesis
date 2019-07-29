using namespace std;
using namespace boost::multiprecision;

class ImpactFinder {
    private:
        Point * M0;
        Point * M1;
        Point * M2;
        Point * M3;
        string filename;
        bool verbose;
        bool debug;
        
        ofstream fout;
        
        int findTDOAs(double (& TDOA)[4], uint1024_t m0, uint1024_t m1,
            uint1024_t m2, uint1024_t m3);
        double firstSustainedHit(uint1024_t m);
        double firstHit(uint1024_t m);
        double crossCorrelation(uint1024_t c, uint1024_t m);

        void findPoint(double TDOA[], int closest);
        std::vector<std::shared_ptr<Point>> DBSCAN(
            std::vector<std::shared_ptr<Point>> intersections,
            double mindist, uint minpts);
        std::vector<std::pair<std::shared_ptr<Point>, int>*> inRange(
            std::vector<std::pair<std::shared_ptr<Point>, int>> & points,
            double maxdist, std::pair<std::shared_ptr<Point>, int> * P);
        std::vector<std::shared_ptr<Point>> SD(
            std::vector<std::shared_ptr<Point>> intersections);
    
    public:
        ImpactFinder(Point * M0, Point * M1, Point * M2, Point * M3,
            string filename, bool verbose, bool debug);
        void dataProcessing();
};

