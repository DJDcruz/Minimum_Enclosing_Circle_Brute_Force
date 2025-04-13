#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include <chrono>
#include <limits>
#include <future>
#include <thread>

using namespace std::chrono;
using namespace std;

struct Point {
    double x, y;
};

struct Circle {
    Point center;
    double radius;
    Circle(Point c, double r) : center(c), radius(r) {}
    Circle() : center({ 0, 0 }), radius(0) {}
};

double distance(const Point& p1, const Point& p2) {
    return std::sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y));
}

Circle computeCircle(const Point& p1, const Point& p2) {
    Point center = { (p1.x + p2.x) / 2, (p1.y + p2.y) / 2 };
    double rSquared = (pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2)) / 2.0;
    double r = sqrt(rSquared); 
    return Circle(center, r);
}

Circle computeCircleFromThreePoints(const Point& p1, const Point& p2, const Point& p3) {
    double A = p1.x * (p2.y - p3.y) + p2.x * (p3.y - p1.y) + p3.x * (p1.y - p2.y);
    if (A == 0) return Circle();

    double Bx = (p1.x * p1.x + p1.y * p1.y) * (p3.y - p2.y) +
        (p2.x * p2.x + p2.y * p2.y) * (p1.y - p3.y) +
        (p3.x * p3.x + p3.y * p3.y) * (p2.y - p1.y);
    double By = (p1.x * p1.x + p1.y * p1.y) * (p2.x - p3.x) +
        (p2.x * p2.x + p2.y * p2.y) * (p3.x - p1.x) +
        (p3.x * p3.x + p3.y * p3.y) * (p1.x - p2.x);
    Point center = { -Bx / (2 * A), -By / (2 * A) };
    double r = distance(center, p1);
    return Circle(center, r);
}

bool checkPointsForPairs(const Circle& c, const std::vector<Point>& points) {
    for (const auto& p : points) {
        if (distance(c.center, p) > c.radius) {
            return false;
        }
    }
    return true;
}

bool checkPointsForTrip(const Circle& c, const std::vector<Point>& points) {
    for (const auto& p : points) {
        if (distance(c.center, p) > c.radius) {
            return false;
        }
    }
    return true;
}

Circle findSmallestCircleSubset(const std::vector<Point>& points, int start, int end) {
    Circle bestCircle({ 0, 0 }, std::numeric_limits<double>::max());

    for (int i = start; i < end; i++) {
        for (int j = i + 1; j < points.size(); j++) {
            Circle c = computeCircle(points[i], points[j]);
            if (checkPointsForPairs(c, points)) {
                if (c.radius < bestCircle.radius) {
                    bestCircle = c;
                }
            }
        }
    }

    return bestCircle;
}

Circle findSmallestCircleTripletSubset(const std::vector<Point>& points, int start, int end) {
    Circle minCircle({ 0, 0 }, std::numeric_limits<double>::max());

    for (int i = start; i < end; i++) {
        for (int j = i + 1; j < points.size(); j++) {
            for (int k = j + 1; k < points.size(); k++) {
                Circle c = computeCircleFromThreePoints(points[i], points[j], points[k]);
                if (c.radius < minCircle.radius && checkPointsForTrip(c, points)) {
                    minCircle = c;
                }
            }
        }
    }

    return minCircle;
}

Circle findSmallestCircleTriplet(const std::vector<Point>& points) {
    int numThreads = 16;
    int pointsPerThread = points.size() / numThreads;
    std::vector<std::future<Circle>> futures;

    for (int i = 0; i < numThreads; i++) {
        int start = i * pointsPerThread;
        int end = (i == numThreads - 1) ? points.size() : (i + 1) * pointsPerThread;
        futures.push_back(std::async(std::launch::async, findSmallestCircleTripletSubset, std::ref(points), start, end));
    }

    Circle minCircle({ 0, 0 }, std::numeric_limits<double>::max());
    for (auto& future : futures) {
        Circle c = future.get();
        if (c.radius < minCircle.radius) {
            minCircle = c;
        }
    }

    return minCircle;
}

Circle findSmallestCircle(const std::vector<Point>& points) {
    int numThreads = std::thread::hardware_concurrency();
    int pointsPerThread = points.size() / numThreads;
    std::vector<std::future<Circle>> futures;

    for (int i = 0; i < numThreads; i++) {
        int start = i * pointsPerThread;
        int end = (i == numThreads - 1) ? points.size() : (i + 1) * pointsPerThread;
        futures.push_back(std::async(std::launch::async, findSmallestCircleSubset, std::ref(points), start, end));
    }

    Circle bestCircle({ 0, 0 }, std::numeric_limits<double>::max());
    for (auto& future : futures) {
        Circle c = future.get();
        if (c.radius < bestCircle.radius) {
            bestCircle = c;
        }
    }

    if (bestCircle.radius == std::numeric_limits<double>::max()) {
        bestCircle = findSmallestCircleTriplet(points);
    }

    return bestCircle;
}

std::vector<Point> readPointsFromFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        exit(1);
    }

    std::vector<Point> points;
    int numPoints;
    file >> numPoints;
    points.reserve(numPoints);

    for (int i = 0; i < numPoints; ++i) {
        double x, y;
        file >> x >> y;
        points.push_back({ x, y });
    }

    file.close();
    return points;
}

void performAnalysis(const std::string& file) {
    std::string filename = "C:/Users/Pixel/OneDrive/Desktop/Min_enclosing_circ_with_task_parallelization/datasets/circle_" + file + ".txt"; //Replace with your path here
    std::vector<Point> points = readPointsFromFile(filename);

    auto start = high_resolution_clock::now();
    Circle resultPair = findSmallestCircle(points);
    auto stopPair = high_resolution_clock::now();
    auto durationPair = duration_cast<milliseconds>(stopPair - start);

    start = high_resolution_clock::now();
    Circle resultTriplet = findSmallestCircleTriplet(points);
    auto stopTriplet = high_resolution_clock::now();
    auto durationTriplet = duration_cast<milliseconds>(stopTriplet - start);

    std::cout << file + " Points Dataset" << std::endl;
    std::cout << "Smallest enclosing circle using pair of points: Center(" << resultPair.center.x << ", " << resultPair.center.y << "), Radius " << resultPair.radius << ", Surface Area: " << M_PI * resultPair.radius * resultPair.radius << " (" << durationPair.count() << " ms)" << endl;
    std::cout << "Smallest enclosing circle using triplet of points: Center(" << resultTriplet.center.x << ", " << resultTriplet.center.y << "), Radius " << resultTriplet.radius << ", Surface Area: " << M_PI * resultTriplet.radius * resultTriplet.radius << " (" << durationTriplet.count() << " ms)" << endl;
    std::cout << "____________________________" << std::endl;
}

int main() {
    std::vector<std::string> datasets = { "10", "20", "40", "100", "400", "800", "1000", "2000", "4000", "20000", "100000", "400000", "800000"};
    for (const auto& dataset : datasets) {
        performAnalysis(dataset);
    }

    std::cout << "Analysis performed for all datasets!" << std::endl;
    return 0;
}
