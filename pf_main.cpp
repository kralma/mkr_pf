/**
 * ICP localization
 * Skeleton code for teaching
 * A3M33MKR
 * Czech Technical University
 * Faculty of Electrical Engineering
 * Intelligent and Mobile Robotics group
 *
 * Authors: Miroslav Kulich, Zdeněk Kasl, Karel Košnar kulich@labe.felk.cvut.cz kosnar@labe.felk.cvut.cz
 *
 * Licence: MIT (see LICENSE file)
 **/

#include<cmath>
#include<cassert>
#include<cstdlib>
#include<fstream>
#include<iostream>
#include <cstdlib>
#include <ctime>
#include <vtkPNGReader.h>
#include <zconf.h>

#include "gui/gui.h"
#include "dataLoader/laserDataLoader.h"
#include "laserSimulator/lasersimulator.h"

#define N_OF_PARTICLES_INITIAL 1000
#define N_OF_PARTICLES_MIN 900
#define N_OF_RANDOM_PARTICLES 100
using namespace std;
using namespace imr;
using namespace gui;
using namespace laserDataLoader;


double toRadians(const double alpha);

//function performing the iterative closest point method - you must to implement it
RobotPosition icp(const Measurement reference, const Measurement actual, RobotPosition previousPosition);

//function convert RobotPosition to Point structure (for drawing)
Point robotPosition2point(const RobotPosition &rp);

//convert the laser scan into vector of points in Cartesian coordinate system using the odometry from the Measurement structure
void calculateRawPoints(RawPoints &rp, Measurement m);

//convert the laser scan into vector of points in Cartesian coordinate system using given robot position p
void calculateRawPoints(RawPoints &rp, Measurement m, RobotPosition p);

//convert point in polar coordinate system  into Cartesian coordinate system (possibly moved by Robot Position)
void polar2cartesian(Point &p, const double &alpha, const double &r, const RobotPosition &p0);

double getDistancePow(Point &p1, Point &p2);
RobotPosition randomizePosition(RobotPosition currPos);
ParticleVector getNewParticles(ParticleVector currentParticles, LaserSimulator simul);

void help(char **argv) {
    std::cout << "\nUsage of the program " << argv[0] + 2 << ":\n"
              << "Parameter [-h or -H] displays this message.\n"
              << "Parameter [-f or -F] specifies path to data."
              << "Parameter [-m or -M] specifies number of measurements taken,\n"
              << "   default number of measurements is 2.\n"
              << std::endl;
}


int main(int argc, char **argv) {
    int nMeasurements = 2;
    char *dataFile;

    // argument count must be greater than three
    // >>> at least source file must be specified
    if (argc < 3) {
        help(argv);
        return EXIT_FAILURE;
    }

    // Parse all console parameters
    for (int i = 0; i < argc; i++) {
        if (argv[i][0] == '-') {
            switch (argv[i][1]) {
                //> HELP
                case 'H' :
                case 'h' :
                    help(argv);
                    break;

                    //> Source file
                case 'F' :
                case 'f' :
                    assert(i + 1 < argc);
                    dataFile = argv[i + 1];
                    break;

                    //> Number of Measurements
                case 'M' :
                case 'm' :
                    assert(i + 1 < argc);
                    assert(atoi(argv[i + 1]) > 1);
                    nMeasurements = atoi(argv[i + 1]);
                    break;

                default :
                    std::cout << "Parameter \033[1;31m" << argv[i] << "\033[0m is not valid!\n"
                              << "Use parameter -h or -H for help." << std::endl;
                    break;
            }
        }
    }
    // All parameters parsed

    std::cout << "Number of measuremetns taken: " << nMeasurements << "\n"
              << "Source file: " << dataFile
              << std::endl;

    // Load data
    LaserDataLoader loader(dataFile, nMeasurements, "FLASER");
    LaserScan scan;
    RawPoints scanPoints;
    Measurement measurement;
    RobotPosition pos;

    // Load map as a grid
    vtkSmartPointer<vtkPNGReader> reader = vtkSmartPointer<vtkPNGReader>::New();
    reader->SetFileName("../data/belgioioso-map5b.png");
    imr::LaserSimulator simul(reader);
    Gui gui(/*initial,tentative,robotPosition2point(loader[0].position),*/reader);
//    std::cout << "X: " << simul.grid2realX(0) << " " << simul.grid2realX(1872) << std::endl;
//    std::cout << "Y: " << simul.grid2realY(0) << " " << simul.grid2realY(5015) << std::endl;
    Particle p;
    ParticleVector particles;
    double x, y, phi;
    for (size_t i = 0; i < N_OF_PARTICLES_INITIAL;) {
        x = (rand() % 3669 - 1696) / 100.0;
        y = (rand() % 9828 - 4325) / 100.0;
        phi = M_PI * 2 * drand48();
//        phi = loader[0].position.phi;
        p.pos = RobotPosition(x, y, phi);
        if (simul.isFeasible(p.pos)) {
            p.weight = 1.0 / N_OF_PARTICLES_INITIAL;
            particles.push_back(p);
            i++;
        }
    }

    clock_t begin = clock();
    RobotPosition prevPos = loader[0].position;
    for (int i = 0; i < nMeasurements; i++) {
       double angle = 8.5/180.0*M_PI;
       pos.x = loader[i].position.x*cos(angle) - loader[i].position.y*sin(angle);
       pos.y = loader[i].position.x*sin(angle) + loader[i].position.y*cos(angle);
       pos.phi = loader[i].position.phi + angle;
//        pos = loader[i].position;
//         std::cout << "POSE " << pos.x << " " << pos.y << std::endl;
        scan = simul.getScan(pos);
        scanPoints = simul.getRawPoints();
        // gui.clearPositionPoints();


        double movX = pos.x - prevPos.x;
        double movY = pos.y - prevPos.y;
        double movD = sqrt(movX * movX + movY * movY);
        double movPhi = pos.phi - prevPos.phi;

        gui.setPosition(robotPosition2point(pos));
        gui.clearMapPoints();
        gui.setPointsToMap(scanPoints, robotPosition2point(pos));
        gui.setParticlePoints(particles);

        for (int k = 0; k < particles.size(); k++) {
            RobotPosition posTmp = particles[k].pos;
            posTmp.phi += movPhi; // ??
            posTmp.x += cos(posTmp.phi) * movD;
            posTmp.y += sin(posTmp.phi) * movD;
            particles[k].pos = posTmp;


                LaserScan scanTmp = simul.getScan(posTmp);
                unsigned long n = scanTmp.size();
                double errors[n];
                double errAvg = 0;
                for (int l = 0; l < n; l++) {
                    double distTmp = abs(scanTmp[l] - scan[l]);
                    errors[l] = distTmp;
                    errAvg += distTmp;
                }
                particles[k].weight = 1/errAvg;
        }

        particles = getNewParticles(particles,simul);
//        usleep(200000);
        gui.startInteractor();
        prevPos = pos;
    }
    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    std::cout << "ELAPSED TIME: " << elapsed_secs << " s" << std::endl;
    gui.startInteractor();
    return EXIT_SUCCESS;
}

/// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
/// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

Point robotPosition2point(const RobotPosition &rp) {
    return Point(rp.x, rp.y);
}

void calculateRawPoints(RawPoints &rp, Measurement m) {
    calculateRawPoints(rp, m, m.position);
}

void calculateRawPoints(RawPoints &rp, Measurement m, RobotPosition pos) {
    const double laserResolution = 0.5; // deg
    const double laserShift = -90.0;
    for (int j = 0; j < m.scan.size(); j++) {
        if (m.scan[j] > 50.0) continue;
        Point p;
        polar2cartesian(p, toRadians(j * laserResolution + laserShift), m.scan[j], pos);
        rp.push_back(p);
    }
}

void polar2cartesian(Point &p, const double &alpha, const double &r, const RobotPosition &p0) {
    p.x = p0.x + r * cos(alpha + p0.phi);
    p.y = p0.y + r * sin(alpha + p0.phi);
}

double toRadians(const double alpha) {
    return (alpha * M_PI) / 180.0;
}

double getDistancePow(Point &p1, Point &p2) {
    return (p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y);
}

ParticleVector getNewParticles(ParticleVector currentParticles, LaserSimulator simul) {
    ParticleVector newParticles;

    int minN = N_OF_PARTICLES_MIN;
    double downSizing = 0.01 * (currentParticles.size() - minN);

    double weightSum = 0;
    for (int i = 0; i < currentParticles.size(); i++) {
        weightSum += currentParticles[i].weight;
    }
    for (int i = 0; i < currentParticles.size() - N_OF_RANDOM_PARTICLES - downSizing;) {
        double t = 0;
        double rnd = drand48() * weightSum;
        for (int j = 0; j < currentParticles.size(); j++) {

            if (currentParticles[j].weight + t < rnd) {
                t += currentParticles[j].weight;
                continue;
            }

            Particle p;
            p.pos = randomizePosition(currentParticles[j].pos);

            if (simul.isFeasible(p.pos)) {
                newParticles.push_back(p);
                i++;
                break;
            }
        }
    }
    for (int i = 0; i < N_OF_RANDOM_PARTICLES;) {
        double x = (rand() % 3669 - 1696) / 100.0;
        double y = (rand() % 9828 - 4325) / 100.0;
        double phi = M_PI * 2 * drand48();
//        phi = loader[0].position.phi;
        Particle p;
        p.pos = RobotPosition(x, y, phi);
        if (simul.isFeasible(p.pos)) {
            p.weight = 1.0 / newParticles.size() + N_OF_RANDOM_PARTICLES - i;
            newParticles.push_back(p);
            i++;
        }
    }

    std::cout<<newParticles.size()<<endl;
    return newParticles;
}

RobotPosition randomizePosition(RobotPosition currPos) {
    double posCoef = 0.1;
    double angCoef = M_PI * 0.4;
    double randAng = M_PI * 2 * drand48() - M_PI;
//    std::cout<< M_2_PI << endl;
    double phi = currPos.phi + angCoef * drand48() - 0.5 * angCoef;
    double x = currPos.x + posCoef * drand48() * cos(randAng);
    double y = currPos.y + posCoef * drand48() * sin(randAng);
    RobotPosition newPos(x, y, phi);
    return newPos;
}