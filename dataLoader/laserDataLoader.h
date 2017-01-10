#ifndef __LASER_DATA_LOADER_H__
#define __LASER_DATA_LOADER_H__

#include<vector>
#include<string>

namespace imr {
  namespace laserDataLoader {
    
    struct RobotPosition{
      double x;
      double y;
      double phi;
      RobotPosition() : x(0.0), y(0.0), phi(0.0) {}
      
      RobotPosition(const double &_x, const double &_y, const double &_phi) : x(_x), y(_y), phi(_phi) {}
    };
    
    typedef std::vector<double> LaserScan;
    
    struct Measurement {
      RobotPosition position;
      LaserScan scan;
    };
    
    typedef std::vector<Measurement> Measurements;
    
    struct Particle {
      RobotPosition pos;
      double weight;
    };
    
    typedef std::vector<Particle> ParticleVector;
    
    
    class LaserDataLoader
    {
    private:
      Measurements measurements;
      
    public:
      LaserDataLoader(const char *_dataFile, const int nMeasurements, const std::string & _keyWord);
      
      const Measurement & operator[](int index) const {
        return measurements[index];
      }
    };
  } // namespae laserDataLoader
} // namespace imr

#endif
