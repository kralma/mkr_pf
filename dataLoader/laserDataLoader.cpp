#include "laserDataLoader.h"

#include<fstream>
#include<sstream>
#include<cstdlib>
#include<iostream>

using namespace imr::laserDataLoader;

LaserDataLoader::LaserDataLoader(const char *_dataFile, const int nMeasurements, const std::string & _keyWord)
{
  int mcnt = 0;
  
  std::ifstream ifs(_dataFile);
  if (ifs.fail()) {
    std::cerr << "File " << _dataFile << " cannot be opened or found." << std::endl;
    exit(EXIT_FAILURE);
  }
  
  while (ifs && mcnt < nMeasurements) {
    std::string line;
    std::string keyword;
    
    // Read line from file input stream
    std::getline(ifs,line,'\n');
    
    // Construct string stream
    std::stringstream iss(line);
    std::getline(iss, keyword, ' ');
    if (keyword.compare(_keyWord) == 0) {
      std::cout << "Loading measurement " << ++mcnt;
      
      std::string temp;
      
      std::getline(iss, temp, ' ');
      const int nmeas = atoi(temp.c_str());
      std::cout << " containing " << nmeas << " measurements." << std::endl;
      
      Measurement newMeasurement;
      // Read all laser measurements
      for (int i=0; i<nmeas; i++) {
        std::getline(iss, temp, ' ');
        // std::cout << temp << std::endl;
        newMeasurement.scan.push_back(atof(temp.c_str()));
      }
      
      // Read position (x,y,theta)
      
//       double angle = 8.0/180.0*M_PI;
//       pos.x = loader[i].position.x*cos(angle) - loader[i].position.y*sin(angle);
//       pos.y = loader[i].position.x*sin(angle) + loader[i].position.y*cos(angle);
//       pos.phi = loader[i].position.phi + angle; 
      
      
      
      std::getline(iss, temp, ' ');
      newMeasurement.position.x=atof(temp.c_str());
      std::getline(iss, temp, ' ');
      newMeasurement.position.y=atof(temp.c_str());
      std::getline(iss, temp, ' ');
      newMeasurement.position.phi=atof(temp.c_str());
      
      measurements.push_back(newMeasurement);
    }
  }
  
  std::cout << "Data loaded ..." << std::endl;
}
