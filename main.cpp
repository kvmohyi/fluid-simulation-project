#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <utility>
#include "Eigen/Eigen"
#ifdef _WIN32
#include <windows.h>
#else
#include <sys/time.h>
#endif

#ifdef OSX
#include <GLUT/glut.h>
#include <OpenGL/glu.h>
#else
#include <GL/glut.h>
#include <GL/glu.h>
#endif

#include <time.h>
#include <math.h>


#define PI 3.14159265  // Should be used from mathlib
inline float sqr(float x) { return x*x; }

using namespace std;
using namespace Eigen;


void loadScene(std::string file) {

  //store variables and set stuff at the end
  int num_grids, radius, num_particles;
  float step_size, volume, density, viscosity, gas_constant;
  Vector3f particle;
  Vector3f gravity;
  std::string fname = "output.bmp";

  std::ifstream inpfile(file.c_str());
  if(!inpfile.is_open()) {
    std::cout << "Unable to open file" << std::endl;
  } else {
    std::string line;
    //MatrixStack mst;

    while(inpfile.good()) {
      std::vector<std::string> splitline;
      std::string buf;

      std::getline(inpfile,line);
      std::stringstream ss(line);

      while (ss >> buf) {
        splitline.push_back(buf);
      }
      //Ignore blank lines
      if(splitline.size() == 0) {
        continue;
      }

      //Ignore comments
      if(splitline[0][0] == '#') {
        continue;
      }

      //Valid commands:
      //size width height
      //  must be first command of file, controls image size
      else if(!splitline[0].compare("step_size")) {
        step_size = atof(splitline[1]);
      }
      //maxdepth depth
      //  max # of bounces for ray (default 5)
      else if(!splitline[0].compare("num_grids")) {
        num_grids = atoi(splitline[1]);
      }
      //output filename
      //  output file to write image to 
      else if(!splitline[0].compare("radius")) {
        radius = atoi(splitline[1]);
      }
      else if(!splitline[0].compare("num_particles")) {
        radius = atoi(splitline[1]);
      }
      else if(!splitline[0].compare("volume")) {
        volume = atof(splitline[1]);
      }
      else if(!splitline[0].compare("density")) {
        density = atof(splitline[1]);
      }
      else if(!splitline[0].compare("viscosity")) {
        viscosity = atoi(splitline[1]);
      }
      else if(!splitline[0].compare("gas_constant")) {
        gas_constant = atoi(splitline[1]);
      }
      else if(!splitline[0].compare("radius")) {
        radius = atoi(splitline[1]);
      }
      else if(!splitline[0].compare("gravity")) {
        gravity = new Vector3f(atof(splitline[1]), atof(splitline[2]), atof(splitline[3]));
      }
      else if(!splitline[0].compare("particle")) {
        particle = new Vector3f(atof(splitline[1]), atof(splitline[2]), atof(splitline[3]));
      }
    }
    inpfile.close();
  }
}

int main(int argc, char* argv[]){
  return 0;
}







