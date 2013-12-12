#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <sstream>

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

#include <FreeImage.h>
#include "simulation.hpp"

#ifdef _WIN32
static DWORD lastTime;
#else
static struct timeval lastTime;
#endif

#define WAIT_KEY true
#define SAVE_IMAGE true

using namespace std;

//****************************************************
// Some Classes
//****************************************************
class Viewport {
  public:
    int w, h; // width and height
};


//****************************************************
// Global Variables
//****************************************************
Viewport viewport;
FluidSimulation* fluidsim = NULL;
bool continueSimulation = false;
int frameNumber = 0;

//****************************************************
// reshape viewport if the window is resized
//****************************************************
void myReshape(int w, int h) {
  viewport.w = w;
  viewport.h = h;

  glViewport(0,0,viewport.w,viewport.h);// sets the rectangle that will be the window
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();                // loading the identity matrix for the screen

  //----------- setting the projection -------------------------
  // glOrtho sets left, right, bottom, top, zNear, zFar of the chord system


  // glOrtho(-1, 1 + (w-400)/200.0 , -1 -(h-400)/200.0, 1, 1, -1); // resize type = add
  // glOrtho(-w/400.0, w/400.0, -h/400.0, h/400.0, 1, -1); // resize type = center

  //glOrtho(-1, 1, -1, 1, 1, -1);    // resize type = stretch
  gluPerspective(50, (float)viewport.w/(float)viewport.h, 0.001, 1000.0);
  //gluLookAt(0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
  //------------------------------------------------------------
}


//****************************************************
// sets the window up
//****************************************************
void initScene(){
  glClearColor(1.0f, 1.0f, 1.0f, 0.0f); // Clear to black, fully transparent
  glShadeModel(GL_SMOOTH);

  GLfloat mat_ambient[] = {0.0f, 0.0f, 1.0f, 1.0f};
  GLfloat mat_diffuse[] = {0.0f, 0.0f, 1.0f, 1.0f};
  GLfloat light_position[] = {0.0f, 10.0f, 0.0f};
  glMaterialfv(GL_FRONT, GL_AMBIENT, mat_ambient);
  glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
  glLightfv(GL_LIGHT0, GL_POSITION, light_position);

  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_DEPTH_TEST);
  //myReshape(viewport.w,viewport.h);
}
void drawCube(RigidBody& rigidBody) {
  glBegin(GL_QUADS);
  float sideLength = fluidsim->worldSize / 2.0;
      // front
  glVertex3f(-1.0 * sideLength, -1.0 * sideLength, sideLength);
  glVertex3f(sideLength, -1.0 * sideLength , sideLength);
  glVertex3f(sideLength, sideLength, sideLength);
  glVertex3f(-1.0 * sideLength, sideLength, sideLength);
      // back
  glVertex3f(-1.0 * sideLength, -1.0 * sideLength, -1.0 * sideLength);
  glVertex3f(sideLength, -1.0 * sideLength , -1.0 * sideLength);
  glVertex3f(sideLength, sideLength, -1.0 * sideLength);
  glVertex3f(-1.0 * sideLength, sideLength,-1.0 * sideLength);
      // right
  glVertex3f(sideLength, -1.0 * sideLength, sideLength);
  glVertex3f(sideLength, -1.0 * sideLength, -1.0f * sideLength);
  glVertex3f(sideLength, sideLength, -1.0f * sideLength);
  glVertex3f(sideLength, sideLength, sideLength);
      // left
  glVertex3f(-1.0 * sideLength, -1.0 * sideLength, sideLength);
  glVertex3f(-1.0 * sideLength, -1.0 * sideLength, -1.0 * sideLength);
  glVertex3f(-1.0 * sideLength, sideLength, -1.0 * sideLength);
  glVertex3f(-1.0 * sideLength, sideLength, sideLength);
      // top
  glVertex3f(sideLength, sideLength, sideLength);
  glVertex3f(-1.0 * sideLength, sideLength, sideLength);
  glVertex3f(-1.0 * sideLength, sideLength, -1.0 * sideLength);
  glVertex3f(sideLength, sideLength, -1.0 * sideLength);

      // bottom
  glVertex3f(sideLength, -1.0 * sideLength, sideLength);
  glVertex3f(-1.0 * sideLength, -1.0 * sideLength, sideLength);
  glVertex3f(-1.0 * sideLength, -1.0 * sideLength, -1.0 * sideLength);
  glVertex3f(sideLength, -1.0 * sideLength, -1.0 * sideLength);
  glEnd();
}

//***************************************************
// function that does the actual drawing
//***************************************************
void myDisplay2D() {
  #if SAVE_IMAGE
    BYTE* pixels = new BYTE[ 3 * viewport.w * viewport.h];
  #endif

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glPushMatrix();
  glTranslatef(0.0f, 0.0f, -2.0f);
  glRotatef(0.0, 1.0, 0.0, 0.0);
  glScalef(1.0f / fluidsim->worldSize, 1.0f / fluidsim->worldSize, 1.0f / fluidsim->worldSize);
  glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
  drawCube(fluidsim->cube);

  vector<vector<Particle> >& gridCells = fluidsim->particleList();

  for (size_t cell = 0; cell < gridCells.size(); cell++) {
    for (size_t i = 0; i < gridCells[cell].size(); i++) {
      #if false
      cout << "Iteration " << fluidsim->numIterations << ", Time " << fluidsim->numIterations * fluidsim->timeStepSize << endl;
      cout << "Position: " << gridCells[cell][i].position.x << ", " << gridCells[cell][i].position.y << ", " << gridCells[cell][i].position.z << endl;
      cout << "Velocity: " << gridCells[cell][i].velocity.x << ", " << gridCells[cell][i].velocity.y << ", " << gridCells[cell][i].velocity.z << endl;
      cout << "Acceleration: " << gridCells[cell][i].acceleration.x << ", " << gridCells[cell][i].acceleration.y << ", " << gridCells[cell][i].acceleration.z << endl;
      cout << endl;
      #endif
      glPushMatrix();
        glTranslatef(gridCells[cell][i].position.x, gridCells[cell][i].position.y, gridCells[cell][i].position.z);
        glutSolidSphere(fluidsim->sphereRadius(gridCells[cell][i]), 20, 20);
      glPopMatrix();
    }
  }

  glPopMatrix();

  #if SAVE_IMAGE
    if (fluidsim->numIterations % 4 == 0) {
      std::stringstream ss;
      if (frameNumber < 10)
        ss << "results/test" << fluidsim->testVersion << "/frame00" << frameNumber++ << ".png";
      else if (frameNumber < 100)
        ss << "results/test" << fluidsim->testVersion << "/frame0" << frameNumber++ << ".png";
      else
        ss << "results/test" << fluidsim->testVersion << "/frame" << frameNumber++ << ".png";
      string file = ss.str();
      cout << file << endl;
      glReadPixels(0, 0, viewport.w, viewport.h, GL_BGR, GL_UNSIGNED_BYTE, pixels);
      FIBITMAP* image = FreeImage_ConvertFromRawBits(pixels, viewport.w, viewport.h, 3 * viewport.w, 24, 0xFF0000, 0x00FF00, 0x0000FF, false);
      FreeImage_Save(FIF_PNG, image, file.c_str(), 0);
      FreeImage_Unload(image);
      delete [] pixels;
    }
  #endif

  glFlush();
  glutSwapBuffers();

  if (fluidsim->numIterations >= 240)
    exit(0);

  fluidsim->elapseTimeGrid();
  //cout << "derp" << endl;
  continueSimulation = false;
}



//****************************************************
// called by glut when there are no messages to handle
//****************************************************
void myFrameMove() {
  //nothing here for now
#ifdef _WIN32
  Sleep(10);                                   //give ~10ms back to OS (so as not to waste the CPU)
#endif
  glutPostRedisplay(); // forces glut to call the display function (myDisplay())
}

void keyPressed (unsigned char key, int x, int y) {  
  continueSimulation = true; 
} 

//****************************************************
// the usual stuff, nothing exciting here
//****************************************************
int main(int argc, char *argv[]) {
  FreeImage_Initialise();
  string file = argv[1];
  fluidsim = new FluidSimulation(file);

  //This initializes glut
  glutInit(&argc, argv);

  //This tells glut to use a double-buffered window with red, green, and blue channels 
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);

  // Initalize theviewport size
  viewport.w = 800;
  viewport.h = 800;

  //The size and position of the window
  glutInitWindowSize(viewport.w, viewport.h);
  glutInitWindowPosition(0, 0);
  glutCreateWindow("Fluid Simulation");

  initScene();                                 // quick function to set up scene

  glutDisplayFunc(myDisplay2D);                  // function to run when its time to draw something
  glutReshapeFunc(myReshape);                  // function to run when the window gets resized
  glutIdleFunc(myFrameMove);                   // function to run when not handling any other task
  glutKeyboardFunc(keyPressed);
  glutMainLoop();                              // infinite loop that will keep drawing and resizing and whatever else

  FreeImage_DeInitialise();

  return 0;
}
