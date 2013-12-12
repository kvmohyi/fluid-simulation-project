#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>

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

#include "simulation.hpp"

#ifdef _WIN32
static DWORD lastTime;
#else
static struct timeval lastTime;
#endif

#define WAIT_KEY true

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


//***************************************************
// function that does the actual drawing
//***************************************************
void myDisplay2D() {
    //if (continueSimulation) {
    if (true) {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glPushMatrix();
      glTranslatef(0.0f, 0.0f, -2.0f);
      glRotatef(10.0, 1.0, 0.0, 0.0);
      glScalef(1.0f / fluidsim->worldSize, 1.0f / fluidsim->worldSize, 1.0f / fluidsim->worldSize);

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

    glFlush();
    glutSwapBuffers();

    fluidsim->elapseTimeGrid();
    //cout << "derp" << endl;
    continueSimulation = false;
  }
}

void drawCube(RigidBody& rigidBody) {
  glBegin(GL_QUADS);
      // front
      glVertex3f(0.0f, 0.0f, 0.0f);
      glVertex3f(1.0f, 0.0f, 0.0f);
      glVertex3f(1.0f, 1.0f, 0.0f);
      glVertex3f(0.0f, 1.0f, 0.0f);
      // back
      glVertex3f(0.0f, 0.0f, -1.0f);
      glVertex3f(1.0f, 0.0f, -1.0f);
      glVertex3f(1.0f, 1.0f, -1.0f);
      glVertex3f(0.0f, 1.0f, -1.0f);
      // right
      glVertex3f(1.0f, 0.0f, 0.0f);
      glVertex3f(1.0f, 0.0f, -1.0f);
      glVertex3f(1.0f, 1.0f, -1.0f);
      glVertex3f(1.0f, 1.0f, 0.0f);
      // left
      glVertex3f(0.0f, 0.0f, 0.0f);
      glVertex3f(0.0f, 0.0f, -1.0f);
      glVertex3f(0.0f, 1.0f, -1.0f);
      glVertex3f(0.0f, 1.0f, 0.0f);
      // top
      glVertex3f(0.0f, 1.0f, 0.0f);
      glVertex3f(1.0f, 1.0f, 0.0f);
      glVertex3f(1.0f, 1.0f, -1.0f);
      glVertex3f(0.0f, 1.0f, -1.0f);
      // bottom
      glVertex3f(0.0f, 0.0f, 0.0f);
      glVertex3f(1.0f, 0.0f, 0.0f);
      glVertex3f(1.0f, 0.0f, -1.0f);
      glVertex3f(0.0f, 0.0f, -1.0f);
  glEnd();
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

  return 0;
}
