#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#include "ray-cast.h"

extern int debug;
int WIDTH;
int HEIGHT;
const int N_SPHERES = 8;
const int N_CUBES = 1;

Light myLights[2];
Sphere mySphere[8];
Cube myCube[1];
Camera myCamera;

void plotPixel(GLdouble x, GLdouble y, Color c){
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	glBegin(GL_QUADS);
	glColor3d(c.red, c.green, c.blue);
	glVertex2d(x, y);
	glVertex2d(x, y+1);
	glVertex2d(x+1, y+1);
	glVertex2d(x+1, y);
	glEnd();

	return;
}

void rayCast () {
	glMatrixMode(GL_MODELVIEW);
	glClear (GL_COLOR_BUFFER_BIT);

	GLdouble i, j;
	GLdouble dW = WIDTH, dH = HEIGHT;
	RayCaster rayzor = newRayCaster(myCamera);

	for (i = 0; i < dW; i++) {
		for (j = 0; j < dH; j++) {
			Line ray = getRay(rayzor, i, j);

			debug = 0;
			Intersection intersected = intersect(ray, myCamera.far, mySphere, N_SPHERES, myCube, N_CUBES);
			if (intersected.i >= 0) {
				GLdouble rayL = lenL(ray);
				GLdouble strength = rayL-intersected.len;
				
				plotPixel(i, j, render(ray, intersected, myLights, 2, mySphere, N_SPHERES, myCube, N_CUBES, strength));
			} else
				plotPixel(i, j, black);
		}
	}

	glFlush();
}

void reshape(int w, int h)
{
	WIDTH = w;
	HEIGHT = h;
	glViewport (0, 0, (GLsizei) w, (GLsizei) h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(0, (GLdouble)WIDTH, 0, (GLdouble)HEIGHT, 0, 100);
}

void keyboard (unsigned char key, int x, int y)
{
	static GLdouble taxa = 10;
	switch (key) {
	  case 't':
		 myLights[0].position.y -= taxa;
		 glutPostRedisplay();
		 break;
	  case 'f':
		 myLights[0].position.x -= taxa;
		 glutPostRedisplay();
		 break;
	  case 'g':
		 myLights[0].position.y += taxa;
		 glutPostRedisplay();
		 break;
	  case 'h':
		 myLights[0].position.x += taxa;
		 glutPostRedisplay();
		 break;

	  case 'r':
		 myLights[0].position.z -= taxa;
		 glutPostRedisplay();
		 break;
	  case 'y':
		 myLights[0].position.z += taxa;
		 glutPostRedisplay();
		 break;

	  case 'v':
		 myLights[0].intensity -= taxa;
		 glutPostRedisplay();
		 break;
	  case 'b':
		 myLights[0].intensity += taxa;
		 glutPostRedisplay();
		 break;

	  case '+':
		 taxa *= 10;
		 break;
	  case '-':
		 taxa /= 10;
		 break;

	  case 'w':
		 myCamera.lookFrom.y -= taxa;
		 glutPostRedisplay();
		 break;
	  case 'a':
		 myCamera.lookFrom.x -= taxa;
		 glutPostRedisplay();
		 break;
	  case 's':
		 myCamera.lookFrom.y += taxa;
		 glutPostRedisplay();
		 break;
	  case 'd':
		 myCamera.lookFrom.x += taxa;
		 glutPostRedisplay();
		 break;

	  case 'q':
		 myCamera.lookFrom.z -= taxa;
		 glutPostRedisplay();
		 break;
	  case 'e':
		 myCamera.lookFrom.z += taxa;
		 glutPostRedisplay();
		 break;
	 case 'z':
		 myCamera.fovY -= taxa;
		 glutPostRedisplay();
		 break;
	  case 'x':
		 myCamera.fovY += taxa;
		 glutPostRedisplay();
		 break;
	  default:
		 break;
	}
	
	char buffer[256];

	sprintf(buffer, "Camera (%f, %f, %f); Light(%f, %f, %f); Rate: %f",
			myCamera.lookFrom.x, myCamera.lookFrom.y, myCamera.lookFrom.z,
			myLights[0].position.x, myLights[0].position.y, myLights[0].position.z,
			taxa);

	glutSetWindowTitle(buffer);
}

void init (int argc, char** argv, GLint x, GLint y) {
	WIDTH = x;
	HEIGHT = y;
	glutInit(&argc, argv);
	glutInitDisplayMode (GLUT_SINGLE | GLUT_RGB);
	glutInitWindowSize (WIDTH, HEIGHT);
	glutInitWindowPosition (100, 100);
	glutCreateWindow ("RayCaster");
	glutDisplayFunc(rayCast);
	glutReshapeFunc(reshape);
	glutKeyboardFunc(keyboard);
	glutMainLoop();
}


int main(int argc, char** argv)
{
	myCamera = newCamera(-20,-20,32, 0,0,0, 0,1,0, 1,300,30);
	//myCamera = newCamera(-12,-25,30, 7,7,0, 0,-1,0, 1,1000,30);
	/*myCamera = newCamera(20,20,0, 5,5,0, 0,1,0, 1,1000,30);
	int i, j;
	for (i = 0; i < 10; ++i)
		for(j = 0; j < 10; ++j)
			mySphere[10*i+j] = newSphere(i,j,0, 1,0, 0.3,100,0.6,0.1, i*0.1,j*0.1,1, 0.4);*/
	mySphere[0] = newSphere( 3, 3, 3, 0.3,0, 0.2,100,0.5,0.3, 1,0,0, 1);
	mySphere[1] = newSphere( 3, 3,-3, 0.3,0, 0.2,100,0.5,0.3, 1,1,0, 1);
	mySphere[2] = newSphere( 3,-3, 3, 0.3,0, 0.2,100,0.5,0.3, 1,0,1, 1);
	mySphere[3] = newSphere( 3,-3,-3, 0.3,0, 0.2,100,0.5,0.3, 0,1,1, 1);
	mySphere[4] = newSphere(-3, 3, 3, 0.3,0, 0.2,100,0.5,0.3, 0,1,0, 1);
	mySphere[5] = newSphere(-3, 3,-3, 0.3,0, 0.2,100,0.5,0.3, 1,1,0, 1);
	mySphere[6] = newSphere(-3,-3, 3, 0.3,0, 0.2,100,0.5,0.3, 1,0,1, 1);
	mySphere[7] = newSphere(-3,-3,-3, 0.3,0, 0.2,100,0.5,0.3, 0,1,1, 1);

	myCube[0] = newCube(0,0,0, 0.99,0, 0.2,100,0.5,0.5, 0.37,0.66,0.34, 4);

	myLights[0] = newLight(5,5,10, 1,1,1, 100);
	myLights[1] = newLight(5,5,10,1,1,1,0);
	init (argc, argv, 250, 250);
	return 0;
}
