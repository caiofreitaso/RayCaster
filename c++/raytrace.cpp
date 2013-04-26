#include <GL/glut.h>
#include <ctime>
#include "raytrace.hpp"

using namespace RayTrace;

World myWorld(0);
Camera* myCamera;
RayTracer* myRay;

void render() {
	static time_t begin, end;
	
	time(&begin);
	myRay->render(myWorld);
	time(&end);

	printf("%.f\n", difftime(end,begin));
}
void reshape(int w, int h)
{
	glViewport (0, 0, w, h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(0, w, 0, h, 0, 100);
	myRay->refresh();
}

void keyboard (unsigned char key, int x, int y)
{
	static GLdouble taxa = 10;
	x=y=x;
	if (key == '+')
		 taxa *= 10;
	else if (key == '-')
		 taxa /= 10;
	else
		switch (key) {
		  case 'w':
			 myCamera->lookFrom.y -= taxa;
			 myRay->changeCamera(*myCamera);
			 glutPostRedisplay();
			 break;
		  case 'a':
			 myCamera->lookFrom.x -= taxa;
			 myRay->changeCamera(*myCamera);
			 glutPostRedisplay();
			 break;
		  case 's':
			 myCamera->lookFrom.y += taxa;
			 myRay->changeCamera(*myCamera);
			 glutPostRedisplay();
			 break;
		  case 'd':
			 myCamera->lookFrom.x += taxa;
			 myRay->changeCamera(*myCamera);
			 glutPostRedisplay();
			 break;
		  case 'q':
			 myCamera->lookFrom.z -= taxa;
			 myRay->changeCamera(*myCamera);
			 glutPostRedisplay();
			 break;
		  case 'e':
			 myCamera->lookFrom.z += taxa;
			 myRay->changeCamera(*myCamera);
			 glutPostRedisplay();
			 break;
		 case 'z':
			 myCamera->fovY -= taxa;
			 myRay->changeCamera(*myCamera);
			 glutPostRedisplay();
			 break;
		  case 'x':
			 myCamera->fovY += taxa;
			 myRay->changeCamera(*myCamera);
			 glutPostRedisplay();
			 break;
		  default:
			 break;
		}
	
	char buffer[256];

	sprintf(buffer, "Camera (%f, %f, %f); Rate: %f", myCamera->lookFrom.x, myCamera->lookFrom.y, myCamera->lookFrom.z, taxa);

	glutSetWindowTitle(buffer);
}

void init (int argc, char** argv, GLint x, GLint y) {
	glutInit(&argc, argv);
	glutInitDisplayMode (GLUT_SINGLE | GLUT_RGB);
	glutInitWindowSize (x, y);
	glutInitWindowPosition (500, 500);
	glutCreateWindow ("RayCaster");
	glutDisplayFunc(render);
	glutReshapeFunc(reshape);
	glutKeyboardFunc(keyboard);
	glutMainLoop();
}


int main(int argc, char** argv)
{
	myCamera = new Camera(Point(0,-2,0),Point(-93,-67,84),Point(0,1,0),1,300,10);
	myRay = new RayTracer(*myCamera, Sampling::circle, 1, 1, 3);
	myRay->changeCamera(*myCamera);
	
	myWorld.add(new Sphere(Point( 3, 3, 3),Point(0,1,0), 2),Material(0, 0.5,50,0.5,0.1,Color(1,1,1)));
	myWorld.add(new Sphere(Point( 3, 3,-3),Point(0,1,0), 2),Material(0, 0.5,50,0.5,0.1,Color(1,1,0)));
	myWorld.add(new Sphere(Point( 3,-3, 3),Point(0,1,0), 2),Material(0, 0.5,50,0.5,0.1,Color(1,0,1)));
	myWorld.add(new Sphere(Point( 3,-3,-3),Point(0,1,0), 2),Material(0, 0.5,50,0.5,0.1,Color(1,0,0)));
	myWorld.add(new Sphere(Point(-3, 3, 3),Point(0,1,0), 2),Material(0, 0.5,50,0.5,0.1,Color(0,1,1)));
	myWorld.add(new Sphere(Point(-3, 3,-3),Point(0,1,0), 2),Material(0, 0.5,50,0.5,0.1,Color(0,1,0)));
	myWorld.add(new Sphere(Point(-3,-3, 3),Point(0,1,0), 2),Material(0, 0.5,50,0.5,0.1,Color(0,0,1)));
	myWorld.add(new Sphere(Point(-3,-3,-3),Point(0,1,0), 2),Material(0, 0.5,50,0.5,0.8,Color(1,0,0)));

//	myWorld.add(new Cube(Point(0,40,0),Point(0,1,0), 80),Material(0, 0.5,50,0.5,0.1,Color(0,0,1)));

//	myWorld.add(new Sphere(Point(0,-20,12),Point(0,1,0), 0.5),Material(0, 0.5,50,0.5,0.8,Color(0,1,0)));
//	myWorld.add(new Sphere(Point(-5,-6,10),Point(0,1,0), 0.5),Material(0, 0.5,50,0.5,0.8,Color(0,1,0)));

	myWorld.add(Light(Point(0,-20,11),Color(1,1,1),250, 0.5));
	myWorld.add(Light(Point(-20,-20,10),Color(1,1,1),500, 0.5));

	init (argc, argv, 200, 200);
	return 0;
}
