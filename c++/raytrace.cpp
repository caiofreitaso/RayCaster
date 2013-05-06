#include <GL/glut.h>
#include <ctime>
#include "raytrace.hpp"

using namespace RayTrace;

World myWorld(0);
Camera* myCamera;
RayData<1,1,1,1> myRay;

void render() {
	static time_t begin, end;
	
	time(&begin);
	render(myRay,myWorld);
	time(&end);

	printf("%.f\n", difftime(end,begin));
}
void reshape(int w, int h)
{
	glViewport (0, 0, w, h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(0, w, 0, h, 0, 100);
	myRay.refresh();
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
				myRay.changeCamera(*myCamera);
				glutPostRedisplay();
				break;
			case 'a':
				myCamera->lookFrom.x -= taxa;
				myRay.changeCamera(*myCamera);
				glutPostRedisplay();
				break;
			case 's':
				myCamera->lookFrom.y += taxa;
				myRay.changeCamera(*myCamera);
				glutPostRedisplay();
				break;
			case 'd':
				myCamera->lookFrom.x += taxa;
				myRay.changeCamera(*myCamera);
				glutPostRedisplay();
				break;
			case 'q':
				myCamera->lookFrom.z -= taxa;
				myRay.changeCamera(*myCamera);
				glutPostRedisplay();
				break;
			case 'e':
				myCamera->lookFrom.z += taxa;
				myRay.changeCamera(*myCamera);
				glutPostRedisplay();
				break;
			case 'z':
				myCamera->fovY -= taxa;
				myRay.changeCamera(*myCamera);
				glutPostRedisplay();
				break;
			case 'x':
				myCamera->fovY += taxa;
				myRay.changeCamera(*myCamera);
				glutPostRedisplay();
				break;
			case 'r':
				myCamera->lensHeight -= taxa;
				myRay.changeCamera(*myCamera);
				glutPostRedisplay();
				break;
			case 'f':
				myCamera->lensHeight += taxa;
				myRay.changeCamera(*myCamera);
				glutPostRedisplay();
				break;
			default:
				break;
		}
	
	char buffer[256];

	sprintf(buffer, "Camera (%f, %f, %f); lens:%f Rate: %f", myCamera->lookFrom.x, myCamera->lookFrom.y, myCamera->lookFrom.z, myCamera->lensHeight, taxa);

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
	//myCamera = new Camera(Point(0,0,0),Point(-30,-40,32),Point(0,1,0),20,3000,30,1);
	myCamera = new Camera(Point(0,0,0),Point(10,-20,10),Point(0,1,0),30,3000,7,0.3);
	myRay.changeCamera(*myCamera);
	
	/*myWorld.add(new Sphere(Point( 3, 3, 3),Point(0,1,0), 2),Material(0, 0.5,50,0.4,0.1,Color(1,1,1)));
	myWorld.add(new Sphere(Point( 3, 3,-3),Point(0,1,0), 2),Material(0, 0.5,50,0.4,0.1,Color(1,1,0)));
	myWorld.add(new Sphere(Point( 3,-3, 3),Point(0,1,0), 2),Material(0, 0.5,50,0.4,0.1,Color(1,0,1)));
	myWorld.add(new Sphere(Point( 3,-3,-3),Point(0,1,0), 2),Material(0, 0.5,50,0.4,0.1,Color(1,0,0)));
	myWorld.add(new SpherePoint(-3, 3, 3),Point(0,1,0), 2),Material(0, 0.5,50,0.4,0.1,Color(0,1,1)));
	myWorld.add(new Sphere(Point(-3, 3,-3),Point(0,1,0), 2),Material(0, 0.5,50,0.4,0.1,Color(0,1,0)));
	myWorld.add(new Sphere(Point(-3,-3, 3),Point(0,1,0), 2),Material(0, 0.5,50,0.4,0.1,Color(0,0,1)));
	myWorld.add(new Sphere(Point(-3,-3,-3),Point(0,1,0), 2),Material(0, 0.5,50,0.4,0.1,Color(1,1,1)));

	myWorld.add(new Cube(Point(0,0,0),Point(0,1,0), 6),Material(0, 0.5,100,0.4,0.1,Color(0,0,1)));*/

	MTRand random;
	for (double i = -1; i < 2; i++)
		for (double j = -1; j < 2; j++)
			myWorld.add(new Cube(Point(2*i,0,2*j),Point(0,1,0), random() * 2),Material(0.4,0.5,0,100,0.1,Color(random(),random(),random())));

//	myWorld.add(new Cube(Point(0,80,0),Point(0,1,0), 160),Material(Color(0.3,0.4,0.5),Color(0.5,0.4,0.3),0,100,0.1,Color(1,0.5,0.5)));
//	myWorld.add(new Cube(Point(0,-1,0),Point(0,1,0), 1),Material(Color(0.5,0.4,0.3),Color(0.3,0.4,0.5),0,100,0.1,Color(0,0,1)));

//	myWorld.add(new Sphere(Point(0,0,0),Point(0,1,0), 1),Material(0, 0.5,50,0.5,0.8,Color(1,1,1)));
//	myWorld.add(new Sphere(Point(-5,-6,10),Point(0,1,0), 0.5),Material(0, 0.5,50,0.5,0.8,Color(0,1,0)));

	myWorld.add(Light(Point(5,-5,10),Color(1,1,1),400, 5));
	myWorld.add(Light(Point(-5,-5,-10),Color(1,1,1),250, 3));

	init (argc, argv, 32, 24);
	return 0;
}
