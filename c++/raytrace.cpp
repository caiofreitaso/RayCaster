#ifndef RAYTRACE_BRUTALMODE
#include <GL/glut.h>
#endif

#include <ctime>
#include "raytrace.hpp"

using namespace RayTrace;

World myWorld(0);
Camera* myCamera;
RayData<1,1,1,0> myRay;

#ifndef RAYTRACE_BRUTALMODE
void render() {
	static time_t begin, end;
	
	time(&begin);
	render(myRay,myWorld);
	time(&end);

	printf("%.f\n", difftime(end,begin));
}
void reshape(int w, int h)
{
	myRay.refresh(w,h);
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
#else
void loop()
{
	static GLdouble taxa = 10;
	while(1) {
		char key = getchar();
		if (key == '1')
			break;
		else if (key == '+')
			taxa *= 10;
		else if (key == '-')
			taxa /= 10;
		switch (key) {
			case 'w':
				myCamera->lookFrom.y -= taxa;
				myRay.changeCamera(*myCamera);
				print(myRay,myWorld);
				printf("1: Quit. Camera (%f, %f, %f); lens:%f Rate: %f\n", myCamera->lookFrom.x, myCamera->lookFrom.y, myCamera->lookFrom.z, myCamera->lensHeight, taxa);
				break;
			case 'a':
				myCamera->lookFrom.x -= taxa;
				myRay.changeCamera(*myCamera);
				print(myRay,myWorld);
				printf("1: Quit. Camera (%f, %f, %f); lens:%f Rate: %f\n", myCamera->lookFrom.x, myCamera->lookFrom.y, myCamera->lookFrom.z, myCamera->lensHeight, taxa);
				break;
			case 's':
				myCamera->lookFrom.y += taxa;
				myRay.changeCamera(*myCamera);
				print(myRay,myWorld);
				printf("1: Quit. Camera (%f, %f, %f); lens:%f Rate: %f\n", myCamera->lookFrom.x, myCamera->lookFrom.y, myCamera->lookFrom.z, myCamera->lensHeight, taxa);
				break;
			case 'd':
				myCamera->lookFrom.x += taxa;
				myRay.changeCamera(*myCamera);
				print(myRay,myWorld);
				printf("1: Quit. Camera (%f, %f, %f); lens:%f Rate: %f\n", myCamera->lookFrom.x, myCamera->lookFrom.y, myCamera->lookFrom.z, myCamera->lensHeight, taxa);
				break;
			case 'q':
				myCamera->lookFrom.z -= taxa;
				myRay.changeCamera(*myCamera);
				print(myRay,myWorld);
				printf("1: Quit. Camera (%f, %f, %f); lens:%f Rate: %f\n", myCamera->lookFrom.x, myCamera->lookFrom.y, myCamera->lookFrom.z, myCamera->lensHeight, taxa);
				break;
			case 'e':
				myCamera->lookFrom.z += taxa;
				myRay.changeCamera(*myCamera);
				print(myRay,myWorld);
				printf("1: Quit. Camera (%f, %f, %f); lens:%f Rate: %f\n", myCamera->lookFrom.x, myCamera->lookFrom.y, myCamera->lookFrom.z, myCamera->lensHeight, taxa);
				break;
			case 'z':
				myCamera->fovY -= taxa;
				myRay.changeCamera(*myCamera);
				print(myRay,myWorld);
				printf("1: Quit. Camera (%f, %f, %f); lens:%f Rate: %f\n", myCamera->lookFrom.x, myCamera->lookFrom.y, myCamera->lookFrom.z, myCamera->lensHeight, taxa);
				break;
			case 'x':
				myCamera->fovY += taxa;
				myRay.changeCamera(*myCamera);
				print(myRay,myWorld);
				printf("1: Quit. Camera (%f, %f, %f); lens:%f Rate: %f\n", myCamera->lookFrom.x, myCamera->lookFrom.y, myCamera->lookFrom.z, myCamera->lensHeight, taxa);
				break;
			case 'r':
				myCamera->lensHeight -= taxa;
				myRay.changeCamera(*myCamera);
				print(myRay,myWorld);
				printf("1: Quit. Camera (%f, %f, %f); lens:%f Rate: %f\n", myCamera->lookFrom.x, myCamera->lookFrom.y, myCamera->lookFrom.z, myCamera->lensHeight, taxa);
				break;
			case 'f':
				myCamera->lensHeight += taxa;
				myRay.changeCamera(*myCamera);
				print(myRay,myWorld);
				printf("1: Quit. Camera (%f, %f, %f); lens:%f Rate: %f\n", myCamera->lookFrom.x, myCamera->lookFrom.y, myCamera->lookFrom.z, myCamera->lensHeight, taxa);
				break;
			default:
				break;
		}
	}
}
#endif

int main(int argc, char** argv)
{
	myCamera = new Camera(Point(0,0,0),Point(-54,-30,14),Point(0,1,0),48,3000,5,0.3);
	//myCamera = new Camera(Point(0,0,0),Point(32,-20,32),Point(0,1,0),48,3000,20,0.3);
	//myCamera = new Camera(Point(0,0,0),Point(10,-20,10),Point(0,1,0),24.5,3000,7,0.3);
	myRay.changeCamera(*myCamera);
	
/*	myWorld.add(new Sphere(Point( 3, 3, 3),Point(0,1,0), 2),Material(0.4,0.5, 0,50,0.1,Color(1,1,1)));
	myWorld.add(new Sphere(Point( 3, 3,-3),Point(0,1,0), 2),Material(0.4,0.5, 0,50,0.1,Color(1,1,0)));
	myWorld.add(new Sphere(Point( 3,-3, 3),Point(0,1,0), 2),Material(0.4,0.5, 0,50,0.1,Color(1,0,1)));
	myWorld.add(new Sphere(Point( 3,-3,-3),Point(0,1,0), 2),Material(0.4,0.5, 0,50,0.1,Color(1,0,0)));
	myWorld.add(new Sphere(Point(-3, 3, 3),Point(0,1,0), 2),Material(0.4,0.5, 0,50,0.1,Color(0,1,1)));
	myWorld.add(new Sphere(Point(-3, 3,-3),Point(0,1,0), 2),Material(0.4,0.5, 0,50,0.1,Color(0,1,0)));
	myWorld.add(new Sphere(Point(-3,-3, 3),Point(0,1,0), 2),Material(0.4,0.5, 0,50,0.1,Color(0,0,1)));
	myWorld.add(new Sphere(Point(-3,-3,-3),Point(0,1,0), 2),Material(0.4,0.5, 0,50,0.1,Color(1,1,1)));

	myWorld.add(new Cube(Point(0,0,0),Point(0,1,0), 6),Material(0.4,0.5, 0,100,0.1,Color(0,0,1)));
*/
/*	MTRand random;
	for (double i = -1; i < 2; i++)
		for (double j = -1; j < 2; j++)
			myWorld.add(new Cube(Point(2*i,0,2*j),Point(0,1,0), random() * 2),Material(0.4,0.5,0,100,0.1,Color(random(),random(),random())));
*/
	myWorld.add(new Cube(Point(0,80,0),Point(0,1,0), 158),Material(0.4,0.1,0,1000,0.4,Color(0,0,1)));
	myWorld.add(new Sphere(Point(0,0,0),Point(0,1,0), 1),Material(0.7,0.3,0,50,0.1,Color(1,1,1)));

//	myWorld.add(new Sphere(Point(0,0,0),Point(0,1,0), 1),Material(0, 0.5,50,0.5,0.8,Color(1,1,1)));
//	myWorld.add(new Sphere(Point(-5,-6,10),Point(0,1,0), 0.5),Material(0, 0.5,50,0.5,0.8,Color(0,1,0)));

	myWorld.add(Light(Point(0,-10,10),Color(1,1,1),200, 2));
	myWorld.add(Light(Point(10,-10,-10),Color(1,1,1),200, 1));

	#ifndef RAYTRACE_BRUTALMODE
	init (argc, argv, 256, 200);
	#else
	myRay.refresh(64,50);
	print(myRay,myWorld);
	loop();
	
	argv[0][0] += argc;
	#endif
	return 0;
}
