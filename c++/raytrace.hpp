#ifndef RAYTRACE_H
#define RAYTRACE_H

#ifndef RAYTRACE_NONPARALLEL
#include <omp.h>
#endif
#include <math.h>
#include <vector>
#include <limits.h>

#include <iostream>
#include <stdlib.h>
#include "MersenneTwister.h"

#ifndef RAYTRACE_BRUTALMODE
#include "GL/gl.h"
#include "GL/glu.h"
#else
typedef double GLdouble;
typedef int GLint;
typedef unsigned GLuint;

void multiply(GLdouble* d, GLdouble* r, GLdouble* l)
{
	d[0] =	r[0]*l[0]  + r[4]*l[1]  + r[8]*l[2]   + r[12]*l[3];
	d[1] =	r[1]*l[0]  + r[5]*l[1]  + r[9]*l[2]   + r[13]*l[3];
	d[2] =	r[2]*l[0]  + r[6]*l[1]  + r[10]*l[2]  + r[14]*l[3];
	d[3] =	r[3]*l[0]  + r[7]*l[1]  + r[11]*l[2]  + r[15]*l[3];

	d[4] =	r[0]*l[4]  + r[4]*l[5]  + r[8]*l[6]   + r[12]*l[7];
	d[5] =	r[1]*l[4]  + r[5]*l[5]  + r[9]*l[6]   + r[13]*l[7];
	d[6] =	r[2]*l[4]  + r[6]*l[5]  + r[10]*l[6]  + r[14]*l[7];
	d[7] =  r[3]*l[4]  + r[7]*l[5]  + r[11]*l[6]  + r[15]*l[7];

	d[8] =	r[0]*l[8]  + r[4]*l[9]  + r[8]*l[10]  + r[12]*l[11];
	d[9] =	r[1]*l[8]  + r[5]*l[9]  + r[9]*l[10]  + r[13]*l[11];
	d[10] =	r[2]*l[8]  + r[6]*l[9]  + r[10]*l[10] + r[14]*l[11];
	d[11] =	r[3]*l[8]  + r[7]*l[9]  + r[11]*l[10] + r[15]*l[11];

	d[12] =	r[0]*l[12] + r[4]*l[13] + r[8]*l[14]  + r[12]*l[15];
	d[13] =	r[1]*l[12] + r[5]*l[13] + r[9]*l[14]  + r[13]*l[15];
	d[14] =	r[2]*l[12] + r[6]*l[13] + r[10]*l[14] + r[14]*l[15];
	d[15] =	r[3]*l[12] + r[7]*l[13] + r[11]*l[14] + r[15]*l[15];
}

void multiply4x1(GLdouble* d, GLdouble* r, GLdouble* l)
{
	d[0] =	r[0]*l[0]  + r[4]*l[1]  + r[8]*l[2]   + r[12]*l[3];
	d[1] =	r[1]*l[0]  + r[5]*l[1]  + r[9]*l[2]   + r[13]*l[3];
	d[2] =	r[2]*l[0]  + r[6]*l[1]  + r[10]*l[2]  + r[14]*l[3];
	d[3] =	r[3]*l[0]  + r[7]*l[1]  + r[11]*l[2]  + r[15]*l[3];
}

GLint invert(GLdouble* d, GLdouble* m)
{
	GLdouble** helper = new GLdouble*[4];
	GLdouble* tmp = 0;
	
	GLdouble val = 0;
	
	GLint i = 0, j = 0, k = 0;

	for(; i < 4; i++) {
		helper[i] = new GLdouble[8];
		for(j = 0; j < 4; j++)
			helper[i][j] = m[j*4+i];
		for(; j < 8; j++)
			helper[i][j] = (4*i+j-4) % 5 == 0;
	}

	for (i = 0; i < 4; i++) {
		for (j = 3; j > i; j--)
			if (fabs(helper[j][i]) > fabs(helper[j-1][i])) {
				tmp = helper[j-1];
				helper[j-1] = helper[j];
				helper[j] = tmp;	
			}
		if (helper[i][i] == 0)
			return 0;

		val = helper[i][i];
		for (j = 0; j < 8; j++)
			helper[i][j] /= val;

		for (k = 0; k < 4; k++) {
			val = helper[k][i];
			for (j = 0; j < 8; j++)
				if (k != i && helper[i][j] != 0)
					helper[k][j] -= helper[i][j] * val;
		}
	}

	for(i = 0; i < 4; i++)
		for(j = 0; j < 4; j++)
			d[j*4+i] = helper[i][j+4];
	return 1;
}

GLint UnProject(GLdouble winx, GLdouble winy, GLdouble winz,
				GLdouble *modelview, GLdouble *projection, GLint *viewport,
				GLdouble *x, GLdouble* y, GLdouble* z)
{
	GLdouble inverse[16], A[16];
	GLdouble in[4], out[4];
	multiply(A, projection, modelview);

	if(!invert(inverse, A))
		return 0;

	in[0]=(winx-(GLdouble)viewport[0])/(GLdouble)viewport[2]*2.0-1.0;
	in[1]=(winy-(GLdouble)viewport[1])/(GLdouble)viewport[3]*2.0-1.0;
	in[2]=2.0*winz-1.0;
	in[3]=1.0;
	
	multiply4x1(out, inverse, in);
	
	if(out[3] == 0)
		return 0;
	out[3] = 1.0/out[3];
	
	*x = out[0]*out[3];
	*y = out[1]*out[3];
	*z = out[2]*out[3];
	
	return 1;
}
#endif

namespace RayTrace {
	struct Point;
	struct Ray;
	struct Color;

	struct Intersection;
	struct World;
	struct RayTracer;
	
	void prerender(RayTracer& rt, World const& world);
	Color shade(RayTracer const& rt, World const& world, Ray& ray, Intersection const& result);
}

inline GLdouble				operator*(RayTrace::Point a, RayTrace::Point b);
inline RayTrace::Point		operator*(GLdouble a, RayTrace::Point b);
inline RayTrace::Point		operator*(RayTrace::Point a, GLdouble b);
inline RayTrace::Point		operator/(RayTrace::Point a, GLdouble b);
inline RayTrace::Point		operator%(RayTrace::Point a, RayTrace::Point b);
inline RayTrace::Point		operator+(RayTrace::Point a, RayTrace::Point b);
inline RayTrace::Point		operator-(RayTrace::Point a, RayTrace::Point b);

inline bool operator==(RayTrace::Point a, RayTrace::Point b);
inline bool operator!=(RayTrace::Point a, RayTrace::Point b);

inline RayTrace::Point&	operator+=(RayTrace::Point& a, RayTrace::Point b);
inline RayTrace::Point&	operator-=(RayTrace::Point& a, RayTrace::Point b);
inline RayTrace::Point&	operator*=(RayTrace::Point& a, GLdouble b);
inline RayTrace::Point&	operator/=(RayTrace::Point& a, GLdouble b);
inline RayTrace::Point&	operator%=(RayTrace::Point& a, RayTrace::Point b);



inline RayTrace::Color		operator+(RayTrace::Color a, RayTrace::Color b);
inline RayTrace::Color		operator-(RayTrace::Color a, RayTrace::Color b);
inline RayTrace::Color		operator*(RayTrace::Color a, RayTrace::Color b);
inline RayTrace::Color		operator*(GLdouble a, RayTrace::Color b);
inline RayTrace::Color		operator*(RayTrace::Color a, GLdouble b);
inline RayTrace::Color		operator/(RayTrace::Color a, RayTrace::Color b);

inline bool operator==(RayTrace::Color a, RayTrace::Color b);
inline bool operator!=(RayTrace::Color a, RayTrace::Color b);

inline RayTrace::Color&	operator+=(RayTrace::Color& a, RayTrace::Color b);
inline RayTrace::Color&	operator*=(RayTrace::Color& a, GLdouble b);

const GLdouble MATH_PI = 3.14159265359;
const GLdouble MATH_2PI = 6.28318530718;

namespace RayTrace {

	const GLdouble PRECISION = 0.0000001;
	const GLdouble SUBPRECISION = 0.01;
	namespace Sampling {
		enum format {
			square,
			circle,
			hexagon
		};

		GLdouble** getPoints(format shape, GLuint count, MTRand& random)
		{
			GLdouble** ret = new GLdouble*[2];
	
			ret[0] = new GLdouble[count];
			ret[1] = new GLdouble[count];

			switch(shape)
			{
				case circle:
					for(GLuint i = 0; i < count; i++) {
						#ifndef RAYTRACE_NONPARALLEL
						#pragma omp critical
						#endif
						{
							ret[0][i] = MATH_2PI * random.rand();
							ret[1][i] = random.rand() + random.rand();
						}
						if (ret[1][i] > 1)
							ret[1][i] = 2 - ret[1][i];
					}
					for(GLuint i = 0; i < count; i++) {
						GLdouble a = ret[0][i];
						ret[0][i] = ret[1][i]*cos(a);
						ret[1][i] *= sin(a);
					}
				default:
					break;
			}

			return ret;
		}
		GLdouble** getPoints(format shape, GLuint count)
		{
			GLdouble** ret = new GLdouble*[2];
	
			ret[0] = new GLdouble[count];
			ret[1] = new GLdouble[count];

			switch(shape)
			{
				case circle:
					for(GLuint i = 0; i < count; i++) {
						ret[0][i] = cos(i*6.28318530718/count);
						ret[1][i] = sin(i*6.28318530718/count);
					}
				default:
					break;
			}

			return ret;
		}
	};

	struct Color
	{
		GLdouble red, green, blue;

		Color (GLdouble r, GLdouble g, GLdouble b):red(r),green(g),blue(b) { }
		Color (GLdouble gray):red(gray),green(gray),blue(gray) { }
		Color ():red(0),green(0),blue(0) { }
	};

	const Color black;

	struct Point
	{
		GLdouble x,y,z;

		Point ():x(0),y(0),z(0) { }
		Point (GLdouble i, GLdouble j, GLdouble k):x(i),y(j),z(k) { }
		
		GLdouble length() { return sqrt(*this * *this); }
		Point unitary() { return *this/length(); }
	};

	const Point origin;

	struct Ray
	{
		Point origin, direction;
		GLdouble strength;

		Ray():origin(RayTrace::origin),direction(RayTrace::origin),strength(0) { }
		Ray(Point a, Point b):origin(a),direction(b),strength(0) { }
		Ray(Point a, Point b, GLdouble str):origin(a),direction(b),strength(str) { }
	};

	struct Line
	{
		Point origin, destiny;

		Line ():origin(RayTrace::origin),destiny(RayTrace::origin) { }
		Line (Point a, Point b):origin(a),destiny(b) { }

		Point toPoint() { return origin - destiny; }
		Ray toRay() { return Ray(origin,direction(),length()); }
		Ray toRay(GLdouble strength) { return Ray(origin,direction(),strength); }
		GLdouble length() { return toPoint().length(); }
		Point direction() { return (destiny-origin).unitary(); }
	};

	struct Intersection
	{
		Point where, normal;
		GLdouble length;
		GLuint index;
		Intersection():length(-1),index(0) { }
	};

	struct Camera
	{
		Point lookAt;
		Point lookFrom;
		Point up;

		GLdouble near;
		GLdouble far;
		GLdouble fovY;

		GLdouble lensHeight;

		Camera(Point at, Point from, Point up,
			   GLdouble near, GLdouble far, GLdouble fovy, GLdouble lensHeight)
		:lookAt(at),lookFrom(from),up(up),near(near),far(far),fovY(fovy),lensHeight(lensHeight) { }

		Camera(Camera const& c)
		:lookAt(c.lookAt),lookFrom(c.lookFrom),up(c.up),near(c.near),far(c.far),
		 fovY(c.fovY),lensHeight(c.lensHeight) { }
	};

	struct Material
	{
		Color diffuse;
		Color specular;
		
		GLdouble reflection;
		GLdouble shinny;
		GLdouble ambient;

		Color color;

		Material(Color diff, Color spec, GLdouble ref, GLdouble shine, GLdouble amb, Color c)
		: specular(spec), diffuse(diff), reflection(ref), shinny(shine), ambient(amb), color(c) { }

		bool operator==(Material const& m)
		{
			return reflection == m.reflection && specular == m.specular && 
			shinny == m.shinny && diffuse == m.diffuse && color == m.color;
		}
	};

	struct Object
	{
		Point position;
		Point up;
		GLint material;
		GLdouble scale;

		Object(Point pos, Point up, GLint material, GLdouble scale)
		: position(pos),up(up),material(material),scale(scale) { }
		virtual Intersection intersect(Ray const&) = 0;
	};

	struct Cube : Object
	{
		Cube (Point pos, Point up, GLdouble side)
		: Object(pos,up,0,side) { }

		Intersection intersect(Ray const& ray)
		{
			Intersection ret;
			ret.length = -1;

			Point vertices[8];
			vertices[0] = position + Point(scale/2,scale/2,-scale/2);
			vertices[1] = position + Point(scale/2,scale/2,scale/2);
			vertices[2] = position + Point(-scale/2,scale/2,-scale/2);
			vertices[3] = position + Point(-scale/2,scale/2,scale/2);
			vertices[4] = position + Point(-scale/2,-scale/2,scale/2);
			vertices[5] = position + Point(scale/2,-scale/2,scale/2);
			vertices[6] = position + Point(scale/2,-scale/2,-scale/2);
			vertices[7] = position + Point(-scale/2,-scale/2,-scale/2);

			intersectPlane(vertices[3],vertices[1],vertices[2], ray, ret); //top
			intersectPlane(vertices[6],vertices[5],vertices[7], ray, ret); //bottom
			intersectPlane(vertices[1],vertices[3],vertices[5], ray, ret); //front
			intersectPlane(vertices[7],vertices[2],vertices[6], ray, ret); //back
			intersectPlane(vertices[0],vertices[1],vertices[6], ray, ret); //right
			intersectPlane(vertices[4],vertices[3],vertices[7], ray, ret); //left
			return ret;
		}

		private:
			void intersectPlane(Point const& pivot, Point const& a, Point const& b,
								Ray const& ray, Intersection& i)
			{
				Point n = Line(pivot,a).direction() % Line(pivot,b).direction();
				GLdouble denominator = ray.direction * n;
				if (denominator != 0) {
					GLdouble numerator = (pivot - ray.origin) * n;
					numerator /= denominator;

					if (numerator < PRECISION)
						return;

					Point p(ray.origin + numerator*ray.direction);
					Point q(p - position);

					GLdouble x = fabs(q.x)-scale/2;
					GLdouble y = fabs(q.y)-scale/2;
					GLdouble z = fabs(q.z)-scale/2;

					if (x > PRECISION || y > PRECISION || z > PRECISION)
						return;

					GLdouble len = Line(p,ray.origin).length();
					if (i.length == -1 || len - i.length < 0) {
						i.where = p;
						i.normal = n;
						i.length = len;
					}
				}
			}
	};

	struct Sphere : Object
	{
		Sphere(Point pos, Point up, GLdouble radius)
		: Object(pos,up,0,radius) { }

		Intersection intersect(Ray const& ray)
		{
			Point oc = ray.origin - position;

			GLdouble b = ray.direction * oc;
			GLdouble c = oc * oc;
			GLdouble delta = b*b - c + scale*scale;

			Intersection ret;
			ret.length = -1;

			if (fabs(delta) >= PRECISION) {
				Point solutions[2];
				GLdouble roots[2];
				GLdouble lengths[2];
				GLint i = 0;

				roots[0] = -b + sqrt(delta);
				roots[1] = -b - sqrt(delta);

				solutions[0] = ray.origin + roots[0]*ray.direction;
				solutions[1] = ray.origin + roots[1]*ray.direction;

				lengths[0] = Line(ray.origin, solutions[0]).length();
				lengths[1] = Line(ray.origin, solutions[1]).length();

				if(roots[0] < PRECISION)
					if(roots[1] < PRECISION)
						return ret;
					else
						i = 1;
				else if (roots[1] > PRECISION && lengths[1] < lengths[0])
					i = 1;

				ret.where = solutions[i];
				ret.length = lengths[i];
				ret.normal = (solutions[i] - position)/scale;
			}

			return ret;
		}
	};

	Point* intersectionPoints(GLuint sampling, Point position, Point where,
							  GLdouble radius, bool random = true)
	{
		static MTRand twister;

		Point* ret = new Point[sampling];
		ret[0] = position;
		if (sampling > 1) {
			Point normal;
			Point x;
			Point y;

			GLdouble** points = 0;

			normal = (position - where).unitary();
			x = Point(1,0,0) * normal == 0 ?
				Point(0,0,1) - (Point(0,0,1)*normal)*normal :
				Point(1,0,0) - (Point(1,0,0)*normal)*normal;
			y = x % normal;

			if (random)
				points = Sampling::getPoints(Sampling::circle, sampling, twister);
			else
				points = Sampling::getPoints(Sampling::circle, sampling);

			for (GLuint i = 1; i < sampling; i++)
				ret[i] = position + points[0][i]*radius*x + points[1][i]*radius*y;

			delete[] points[0];
			delete[] points[1];
			delete[] points;
		}

		return ret;
	}

	inline void intersectionPoints(Point* ret, GLuint sampling, Point position, Point where,
								   GLdouble radius, bool random = true)
	{
		static MTRand twister;
		ret[0] = position;
		if (sampling > 1) {
			Point normal;
			Point x;
			Point y;

			GLdouble** points = 0;

			normal = (position - where).unitary();
			x = Point(1,0,0) * normal == 0 ?
				Point(0,0,1) - (Point(0,0,1)*normal)*normal :
				Point(1,0,0) - (Point(1,0,0)*normal)*normal;
			y = x % normal;

			if (random)
				points = Sampling::getPoints(Sampling::circle, sampling-1, twister);
			else
				points = Sampling::getPoints(Sampling::circle, sampling-1);

			for (GLuint i = 1; i < sampling; i++)
				ret[i] = position + points[0][i-1]*radius*x + points[1][i-1]*radius*y;

			delete[] points[0];
			delete[] points[1];
			delete[] points;
		}
	}

	struct Light
	{
		Point position;
		Color color;
		GLdouble intensity;
		GLdouble radius;
		bool cube;

		Light(Point p, Color c, GLdouble i, GLdouble r, bool cube = false)
		: position(p),color(c),intensity(i),radius(r),cube(cube) { }

		Point* intersectionPoints(GLuint sampling, Point where) const
		{ return RayTrace::intersectionPoints(sampling,position,where,radius); }
	};

	struct World {
		std::vector<Object*> objects;
		std::vector<Material> materials;
		std::vector<Light> lights;
		GLdouble ambientIntensity;
		
		World(GLdouble light):ambientIntensity(light) { }
		void add(Light const& l) { lights.push_back(l); }
		void add(Object* const& obj, Material const& m)
		{
			objects.push_back(obj);
			for (GLuint i = 0; i < materials.size(); i++)
				if (materials[i] == m) {
					objects[objects.size()-1]->material = i;
					return;
				}
			materials.push_back(m);
			objects[objects.size()-1]->material = materials.size() - 1;
		}
		Intersection intersect(Ray const& ray) const
		{
			Intersection ret, tmp;
			GLdouble distance = ray.strength;
			for(GLuint i = 0; i < objects.size(); i++) {
				tmp = objects[i]->intersect(ray);
				if (tmp.length >= 0)
					if (tmp.length < distance) {
						distance = tmp.length;
						ret = tmp;
						ret.index = i;
					}
			}

			return ret;
		}
	};

	#ifndef RAYTRACE_BRUTALMODE
	void plot(Color c, GLdouble x, GLdouble y)
	{
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();

		glBegin(GL_QUADS);
		glColor3d(c.red, c.green, c.blue);
		glVertex2d(x, y);
		glVertex2d(x, y+1);
		glVertex2d(x+1, y+1);
		glVertex2d(x+1, y);
		glEnd();
	}
	#endif

	template<GLuint antialias, GLuint depthRays, GLuint shadows, GLuint interreflections>
	struct RayData
	{
		GLdouble compensation;
		GLdouble shadows_compensation;
		GLdouble interreflections_compensation;

		GLdouble modelview[16], projection[16];
		GLint viewport[4];

		Color** buffer;

		Camera camera;

		bool changed;

		RayData()
		: compensation(1/((GLdouble)(antialias * depthRays))),
		  shadows_compensation(1/((GLdouble)shadows)),
		  interreflections_compensation(1/((GLdouble)interreflections)),
		  buffer(0),camera(Camera(origin,origin,origin,0,0,0,0)),changed(false)
		{ }
		RayData(Camera c)
		: compensation(1/((GLdouble)(antialias * depthRays))),
		  shadows_compensation(1/((GLdouble)shadows)),
		  interreflections_compensation(1/((GLdouble)interreflections)),
		  buffer(0),camera(c),changed(true)
		{
			init();
		}

		void refreshCamera()
		{
			changed = true;
			#ifndef RAYTRACE_BRUTALMODE
			static GLdouble m[16], p[16];
			glGetDoublev (GL_MODELVIEW_MATRIX, m);
			glGetDoublev (GL_PROJECTION_MATRIX, p);

			glMatrixMode(GL_MODELVIEW);
			glLoadIdentity();
			gluLookAt(camera.lookFrom.x, camera.lookFrom.y, camera.lookFrom.z,
					  camera.lookAt.x, camera.lookAt.y, camera.lookAt.z,
					  camera.up.x, camera.up.y, camera.up.z);

			glMatrixMode(GL_PROJECTION);
			glLoadIdentity();
			gluPerspective(camera.fovY, (GLdouble) viewport[2]/viewport[3], camera.near, camera.far);
			
			glGetDoublev (GL_MODELVIEW_MATRIX, modelview);
			glGetDoublev (GL_PROJECTION_MATRIX, projection);

			glMatrixMode(GL_MODELVIEW);
			glLoadIdentity();
			glMultMatrixd(m);

			glMatrixMode(GL_PROJECTION);
			glLoadIdentity();
			glMultMatrixd(p);
			#else
			Point x,y,z;
			z = (camera.lookFrom - camera.lookAt).unitary();
			x = (camera.up % z).unitary();
			y = z % x;
			modelview[0] = x.x;
			modelview[1] = y.x;
			modelview[2] = z.x;
			modelview[3] = 0;

			modelview[4] = x.y;
			modelview[5] = y.y;
			modelview[6] = z.y;
			modelview[7] = 0;

			modelview[8] = x.z;
			modelview[9] = y.z;
			modelview[10] = z.z;
			modelview[11] = 0;

			modelview[12] = 0;
			modelview[13] = 0;
			modelview[14] = -(camera.lookFrom - camera.lookAt).length();
			modelview[15] = 1;

			GLdouble tangent = tanf(camera.fovY/2 * (MATH_PI/180));
			GLdouble height = camera.near * tangent;
			GLdouble width = height * ((GLdouble) viewport[2]/viewport[3]);

			projection[0] = camera.near/width;
			projection[1] = 0;
			projection[2] = 0;
			projection[3] = 0;

			projection[4] = 0;
			projection[5] = camera.near/height;
			projection[6] = 0;
			projection[7] = 0;

			projection[8] = 0;
			projection[9] = 0;
			projection[10] = -(camera.far + camera.near)/(camera.far - camera.near);
			projection[11] = -1;

			projection[12] = 0;
			projection[13] = 0;
			projection[14] = -(2 * camera.far * camera.near)/(camera.far - camera.near);
			projection[15] = 0;
			#endif
		}

		void changeCamera(Camera& c)
		{
			camera = c;
			refreshCamera();
		}

		void init(GLint width = 0, GLint height = 0)
		{
			#ifndef RAYTRACE_BRUTALMODE
			if (width) {
				glViewport (0, 0, width, height);
				glMatrixMode(GL_PROJECTION);
				glLoadIdentity();
				glOrtho(0, width, 0, height, 0, 100);
			}
			glGetIntegerv (GL_VIEWPORT, viewport);
			#else
			viewport[0] = viewport[1] = 0;
			viewport[2] = width;
			viewport[3] = height;
			#endif
			refreshCamera();

			buffer = new Color*[viewport[2]];
			for (GLint i = 0; i < viewport[2]; i++)
				buffer[i] = new Color[viewport[3]];

			#ifndef RAYTRACE_BRUTALMODE
			glMatrixMode(GL_MODELVIEW);
			#endif
		}

		void refresh(GLint width = 0, GLint height = 0)
		{
			for (GLint i = 0; i < viewport[2]; i++)
				delete[] buffer[i];
			delete[] buffer;
			init(width,height);
		}
	};

	#ifdef RAYTRACE_CACHE
	struct RayCache
	{
		std::vector<Point> intersections;
		std::vector<Color> colors;
	};
	#endif


	template<GLuint AA, GLuint D, GLuint S, GLuint I>
	inline Ray getRay(RayData<AA,D,S,I>& data, GLdouble x, GLdouble y)
	{
		Point end;
		GLdouble di = x, dj = data.viewport[3];
		#ifndef RAYTRACE_BRUTALMODE
		dj -= y+1;
		gluUnProject(di, dj, 0, data.modelview, data.projection, data.viewport,
					 &end.x, &end.y, &end.z);
		#else
		UnProject(di, dj, 0, data.modelview, data.projection, data.viewport,
				  &end.x, &end.y, &end.z);
		#endif
		
		return Line(data.camera.lookFrom,end).toRay(data.camera.far);
	}

	#ifndef RAYTRACE_BRUTALMODE
	template<GLuint AA, GLuint D, GLuint S, GLuint I>
	void render(RayData<AA,D,S,I>& data, World const& world)
	{
		glMatrixMode(GL_MODELVIEW);
		glClear (GL_COLOR_BUFFER_BIT);
		
		if (data.changed)
			prerender(data,world);
		for (GLint i = 0; i < data.viewport[2]; i++)
			for (GLint j = 0; j < data.viewport[3]; j++)
				plot(data.buffer[i][j],i,j);

		glFlush();
	}
	#endif

	template<GLuint AA, GLuint D, GLuint S, GLuint I>
	void print(RayData<AA,D,S,I>& data, World const& world)
	{
		static char ASCII_ART[65] = " `-.'~,\";!+<r?/|)icvlI1}]unxzjftLCJUYXZO0Qoahkbdpqwm*WMB8&%$=#@";

		if (data.changed)
			prerender(data,world);
		if (!system("clear"))
			for (GLint i = (data.viewport[3]-1)/2; i >= 0; i--) {
				for (GLint j = 0; j < data.viewport[2]; j++) {
					Color target = (data.buffer[j][2*i] + data.buffer[j][2*i+1])/2;
					GLdouble luma = 0.2126 * target.red + 0.7152 * target.green + 0.0722 * target.blue;
					GLint index = luma >= 1 ? 63 : luma/0.015625;


					std::cout << ASCII_ART[index];
				}
				std::cout << "\n";
			}
		std::cout << data.viewport[2] << "x"<<data.viewport[3]<<"\n";
	}

	template<GLuint AA, GLuint D, GLuint S, GLuint I>
	void prerender(RayData<AA,D,S,I>& data, World const& world)
	{
		static Point end;
		static Point depth[D];
		static Ray ray;
		static Color tmp;
		static Intersection intersected;
		#ifdef RAYTRACE_CACHE
		static RayCache cache;
		#endif

		static bool done = false;
		static GLdouble** circle = Sampling::getPoints(Sampling::circle,AA);
		if (!done) {
			for (GLuint i = 0; i < AA; i++) {
				circle[0][i] *= 0.5;
				circle[1][i] *= 0.5;
			}
			done = true;
		}
		
		#ifndef RAYTRACE_CACHE
		#ifndef RAYTRACE_NONPARALLEL
		#pragma omp parallel for schedule(static) private(tmp,ray,intersected,end,depth)
		#endif
		#endif
		for (GLint i = 0; i < data.viewport[2]; i++)
			for (GLint j = 0; j < data.viewport[3]; j++) {
				tmp = black;
				for (GLuint k = 0; k < AA; k++) {
					#ifndef RAYTRACE_BRUTALMODE
					gluUnProject(i + circle[0][k], data.viewport[3]-j-1+circle[1][k],
								 0, data.modelview, data.projection, data.viewport,
								 &end.x, &end.y, &end.z);
					#else
					UnProject(i + circle[0][k], data.viewport[3]-j-1+circle[1][k],
							  0, data.modelview, data.projection, data.viewport,
							  &end.x, &end.y, &end.z);
					#endif

					intersectionPoints(depth,D,data.camera.lookFrom,end,
									   data.camera.lensHeight,false);
					for (GLuint r = 0; r < D; r++) {
						ray = Line(depth[r],end).toRay(data.camera.far);
						intersected = world.intersect(ray);
						if (intersected.length > PRECISION)
							#ifdef RAYTRACE_CACHE
							tmp += propagateRay(cache,data,world,ray,intersected);
							#else
							tmp += propagateRay(data,world,ray,intersected);
							#endif
					}
				}
				tmp *= data.compensation;
				data.buffer[i][j] = tmp;
			}
		data.changed = false;
	}

	template<GLuint AA, GLuint D, GLuint S, GLuint I>
	#ifdef RAYTRACE_CACHE
	Color propagateRay(RayCache& cache, RayData<AA,D,S,I>& data, World const& world, Ray ray, Intersection result)
	#else
	Color propagateRay(RayData<AA,D,S,I>& data, World const& world, Ray ray, Intersection result)
	#endif
	{
		Color ret;

		#ifdef RAYTRACE_CACHE
		GLdouble found = 0;
		for (GLuint k = 0; k < cache.intersections.size(); k++)
			if (fabs(result.where.x - cache.intersections[k].x) < PRECISION &&
				fabs(result.where.y - cache.intersections[k].y) < PRECISION &&
				fabs(result.where.z - cache.intersections[k].z) < PRECISION) {
				found++;
				ret += cache.colors[k];
			}
		if (found > 0) {
			ret *= 1/found;
			return ret;
		}
		#endif

		Material material = world.materials[world.objects[result.index]->material];
		ret = material.color * world.ambientIntensity * material.ambient;
		
		GLdouble str,str2;

		Color tmp;
		Line tmpLine;
		Ray tmpRay;
		Intersection tmpIntsc;
		
		Point points[S > I ? S : I];

		if (material.reflection > 0) {
			Point origin = ray.origin.unitary();
			tmpRay = Line(result.where,
						  result.where + ray.strength *
						  ((2*(result.normal*origin))*result.normal - origin)).toRay(ray.strength);

			tmpIntsc = world.intersect(tmpRay);

			tmpRay.strength -= tmpIntsc.length;
			tmpRay.strength *= material.reflection;
			
			ret *= 1 - material.reflection;

			if (tmpRay.strength/data.camera.far > SUBPRECISION)
				if (tmpIntsc.length > 0)
					if (tmpIntsc.where != result.where)
						#ifdef RAYTRACE_CACHE
						ret += (material.reflection) * propagateRay(cache,data,world,tmpRay,tmpIntsc);
						#else
						ret += (material.reflection) * propagateRay(data,world,tmpRay,tmpIntsc);
						#endif

		}

		for (GLuint i = 0; i < world.lights.size(); i++) {
			str = str2 = 0;

			intersectionPoints(points, S, world.lights[i].position,
							   result.where, world.lights[i].radius);
			#ifndef RAYTRACE_CACHE
			#ifndef RAYTRACE_NONPARALLEL
			#pragma omp parallel for schedule(static) private(tmpRay,tmpIntsc) reduction(+:str,str2)
			#endif
			#endif
			for (GLuint k = 0; k < S; k++) {
				Point light = points[k] - result.where;
				tmpRay = Line(points[k],result.where).toRay(world.lights[i].intensity);

				tmpIntsc = world.intersect(tmpRay);
				
				if (tmpIntsc.where == result.where)
				{
					GLdouble iLight = world.lights[i].intensity / (light.length() * light.length());

					light = light.unitary();

					GLdouble NL = result.normal * light;
					Point reflectedLight = (2*NL)*result.normal - light;
					GLdouble phi = (reflectedLight * ray.origin) / 
								   (reflectedLight.length() * ray.origin.length());

					if (phi > 0)
						str2 += pow(phi, material.shinny) * iLight;

					if (material.reflection < 1) {
						GLdouble tmpStr = NL * iLight;
						if (tmpStr > 0)
							str += tmpStr;
					}
				}
			}
			str *= data.shadows_compensation;
			str2 *= data.shadows_compensation;
			ret += ((str * (1 - material.reflection)) * (material.diffuse * world.lights[i].color)) *
					material.color + (str2 * (world.lights[i].color * material.specular));
		}

		for(GLuint i = 0; i < world.objects.size(); i++)
			if (i != result.index) {
				tmp = black;
				intersectionPoints(points, I, world.objects[i]->position,
								   result.where, world.objects[i]->scale,false);
				for(GLuint j = 0; j < I; j++) {
					tmpRay = Line(result.where,points[j]).toRay(ray.strength);
					tmpIntsc = world.intersect(tmpRay);

					tmpRay.strength -= tmpIntsc.length;
					tmpRay.strength *= fmax(fmax(material.diffuse.blue,material.diffuse.red),
											material.diffuse.green);
					
					if (tmpRay.strength/data.camera.far > SUBPRECISION) {
						Point there = tmpIntsc.where - result.where;
						str = result.normal * there.unitary();
						if (str > 0) {
							//str /= there.length() * there.length();
							Color back;
							#ifdef RAYTRACE_CACHE
							found = 0;
							for (GLuint k = 0; k < cache.intersections.size(); k++)
								if (fabs(tmpIntsc.where.x - cache.intersections[k].x) < SUBPRECISION &&
									fabs(tmpIntsc.where.y - cache.intersections[k].y) < SUBPRECISION &&
									fabs(tmpIntsc.where.z - cache.intersections[k].z) < SUBPRECISION) {
									found++;
									back += cache.colors[k];
								}
							back *= 1/found;
							if (found == 0)
								back = propagateRay(cache,data,world,tmpRay,tmpIntsc);
							#else
								back = propagateRay(data,world,tmpRay,tmpIntsc);
							#endif
							tmp += str * back;
						}
					}
				}
				tmp *= data.interreflections_compensation;
				ret += material.diffuse * material.color * tmp;
			}
		
		#ifdef RAYTRACE_CACHE
		cache.intersections.push_back(result.where);
		cache.colors.push_back(ret);
		#endif

		ret *= (data.camera.far-result.length)/data.camera.far;

		return ret;
	}

	template<GLuint AA, GLuint D, GLuint S>
	#ifdef RAYTRACE_CACHE
	Color propagateRay(RayCache& cache, RayData<AA,D,S,0>& data, World const& world, Ray ray, Intersection result)
	#else
	Color propagateRay(RayData<AA,D,S,0>& data, World const& world, Ray ray, Intersection result)
	#endif
	{
		Color ret;

		Material material = world.materials[world.objects[result.index]->material];
		ret = material.color * world.ambientIntensity * material.ambient;
		
		GLdouble str,str2;

		Color tmp;
		Line tmpLine;
		Ray tmpRay;
		Intersection tmpIntsc;
		
		Point points[S];

		if (material.reflection > 0) {
			Point origin = ray.origin.unitary();
			tmpRay = Line(result.where,
						  result.where + ray.strength *
						  ((2*(result.normal*origin))*result.normal - origin)).toRay(ray.strength);

			tmpIntsc = world.intersect(tmpRay);

			tmpRay.strength -= tmpIntsc.length;
			tmpRay.strength *= material.reflection;
			
			ret *= 1 - material.reflection;

			if (tmpRay.strength/data.camera.far > SUBPRECISION)
				if (tmpIntsc.length > 0)
					if (tmpIntsc.where != result.where)
						#ifdef RAYTRACE_CACHE
						ret += (material.reflection) * propagateRay(cache,data,world,tmpRay,tmpIntsc);
						#else
						ret += (material.reflection) * propagateRay(data,world,tmpRay,tmpIntsc);
						#endif
		}

		for (GLuint i = 0; i < world.lights.size(); i++) {
			str = str2 = 0;

			intersectionPoints(points, S, world.lights[i].position,
							   result.where, world.lights[i].radius);
			for (GLuint k = 0; k < S; k++) {
				Point light = points[k] - result.where;
				Ray shadow = Line(points[k],result.where).toRay(world.lights[i].intensity);

				Intersection isShadow = world.intersect(shadow);
				
				if (isShadow.where == result.where)
				{
					GLdouble iLight = world.lights[i].intensity / (light.length() * light.length());

					light = light.unitary();

					GLdouble NL = result.normal * light;
					Point reflectedLight = (2*NL)*result.normal - light;
					GLdouble phi = (reflectedLight * ray.origin) / 
								   (reflectedLight.length() * ray.origin.length());

					if (phi > 0)
						str2 += pow(phi, material.shinny) * iLight;

					if (material.reflection < 1) {
						GLdouble tmpStr = NL * iLight;
						if (tmpStr > 0)
							str += tmpStr;
					}
				}
			}
			str *= data.shadows_compensation;
			str2 *= data.shadows_compensation;
			ret += ((str * (1 - material.reflection)) * material.diffuse * world.lights[i].color) *
					material.color + (str2 * world.lights[i].color * material.specular);
		}

		#ifdef RAYTRACE_CACHE
		#ifndef RAYTRACE_NONPARALLEL
		#pragma omp critical
		#endif
		{
			cache.intersections.push_back(result.where);
			cache.colors.push_back(ret);
		}
		#endif

		ret *= (data.camera.far-result.length)/data.camera.far;

		return ret;
	}
}

inline GLdouble				operator*(RayTrace::Point a, RayTrace::Point b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline RayTrace::Point		operator*(GLdouble a, RayTrace::Point b)
{
	b.x *= a;
	b.y *= a;
	b.z *= a;
	return b;
}
inline RayTrace::Point		operator*(RayTrace::Point a, GLdouble b)
{
	a.x *= b;
	a.y *= b;
	a.z *= b;
	return a;
}
inline RayTrace::Point		operator/(RayTrace::Point a, GLdouble b)
{
	a.x /= b;
	a.y /= b;
	a.z /= b;
	return a;
}
inline RayTrace::Point		operator%(RayTrace::Point a, RayTrace::Point b) { return RayTrace::Point(a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x); }
inline RayTrace::Point		operator+(RayTrace::Point a, RayTrace::Point b)
{
	a.x += b.x;
	a.y += b.y;
	a.z += b.z;
	return a;
}
inline RayTrace::Point		operator-(RayTrace::Point a, RayTrace::Point b)
{
	a.x -= b.x;
	a.y -= b.y;
	a.z -= b.z;
	return a;
}

inline bool operator==(RayTrace::Point a, RayTrace::Point b)
{
	return fabs(a.x - b.x) < RayTrace::PRECISION &&
		   fabs(a.y - b.y) < RayTrace::PRECISION &&
		   fabs(a.z - b.z) < RayTrace::PRECISION;
}
inline bool operator!=(RayTrace::Point a, RayTrace::Point b) { return !(a == b); }

inline RayTrace::Point&	operator+=(RayTrace::Point& a, RayTrace::Point b)
{
	a.x += b.x;
	a.y += b.y;
	a.z += b.z;
	return a;
}
inline RayTrace::Point&	operator-=(RayTrace::Point& a, RayTrace::Point b)
{
	a.x -= b.x;
	a.y -= b.y;
	a.z -= b.z;
	return a;
}
inline RayTrace::Point&	operator*=(RayTrace::Point& a, GLdouble b)
{
	a.x *= b;
	a.y *= b;
	a.z *= b;
	return a;
}
inline RayTrace::Point&	operator/=(RayTrace::Point& a, GLdouble b)
{
	a.x /= b;
	a.y /= b;
	a.z /= b;
	return a;
}
inline RayTrace::Point&	operator%=(RayTrace::Point& a, RayTrace::Point b)
{
	a = a % b;
	return a;
}



inline RayTrace::Color		operator+(RayTrace::Color a, RayTrace::Color b)
{
	a.red += b.red;
	a.green += b.green;
	a.blue += b.blue;
	return a;
}
inline RayTrace::Color		operator-(RayTrace::Color a, RayTrace::Color b)
{
	a.red -= b.red;
	a.green -= b.green;
	a.blue -= b.blue;
	return a;
}
inline RayTrace::Color		operator*(RayTrace::Color a, RayTrace::Color b)
{
	a.red *= b.red;
	a.green *= b.green;
	a.blue *= b.blue;
	return a;
}
inline RayTrace::Color		operator*(GLdouble a, RayTrace::Color b)
{
	b.red *= a;
	b.green *= a;
	b.blue *= a;
	return b;
}
inline RayTrace::Color		operator*(RayTrace::Color a, GLdouble b)
{
	a.red *= b;
	a.green *= b;
	a.blue *= b;
	return a;
}
inline RayTrace::Color		operator/(RayTrace::Color a, RayTrace::Color b)
{
	a.red /= b.red;
	a.green /= b.green;
	a.blue /= b.blue;
	return a;
}

inline RayTrace::Color&	operator+=(RayTrace::Color& a, RayTrace::Color b)
{
	a.red += b.red;
	a.green += b.green;
	a.blue += b.blue;
	return a;
}
inline RayTrace::Color&	operator*=(RayTrace::Color& a, GLdouble b)
{
	a.red *= b;
	a.green *= b;
	a.blue *= b;
	return a;
}

inline bool operator==(RayTrace::Color a, RayTrace::Color b) { return a.red == b.red && a.green == b.green && a.blue == b.blue; }
inline bool operator!=(RayTrace::Color a, RayTrace::Color b) { return !(a==b); }

#endif

