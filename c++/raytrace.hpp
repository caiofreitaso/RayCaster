#ifndef RAYTRACE_H
#define RAYTRACE_H

#include "GL/gl.h"
#include "GL/glu.h"
#include <math.h>
#include <vector>
#include <stdio.h>
#include <stdlib.h>

namespace RayTrace {
	struct Point;
	struct Color;
}

GLdouble			operator*(RayTrace::Point a, RayTrace::Point b);
RayTrace::Point		operator*(GLdouble a, RayTrace::Point b);
RayTrace::Point		operator*(RayTrace::Point a, GLdouble b);
RayTrace::Point		operator/(RayTrace::Point a, GLdouble b);
RayTrace::Point		operator%(RayTrace::Point a, RayTrace::Point b);
RayTrace::Point		operator+(RayTrace::Point a, RayTrace::Point b);
RayTrace::Point		operator-(RayTrace::Point a, RayTrace::Point b);

bool operator==(RayTrace::Point a, RayTrace::Point b);
bool operator!=(RayTrace::Point a, RayTrace::Point b);

RayTrace::Point&	operator+=(RayTrace::Point& a, RayTrace::Point& b);
RayTrace::Point&	operator-=(RayTrace::Point& a, RayTrace::Point& b);
RayTrace::Point&	operator*=(RayTrace::Point& a, GLdouble b);
RayTrace::Point&	operator/=(RayTrace::Point& a, GLdouble b);
RayTrace::Point&	operator%=(RayTrace::Point& a, RayTrace::Point& b);



RayTrace::Color		operator+(RayTrace::Color a, RayTrace::Color b);
RayTrace::Color		operator-(RayTrace::Color a, RayTrace::Color b);
RayTrace::Color		operator*(RayTrace::Color a, RayTrace::Color b);
RayTrace::Color		operator*(GLdouble a, RayTrace::Color b);
RayTrace::Color		operator*(RayTrace::Color a, GLdouble b);
RayTrace::Color		operator/(RayTrace::Color a, RayTrace::Color b);

bool operator==(RayTrace::Color a, RayTrace::Color b);
bool operator!=(RayTrace::Color a, RayTrace::Color b);

RayTrace::Color&	operator+=(RayTrace::Color& a, RayTrace::Color b);
RayTrace::Color&	operator*=(RayTrace::Color& a, GLdouble b);

namespace RayTrace {

	const GLdouble PRECISION = 0.0000001;
	namespace Sampling {
		enum format {
			square,
			circle,
			hexagon
		};

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

		const GLdouble square_x[] = { 0, -0.5, 0.5, 0.5, -0.5, 0, 0.5, 0, -0.5 };
		const GLdouble square_y[] = { 0, 0.5, 0.5, -0.5, -0.5, 0.5, 0, -0.5, 0 };
		const GLdouble circle_x[] = { 0, -0.3535533906, 0.3535533906,
									   0.3535533906, -0.3535533906,
									  0, 0.5, 0, -0.5 };
		const GLdouble circle_y[] = { 0, 0.3535533906, 0.3535533906,
									  -0.3535533906, -0.3535533906,
									  0.5, 0, -0.5, 0 };
		const GLdouble hexagon_x[] = { 0, -0.4330127019, 0.4330127019,
									   0.4330127019, -0.4330127019,
									   0, 0 };
		const GLdouble hexagon_y[] = { 0, 0.25, 0.25, -0.25, -0.25, 0.5, -0.5 };
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

		Ray(Point a, Point b):origin(a),direction(b),strength(0) { }
		Ray(Point a, Point b, GLdouble str):origin(a),direction(b),strength(str) { }
	};

	struct Line
	{
		Point origin, destiny;

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
		GLint index;
		Intersection():length(-1),index(-1) { }
	};

	struct Camera
	{
		Point lookAt;
		Point lookFrom;
		Point up;

		GLdouble near;
		GLdouble far;
		GLdouble fovY;

		Camera(Point at, Point from, Point up, GLdouble near, GLdouble far, GLdouble fovy)
		:lookAt(at),lookFrom(from),up(up),near(near),far(far),fovY(fovy) { }
		Camera(Camera const& c)
		:lookAt(c.lookAt),lookFrom(c.lookFrom),up(c.up),near(c.near),far(c.far),fovY(c.fovY) { }
	};

	struct Material
	{
		GLdouble reflection;

		GLdouble specular;
		GLdouble shinny;
		GLdouble diffuse;
		GLdouble ambient;

		Color color;

		Material(GLdouble ref, GLdouble spec, GLdouble shine, GLdouble diff, GLdouble amb, Color c)
		: reflection(ref), specular(spec), shinny(shine), diffuse(diff), ambient(amb), color(c) { }

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
			void intersectPlane(Point const& pivot, Point const& a, Point const& b, Ray const& ray, Intersection& i)
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
				else if (lengths[1] < lengths[0])
					i = 1;

				ret.where = solutions[i];
				ret.length = lengths[i];
				ret.normal = (solutions[i] - position)/scale;
			}

			return ret;
		}
	};

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
		{
			Point* ret = new Point[sampling];

			ret[0] = position;

			if (sampling) {
				static Point normal;
				static Point x;
				static Point y;

				normal = (position - where).unitary();
				x = Point(1,0,0) * normal == 0 ?
					Point(0,0,1) - (Point(0,0,1)*normal)*normal :
					Point(1,0,0) - (Point(1,0,0)*normal)*normal;
				y = x % normal;

				if (cube)
					for (GLuint i = 0; i < sampling; i++)
						ret[i] = position + 2*Sampling::square_x[i+1]*radius*x + 2*Sampling::square_y[i+1]*radius*y;
				else {
					static GLuint n_points = 0;
					static GLdouble** points = 0;
					if (n_points != sampling) {
						n_points = sampling;
						points = Sampling::getPoints(Sampling::circle, sampling);
					}
					for (GLuint i = 1; i < sampling; i++)
						ret[i] = position + points[0][i]*radius*x + points[1][i]*radius*y;
				}
			}
			return ret;
		}
	};

	Point* intersectionPoints(GLuint sampling, Point position, Point where, GLdouble radius)
	{
		Point* ret = new Point[sampling];
		ret[0] = position;
		if (sampling > 1) {
			static Point normal;
			static Point x;
			static Point y;
			static GLdouble** points = 0;

			normal = (position - where).unitary();
			x = Point(1,0,0) * normal == 0 ?
				Point(0,0,1) - (Point(0,0,1)*normal)*normal :
				Point(1,0,0) - (Point(1,0,0)*normal)*normal;
			y = x % normal;

			points = Sampling::getPoints(Sampling::circle, sampling);
			for (GLuint i = 1; i < sampling; i++)
				ret[i] = position + points[0][i]*radius*x + points[1][i]*radius*y;

			delete[] points[0];
			delete[] points[1];
			delete[] points;
		}

		return ret;
	}

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
			ret.index = -1;
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

	struct RayTracer
	{
		GLdouble modelview[16], projection[16];
		GLint viewport[4];

		Camera camera;
		Sampling::format format;
		GLint sampling;
		GLint softShadows;
		
		Color** buffer;

		RayTracer(Camera c, Sampling::format f, GLint sampling, GLint softShadows)
		: camera(c),format(f),sampling(sampling),softShadows(softShadows)
		{
			init();
		}
		~RayTracer()
		{
			for (GLint i = 0; i < viewport[2]; i++)
				delete[] buffer[i];
			delete[] buffer;
		}

		void changeCamera(Camera c)
		{
			camera = c;
			refreshCamera();
		}
		void render(World const& world)
		{
			glMatrixMode(GL_MODELVIEW);
			glClear (GL_COLOR_BUFFER_BIT);
			
			if (changed)
				prerender(world);
			for (GLint i = 0; i < viewport[2]; i++)
				for (GLint j = 0; j < viewport[3]; j++)
					plot(buffer[i][j],i,j);

			glFlush();
			printf("OK\n");
		}

		void prerender(World const& world)
		{
			GLdouble compensation = sampling;
			compensation = 1/compensation;

			static Ray ray(origin,origin);
			static Intersection intersected;

			switch(format)
			{
				case Sampling::circle:
					for (GLint i = 0; i < viewport[2]; i++)
						for (GLint j = 0; j < viewport[3]; j++) {
							buffer[i][j] = black;
							for (GLint k = 0; k < sampling; k++) {
								ray = getRay(i + Sampling::circle_x[k],j + Sampling::circle_y[k]);
								intersected = world.intersect(ray);
								if (intersected.index >= 0)
									buffer[i][j] += shade(world,ray,intersected);
							}
							buffer[i][j] *= compensation;
						}
					break;
				case Sampling::square:
					for (GLint i = 0; i < viewport[2]; i++)
						for (GLint j = 0; j < viewport[3]; j++) {
							buffer[i][j] = black;
							for (GLint k = 0; k < sampling; k++) {
								ray = getRay(i + Sampling::square_x[k],j + Sampling::square_y[k]);
								intersected = world.intersect(ray);
								if (intersected.index >= 0)
									buffer[i][j] += shade(world,ray,intersected);
							}
							buffer[i][j] *= compensation;
						}
					break;
				case Sampling::hexagon:
					for (GLint i = 0; i < viewport[2]; i++)
						for (GLint j = 0; j < viewport[3]; j++) {
							buffer[i][j] = black;
							for (GLint k = 0; k < sampling; k++) {
								ray = getRay(i + Sampling::hexagon_x[k],j + Sampling::hexagon_y[k]);
								intersected = world.intersect(ray);
								if (intersected.index >= 0)
									buffer[i][j] += shade(world,ray,intersected);
							}
							buffer[i][j] *= compensation;
						}
					break;
			}

			
			changed = false;
		}
		void refresh()
		{
			for (GLint i = 0; i < viewport[2]; i++)
				delete[] buffer[i];
			delete[] buffer;
			init();
		}
		private:
			bool changed;
			Ray getRay(GLdouble x, GLdouble y) const
			{
				static GLdouble dY = 0;
				static Point end;
				static Ray ret(origin,origin);

				dY = viewport[3];
				dY -= y+1;
				gluUnProject(x, dY, 1, modelview, projection, viewport, &end.x, &end.y, &end.z);
				//Point front;
				//gluUnProject(x, dY, 0, modelview, projection, viewport, &front.x, &front.y, &front.z);
				//printf("RAY: %f,%f,%f -> %f,%f,%f\n", camera.lookFrom.x,camera.lookFrom.y,camera.lookFrom.z,end.x,end.y,end.z);
				ret = Line(camera.lookFrom,end).toRay();
				ret.strength = camera.far;
				return ret;
			}
			void plot(Color const& c, GLdouble x, GLdouble y) const
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
			void refreshCamera()
			{
				changed = true;
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
			}
			void init()
			{
				glGetIntegerv (GL_VIEWPORT, viewport);
				refreshCamera();
				buffer = new Color*[viewport[2]];
				for (GLint i = 0; i < viewport[2]; i++)
					buffer[i] = new Color[viewport[3]];
				
				glMatrixMode(GL_MODELVIEW);
			}
			Color shade(World const& world, Ray& ray, Intersection const& result) const
			{
				Material material = world.materials[world.objects[result.index]->material];
				Color ret = material.color * world.ambientIntensity * material.ambient;
				
				GLdouble compensation = softShadows;
				compensation = 1/compensation;

				static GLdouble str;

				static Color tmp;
				static Line tmpLine(origin,origin);
				static Ray tmpRay(origin,origin);
				static Intersection tmpIntsc;
				
				static Point* points = 0;

				if (material.reflection > 0) {
					Point origin = ray.origin.unitary();
					//printf("%f %f %f -> %f %f %f\n", l.origin.x, l.origin.y, l.origin.z, l.destiny.x, l.destiny.y, l.destiny.z);
					tmpRay = Line(result.where, result.where + ray.strength * ((2*(result.normal*origin)) * result.normal - origin)).toRay(ray.strength);

					tmpIntsc = world.intersect(tmpRay);

					tmpRay.strength -= tmpIntsc.length;
					tmpRay.strength *= material.reflection;
					
					ret *= 1 - material.reflection;

					if (tmpRay.strength > PRECISION)
						if (tmpIntsc.index >= 0)
							if (tmpIntsc.where != result.where)
								ret += (material.reflection) * shade(world,tmpRay,tmpIntsc);
				}

				for(GLuint i = 0; i < world.objects.size(); i++) {
					tmp = black;
					points = intersectionPoints(4, world.objects[i]->position,
						result.where,world.objects[i]->scale);
					for(GLuint j = 0; j < 4; j++) {
						tmpLine = Line(result.where,points[j]); 
						tmpRay = tmpLine.toRay(ray.strength);
						tmpIntsc = world.intersect(tmpRay);

						tmpRay.strength -= tmpIntsc.length;
						tmpRay.strength *= material.specular*material.ambient;
						
						//if (tmpRay.strength > PRECISION)
							if (tmpIntsc.index >= 0)
								if (tmpIntsc.where != result.where)
									tmp += material.diffuse * shade(world,tmpRay,tmpIntsc);
					}
					tmp *= 0.25;
					ret += tmp;
				}
				
				for (GLuint i = 0; i < world.lights.size(); i++) {
					tmp = black;

					points = world.lights[i].intersectionPoints(softShadows,result.where);
					for (GLint k = 0; k < softShadows; k++) {
						Point light = points[k] - result.where;
						Ray shadow = Line(points[k], result.where).toRay(world.lights[i].intensity);

						Intersection isShadow = world.intersect(shadow);
						
						if (isShadow.where == result.where)
						{
							GLdouble iLight = world.lights[i].intensity / (light.length() * light.length());

							light = light.unitary();

							GLdouble NL = result.normal * light;
							Point reflectedLight = (2*NL)*result.normal - light;
							GLdouble phi = (reflectedLight * ray.origin)/(reflectedLight.length() * ray.origin.length());

							if (phi > 0)
							{
								str = material.specular;
								str *= pow(phi, material.shinny);
								str *= iLight;
								tmp += str * world.lights[i].color;
							}
							if (material.reflection < 1) {
								str = material.diffuse * (1 - material.reflection);
								str *= NL;
								str *= iLight;
								if (str < 0) str = 0;
								tmp += material.color * (str * world.lights[i].color);
							}
						}
					}
					tmp *= compensation;
					ret += tmp;
					delete[] points;
				}

				return ret;
			}
	};

}

GLdouble			operator*(RayTrace::Point a, RayTrace::Point b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
RayTrace::Point		operator*(GLdouble a, RayTrace::Point b)
{
	b.x *= a;
	b.y *= a;
	b.z *= a;
	return b;
}
RayTrace::Point		operator*(RayTrace::Point a, GLdouble b)
{
	a.x *= b;
	a.y *= b;
	a.z *= b;
	return a;
}
RayTrace::Point		operator/(RayTrace::Point a, GLdouble b)
{
	a.x /= b;
	a.y /= b;
	a.z /= b;
	return a;
}
RayTrace::Point		operator%(RayTrace::Point a, RayTrace::Point b) { return RayTrace::Point(a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x); }
RayTrace::Point		operator+(RayTrace::Point a, RayTrace::Point b)
{
	a.x += b.x;
	a.y += b.y;
	a.z += b.z;
	return a;
}
RayTrace::Point		operator-(RayTrace::Point a, RayTrace::Point b)
{
	a.x -= b.x;
	a.y -= b.y;
	a.z -= b.z;
	return a;
}

bool operator==(RayTrace::Point a, RayTrace::Point b)
{
	return fabs(a.x - b.x) < RayTrace::PRECISION &&
		   fabs(a.y - b.y) < RayTrace::PRECISION &&
		   fabs(a.z - b.z) < RayTrace::PRECISION;
}
bool operator!=(RayTrace::Point a, RayTrace::Point b) { return !(a == b); }

RayTrace::Point&	operator+=(RayTrace::Point& a, RayTrace::Point& b)
{
	a.x += b.x;
	a.y += b.y;
	a.z += b.z;
	return a;
}
RayTrace::Point&	operator-=(RayTrace::Point& a, RayTrace::Point& b)
{
	a.x -= b.x;
	a.y -= b.y;
	a.z -= b.z;
	return a;
}
RayTrace::Point&	operator*=(RayTrace::Point& a, GLdouble b)
{
	a.x *= b;
	a.y *= b;
	a.z *= b;
	return a;
}
RayTrace::Point&	operator/=(RayTrace::Point& a, GLdouble b)
{
	a.x /= b;
	a.y /= b;
	a.z /= b;
	return a;
}
RayTrace::Point&	operator%=(RayTrace::Point& a, RayTrace::Point& b)
{
	a = a % b;
	return a;
}



RayTrace::Color		operator+(RayTrace::Color a, RayTrace::Color b)
{
	a.red += b.red;
	a.green += b.green;
	a.blue += b.blue;
	return a;
}
RayTrace::Color		operator-(RayTrace::Color a, RayTrace::Color b)
{
	a.red -= b.red;
	a.green -= b.green;
	a.blue -= b.blue;
	return a;
}
RayTrace::Color		operator*(RayTrace::Color a, RayTrace::Color b)
{
	a.red *= b.red;
	a.green *= b.green;
	a.blue *= b.blue;
	return a;
}
RayTrace::Color		operator*(GLdouble a, RayTrace::Color b)
{
	b.red *= a;
	b.green *= a;
	b.blue *= a;
	return b;
}
RayTrace::Color		operator*(RayTrace::Color a, GLdouble b)
{
	a.red *= b;
	a.green *= b;
	a.blue *= b;
	return a;
}
RayTrace::Color		operator/(RayTrace::Color a, RayTrace::Color b)
{
	a.red /= b.red;
	a.green /= b.green;
	a.blue /= b.blue;
	return a;
}

RayTrace::Color&	operator+=(RayTrace::Color& a, RayTrace::Color b)
{
	a.red += b.red;
	a.green += b.green;
	a.blue += b.blue;
	return a;
}
RayTrace::Color&	operator*=(RayTrace::Color& a, GLdouble b)
{
	a.red *= b;
	a.green *= b;
	a.blue *= b;
	return a;
}

bool operator==(RayTrace::Color a, RayTrace::Color b) { return a.red == b.red && a.green == b.green && a.blue == b.blue; }
bool operator!=(RayTrace::Color a, RayTrace::Color b) { return !(a==b); }

#endif

