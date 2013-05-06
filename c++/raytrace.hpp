#ifndef RAYTRACE_H
#define RAYTRACE_H

#include "GL/gl.h"
#include "GL/glu.h"
#include <math.h>
#include <vector>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include "MersenneTwister.h"

#ifndef RAYTRACE_NONPARALLEL
#include <omp.h>
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

namespace RayTrace {

	const GLdouble PRECISION = 0.0000001;
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
							ret[0][i] = 6.28318530718 * random.rand();
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

		const GLdouble square_x[] = { 0, -0.5, 0.5,  0.5, -0.5,
										  0,   0.5,  0,   -0.5 };
		const GLdouble square_y[] = { 0,  0.5, 0.5, -0.5, -0.5,
										  0.5, 0,   -0.5,  0 };
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
				points = Sampling::getPoints(Sampling::circle, sampling, twister);
			else
				points = Sampling::getPoints(Sampling::circle, sampling);

			for (GLuint i = 1; i < sampling; i++)
				ret[i] = position + points[0][i]*radius*x + points[1][i]*radius*y;

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

	template<GLuint antialias, GLuint depthRays, GLuint shadows, GLuint interreflections>
	struct RayData
	{
		GLdouble compensation;
		GLdouble shadows_compensation;
		GLdouble interreflections_compensation;

		#ifdef RAYTRACE_CACHE
		std::vector<Point> intersections;
		std::vector<Color> colors;
		#endif

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

		void changeCamera(Camera& c)
		{
			camera = c;
			refreshCamera();
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

		void refresh()
		{
			for (GLint i = 0; i < viewport[2]; i++)
				delete[] buffer[i];
			delete[] buffer;
			init();
		}
	};

	typedef RayData<1,1,1,0> SimpleRT;

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

	template<GLuint AA, GLuint D, GLuint S, GLuint I>
	void prerender(RayData<AA,D,S,I>& data, World const& world)
	{
		static Point end;
		static Point depth[D];
		static Ray ray;
		static Color tmp;
		static Intersection intersected;

		#ifndef RAYTRACE_NONPARALLEL
		#pragma omp parallel for schedule(static) private(tmp,ray,intersected,end,depth)
		#endif
		for (GLint i = 0; i < data.viewport[2]; i++)
			for (GLint j = 0; j < data.viewport[3]; j++) {
				tmp = black;
				for (GLuint k = 0; k < AA; k++) {
					gluUnProject(i+Sampling::circle_x[k],data.viewport[3]-j-1+Sampling::circle_y[k],
								 0, data.modelview, data.projection, data.viewport,
								 &end.x, &end.y, &end.z);

					intersectionPoints(depth,D,data.camera.lookFrom,end,
									   data.camera.lensHeight,false);
					for (GLuint r = 0; r < D; r++) {
						ray = Line(depth[r],end).toRay(data.camera.far);
						intersected = world.intersect(ray);
						if (intersected.length > PRECISION)
							tmp += propagateRay(data,world,ray,intersected);
					}
					tmp *= data.compensation;
				}
				data.buffer[i][j] = tmp;
			}
		data.changed = false;
	}

	template<GLuint AA, GLuint D, GLuint S, GLuint I>
	Color propagateRay(RayData<AA,D,S,I>& data, World const& world, Ray ray, Intersection result)
	{
		Color ret;
		#ifdef RAYTRACE_CACHE
		bool found = false;
		#ifndef RAYTRACE_NONPARALLEL
		#pragma omp critical
		#endif
		{
			for (GLuint i = 0; i < data.intersections.size() && !found; i++)
				if (result.where == data.intersections[i]) {
					ret = data.colors[i];
					found = true;
				}
		}
		if (found)
			return ret;
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

			if (tmpRay.strength/data.camera.far > 0.01)
				if (tmpIntsc.length > 0)
					if (tmpIntsc.where != result.where)
						ret += (material.reflection) * propagateRay(data,world,tmpRay,tmpIntsc);
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

		for(GLuint i = 0; i < world.objects.size(); i++)
			if (i != result.index) {
				tmp = black;
				intersectionPoints(points, I, world.objects[i]->position,
								   result.where, world.objects[i]->scale,false);
				for(GLuint j = 0; j < I; j++) {
					tmpLine = Line(result.where,points[j]); 
					tmpRay = tmpLine.toRay(ray.strength);
					tmpIntsc = world.intersect(tmpRay);

					tmpRay.strength -= tmpIntsc.length;
					tmpRay.strength *= fmax(fmax(material.diffuse.blue,material.diffuse.red),
											material.diffuse.green);
					
					if (tmpRay.strength/data.camera.far > 0.01) {
						str = result.normal * (tmpIntsc.where - result.where).unitary();
						if (str > 0) {
							#ifdef RAYTRACE_CACHE
							found = false;
							#ifndef RAYTRACE_NONPARALLEL
							#pragma omp critical
							#endif
							{
								for (GLuint k = 0; k < data.intersections.size() && !found; k++)
									if (fabs(result.where.x - data.intersections[k].x) < 0.001 &&
										fabs(result.where.y - data.intersections[k].y) < 0.001 &&
										fabs(result.where.z - data.intersections[k].z) < 0.001) {
										found = true;
										tmp += str * data.colors[k];
									}
							}
							if (!found)
							#endif
								tmp += str * propagateRay(data,world,tmpRay,tmpIntsc);
						}
					}
				}
				tmp *= data.interreflections_compensation;
				ret += material.diffuse * (ret * tmp);
			}

		#ifdef RAYTRACE_CACHE
		#ifndef RAYTRACE_NONPARALLEL
		#pragma omp critical
		#endif
		{
			data.intersections.push_back(result.where);
			data.colors.push_back(ret);
		}
		#endif

		ret *= (data.camera.far-result.length)/data.camera.far;

		return ret;
	}

	template<GLuint AA, GLuint D, GLuint S>
	Color propagateRay(RayData<AA,D,S,0>& data, World const& world, Ray ray, Intersection result)
	{
		Color ret;
		#ifdef RAYTRACE_CACHE
		bool found = false;
		for (GLuint i = 0; i < data.intersections.size() && !found; i++)
		#ifndef RAYTRACE_NONPARALLEL
		#pragma omp critical
		#endif
		{
			if (result.where == data.intersections[i]) {
				ret = data.colors[i];
				found = true;
			}
		}
		if (found)
			return ret;
		#endif

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

			if (tmpRay.strength/data.camera.far > 0.01)
				if (tmpIntsc.length > 0)
					if (tmpIntsc.where != result.where)
						ret += (material.reflection) * propagateRay(data,world,tmpRay,tmpIntsc);
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
			data.intersections.push_back(result.where);
			data.colors.push_back(ret);
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

