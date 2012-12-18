#include <math.h>
#include <stdio.h>

int debug = 0;
GLdouble RAYCASTER_PRECISION = 0.000000001;

typedef struct {
	GLdouble red;
	GLdouble green;
	GLdouble blue;
} Color;

const Color red   = { .red = 1, .green = 0, .blue = 0 };
const Color green = { .red = 0, .green = 1, .blue = 0 };
const Color blue  = { .red = 0, .green = 0, .blue = 1 };
const Color black = { .red = 0, .green = 0, .blue = 0 };
const Color white = { .red = 1, .green = 1, .blue = 1 };

typedef struct {
	GLdouble x, y, z;
} Point;

typedef struct {
	Point origin, destiny;
} Line;

Point point(GLdouble x, GLdouble y, GLdouble z);
GLdouble dot(Point a, Point b);
GLdouble len(Point a);
GLdouble len2p(Point a, Point b);
GLdouble lenL(Line a);
Point mul(GLdouble a, Point b);
Point dv(Point a, GLdouble b);
Point cross(Point a, Point b);
Point add(Point a, Point b);
Point sub(Point a, Point b);
GLint eq(Point a, Point b);
void printP(Point a);

Point toPoint(Line a);
Point direction(Line a);
Line line(Point a, Point b);

Color addC(Color a, Color b);
Color muldC(GLdouble a, Color b);
Color mul2C(Color a, Color b);
void printC(Color a);

typedef struct {
	GLdouble modelview[16], projection[16];
	GLint viewport[4];
	Point camera;
} RayCaster;

typedef struct {
	GLdouble reflection;
	GLdouble refraction;
	
	GLdouble specular;
	GLdouble shiny;
	GLdouble diffuse;
	GLdouble ambient;

	Color color;
} Material;

typedef struct {
	Point position;
	Material material;
	GLdouble radius;
} Sphere;

typedef struct {
	Point position;
	Material material;
	GLdouble side;
} Cube;

typedef struct {
	Point lookFrom;
	Point lookAt;
	Point up;

	GLdouble near;
	GLdouble far;
	GLdouble fovY;
} Camera;

typedef struct {
	Point position;
	Color color;
	GLdouble intensity;
} Light;

struct intersection {
	Point p, normal;
	GLdouble len;
};

typedef struct {
	Point p, normal;
	GLdouble len;
	GLint i;
	GLint type;
} Intersection;

Camera newCamera(GLdouble lookFromX, GLdouble lookFromY, GLdouble lookFromZ,
				 GLdouble lookAtX, GLdouble lookAtY, GLdouble lookAtZ,
				 GLdouble upX, GLdouble upY, GLdouble upZ,
				 GLdouble near, GLdouble far, GLdouble fovY)
{
	Camera ret;

	ret.lookFrom.x = lookFromX;
	ret.lookFrom.y = lookFromY;
	ret.lookFrom.z = lookFromZ;

	ret.lookAt.x = lookAtX;
	ret.lookAt.y = lookAtY;
	ret.lookAt.z = lookAtZ;

	ret.up.x = upX;
	ret.up.y = upY;
	ret.up.z = upZ;

	ret.near = near;
	ret.far = far;
	ret.fovY = fovY;

	return ret;
}

Light newLight(GLdouble x, GLdouble y, GLdouble z,
			   GLdouble r, GLdouble g, GLdouble b,
			   GLdouble i)
{
	Light ret;
	ret.position.x = x;
	ret.position.y = y;
	ret.position.z = z;

	ret.color.red = r;
	ret.color.green = g;
	ret.color.blue = b;

	ret.intensity = i;

	return ret;
}

Sphere newSphere(GLdouble x, GLdouble y, GLdouble z,
				 GLdouble reflection, GLdouble refraction,
				 GLdouble specular, GLdouble shiny, GLdouble diffuse, GLdouble ambient,
				 GLdouble red, GLdouble green, GLdouble blue,
				 GLdouble radius)
{
	Sphere ret;

	ret.position.x = x;
	ret.position.y = y;
	ret.position.z = z;

	ret.material.refraction = refraction;
	ret.material.reflection = reflection;

	ret.material.specular = specular;
	ret.material.shiny = shiny;
	ret.material.diffuse = diffuse;
	ret.material.ambient = ambient;

	ret.material.color.red   = red;
	ret.material.color.green = green;
	ret.material.color.blue  = blue;

	ret.radius = radius;

	return ret;
}

Cube newCube(GLdouble x, GLdouble y, GLdouble z,
			 GLdouble reflection, GLdouble refraction,
			 GLdouble specular, GLdouble shiny, GLdouble diffuse, GLdouble ambient,
			 GLdouble red, GLdouble green, GLdouble blue,
			 GLdouble side)
{
	Cube ret;

	ret.position.x = x;
	ret.position.y = y;
	ret.position.z = z;

	ret.material.refraction = refraction;
	ret.material.reflection = reflection;

	ret.material.specular = specular;
	ret.material.shiny = shiny;
	ret.material.diffuse = diffuse;
	ret.material.ambient = ambient;

	ret.material.color.red   = red;
	ret.material.color.green = green;
	ret.material.color.blue  = blue;

	ret.side = side;

	return ret;
}

RayCaster newRayCaster(Camera c) {
	GLdouble m[16], p[16];
	glGetDoublev (GL_MODELVIEW_MATRIX, m);
	glGetDoublev (GL_PROJECTION_MATRIX, p);

	RayCaster ret;
	glGetIntegerv (GL_VIEWPORT, ret.viewport);
	ret.camera = c.lookFrom;

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(c.lookFrom.x, c.lookFrom.y, c.lookFrom.z,
			  c.lookAt.x, c.lookAt.y, c.lookAt.z,
			  c.up.x, c.up.y, c.up.z);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(c.fovY, (GLdouble) ret.viewport[2]/ret.viewport[3], c.near, c.far);
	
	glGetDoublev (GL_MODELVIEW_MATRIX, ret.modelview);
	glGetDoublev (GL_PROJECTION_MATRIX, ret.projection);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glMultMatrixd(m);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glMultMatrixd(p);

	return ret;
}

Line getRay(RayCaster caster, GLdouble x, GLdouble y) {
	static Line ret;
	static GLdouble dY;
	dY = (GLdouble)caster.viewport[3] - y - 1;
	ret.origin = caster.camera;
	//gluUnProject (x, dY, 0.0, caster.modelview, caster.projection, caster.viewport, &ret.origin.x, &ret.origin.y, &ret.origin.z);
	gluUnProject (x, dY, 1, caster.modelview, caster.projection, caster.viewport, &ret.destiny.x, &ret.destiny.y, &ret.destiny.z);
	return ret;
}

struct intersection intersectSphere(Point l, Point origin, Sphere object) {
	Point oc = sub(origin,object.position);
	
	GLdouble b = dot(l,oc);
	GLdouble c = dot(oc,oc);

	GLdouble delta = b*b - c + object.radius*object.radius;
	struct intersection ret;
	ret.len = -1;

	if (delta < 0)
		return ret;
	else {
		Point solutions[2];

		GLdouble roots[2];
		roots[0] = -b + sqrt(delta);
		roots[1] = -b - sqrt(delta);

		solutions[0] = add(origin, mul(roots[0],l));
		solutions[1] = add(origin, mul(roots[1],l));

		GLdouble lengths[2];
		lengths[0] = len2p(origin, solutions[0]);
		lengths[1] = len2p(origin, solutions[1]);

		GLubyte index = 0;

		if (roots[0] < 0)
			if (roots[1] < 0)
				return ret;
			else
				index = 1;
		else
			index = 0;


		if (lengths[1] < lengths[0])
			index = 1;

		ret.p = solutions[index];
		ret.len = lengths[index];

		return ret;
	}
}
Intersection intersectAllSpheres(Line ray, GLdouble strength, Sphere objects[], GLint n_objects) {
	Intersection intersected;
	intersected.len = -1;
	intersected.i = -1;
	intersected.type = 0;

	if (n_objects > 0) {
		Point l = direction(ray);
		GLdouble distance = strength;

		struct intersection test;
		test.len = -1;
		
		GLint final = -1, k;
		for (k = 0; k < n_objects; ++k) {
			test = intersectSphere(l, ray.origin, objects[k]);
			if (test.len > 0)
				if (test.len < distance) {
					distance = test.len;
					intersected.len = test.len;
					intersected.p = test.p;
					final = k;
				}
		}

		intersected.i = final;
	}

	return intersected;
}

void intersectPlane(GLdouble half, Point center, Point pivot, Point a, Point b, Point origin, Point l, struct intersection* ret) {
	Point n = cross(direction(line(pivot,a)),direction(line(pivot,b)));
	GLdouble denominator = dot(l,n);
	if (denominator != 0) {
		GLdouble numerator = dot(sub(pivot,origin), n);
		numerator /= denominator;
		
		if (numerator < RAYCASTER_PRECISION)
			return;

		Point p = add(origin, mul(numerator, l));
		Point q = sub(p, center);

		GLdouble x = fabs(q.x)-half;
		GLdouble y = fabs(q.y)-half;
		GLdouble z = fabs(q.z)-half;
		
		if (x > RAYCASTER_PRECISION || y > RAYCASTER_PRECISION || z > RAYCASTER_PRECISION)
			return;

		GLdouble len = len2p(p, origin);
		if (ret->len == -1 || len - ret->len < 0) {
			ret->p = p;
			ret->normal = n;
			ret->len = len;
		}
	}
}
struct intersection intersectCube(Point l, Point origin, Cube object) {
	GLdouble half = object.side * 0.5;

	Point a, b, c, d, e, f, g, h;
	a = add(object.position, point(half, half, -half));
	b = add(object.position, point(half, half, half));
	c = add(object.position, point(-half, half, -half));
	d = add(object.position, point(-half, half, half));
	e = add(object.position, point(-half, -half, half));
	f = add(object.position, point(half, -half, half));
	g = add(object.position, point(half, -half, -half));
	h = add(object.position, point(-half, -half, -half));

	struct intersection ret;
	ret.len = -1;
	
	//TOP
	intersectPlane(half, object.position, d,b,c, origin,l, &ret);
	//BOTTOM
	intersectPlane(half, object.position, g,f,h, origin,l, &ret);
	//FRONT
	intersectPlane(half, object.position, b,d,f, origin,l, &ret);
	//BACK
	intersectPlane(half, object.position, h,c,g, origin,l, &ret);
	//RIGHT
	intersectPlane(half, object.position, a,b,g, origin,l, &ret);
	//LEFT
	intersectPlane(half, object.position, e,d,h, origin,l, &ret);
	
	return ret;
}
Intersection intersectAllCubes(Line ray, GLdouble strength, Cube objects[], GLint n_objects) {
	Intersection intersected;
	intersected.i = -1;
	intersected.len = -1;
	intersected.type = 1;

	if (n_objects > 0) {
		Point l = direction(ray);
		GLdouble distance = strength;

		struct intersection test;
		test.len = -1;
		
		GLint k;
		for (k = 0; k < n_objects; ++k) {
			test = intersectCube(l, ray.origin, objects[k]);
			if (test.len > 0)
				if (test.len < distance) {
					distance = test.len;
					intersected.len = test.len;
					intersected.normal = test.normal;
					intersected.p = test.p;
					intersected.i = k;
				}
		}
	}

	return intersected;
}

Intersection intersect(Line ray, GLdouble strength, Sphere spheres[], GLint n_spheres, Cube cubes[], GLint n_cubes) {
	Intersection cube = intersectAllCubes(ray, strength, cubes, n_cubes),
				 sphere = intersectAllSpheres(ray, strength, spheres, n_spheres);

	if (cube.len == -1)
		return sphere;
	if (sphere.len == -1)
		return cube;

	if (cube.len < sphere.len)
		return cube;
	return sphere;
}
Color render(Line ray, Intersection result, Light sources[], GLint n_sources, Sphere spheres[], GLint n_spheres, Cube cubes[], GLint n_cubes, GLdouble strength) {
	Point normal = result.type ? result.normal : dv(sub(result.p, spheres[result.i].position), spheres[result.i].radius);
	GLdouble originLen = len(ray.origin);
	
	Material material = result.type ? cubes[result.i].material : spheres[result.i].material;
	Color ret = muldC(material.ambient, material.color);

	if (material.reflection > 0) {
		Point camera = dv(ray.origin, originLen);

		Line reflectedRay = line(result.p, add(result.p, mul(strength, sub(mul(2*dot(normal,camera),normal),camera))));
		Intersection reflection = intersect(reflectedRay, strength, spheres, n_spheres, cubes, n_cubes);

		GLdouble str = strength-reflection.len;
		str *= material.reflection;
		
		if (reflection.i >= 0) {
			if (!eq(reflection.p,result.p))
			{
				ret = addC(muldC(1-material.reflection,ret),
						   muldC(material.reflection,
						   	     render(reflectedRay, reflection, sources, n_sources,
						   	     	    spheres, n_spheres, cubes, n_cubes, str))
					//green)
					);
			} else
				ret = muldC(1-material.reflection,ret);
		} else
			ret = muldC(1-material.reflection,ret);
	}

	if (material.refraction > 0) {
		Point rayDirection = direction(ray);
		GLdouble cosTi = -dot(normal,rayDirection);

		Point transmitted = add(mul(material.refraction, rayDirection),
							    mul(material.refraction * cosTi -
							    	sqrt(1 + material.refraction * material.refraction *
							    		 (cosTi * cosTi - 1)), normal));

		GLdouble cosTt = dot(normal, transmitted);
		
		if (fabs(cosTt) > RAYCASTER_PRECISION) {
			Line refractedRay = line(result.p, add(result.p, mul(strength, transmitted)));
			Intersection refraction = intersect(refractedRay, strength, spheres, n_spheres, cubes, n_cubes);

			GLdouble str = strength-refraction.len;
			str *= material.refraction;

			if (refraction.i >= 0) {
				if (!eq(refraction.p,result.p))
				{
					ret = addC(muldC(1-material.refraction,ret),
							   muldC(material.refraction,
							   	     render(refractedRay, refraction, sources, n_sources,
							   	     	    spheres, n_spheres, cubes, n_cubes, str))
					//	green)
						);
				} else
					ret = muldC(1-material.refraction,ret);
			} else
				ret = muldC(1-material.refraction,ret);		
		}
	}

	GLint i;
	for (i = 0; i < n_sources; i++) {
		GLdouble str;
		Color tmp = black;

		Point light = sub(sources[i].position, result.p);
		Line shadow = line(sources[i].position, result.p);

		Intersection isShadow = intersect(shadow, sources[i].intensity, spheres, n_spheres, cubes, n_cubes);
		
		if (eq(isShadow.p,result.p))
		{
			GLdouble lightL = len(light);
			GLdouble iLight = sources[i].intensity / (lightL * lightL);

			light = dv(light, lightL);

			GLdouble NL = dot(normal, light);
			Point reflectedLight = sub(mul(2*NL,normal),light);
			GLdouble phi = dot(reflectedLight, ray.origin)/(len(reflectedLight)*originLen);

			if (phi > 0)
			{
				str = material.specular;
				str *= pow(phi, material.shiny);
				str *= iLight;
				tmp = addC(tmp, muldC(str,sources[i].color));
			}

			if (material.reflection < 1) {
				str = material.diffuse * (1-material.reflection);
				str *= NL;
				str *= iLight;
				if (str < 0) str = 0;
				tmp = addC(tmp, mul2C(material.color, muldC(str,sources[i].color)));
			}
		}
		
		ret = addC(ret, tmp);
	}
	ret = muldC(strength/lenL(ray), ret);

	return ret;
}

Point point(GLdouble x, GLdouble y, GLdouble z) {
	Point r;
	r.x = x;
	r.y = y;
	r.z = z;
	return r;
}
GLdouble dot(Point a, Point b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
GLdouble len(Point a) { return sqrt(dot(a,a)); }
GLdouble len2p(Point a, Point b) {
	a.x = fabs(a.x - b.x);
	a.y = fabs(a.y - b.y);
	a.z = fabs(a.z - b.z);
	return sqrt(dot(a,a));
}
GLdouble lenL(Line a) { return len(toPoint(a)); }
Point mul(GLdouble a, Point b) {
	b.x *= a;
	b.y *= a;
	b.z *= a;
	return b;
}
Point cross(Point a, Point b) {
	Point ret;
	ret.x = a.y*b.z - a.z*b.y;
	ret.y = a.z*b.x - a.x*b.z;
	ret.z = a.x*b.y - a.y*b.x;
	return ret;
}
Point add(Point a, Point b) {
	a.x += b.x;
	a.y += b.y;
	a.z += b.z;
	return a;
}
Point sub(Point a, Point b) {
	a.x -= b.x;
	a.y -= b.y;
	a.z -= b.z;
	return a;
}
Point dv(Point a, GLdouble b) {
	a.x /= b;
	a.y /= b;
	a.z /= b;
	return a;
}
GLint eq(Point a, Point b) {
	return fabs(a.x - b.x) < RAYCASTER_PRECISION &&
		   fabs(a.y - b.y) < RAYCASTER_PRECISION &&
		   fabs(a.z - b.z) < RAYCASTER_PRECISION;
}
void printP(Point a) { printf("%.15f, %.15f, %.15f\n", a.x, a.y, a.z); }

Point toPoint(Line a) {
	Point p;
	p.x = fabs(a.origin.x - a.destiny.x);
	p.y = fabs(a.origin.y - a.destiny.y);
	p.z = fabs(a.origin.z - a.destiny.z);
	return p;
}
Point direction(Line a) {
	Point ret = sub(a.destiny, a.origin);
	ret = dv(ret, len(ret));
	return ret;
}
Line line(Point a, Point b) {
	Line ret;
	ret.origin = a;
	ret.destiny = b;
	return ret;
}

Color gray(GLdouble a) {
	Color r;
	r.red = r.green = r.blue = a;
	return r;
}
Color addC(Color a, Color b) {
	GLdouble tmp = a.red + b.red;
	a.red = tmp > 1 ? 1 : tmp;

	tmp = a.green + b.green;
	a.green = tmp > 1 ? 1 : tmp;

	tmp = a.blue + b.blue;
	a.blue = tmp > 1 ? 1 : tmp;

	return a;
}
Color muldC(GLdouble a, Color b) {
	b.red *= a;
	b.green *= a;
	b.blue *= a;
	return b;
}
Color mul2C(Color a, Color b) {
	b.red *= a.red;
	b.green *= a.green;
	b.blue *= a.blue;
	return b;
}
void printC(Color a) { printf("RGB: %f %f %f", a.red, a.green, a.blue); }