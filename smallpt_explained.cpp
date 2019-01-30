/*
	Human readable variant of
	smallpt, a Path Tracer by Kevin Beason, 2008
	modified by Vassillen Chizhov, 2019
*/

// Make : g++ -O3 -fopenmp smallpt.cpp -o smallpt
// Remove "-fopenmp" for g++ version < 4.2

#include <math.h>   
#include <stdlib.h> 
#include <stdio.h>  

#include <random>
//#include <limits>


/*
	Changes:

	modified prng to Mersenne twister
	added depth limit to avoid infinite recursion on 'bad' scenes (introduces bias)
	modified the probability for the Russian roulete to use luminance rather than max of the 3 components
	minor changes to the Vec class (now vec3)
	moved eps (and the check) from the Sphere class to the scene intersect function
	modified Sphere::intersect() to return negative if there is no intersection (saves one ?: evaluation)
	changed the arbitrary large number in intersect to std::numeric_limits<double>::infinity()
	added a function returning the normal of the sphere at some point (useful if a base class Surface is made)
*/

#define MAX_DEPTH 10

double random()
{
	static thread_local std::mt19937 prng{ std::random_device{}() };
	static thread_local std::uniform_real_distribution<double> dist(0.0f, 1.0f);
	return dist(prng);
}

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif



/////////////////////////////////////////////////////////////////////////////////////
//										VEC3									   //
/////////////////////////////////////////////////////////////////////////////////////
struct vec3 {
	// convenience
	union
	{
		struct
		{
			double x, y, z;
		};
		struct
		{
			double r, g, b;
		};
		double e[3];
	};
	// changed ctors
	vec3() = default;
	vec3(double x, double y, double z) : x(x), y(y), z(z) {}
	vec3(double s) : x(s), y(s), z(s) {}

	// added unary - and +=, -= , *=, /=
	vec3 operator-() const { return vec3(-x, -y, -z); }
	vec3& operator+=(const vec3& b)
	{
		x += b.x;
		y += b.y;
		z += b.z;
		return *this;
	}
	vec3& operator-=(const vec3& b)
	{
		x -= b.x;
		y -= b.y;
		z -= b.z;
		return *this;
	}
	vec3& operator*=(const vec3& b)
	{
		x *= b.x;
		y *= b.y;
		z *= b.z;
		return *this;
	}
	vec3& operator*=(double b)
	{
		x *= b;
		y *= b;
		z *= b;
		return *this;
	}
	vec3& operator/=(const vec3& b)
	{
		x /= b.x;
		y /= b.y;
		z /= b.z;
		return *this;
	}
	vec3& operator/=(double b)
	{
		double invB = 1 / b;
		x *= invB;
		y *= invB;
		z *= invB;
		return *this;
	}
};

// took out the binary operators to emphasize commutativity
vec3 operator+(const vec3& a, const vec3 &b) { return vec3(a.x + b.x, a.y + b.y, a.z + b.z); }
vec3 operator-(const vec3& a, const vec3 &b) { return vec3(a.x - b.x, a.y - b.y, a.z - b.z); }
vec3 operator*(const vec3& a, const vec3 &b) { return vec3(a.x*b.x, a.y*b.y, a.z*b.z); }
vec3 operator*(const vec3& a, double b) { return vec3(a.x*b, a.y*b, a.z*b); }
// provided right side scalar multiplication
vec3 operator*(double a, const vec3& b) { return vec3(a*b.x, a*b.y, a*b.z); }
// added / operator
vec3 operator/(const vec3& a, const vec3 &b) { return vec3(a.x / b.x, a.y / b.y, a.z / b.z); }
vec3 operator/(const vec3& a, double b) { double invB = 1 / b; return vec3(a.x*invB, a.y*invB, a.z*invB); }
vec3 operator/(double a, const vec3& b) { return vec3(a / b.x, a / b.y, a / b.z); }
double dot(const vec3& a, const vec3 &b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
// renamed % to cross
vec3 cross(const vec3& a, const vec3& b) { return vec3(a.y*b.z - a.z * b.y, a.z*b.x - a.x * b.z, a.x*b.y - a.y * b.x); }
// removed .norm() and made normalize()
vec3 normalize(const vec3& v) { return v / sqrt(dot(v, v)); }


struct Ray 
{ 
	vec3 o, d; 
	Ray(vec3 o, vec3 d) : o(o), d(d) {} 
};

enum Refl_t { DIFF, SPEC, REFR };  // material types, used in radiance()

struct Sphere {
	double rad;       // radius
	vec3 p, e, c;      // position, emission, color
	Refl_t refl;      // reflection type (DIFFuse, SPECular, REFRactive)
	Sphere(double rad, vec3 p, vec3 e, vec3 c, Refl_t refl) :
		rad(rad), p(p), e(e), c(c), refl(refl) {}
	double intersect(const Ray &r) const
	{
		// returns distance, negative if nohit
		// |r.o+t*r.d-p|^2 = rad^2
		// |r.d|^2*t^2 - 2*<r.d,p-r.o> + |p-r.o|^2 - rad^2 = 0
		// A = |r.d|^2 = 1 (r.d normalized), B = <r.d,p-r.o>, C = |p-r.o|^2 - rad^2, det = B^2 - C
		// t1,2 = B -+ sqrt(det) (as long as det>=0)
		vec3 op = p - r.o;
		// moved eps to the scene intersect routine, and modified sphere intersect to return negative if no intersection
		double t;
		double b = dot(op, r.d);
		double det = b * b - dot(op, op) + rad * rad;
		if (det < 0) return -1; else det = sqrt(det);
		// if t<0 -> no intersection
		return (t = b - det) > 0 ? t : b + det;
	}

	vec3 normal(const vec3& v) const
	{
		return (v - p) / rad;
	}
};


Sphere spheres[] = {//Scene: radius, position, emission, color, material
  Sphere(1e5, vec3(1e5 + 1,40.8,81.6), vec3(0),vec3(.75,.25,.25),DIFF),//Left
  Sphere(1e5, vec3(-1e5 + 99,40.8,81.6),vec3(0),vec3(.25,.25,.75),DIFF),//Rght
  Sphere(1e5, vec3(50,40.8, 1e5),     vec3(0),vec3(.75,.75,.75),DIFF),//Back
  Sphere(1e5, vec3(50,40.8,-1e5 + 170), vec3(0),vec3(0),           DIFF),//Frnt
  Sphere(1e5, vec3(50, 1e5, 81.6),    vec3(0),vec3(.75,.75,.75),DIFF),//Botm
  Sphere(1e5, vec3(50,-1e5 + 81.6,81.6),vec3(0),vec3(.75,.75,.75),DIFF),//Top
  Sphere(16.5,vec3(27,16.5,47),       vec3(0),vec3(1,1,1)*.999, SPEC),//Mirr
  Sphere(16.5,vec3(73,16.5,78),       vec3(0),vec3(1,1,1)*.999, REFR),//Glas
  Sphere(600, vec3(50,681.6 - .27,81.6),vec3(12,12,12),  vec3(0), DIFF) //Lite
};

inline double clamp(double x) { return x < 0 ? 0 : x>1 ? 1 : x; }

inline int toInt(double x) { return int(pow(clamp(x), 1 / 2.2) * 255 + .5); }


inline bool intersect(const Ray &r, double &t, int &id) {
	double n = sizeof(spheres) / sizeof(Sphere);
	double d;
	double inf = t = std::numeric_limits<double>::infinity();
	double eps = 1e-4;
	for (int i = int(n); i--;) 
		if ((d = spheres[i].intersect(r))>eps && d < t) 
		{ 
			t = d; 
			id = i; 
		}
	return t < inf;
}


/////////////////////////////////////////////////////////////////////////////////////
//									RAYTRACE									   //
/////////////////////////////////////////////////////////////////////////////////////
// https://en.wikipedia.org/wiki/Grayscale#Luma_coding_in_video_systems
double luma(const vec3& color)
{
	return dot(color, vec3(0.2126f, 0.7152f, 0.0722f));
}
vec3 radiance(const Ray &r, int depth) {
	// Limit max depth (or you'll run into a stackoverflow on some scenes)
	if (depth > MAX_DEPTH) return vec3(0);


	double t;                                 // distance to intersection
	int id = 0;                               // id of intersected object
	if (!intersect(r, t, id)) return vec3(0); // if miss, return black


	const Sphere &obj = spheres[id];        // the hit object
	// intersection point
	vec3 x = r.o + r.d*t;
	// changed to use a generic method to compute the normal
	vec3 n = obj.normal(x);
	// find the correct facing normal (if the ray is on the 'inside' when hitting obj, then flip it)
	vec3 nl = dot(n, r.d) < 0 ? n : n * -1;
	// albedo
	vec3 albedo = obj.c;



	// Russian Roulette:
	// probability to continue the ray (the less reflective the material, the lower)
	// modified to use luma rather than max(r, max(g, b))
	double russianRouletteProb = luma(albedo); 
	// Apply Russian Roulette only after the 5th bounce, and increment depth by 1
	if (++depth > 5) 
		if (random() < russianRouletteProb) 
			albedo /= russianRouletteProb; // boost the ray to compesnate for the probability of terminating it
		else 
			return obj.e; // terminate the ray




	if (obj.refl == DIFF) {                  // Ideal DIFFUSE reflection
		// generate cosine weighted points on the upper hemisphere through inverse transform sampling:
		// pdf = cos(theta) * sin(theta) / pi
		// integrating out phi ([0,2pi]), and then integrating over theta yields the marginal CDF for theta: 
		// (1-cos^2(theta))/2 = r2 -> cos(theta) = sqrt(1-r2) ~ sqrt(r2)
		// the last equivalence follows from 1-r2, and r2 being identically distributed
		// since the joint pdf is separable we don't need to use the conditional distribution and can integrate out theta ([0,pi/2])
		// and then integrate over phi: phi / (2*pi) = r1 -> phi = 2 * pi * r1
		double phi = 2 * M_PI*random();
		double r2 = random();
		double sinTheta = sqrt(r2);
		double cosTheta = sqrt(1 - r2);

		// use correct facing normal for the central vector of the hemisphere
		vec3 w = nl;
		// build an orthonormal basis from it
		vec3 u = normalize(cross(fabs(w.x) > .1 ? vec3(0, 1, 0) : vec3(1, 0, 0), w));
		vec3 v = cross(w, u);
		// transform the generated points (in spherical coordinates) into the orthonormal basis defined by the correct facing normal
		// use the sampled point as a direction for the bounce
		vec3 d = normalize(u*cos(phi)*sinTheta + v * sin(phi)*sinTheta + w * cosTheta);

		// From the estimator of the rendering equation (integrand / pdf):
		// radiance += emitted + brdf * radiance * cos(theta)/cos(theta) * PI
		// note that the cosine term in the integrand gets canceled by the cosine weighted pdf from which we sample the random direction
		// we assume that PI was already encoded in the brdf
		return obj.e + albedo*radiance(Ray(x, d), depth);
	}
	else if (obj.refl == SPEC)            // Ideal SPECULAR reflection
	{
		// reflection around the normal:
		//v___| n
		// \  |
		//  \ | -dot(d,n)*n
		// d \| 

		// v ___|___ v
		//   \  |  /
		//  d \ | / r
		//     \|/
		//

		// v = d - dot(d,n)*n
		// r = -d + 2*v = -d + 2*d -2*dot(d,n)*n = d - 2*dot(d,n)*n
		// r = r.d - 2.0*dot(r.d,n)*n;

		// reflect the ray around the normal
		vec3 refl = r.d - n * 2 * dot(n, r.d);

		// rendering equation estimator (integrand/pdf) for mirror brdf
		// mirror brdf = albedo * delta(theta_refl-theta_in) * delta(phi_refl - phi_in +- pi)  / (cos(theta_in) * sin(theta_in))
		return obj.e + albedo * radiance(Ray(x, refl), depth);
	}
	
	// refraction in the plane defined by d and n:
	// sin(theta1)*eta1 = sin(theta2)*eta2
	// sin(theta2) = sin(theta1) * eta1 / eta2
	// if eta1/eta2>1, it is possible that |sin(theta1)*eta1/eta2|>1
	// there is no angle theta1 that satisifies this -> total internal reflection
	// otherwise:
	//\  |   eta1
	// \ |      
	//d \| n          theta1 = angle(-d,n)
	//---------
	//   ||   eta2    theta2 = angle(r,-n)
	//   | |
	//-n |  | r

	// r = cos(theta2)*(-n) + sin(theta2)*perp(n)
	// perp(n) = d - dot(d,n)*n / |d - dot(d,n)*n| = (d - cos(theta1)*n)/sin(theta1)
	// cos(theta2) = sqrt(1-sin^2(theta2)) = sqrt(1-(eta1/eta2*sin(theta1))^2)
	// sin(theta2) = eta1/eta2*sin(theta1)
	// r = cos(theta2)*(-n) + eta1/eta2*sin(theta1)/sin(theta1)* (d - cos(theta1)*n)
	// r = -sqrt(1-(eta1/eta2*sin(theta1))^2)*n + eta1/eta2*(d-cos(theta1)*n)
	// r = eta1/eta2*d + (eta1/eta2*dot(d,n)-sqrt(1-(eta1/eta2)^2*(1-dot(d,n)^2)))*n


	Ray reflRay(x, r.d - n * 2 * dot(n,r.d));     // Ideal dielectric REFRACTION
	bool into = dot(n,nl) > 0;                    // Ray from outside going in?
	// indices of refraction
	double nc = 1;
	double nt = 1.5;
	// indices of refraction ratio
	double nnt = into ? nc / nt : nt / nc;
	// cosTheta
	double cosTheta = dot(r.d, nl);
	double cosTheta2Sqr;

	if ((cosTheta2Sqr = 1 - nnt * nnt*(1 - cosTheta * cosTheta)) < 0)    // Total internal reflection
		return obj.e + albedo*radiance(reflRay, depth);

	// refracted ray direction
	vec3 tdir = normalize(r.d*nnt - n * ((into ? 1 : -1)*(cosTheta*cosTheta + sqrt(cosTheta2Sqr))));

	// Schlick's Fresnel approximation:  Schlick, Christophe, An Inexpensive BDRF Model for Physically based Rendering,
	double a = nt - nc;
	double b = nt + nc;
	double R0 = a * a / (b*b);
	double cosTheta2 = dot(tdir, n);
	double c = 1 - (into ? -cosTheta : cosTheta2);
	double Re = R0 + (1 - R0)*c*c*c*c*c; // reflection weight
	double Tr = 1 - Re; // refraction weight

	// Russian roulette probability (for reflection)
	double P = .25 + .5*Re;
	// reflection weight boosted by the russian roulette probability
	double RP = Re / P;
	// refraction weight boosted by the russian roulette probability
	double TP = Tr / (1 - P);

	// splitting below depth 3
	if (depth < 3)
	{
		return obj.e + albedo*(radiance(reflRay, depth)*Re + radiance(Ray(x, tdir), depth)*Tr); // weighted relfection + refraction based on fresnel
	}
	else
	{
		// Russian roulette decision (between reflected and refracted ray)
		if (random() < P)
			return obj.e + albedo * radiance(reflRay, depth)*RP; //reflect
		else
			return obj.e + albedo * radiance(Ray(x, tdir), depth)*TP; //refract
	}
		
}
int main(int argc, char *argv[]) {
	int w = 1024;
	int h = 768;
	int samps = argc == 2 ? atoi(argv[1]) / 4 : 1; // # samples

	// camera position and forward vector (|forward| = 1)
	Ray cam(vec3(50, 52, 295.6), normalize(vec3(0, -0.042612, -1)));


	double aspectRatio = double(w) / double(h);
	// vertical field of view scaling:
	// tan(vfov/2) = fovScale / (2 * |forward|) -> fovScale = 2 * tan(vfov/2)
	double vfov = 0.502643; // ~ 58 degree vertical opening angle
	// axis scale in order to achieve the vfov above
	double fovScale = 2 * tan(0.5*vfov); // fovScale == 0.5135


	// horizontal camera axis, note that it is orthogonal to the forward vector by construction
	vec3 cx = vec3(aspectRatio, 0, 0)*fovScale;
	// vertical camera axis, orthogonal to both the horiztonal and forward axis due to the cross product
	vec3 cy = normalize(cross(cx, cam.d))*fovScale;


	vec3 r; //a temporary variable for the radiance, will be made private in the openMP loop

	// the image array
	vec3* c = new vec3[w*h];
	// constructor doesn't init vecs to 0 anymore, so memset the array to 0
	memset(c, 0, sizeof(vec3)*w*h);

// dynamic scheduling since work done for different rays can vary
#pragma omp parallel for schedule(dynamic, 1) private(r)       // OpenMP

	// Loop over the image rows
	for (int y = 0; y < h; ++y) 
	{    
		// print progress: samples and progress (percent)
		fprintf(stderr, "\rRendering (%d spp) %5.2f%%", samps * 4, 100.0*y / (h - 1));

		// Loop over the image columns
		for (unsigned short x = 0; x < w; ++x)
		{
			// 2x2 subpixel rows
			for (int sy = 0, i = (h - y - 1)*w + x; sy < 2; ++sy)
			{
				// 2x2 subpixel cols
				for (int sx = 0; sx < 2; sx++, r = vec3(0))
				{
					// samples per subpixel
					for (int s = 0; s < samps; s++)
					{
						// Bartlett window / tent function / triangular function random sampling
						// f(x) = 1 - |x|, |x|<=1, else f(x) = 0
						// Let F(x) be the indefinite integral of f(x): x in [-1,0], F(x) = x + x^2/2 + C; x in [0,1] F(x) = x - x^2/2 + C
						// We can use the inverse transform sampling method to sample from the tent pdf, we'll first split it into 2 pdfs:
						// in [-1,0] and [0,1], since both intervals have equal probability (1/2) 
						// Then for r1<=0.5 sample from the first, for r1>0.5 from the second

						// set r1' = 2 * r1

						// renormalize pdfs: f'(x) = 2 * f(x)
						// r1'<=1 -> 
						// -1<= x <= 0, F(x) - F(-1) = 2x + x^2 + 2 - 1 = (x+1)^2 = r1', then:
						//  x = sqrt(r1')-1

						// r1' > 1 ->
						// 0<=x<=1, F(x) - F(0) = 2x - x^2 = -(x^2 - 2x + 1) +1 = -(x-1)^2 + 1 = r1' - 1
						// x = 1 - sqrt(2-r1'), note that the negative root is taken since x<=1
						double r1 = 2 * random();
						double dx = r1 < 1 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);
						// the joint pdf is separable, so similarly for y:
						double r2 = 2 * random();
						double dy = r2 < 1 ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);

						// primary ray direction
						vec3 d = cx * (((sx + .5 + dx) / 2 + x) / w - .5) +
							cy * (((sy + .5 + dy) / 2 + y) / h - .5) + cam.d;
						
						// start the ray from the inside of the cornell 'box'
						// and accumulate the radiance per subpixel
						r = r + radiance(Ray(cam.o + d * 140, normalize(d)), 0)*(1. / samps);
					}
					// accumulate subpixels
					c[i] = c[i] + vec3(clamp(r.x), clamp(r.y), clamp(r.z))*.25;
				}
			}
		}
	}
	FILE *f = fopen("image.ppm", "w");         // Write image to PPM file.
	fprintf(f, "P3\n%d %d\n%d\n", w, h, 255);
	for (int i = 0; i < w*h; i++)
		fprintf(f, "%d %d %d ", toInt(c[i].x), toInt(c[i].y), toInt(c[i].z));
}
