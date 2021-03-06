//
// Flareon: lens flare shader program
//
// Adaptation and implementation of http://resources.mpi-inf.mpg.de/lensflareRendering/
//

#version 330 core

#shader vertex
#shader geometry
#shader fragment

#define PI 3.141592654

// max number of interfaces in the lens configurations
#define MAX_INTERFACES 50

// max number of ghosts / bounce sequences
#define MAX_GHOSTS 200

// max number of wavelengths to render at
#define MAX_WAVELENGTHS 8

// x/y offsets for secondary rays (used for area factor calculation)
#define SECONDARY_RAY_OFFSET 0.0001

uniform mat4 proj_matrix;
uniform float lens_scale;
uniform vec3 light_norm;
uniform uint num_quads;

// an optical interface
// fields are layed out for 8-float size with std140
struct LensInterface {
	// sphere radius
	// +ve => front-convex, -ve => front-concave
	// 0 => this isnt a lens
	float sr;
	// aperture radius
	// 0 => sensor plane
	float ar;
	// z-position of this surface on optical axis
	float z;
	// thickness of antireflective coating
	float d1;
	// refractive indices
	// n.x = refractive index of front material
	// n.y = refractive index of anti-reflective coating
	// n.z = refractive index of back material
	vec3 n;
};

struct Ray {
	vec3 pos;
	vec3 dir;
	// entrance plane area : sensor plane area
	float area_factor;
	// (aperture tex coord, r_rel, irradiance factor)
	vec4 tex;
};

struct Intersection {
	vec3 pos;
	vec3 norm;
	float theta;
	bool hit;
	//bool inverted;
};

// lens configuration
layout(std140) uniform InterfacesBlock {
	uint num_interfaces;
	uint aperture_index;
	LensInterface interfaces[MAX_INTERFACES];
};

// the interfaces to bounce off
layout(std140) uniform BouncesBlock {
	uvec2 bounces[MAX_GHOSTS];
};

// wavelengths
uniform uint num_wavelengths;
// rgb: weights for RGB conversion, w: wavelength
// TODO use the sensitivity texture?
uniform vec4 wavelengths[MAX_WAVELENGTHS];

// min/max frft orders (texture boundaries)
uniform float frft_min;
uniform float frft_max;

// textures
uniform sampler3D sampler_frft;

// TODO starburst texture

Intersection intersect_plane(Ray ray, LensInterface li) {
	Intersection isect;
	isect.pos = ray.pos + ray.dir * ((li.z - ray.pos.z) / ray.dir.z);
	isect.norm = vec3(0.0, 0.0, mix(1.0, -1.0, ray.dir.z > 0.0));
	isect.theta = 0.0; // meaningless
	isect.hit = true;
	//isect.inverted = false;
	return isect;
}

Intersection intersect_sphere(Ray ray, LensInterface li) {
	Intersection isect;
	vec3 centre = vec3(0.0, 0.0, li.z + li.sr);
	vec3 d = ray.pos - centre;
	float b = dot(d, ray.dir);
	float c = dot(d, d) - li.sr * li.sr;
	float b2_c = b * b - c;
	if (b2_c < 0.0) {
		// no intersection
		isect.hit = false;
		return isect;
	}

	float sgn = mix(-1.0, 1.0, li.sr * ray.dir.z < 0.0);
	float t = sqrt(b2_c) * sgn - b;
	isect.pos = ray.dir * t + ray.pos;
	isect.norm = normalize(isect.pos - centre);
	isect.norm *= mix(1.0, -1.0, dot(isect.norm, ray.dir) > 0.0);
	isect.theta = acos(dot(-ray.dir, isect.norm));
	isect.hit = true;
	//isect.inverted = t < 0.0;

	return isect;
}

float fresnel_ar(float theta0, float wavelen, float d1, vec3 n) {
	// TODO - how does this work?
	
	// refraction angles in coating and 2nd medium
	float theta1 = asin(sin(theta0) * n.x / n.y);
	float theta2 = asin(sin(theta0) * n.x / n.z);

	// amplitude for outer reflection / transmission on topmost interface
	float rs01 = -sin(theta0 - theta1) / sin(theta0 + theta1);
	float rp01 = tan(theta0 - theta1) / tan(theta0 + theta1);
	float ts01 = 2.0 * sin(theta1) * cos(theta0) / sin(theta0 + theta1);
	float tp01 = ts01 * cos(theta0 - theta1);

	// amplitude for inner reflection
	float rs12 = -sin(theta1 - theta2) / sin(theta1 + theta2);
	float rp12 = tan(theta1 - theta2) / tan(theta1 + theta2);

	// after passing through first surface twice:
	// 2 transmissions and 1 reflection
	float ris = ts01 * ts01 * rs12;
	float rip = tp01 * tp01 * rp12;

	// phase difference between outer and inner reflections
	float dy = d1 * n.y;
	float dx = tan(theta1) * dy;
	float delay = sqrt(dx * dx + dy * dy);
	float rel_phase = 4.0 * PI / wavelen * (delay - dx * sin(theta0));

	// add up sines of different phase and amplitude
	float out_s2 = rs01 * rs01 + ris * ris + 2.0 * rs01 * ris * cos(rel_phase);
	float out_p2 = rp01 * rp01 + rip * rip + 2.0 * rp01 * rip * cos(rel_phase);

	return (out_s2 + out_p2) * 0.5;
}

// main ray tracing function
Ray trace(uint gid, Ray ray, float wavelen) {
	
	// trace 2 secondary rays so we can calculate the area factor
	Ray ray1 = ray;
	Ray ray2 = ray;
	ray1.pos.x += SECONDARY_RAY_OFFSET;
	ray2.pos.y += SECONDARY_RAY_OFFSET;

	// initial area of micro-beam (apart from the factor of 0.5)
	float a0 = length(cross(ray1.pos - ray.pos, ray2.pos - ray.pos));

	// take into account the lambertian cosine term on the entrance plane
	ray.tex.a *= abs(ray.dir.z);

	// bounce sequence
	ivec3 bounce = ivec3(ivec2(bounces[gid]), -1);
	// interface index delta
	int delta = 1;
	// ray-tracing stage
	uint stage = 0u;
	// interface index
	int i = 0;

	// iterate!
	for (; i < int(num_interfaces); i += delta) {
		
		// current interface
		LensInterface li = interfaces[i];

		// does this sequence reflect off this interface?
		bool should_reflect = i == bounce[stage];
		delta = should_reflect ? -delta : delta;
		stage += uint(should_reflect);

		// test intersection, record aperture tex coord
		Intersection isect, isect1, isect2;
		if (abs(li.sr) > 0.0) {
			// spherical lens
			isect = intersect_sphere(ray, li);
			isect1 = intersect_sphere(ray1, li);
			isect2 = intersect_sphere(ray2, li);
		} else {
			isect = intersect_plane(ray, li);
			isect1 = intersect_plane(ray1, li);
			isect2 = intersect_plane(ray2, li);
			// record coord if its the aperture, otherwise do nothing
			// if this happens more than once, only the last intersection is used
			ray.tex.xy = mix(ray.tex.xy, isect.pos.xy / li.ar, bvec2(i == int(aperture_index)));
		}

		// exit on miss
		if (!isect.hit) break;
		
		// record max relative radius
		ray.tex.z = max(ray.tex.z, mix(0.0, length(isect.pos.xy) / li.ar, li.ar > 0.0));

		// update ray direction and position
		//ray.dir = normalize(isect.pos - ray.pos); // why do i need this?
		//ray.dir *= mix(1.0, -1.0, isect.inverted); // correct inverted intersection

		ray.pos = isect.pos;
		ray1.pos = isect1.pos;
		ray2.pos = isect2.pos;

		// reflection / refraction
		if (abs(li.sr) > 0.0) {
			
			// swap order of refractive indices if ray going in 'reverse'
			vec3 n = mix(li.n.xyz, li.n.zyx, bvec3(ray.dir.z < 0.0));

			// TODO temporary hack for wavelength-dependent refractive index
			n.xz -= (wavelen - 400e-9) * 1.6e5;

			if (should_reflect) {
				// reflection with AR coating
				ray.dir = reflect(ray.dir, isect.norm);
				ray1.dir = reflect(ray1.dir, isect1.norm);
				ray2.dir = reflect(ray2.dir, isect2.norm);
				float a = fresnel_ar(isect.theta, wavelen, li.d1, n);
				a = mix(a, 0.0, isnan(a) || isinf(a));
				ray.tex.a *= a;
			} else {
				// refraction
				ray.dir = refract(ray.dir, isect.norm, n.x / n.z);
				ray1.dir = refract(ray1.dir, isect1.norm, n.x / n.z);
				ray2.dir = refract(ray2.dir, isect2.norm, n.x / n.z);
				// test for total internal reflection
				if (all(equal(ray.dir, vec3(0.0)))) break;
			}
		}
	}

	// final area of micro-beam
	float a1 = length(cross(ray1.pos - ray.pos, ray2.pos - ray.pos));

	// area factor
	// a1 can contain nans (untested TIR i guess)
	ray.area_factor = mix(a0 / a1, 0.0, isnan(a1));

	// ignore early exits
	ray.tex.a = mix(ray.tex.a, 0.0, i < int(num_interfaces));

	// ignore nans
	//ray.tex.a = mix(ray.tex.a, 0.0, isnan(ray.tex.a));

	return ray;
}


// vertex shader
#ifdef _VERTEX_

layout(location = 0) in vec2 pos_p;

out VertexData {
	flat uint id;
	flat uint wid;
	flat Ray ray;
} vertex_out;

void main() {
	vertex_out.id = uint(gl_VertexID);
	// id 0 is reserved for 'not a vertex', so early exit
	if (gl_VertexID == 0) return;
	// entrance ray
	Ray ray;
	ray.pos = vec3(lens_scale * pos_p, interfaces[0].z);
	ray.dir = -light_norm;
	ray.tex = vec4(vec3(0.0), 1.0);
	// get ghost id and wavelength id from instance id
	uint gid = uint(gl_InstanceID) / num_wavelengths;
	vertex_out.wid = uint(gl_InstanceID) % num_wavelengths;
	// trace!
	vertex_out.ray = trace(gid, ray, wavelengths[vertex_out.wid].w);
}

#endif

// geometry shader
#ifdef _GEOMETRY_

// probably dont need the adjacency info
layout(triangles_adjacency) in;
layout(triangle_strip, max_vertices = 3) out;

in VertexData {
	flat uint id;
	flat uint wid;
	flat Ray ray;
} vertex_in[];

out VertexData {
	noperspective vec4 tex;
	// flat is ok, this will be constant across any instance
	flat uint wid;
} vertex_out;

void main() {
	Ray r0 = vertex_in[0].ray;
	Ray r2 = vertex_in[2].ray;
	Ray r4 = vertex_in[4].ray;

	// area factor for real triangle
	float a0 = 4.0 * lens_scale * lens_scale / float(num_quads);
	float a1 = length(cross(r2.pos - r0.pos, r4.pos - r0.pos));
	float af = mix(a0 / a1, 0.0, isnan(a1));

	// average area factor from rays
	float aaf = (r0.area_factor + r2.area_factor + r4.area_factor) / 3.0;

	// speculative correction
	// TODO idk
	float maf = mix(9001.0, af, aaf > 2.0 * af);

	r0.tex.a *= min(r0.area_factor, maf);
	r2.tex.a *= min(r2.area_factor, maf);
	r4.tex.a *= min(r4.area_factor, maf);

	// set wavelength id once
	vertex_out.wid = vertex_in[0].wid;
	// set ray data
	vertex_out.tex = r0.tex;
	gl_Position = proj_matrix * vec4(r0.pos.xy, 0.0, 1.0);
	EmitVertex();
	vertex_out.tex = r2.tex;
	gl_Position = proj_matrix * vec4(r2.pos.xy, 0.0, 1.0);
	EmitVertex();
	vertex_out.tex = r4.tex;
	gl_Position = proj_matrix * vec4(r4.pos.xy, 0.0, 1.0);
	EmitVertex();
	EndPrimitive();
}

#endif

// fragment shader
#ifdef _FRAGMENT_

in VertexData {
	noperspective vec4 tex;
	flat uint wid;
} vertex_in;

out vec4 frag_color;

void main() {
	if (vertex_in.tex.z > 1.0) discard;
	
	//frag_color = vec4(0.2 * (vertex_in.tex.xy * 0.5 + 0.5), 0.0, 1.0);
	
	// TODO f-stop
	float a = 0.15 * wavelengths[vertex_in.wid].w / 400e-9 * (5.0 / 18.0);
	
	float frft = texture(sampler_frft,
		vec3(
			vertex_in.tex.xy * 0.5 + 0.5,
			(a - frft_min) / (frft_max - frft_min)
		)
	).r;
	
	frag_color = vec4(1000.0 * vertex_in.tex.a * wavelengths[vertex_in.wid].rgb * frft, 1.0);
	
	//frag_color = vec4(vertex_in.tex.xy, 0.0, 1.0);
	
	// TODO sensible irradiance ?
}

#endif































