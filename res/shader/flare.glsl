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

uniform mat4 proj_matrix;
uniform float lens_scale;

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
	// (aperture tex coord, r_rel, intensity)
	vec4 tex;
};

struct Intersection {
	vec3 pos;
	vec3 norm;
	float theta;
	bool hit;
	bool inverted;
};

// lens configuration
layout(std140) uniform InterfacesBlock {
	uint num_interfaces;
	LensInterface interfaces[MAX_INTERFACES];
};

// the interfaces to bounce off
layout(std140) uniform BouncesBlock {
	uvec2 bounces[MAX_GHOSTS];
};

Intersection intersect_plane(Ray ray, LensInterface li) {
	Intersection isect;
	isect.pos = ray.pos + ray.dir * ((li.z - ray.pos.z) / ray.dir.z);
	isect.norm = vec3(0.0, 0.0, mix(1.0, -1.0, ray.dir.z > 0.0));
	isect.theta = 0.0; // meaningless
	isect.hit = true;
	isect.inverted = false;
	return isect;
}

Intersection intersect_sphere(Ray ray, LensInterface li) {
	Intersection isect;
	vec3 centre = vec3(0.0, 0.0, li.z - li.sr);
	vec3 d = ray.pos - centre;
	float b = dot(d, ray.dir);
	float c = dot(d, d) - li.sr * li.sr;
	float b2_c = b * b - c;
	if (b2_c < 0.0) {
		// no intersection
		isect.hit = false;
		return isect;
	}

	float sgn = mix(-1.0, 1.0, li.sr * ray.dir.z > 0.0);
	float t = sqrt(b2_c) * sgn - b;
	isect.pos = ray.dir * t + ray.pos;
	isect.norm = normalize(isect.pos - centre);
	isect.norm *= mix(1.0, -1.0, dot(isect.norm, ray.dir) > 0.0);
	isect.theta = acos(dot(-ray.dir, isect.norm));
	isect.hit = true;
	isect.inverted = t < 0.0;

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
	float rel_phase = 2.0 * PI / wavelen * (delay - dx * sin(theta0));

	// add up sines of different phase and amplitude
	float out_s2 = rs01 * rs01 + ris * ris + 2.0 * rs01 * ris * cos(rel_phase);
	float out_p2 = rp01 * rp01 + rip * rip + 2.0 * rp01 * rip * cos(rel_phase);

	return (out_s2 + out_p2) * 0.5;
}

// main ray tracing function
Ray trace(uint gid, Ray ray, float wavelen) {
	
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
		delta *= should_reflect ? -1 : 1;
		stage += uint(should_reflect);

		// test intersection, record max relative radius or aperture tex coord
		Intersection isect;
		if (li.sr > 0.0) {
			// spherical lens
			isect = intersect_sphere(ray, li);
			ray.tex.z = max(ray.tex.z, length(isect.pos.xy) / li.ar);
		} else {
			isect = intersect_plane(ray, li);
			// assume its the aperture (if radius > 0)
			ray.tex.xy = mix(ray.tex.xy, isect.pos.xy / li.ar, bvec2(li.ar > 0.0));
		}

		// exit on miss
		if (!isect.hit) break;

		// update ray direction and position
		ray.dir = normalize(isect.pos - ray.pos); // why do i need this?
		ray.dir *= mix(1.0, -1.0, isect.inverted); // correct inverted intersection
		ray.pos = isect.pos;

		// reflection / refraction
		if (li.sr > 0.0) {
			
			// swap order of refractive indices if ray going in 'reverse'
			vec3 n = mix(li.n.xyz, li.n.zyx, bvec3(ray.dir.z > 0.0));

			if (should_reflect) {
				// reflection with AR coating
				ray.dir = reflect(ray.dir, isect.norm);
				ray.tex.a *= fresnel_ar(isect.theta, wavelen, li.d1, n);
			} else {
				// refraction
				ray.dir = refract(ray.dir, isect.norm, n.x / n.z);
				// test for total internal reflection
				if (all(equal(ray.dir, vec3(0.0)))) break;
			}
		}
	}

	if (i < int(num_interfaces)) {
		// ignore early exits
		ray.tex.a = 0.0;
	}

	return ray;
}


// vertex shader
#ifdef _VERTEX_

layout(location = 0) in vec2 pos_p;

out VertexData {
	flat uint id;
	flat Ray ray;
} vertex_out;

void main() {
	vertex_out.id = uint(gl_VertexID);
	// id 0 is reserved for 'not a vertex', so early exit
	if (gl_VertexID == 0) return;
	// entrance ray
	Ray r0;
	r0.pos = vec3(lens_scale * pos_p, interfaces[0].z + 0.001);
	r0.dir = vec3(0.0, 0.0, -1.0);
	r0.tex = vec4(vec3(0.0), 1.0);
	// TODO wavelength? light direction?
	vertex_out.ray = trace(uint(gl_InstanceID), r0, 550e-9);
}

#endif

// geometry shader
#ifdef _GEOMETRY_

layout(triangles_adjacency) in;
layout(triangle_strip, max_vertices = 3) out;

in VertexData {
	flat uint id;
	flat Ray ray;
} vertex_in[];

out VertexData {
	noperspective vec4 tex;
} vertex_out;

void main() {
	Ray r0 = vertex_in[0].ray;
	Ray r2 = vertex_in[2].ray;
	Ray r4 = vertex_in[4].ray;
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
} vertex_in;

out vec4 frag_color;

void main() {
	if (vertex_in.tex.z > 1.0) discard;

	frag_color = vec4(0.4 * (vertex_in.tex.xy * 0.5 + 0.5), 0.0, 1.0);
	//frag_color = length(vertex_in.tex.xy) < 0.2 ? vec4(10.0 * vertex_in.tex.a, 0.0, 0.0, 1.0) : vec4(0.0);
}

#endif































