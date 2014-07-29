//
// Flareon: lens lare shader program
//

#version 330 core

#shader vertex
#shader geometry
#shader fragment

#define PI 3.141592654

// max number of interfaces in the lens configurations
#define MAX_INTERFACES 50

// max number of 2-bounce sequences (i.e. ghosts)
#define MAX_GHOSTS 100

struct LensInterface {
	// sphere radius
	// +ve => front-convex, -ve => front-concave
	// 0 => this isnt a lens
	float sr;
	// z-position of this surface on optical axis
	float z;
	// refractive indices
	// n.x = refractive index of front material
	// n.y = refractive index of anti-reflective coating
	// n.z = refractive index of back material
	vec3 n;
	// aperture radius
	// 0 => sensor plane
	float ap;
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

// TODO buffer objects instead?

// lens configuration
uniform LensInterface interfaces[MAX_INTERFACES];

// bounce sequences
uniform uvec2 bounces[MAX_GHOSTS];

// number of lens interfaces
uniform uint num_interfaces;

// number of ghosts / bounce sequences
uniform uint num_ghosts;

// main ray tracing function
Ray trace(uint gid, Ray r, float lambda) {
	
	return r;
}


// vertex shader
#ifdef _VERTEX_

layout(location = 0) in vec2 pos0;

out VertexData {
	flat int id;
	flat Ray rays[MAX_GHOSTS];
} vertex_out;

void main() {
	vertex_out.id = gl_VertexID;
	// id 0 is reserved for 'not a vertex', so early exit
	if (gl_VertexID == 0) return;
	// entrance ray
	Ray r0;
	r0.pos = vec3(pos0, interfaces[0].z + 0.001);
	r0.dir = vec3(0.0, 0.0, -1.0);
	r0.tex = vec4(0.0);
	// trace each ghost
	for (uint gid = 0u; gid < num_ghosts; gid++) {
		// TODO wavelength?
		vertex_out.rays[gid] = trace(gid, r0, 550e-9);
	}
}

#endif

// geometry shader
#ifdef _GEOMETRY_

layout(triangles_adjacency) in;

// TODO apparently i cant do this
layout(triangle_strip, max_vertices = 3 * MAX_GHOSTS) out;

in VertexData {
	flat int id;
	flat Ray rays[MAX_GHOSTS];
} vertex_in[];

out VertexData {
	vec4 tex;
} vertex_out;

void main() {
	for (uint gid = 0u; gid < num_ghosts; gid++) {
		Ray r0 = vertex_in[0].rays[gid];
		Ray r2 = vertex_in[2].rays[gid];
		Ray r4 = vertex_in[4].rays[gid];
		vertex_out.tex = r0.tex;
		gl_Position = vec4(r0.pos.xy, 0.0, 1.0);
		EmitVertex();
		vertex_out.tex = r1.tex;
		gl_Position = vec4(r1.pos.xy, 0.0, 1.0);
		EmitVertex();
		vertex_out.tex = r1.tex;
		gl_Position = vec4(r1.pos.xy, 0.0, 1.0);
		EmitVertex();
		EndPrimitive();
	}
}

#endif

// fragment shader
#ifdef _FRAGMENT_

in VertexData {
	vec4 tex;
} vertex_in;

out vec4 frag_color;

void main() {
	frag_color = vec4(1.0, 0.0, 0.0, 1.0);
}

#endif































