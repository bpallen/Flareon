
#version 330 core

#shader vertex
#shader geometry
#shader fragment

#define PI 3.141592654

uniform uint sides;
uniform float radius;



#ifdef _VERTEX_

void main() { }

#endif

#ifdef _GEOMETRY_

layout(points) in;
layout(triangle_strip, max_vertices = 100) out;

void main() {
	
	// this doesnt do any aspect ratio correction
	
	vec4 points[20];
	
	float theta = 0.093;
	for (uint i = 0u; i < sides; i++) {
		float x = radius * cos(theta);
		float y = radius * sin(theta);
		points[i] = vec4(x, y, 0.0, 1.0);
		theta += (2.0 * PI) / float(sides);
	}
	
	for (uint i = 0u; i < sides; i++) {
		gl_Position = vec4(vec3(0.0), 1.0);
		EmitVertex();
		gl_Position = points[i];
		EmitVertex();
		gl_Position = points[(i + 1u) % (sides)];
		EmitVertex();
		EndPrimitive();
	}
	
}

#endif

#ifdef _FRAGMENT_

out vec4 frag_color;

void main() {
	frag_color = vec4(1.0);
}

#endif
