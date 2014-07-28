//
// Flareon: lens lare shader program
//

#version 330 core

#shader vertex
#shader geometry
#shader fragment

// vertex shader
#ifdef _VERTEX_

layout(location = 0) in vec2 pos_p;

out VertexData {
	vec4 tex;
} vertex_out;

void main() {
	gl_Position = vec4(pos_p, 0, 1);
}

#endif

// geometry shader
#ifdef _GEOMETRY_

layout(triangles_adjacency) in;
layout(triangle_strip, max_vertices = 3) out;

in VertexData {
	vec4 tex;
} vertex_in[];

out VertexData {
	vec4 tex;
} vertex_out;

void main() {
	gl_Position = gl_in[0].gl_Position;
	EmitVertex();
	gl_Position = gl_in[2].gl_Position;
	EmitVertex();
	gl_Position = gl_in[4].gl_Position;
	EmitVertex();
	EndPrimitive();
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