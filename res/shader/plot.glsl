
#version 330 core

#shader vertex
#shader geometry
#shader fragment

uniform sampler1D sampler_data;
uniform float plot_min, plot_max;

#ifdef _VERTEX_

flat out int id;

void main() {
	id = gl_InstanceID;
}

#endif

#ifdef _GEOMETRY_

layout(points) in;
layout(line_strip, max_vertices = 2) out;

flat in int id[];

out vec3 color;

void main() {
	
	int size = textureSize(sampler_data, 0);

	// generate line
	for (int i = id[0]; i < id[0] + 2; i++) {
		color = vec3(0.0, 0.0, 1.0);
		float v = texelFetch(sampler_data, i, 0);
		float x = (float(i) + 0.5) / float(size) * 2.0 - 1.0;
		float y = (v - plot_min) / (plot_max - plot_min) * 2.0 - 1.0;
		gl_Position = vec4(x, y, 0.0, 1.0);
		EmitVertex();
	}
	EndPrimitive();

}

#endif

#ifdef _FRAGMENT_

in vec3 color;

out vec4 frag_color;

void main() {
	frag_color = vec4(color, 1.0);
}

#endif