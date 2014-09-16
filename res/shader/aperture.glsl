
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

float rand(vec2 seed) {
	return fract(sin(dot(seed + 0.317, vec2(12.9898, 78.233))) * 43758.5453);
}

void main() {
	
	// this produces fairly random (though somewhat grid-aligned) smoothed circular points
	
	float spot = 1.0;
	
	vec2 fc0 = floor(gl_FragCoord.xy / spot);
	vec2 fc1 = fc0 * spot + spot * 0.5 + 0.5;
	
	float junk = rand(fc0);
	
	// frag_color = vec4(mix(1.0, 1.0 - exp(-(2.0 / spot) * pow(distance(fc1, gl_FragCoord.xy), 2.0)), junk > 0.9));
	
	float d = distance(fc1, gl_FragCoord.xy);
	
	frag_color = vec4(mix(1.0, 1.0 - 0.9 * exp(-(2.0 / spot) * pow(d / spot * 5.0, spot / 2.0)), junk > 0.87));
}

#endif

















