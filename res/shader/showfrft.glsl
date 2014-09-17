
#version 330 core

#include "fullscreen.glsl"

uniform sampler3D sampler_frft;
uniform float exposure;
uniform float tz;

#ifdef _FRAGMENT_

out vec4 frag_color;

vec3 hdr(vec3 e) {
	return vec3(1.0) - exp(-exposure * e);
}

void main() {
	
	float a = texture(sampler_frft, vec3(texCoord, tz)).r;

	float p = a * a;

	// display hdr'd power spectrum
	frag_color = vec4(hdr(vec3(p)), 1.0);

}

#endif
