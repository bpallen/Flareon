
#version 330 core

#include "fullscreen.glsl"

uniform sampler2D sampler_hdr;

uniform float exposure;

#ifdef _FRAGMENT_

out vec4 frag_color;

vec3 hdr(vec3 e) {
	return vec3(1.0) - exp(-exposure * e);
}

void main() {
	frag_color = vec4(hdr(texture(sampler_hdr, texCoord).rgb), 1.0);
}

#endif
