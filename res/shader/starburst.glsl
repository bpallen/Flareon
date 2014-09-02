
#version 330 core

#include "fullscreen.glsl"

uniform float exposure;

// FFT of the aperture
uniform sampler2D sampler_fourier;

// RGB sensitivities
uniform sampler1D sampler_rgb;

// min / max wavelengths for the sensitivity table
uniform float wavelen_min;
uniform float wavelen_max;

vec3 rgb(float wavelen) {
	float a = (wavelen - wavelen_min) / (wavelen_max - wavelen_min);
	return texture(sampler_rgb, a).rgb;
}

vec3 hdr(vec3 e) {
	return vec3(1.0) - exp(-exposure * e);
}

#ifdef _FRAGMENT_

out vec4 frag_color;

void main() {
	
	// wavelength for base fourier image
	const float w0 = 575e-9;

	const uint samples = 200u;

	vec3 f = vec3(0.0);

	vec3 rgb_t = vec3(0.0);

	for (uint i = 0u; i < samples; i++) {
		float wi = wavelen_min + (float(i) + 0.5) / float(samples) * (wavelen_max - wavelen_min);
		float scale = w0 / wi;
		float a = texture(sampler_fourier, clamp((texCoord - 0.5) * scale, -0.5, 0.5)).r;
		rgb_t += rgb(wi);
		f += rgb(wi) * a * a;
	}

	// f /= float(samples);
	f /= rgb_t;

	frag_color = vec4(f * exposure, 1.0);

}

#endif
