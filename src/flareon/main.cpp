

#include <cstdio>
#include <cstdlib>
#include <cassert>

#include <iostream>
#include <chrono>
#include <string>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <thread>
#include <memory>
#include <complex>

#include <GFFT/GFFT.h>

#include "Initial3D.hpp"
#include "Flareon.hpp"
#include "Log.hpp"
#include "Window.hpp"
#include "Shader.hpp"
#include "Chrono.hpp"

using namespace std;
using namespace ambition;
using namespace initial3d;

ShaderManager *shaderman;

vec3d light_norm = vec3d::k(-1);

double exposure = 1.0;

// uniform buffer objects for lens system and precomputed bounce sequences
GLuint ubo_lens = 0;
GLuint ubo_bounce = 0;

// RGB wavelength sensitivities
GLuint tex_rgb = 0;
double wavelen_min = 9001, wavelen_max = 0;

// framebuffer and textures for deferred hdr
GLuint fbo_hdr = 0;
GLuint tex_hdr = 0;

// framebuffer and textures for aperture texture synthesis
GLuint fbo_ap = 0;
GLuint tex_ap = 0;
GLuint tex_ap_fft = 0;
GLuint tex_ap_frft = 0;
GLuint tex_star = 0;
GLuint tex_ghost = 0;

// wavelength to optimize antireflective coating for
double wavelen_ar = 550e-9;

// rgb weights and rendering wavelengths
float wavelengths[] = {
	1, 0, 0, 610e-9,
	0, 1, 0, 555e-9,
	0, 0, 1, 465e-9
};

unsigned num_wavelengths = 3;

struct lens_interface {
	// sphere radius
	// +ve => front-convex, -ve => front-concave
	// 0 => this isnt a lens
	double sr;
	// aperture radius
	// 0 => sensor plane
	double ar;
	// z-position of this surface on optical axis
	double z;
	// thickness of antireflective coating
	double d1;
	// refractive indices
	// n0 = refractive index of front material
	// n1 = refractive index of anti-reflective coating
	// n2 = refractive index of back material
	double n0, n1, n2;
};

// lens interfaces, in order from entry to sensor
vector<lens_interface> interfaces;
unsigned num_ghosts = 0;


void draw_fullscreen() {
	static GLuint vao = 0;
	if (vao == 0) {
		glGenVertexArrays(1, &vao);
	}
	glBindVertexArray(vao);
	glDrawArrays(GL_POINTS, 0, 1);
	glBindVertexArray(0);
}

void draw_aperture() {
	GLuint prog_ap = shaderman->getProgram("aperture.glsl");
	glUseProgram(prog_ap);
	glUniform1ui(glGetUniformLocation(prog_ap, "sides"), 8);
	glUniform1f(glGetUniformLocation(prog_ap, "radius"), 0.5);
	draw_fullscreen();
}

void draw_fourier(GLuint tex) {

	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, tex);

	GLuint prog = shaderman->getProgram("showfourier.glsl");
	glUseProgram(prog);

	glUniform1f(glGetUniformLocation(prog, "exposure"), exposure);
	glUniform1i(glGetUniformLocation(prog, "sampler_fourier"), 0);

	draw_fullscreen();

}

void draw_starburst() {

	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, tex_ap_fft);
	glActiveTexture(GL_TEXTURE1);
	glBindTexture(GL_TEXTURE_1D, tex_rgb);

	GLuint prog = shaderman->getProgram("starburst.glsl");
	glUseProgram(prog);

	glUniform1f(glGetUniformLocation(prog, "exposure"), exposure);
	glUniform1i(glGetUniformLocation(prog, "sampler_fourier"), 0);
	glUniform1i(glGetUniformLocation(prog, "sampler_rgb"), 1);
	glUniform1f(glGetUniformLocation(prog, "wavelen_min"), wavelen_min);
	glUniform1f(glGetUniformLocation(prog, "wavelen_max"), wavelen_max);

	draw_fullscreen();

}

// in-place transpose (square)
// input data is interleaved real and imaginary parts (ie 2 elements per sample)
void transpose(unsigned size, double *data) {
	for (unsigned i = 0; i < size; i++) {
		for (unsigned j = i + 1; j < size; j++) {
			unsigned k0 = 2 * size * i + 2 * j;
			unsigned k1 = 2 * size * j + 2 * i;
			swap(data[k0 + 0], data[k1 + 0]);
			swap(data[k0 + 1], data[k1 + 1]);
		}
	}
}

// in-place 1D FFT on multiple datasets (residing sequentially in memory)
// input data is interleaved real and imaginary parts (ie 2 elements per sample)
// result will need 'fftshifting'
void fft(unsigned size, unsigned count, double *data) {
	using namespace gfft;

	assert((size & (size - 1)) == 0 && "FFT: size must be a non-zero power of two");

	// initialization of the object factory
	Loki::Factory<AbstractFFT<double>,unsigned int> gfft_factory;
	FactoryInit<GFFTList<GFFT, 1, 27>::Result>::apply(gfft_factory);

	// power-of-two for data length
	unsigned p = 0;
	// assume size >= 1 (guaranteed by the above assertion)
	for (unsigned i = size >> 1; i > 0; ) {
		i >>= 1;
		p++;
	}

	// create an instance of the GFFT
	auto gfft = unique_ptr<AbstractFFT<double>>(gfft_factory.CreateObject(p));

	// run the FFTs
	for (unsigned i = 0; i < count; i++) {
		gfft->fft(data + 2 * size * i);
	}

}

// in-place 2D FFT (square)
// input data is interleaved real and imaginary parts (ie 2 elements per sample)
void fft2(unsigned size, double *data) {

	auto time0 = really_high_resolution_clock::now();

	// use seperability for 2D
	fft(size, size, data);
	transpose(size, data);
	fft(size, size, data);
	transpose(size, data);

	double dt = chrono::duration_cast<chrono::duration<double>>(really_high_resolution_clock::now() - time0).count();

	log("FFT2") << "size=" << size << ", took " << dt << "s";
}

void load_rgb_sensitivities() {

	if (tex_rgb) return;

	// parse file
	vector<float> rgb;
	ifstream ifs("./res/cones.txt");
	if (!ifs.good()) {
		throw runtime_error("unable to open file 'cones.txt'");
	}
	while (ifs.good()) {
		string line;
		getline(ifs, line);
		istringstream iss(line);
		string tok0;
		iss >> tok0;
		if (iss.fail()) continue;
		if (tok0[0] == '#') continue;
		iss.seekg(0);
		double wl, r, g, b;
		iss >> wl >> r >> g >> b;
		if (iss.fail()) {
			log("RGB").warning() << "Encountered invalid line '" << line << "'";
		} else {
			wl *= 1e-9;
			r = math::pow(10.0, r);
			g = math::pow(10.0, g);
			b = math::pow(10.0, b);
			wavelen_min = math::min(wavelen_min, wl);
			wavelen_max = math::max(wavelen_max, wl);
			// lines are assumed to be in increasing order of equally-spaced-wavelengths
			rgb.push_back(r);
			rgb.push_back(g);
			rgb.push_back(b);
			rgb.push_back(0);
		}
	}

	log("RGB") << "Read " << (rgb.size() / 4) << " sensitivity entries";
	log("RGB") << "Min wavelength: " << wavelen_min;
	log("RGB") << "Max wavelength: " << wavelen_max;

	glGenTextures(1, &tex_rgb);
	glBindTexture(GL_TEXTURE_1D, tex_rgb);
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexImage1D(GL_TEXTURE_1D, 0, GL_RGBA16F, rgb.size() / 4, 0, GL_RGBA, GL_FLOAT, &rgb[0]);

}

void make_textures() {
	
	static const unsigned tex_size = 1024;
	
	load_rgb_sensitivities();

	if (!fbo_ap) glGenFramebuffers(1, &fbo_ap);
	glBindFramebuffer(GL_FRAMEBUFFER, fbo_ap);

	if (!tex_ap) {
		glGenTextures(1, &tex_ap);
		glBindTexture(GL_TEXTURE_2D, tex_ap);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, tex_size, tex_size, 0, GL_RGBA, GL_FLOAT, nullptr);
	}
	
	if (!tex_ap_fft) {
		glGenTextures(1, &tex_ap_fft);
		glBindTexture(GL_TEXTURE_2D, tex_ap_fft);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	}
	
	if (!tex_ap_frft) {
		glGenTextures(1, &tex_ap_frft);
		glBindTexture(GL_TEXTURE_2D, tex_ap_frft);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	}
	
	if (!tex_star) {
		glGenTextures(1, &tex_star);
		glBindTexture(GL_TEXTURE_2D, tex_star);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
		// TODO 32f for the moment
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, tex_size, tex_size, 0, GL_RGBA, GL_FLOAT, nullptr);
	}
	
	if (!tex_ghost) {
		glGenTextures(1, &tex_ghost);
		glBindTexture(GL_TEXTURE_2D, tex_ghost);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16F, tex_size, tex_size, 0, GL_RGBA, GL_FLOAT, nullptr);
	}
	
	// draw the aperture transmission function
	glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, tex_ap, 0);
	glDrawBuffer(GL_COLOR_ATTACHMENT0);
	glViewport(0, 0, tex_size, tex_size);
	glClear(GL_COLOR_BUFFER_BIT);
	draw_aperture();
	checkGL();
	
	// read back aperture
	vector<float> ap_raw(tex_size * tex_size);
	glReadBuffer(GL_COLOR_ATTACHMENT0);
	glReadPixels(0, 0, tex_size, tex_size, GL_RED, GL_FLOAT, &ap_raw[0]);
	
	// unpack aperture for FFT
	vector<double> ap_temp(2 * tex_size * tex_size);
	for (unsigned i = 0; i < ap_raw.size(); i++) {
		// flattened complex numbers
		ap_temp[2 * i + 0] = ap_raw[i];
		ap_temp[2 * i + 1] = 0;
	}

	// FFT
	fft2(tex_size, &ap_temp[0]);
	
	// repack aperture FFT
	vector<float> ap_fft(2 * tex_size * tex_size);
	for (unsigned i = 0; i < ap_raw.size(); i++) {
		// flattened complex numbers - as amplitude and phase
		// this interpolates a lot better than real/imag
		auto z = complex<double>(ap_temp[2 * i + 0], ap_temp[2 * i + 1]);
		ap_fft[2 * i + 0] = abs(z);
		ap_fft[2 * i + 1] = arg(z);
	}

	// upload aperture FFT texture
	glBindTexture(GL_TEXTURE_2D, tex_ap_fft);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RG32F, tex_size, tex_size, 0, GL_RG, GL_FLOAT, &ap_fft[0]);
	glGenerateMipmap(GL_TEXTURE_2D);
	checkGL();
	
	// draw the starburst texture
	glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, tex_star, 0);
	glDrawBuffer(GL_COLOR_ATTACHMENT0);
	glViewport(0, 0, tex_size, tex_size);
	glClear(GL_COLOR_BUFFER_BIT);
	// TODO cleaner...
	exposure = 1.0;
	draw_starburst();
	checkGL();
	vector<float> star_temp(4 * tex_size * tex_size);
	glReadBuffer(GL_COLOR_ATTACHMENT0);
	glReadPixels(0, 0, tex_size, tex_size, GL_RGBA, GL_FLOAT, &star_temp[0]);
	
	// TODO save to exr for andrew
	
	// cleanup
	checkGL();
	glFinish();
	glUseProgram(0);
	glBindTexture(GL_TEXTURE_2D, 0);
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
	
}

void precompute(const string &path) {
	// parse the lens specification file
	ifstream ifs(path);
	if (!ifs.is_open()) {
		throw runtime_error("unable to open file '" + path + "'");
	}
	
	double n_last = 1.0;
	double dz_last = 0.0;
	while (ifs.good()) {
		string line;
		getline(ifs, line);
		istringstream iss(line);
		string t0;
		iss >> t0;
		if (iss.fail()) continue;
		if (t0[0] == '#') continue;
		iss.seekg(0);
		double sr, dz, n, ar;
		iss >> sr >> dz >> n >> ar;
		if (iss.fail()) {
			log("LensData").warning() << "Encountered invalid line '" << line << "' in file '" << path << "'";
		} else {
			// correct mm -> m
			sr /= 1000.0;
			dz /= 1000.0;
			ar /= 1000.0;
			// correct aperture diameter -> radius
			ar *= 0.5;
			// create interface
			lens_interface li;
			li.sr = sr;
			li.z = 0;
			li.n0 = n_last;
			li.n2 = n;
			li.ar = ar;
			// compute characteristics of antireflective coating
			li.n1 = math::max(math::sqrt(li.n0 * li.n2), 1.38);
			li.d1 = wavelen_ar / 4.0 / li.n1;
			// move absolute positions of already-loaded interfaces (because they're specified relatively, front to back)
			for (auto it = interfaces.begin(); it != interfaces.end(); it++) {
				it->z -= dz_last;
			}
			// update parser state
			n_last = n;
			dz_last = dz;
			// add to interfaces list
			interfaces.push_back(li);
		}
	}
	
	if (interfaces.size() == 0) {
		log("LensData").error() << "Lens data is empty";
		throw runtime_error("lens data is empty");
	}
	
	if (interfaces.back().ar > 0) {
		log("LensData").warning() << "Lens configuration has no sensor plane";
	}

	log("LensData") << interfaces.size() << " lens interfaces loaded";
	for (auto li : interfaces) {
		log("LensData") << "z=" << li.z;
	}
	
	// find aperture interface
	unsigned aperture_index = 0;
	for (unsigned i = 0; i < interfaces.size(); i++) {
		const lens_interface &li = interfaces[i];
		if (li.sr == 0.0 && li.ar > 0) {
			aperture_index = i;
		}
	}

	log("LensData") << "aperture at index " << aperture_index;

	// create uniform block for lens configuration
	// this relies on the std140 layout of the uniform block
	vector<GLuint> lens_block_words(4 + interfaces.size() * 8);
	lens_block_words[0] = interfaces.size();
	lens_block_words[1] = aperture_index;
	for (unsigned i = 0; i < interfaces.size(); i++) {
		unsigned j = 4 + i * 8;
		// sr, ar, z, d1
		reinterpret_cast<float &>(lens_block_words[j + 0]) = interfaces[i].sr;
		reinterpret_cast<float &>(lens_block_words[j + 1]) = interfaces[i].ar;
		reinterpret_cast<float &>(lens_block_words[j + 2]) = interfaces[i].z;
		reinterpret_cast<float &>(lens_block_words[j + 3]) = interfaces[i].d1;
		// n0, n1, n2
		reinterpret_cast<float &>(lens_block_words[j + 4]) = interfaces[i].n0;
		reinterpret_cast<float &>(lens_block_words[j + 5]) = interfaces[i].n1;
		reinterpret_cast<float &>(lens_block_words[j + 6]) = interfaces[i].n2;
	}

	// upload lens uniform buffer
	if (ubo_lens) glDeleteBuffers(1, &ubo_lens);
	glGenBuffers(1, &ubo_lens);
	glBindBuffer(GL_UNIFORM_BUFFER, ubo_lens);
	glBufferData(GL_UNIFORM_BUFFER, lens_block_words.size() * sizeof(GLuint), &lens_block_words[0], GL_STATIC_DRAW);

	// create uniform block for bounce enumeration (and enumerate bounces)
	// this relies on the std140 layout of the uniform block
	num_ghosts = 0;
	vector<GLuint> bounce_block_words;
	for (unsigned i = 1; i < interfaces.size(); i++) {
		const lens_interface &li0 = interfaces[i];
		if (math::abs(li0.sr) > 0 && li0.ar > 0) {
			for (unsigned j = 0; j < i; j++) {
				const lens_interface &li1 = interfaces[j];
				if (math::abs(li1.sr) > 0 && li1.ar > 0) {
					bounce_block_words.push_back(i);
					bounce_block_words.push_back(j);
					// need to pad to vec4 alignment
					bounce_block_words.push_back(0);
					bounce_block_words.push_back(0);
					num_ghosts++;
				}
			}
		}
	}

	log("LensData") << "enumerated " << num_ghosts << " ghosts";
	
	// upload bounce uniform buffer
	if (ubo_bounce) glDeleteBuffers(1, &ubo_bounce);
	glGenBuffers(1, &ubo_bounce);
	glBindBuffer(GL_UNIFORM_BUFFER, ubo_bounce);
	glBufferData(GL_UNIFORM_BUFFER, bounce_block_words.size() * sizeof(GLuint), &bounce_block_words[0], GL_STATIC_DRAW);

	// cleanup
	glBindBuffer(GL_UNIFORM_BUFFER, 0);
	
	// texture synthesis
	make_textures();
	
}

void init_fbo_hdr(const size2i &size) {
	
	if (fbo_hdr) glDeleteFramebuffers(1, &fbo_hdr);
	glGenFramebuffers(1, &fbo_hdr);
	
	glBindFramebuffer(GL_DRAW_FRAMEBUFFER, fbo_hdr);
	
	if (tex_hdr) glDeleteTextures(1, &tex_hdr);
	glGenTextures(1, &tex_hdr);
	
	glActiveTexture(GL_TEXTURE0);
	
	glBindTexture(GL_TEXTURE_2D, tex_hdr);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16F, size.w, size.h, 0, GL_RGBA, GL_FLOAT, nullptr);
	glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, tex_hdr, 0);
	checkGL();
	
	GLenum bufs[] = { GL_COLOR_ATTACHMENT0 };
	glDrawBuffers(1, bufs);
	
	checkFB();

	glBindTexture(GL_TEXTURE_2D, 0);
	
}

// width and height are in terms of on-screen triangles
template <unsigned Width, unsigned Height>
void draw_fullscreen_grid_adjacency_border_instanced(unsigned instances) {
	static_assert(Width >= 1 && Height >= 1, "grid must be at least 1x1");
	static GLuint vao = 0;
	if (vao == 0) {
		glGenVertexArrays(1, &vao);
		glBindVertexArray(vao);
		GLuint vbo, ibo;
		glGenBuffers(1, &vbo);
		glGenBuffers(1, &ibo);
		glBindBuffer(GL_ARRAY_BUFFER, vbo); // this doesnt stick to the vao
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo); // this does
		// generate x and y components
		vector<float> x_vals, y_vals;
		double dx = 2.0 / Width, dy = 2.0 / Height;
		// bottom and left off-screen vertices
		x_vals.push_back(-1.0 - dx);
		y_vals.push_back(-1.0 - dy);
		// on-screen vertices
		for (unsigned i = 0; i < Width; i++) {
			x_vals.push_back(-1.0 + dx * i);
		}
		for (unsigned i = 0; i < Height; i++) {
			y_vals.push_back(-1.0 + dy * i);
		}
		// force the last on-screen values of x/y to be correct
		x_vals.push_back(1.0);
		y_vals.push_back(1.0);
		// top and right off-screen vertices
		x_vals.push_back(1.0 + dx);
		y_vals.push_back(1.0 + dy);
		// generate actual points, index 0 is reserved for 'not a vertex'
		vector<float> points { 0.0f, 0.0f };
		for (unsigned i = 0; i < x_vals.size(); i++) {
			for (unsigned j = 0; j < y_vals.size(); j++) {
				points.push_back(x_vals[i]);
				points.push_back(y_vals[j]);
			}
		}
		// generate triangle vertex indices
		vector<GLuint> indices;
		auto push_index = [&](int i, int j) {
			// index 0 is reserved for 'not a vertex'
			int ii = 0;
			if (i >= 0 && i <= Width + 2 && j >= 0 && j <= Height + 2) {
				ii = 1 + y_vals.size() * i + j;
			}
			indices.push_back(unsigned(ii));
		};
		for (unsigned i = 0; i < Width + 2; i++) {
			for (unsigned j = 0; j < Height + 2; j++) {
				// first tri (:), p1=(i,j)
				// 6---5---4 //
				//   \ |:\ | //
				//     1---3 //
				//       \ | //
				//         2 //
				push_index(i, j);
				push_index((i + 1), j - 1);
				push_index((i + 1), j);
				push_index((i + 1), j + 1);
				push_index(i, j + 1);
				push_index((i - 1), j + 1);
				// second tri (:), p6=(i,j)
				// 4         //
				// | \       //
				// 5---3     //
				// | \:| \   //
				// 6---1---2 //
				push_index((i + 1), j);
				push_index((i + 2), j);
				push_index((i + 1), j + 1);
				push_index(i, j + 2);
				push_index(i, j + 1);
				push_index(i, j);
			}
		}
		// upload data
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(GLuint), &indices[0], GL_STATIC_DRAW);
		glBufferData(GL_ARRAY_BUFFER, points.size() * sizeof(float), &points[0], GL_STATIC_DRAW);
		glEnableVertexAttribArray(0);
		glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 0, nullptr);
		glBindBuffer(GL_ARRAY_BUFFER, 0); // cleanup
	}
	glBindVertexArray(vao);
	glDrawElementsInstanced(GL_TRIANGLES_ADJACENCY, (Width + 2) * (Height + 2) * 12, GL_UNSIGNED_INT, nullptr, instances);
	glBindVertexArray(0);
}

void display(const size2i &size) {
	if (size.w < 1 || size.h < 1) return;

	// height of the sensor
	// temp - set to fraction of diameter of first interface's aperture
	double sensor_h = interfaces[0].ar * 2;

	// orthographic projection, just to scale things properly
	mat4d proj = mat4d::scale(double(size.h) / double(size.w), 1.0, 1.0) * mat4d::scale(2.0 / sensor_h);

	// composite as HDR
	glBindFramebuffer(GL_DRAW_FRAMEBUFFER, fbo_hdr);
	glViewport(0, 0, size.w, size.h);
	glClear(GL_COLOR_BUFFER_BIT);

	GLuint prog_flare = shaderman->getProgram("flare.glsl");

	glUseProgram(prog_flare);
	
	// bind uniform buffer objects
	glBindBufferBase(GL_UNIFORM_BUFFER, 0, ubo_lens);
	glBindBufferBase(GL_UNIFORM_BUFFER, 1, ubo_bounce);
	glUniformBlockBinding(prog_flare, glGetUniformBlockIndex(prog_flare, "InterfacesBlock"), 0);
	glUniformBlockBinding(prog_flare, glGetUniformBlockIndex(prog_flare, "BouncesBlock"), 1);

	// set other uniforms
	glUniformMatrix4fv(glGetUniformLocation(prog_flare, "proj_matrix"), 1, true, mat4f(proj));
	glUniform1f(glGetUniformLocation(prog_flare, "lens_scale"), interfaces[0].ar);
	glUniform3fv(glGetUniformLocation(prog_flare, "light_norm"), 1, vec3f(light_norm));
	glUniform1ui(glGetUniformLocation(prog_flare, "num_wavelengths"), num_wavelengths);
	glUniform4fv(glGetUniformLocation(prog_flare, "wavelengths"), num_wavelengths, wavelengths);

	checkGL();
	
	// temp
	//glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	
	// instance the ghosts (batched by complexity == grid resolution)
	// blend ghosts together
	glDisable(GL_CULL_FACE);
	glDisable(GL_DEPTH_TEST);
	glEnable(GL_BLEND);
	glBlendEquation(GL_FUNC_ADD);
	glBlendFunc(GL_ONE, GL_ONE);
	glUniform1ui(glGetUniformLocation(prog_flare, "num_quads"), 64 * 64);
	//draw_fullscreen_grid_adjacency_border_instanced<64, 64>(num_ghosts * num_wavelengths);
	glDisable(GL_BLEND);
	checkGL();
	
	glFinish();
	
	// apply HDR tonemap
	glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
	glViewport(0, 0, size.w, size.h);
	
	GLuint prog_hdr = shaderman->getProgram("hdr.glsl");
	
	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, tex_hdr);
	
	glUseProgram(prog_hdr);
	
	glUniform1i(glGetUniformLocation(prog_hdr, "sampler_hdr"), 0);
	glUniform1f(glGetUniformLocation(prog_hdr, "exposure"), exposure);
	
	//draw_fullscreen();
	
	draw_starburst();
	//draw_aperture();
	//draw_fourier(tex_ap_fft);

	checkGL();

	glFinish();

}

int main(int argc, char *argv[]) {
	
	if (argc < 2) {
		log().error() << "Usage: flareon <lensdata>";
		exit(1);
	}
	
	Window *win = createWindow().size(1024, 1024).title("I Choose You, Flareon!").visible(true).debug(true);
	win->makeContextCurrent();

	shaderman = new ShaderManager("./res/shader");
	
	win->onResize.attach([](const window_size_event &e) {
		if (e.size.w > 0 && e.size.h > 0) {
			init_fbo_hdr(e.size);
		}
		return false;
	});

	win->onKeyPress.attach([](const key_event &e) {
		double rot_angle = math::pi() / 180;
		quatd rot = quatd::one();
		if (e.key == GLFW_KEY_UP) rot = quatd::axisangle(vec3d::i(), rot_angle);
		if (e.key == GLFW_KEY_DOWN) rot = quatd::axisangle(vec3d::i(), -rot_angle);
		if (e.key == GLFW_KEY_RIGHT) rot = quatd::axisangle(vec3d::j(), -rot_angle);
		if (e.key == GLFW_KEY_LEFT) rot = quatd::axisangle(vec3d::j(), rot_angle);
		light_norm = ~(rot * light_norm);

		if (e.key == GLFW_KEY_EQUAL) exposure *= 1.2;
		if (e.key == GLFW_KEY_MINUS) exposure /= 1.2;

		return false;
	});

	// lens-specific precomputation
	precompute(argv[1]);
	
	auto time_fps = chrono::steady_clock::now();
	unsigned fps = 0;
	
	init_fbo_hdr(win->size());
	
	while (!win->shouldClose()) {
		auto now = chrono::steady_clock::now();
		if (now - time_fps > chrono::seconds(1)) {
			time_fps = now;
			char cbuf[100];
			sprintf(cbuf, "Flareon [%d FPS @%dx%d]", fps, win->width(), win->height());
			win->title(cbuf);
			fps = 0;
		}

		glfwPollEvents();

		display(win->size());

		glFinish();
		win->swapBuffers();
		fps++;
	}

	delete win;
	glfwTerminate();

	return 0;
}






























