

#include <cstdio>
#include <cstdlib>

#include <iostream>
#include <chrono>
#include <string>
#include <fstream>
#include <sstream>
#include <stdexcept>

#include "Initial3D.hpp"
#include "Flareon.hpp"
#include "Log.hpp"
#include "Window.hpp"
#include "Shader.hpp"

using namespace std;
using namespace ambition;
using namespace initial3d;

ShaderManager *shaderman;

// uniform buffer objects for lens system and precomputed bounce sequences
GLuint ubo_lens = 0;
GLuint ubo_bounce = 0;

// wavelength to optimize antireflective coating for
double wavelen_ar = 500e-9;

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
				it->z += dz_last;
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
	
	// create uniform block for lens configuration
	// this relies on the std140 layout of the uniform block
	vector<GLuint> lens_block_words(4 + interfaces.size() * 8);
	lens_block_words[0] = interfaces.size();
	// TODO properly, hard coded aperture index for now
	lens_block_words[1] = 5;
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
	if (!ubo_lens) glDeleteBuffers(1, &ubo_lens);
	glGenBuffers(1, &ubo_lens);
	glBindBuffer(GL_UNIFORM_BUFFER, ubo_lens);
	glBufferData(GL_UNIFORM_BUFFER, lens_block_words.size() * sizeof(GLuint), &lens_block_words[0], GL_STATIC_DRAW);

	// create uniform block for bounce enumeration
	// TODO properly, hard coded single ghost for now
	vector<GLuint> bounce_block_words(4);
	bounce_block_words[0] = 7;
	bounce_block_words[1] = 0;
	num_ghosts = 1;

	// upload bounce uniform buffer
	if (!ubo_bounce) glDeleteBuffers(1, &ubo_bounce);
	glGenBuffers(1, &ubo_bounce);
	glBindBuffer(GL_UNIFORM_BUFFER, ubo_bounce);
	glBufferData(GL_UNIFORM_BUFFER, bounce_block_words.size() * sizeof(GLuint), &bounce_block_words[0], GL_STATIC_DRAW);

	// cleanup
	glBindBuffer(GL_UNIFORM_BUFFER, 0);
	
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

void display(const size2i &sz) {

	// height of the sensor
	// temp - set to fraction of diameter of first interface's aperture
	double sensor_h = interfaces[0].ar * 2;

	// orthographic projection, just to scale things properly
	mat4d proj = mat4d::scale(double(sz.h) / double(sz.w), 1.0, 1.0) * mat4d::scale(2.0 / sensor_h);

	glViewport(0, 0, sz.w, sz.h);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

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

	checkGL();
	
	// temp
	//glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	
	// instance the ghosts (batched by complexity == grid resolution)
	// do multiple wavelengths at once
	// blend ghosts together
	glDisable(GL_CULL_FACE);
	glDisable(GL_DEPTH_TEST);
	glEnable(GL_BLEND);
	glBlendEquation(GL_FUNC_ADD);
	glBlendFunc(GL_ONE, GL_ONE);
	draw_fullscreen_grid_adjacency_border_instanced<512, 512>(1);
	glDisable(GL_BLEND);
	checkGL();

	glUseProgram(0);

	glFinish();

}

int main(int argc, char *argv[]) {
	
	if (argc < 2) {
		log().error() << "Usage: flareon <lensdata>";
		exit(1);
	}
	
	Window *win = createWindow().size(512, 512).title("I Choose You, Flareon!").visible(true);
	win->makeContextCurrent();

	shaderman = new ShaderManager("./res/shader");
	
	win->onResize.attach([](const window_size_event &e) {

		return false;
	});

	// lens-specific precomputation
	precompute(argv[1]);
	
	auto time_fps = chrono::steady_clock::now();
	unsigned fps = 0;
	
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






























