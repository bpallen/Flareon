

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

struct lens_interface {
	// sphere radius
	// +ve => front-convex, -ve => front-concave
	// 0 => this isnt a lens
	double sr;
	// z-position of this surface on optical axis
	double z;
	// refractive indices in front of and behind this surface
	double n0, n2;
	// aperture radius
	// 0 => sensor plane
	double ap;
};

// lens interfaces, in order from entry to sensor
vector<lens_interface> interfaces;

void precompute(const string &path) {
	// parse the lens specification file
	ifstream ifs(path);
	if (!ifs.is_open()) {
		throw runtime_error("unable to open file '" + path + "'");
	}
	
	double n_last = 1.0;
	while (ifs.good()) {
		string line;
		getline(ifs, line);
		istringstream iss(line);
		string t0;
		iss >> t0;
		if (!iss.good()) continue;
		if (t0[0] == '#') continue;
		iss.seekg(0);
		double sr, dz, n, ap;
		iss >> sr >> dz >> n >> ap;
		if (iss.fail()) {
			log("LensData").warning() << "Encountered invalid line '" << line << "' in file '" << path << "'";
		} else {
			lens_interface li;
			li.sr = sr;
			li.z = 0;
			li.n0 = n_last;
			li.n2 = n;
			li.ap = ap;
			n_last = n;
			for (auto it = interfaces.begin(); it != interfaces.end(); it++) {
				it->z += dz;
			}
			interfaces.push_back(li);
		}
	}
	
	if (interfaces.size() == 0) {
		log("LensData").error() << "Lens data is empty";
		throw runtime_error("lens data is empty");
	}
	
	if (interfaces.back().ap > 0) {
		log("LensData").warning() << "Lens configuration has no sensor plane";
	}
	
	// ???
	
}

void upload_uniforms(GLuint prog) {
	glUniform1ui(glGetUniformLocation(prog, "num_interfaces"), interfaces.size());
	glUniform1ui(glGetUniformLocation(prog, "num_ghosts"), 1);
	for (unsigned i = 0; i < interfaces.size(); i++) {
		// TODO
	}
}

// width and height are in terms of on-screen triangles
template <unsigned Width, unsigned Height>
void draw_fullscreen_grid_adjacency() {
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
			if (i >= 0 && i <= Width + 1 && j >= 0 && j <= Height + 1) {
				ii = 1 + y_vals.size() * i + j;
			}
			indices.push_back(unsigned(ii));
		};
		for (unsigned i = 1; i < Width + 1; i++) {
			for (unsigned j = 1; j < Height + 1; j++) {
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
	glDrawElements(GL_TRIANGLES_ADJACENCY, Width * Height * 12, GL_UNSIGNED_INT, nullptr);
	glBindVertexArray(0);
}

void display(const size2i &sz) {

	glViewport(0, 0, sz.w, sz.h);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	GLuint prog_flare = shaderman->getProgram("flare.glsl");

	glUseProgram(prog_flare);
	upload_uniforms(prog_flare);

	checkGL();
	
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	
	// one pass per wavelength and blend?
	
	draw_fullscreen_grid_adjacency<32, 32>();

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

	auto time_fps = chrono::steady_clock::now();
	unsigned fps = 0;
	
	// lens-specific precomputation
	precompute(argv[1]);
	
	while (!win->shouldClose()) {
		auto now = chrono::steady_clock::now();
		if (now - time_fps > chrono::seconds(1)) {
			time_fps = now;
			char cbuf[100];
			sprintf(cbuf, "Flareon [%d FPS @%dx%d", fps, win->width(), win->height());
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






























