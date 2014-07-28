


#include <iostream>

#include "Initial3D.hpp"
#include "Flareon.hpp"
#include "Log.hpp"
#include "Window.hpp"

using namespace std;
using namespace ambition;
using namespace initial3d;

void display() {

}

int main(int argc, char *argv[]) {

	Window *win = createWindow().size(512, 512).title("I Choose You, Flareon!").visible(true);
	win->makeContextCurrent();

	win->onResize.attach([](const window_size_event &e) {

		return false;
	});

	while (!win->shouldClose()) {
		glfwPollEvents();

		display();

		glFinish();
		win->swapBuffers();
	}

	delete win;
	glfwTerminate();

	return 0;
}