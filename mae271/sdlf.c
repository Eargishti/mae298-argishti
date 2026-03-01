#include <SDL2/SDL.h>
#include <SDL2/SDL_events.h>
#include <SDL2/SDL_render.h>
#include <SDL2/SDL_video.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

double dt = 0.0001;
void endprocesses(SDL_Renderer *renderer, SDL_Window *window1) {
  SDL_DestroyRenderer(renderer);
  SDL_DestroyWindow(window1);
  SDL_Quit();
};

#define WINDOW_WIDTH 900
#define WINDOW_HEIGHT 612
#define WINDOW_x 300
#define WINDOW_y 300

const double k = 1;

int main() {
  int running = 1;
  SDL_Window *window1 = NULL;
  SDL_Renderer *renderer = NULL;
  if (SDL_Init(SDL_INIT_VIDEO) != 0) {
    printf("failed oti nisfdnin");
    return (1);
  };
  window1 = SDL_CreateWindow("Simulation", WINDOW_x, WINDOW_y, WINDOW_WIDTH,
                             WINDOW_HEIGHT, 0);
  if (!window1) {
    printf("failed to open window");
    return (2);
  };
  renderer = SDL_CreateRenderer(window1, -1, 0);
  if (!renderer) {
    printf("failed to open renderer");
    return (3);
  };
  SDL_SetRenderDrawColor(renderer, 0, 255, 0, 255);
  double theta = 0;
  double r = 100;
  int jmax = 1000;
  double p1[3] = {0.5,0.23,0.19};
  double p2[3] = {0.5,5.23,2.9};
  double L = sqrt( pow(p2[0] - p1[0],2.0) + pow(p2[1] - p1[1],2.0) + pow(p2[2] - p1[2],2.0) ); 
  double T = 0.3;
  double camerapos[3] = {-10, 0, 0};
  double cameraFOV = 60;





   for (int j1 = 0; j1 <= jmax; j1++) {
    theta = ((double)j1 / (double)jmax) * 2.0 * M_PI;
    double tempx = -cos(theta);
    double tempy = sin(theta);
    SDL_RenderDrawPointF(renderer, 450 + (int)(r * tempx),
                        300 + (int)(r * tempy));
  }; 
  SDL_RenderPresent(renderer);
  SDL_Event event;
  while (running) {
    while (SDL_PollEvent(&event)) {
      if (event.type == SDL_QUIT) {
        running = 0;
      }
    }
    SDL_Delay(16); // prevents 100% CPU usage
  }
  endprocesses(renderer, window1);
};
