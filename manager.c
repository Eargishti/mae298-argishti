#include <X11/X.h>
#include <X11/XKBlib.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <err.h>
#include <stdio.h>

static Display *dpy;

static int scr;

static Window root;

#define POSX 200
#define POSY 200
#define WIDTH 400
#define HEIGHT 400
#define BORDER 15

int main() {

  Window win1;
  XEvent ev;

  if ((dpy = XOpenDisplay(NULL)) == NULL) {

    err(1, "Coudln't open display");
  };

  /*Get default screen and root window */

  scr = DefaultScreen(dpy);

  root = RootWindow(dpy, scr);

  /*Creating our simple window */
  win1 = XCreateSimpleWindow(dpy, root, POSX, POSY, WIDTH, HEIGHT, BORDER,
                             BlackPixel(dpy, scr), WhitePixel(dpy, scr));
  /*Map our window to display server */

  XMapWindow(dpy, win1);

  while (XNextEvent(dpy, &ev) == 0) {
  };

  /*unmap our simple window */
  XUnmapWindow(dpy, win1);

  XDestroyWindow(dpy, win1);

  // Close Connection wiht display server

  XCloseDisplay(dpy);

  printf("");

  return 0;
};
