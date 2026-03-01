#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define NON_LINEAR 1


int main(int argc, char *argv[]) {
  FILE *out;
  out = fopen("data.txt", "w");
  printf("usage: ./simulation [u] [A]\n u is the speed in mph\n A is the height in inches\n");
  //car weighs 3500 lbs
  double g, m_t, msmus, m_us, m_s, k_s, b_s, k_t, A, d, u, G,  qs_0, v_s, v_us, v_rel, x, y, dy, v_in, dq_t, F_s, F_d, N, dp_s, dp_us, a_s, a_us, dq_s;

  g = 9.81;
  m_t = 3500 / 2.2; // total mass
  msmus = 5;
  m_us = m_t / (1 + msmus);
  m_s = m_t - m_us;
  k_s = m_s * pow(2 * M_PI * 1, 2.0);
  b_s = 2 * .5 * sqrt(k_s * m_s);
  k_t = m_us * pow(2 * M_PI * 8, 2.0);
  double q_s;


#if NON_LINEAR == 1
   q_s = 3 * m_s * g / k_s;
   qs_0 = 3 * m_s * g / k_s;
   G = pow((k_s / pow((3 * m_s * g), 2.0 / 3.0)), 3.0);
   G = m_s * g / pow(qs_0, 3.0);
   //y =  ( A / 2 ) (1 - cos ( 2 pi x / d) )
#else
   q_s = m_s * g / k_s;
   qs_0 = m_s * g / k_s;
#endif


  A = 2 * 0.0254;
  d = 0.5;
  u = 20 * 0.447038888889;

  char *endptr;
  if (argc >= 3){
  u = strtod(argv[1], &endptr);
  A = strtod(argv[2], &endptr);
  u*= 0.44703888;
  A*=0.0254;
  };
  x = 0.0;



  double t = 0.0;
  double t_f = 1.0;
  double dt = 0.0001;
  double p_s = 0.0, p_us = 0.0, q_t = m_t * g / k_t;
  
  FILE *sim = fopen("sim.txt", "w");
  FILE *SprungMassAccel = fopen("SprunMassAccel.txt", "w");
  FILE *RelSuspDisp = fopen("RelSuspDisp.txt", "w");
  FILE *TireForce = fopen("TireForce.txt", "w");

  while (t < t_f){
  v_s = p_s / m_s;
  v_us = p_us / m_us;
  v_rel = v_us - v_s;
  x = u * t;
  y = 0.5 * A * (1 - cos ( 2 * M_PI * x / d)) * (x < d);
  dy = 0.5 * A * (2 *M_PI/d) * sin( 2 * M_PI * x / d) * (x < d);
  v_in = u * dy;
  dq_s = v_rel;
  dq_t = v_in - v_us;
  F_s = k_s * q_s *(NON_LINEAR == 0) + G * pow(q_s, 3) * (NON_LINEAR == 1);
  F_d = b_s * v_rel;
  N = (k_t * q_t) * (q_t > 0);
  dp_s = -m_s * g + F_s + F_d;
  dp_us = -m_us * g + N - F_s - F_d;
  a_s = dp_s / m_s;
  a_us = dp_us / m_us;

  fprintf(SprungMassAccel, "%.6lf\t%.6lf\n", t, a_s);
  fprintf(RelSuspDisp, "%.6lf\t%.6lf\n", t, q_s - qs_0);
  fprintf(TireForce, "%.6lf\t%.6lf\n", t, N);

  //plot sprung mass acceleration, q_s, tire force N
  q_s = q_s + dq_s * dt;
  q_t = q_t + dq_t * dt;
  p_s = p_s + dp_s * dt;
  p_us = p_us + dp_us * dt;

  t+= dt;
  };


};
