/* Implementation of Leap-Frog Method (Mid-point Euler Method)
* The simulation solves a simple, normalized second order differential equation:
* \ddot{y} = -\kappa*sin(y) - \alpha*\dot{y} + \eta, y(0) = 0.4, \dot{y} = 0.4
* where \kappa = 1 and \alpha = 0.05
*/


#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/times.h>
#include<sys/param.h>

#define BUF1 32

int main()
{
   /* Declare variables */
   int     n_dt, n;
   double  fn, t0, tn, y0, yn, ym, kappa, eta;
   double  fm, alpha,  v0, vn, vm, a, b;
   double  t_start, t_stop, dt;
   char    *rundat = "Dat";
   char    datfile[BUF1];
   FILE    *f_dat;

   /* Declaration of right-hand side of differential equation */
   extern double func();

   /* Open the data file Dat  -- leave it as is unless you know what you are doing */
   sprintf(datfile, "%s", rundat);
   if ((f_dat = fopen(datfile, "w")) == NULL) {
      fprintf(stdout, "Can't open Dat file\n");
      fflush(stdout);
      exit(1);
   }

   /* Set system parameters */
   eta   = 0.5;
   kappa = 1.0;
   alpha = 0.05;

   t_start = 0.0;             /* initial time */
   t_stop  = 100.0;           /* final time */
   n_dt    = 1000;            /* number of time steps */
   dt = (t_stop-t_start)/n_dt;   /* time step */
   a = (1.0-0.5*alpha*dt)/(1.0+0.5*alpha*dt);
   b =                1.0/(1.0+0.5*alpha*dt);

   /* Initial condition */
   t0 = t_start; /* initial time */
   y0 = 0.4;     /* initial value */
   v0 = 0.4;     /* initial value */

   int flag = 1;

   double prev = 0;
   double curr = 0;

   /* now dump the initial value and the function to the Dat file */
   fprintf(f_dat, "%8d %.15e %.15e %.15e\n", 0, t0, y0, v0); fflush(f_dat);

   /* Prepare for time-step loop */
   yn = y0;
   fn = func(t0, y0, kappa, eta);
   vn = (a/b)*v0 + (dt*0.5)*fn;

   /* enter the discrete-time loop of the Euler map */
   for (n = 0; n < n_dt; n++) 
   {
      ym = yn + dt*vn; /* y updated */
      tn = t0+n*dt; /*time updated */
      fm = func(tn + dt, ym, kappa, eta); /* new f value for the new (n+1) time */
      vm = a*vn + b*dt*fm; /* v updated */

      fprintf(f_dat, "%8d %.15e %.15e %.15e\n", n+1, tn+dt, ym, vn); fflush(f_dat);

      yn = ym; /* Prepare for next n */
      vn = vm; /* Prepare for next n */
   }
   fclose(f_dat);

   /* Print final values to screen */
   fprintf(stdout, "%8d %.15e %.15e %.15e\n", n, tn+dt, ym, vm);
   fflush(stdout);
}


/* the right-hand side function f */
double func(t,y,c, eta)
double t, y, c, eta;
{
   double f;

   f = -c*sin(y)+eta;

   return f;
}