/*=================================================================================
	Name		: Gravity Forward Modeling
	Version		: 1.0
	Programmer	: Ir. Indra Gunawan, M.Sc.
	Institution	: Geophysical Engineering, Institut Teknologi Bandung
=================================================================================*/

#ifndef GFM_H
#define GFM_H
 
#include <gsl/gsl_integration.h>
#include <boost/format.hpp>
#include "vikidtype.h"

#define GRAV 6.67384E-11

double rijk(double ai, double bj, double jk)
{
   double res = sqrt(pow(ai,2) + pow(bj,2) + pow(jk,2));
   return res;
}


double sorokinForward(Body* bodies, int nBodies, XYZ station)
{
   double res = 0;
   double a[2], b[2], c[2], x[2], y[2], z[2];
   double rho;
   for (int n=0; n<nBodies; n++)
   {
      double result =0;
      a[0] = bodies[n].x;
      b[0] = bodies[n].y;
      c[0] = bodies[n].z;
      a[1] = bodies[n].x + bodies[n].dx;
      b[1] = bodies[n].y + bodies[n].dy;
      c[1] = bodies[n].z + bodies[n].dz;
      for (int i=0; i<2; i++)
      {
         for (int j=0; j<2; j++)
         {
            for (int k=0; k<2; k++)
            {
               x[i] = station.x - a[i];
               y[j] = station.y - b[j];
               z[k] = station.z - c[k];
               double distance = rijk(x[i], y[j], z[k]);
               double nu = pow(-1,(i+1)) * pow(-1,(j+1)) * pow(-1,(k+1));
               result = result + nu * (x[i] * log(y[j] + distance) + y[j] * (log(x[i] + distance)) + z[k] * atan2((z[k] * distance),(x[i] * y[j])));
            }
         }
      }
      res = res + (-1 * GRAV * bodies[n].rho * result); 
   }
   res = res * 1E5;
   return res;   
}

double plouffForward(Body* bodies, int nBodies, XYZ station)
{
   double res = 0;
   double a[2], b[2], c[2], x[2], y[2], z[2];
   double rho;
   for (int n=0; n<nBodies; n++)
   {
      double result =0;
      a[0] = bodies[n].x;
      b[0] = bodies[n].y;
      c[0] = bodies[n].z;
      a[1] = bodies[n].x + bodies[n].dx;
      b[1] = bodies[n].y + bodies[n].dy;
      c[1] = bodies[n].z + bodies[n].dz;
      for (int i=0; i<2; i++)
      {
         for (int j=0; j<2; j++)
         {
            for (int k=0; k<2; k++)
            {
               x[i] = station.x - a[i];
               y[j] = station.y - b[j];
               z[k] = station.z - c[k];
               double distance = rijk((x[i]), (y[j]), (z[k]));
               double nu = pow(-1,(i+1)) * pow(-1,(j+1)) * pow(-1,(k+1));
               result = result + nu * (x[i] * log(y[j] + distance) + y[j] * (log(x[i] + distance)) - z[k] * atan2((x[i] * y[j]),(z[k] * distance)));
            }
         }
      }
      res = res + (-1 * GRAV * bodies[n].rho * result); 
   }
   res = res * 1E5;
   return res;   
}

double nagyForward(Body* bodies, int nBodies, XYZ station)
{
   double res = 0;
   double a[2], b[2], c[2], x[2], y[2], z[2];
   double rho;
   for (int n=0; n<nBodies; n++)
   {
      double result =0;
      a[0] = bodies[n].x;
      b[0] = bodies[n].y;
      c[0] = bodies[n].z;
      a[1] = bodies[n].x + bodies[n].dx;
      b[1] = bodies[n].y + bodies[n].dy;
      c[1] = bodies[n].z + bodies[n].dz;
      for (int i=0; i<2; i++)
      {
         for (int j=0; j<2; j++)
         {
            for (int k=0; k<2; k++)
            {
               x[i] = station.x - a[i];
               y[j] = station.y - b[j];
               z[k] = station.z - c[k];
               double distance = rijk(x[i], y[j], z[k]);
               double nu = pow(-1,(i+1)) * pow(-1,(j+1)) * pow(-1,(k+1));
               result = result + nu * (x[i] * log(y[j] + distance) + y[j] * (log(x[i] + distance)) - z[k] * asin((pow(y[j],2) + pow(z[k],2) + y[j] * distance)/((y[j] + distance)* sqrt(pow(y[j],2) + pow(z[k],2)))));
            }
         }
      }
      res = res + (-1 * GRAV * bodies[n].rho * result); 
   }
   res = res * 1E5;
   return res;   
}

double okabeForward(Body* bodies, int nBodies, XYZ station)
{
   double res = 0;
   double a[2], b[2], c[2], x[2], y[2], z[2];
   double rho;
   for (int n=0; n<nBodies; n++)
   {
      double result =0;
      a[0] = bodies[n].x;
      b[0] = bodies[n].y;
      c[0] = bodies[n].z;
      a[1] = bodies[n].x + bodies[n].dx;
      b[1] = bodies[n].y + bodies[n].dy;
      c[1] = bodies[n].z + bodies[n].dz;
      for (int i=0; i<2; i++)
      {
         for (int j=0; j<2; j++)
         {
            for (int k=0; k<2; k++)
            {
               x[i] = station.x - a[i];
               y[j] = station.y - b[j];
               z[k] = station.z - c[k];
               double distance = rijk(x[i], y[j], z[k]);
               double nu = pow(-1,(i+1)) * pow(-1,(j+1)) * pow(-1,(k+1));
               result = result + nu * (x[i] * log(y[j] + distance) + y[j] * (log(x[i] + distance)) + 2 * z[k] * atan2((x[i] + y[j] + distance),(z[k])));
            }
         }
      }
      res = res + (-1 * GRAV * bodies[n].rho * result); 
   }
   res = res * 1E5;
   return res;   
}

#endif
