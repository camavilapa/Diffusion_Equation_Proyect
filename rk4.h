#pragma once
#include <iostream>

template <class deriv_t, class system_t, class printer_t>
void integrate_rk4(deriv_t fderiv, system_t & s, double tinit, double tend, 
                    double dt, printer_t writer)
{
    // vector to store derivs
    system_t dsdt(s.size()); 
    system_t k1(s.size());
    system_t k2(s.size());
    system_t k3(s.size());
    system_t k4(s.size());
    // time loop
    for(double t = tinit; t < tend; t = t + dt) { // NOTE: Last time step not necessarily tf
       
        //escribir el estado
        writer(s, t);

        fderiv(s, dsdt, t);
        k1 = dsdt;
        fderiv(s + dt*k1/2, dsdt, t+dt/2);
        k2 = dsdt;
        fderiv(s + dt*k2/2, dsdt, t+dt/2);
        k3 = dsdt;
        fderiv(s + dt*k3, dsdt, t+dt);
        k4 = dsdt;

        // compute new state. NOTE: Not using components, assuming valarray or similar
        s = s + dt*(k1 + 2*k2 + 2*k3 + k4)/6.0;



        

      }
}
