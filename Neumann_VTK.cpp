#include <iostream>
#include <vector>
#include <cmath>
#include <print>
#include <eigen3/Eigen/Core> // Estructura Core de Eigen para vectores dinámicos
#include "rk4.h"
#include "vtk_exporter.h"

void set_initial_conditions(Eigen::VectorXd &s);
void fderiv(const Eigen::VectorXd &s, Eigen::VectorXd &dsdt, double t);
void writer(const Eigen::VectorXd &s, double t);

/*
Compílese con g++ -std=c++23 Neumann_VTK.cpp vtk_exporter.cpp -o neumann.x
Cada paso de tiempo se exporta como un archivo .vtk en la carpeta "output_vtk/",
listo para abrirse en ParaView (File > Open, seleccionar toda la serie frame_*.vtk).

*/

const int Nx = 10;  // Número de nodos en X (incluyendo fronteras)
const int Ny = 10;  // Número de nodos en Y (incluyendo fronteras)
const double Lx = 1.0;
const double Ly = 1.0;

const double dx = Lx / (Nx);
const double dy = Ly / (Ny);

// Condiciones de frontera térmicas en Celsius (°C) o Kelvin (K)
const double T_izq = 0.0; 
const double T_der = 0.0;
const double T_inf = 0.0;
const double T_sup = 100.0;
const double T_inicial = 0.0; // Interior de la placa

// Vector global para almacenar la difusividad térmica espacial alpha(x,y)
Eigen::VectorXd alpha_grid(Nx * Ny);

int main() {
    int total_nodes = Nx * Ny;
    Eigen::VectorXd T(total_nodes);

    // se define una difusividad constante 0.1 para toda la placa
    alpha_grid.setConstant(0.1); 

    // condiciones iniciales
    set_initial_conditions(T);

    // criterio de estabilidad de Courant-Friedrichs-Lewy en 2D

    double max_alpha = alpha_grid.maxCoeff();
    double dt_cfl = 1.0 / (2.0 * max_alpha * ((1.0 / (dx * dx)) + (1.0 / (dy * dy))));
    
    double dt = dt_cfl * 0.8; 
    
    double t_init = 0.0;
    double t_end = 2.0; // Tiempo suficiente para visualizar el comportamiento transitorio

    // Ejecutar la integración numérica con el algoritmo RK4
    integrate_rk4(fderiv, T, t_init, t_end, dt, writer);

    return 0;
}
void set_initial_conditions(Eigen::VectorXd &s) {
    for (int i = 0; i < Ny; ++i) {
        for (int j = 0; j < Nx; ++j) {
            int k = i * Nx + j; // Aplanamiento lineal de la malla indexada

            // Temperatura en las fronteras, condiciones iniciales
            if (i == 0)       s[k] = T_sup;
            else if (i == Ny - 1)  s[k] = T_inf;
            else if (j == 0)       s[k] = T_izq;
            else if (j == Nx - 1)  s[k] = T_der;
            // Interior de la placa
            else                   s[k] = T_inicial;
        }
    }
}

void fderiv(const Eigen::VectorXd &s, Eigen::VectorXd &dsdt, double t) {
    // Inicializamos todas las derivadas en cero.
    dsdt.setZero();

    for (int i = 0; i < Ny; ++i) {
        for (int j = 0; j < Nx; ++j) {
            int k = i * Nx + j;

            // Índices lineales de los 4 vecinos contiguos
            int k_right = i * Nx + j + 1;
            int k_left  = i * Nx + j - 1;
            int k_up    = (i-1) * Nx + j;
            int k_down  = (i+1) * Nx + j;
            
            if(i == 0)  k_up=k_down;
            if(i == Ny-1)  k_down=k_up;
            if(j == 0)  k_left=k_right;
            if(j == Nx-1)  k_right=k_left;

            // Difusividades promedio en las interfaces 

            double alpha_right = 0.5 * (alpha_grid[k] + alpha_grid[k_right]);
            double alpha_left  = 0.5 * (alpha_grid[k] + alpha_grid[k_left]);
            double alpha_up    = 0.5 * (alpha_grid[k] + alpha_grid[k_up]);
            double alpha_down  = 0.5 * (alpha_grid[k] + alpha_grid[k_down]);

            // Flujos aproximados por diferencias finitas de segundo orden
            double flux_x = (alpha_right * (s[k_right] - s[k]) - alpha_left * (s[k] - s[k_left])) / (dx * dx);
            double flux_y = (alpha_up * (s[k_up] - s[k]) - alpha_down * (s[k] - s[k_down])) / (dy * dy);

            // Ecuación diferencial ordinaria acoplada para el nodo k
            dsdt[k] = flux_x + flux_y;
        }
    }
}

void writer(const Eigen::VectorXd &s, double t) {
    // Exportador VTK: se crea una única vez (static) y se reutiliza en cada
    // llamada para ir numerando los frames (frame_000000.vtk, frame_000001.vtk, ...)
    static VTKExporter vtk_exporter("output_Neumann", Nx, Ny, dx, dy);

    vtk_exporter.write_frame(s, t, "temperatura");
}


