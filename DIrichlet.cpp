#include <iostream>
#include <vector>
#include <cmath>
#include <print>
#include <eigen3/Eigen/Core> // Estructura Core de Eigen para vectores dinámicos
#include "rk4.h"

void set_initial_conditions(Eigen::VectorXd &s);
void fderiv(const Eigen::VectorXd &s, Eigen::VectorXd &dsdt, double t);
void writer(const Eigen::VectorXd &s, double t);

/*
En este primer avance se utilizó rk4 para simular matemáticamente la difusión del calor en una placa homogénea con
difusividad constante 0.1 cuando el borde izquierdo se encuentra a 100°C y el resto de la placa a 0°C
*/

const int Nx = 5;  // Número de nodos en X (incluyendo fronteras)
const int Ny = 5;  // Número de nodos en Y (incluyendo fronteras)
const double Lx = 1.0;
const double Ly = 1.0;

const double dx = Lx / (Nx);
const double dy = Ly / (Ny);

// Condiciones de frontera térmicas en Celsius (°C) o Kelvin (K)
const double T_izq = 100.0; 
const double T_der = 0.0;
const double T_inf = 0.0;
const double T_sup = 0.0;
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

    // Iteramos únicamente sobre los nodos internos dejando intactos los bordes
    for (int i = 1; i < Ny - 1; ++i) {
        for (int j = 1; j < Nx - 1; ++j) {
            int k = i * Nx + j;

            // Índices lineales de los 4 vecinos contiguos
            int k_right = i * Nx + j + 1;
            int k_left  = i * Nx + j - 1;
            int k_up    = (i-1) * Nx + j;
            int k_down  = (i+1) * Nx + j;

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
    // Imprimir datos con esta forma [Tiempo] [T_0] [T_1] ... [T_n] separado por tabs
    std::print("{:.2f}",t);
    //imprimos el vector estado s
    for (int i = 0; i < s.size(); ++i) {
        std::print("\t{:.2f}",s[i]);
    }
    std::println("");
}


