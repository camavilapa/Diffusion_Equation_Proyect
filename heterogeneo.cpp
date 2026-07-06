#include <iostream>
#include <vector>
#include <cmath>
#include <eigen3/Eigen/Core> // Usamos la estructura Core de Eigen para vectores dinámicos
#include "rk4.h"

void writer(const Eigen::VectorXd &s, double t);
void set_initial_conditions(Eigen::VectorXd &s);
void fderiv(const Eigen::VectorXd &s, Eigen::VectorXd &dsdt, double t);

// ==========================================
// CONSTANTES COORDENADAS Y CONFIGURACIÓN DE MALLA
// ==========================================
const int Nx = 20;  // Número de nodos en X (incluyendo fronteras)
const int Ny = 20;  // Número de nodos en Y (incluyendo fronteras)
const double Lx = 1.0;
const double Ly = 1.0;

const double dx = Lx / (Nx - 1);
const double dy = Ly / (Ny - 1);

// Condiciones de frontera térmicas (Dirichlet)
const double T_izq = 100.0;
const double T_der = 50.0;
const double T_inf = 0.0;
const double T_sup = 20.0;
const double T_inicial = 25.0; // Interior de la placa

// Vector global para almacenar la difusividad térmica espacial alpha(x,y)
Eigen::VectorXd alpha_grid(Nx * Ny);

// ==========================================
// FUNCIÓN DE CONDICIONES INICIALES
// ==========================================
void set_initial_conditions(Eigen::VectorXd &s) {
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            int k = i * Ny + j; // Aplanamiento lineal de la malla indexada

            // Imponer valores en las fronteras de la placa
            if (i == 0)       s[k] = T_izq;
            else if (i == Nx - 1)  s[k] = T_der;
            else if (j == 0)       s[k] = T_inf;
            else if (j == Ny - 1)  s[k] = T_sup;
            // Interior de la placa
            else                   s[k] = T_inicial;
        }
    }
}

// ==========================================
// FUNCIÓN DERIVADA (DIFUSIÓN DEL CALOR EN 2D)
// ==========================================
void fderiv(const Eigen::VectorXd &s, Eigen::VectorXd &dsdt, double t) {
    // Inicializamos todas las derivadas en cero.
    // Esto asegura automáticamente que los nodos frontera mantengan derivada 0 y queden fijos.
    dsdt.setZero();

    // Iteramos únicamente sobre los nodos internos (dejando intactos los bordes)
    for (int i = 1; i < Nx - 1; ++i) {
        for (int j = 1; j < Ny - 1; ++j) {
            int k = i * Ny + j;

            // Índices lineales de los 4 vecinos contiguos
            int k_right = (i + 1) * Ny + j;
            int k_left  = (i - 1) * Ny + j;
            int k_up    = i * Ny + (j + 1);
            int k_down  = i * Ny + (j - 1);

            // Difusividades promedio en las interfaces (Forma Conservativa Heterogénea)
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

// ==========================================
// FUNCIÓN IMPRESORA (WRITER)
// ==========================================
void writer(const Eigen::VectorXd &s, double t) {
    // Imprime en consola: [Tiempo] [T_0] [T_1] ... [T_n] separado por espacios
    // Ideal para ser redirigido a un archivo .txt y leído fácilmente en Python con NumPy
    std::cout << t;
    for (int i = 0; i < s.size(); ++i) {
        std::cout << " " << s[i];
    }
    std::cout << "\n";
}

// ==========================================
// PROGRAMA PRINCIPAL
// ==========================================
int main() {
    int total_nodes = Nx * Ny;
    Eigen::VectorXd T(total_nodes);

    // --- Definición del Material (Ejemplo Heterogéneo Inicial) ---
    // Inicializamos con un valor base homogéneo
    alpha_grid.setConstant(0.1); 
    
    // Agregamos una imperfección / barrera aislante en el centro de la placa (Caso Heterogéneo)
    for (int i = Nx / 3; i < 2 * Nx / 3; ++i) {
        for (int j = Ny / 3; j < 2 * Ny / 3; ++j) {
            alpha_grid[i * Ny + j] = 0.005; // Difusividad críticamente baja
        }
    }

    // Configurar el estado de temperaturas inicial
    set_initial_conditions(T);

    // --- Control de Estabilidad Numérica (Criterio CFL de 2D) ---
    double max_alpha = alpha_grid.maxCoeff();
    double dt_cfl = 1.0 / (2.0 * max_alpha * ((1.0 / (dx * dx)) + (1.0 / (dy * dy))));
    
    // Tomamos un paso de tiempo seguro, por debajo del límite crítico de convergencia
    double dt = dt_cfl * 0.8; 
    
    double t_init = 0.0;
    double t_end = 2.0; // Tiempo suficiente para visualizar el comportamiento transitorio

    // Ejecutar la integración numérica con el algoritmo RK4 plantilla
    integrate_rk4(fderiv, T, t_init, t_end, dt, writer);

    return 0;
}