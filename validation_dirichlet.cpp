#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <numbers>
#include <eigen3/Eigen/Core>
#include "rk4.h"
#include "vtk_exporter.h"

// Compilar con: g++ -std=c++23 validacion_Dirichlet.cpp vtk_exporter.cpp -o validacion.x

// --- PARÁMETROS CONFIGURABLES DE LA SIMULACIÓN ---
const int Nx = 50;  // Subir a 50 o 100 para alta resolución
const int Ny = 50;
const double Lx = 1.0;
const double Ly = 1.0;

const double dx = Lx / (Nx - 1);
const double dy = Ly / (Ny - 1);

const double T_izq = 0.0; 
const double T_der = 0.0;
const double T_inf = 0.0;
const double T_sup = 100.0;
const double T_inicial = 0.0;

Eigen::VectorXd alpha_grid(Nx * Ny);
std::ofstream archivo_error("error_vs_tiempo.csv");

// --- SOLUCIÓN ANALÍTICA TRANSITORIA POR SERIES DE FOURIER ---
double T_analitica_transitoria(double x, double y, double t, double alpha, int N_terms = 40) {
    const double pi = std::numbers::pi;
    
    // 1. Componente Estacionaria
    double T_est = 0.0;
    for (int n = 1; n <= N_terms; n += 2) {
        double factor = (4.0 * T_sup) / (n * pi);
        double num_y = std::sinh((n * pi * y) / Lx);
        double den_y = std::sinh((n * pi * Ly) / Lx);
        T_est += factor * (num_y / den_y) * std::sin((n * pi * x) / Lx);
    }

    if (t <= 0.0) return 0.0;

    // 2. Componente Transitoria
    double T_trans = 0.0;
    for (int m = 1; m <= N_terms; m += 2) {
        for (int n = 1; n <= N_terms; n += 2) {
            double C_mn = (8.0 * T_sup) / (n * m * pi * pi);
            double lambda_sq = pi * pi * ((n * n) / (Lx * Lx) + (m * m) / (Ly * Ly));
            double spatial = std::sin((n * pi * x) / Lx) * std::sin((m * pi * y) / Ly);
            double temporal = std::exp(-alpha * lambda_sq * t);

            T_trans += C_mn * spatial * temporal;
        }
    }
    return T_est - T_trans;
}

void set_initial_conditions(Eigen::VectorXd &s) {
    for (int i = 0; i < Ny; ++i) {
        for (int j = 0; j < Nx; ++j) {
            int k = i * Nx + j;
            if (i == 0)            s[k] = T_sup;
            else if (i == Ny - 1)  s[k] = T_inf;
            else if (j == 0)       s[k] = T_izq;
            else if (j == Nx - 1)  s[k] = T_der;
            else                   s[k] = T_inicial;
        }
    }
}

// Derivada limpia: Solo nodos interiores (Dirichlet Puro)
void fderiv(const Eigen::VectorXd &s, Eigen::VectorXd &dsdt, double t) {
    dsdt.setZero(); // Mantiene bordes fijos en dsdt = 0

    for (int i = 1; i < Ny - 1; ++i) {
        for (int j = 1; j < Nx - 1; ++j) {
            int k = i * Nx + j;

            int k_right = i * Nx + j + 1;
            int k_left  = i * Nx + j - 1;
            int k_up    = (i - 1) * Nx + j;
            int k_down  = (i + 1) * Nx + j;

            double alpha_right = 0.5 * (alpha_grid[k] + alpha_grid[k_right]);
            double alpha_left  = 0.5 * (alpha_grid[k] + alpha_grid[k_left]);
            double alpha_up    = 0.5 * (alpha_grid[k] + alpha_grid[k_up]);
            double alpha_down  = 0.5 * (alpha_grid[k] + alpha_grid[k_down]);

            double flux_x = (alpha_right * (s[k_right] - s[k]) - alpha_left * (s[k] - s[k_left])) / (dx * dx);
            double flux_y = (alpha_up * (s[k_up] - s[k]) - alpha_down * (s[k] - s[k_down])) / (dy * dy);

            dsdt[k] = flux_x + flux_y;
        }
    }
}

void writer(const Eigen::VectorXd &s, double t) {
    static VTKExporter vtk_exporter("output_validacion", Nx, Ny, dx, dy);
    
    // Exportar VTK cada N iteraciones para no saturar memoria
    static int frame = 0;
    if (frame % 5 == 0) {
        vtk_exporter.write_frame(s, t, "temperatura");
    }
    frame++;

    // Calcular y guardar errores en el CSV
    double suma_error_sq = 0.0;
    double error_max = 0.0;
    int nodos_interiores = 0;
    double alpha = alpha_grid.maxCoeff();

    for (int i = 1; i < Ny - 1; ++i) {
        double y = Ly - i * dy; // Coordenada Y
        for (int j = 1; j < Nx - 1; ++j) {
            double x = j * dx;  // Coordenada X
            int k = i * Nx + j;

            double T_exacta = T_analitica_transitoria(x, y, t, alpha);
            double diff = std::abs(s[k] - T_exacta);

            if (diff > error_max) error_max = diff;
            suma_error_sq += diff * diff;
            nodos_interiores++;
        }
    }

    double error_rms = std::sqrt(suma_error_sq / nodos_interiores);
    archivo_error << t << "," << error_rms << "," << error_max << "\n";
}

int main() {
    int total_nodes = Nx * Ny;
    Eigen::VectorXd T(total_nodes);

    alpha_grid.setConstant(0.1); 
    set_initial_conditions(T);

    double max_alpha = alpha_grid.maxCoeff();
    double dt_cfl = 1.0 / (2.0 * max_alpha * ((1.0 / (dx * dx)) + (1.0 / (dy * dy))));
    double dt = dt_cfl * 0.8; 
    
    double t_init = 0.0;
    double t_end = 2.0;

    archivo_error << "tiempo,error_rms,error_max\n";

    integrate_rk4(fderiv, T, t_init, t_end, dt, writer);

    archivo_error.close();
    std::cout << "Simulación finalizada. Datos guardados en 'error_vs_tiempo.csv'\n";
    return 0;
}
